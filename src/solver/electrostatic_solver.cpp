/**
 * @file electrostatic_solver.cpp
 * @brief 求解器层 - 静电场求解器实现文件
 * @details 实现ElectrostaticSolver类的完整求解流程：
 *          setup() → solve() → postProcess()
 *          包含前处理检查、DOF分配、矩阵装配、边界条件施加、
 *          线性求解、解向量回扩、后处理计算等核心步骤。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "electrostatic_solver.hpp"
#include "logger_factory.hpp"
#include "em_enums.hpp"

#include <algorithm>
#include <set>

namespace solver {

// ==================== 构造与析构 ====================

ElectrostaticSolver::ElectrostaticSolver()
    : dim_type_(tool::DimType::D3)
    , setup_done_(false)
    , solve_done_(false)
    , dof_manager_(nullptr)
    , assembly_(nullptr)
    , linear_solver_(nullptr)
    , mesh_data_ptr_(nullptr)
{
}

ElectrostaticSolver::ElectrostaticSolver(tool::DimType dim_type)
    : dim_type_(dim_type)
    , setup_done_(false)
    , solve_done_(false)
    , dof_manager_(nullptr)
    , assembly_(nullptr)
    , linear_solver_(nullptr)
    , mesh_data_ptr_(nullptr)
{
}

// ==================== PhysicsField接口实现 ====================

bool ElectrostaticSolver::setup(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    const std::vector<tool::Excitation>& excitations
) {
    FEEM_INFO("ElectrostaticSolver::setup 开始初始化静电场求解器");

    /* 1. 缓存输入数据 */
    mesh_data_ptr_ = &mesh_data;
    materials_ = materials;

    /* 2. 前处理检查 */
    if (!preCheck()) {
        FEEM_ERROR("ElectrostaticSolver::setup 前处理检查失败");
        return false;
    }

    /* 3. 构建材料映射 */
    buildMaterialMapping(materials);

    /* 4. 处理激励：电压源→Dirichlet BC，体电荷密度→源项 */
    std::vector<tool::Boundary> excitation_bcs =
        excitation_manager_.processExcitations(excitations, mesh_data);

    /* 5. 合并原始边界条件与激励生成的等效边界条件 */
    std::vector<tool::Boundary> all_boundaries = boundaries;
    all_boundaries.insert(
        all_boundaries.end(),
        excitation_bcs.begin(),
        excitation_bcs.end()
    );

    /* 6. DOF分配 */
    if (!allocateDOFs()) {
        FEEM_ERROR("ElectrostaticSolver::setup DOF分配失败");
        return false;
    }

    /* 7. 全局矩阵装配 */
    if (!assembleSystem()) {
        FEEM_ERROR("ElectrostaticSolver::setup 矩阵装配失败");
        return false;
    }

    /* 8. 施加边界条件和激励 */
    if (!applyBoundaryAndExcitation(all_boundaries, excitations)) {
        FEEM_ERROR("ElectrostaticSolver::setup 边界条件施加失败");
        return false;
    }

    setup_done_ = true;
    FEEM_INFO("ElectrostaticSolver::setup 初始化完成");
    return true;
}

bool ElectrostaticSolver::solve() {
    if (!setup_done_) {
        FEEM_ERROR("ElectrostaticSolver::solve setup未完成，无法求解");
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::solve 开始线性求解");

    /* 1. 创建线性求解器（静电场K为SPD矩阵） */
    linear_solver_ = createLinearSolver();
    if (!linear_solver_) {
        FEEM_ERROR("ElectrostaticSolver::solve 线性求解器创建失败");
        return false;
    }

    /* 2. 获取刚度矩阵和修正后的右端项 */
    const numeric::CsrMatrix<double>& K_csr = assembly_->getStiffnessMatrix();
    /* f_free_在applyBoundaryAndExcitation中存储了施加边界条件后的右端项向量 */
    const Eigen::VectorXd& F = f_free_;

    /* 3. 设置系数矩阵（触发分解） */
    linear_solver_->set_matrix(K_csr);

    /* 4. 求解线性系统 K·phi_free = F */
    numeric::SolverResult result = linear_solver_->solve(F);

    /* 5. 检查求解状态 */
    if (result.status != numeric::SolverStatus::SUCCESS) {
        FEEM_ERROR("ElectrostaticSolver::solve 线性求解失败, 状态: {}",
                   static_cast<int>(result.status));
        if (!result.error_msg.empty()) {
            FEEM_ERROR("  错误信息: {}", result.error_msg);
        }
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::solve 线性求解成功, 耗时: {:.2f} ms",
              result.solve_time_ms);

    /* 6. 保存自由DOF解向量 */
    phi_free_ = result.x;

    /* 7. 回扩到完整解向量（包含约束DOF） */
    phi_full_ = bc_manager_.expandSolution(phi_free_, *dof_manager_);

    solve_done_ = true;
    FEEM_INFO("ElectrostaticSolver::solve 求解完成, 自由DOF数: {}, 总DOF数: {}",
              phi_free_.size(), phi_full_.size());
    return true;
}

bool ElectrostaticSolver::postProcess() {
    if (!solve_done_) {
        FEEM_ERROR("ElectrostaticSolver::postProcess solve未完成，无法后处理");
        return false;
    }

    if (!mesh_data_ptr_) {
        FEEM_ERROR("ElectrostaticSolver::postProcess 网格数据无效");
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::postProcess 开始后处理计算");

    /* 1. 设置节点电位 */
    field_data_.setNodalPotential(phi_full_, static_cast<int>(phi_full_.size()));

    /* 2. 计算电场强度 E = -∇φ */
    field_data_.computeElectricField(*mesh_data_ptr_, materials_);

    /* 3. 计算电位移矢量 D = εE */
    field_data_.computeElectricDisplacement(*mesh_data_ptr_, materials_);

    /* 4. 计算静电能量 W = 0.5 * φ^T * K * φ */
    const numeric::CsrMatrix<double>& K_csr = assembly_->getStiffnessMatrix();
    double energy = field_data_.computeElectrostaticEnergy(K_csr, phi_free_);
    FEEM_INFO("ElectrostaticSolver::postProcess 静电能量: {:.6e} J", energy);

    FEEM_INFO("ElectrostaticSolver::postProcess 后处理完成");
    return true;
}

const FieldData& ElectrostaticSolver::getFieldData() const {
    return field_data_;
}

tool::SimulationType ElectrostaticSolver::getSimulationType() const {
    return tool::SimulationType::ELECTROSTATIC;
}

std::string ElectrostaticSolver::getSolverName() const {
    return "ElectrostaticSolver";
}

void ElectrostaticSolver::clear() {
    setup_done_ = false;
    solve_done_ = false;

    dof_manager_.reset();
    assembly_.reset();
    linear_solver_.reset();
    bc_manager_.clear();
    excitation_manager_.clear();
    field_data_.clear();

    phi_free_.resize(0);
    phi_full_.resize(0);
    f_free_.resize(0);
    elem_l2g_.clear();
    material_list_.clear();
    materials_.clear();
    mesh_data_ptr_ = nullptr;
}

// ==================== 静电场专用接口 ====================

void ElectrostaticSolver::setDimType(tool::DimType dim_type) {
    dim_type_ = dim_type;
}

tool::DimType ElectrostaticSolver::getDimType() const {
    return dim_type_;
}

const fe_em::EMDOFManager& ElectrostaticSolver::getDOFManager() const {
    return *dof_manager_;
}

const numeric::CsrMatrix<double>& ElectrostaticSolver::getStiffnessMatrix() const {
    return assembly_->getStiffnessMatrix();
}

const Eigen::VectorXd& ElectrostaticSolver::getSourceVector() const {
    return assembly_->getSourceVector();
}

const Eigen::VectorXd& ElectrostaticSolver::getFullSolution() const {
    return phi_full_;
}

// ==================== 私有方法实现 ====================

bool ElectrostaticSolver::preCheck() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("ElectrostaticSolver::preCheck 网格数据指针为空");
        return false;
    }

    const fe_em::EMMeshData& mesh = *mesh_data_ptr_;

    /* 检查1：网格非空 */
    if (mesh.isEmpty()) {
        FEEM_ERROR("ElectrostaticSolver::preCheck 网格数据为空（无节点和单元）");
        return false;
    }
    if (mesh.nodes.empty()) {
        FEEM_ERROR("ElectrostaticSolver::preCheck 网格节点列表为空");
        return false;
    }
    if (mesh.elements.empty()) {
        FEEM_ERROR("ElectrostaticSolver::preCheck 网格单元列表为空");
        return false;
    }

    /* 检查2：所有单元的material_id有效 */
    std::set<int> material_ids;
    for (const auto& elem : mesh.elements) {
        if (elem.material_id < 0) {
            FEEM_ERROR("ElectrostaticSolver::preCheck 单元{}的material_id={}无效（负数）",
                       elem.id, elem.material_id);
            return false;
        }
        material_ids.insert(elem.material_id);
    }

    /* 检查3：材料映射覆盖所有单元 */
    for (int mid : material_ids) {
        if (materials_.find(mid) == materials_.end()) {
            FEEM_WARN("ElectrostaticSolver::preCheck material_id={}未在材料映射中找到，"
                      "将使用真空默认参数", mid);
        }
    }

    /* 检查4：至少存在一个Dirichlet边界条件 */
    bool has_dirichlet = false;
    for (const auto& marker : mesh.boundary_markers) {
        if (marker.bnd_type == tool::BndType::DIRICHLET ||
            marker.bnd_type == tool::BndType::PERFECT_E ||
            marker.bnd_type == tool::BndType::ODD_SYMMETRY) {
            has_dirichlet = true;
            break;
        }
    }
    if (!has_dirichlet) {
        FEEM_WARN("ElectrostaticSolver::preCheck 未检测到Dirichlet边界条件，"
                  "系统矩阵可能奇异导致求解失败");
    }

    FEEM_INFO("ElectrostaticSolver::preCheck 前处理检查通过: "
              "节点数={}, 单元数={}, 材料数={}, 边界标记数={}",
              mesh.getNodeCount(), mesh.getElementCount(),
              material_ids.size(), mesh.getBoundaryMarkerCount());
    return true;
}

void ElectrostaticSolver::buildMaterialMapping(
    const std::map<int, numeric::MaterialProperties>& materials
) {
    /* 确定material_id的最大值，构建按material_id索引的向量 */
    int max_material_id = 0;
    for (const auto& [mid, _] : materials) {
        max_material_id = std::max(max_material_id, mid);
    }

    /* 同时检查网格中引用的material_id */
    if (mesh_data_ptr_) {
        for (const auto& elem : mesh_data_ptr_->elements) {
            max_material_id = std::max(max_material_id, elem.material_id);
        }
    }

    /* 构建向量：material_list_[material_id] = MaterialProperties */
    material_list_.resize(static_cast<size_t>(max_material_id) + 1);
    for (const auto& [mid, props] : materials) {
        material_list_[static_cast<size_t>(mid)] = props;
    }

    FEEM_DEBUG("ElectrostaticSolver::buildMaterialMapping 材料映射构建完成, "
               "最大material_id={}, 材料数={}",
               max_material_id, materials.size());
}

bool ElectrostaticSolver::allocateDOFs() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("ElectrostaticSolver::allocateDOFs 网格数据无效");
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::allocateDOFs 开始DOF分配");

    /* 静电场为纯标量场（SCALAR_ONLY），无需棱边DOF */
    std::vector<std::vector<int>> empty_edge_mapping;
    dof_manager_ = std::make_unique<fe_em::EMDOFManager>(
        *mesh_data_ptr_, empty_edge_mapping, nullptr
    );

    /* 执行完整DOF编号流程 */
    dof_manager_->build();

    /* 获取DOF映射表 */
    elem_l2g_ = dof_manager_->getElemLocalToGlobal();

    int num_free = dof_manager_->getNumFreeDOFs();
    FEEM_INFO("ElectrostaticSolver::allocateDOFs DOF分配完成, "
              "自由DOF数={}", num_free);
    return true;
}

bool ElectrostaticSolver::assembleSystem() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("ElectrostaticSolver::assembleSystem 网格数据无效");
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::assembleSystem 开始全局矩阵装配");

    assembly_ = std::make_unique<numeric::EMAssembly>();

    /* 静电场为静态问题，is_transient=false */
    bool success = assembly_->assemble(
        *mesh_data_ptr_,
        elem_l2g_,
        material_list_,
        false
    );

    if (!success) {
        numeric::AssemblyStats stats = assembly_->getStats();
        FEEM_ERROR("ElectrostaticSolver::assembleSystem 矩阵装配失败: {}",
                   stats.error_message);
        return false;
    }

    numeric::AssemblyStats stats = assembly_->getStats();
    FEEM_INFO("ElectrostaticSolver::assembleSystem 矩阵装配完成: "
              "耗时={:.2f}ms, 自由DOF={}, 约束DOF={}, "
              "K非零元素={}, K对称={}",
              stats.assembly_time_ms, stats.free_dofs, stats.constrained_dofs,
              stats.nnz_stiffness, stats.is_symmetric_k ? "是" : "否");
    return true;
}

bool ElectrostaticSolver::applyBoundaryAndExcitation(
    const std::vector<tool::Boundary>& boundaries,
    const std::vector<tool::Excitation>& excitations
) {
    if (!mesh_data_ptr_ || !dof_manager_) {
        FEEM_ERROR("ElectrostaticSolver::applyBoundaryAndExcitation "
                   "网格数据或DOF管理器无效");
        return false;
    }

    FEEM_INFO("ElectrostaticSolver::applyBoundaryAndExcitation "
              "开始施加边界条件和激励");

    /* excitations已在setup()中通过excitation_manager_处理，
       此处保留参数以备后续扩展（如非线性迭代中重新施加激励） */
    (void)excitations;

    /* 1. 提取Dirichlet约束DOF */
    if (!bc_manager_.extractDirichletDOFs(boundaries, *mesh_data_ptr_, *dof_manager_)) {
        FEEM_ERROR("ElectrostaticSolver::applyBoundaryAndExcitation "
                   "Dirichlet DOF提取失败");
        return false;
    }

    if (!bc_manager_.hasDirichletBC()) {
        FEEM_WARN("ElectrostaticSolver::applyBoundaryAndExcitation "
                  "未检测到Dirichlet边界条件，系统可能奇异");
    }

    /* 2. 获取右端项向量的可修改副本 */
    Eigen::VectorXd rhs = assembly_->getSourceVector();

    /* 3. 施加Dirichlet边界条件修正（消去法：F_free -= K_fc * phi_c） */
    bc_manager_.applyDirichletToRHS(
        rhs, *mesh_data_ptr_, elem_l2g_, material_list_
    );

    /* 4. 施加体电荷密度源项 */
    const std::map<int, double>& charge_density =
        excitation_manager_.getVolumeChargeDensity();
    if (!charge_density.empty()) {
        excitation_manager_.applyVolumeChargeSource(
            rhs, *mesh_data_ptr_, elem_l2g_, charge_density
        );
    }

    /* 5. 保存修正后的右端项向量供solve()使用 */
    f_free_ = rhs;

    FEEM_INFO("ElectrostaticSolver::applyBoundaryAndExcitation "
              "边界条件施加完成, Dirichlet约束DOF数={}",
              bc_manager_.getDirichletDOFs().size());
    return true;
}

std::unique_ptr<numeric::EMLinearSolverBase>
ElectrostaticSolver::createLinearSolver() {
    /* 静电场刚度矩阵K为对称正定(SPD)矩阵，
       使用EMSolverFactory创建专用的对称直接求解器（Cholesky分解） */
    FEEM_INFO("ElectrostaticSolver::createLinearSolver "
              "创建SPD对称直接求解器（Cholesky分解）");
    return numeric::EMSolverFactory::create_solver_for_electrostatic();
}

} // namespace solver
