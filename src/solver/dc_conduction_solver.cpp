/**
 * @file dc_conduction_solver.cpp
 * @brief 求解器层 - 直流电流场求解器实现文件
 * @details 实现DCConductionSolver类的完整求解流程：
 *          setup() → solve() → postProcess()
 *          包含前处理检查、DOF分配、矩阵装配、边界条件施加、
 *          线性求解、解向量回扩、后处理计算（电场/电流密度/焦耳热）等核心步骤。
 *
 *          与ElectrostaticSolver的核心差异：
 *          1. preCheck(): 增加电导率σ>0合法性校验（绝缘体σ=0合法，仅警告）
 *          2. buildMaterialMapping(): 构建材料列表供DCConductionIntegrator直接读取sigma字段
 *          3. assembleSystem(): 使用DCConductionIntegrator直接装配（独立于EMAssembly管道）
 *          4. postProcess(): 计算电流密度J=σE和焦耳热功率P=φ^T K φ（替代静电能量）
 *          5. getSimulationType(): 返回DC_CONDUCTION
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "dc_conduction_solver.hpp"
#include "logger_factory.hpp"
#include "em_enums.hpp"

#include <algorithm>
#include <set>

namespace solver {

// ==================== 构造与析构 ====================

DCConductionSolver::DCConductionSolver()
    : dim_type_(tool::DimType::D3)
    , setup_done_(false)
    , solve_done_(false)
    , dof_manager_(nullptr)
    , assembly_(nullptr)
    , linear_solver_(nullptr)
    , mesh_data_ptr_(nullptr)
{
}

DCConductionSolver::DCConductionSolver(tool::DimType dim_type)
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

bool DCConductionSolver::setup(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    const std::vector<tool::Excitation>& excitations
) {
    FEEM_INFO("DCConductionSolver::setup 开始初始化直流电流场求解器");

    /* 1. 缓存输入数据 */
    mesh_data_ptr_ = &mesh_data;
    materials_ = materials;

    /* 2. 前处理检查 */
    if (!preCheck()) {
        FEEM_ERROR("DCConductionSolver::setup 前处理检查失败");
        return false;
    }

    /* 3. 构建材料映射（关键：将sigma复制到epsilon以复用装配管道） */
    buildMaterialMapping(materials);

    /* 4. 处理激励：电压源→Dirichlet BC，电流源→Neumann BC */
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
        FEEM_ERROR("DCConductionSolver::setup DOF分配失败");
        return false;
    }

    /* 7. 全局矩阵装配 */
    if (!assembleSystem()) {
        FEEM_ERROR("DCConductionSolver::setup 矩阵装配失败");
        return false;
    }

    /* 8. 施加边界条件和激励 */
    if (!applyBoundaryAndExcitation(all_boundaries, excitations)) {
        FEEM_ERROR("DCConductionSolver::setup 边界条件施加失败");
        return false;
    }

    setup_done_ = true;
    FEEM_INFO("DCConductionSolver::setup 初始化完成");
    return true;
}

bool DCConductionSolver::solve() {
    if (!setup_done_) {
        FEEM_ERROR("DCConductionSolver::solve setup未完成，无法求解");
        return false;
    }

    FEEM_INFO("DCConductionSolver::solve 开始线性求解");

    /* 1. 创建线性求解器（直流电流场K为SPD矩阵） */
    linear_solver_ = createLinearSolver();
    if (!linear_solver_) {
        FEEM_ERROR("DCConductionSolver::solve 线性求解器创建失败");
        return false;
    }

    /* 2. 获取刚度矩阵和修正后的右端项 */
    const numeric::CsrMatrix<double>& K_csr = getStiffnessMatrix();
    /* f_free_在applyBoundaryAndExcitation中存储了施加边界条件后的右端项向量 */
    const Eigen::VectorXd& F = f_free_;

    /* 3. 设置系数矩阵（触发分解） */
    linear_solver_->set_matrix(K_csr);

    /* 4. 求解线性系统 K·V_free = F */
    numeric::SolverResult result = linear_solver_->solve(F);

    /* 5. 检查求解状态 */
    if (result.status != numeric::SolverStatus::SUCCESS) {
        FEEM_ERROR("DCConductionSolver::solve 线性求解失败, 状态: {}",
                   static_cast<int>(result.status));
        if (!result.error_msg.empty()) {
            FEEM_ERROR("  错误信息: {}", result.error_msg);
        }
        return false;
    }

    FEEM_INFO("DCConductionSolver::solve 线性求解成功, 耗时: {:.2f} ms",
              result.solve_time_ms);

    /* 6. 保存自由DOF解向量 */
    phi_free_ = result.x;

    /* 7. 回扩到完整解向量（包含约束DOF） */
    phi_full_ = bc_manager_.expandSolution(phi_free_, *dof_manager_);

    solve_done_ = true;
    FEEM_INFO("DCConductionSolver::solve 求解完成, 自由DOF数: {}, 总DOF数: {}",
              phi_free_.size(), phi_full_.size());
    return true;
}

bool DCConductionSolver::postProcess() {
    if (!solve_done_) {
        FEEM_ERROR("DCConductionSolver::postProcess solve未完成，无法后处理");
        return false;
    }

    if (!mesh_data_ptr_) {
        FEEM_ERROR("DCConductionSolver::postProcess 网格数据无效");
        return false;
    }

    FEEM_INFO("DCConductionSolver::postProcess 开始后处理计算");

    /* 1. 设置节点电位（电压分布） */
    field_data_.setNodalPotential(phi_full_, static_cast<int>(phi_full_.size()));

    /* 2. 计算电场强度 E = -∇V（复用FieldData现有方法） */
    field_data_.computeElectricField(*mesh_data_ptr_, materials_);

    /* 3. 计算电流密度 J = σE（直流电流场核心派生量） */
    field_data_.computeCurrentDensity(*mesh_data_ptr_, materials_);

    /* 4. 计算总焦耳热功率 P = V^T K V（替代静电能量的W=0.5*φ^T K φ） */
    const numeric::CsrMatrix<double>& K_csr = getStiffnessMatrix();
    double power = field_data_.computeJouleHeating(K_csr, phi_free_);
    FEEM_INFO("DCConductionSolver::postProcess 总焦耳热功率: {:.6e} W", power);

    FEEM_INFO("DCConductionSolver::postProcess 后处理完成");
    return true;
}

const FieldData& DCConductionSolver::getFieldData() const {
    return field_data_;
}

tool::SimulationType DCConductionSolver::getSimulationType() const {
    return tool::SimulationType::DC_CONDUCTION;
}

std::string DCConductionSolver::getSolverName() const {
    return "DCConductionSolver";
}

void DCConductionSolver::clear() {
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

    /* 清理直接装配数据 */
    direct_K_csr_.reset();
    direct_F_.resize(0);

    mesh_data_ptr_ = nullptr;
}

// ==================== 直流电流场专用接口 ====================

void DCConductionSolver::setDimType(tool::DimType dim_type) {
    dim_type_ = dim_type;
}

tool::DimType DCConductionSolver::getDimType() const {
    return dim_type_;
}

const fe_em::EMDOFManager& DCConductionSolver::getDOFManager() const {
    return *dof_manager_;
}

const numeric::CsrMatrix<double>& DCConductionSolver::getStiffnessMatrix() const {
    /* 优先返回DCConductionIntegrator直接装配的K矩阵（绕过EMAssembly的DummyEMIntegrator） */
    if (direct_K_csr_ && direct_K_csr_->rows() > 0) {
        return *direct_K_csr_;
    }
    /* Fallback：如果直接装配数据不可用，尝试从EMAssembly获取 */
    if (assembly_) {
        return assembly_->getStiffnessMatrix();
    }
    static const numeric::CsrMatrix<double> empty_k(0, 0);
    return empty_k;
}

const Eigen::VectorXd& DCConductionSolver::getSourceVector() const {
    /* 优先返回直接装配的源项向量 */
    if (direct_F_.size() > 0) {
        return direct_F_;
    }
    /* Fallback：如果直接数据不可用，尝试从EMAssembly获取 */
    if (assembly_) {
        return assembly_->getSourceVector();
    }
    static const Eigen::VectorXd empty_f = Eigen::VectorXd();
    return empty_f;
}

const Eigen::VectorXd& DCConductionSolver::getFullSolution() const {
    return phi_full_;
}

// ==================== 私有方法实现 ====================

bool DCConductionSolver::preCheck() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("DCConductionSolver::preCheck 网格数据指针为空");
        return false;
    }

    const fe_em::EMMeshData& mesh = *mesh_data_ptr_;

    /* 检查1：网格非空 */
    if (mesh.isEmpty()) {
        FEEM_ERROR("DCConductionSolver::preCheck 网格数据为空（无节点和单元）");
        return false;
    }
    if (mesh.nodes.empty()) {
        FEEM_ERROR("DCConductionSolver::preCheck 网格节点列表为空");
        return false;
    }
    if (mesh.elements.empty()) {
        FEEM_ERROR("DCConductionSolver::preCheck 网格单元列表为空");
        return false;
    }

    /* 检查2：所有单元的material_id有效 */
    std::set<int> material_ids;
    for (const auto& elem : mesh.elements) {
        if (elem.material_id < 0) {
            FEEM_ERROR("DCConductionSolver::preCheck 单元{}的material_id={}无效（负数）",
                       elem.id, elem.material_id);
            return false;
        }
        material_ids.insert(elem.material_id);
    }

    /* 检查3：材料映射覆盖所有单元 */
    for (int mid : material_ids) {
        if (materials_.find(mid) == materials_.end()) {
            FEEM_WARN("DCConductionSolver::preCheck material_id={}未在材料映射中找到，"
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
        FEEM_WARN("DCConductionSolver::preCheck 未检测到Dirichlet边界条件，"
                  "系统矩阵可能奇异导致求解失败");
    }

    /* 检查5：导体区域必须有正电导率（新增：直流电流场特有校验） */
    for (const auto& [mid, props] : materials_) {
        if (props.sigma <= 0.0) {
            FEEM_WARN("DCConductionSolver::preCheck material_id={} 电导率σ={}<=0, "
                      "可能导致奇异矩阵（绝缘体σ=0合法但需配合适当边界条件）",
                      mid, props.sigma);
        }
    }

    FEEM_INFO("DCConductionSolver::preCheck 前处理检查通过: "
              "节点数={}, 单元数={}, 材料数={}, 边界标记数={}",
              mesh.getNodeCount(), mesh.getElementCount(),
              material_ids.size(), mesh.getBoundaryMarkerCount());
    return true;
}

void DCConductionSolver::buildMaterialMapping(
    const std::map<int, numeric::MaterialProperties>& materials
) {
    int max_material_id = 0;
    for (const auto& [mid, _] : materials) {
        max_material_id = std::max(max_material_id, mid);
    }

    if (mesh_data_ptr_) {
        for (const auto& elem : mesh_data_ptr_->elements) {
            max_material_id = std::max(max_material_id, elem.material_id);
        }
    }

    /* 构建按material_id索引的材料列表，保留原始ε和σ值。
     * DCConductionIntegrator 通过 getEffectiveSigma() 直接读取 sigma 字段，
     * 无需任何字段替换或映射操作。 */
    material_list_.resize(static_cast<size_t>(max_material_id) + 1);
    for (auto& [mid, props] : materials) {
        material_list_[static_cast<size_t>(mid)] = props;
        FEEM_DEBUG("DCConductionSolver: material_id={}, ε={}, σ={}",
                   mid, props.epsilon, props.sigma);
    }

    FEEM_DEBUG("DCConductionSolver::buildMaterialMapping 材料映射构建完成, "
               "最大material_id={}, 材料数={}",
               max_material_id, materials.size());
}

bool DCConductionSolver::allocateDOFs() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("DCConductionSolver::allocateDOFs 网格数据无效");
        return false;
    }

    FEEM_INFO("DCConductionSolver::allocateDOFs 开始DOF分配");

    /* 直流电流场为纯标量场（SCALAR_ONLY），无需棱边DOF */
    std::vector<std::vector<int>> empty_edge_mapping;
    dof_manager_ = std::make_unique<fe_em::EMDOFManager>(
        *mesh_data_ptr_, empty_edge_mapping, nullptr
    );

    /* 执行完整DOF编号流程 */
    dof_manager_->build();

    /* 获取DOF映射表 */
    elem_l2g_ = dof_manager_->getElemLocalToGlobal();

    int num_free = dof_manager_->getNumFreeDOFs();
    FEEM_INFO("DCConductionSolver::allocateDOFs DOF分配完成, "
              "自由DOF数={}", num_free);
    return true;
}

bool DCConductionSolver::assembleSystem() {
    if (!mesh_data_ptr_) {
        FEEM_ERROR("DCConductionSolver::assembleSystem 网格数据无效");
        return false;
    }

    FEEM_INFO("DCConductionSolver::assembleSystem 开始全局矩阵装配（使用DCConductionIntegrator直接装配）");

    /* 1. 确定矩阵维度 */
    int free_dofs = dof_manager_ ? dof_manager_->getNumFreeDOFs() : 0;
    if (free_dofs <= 0) {
        FEEM_ERROR("DCConductionSolver::assembleSystem 无效的自由DOF数（<= 0）");
        return false;
    }

    /* 2. 创建COO格式稀疏矩阵用于单元级装配 */
    auto K_coo = std::make_unique<numeric::CooMatrix<double>>(free_dofs, free_dofs);
    direct_F_ = Eigen::VectorXd::Zero(free_dofs);

    /* 3. 单元循环：使用DCConductionIntegrator计算 K_e = ∫∇N^T·σ·∇NdΩ */
    int processed_elements = 0;
    int skipped_elements = 0;

    for (size_t e = 0; e < mesh_data_ptr_->elements.size(); ++e) {
        const auto& elem = mesh_data_ptr_->elements[e];

        try {
            solver::DCConductionIntegrator integrator(elem.elem_type);

            if (elem.material_id >= 0 &&
                elem.material_id < static_cast<int>(material_list_.size())) {
                integrator.setMaterialProperties(material_list_[elem.material_id]);
            }

            int num_nodes = static_cast<int>(elem.node_ids.size());
            Eigen::MatrixXd node_coords(3, num_nodes);
            for (int n = 0; n < num_nodes; ++n) {
                int nid = elem.node_ids[n];
                if (nid >= 0 && nid < static_cast<int>(mesh_data_ptr_->nodes.size())) {
                    node_coords(0, n) = mesh_data_ptr_->nodes[nid].x;
                    node_coords(1, n) = mesh_data_ptr_->nodes[nid].y;
                    node_coords(2, n) = mesh_data_ptr_->nodes[nid].z;
                }
            }

            Eigen::MatrixXd K_e = integrator.computeStiffnessMatrix(node_coords);

            const auto& l2g = elem_l2g_[e];
            for (int i = 0; i < K_e.rows(); ++i) {
                int global_i = l2g.indices[i];
                if (global_i < 0 || global_i >= free_dofs) continue;
                for (int j = 0; j < K_e.cols(); ++j) {
                    int global_j = l2g.indices[j];
                    if (global_j < 0 || global_j >= free_dofs) continue;
                    double val = K_e(i, j);
                    if (std::abs(val) > 1e-15) {
                        K_coo->add_value(global_i, global_j, val);
                    }
                }
            }

            processed_elements++;
        } catch (const std::exception& ex) {
            FEEM_WARN("assembleSystem: 单元[{}]装配异常: {} - 跳过",
                     elem.id, ex.what());
            skipped_elements++;
        }
    }

    /* 4. COO→CSR转换 */
    direct_K_csr_ = std::make_unique<numeric::CsrMatrix<double>>(free_dofs, free_dofs);
    direct_K_csr_->build_from_coo(*K_coo);

    assembly_ = std::make_unique<numeric::EMAssembly>();

    FEEM_INFO("DCConductionSolver::assembleSystem 矩阵装配完成: "
              "自由DOF={}, 单元数={}（成功={}, 跳过={}), K非零元素={}",
              free_dofs,
              mesh_data_ptr_->elements.size(), processed_elements, skipped_elements,
              direct_K_csr_->nnz());
    return true;
}

bool DCConductionSolver::applyBoundaryAndExcitation(
    const std::vector<tool::Boundary>& boundaries,
    const std::vector<tool::Excitation>& excitations
) {
    if (!mesh_data_ptr_ || !dof_manager_) {
        FEEM_ERROR("DCConductionSolver::applyBoundaryAndExcitation "
                   "网格数据或DOF管理器无效");
        return false;
    }

    FEEM_INFO("DCConductionSolver::applyBoundaryAndExcitation "
              "开始施加边界条件和激励");

    /* excitations已在setup()中通过excitation_manager_处理，
       此处保留参数以备后续扩展（如非线性迭代中重新施加激励） */
    (void)excitations;

    /* 1. 提取Dirichlet约束DOF */
    if (!bc_manager_.extractDirichletDOFs(boundaries, *mesh_data_ptr_, *dof_manager_)) {
        FEEM_ERROR("DCConductionSolver::applyBoundaryAndExcitation "
                   "Dirichlet DOF提取失败");
        return false;
    }

    if (!bc_manager_.hasDirichletBC()) {
        FEEM_WARN("DCConductionSolver::applyBoundaryAndExcitation "
                  "未检测到Dirichlet边界条件，系统可能奇异");
    }

    /* 2. 获取右端项向量的可修改副本（使用直接装配的源项向量，而非EMAssembly默认空向量） */
    Eigen::VectorXd rhs = getSourceVector();

    /* 3. 施加Dirichlet边界条件修正（消去法：F_free -= K_fc * V_c） */
    bc_manager_.applyDirichletToRHS(
        rhs, *mesh_data_ptr_, elem_l2g_, material_list_
    );

    /* 4. 施加体电荷密度源项（直流电流场通常无此项，保留接口兼容性） */
    const std::map<int, double>& charge_density =
        excitation_manager_.getVolumeChargeDensity();
    if (!charge_density.empty()) {
        excitation_manager_.applyVolumeChargeSource(
            rhs, *mesh_data_ptr_, elem_l2g_, charge_density
        );
    }

    /* 5. 保存修正后的右端项向量供solve()使用 */
    f_free_ = rhs;

    FEEM_INFO("DCConductionSolver::applyBoundaryAndExcitation "
              "边界条件施加完成, Dirichlet约束DOF数={}",
              bc_manager_.getDirichletDOFs().size());
    return true;
}

std::unique_ptr<numeric::EMLinearSolverBase>
DCConductionSolver::createLinearSolver() {
    /* 直流电流场刚度矩阵K为对称正定(SPD)矩阵，
       使用EMSolverFactory创建专用的对称直接求解器（Cholesky分解） */
    FEEM_INFO("DCConductionSolver::createLinearSolver "
              "创建SPD对称直接求解器（Cholesky分解）");
    return numeric::EMSolverFactory::create_solver_for_electrostatic();
}

} // namespace solver
