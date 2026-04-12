/**
 * @file solver_scheduler.cpp
 * @brief 求解器层 - 求解器调度器实现
 * @details 实现求解流程的统一调度与管理，
 *          连接应用层(SolverApp)与物理场求解器的桥梁。
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "solver_scheduler.hpp"
#include "electrostatic_solver.hpp"
#include "dc_conduction_solver.hpp"
#include "project_manager.hpp"
#include "logger_factory.hpp"
#include "math_constants.hpp"

namespace solver {

// ============================================================
// 构造/析构
// ============================================================

SolverScheduler::SolverScheduler()
    : physics_field_(nullptr)
    , dim_type_(tool::DimType::D3)
    , sim_type_(tool::SimulationType::ELECTROSTATIC)
{
    FEEM_DEBUG("SolverScheduler 构造完成");
}

SolverScheduler::~SolverScheduler() {
    clear();
}

// ============================================================
// 初始化
// ============================================================

bool SolverScheduler::initialize(tool::ProjectManager& pm) {
    FEEM_INFO("========================================");
    FEEM_INFO("SolverScheduler 初始化");
    FEEM_INFO("========================================");

    // 检查网格数据是否存在
    if (!pm.hasEMMeshData()) {
        FEEM_ERROR("ProjectManager中没有网格拓扑数据，请先导入或生成网格");
        return false;
    }

    auto* mesh_data = pm.getEMMeshData();
    if (!mesh_data || mesh_data->isEmpty()) {
        FEEM_ERROR("网格拓扑数据为空");
        return false;
    }

    // 提取维度类型和求解类型
    dim_type_ = pm.getDesignType();

    // 优先从SolutionSetup获取求解类型，其次使用ProjectManager默认值
    sim_type_ = pm.getSolutionType();
    const auto& setups = pm.getAllSolutionSetups();
    if (!setups.empty()) {
        auto first_setup = setups.begin()->second;
        if (first_setup) {
            sim_type_ = first_setup->getSolutionType();
        }
    }

    // 打印初始化摘要
    FEEM_INFO("  维度类型: {}", static_cast<int>(dim_type_));
    FEEM_INFO("  求解类型: {}", static_cast<int>(sim_type_));
    FEEM_INFO("  网格节点数: {}", mesh_data->getNodeCount());
    FEEM_INFO("  网格单元数: {}", mesh_data->getElementCount());
    FEEM_INFO("  材料数量: {}", pm.getAllMaterials().size());
    FEEM_INFO("  边界条件数量: {}", pm.getAllBoundaries().size());
    FEEM_INFO("  激励源数量: {}", pm.getAllExcitations().size());

    FEEM_INFO("SolverScheduler 初始化成功");
    return true;
}

// ============================================================
// 执行求解
// ============================================================

bool SolverScheduler::run() {
    if (!physics_field_) {
        FEEM_INFO("未设置自定义物理场求解器，自动创建...");

        physics_field_ = createPhysicsField(sim_type_, dim_type_);
        if (!physics_field_) {
            FEEM_ERROR("无法创建物理场求解器（类型: SimulationType={}, DimType={})",
                       static_cast<int>(sim_type_), static_cast<int>(dim_type_));
            return false;
        }
        FEEM_INFO("自动创建物理场求解器: {}", physics_field_->getSolverName());
    }

    auto& pm = tool::ProjectManager::getInstance();

    // 获取网格数据指针
    auto* mesh_data = pm.getEMMeshData();
    if (!mesh_data) {
        FEEM_ERROR("网格数据指针为空");
        return false;
    }

    // 构建材料参数映射表
    auto materials = buildMaterialMap(pm);
    FEEM_INFO("材料映射表构建完成，包含 {} 种材料", materials.size());

    // 构建边界条件列表
    auto boundaries = buildBoundaryList(pm);
    FEEM_INFO("边界条件列表构建完成，包含 {} 个边界条件", boundaries.size());

    // 构建激励源列表
    auto excitations = buildExcitationList(pm);
    FEEM_INFO("激励源列表构建完成，包含 {} 个激励源", excitations.size());

    // ========== 阶段1：Setup（前处理+矩阵装配）==========
    FEEM_INFO("");
    FEEM_INFO("[阶段1] 物理场求解器 Setup...");
    if (!physics_field_->setup(*mesh_data, materials, boundaries, excitations)) {
        FEEM_ERROR("物理场求解器 Setup 失败");
        return false;
    }
    FEEM_INFO("[阶段1] Setup 完成");

    // ========== 阶段2：Solve（线性求解）==========
    FEEM_INFO("[阶段2] 物理场求解器 Solve...");
    if (!physics_field_->solve()) {
        FEEM_ERROR("物理场求解器 Solve 失败");
        return false;
    }
    FEEM_INFO("[阶段2] Solve 完成");

    // ========== 阶段3：PostProcess（后处理）==========
    FEEM_INFO("[阶段3] 物理场求解器 PostProcess...");
    if (!physics_field_->postProcess()) {
        FEEM_ERROR("物理场求解器 PostProcess 失败");
        return false;
    }
    FEEM_INFO("[阶段3] PostProcess 完成");

    // 缓存场数据引用
    field_data_ = physics_field_->getFieldData();

    FEEM_INFO("");
    FEEM_INFO("========================================");
    FEEM_INFO("求解流程执行成功！");
    FEEM_INFO("  求解器: {}", physics_field_->getSolverName());
    FEEM_INFO("  场数据状态: {}", field_data_.hasData() ? "有效" : "无效");
    FEEM_INFO("========================================");

    return true;
}

// ============================================================
// 设置自定义物理场求解器
// ============================================================

void SolverScheduler::setPhysicsField(std::unique_ptr<PhysicsField> field) {
    physics_field_ = std::move(field);
    if (physics_field_) {
        FEEM_INFO("已设置自定义物理场求解器: {}", physics_field_->getSolverName());
    }
}

// ============================================================
// 获取场数据
// ============================================================

const FieldData& SolverScheduler::getFieldData() const {
    return field_data_;
}

const PhysicsField* SolverScheduler::getPhysicsField() const {
    return physics_field_.get();
}

// ============================================================
// 清空资源
// ============================================================

void SolverScheduler::clear() {
    if (physics_field_) {
        physics_field_->clear();
        physics_field_.reset();
    }
    field_data_.clear();
    dim_type_ = tool::DimType::D3;
    sim_type_ = tool::SimulationType::ELECTROSTATIC;
    FEEM_DEBUG("SolverScheduler 资源已释放");
}

// ============================================================
// 私有方法：创建物理场求解器（工厂模式）
// ============================================================

std::unique_ptr<PhysicsField> SolverScheduler::createPhysicsField(
    tool::SimulationType sim_type,
    tool::DimType dim_type
) const {
    switch (sim_type) {
        case tool::SimulationType::ELECTROSTATIC:
            FEEM_INFO("创建静电场求解器 (ElectrostaticSolver)");
            return std::make_unique<ElectrostaticSolver>(dim_type);

        case tool::SimulationType::DC_CONDUCTION:
            FEEM_INFO("创建直流电流场求解器 (DCConductionSolver)");
            return std::make_unique<DCConductionSolver>(dim_type);

        case tool::SimulationType::MAGNETOSTATIC:
            FEEM_ERROR("静磁场求解器 (MagnetostaticSolver) 尚未实现");
            return nullptr;

        case tool::SimulationType::EDDYCURRENT:
            FEEM_ERROR("涡流场求解器 (EddyCurrentSolver) 尚未实现");
            return nullptr;

        case tool::SimulationType::TRANSIENT:
            FEEM_ERROR("瞬态场求解器 (TransientSolver) 尚未实现");
            return nullptr;

        default:
            FEEM_ERROR("不支持的求解类型: {}", static_cast<int>(sim_type));
            return nullptr;
    }
}

// ============================================================
// 私有方法：构建材料参数映射表
// ============================================================

std::map<int, numeric::MaterialProperties> SolverScheduler::buildMaterialMap(
    tool::ProjectManager& pm
) {
    std::map<int, numeric::MaterialProperties> material_map;

    const auto& materials = pm.getAllMaterials();
    for (const auto& [name, mat_ptr] : materials) {
        if (!mat_ptr) continue;

        // 以Material的ID作为映射键（转换为int）
        int material_id = static_cast<int>(mat_ptr->getID());

        // 构建MaterialProperties结构体
        numeric::MaterialProperties props;

        // 介电常数 ε = ε_r × ε₀
        double permittivity = mat_ptr->getPermittivity();
        if (permittivity > 0.0) {
            props.epsilon = permittivity * numeric::EPSILON0;
        } else {
            // 默认使用真空介电常数
            props.epsilon = numeric::EPSILON0;
            FEEM_DEBUG("材料 [{}] 未指定介电常数，使用真空值 ε₀",
                       mat_ptr->getName());
        }

        // 磁导率 μ = μ_r × μ₀
        props.mu = mat_ptr->getRelativePermeability() * numeric::MU0;

        // 电导率 σ
        props.sigma = mat_ptr->getConductivity();

        material_map[material_id] = props;

        FEEM_DEBUG("材料映射: ID={}, name={}, ε={:.4e}, μ={:.4e}, σ={:.4e}",
                   material_id, mat_ptr->getName(),
                   props.epsilon, props.mu, props.sigma);
    }

    return material_map;
}

// ============================================================
// 私有方法：构建边界条件列表
// ============================================================

std::vector<tool::Boundary> SolverScheduler::buildBoundaryList(
    tool::ProjectManager& pm
) {
    std::vector<tool::Boundary> boundary_list;

    const auto& boundaries = pm.getAllBoundaries();
    boundary_list.reserve(boundaries.size());

    for (const auto& [name, bc_ptr] : boundaries) {
        if (bc_ptr) {
            boundary_list.push_back(*bc_ptr);
        }
    }

    return boundary_list;
}

// ============================================================
// 私有方法：构建激励源列表
// ============================================================

std::vector<tool::Excitation> SolverScheduler::buildExcitationList(
    tool::ProjectManager& pm
) {
    std::vector<tool::Excitation> excitation_list;

    const auto& excitations = pm.getAllExcitations();
    excitation_list.reserve(excitations.size());

    for (const auto& [name, exc_ptr] : excitations) {
        if (exc_ptr) {
            excitation_list.push_back(*exc_ptr);
        }
    }

    return excitation_list;
}

} // namespace solver
