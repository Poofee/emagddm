/**
 * @file solver_app.cpp
 * @brief 应用层 - 求解器应用模块实现
 * @details 核心流程：输入文件 → ProjectManager（内部自动选择加载器） → SolverScheduler → 求解
 * @author Poofee
 * @date 2026-04-04
 * @version 6.0
 */

#include "solver_app.hpp"
#include "solver_scheduler.hpp"
#include "field_data.hpp"
#include "logger_factory.hpp"

namespace app {

SolverApp::SolverApp() : initialized_(false) {
}

SolverApp::~SolverApp() {
    cleanup();
}

bool SolverApp::initialize(const std::string& input_file) {
    if (input_file.empty()) {
        FEEM_ERROR("输入文件路径为空");
        return false;
    }

    input_file_ = input_file;

    auto& pm = tool::ProjectManager::getInstance();

    FEEM_INFO("========================================");
    FEEM_INFO("正在打开项目文件: {}", input_file);
    FEEM_INFO("========================================");

    if (!pm.openProject(input_file)) {
        FEEM_ERROR("项目打开失败: {}", pm.getLastError());
        return false;
    }

    initialized_ = true;

    FEEM_INFO("");
    FEEM_INFO("========================================");
    FEEM_INFO("数据加载完成摘要:");
    FEEM_INFO("  项目名称: {}", pm.getProjectName());
    FEEM_INFO("  材料数量: {}", pm.getAllMaterials().size());
    FEEM_INFO("  边界条件数量: {}", pm.getAllBoundaries().size());
    FEEM_INFO("  激励源数量: {}", pm.getAllExcitations().size());
    FEEM_INFO("  求解设置数量: {}", pm.getAllSolutionSetups().size());

    if (pm.hasEMMeshData()) {
        auto* mesh = pm.getEMMeshData();
        if (mesh) {
            FEEM_INFO("  网格节点数: {}", mesh->getNodeCount());
            FEEM_INFO("  网格单元数: {}", mesh->getElementCount());
        }
    } else {
        FEEM_WARN("  网格拓扑数据: 未加载");
    }

    FEEM_INFO("========================================");

    return true;
}

bool SolverApp::run() {
    if (!initialized_) {
        FEEM_ERROR("求解器应用未初始化，无法运行");
        return false;
    }

    FEEM_INFO("========================================");
    FEEM_INFO("开始电磁场求解");
    FEEM_INFO("========================================");

    // 创建求解器调度器
    solver::SolverScheduler scheduler;

    // 初始化调度器（从ProjectManager提取全部求解数据）
    if (!scheduler.initialize(tool::ProjectManager::getInstance())) {
        FEEM_ERROR("求解器调度器初始化失败");
        return false;
    }

    // 执行完整求解流程（setup → solve → postProcess）
    if (!scheduler.run()) {
        FEEM_ERROR("求解流程执行失败");
        return false;
    }

    // 获取求解结果
    const auto& field_data = scheduler.getFieldData();

    // 输出结果摘要
    FEEM_INFO("");
    FEEM_INFO("[步骤4] 求解结果输出:");

    if (field_data.hasData()) {
        const auto& phi = field_data.getNodalPotential();
        FEEM_INFO("  节点电位范围: [{:.6e}, {:.6e}] V",
                  phi.minCoeff(), phi.maxCoeff());

        if (!field_data.getElementElectricField().empty()) {
            double max_E = 0.0;
            for (const auto& E : field_data.getElementElectricField()) {
                double E_mag = E.norm();
                if (E_mag > max_E) max_E = E_mag;
            }
            FEEM_INFO("  最大电场强度: {:.6e} V/m", max_E);
        }

        // TODO: 根据SolutionSetup中的输出配置导出VTK/CSV文件
        // field_data.exportVTK(output_path, *mesh_data);
        // field_data.exportCSV(output_path, *mesh_data);
    } else {
        FEEM_WARN("  场数据为空，无有效结果可输出");
    }

    FEEM_INFO("");
    FEEM_INFO("========================================");
    FEEM_INFO("求解完成！");
    FEEM_INFO("========================================");

    return true;
}

void SolverApp::cleanup() {
    if (initialized_) {
        FEEM_INFO("清理求解器应用资源");
        tool::ProjectManager::getInstance().closeProject();
        initialized_ = false;
    }
}

} // namespace app
