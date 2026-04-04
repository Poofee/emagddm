/**
 * @file solver_app.cpp
 * @brief 应用层 - 求解器应用模块实现
 * @details 核心流程：输入文件 → ProjectManager（内部自动选择加载器） → 求解
 * @author Poofee
 * @date 2026-04-04
 * @version 5.0
 */

#include "solver_app.hpp"
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
    FEEM_INFO("========================================");

    return true;
}

bool SolverApp::run() {
    if (!initialized_) {
        FEEM_ERROR("求解器应用未初始化，无法运行");
        return false;
    }

    FEEM_INFO("========================================");
    FEEM_INFO("开始FETI-DP电磁场求解");
    FEEM_INFO("========================================");

    // 步骤1：根据已加载的项目数据进行网格生成/导入
    FEEM_INFO("");
    FEEM_INFO("[步骤1] 网格生成/导入...");
    // TODO: 使用 project_manager_ 中的几何和网格设置生成有限元网格

    // 步骤2：组装有限元矩阵
    FEEM_INFO("[步骤2] 组装刚度矩阵...");
    // TODO: 遍历 project_manager_ 的材料、边界条件等，组装全局刚度矩阵 K 和载荷向量 F

    // 步骤3：执行FETI-DP迭代求解
    FEEM_INFO("[步骤3] FETI-DP迭代求解...");
    // TODO: 子域分解、界面预处理、共轭梯度迭代

    // 步骤4：输出求解结果
    FEEM_INFO("[步骤4] 输出求解结果...");
    // TODO: 计算场量（磁通密度B、磁场强度H等），写入结果文件

    FEEM_INFO("");
    FEEM_INFO("========================================");
    FEEM_INFO("求解完成（核心求解流程待实现）");
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
