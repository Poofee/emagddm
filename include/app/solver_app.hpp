/**
 * @file solver_app.hpp
 * @brief 应用层 - 求解器应用模块头文件
 * @details 提供FETI-DP求解器的应用接口
 *          核心流程：输入文件 → ProjectManager(内部自动选择加载器) → 求解
 * @author Poofee
 * @date 2026-04-04
 * @version 5.0
 */

#pragma once

#include <string>
#include "tool/project_manager.hpp"

namespace app {

/**
 * @class SolverApp
 * @brief 求解器应用类
 * 
 * 核心设计：
 * - 格式无关：SolverApp 不关心文件格式，只调用 ProjectManager::openProject()
 * - 加载器封装：IProjectLoader 是 ProjectManager 的内部实现细节
 * - 职责清晰：SolverApp 只负责初始化日志、调用 openProject、执行求解流程
 */
class SolverApp {
public:
    SolverApp();
    ~SolverApp();

    SolverApp(const SolverApp&) = delete;
    SolverApp& operator=(const SolverApp&) = delete;

    /**
     * @brief 初始化求解器应用
     * 
     * 将输入文件交给 ProjectManager，由其内部自动选择对应的格式加载器。
     * 
     * @param input_file 输入文件路径（.aedt/.json/.xml）
     * @return bool 初始化成功返回true
     */
    bool initialize(const std::string& input_file);

    /**
     * @brief 运行求解器
     * 
     * 使用已加载到 ProjectManager 中的数据执行完整求解流程：
     * 1. 网格生成/导入
     * 2. 组装有限元矩阵
     * 3. FETI-DP迭代求解
     * 4. 输出结果
     * 
     * @return bool 求解成功返回true
     */
    bool run();

    /**
     * @brief 清理资源
     */
    void cleanup();

    /**
     * @brief 获取项目管理器引用（用于访问加载的数据）
     * @return ProjectManager& 项目管理器引用
     */
    tool::ProjectManager& getProjectManager() { return tool::ProjectManager::getInstance(); }

private:
    std::string input_file_;
    bool initialized_;
};

} // namespace app
