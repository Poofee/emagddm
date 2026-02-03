/**
 * @file solver_app.hpp
 * @brief 应用层 - 求解器应用模块头文件
 * @details 提供FETI-DP求解器的应用接口，整合各模块功能
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include <string>

namespace app {

/**
 * @class SolverApp
 * @brief 求解器应用类
 */
class SolverApp {
public:
    /**
     * @brief 构造函数
     */
    SolverApp();

    /**
     * @brief 析构函数
     */
    ~SolverApp();

    /**
     * @brief 初始化求解器应用
     * @param config_file 配置文件路径
     * @return bool 初始化成功返回true，失败返回false
     */
    bool initialize(const std::string& config_file);

    /**
     * @brief 运行求解器
     * @return bool 求解成功返回true，失败返回false
     */
    bool run();

    /**
     * @brief 清理资源
     */
    void cleanup();

private:
    bool initialized_;
};

} // namespace app