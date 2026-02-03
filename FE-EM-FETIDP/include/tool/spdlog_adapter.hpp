/**
 * @file spdlog_adapter.hpp
 * @brief 基础工具层 - spdlog适配器模块头文件
 * @details 实现spdlog日志库的适配器，封装spdlog接口为统一的ILogger接口
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include "log_interface.hpp"

namespace tool {

/**
 * @class SpdlogAdapter
 * @brief spdlog适配器类
 * @details 将spdlog库的功能适配到统一的ILogger接口
 */
class SpdlogAdapter : public ILogger {
public:
    /**
     * @brief 构造函数
     */
    SpdlogAdapter();

    /**
     * @brief 析构函数
     */
    ~SpdlogAdapter() override;

    /**
     * @brief 初始化日志系统
     * @param log_file 日志文件路径
     * @param console_output 是否输出到控制台
     * @return bool 初始化成功返回true，失败返回false
     */
    bool initialize(const std::string& log_file, bool console_output = true) override;

    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    void setLevel(LogLevel level) override;

    /**
     * @brief 输出跟踪日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void trace(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 输出调试日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void debug(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 输出信息日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void info(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 输出警告日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void warn(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 输出错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void err(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 输出严重错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void critical(const std::string& message, const std::string& module = "") override;

    /**
     * @brief 刷新日志缓冲区
     */
    void flush() override;

    /**
     * @brief 检查是否已初始化
     * @return bool 已初始化返回true，否则返回false
     */
    bool isInitialized() const override;

private:
    class Impl;  ///< 实现类前向声明
    std::unique_ptr<Impl> pimpl_;  ///< PIMPL模式实现
};

} // namespace tool