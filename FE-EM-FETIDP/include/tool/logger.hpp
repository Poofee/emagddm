/**
 * @file logger.hpp
 * @brief 基础工具层 - 日志管理模块头文件
 * @details 提供统一的日志接口，基于抽象接口支持多种日志库实现
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#pragma once

#include "log_interface.hpp"
#include <string>
#include <memory>

namespace tool {

/**
 * @class Logger
 * @brief 日志管理类（单例模式）
 * @details 基于抽象日志接口，支持多种日志库实现，提供统一的日志操作接口
 */
class Logger {
public:
    /**
     * @brief 获取日志实例（单例模式）
     * @return Logger& 日志实例引用
     */
    static Logger& getInstance();

    /**
     * @brief 初始化日志系统
     * @param log_file 日志文件路径
     * @param console_output 是否输出到控制台
     * @param logger_type 日志库类型（默认spdlog）
     * @return bool 初始化成功返回true，失败返回false
     */
    bool initialize(const std::string& log_file, 
                   bool console_output = true,
                   LoggerType logger_type = LoggerType::SPDLOG);

    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    void setLevel(LogLevel level);

    /**
     * @brief 输出跟踪日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void trace(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出调试日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void debug(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出信息日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void info(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出警告日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void warn(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void err(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出严重错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void critical(const std::string& message, const std::string& module = "");

    /**
     * @brief 刷新日志缓冲区
     */
    void flush();

    /**
     * @brief 检查是否已初始化
     * @return bool 已初始化返回true，否则返回false
     */
    bool isInitialized() const;

    /**
     * @brief 获取当前使用的日志库类型
     * @return LoggerType 日志库类型
     */
    LoggerType getLoggerType() const;

private:
    Logger();
    ~Logger();

    std::unique_ptr<ILogger> logger_impl_;
    LoggerType logger_type_;
    bool initialized_;
};

} // namespace tool