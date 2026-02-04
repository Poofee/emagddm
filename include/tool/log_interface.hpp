/**
 * @file log_interface.hpp
 * @brief 基础工具层 - 日志抽象接口模块头文件
 * @details 提供统一的日志接口抽象，支持多种日志库实现（spdlog、glog等）
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include <string>
#include <memory>

namespace tool {

/**
 * @enum LogLevel
 * @brief 日志级别枚举
 */
enum class LogLevel {
    TRACE,    ///< 跟踪级别，用于详细跟踪信息
    DEBUG,    ///< 调试级别，用于开发调试
    INFO,     ///< 信息级别，用于正常运行信息
    WARN,     ///< 警告级别，用于非致命错误
    ERR,      ///< 错误级别，用于致命错误（避免与Windows ERROR宏冲突）
    CRITICAL  ///< 严重级别，用于系统级严重错误
};

/**
 * @class ILogger
 * @brief 日志接口抽象类
 * @details 定义统一的日志操作接口，支持多种日志库实现
 */
class ILogger {
public:
    virtual ~ILogger() = default;

    /**
     * @brief 初始化日志系统
     * @param log_file 日志文件路径
     * @param console_output 是否输出到控制台
     * @return bool 初始化成功返回true，失败返回false
     */
    virtual bool initialize(const std::string& log_file, bool console_output = true) = 0;

    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    virtual void setLevel(LogLevel level) = 0;

    /**
     * @brief 输出跟踪日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void trace(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 输出调试日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void debug(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 输出信息日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void info(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 输出警告日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void warn(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 输出错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void err(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 输出严重错误日志
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    virtual void critical(const std::string& message, const std::string& module = "") = 0;

    /**
     * @brief 刷新日志缓冲区
     */
    virtual void flush() = 0;

    /**
     * @brief 检查是否已初始化
     * @return bool 已初始化返回true，否则返回false
     */
    virtual bool isInitialized() const = 0;
};

/**
 * @enum LoggerType
 * @brief 日志库类型枚举
 */
enum class LoggerType {
    SPDLOG,  ///< spdlog日志库
    GLOG,    ///< Google glog日志库
    CUSTOM   ///< 自定义日志库
};

/**
 * @brief 创建日志器实例
 * @param type 日志库类型
 * @return std::unique_ptr<ILogger> 日志器实例指针
 */
std::unique_ptr<ILogger> createLogger(LoggerType type = LoggerType::SPDLOG);

/**
 * @brief 获取日志级别字符串表示
 * @param level 日志级别
 * @return std::string 日志级别字符串
 */
std::string logLevelToString(LogLevel level);

/**
 * @brief 从字符串解析日志级别
 * @param level_str 日志级别字符串
 * @return LogLevel 日志级别
 */
LogLevel logLevelFromString(const std::string& level_str);

} // namespace tool