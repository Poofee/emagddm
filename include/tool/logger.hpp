/**
 * @file logger.hpp
 * @brief 基础工具层 - 日志模块头文件
 * @details 提供统一的日志接口，支持多种日志库实现（spdlog、glog等）
 *          支持通过宏简化调用方式，支持格式化参数
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#pragma once

#include "log_interface.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

namespace tool {

/**
 * @class Logger
 * @brief 日志器类
 * @details 封装日志操作，提供统一的日志输出接口，支持格式化参数
 *          使用单例模式管理全局唯一日志器实例
 */
class Logger {
public:
    /**
     * @brief 构造函数
     * @param type 日志库类型，默认使用SPDLOG
     */
    explicit Logger(LoggerType type = LoggerType::SPDLOG);
    
    /**
     * @brief 析构函数，负责释放日志资源
     */
    ~Logger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    /**
     * @brief 获取全局唯一日志器实例（单例模式）
     * @return Logger& 全局日志器引用
     */
    static Logger& getInstance();

    /**
     * @brief 初始化日志系统
     * @param log_file 日志文件路径，为空则仅输出到控制台
     * @param console_output 是否同时输出到控制台
     * @return bool 初始化成功返回true
     */
    bool initialize(const std::string& log_file = "", bool console_output = true);
    
    /**
     * @brief 设置日志级别
     * @param level 目标日志级别
     */
    void setLevel(LogLevel level);
    
    /**
     * @brief 输出跟踪级别日志（无格式化）
     * @param message 日志消息
     * @param module 模块名称（可选）
     */
    void trace(const std::string& message, const std::string& module = "");
    
    /** 
     * @brief 输出调试级别日志（无格式化）
     */
    void debug(const std::string& message, const std::string& module = "");
    
    /** 
     * @brief 输出信息级别日志（无格式化）
     */
    void info(const std::string& message, const std::string& module = "");
    
    /** 
     * @brief 输出警告级别日志（无格式化）
     */
    void warn(const std::string& message, const std::string& module = "");
    
    /** 
     * @brief 输出错误级别日志（无格式化）
     */
    void err(const std::string& message, const std::string& module = "");
    
    /** 
     * @brief 输出严重错误级别日志（无格式化）
     */
    void critical(const std::string& message, const std::string& module = "");

    /**
     * @brief 输出跟踪级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串，使用{}作为占位符
     * @param args 要填充的参数列表
     * 
     * @code
     * traceFmt("处理文件: {}, 大小: {}", filename, size);
     * @endcode
     */
    template<typename... Args>
    void traceFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->trace(msg);
        }
    }

    /**
     * @brief 输出调试级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串
     * @param args 参数列表
     */
    template<typename... Args>
    void debugFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->debug(msg);
        }
    }

    /**
     * @brief 输出信息级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串
     * @param args 参数列表
     */
    template<typename... Args>
    void infoFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->info(msg);
        }
    }

    /**
     * @brief 输出警告级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串
     * @param args 参数列表
     */
    template<typename... Args>
    void warnFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->warn(msg);
        }
    }

    /**
     * @brief 输出错误级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串
     * @param args 参数列表
     */
    template<typename... Args>
    void errFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->err(msg);
        }
    }

    /**
     * @brief 输出严重错误级别日志（支持格式化参数）
     * @tparam Args 可变参数类型
     * @param fmt_str 格式字符串
     * @param args 参数列表
     */
    template<typename... Args>
    void criticalFmt(const std::string& fmt_str, Args&&... args) {
        if (logger_impl_ && initialized_) {
            auto msg = spdlog::fmt_lib::format(fmt_str, std::forward<Args>(args)...);
            logger_impl_->critical(msg);
        }
    }

    /**
     * @brief 刷新日志缓冲区，确保所有日志写入目标
     */
    void flush();
    
    /**
     * @brief 检查日志器是否已初始化
     * @return bool 已初始化返回true
     */
    bool isInitialized() const;
    
    /**
     * @brief 获取当前使用的日志库类型
     * @return LoggerType 日志库类型枚举值
     */
    LoggerType getLoggerType() const;

private:
    static Logger instance_;  ///< 全局唯一实例（单例模式）
    std::unique_ptr<ILogger> logger_impl_;  ///< 日志实现接口指针
    LoggerType logger_type_;  ///< 当前使用的日志库类型
    bool initialized_;  ///< 是否已初始化标志
};

} // namespace tool

/**
 * @def FEEM_TRACE(...)
 * @brief 输出跟踪级别日志的宏（支持格式化参数）
 * @param ... 格式字符串和参数列表，使用{}占位符
 * 
 * @code
 * FEEM_TRACE("变量x的值为: {}", x);
 * @endcode
 */
#define FEEM_TRACE(...)   tool::LoggerFactory::getDefaultLogger().traceFmt(__VA_ARGS__)

/**
 * @def FEEM_DEBUG(...)
 * @brief 输出调试级别日志的宏（支持格式化参数）
 */
#define FEEM_DEBUG(...)   tool::LoggerFactory::getDefaultLogger().debugFmt(__VA_ARGS__)

/**
 * @def FEEM_INFO(...)
 * @brief 输出信息级别日志的宏（支持格式化参数）
 */
#define FEEM_INFO(...)    tool::LoggerFactory::getDefaultLogger().infoFmt(__VA_ARGS__)

/**
 * @def FEEM_WARN(...)
 * @brief 输出警告级别日志的宏（支持格式化参数）
 */
#define FEEM_WARN(...)    tool::LoggerFactory::getDefaultLogger().warnFmt(__VA_ARGS__)

/**
 * @def FEEM_ERROR(...)
 * @brief 输出错误级别日志的宏（支持格式化参数）
 */
#define FEEM_ERROR(...)   tool::LoggerFactory::getDefaultLogger().errFmt(__VA_ARGS__)

/**
 * @def FEEM_CRITICAL(...)
 * @brief 输出严重错误级别日志的宏（支持格式化参数）
 */
#define FEEM_CRITICAL(...) tool::LoggerFactory::getDefaultLogger().criticalFmt(__VA_ARGS__)

/**
 * @def FEEM_LOG_INIT(file_path, console_output)
 * @brief 初始化默认日志系统的便捷宏
 * @param file_path 日志文件路径（空字符串表示不写文件）
 * @param console_output 是否同时输出到控制台
 */
#define FEEM_LOG_INIT(file_path, console_output) \
    tool::LoggerFactory::initializeDefaultLogger(file_path, console_output)

/**
 * @def FEEM_SET_LEVEL(level)
 * @brief 设置默认日志级别的便捷宏
 * @param level 日志级别（LogLevel枚举值）
 */
#define FEEM_SET_LEVEL(level) \
    tool::LoggerFactory::getDefaultLogger().setLevel(level)
