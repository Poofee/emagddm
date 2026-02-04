/**
 * @file logger_factory.hpp
 * @brief 基础工具层 - 日志工厂类
 * @details 提供默认日志器获取和日志宏定义
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include "logger.hpp"

namespace tool {

/**
 * @class LoggerFactory
 * @brief 日志工厂类，提供默认日志器获取
 */
class LoggerFactory {
public:
    /**
     * @brief 获取默认日志器
     * @return Logger& 默认日志器引用
     */
    static Logger& getDefaultLogger() {
        return Logger::getInstance();
    }
    
    /**
     * @brief 初始化默认日志器
     * @param log_file 日志文件路径
     * @param console_output 是否输出到控制台
     * @param logger_type 日志库类型
     * @return bool 初始化是否成功
     */
    static bool initializeDefaultLogger(const std::string& log_file = "",
                                       bool console_output = true,
                                       LoggerType logger_type = LoggerType::SPDLOG) {
        return getDefaultLogger().initialize(log_file, console_output, logger_type);
    }
    
    /**
     * @brief 设置默认日志器级别
     * @param level 日志级别
     */
    static void setDefaultLoggerLevel(LogLevel level) {
        getDefaultLogger().setLevel(level);
    }
    
    /**
     * @brief 检查默认日志器是否已初始化
     * @return bool 是否已初始化
     */
    static bool isDefaultLoggerInitialized() {
        return getDefaultLogger().isInitialized();
    }
};

} // namespace tool

/**
 * @def FEEM_TRACE(...)
 * @brief 输出跟踪级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_TRACE(...) tool::LoggerFactory::getDefaultLogger().trace(__VA_ARGS__)

/**
 * @def FEEM_DEBUG(...)
 * @brief 输出调试级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_DEBUG(...) tool::LoggerFactory::getDefaultLogger().debug(__VA_ARGS__)

/**
 * @def FEEM_INFO(...)
 * @brief 输出信息级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_INFO(...) tool::LoggerFactory::getDefaultLogger().info(__VA_ARGS__)

/**
 * @def FEEM_WARN(...)
 * @brief 输出警告级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_WARN(...) tool::LoggerFactory::getDefaultLogger().warn(__VA_ARGS__)

/**
 * @def FEEM_ERROR(...)
 * @brief 输出错误级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_ERROR(...) tool::LoggerFactory::getDefaultLogger().err(__VA_ARGS__)

/**
 * @def FEEM_CRITICAL(...)
 * @brief 输出严重错误级别日志的宏
 * @param ... 日志消息和模块名称（可选）
 */
#define FEEM_CRITICAL(...) tool::LoggerFactory::getDefaultLogger().critical(__VA_ARGS__)