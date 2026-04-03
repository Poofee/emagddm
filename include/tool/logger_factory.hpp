/**
 * @file logger_factory.hpp
 * @brief 基础工具层 - 日志工厂类
 * @details 提供默认日志器获取和日志宏定义
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
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
    static Logger& getDefaultLogger();
    
    static bool initializeDefaultLogger(const std::string& log_file = "",
                                       bool console_output = true,
                                       LoggerType logger_type = LoggerType::SPDLOG);
    
    static void setDefaultLoggerLevel(LogLevel level);
    
    static bool isDefaultLoggerInitialized();
};

} // namespace tool
