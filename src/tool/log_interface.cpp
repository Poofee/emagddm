/**
 * @file log_interface.cpp
 * @brief 基础工具层 - 日志抽象接口模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "log_interface.hpp"
#include "spdlog_adapter.hpp"
#include <stdexcept>

namespace tool {

std::unique_ptr<ILogger> createLogger(LoggerType type) {
    switch (type) {
        case LoggerType::SPDLOG:
            return std::make_unique<SpdlogAdapter>();
        case LoggerType::GLOG:
            // TODO: 实现glog适配器
            throw std::runtime_error("GLOG logger not implemented yet");
        case LoggerType::CUSTOM:
            // TODO: 实现自定义日志器
            throw std::runtime_error("Custom logger not implemented yet");
        default:
            throw std::runtime_error("Unknown logger type");
    }
}

std::string logLevelToString(LogLevel level) {
    switch (level) {
        case LogLevel::TRACE:
            return "TRACE";
        case LogLevel::DEBUG:
            return "DEBUG";
        case LogLevel::INFO:
            return "INFO";
        case LogLevel::WARN:
            return "WARN";
        case LogLevel::ERR:
            return "ERROR";
        case LogLevel::CRITICAL:
            return "CRITICAL";
        default:
            return "UNKNOWN";
    }
}

LogLevel logLevelFromString(const std::string& level_str) {
    if (level_str == "TRACE" || level_str == "trace") {
        return LogLevel::TRACE;
    } else if (level_str == "DEBUG" || level_str == "debug") {
        return LogLevel::DEBUG;
    } else if (level_str == "INFO" || level_str == "info") {
        return LogLevel::INFO;
    } else if (level_str == "WARN" || level_str == "warn") {
        return LogLevel::WARN;
    } else if (level_str == "ERROR" || level_str == "error") {
        return LogLevel::ERR;
    } else if (level_str == "CRITICAL" || level_str == "critical") {
        return LogLevel::CRITICAL;
    } else {
        throw std::runtime_error("Unknown log level: " + level_str);
    }
}

} // namespace tool