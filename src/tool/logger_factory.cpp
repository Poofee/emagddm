/**
 * @file logger_factory.cpp
 * @brief 日志工厂类实现
 * @details 实现LoggerFactory的静态方法，提供全局默认日志器的访问接口
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#include "logger_factory.hpp"

namespace tool {

/**
 * @brief 获取全局默认日志器实例
 * 
 * 通过Logger单例模式获取唯一的日志器对象，
 * 所有模块共享同一个日志器，确保日志输出统一管理
 * 
 * @return Logger& 全局日志器引用
 */
Logger& LoggerFactory::getDefaultLogger() {
    return Logger::getInstance();
}

/**
 * @brief 初始化默认日志器
 * 
 * 创建并配置全局日志器，设置输出目标（文件/控制台）。
 * 如果日志器已初始化则直接返回成功。
 * 
 * @param log_file 日志文件路径（空字符串表示不写文件）
 * @param console_output 是否同时输出到控制台
 * @param logger_type 日志库类型（当前保留参数，仅支持SPDLOG）
 * @return bool 初始化成功返回true
 */
bool LoggerFactory::initializeDefaultLogger(const std::string& log_file,
                                           bool console_output,
                                           LoggerType logger_type) {
    (void)logger_type;  // 当前版本固定使用SPDLOG，预留扩展接口
    return getDefaultLogger().initialize(log_file, console_output);
}

/**
 * @brief 设置默认日志器的日志级别
 * @param level 目标日志级别（TRACE/DEBUG/INFO/WARN/ERROR/CRITICAL）
 */
void LoggerFactory::setDefaultLoggerLevel(LogLevel level) {
    getDefaultLogger().setLevel(level);
}

/**
 * @brief 检查默认日志器是否已初始化
 * @return bool 已初始化返回true
 */
bool LoggerFactory::isDefaultLoggerInitialized() {
    return getDefaultLogger().isInitialized();
}

} // namespace tool
