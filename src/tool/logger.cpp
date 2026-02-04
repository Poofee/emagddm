/**
 * @file logger.cpp
 * @brief 基础工具层 - 日志管理模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#include "logger.hpp"

namespace tool {

Logger::Logger() : logger_type_(LoggerType::SPDLOG), initialized_(false) {
}

Logger::~Logger() {
    if (logger_impl_) {
        logger_impl_->flush();
    }
}

Logger& Logger::getInstance() {
    static Logger instance;
    return instance;
}

bool Logger::initialize(const std::string& log_file, 
                       bool console_output,
                       LoggerType logger_type) {
    if (initialized_) {
        // 已经初始化，先清理
        logger_impl_.reset();
        initialized_ = false;
    }
    
    try {
        logger_impl_ = createLogger(logger_type);
        logger_type_ = logger_type;
        
        if (logger_impl_->initialize(log_file, console_output)) {
            initialized_ = true;
            return true;
        }
    }
    catch (const std::exception& e) {
        // 初始化失败，回退到默认实现
        logger_impl_.reset();
    }
    
    return false;
}

void Logger::setLevel(LogLevel level) {
    if (logger_impl_ && initialized_) {
        logger_impl_->setLevel(level);
    }
}

void Logger::trace(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->trace(message, module);
    }
}

void Logger::debug(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->debug(message, module);
    }
}

void Logger::info(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->info(message, module);
    }
}

void Logger::warn(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->warn(message, module);
    }
}

void Logger::err(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->err(message, module);
    }
}

void Logger::critical(const std::string& message, const std::string& module) {
    if (logger_impl_ && initialized_) {
        logger_impl_->critical(message, module);
    }
}

void Logger::flush() {
    if (logger_impl_ && initialized_) {
        logger_impl_->flush();
    }
}

bool Logger::isInitialized() const {
    return initialized_ && logger_impl_ && logger_impl_->isInitialized();
}

LoggerType Logger::getLoggerType() const {
    return logger_type_;
}

} // namespace tool