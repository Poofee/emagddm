/**
 * @file spdlog_adapter.cpp
 * @brief 基础工具层 - spdlog适配器模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "spdlog_adapter.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <memory>
#include <vector>

namespace tool {

// PIMPL实现类
class SpdlogAdapter::Impl {
public:
    Impl() : logger_(nullptr), initialized_(false) {}
    
    ~Impl() {
        if (logger_) {
            spdlog::drop(logger_->name());
        }
    }
    
    bool initialize(const std::string& log_file, bool console_output) {
        try {
            std::vector<spdlog::sink_ptr> sinks;
            
            // 控制台输出
            if (console_output) {
                auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
                console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
                sinks.push_back(console_sink);
            }
            
            // 文件输出
            if (!log_file.empty()) {
                auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
                file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] [%n] %v");
                sinks.push_back(file_sink);
            }
            
            // 创建日志器
            logger_ = std::make_shared<spdlog::logger>("fetidp", begin(sinks), end(sinks));
            spdlog::register_logger(logger_);
            
            // 设置默认级别
            logger_->set_level(spdlog::level::info);
            logger_->flush_on(spdlog::level::warn);
            
            initialized_ = true;
            return true;
        }
        catch (const spdlog::spdlog_ex&) {
            // 如果初始化失败，尝试创建简单的控制台日志器
            try {
                logger_ = spdlog::stdout_color_mt("fetidp");
                logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
                logger_->set_level(spdlog::level::info);
                initialized_ = true;
                return true;
            }
            catch (...) {
                return false;
            }
        }
    }
    
    void setLevel(LogLevel level) {
        if (!logger_) return;
        
        switch (level) {
            case LogLevel::TRACE:
                logger_->set_level(spdlog::level::trace);
                break;
            case LogLevel::DEBUG:
                logger_->set_level(spdlog::level::debug);
                break;
            case LogLevel::INFO:
                logger_->set_level(spdlog::level::info);
                break;
            case LogLevel::WARN:
                logger_->set_level(spdlog::level::warn);
                break;
            case LogLevel::ERR:
                logger_->set_level(spdlog::level::err);
                break;
            case LogLevel::CRITICAL:
                logger_->set_level(spdlog::level::critical);
                break;
            default:
                logger_->set_level(spdlog::level::info);
                break;
        }
    }
    
    void trace(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->trace(message);
            } else {
                logger_->trace("[{}] {}", module, message);
            }
        }
    }
    
    void debug(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->debug(message);
            } else {
                logger_->debug("[{}] {}", module, message);
            }
        }
    }
    
    void info(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->info(message);
            } else {
                logger_->info("[{}] {}", module, message);
            }
        }
    }
    
    void warn(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->warn(message);
            } else {
                logger_->warn("[{}] {}", module, message);
            }
        }
    }
    
    void err(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->error(message);
            } else {
                logger_->error("[" + module + "] " + message);
            }
        }
    }
    
    void critical(const std::string& message, const std::string& module) {
        if (logger_ && initialized_) {
            if (module.empty()) {
                logger_->critical(message);
            } else {
                logger_->critical("[{}] {}", module, message);
            }
        }
    }
    
    void flush() {
        if (logger_ && initialized_) {
            logger_->flush();
        }
    }
    
    bool isInitialized() const {
        return initialized_ && (logger_ != nullptr);
    }
    
private:
    std::shared_ptr<spdlog::logger> logger_;
    bool initialized_;
};

// SpdlogAdapter实现
SpdlogAdapter::SpdlogAdapter() : pimpl_(std::make_unique<Impl>()) {}

SpdlogAdapter::~SpdlogAdapter() = default;

bool SpdlogAdapter::initialize(const std::string& log_file, bool console_output) {
    return pimpl_->initialize(log_file, console_output);
}

void SpdlogAdapter::setLevel(LogLevel level) {
    pimpl_->setLevel(level);
}

void SpdlogAdapter::trace(const std::string& message, const std::string& module) {
    pimpl_->trace(message, module);
}

void SpdlogAdapter::debug(const std::string& message, const std::string& module) {
    pimpl_->debug(message, module);
}

void SpdlogAdapter::info(const std::string& message, const std::string& module) {
    pimpl_->info(message, module);
}

void SpdlogAdapter::warn(const std::string& message, const std::string& module) {
    pimpl_->warn(message, module);
}

void SpdlogAdapter::critical(const std::string& message, const std::string& module) {
    pimpl_->critical(message, module);
}

void SpdlogAdapter::err(const std::string& message, const std::string& module) {
    pimpl_->err(message, module);
}

void SpdlogAdapter::flush() {
    pimpl_->flush();
}

bool SpdlogAdapter::isInitialized() const {
    return pimpl_->isInitialized();
}

} // namespace tool