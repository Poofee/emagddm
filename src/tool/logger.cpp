/**
 * @file logger.cpp
 * @brief 日志模块实现
 * @details 实现Logger类的所有方法，包括单例管理、日志输出、格式化等
 *          使用PIMPL模式封装具体日志库实现（当前为spdlog）
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#include "logger.hpp"
#include "spdlog_adapter.hpp"

namespace tool {

// 全局唯一实例定义（单例模式）
Logger Logger::instance_;

/**
 * @brief 构造函数，初始化成员变量默认值
 * @param type 日志库类型，默认使用SPDLOG
 */
Logger::Logger(LoggerType type) : logger_type_(type), initialized_(false) {
    (void)type;  // 预留扩展接口，当前固定使用SPDLOG
}

/**
 * @brief 析构函数，安全释放日志资源
 * 
 * 使用try-catch确保析构过程不会抛出异常，
 * 避免程序退出时因日志系统销毁顺序问题导致崩溃
 */
Logger::~Logger() {
    try {
        if (logger_impl_ && initialized_) {
            logger_impl_.reset();
        }
        initialized_ = false;
    } catch (...) {
    }
}

/**
 * @brief 获取全局单例实例引用
 * @return Logger& 唯一的Logger对象引用
 */
Logger& Logger::getInstance() {
    return instance_;
}

/**
 * @brief 初始化日志系统
 * 
 * 创建具体的日志库实现对象并配置输出目标。
 * 支持重复调用检查，已初始化则直接返回成功。
 * 
 * @param log_file 日志文件路径（空字符串表示仅控制台输出）
 * @param console_output 是否同时输出到控制台
 * @return bool 初始化成功返回true，失败返回false
 */
bool Logger::initialize(const std::string& log_file, bool console_output) {
    // 防止重复初始化
    if (initialized_) {
        return true;
    }

    try {
        // 创建日志库适配器实例
        logger_impl_ = createLogger(logger_type_);
        if (!logger_impl_) {
            return false;
        }

        // 调用底层接口完成初始化配置
        if (!logger_impl_->initialize(log_file, console_output)) {
            logger_impl_.reset();  // 初始化失败时释放资源
            return false;
        }

        initialized_ = true;
        return true;
    } 
    // 捕获标准异常（如文件创建失败、权限不足等）
    catch (const std::exception& e) {
        (void)e;  // 静默处理错误信息，避免日志系统自身出错导致递归
        logger_impl_.reset();
        initialized_ = false;
        return false;
    } 
    // 捕获未知异常
    catch (...) {
        logger_impl_.reset();
        initialized_ = false;
        return false;
    }
}

/**
 * @brief 设置日志过滤级别
 * 
 * 只有大于等于该级别的日志才会被输出，
 * 用于控制日志输出的详细程度
 * 
 * @param level 目标日志级别
 */
void Logger::setLevel(LogLevel level) {
    if (logger_impl_) {
        logger_impl_->setLevel(level);
    }
}

void Logger::trace(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->trace(message, module);
    }
}

void Logger::debug(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->debug(message, module);
    }
}

void Logger::info(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->info(message, module);
    }
}

void Logger::warn(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->warn(message, module);
    }
}

void Logger::err(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->err(message, module);
    }
}

void Logger::critical(const std::string& message, const std::string& module) {
    if (logger_impl_) {
        logger_impl_->critical(message, module);
    }
}

/**
 * @brief 强制刷新日志缓冲区
 * 
 * 将缓冲区中的待写入日志立即刷新到输出目标（文件/控制台），
 * 通常在程序退出或关键操作点调用以确保日志完整性
 */
void Logger::flush() {
    if (logger_impl_) {
        logger_impl_->flush();
    }
}

/**
 * @brief 检查日志器是否已完成初始化
 * @return bool 已初始化返回true
 */
bool Logger::isInitialized() const {
    return initialized_;
}

/**
 * @brief 获取当前使用的日志库类型
 * @return LoggerType 日志库类型枚举值
 */
LoggerType Logger::getLoggerType() const {
    return logger_type_;
}

} // namespace tool
