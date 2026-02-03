/**
 * @file test_logger.cpp
 * @brief 日志功能测试程序
 * @details 测试新设计的日志抽象接口和spdlog适配器
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#include "tool/logger_factory.hpp"
#include <iostream>
#include <thread>
#include <chrono>

/**
 * @brief 测试基本日志功能
 */
void testBasicLogging() {
    std::cout << "=== 测试基本日志功能 ===" << std::endl;
    
    // 初始化日志系统
    if (!tool::LoggerFactory::initializeDefaultLogger("output/log/test.log", true)) {
        std::cerr << "日志系统初始化失败" << std::endl;
        return;
    }
    
    std::cout << "日志系统初始化成功，使用日志库: " 
              << static_cast<int>(tool::LoggerFactory::getDefaultLogger().getLoggerType()) << std::endl;
    
    // 测试不同级别的日志
    FEEM_DEBUG("这是一条调试日志", "test");
    FEEM_INFO("这是一条信息日志", "test");
    FEEM_WARN("这是一条警告日志", "test");
    FEEM_ERROR("这是一条错误日志", "test");
    
    // 测试不带模块名的日志
    FEEM_INFO("这是不带模块名的信息日志");
    
    std::cout << "基本日志功能测试完成" << std::endl;
}

/**
 * @brief 测试日志级别设置
 */
void testLogLevel() {
    std::cout << "\n=== 测试日志级别设置 ===" << std::endl;
    
    // 设置为DEBUG级别
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::DEBUG);
    FEEM_DEBUG("DEBUG级别下可以看到这条日志", "level_test");
    
    // 设置为INFO级别
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::INFO);
    FEEM_DEBUG("INFO级别下不应该看到这条DEBUG日志", "level_test");
    FEEM_INFO("INFO级别下可以看到这条日志", "level_test");
    
    // 设置为WARN级别
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::WARN);
    FEEM_INFO("WARN级别下不应该看到这条INFO日志", "level_test");
    FEEM_WARN("WARN级别下可以看到这条日志", "level_test");
    
    // 设置为ERR级别
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::ERR);
    FEEM_WARN("ERR级别下不应该看到这条WARN日志", "level_test");
    FEEM_ERROR("ERR级别下可以看到这条日志", "level_test");
    
    // 测试TRACE和CRITICAL级别
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::TRACE);
    FEEM_TRACE("TRACE级别下可以看到这条日志", "level_test");
    
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::CRITICAL);
    FEEM_ERROR("CRITICAL级别下不应该看到这条ERR日志", "level_test");
    FEEM_CRITICAL("CRITICAL级别下可以看到这条日志", "level_test");
    
    std::cout << "日志级别设置测试完成" << std::endl;
}

/**
 * @brief 测试多线程日志
 */
void testMultiThreadLogging() {
    std::cout << "\n=== 测试多线程日志 ===" << std::endl;
    
    // 创建多个线程同时写日志
    auto logTask = [](int thread_id) {
        for (int i = 0; i < 5; ++i) {
            FEEM_INFO("线程 " + std::to_string(thread_id) + " 的第 " + 
                      std::to_string(i) + " 条日志", "thread_test");
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    };
    
    std::thread t1(logTask, 1);
    std::thread t2(logTask, 2);
    std::thread t3(logTask, 3);
    
    t1.join();
    t2.join();
    t3.join();
    
    std::cout << "多线程日志测试完成" << std::endl;
}

/**
 * @brief 测试日志刷新
 */
void testLogFlush() {
    std::cout << "\n=== 测试日志刷新 ===" << std::endl;
    
    FEEM_INFO("刷新前的日志", "flush_test");
    tool::LoggerFactory::getDefaultLogger().flush();
    FEEM_INFO("刷新后的日志", "flush_test");
    
    std::cout << "日志刷新测试完成" << std::endl;
}

/**
 * @brief 主函数
 */
int main() {
    std::cout << "开始测试日志功能..." << std::endl;
    
    try {
        testBasicLogging();
        testLogLevel();
        testMultiThreadLogging();
        testLogFlush();
        
        std::cout << "\n=== 所有测试完成 ===" << std::endl;
        std::cout << "请检查 output/log/test.log 文件查看日志输出" << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "测试过程中发生异常: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "测试过程中发生未知异常" << std::endl;
        return 1;
    }
}