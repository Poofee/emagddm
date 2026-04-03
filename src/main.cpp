/**
 * @file main.cpp
 * @brief 驱动层主函数
 * @details FETI-DP求解器的主入口程序
 *          main函数只负责：命令行参数解析、创建SolverApp、调用run()
 *          所有业务逻辑都在SolverApp内部处理
 * @author Poofee
 * @date 2026-XX-XX
 * @version 2.0
 */

#include <iostream>
#include "solver_app.hpp"
#include "logger_factory.hpp"

/**
 * @brief 显示程序使用帮助信息
 */
void showUsage() {
    std::cout << "FETI-DP电磁有限元求解器" << std::endl;
    std::cout << "用法: fetidp_solver <输入文件>" << std::endl;
    std::cout << std::endl;
    std::cout << "支持的输入文件:" << std::endl;
    std::cout << "  *.json   - JSON配置文件（启动求解模式）" << std::endl;
    std::cout << "  *.aedt   - ANSYS Maxwell项目文件（启动解析模式）" << std::endl;
    std::cout << std::endl;
    std::cout << "示例:" << std::endl;
    std::cout << "  fetidp_solver config/motor2d_steady.json" << std::endl;
    std::cout << "  fetidp_solver docs/project/Project49.aedt" << std::endl;
}

/**
 * @brief 程序主入口函数
 * 
 * 职责：
 * 1. 初始化日志系统
 * 2. 解析命令行参数
 * 3. 创建SolverApp实例
 * 4. 调用initialize()和run()
 * 5. 返回退出码
 * 
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return 程序退出码（0成功，非0失败）
 */
int main(int argc, char* argv[]) {
    // 初始化日志系统
    tool::LoggerFactory::initializeDefaultLogger("", true);
    
    int exit_code = 0;
    
    try {
        // 至少需要一个命令行参数（输入文件路径）
        if (argc < 2) {
            showUsage();
            return 1;
        }
        
        // 处理--help选项
        std::string arg1 = argv[1];
        if (arg1 == "--help" || arg1 == "-h") {
            showUsage();
            return 0;
        }
        
        // 创建求解器应用实例
        app::SolverApp solver;
        
        // 初始化（自动根据文件扩展名判断运行模式）
        if (!solver.initialize(arg1)) {
            FEEM_ERROR("求解器初始化失败");
            return 1;
        }
        
        // 运行求解器（根据模式执行对应逻辑）
        if (!solver.run()) {
            FEEM_ERROR("求解器运行失败");
            exit_code = 1;
        } else {
            FEEM_INFO("求解器运行成功完成");
            exit_code = 0;
        }
        
    } catch (const std::exception& e) {
        FEEM_ERROR("程序异常: {}", e.what());
        exit_code = 1;
    } catch (...) {
        FEEM_ERROR("未知异常");
        exit_code = 1;
    }
    
    // 刷新日志缓冲区
    tool::LoggerFactory::getDefaultLogger().flush();
    
    return exit_code;
}
