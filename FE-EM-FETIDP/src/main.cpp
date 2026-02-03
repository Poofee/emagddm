/**
 * @file main.cpp
 * @brief 驱动层主函数
 * @details FETI-DP求解器的主入口程序，负责命令行参数解析、应用初始化和运行
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include <iostream>
#include <string>
#include "solver_app.hpp"

/**
 * @brief 显示使用帮助
 */
void showUsage() {
    std::cout << "FETI-DP电磁有限元求解器" << std::endl;
    std::cout << "用法: fetidp_solver <配置文件路径>" << std::endl;
    std::cout << "示例: fetidp_solver config/motor2d_steady.json" << std::endl;
}

/**
 * @brief 主函数
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return int 程序退出码（0表示成功，非0表示失败）
 */
int main(int argc, char* argv[]) {
    // 检查命令行参数
    if (argc != 2) {
        showUsage();
        return 1;
    }

    std::string config_file = argv[1];
    
    // 创建求解器应用实例
    app::SolverApp solver;
    
    // 初始化求解器
    if (!solver.initialize(config_file)) {
        std::cerr << "求解器初始化失败" << std::endl;
        return 1;
    }
    
    // 运行求解器
    if (!solver.run()) {
        std::cerr << "求解器运行失败" << std::endl;
        return 1;
    }
    
    std::cout << "求解器运行成功完成" << std::endl;
    
    return 0;
}