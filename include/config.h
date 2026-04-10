// 项目配置头文件
// 定义编译期开关宏，用于控制MPI、OpenMP功能和跨平台适配

#ifndef FE_EM_CONFIG_H
#define FE_EM_CONFIG_H

// ============================================================
// 平台检测宏（由CMake根据目标操作系统自动定义）
// 使用方法：在C++源码中通过 #ifdef 进行条件编译
// ============================================================

// 操作系统类型（互斥，同一时间只有一个为1）
// - PLATFORM_WINDOWS: Windows系统（使用Win32/Win64 API）
// - PLATFORM_LINUX:   Linux系统（使用POSIX接口）
// - PLATFORM_MACOS:   macOS系统（POSIX + Darwin扩展）

#ifndef PLATFORM_WINDOWS
#define PLATFORM_WINDOWS 0
#endif

#ifndef PLATFORM_LINUX
#define PLATFORM_LINUX 0
#endif

#ifndef PLATFORM_MACOS
#define PLATFORM_MACOS 0
#endif

// 平台特定特性检测辅助宏（可选，按需启用）
#if PLATFORM_WINDOWS
    // Windows特有：使用Windows API头文件
    // 示例：#include <windows.h>
    #define USE_WIN32_API 1
    #define USE_POSIX_API 0
#elif PLATFORM_LINUX || PLATFORM_MACOS
    // Unix系系统：使用POSIX标准接口
    #define USE_WIN32_API 0
    #define USE_POSIX_API 1

    #if PLATFORM_MACOS
        // macOS特有：Darwin内核扩展
        #define USE_DARWIN_EXTENSIONS 1
    #else
        #define USE_DARWIN_EXTENSIONS 0
    #endif
#else
    // 未知平台：默认假设为类Unix系统（保守策略）
    #define USE_WIN32_API 0
    #define USE_POSIX_API 1
    #define USE_DARWIN_EXTENSIONS 0
#endif

// ============================================================
// MPI编译开关
// 0=关闭MPI，1=开启MPI
#ifndef USE_MPI
#define USE_MPI 0
#endif

// OpenMP编译开关
// 0=关闭OpenMP，1=开启OpenMP
#ifndef USE_OPENMP
#define USE_OPENMP 0
#endif

// MPI根目录（由CMake自动检测并定义）
#ifndef MPI_ROOT
#define MPI_ROOT ""
#endif

// ============================================================
// 编译模式自动推导
// 根据MPI和OpenMP的组合确定运行模式
// ============================================================

#if USE_MPI && USE_OPENMP
#define HYBRID_MODE_MPI_OMP 1
#define SERIAL_MODE 0
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 0
#elif USE_MPI
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 0
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 1
#elif USE_OPENMP
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 0
#define PURE_OMP_MODE 1
#define PURE_MPI_MODE 0
#else
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 1
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 0
#endif

// 模式描述字符串（用于日志输出和调试信息）
#if HYBRID_MODE_MPI_OMP
#define COMPILATION_MODE "MPI+OpenMP混合模式"
#elif PURE_MPI_MODE
#define COMPILATION_MODE "纯MPI模式"
#elif PURE_OMP_MODE
#define COMPILATION_MODE "纯OpenMP模式"
#else
#define COMPILATION_MODE "串行模式"
#endif

// ============================================================
// 跨平台工具函数示例（供参考）
// 在实际源文件中使用时，请参考以下模式：
//
// #include "config.h"
//
// #if USE_WIN32_API
//     #include <windows.h>
//     // Windows实现
// #elif USE_POSIX_API
//     #include <unistd.h>
//     #include <sys/stat.h>
//     // POSIX实现
// #endif
//
// void platform_specific_function() {
//     #if USE_WIN32_API
//         // Windows代码路径
//         Sleep(1000);  // Windows睡眠函数（毫秒）
//     #elif USE_POSIX_API
//         // Unix代码路径
//         usleep(1000000);  // POSIX睡眠函数（微秒）
//     #endif
// }
// ============================================================

#endif // FE_EM_CONFIG_H