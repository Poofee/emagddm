/**
 * @file em_solver_backends.hpp
 * @brief 核心数值层 - 直接求解器后端条件编译配置
 * @details 基于CMake编译期宏（HAVE_SUPERLU/HAVE_MUMPS）配置第三方直接求解后端的可用性，
 *          提供统一的宏接口用于条件编译控制。此文件应在所有直接求解器实现文件中首先包含。
 *
 * @par 编译配置说明：
 * - SuperLU支持：在CMakeLists.txt中通过 -DHAVE_SUPERLU 启用，需安装SuperLU库
 * - MUMPS支持：在CMakeLists.txt中通过 -DHAVE_MUMPS 启用，需安装MUMPS库（含MPI依赖）
 * - 默认仅启用Eigen内置求解器（无需额外依赖）
 *
 * @par 使用示例：
 * @code
 * #include "em_solver_backends.hpp"
 *
 * #if EM_SOLVER_HAS_SUPERLU
 *     // 使用SuperLU特定代码
 * #endif
 *
 * #if EM_SOLVER_HAS_MUMPS
 *     // 使用MUMPS特定代码
 * #endif
 * @endcode
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

// SuperLU 后端支持检测（在CMakeLists.txt中通过 target_compile_definitions 启用）
#ifdef HAVE_SUPERLU
    #define EM_SOLVER_HAS_SUPERLU 1
    // 包含SuperLU双精度头文件（根据实际安装路径通过CMake的target_include_directories提供）
    #include "slu_ddefs.h"
#else
    #define EM_SOLVER_HAS_SUPERLU 0
#endif

// MUMPS 后端支持检测（在CMakeLists.txt中通过 target_compile_definitions 启用）
#ifdef HAVE_MUMPS
    #define EM_SOLVER_HAS_MUMPS 1
    // 包含MUMPS双精度C接口头文件（根据实际安装路径通过CMake的target_include_directories提供）
    #include "dmumps_c.h"
#else
    #define EM_SOLVER_HAS_MUMPS 0
#endif
