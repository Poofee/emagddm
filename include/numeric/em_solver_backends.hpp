/**
 * @file em_solver_backends.hpp
 * @brief 核心数值层 - 直接求解器后端条件编译配置
 * @details 基于 CMake 编译期宏（HAVE_SUPERLU/HAVE_MUMPS）配置第三方直接求解后端的可用性，
 *          提供统一的宏接口用于条件编译控制。此文件应在所有直接求解器实现文件中首先包含。
 *
 * @par 编译配置说明：
 * - SuperLU_MT 支持：在 CMakeLists.txt 中通过 USE_SUPERLU_MT=ON 启用，使用 OpenMP 模式（避免 pthread 问题）
 * - MUMPS 支持：在 CMakeLists.txt 中通过 -DHAVE_MUMPS 启用，需安装 MUMPS 库（含 MPI 依赖）
 * - 默认仅启用 Eigen 内置求解器（无需额外依赖）
 *
 * @par 使用示例：
 * @code
 * #include "em_solver_backends.hpp"
 *
 * #ifdef HAVE_SUPERLU
 *     // 使用 SuperLU_MT 特定代码
 * #endif
 *
 * #ifdef HAVE_MUMPS
 *     // 使用 MUMPS 特定代码
 * #endif
 * @endcode
 *
 * @author Poofee
 * @date 2026-04-08
 * @version 1.3 (OpenMP mode - no pthread required)
 */

#pragma once

// SuperLU_MT 后端支持检测（在 CMakeLists.txt 中通过 target_compile_definitions 启用）
#ifdef HAVE_SUPERLU
    // 包含 SuperLU_MT 双精度头文件（OpenMP 并行 LU 分解）
    // slu_mt_ddefs.h 内部会依次包含：slu_mt_machines.h, slu_mt_Cnames.h, supermatrix.h,
    // slu_mt_util.h, pxgstrf_synch.h，提供完整的 SuperLU_MT API 声明
    #include "slu_mt_ddefs.h"
#endif

// MUMPS 后端支持检测（在 CMakeLists.txt 中通过 target_compile_definitions 启用）
#ifdef HAVE_MUMPS
    #include "dmumps_c.h"
    #include "zmumps_c.h"
#endif
