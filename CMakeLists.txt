# FETI-DP电磁有限元求解器 - 根CMake配置文件
# 负责全局编译选项、第三方库依赖、子目录包含、跨平台适配

cmake_minimum_required(VERSION 3.20)
project(FE_EM_FETIDP VERSION 1.0.0 LANGUAGES CXX C)

# ============================================================
# 跨平台操作系统检测与平台变量定义
# 使用CMake内置变量识别目标操作系统，为后续条件编译提供基础
# ============================================================

message(STATUS "==========================================")
message(STATUS "  跨平台编译环境检测")
message(STATUS "==========================================")
message(STATUS "  操作系统: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  系统版本: ${CMAKE_SYSTEM_VERSION}")
message(STATUS "  处理器架构: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  C++ 编译器: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# 操作系统类型检测（用于条件编译）
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(IS_WINDOWS TRUE)
    set(IS_LINUX FALSE)
    set(IS_MACOS FALSE)
    message(STATUS "  平台标识: Windows (Win32/Win64 API)")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(IS_WINDOWS FALSE)
    set(IS_LINUX TRUE)
    set(IS_MACOS FALSE)
    message(STATUS "  平台标识: Linux (POSIX接口)")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(IS_WINDOWS FALSE)
    set(IS_LINUX FALSE)
    set(IS_MACOS TRUE)
    message(STATUS "  平台标识: macOS (POSIX + Darwin API)")
else()
    message(WARNING "  未知的操作系统: ${CMAKE_SYSTEM_NAME}，将使用通用Unix配置")
    set(IS_WINDOWS FALSE)
    set(IS_LINUX TRUE)  # 默认使用Linux配置作为fallback
    set(IS_MACOS FALSE)
endif()

# 编译器家族检测（MSVC / GCC / Clang）
if(MSVC)
    set(COMPILER_FAMILY "MSVC")
    message(STATUS "  编译器家族: MSVC (Visual Studio)")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(COMPILER_FAMILY "GCC")
    message(STATUS "  编译器家族: GCC (GNU Compiler Collection)")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(COMPILER_FAMILY "Clang")
    message(STATUS "  编译器家族: Clang (LLVM/Apple Clang)")
else()
    set(COMPILER_FAMILY "UNKNOWN")
    message(WARNING "  未知的编译器: ${CMAKE_CXX_COMPILER_ID}")
endif()

# 架构检测（32位 vs 64位）
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(IS_64BIT TRUE)
    message(STATUS "  架构: 64位")
else()
    set(IS_64BIT FALSE)
    message(STATUS "  架构: 32位")
endif()

message(STATUS "")

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI/OpenMP 编译选项
option(USE_MPI "Enable MPI support" OFF)
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_SUPERLU_MT "Enable SuperLU_MT support" ON)

# ============================================================
# 平台特定的全局编译选项配置
# 根据不同操作系统和编译器设置警告级别、优化选项等
# ============================================================

message(STATUS "==========================================")
message(STATUS "  编译选项配置 (${COMPILER_FAMILY})")
message(STATUS "==========================================")

if(MSVC)
    # Windows/MSVC 特定编译选项
    # 原因说明：
    # - /W4: 启用Level 4警告，捕获更多潜在问题
    # - /utf-8: 强制源文件和执行字符集使用UTF-8，避免中文路径问题
    # - /wd4100: 禁用未引用形式参数警告（接口兼容性需要）
    # - /wd4101: 禁用未引用局部变量警告（调试代码常见）
    add_compile_options(/W4 /utf-8)
    add_compile_options(/wd4100 /wd4101)

    # MSVC安全特性：运行时类型检查和缓冲区安全
    add_compile_options(/RTC1)  # 运行时错误检查（堆栈初始化+未初始化变量检测）

    # OpenMP支持（MSVC原生支持，无需外部查找）
    if(USE_OPENMP)
        add_compile_options(/openmp)
        message(STATUS "  [Windows] OpenMP 已启用 (/openmp)")
    endif()

    # 预处理器定义：标识Windows平台
    add_compile_definitions(_WIN32_WINNT=0x0601)  # 目标Windows 7及以上版本
    add_compile_definitions(WIN32_LEAN_AND_MEAN)   # 减少Windows头文件包含范围，加速编译
    add_compile_definitions(NOMINMAX)              # 避免min/max宏与std::min/max冲突

    message(STATUS "  [Windows] MSVC编译选项已配置")

elseif(COMPILER_FAMILY STREQUAL "GCC")
    # Linux/GCC 特定编译选项
    # 原因说明：
    # - -Wall: 启用所有常用警告
    # - -Wextra: 启用额外警告（未使用参数等）
    # - -Wpedantic: 严格遵循ISO C++标准
    # - -fPIC: 生成位置无关代码（共享库必需）
    add_compile_options(-Wall -Wextra -Werror -pedantic)
    add_compile_options(-fPIC)

    # GCC调试信息格式（支持更好的调试体验）
    add_compile_options(-g3)  # 最大级别的调试信息

    # OpenMP支持（需要find_package定位OpenMP库）
    if(USE_OPENMP)
        find_package(OpenMP REQUIRED)
        if(OpenMP_CXX_FOUND)
            add_compile_options(${OpenMP_CXX_FLAGS})
            message(STATUS "  [Linux/GCC] OpenMP 已启用 (版本: ${OpenMP_CXX_VERSION})")
        else()
            message(FATAL_ERROR "  [Linux/GCC] OpenMP 未找到但被要求启用")
        endif()
    endif()

    # POSIX线程库（Linux下OpenMP和某些系统调用需要）
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    message(STATUS "  [Linux/GCC] GCC编译选项已配置")

elseif(COMPILER_FAMILY STREQUAL "Clang")
    # macOS/Clang 或 Linux/Clang 特定编译选项
    # 原因说明：
    # - Clang与GCC大部分标志兼容，但有一些独特优化选项
    # - -fcolor-diagnostics: 彩色诊断输出（提升可读性）
    # - -fPIC: macOS也推荐使用位置无关代码
    add_compile_options(-Wall -Wextra -Werror -pedantic)
    add_compile_options(-fPIC -fcolor-diagnostics)

    # Clang特有的静态分析友好选项
    add_compile_options(-Wdocumentation)  # 检查文档注释一致性

    # OpenMP支持
    if(USE_OPENMP)
        find_package(OpenMP REQUIRED)
        if(OpenMP_CXX_FOUND)
            add_compile_options(${OpenMP_CXX_FLAGS})
            message(STATUS "  [macOS/Clang] OpenMP 已启用 (版本: ${OpenMP_CXX_VERSION})")
        else()
            message(FATAL_ERROR "  [macOS/Clang] OpenMP 未找到但被要求启用")
        endif()
    endif()

    # macOS特定：线程库处理
    if(IS_MACOS)
        set(THREADS_PREFER_PTHREAD_FLAG ON)
        find_package(Threads REQUIRED)
        message(STATUS "  [macOS/Clang] macOS Clang编译选项已配置")
    else()
        find_package(Threads REQUIRED)
        message(STATUS "  [Linux/Clang] Linux Clang编译选项已配置")
    endif()

else()
    # 未知编译器的fallback配置（保守策略）
    message(WARNING "  使用保守的通用编译选项")
    add_compile_options(-Wall)
    if(USE_OPENMP)
        find_package(OpenMP)
        if(OpenMP_FOUND)
            add_compile_options(${OpenMP_CXX_FLAGS})
        endif()
    endif()
endif()

message(STATUS "")

# ============================================================
# MPI支持（跨平台配置）
# Windows: 使用Microsoft MPI SDK（预安装在Visual Studio中）
# Linux/macOS: 使用系统安装的OpenMPI或MPICH
# ============================================================

if(USE_MPI)
    message(STATUS "==========================================")
    message(STATUS "  MPI 配置 (跨平台)")
    message(STATUS "==========================================")

    if(IS_WINDOWS)
        # Windows平台：使用Microsoft MPI (MS-MPI)
        # 原因说明：
        # - MS-MPI是Windows上最常用的MPI实现，与MSVC完美兼容
        # - 通常安装在"Program Files"下的Microsoft SDKs目录
        # - 链接库为msmpi.lib（动态链接）或msmpifec.lib（静态链接）
        set(MPI_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/Microsoft SDKs/MPI/Include")

        if(EXISTS "${MPI_INCLUDE_PATH}/mpi.h")
            message(STATUS "  [Windows] 找到本地MPI库: ${MPI_INCLUDE_PATH}")
            add_definitions(-DUSE_MPI=1)
            include_directories(${MPI_INCLUDE_PATH})

            # Windows下链接MS-MPI运行时库
            set(MPI_LIBRARIES "msmpi.lib")
            link_libraries(${MPI_LIBRARIES})
            message(STATUS "  [Windows] MPI库: msmpi.lib (Microsoft MPI)")
        else()
            message(FATAL_ERROR "  [Windows] MPI头文件未找到: ${MPI_INCLUDE_PATH}/mpi.h")
        endif()

    elseif(IS_LINUX OR IS_MACOS)
        # Linux/macOS平台：使用CMake模块自动查找系统MPI
        # 原因说明：
        # - Linux发行版通常通过包管理器安装OpenMPI或MPICH
        # - macOS可通过Homebrew安装：brew install open-mpi
        # - CMake的FindMPI模块能自动检测常见安装位置
        find_package(MPI REQUIRED)

        if(MPI_FOUND)
            message(STATUS "  [Unix] MPI编译器: ${MPI_CXX_COMPILER}")
            message(STATUS "  [Unix] MPI版本: ${MPI_VERSION}")
            message(STATUS "  [Unix] MPI包含目录: ${MPI_CXX_INCLUDE_DIRS}")
            add_definitions(-DUSE_MPI=1)

            # 将MPI头文件和库路径添加到全局范围
            include_directories(${MPI_CXX_INCLUDE_DIRS})
            link_libraries(${MPI_CXX_LIBRARIES})

            # 设置MPI_ROOT供config.h使用
            if(MPI_CXX_INCLUDE_DIRS)
                list(GET MPI_CXX_INCLUDE_DIRS 0 MPI_FIRST_INC_DIR)
                get_filename_component(MPI_ROOT "${MPI_FIRST_INC_DIR}" DIRECTORY)
                add_definitions("-DMPI_ROOT=\"${MPI_ROOT}\"")
                message(STATUS "  [Unix] MPI根目录: ${MPI_ROOT}")
            endif()
        else()
            message(FATAL_ERROR "  [Unix] 未找到MPI实现。请安装OpenMPI或MPICH后重试。")
        endif()

    else()
        message(FATAL_ERROR "  不支持的操作系统，无法配置MPI")
    endif()

    message(STATUS "")
else()
    add_definitions(-DUSE_MPI=0)
    add_definitions(-DMPI_ROOT="")
endif()

# OpenMP编译定义（全局宏，供C++源代码使用）
if(USE_OPENMP)
    add_definitions(-DUSE_OPENMP=1)
else()
    add_definitions(-DUSE_OPENMP=0)
endif()

# ============================================================
# 平台特定的输出目录配置
# 原因说明：
# - Windows: 输出到bin/lib目录，符合Visual Studio惯例
# - Linux/macOS: 可选支持标准Unix目录结构（usr/local）
# - 统一使用CMAKE_BINARY_DIR确保构建产物在build树内
# ============================================================

message(STATUS "==========================================")
message(STATUS "  输出目录配置")
message(STATUS "==========================================")

# 基础输出目录（所有平台通用）
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 平台特定的子目录组织（可选，用于区分调试/发布模式）
if(IS_WINDOWS)
    # Windows下按配置类型分子目录（Debug/Release/RelWithDebInfo）
    foreach(outputconfig IN ITEMS CMAKE_RUNTIME_OUTPUT_DIRECTORY CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set(${outputconfig}_DEBUG "${${outputconfig}}/Debug")
        set(${outputconfig}_RELEASE "${${outputconfig}}/Release")
        set(${outputconfig}_RELWITHDEBINFO "${${outputconfig}}/RelWithDebInfo")
        set(${outputconfig}_MINSIZEREL "${${outputconfig}}/MinSizeRel")
    endforeach()
    message(STATUS "  [Windows] 输出目录: ${CMAKE_BINARY_DIR}/bin/<Configuration>")
else()
    # Unix系统通常不区分子目录，保持简洁
    message(STATUS "  [Unix] 输出目录:")
    message(STATUS "    可执行文件 -> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    message(STATUS "    共享库     -> ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    message(STATUS "    静态库     -> ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
endif()

message(STATUS "")

# ============================================================
# 第三方库配置（必须在test之前，因为test依赖这些库的编译定义）
# ============================================================

# 单头文件库配置
include_directories(
    ${CMAKE_SOURCE_DIR}/lib/eigen
    ${CMAKE_SOURCE_DIR}/lib/pugixml
    ${CMAKE_SOURCE_DIR}/lib/spdlog/include
    ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

# spdlog库配置（需要编译部分源文件，含bundled fmt库）
add_library(spdlog_lib
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/spdlog.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/cfg.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/color_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/file_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/stdout_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/async.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/bundled_fmtlib_format.cpp
)

target_compile_definitions(spdlog_lib PRIVATE SPDLOG_COMPILED_LIB)

target_include_directories(spdlog_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

# pugixml库配置
add_library(pugixml_lib
    ${CMAKE_SOURCE_DIR}/lib/pugixml/src/pugixml.cpp
)

target_include_directories(pugixml_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

# 设置编译选项（跨平台适配）
if(MSVC)
    target_compile_options(spdlog_lib PRIVATE /W4)
    target_compile_options(spdlog_lib PRIVATE /O2)
elseif(COMPILER_FAMILY STREQUAL "GCC")
    target_compile_options(spdlog_lib PRIVATE -Wall -Wextra)
    target_compile_options(spdlog_lib PRIVATE -O2)
elseif(COMPILER_FAMILY STREQUAL "Clang")
    target_compile_options(spdlog_lib PRIVATE -Wall -Wextra)
    target_compile_options(spdlog_lib PRIVATE -O2)
else()
    target_compile_options(spdlog_lib PRIVATE -Wall)
endif()

# ============================================================
# SuperLU_MT 库配置（多线程稀疏 LU 分解直接求解器）
# 使用 OpenMP 模式，避免 Windows pthread 兼容性问题
# ============================================================
if(USE_SUPERLU_MT)
    message(STATUS "Configuring SuperLU_MT with OpenMP...")

    # SuperLU_MT 编译选项设置
    set(PLAT "_OPENMP" CACHE STRING "SuperLU_MT threading flavor (OpenMP)" FORCE)

    # 禁用内置 CBLAS，使用外部 OpenBLAS（避免 35 个重复 BLAS 符号错误）
    # 条件：仅当 OpenBLAS 启用时才禁用内部 CBLAS
    if(USE_OPENBLAS)
        set(enable_internal_blaslib OFF CACHE BOOL "Use external OpenBLAS for SuperLU_MT" FORCE)
        # 直接计算 OpenBLAS 静态库路径（不依赖 target，因为 OpenBLAS 在后面配置）
        set(_OPENBLAS_LIB_PATH "${CMAKE_BINARY_DIR}/openblas_build/lib/libopenblas.a")
        if(EXISTS "${_OPENBLAS_LIB_PATH}")
            set(TPL_BLAS_LIBRARIES "${_OPENBLAS_LIB_PATH}" CACHE FILEPATH "External BLAS for SuperLU_MT" FORCE)
            message(STATUS "  SuperLU_MT: 使用外部 OpenBLAS (${_OPENBLAS_LIB_PATH})")
        else()
            # 回退到内置 CBLAS（OpenBLAS 尚未编译）
            set(enable_internal_blaslib ON CACHE BOOL "Build internal CBLAS for SuperLU_MT" FORCE)
            message(STATUS "  SuperLU_MT: OpenBLAS 未找到，回退到内置 CBLAS")
        endif()
    else()
        set(enable_internal_blaslib ON CACHE BOOL "Build internal CBLAS for SuperLU_MT" FORCE)
        message(STATUS "  SuperLU_MT: 使用内置 CBLAS")
    endif()

    set(enable_double ON CACHE BOOL "Enable double precision" FORCE)
    set(enable_single OFF CACHE BOOL "Disable single precision" FORCE)
    set(enable_complex OFF CACHE BOOL "Disable complex precision" FORCE)
    set(enable_complex16 OFF CACHE BOOL "Disable complex16 precision" FORCE)
    set(enable_matlabmex OFF CACHE BOOL "Disable MATLAB mex" FORCE)
    set(enable_examples OFF CACHE BOOL "Disable examples" FORCE)
    set(enable_tests OFF CACHE BOOL "Disable SuperLU_MT tests" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static library" FORCE)

    # 启用 OpenMP
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        message(STATUS "  OpenMP found: ${OpenMP_C_VERSION}")
    else()
        message(FATAL_ERROR "OpenMP is required for SuperLU_MT but not found")
    endif()

    # 添加 SuperLU_MT 子目录
    add_subdirectory(${CMAKE_SOURCE_DIR}/lib/superlu_mt ${CMAKE_BINARY_DIR}/superlu_mt)

    # 获取 SuperLU_MT 库目标名称（根据 PLAT 变量确定）
    set(SUPERLU_MT_TARGET superlu_mt${PLAT})

    # 检查目标是否存在
    if(TARGET ${SUPERLU_MT_TARGET})
        message(STATUS "  SuperLU_MT target found: ${SUPERLU_MT_TARGET}")

        # 添加 SuperLU_MT 头文件路径到全局包含目录
        include_directories(${CMAKE_SOURCE_DIR}/lib/superlu_mt/SRC)

        # 为 numeric_lib 定义 HAVE_SUPERLU 宏并链接 SuperLU_MT
        # 注意：在 numeric_lib 目标创建后（src/CMakeLists.txt 中）会执行链接操作
        set(SUPERLU_MT_AVAILABLE TRUE)
        message(STATUS "  SuperLU_MT configured successfully (OpenMP mode)")
    else()
        message(WARNING "  SuperLU_MT target '${SUPERLU_MT_TARGET}' not found, disabling")
        set(SUPERLU_MT_AVAILABLE FALSE)
    endif()
else()
    set(SUPERLU_MT_AVAILABLE FALSE)
    message(STATUS "SuperLU_MT support disabled by user option")
endif()

# ============================================================
# OpenBLAS 配置（高性能 BLAS 库，MUMPS 依赖）
# 跨平台支持：
# - Windows:   使用预编译的静态库
# - Linux:     使用系统包管理器或从源码编译
# - macOS:     从源码编译（针对 Apple Silicon 优化）
# ============================================================
option(USE_OPENBLAS "从源码编译 OpenBLAS (MUMPS依赖)" OFF)

if(USE_OPENBLAS)
    message(STATUS "==========================================")
    message(STATUS "  OpenBLAS 配置 (跨平台)")
    message(STATUS "==========================================")

    set(OPENBLAS_ROOT "${CMAKE_SOURCE_DIR}/lib/OpenBLAS")
    set(OPENBLAS_BUILD_DIR "${CMAKE_BINARY_DIR}/openblas_build")

    if(IS_MACOS AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        # ============================================================
        # macOS Apple Silicon (arm64) 特殊配置
        # 使用 OpenBLAS 自带的 Makefile.arm64 进行优化编译
        # ============================================================
        message(STATUS "  [macOS/arm64] 检测到 Apple Silicon，使用优化的 ARM64 构建")

        # 检查是否已编译
        set(OPENBLAS_STATIC_LIB "${OPENBLAS_BUILD_DIR}/lib/libopenblas.a")

        if(NOT EXISTS ${OPENBLAS_STATIC_LIB})
            message(STATUS "  [macOS] 开始编译 OpenBLAS (ARM64 优化版本)...")

            # 创建构建目录
            file(MAKE_DIRECTORY ${OPENBLAS_BUILD_DIR})

            # 设置编译参数（针对 Apple Silicon M1/M2/M3/M4 优化）
            find_program(GFORTRAN_EXECUTABLE gfortran)
            if(NOT GFORTRAN_EXECUTABLE)
                message(FATAL_ERROR "OpenBLAS 需要 Fortran 编译器 (gfortran)。请安装: brew install gcc")
            endif()

            # 获取 CPU 核心数用于并行编译
            execute_process(
                COMMAND sysctl -n hw.ncpu
                OUTPUT_VARIABLE NUM_CORES
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            set(OPENBLAS_MAKE_ARGS
                TARGET=ARMV8
                CC=${CMAKE_C_COMPILER}
                FC=${GFORTRAN_EXECUTABLE}
                PREFIX=${OPENBLAS_BUILD_DIR}
                NO_LAPACK=0
                NO_LAPACKE=0
                BUILD_LAPACK_DEPRECATED=1
                USE_OPENMP=0
                USE_THREAD=0
                NUM_THREADS=${NUM_CORES}
                DYNAMIC_ARCH=0
            )

            # 执行 make 编译（使用多核并行编译）
            execute_process(
                COMMAND make -j${NUM_CORES} ${OPENBLAS_MAKE_ARGS} all
                WORKING_DIRECTORY ${OPENBLAS_ROOT}
                RESULT_VARIABLE OPENBLAS_MAKE_RESULT
                OUTPUT_VARIABLE OPENBLAS_MAKE_OUTPUT
                ERROR_VARIABLE OPENBLAS_MAKE_ERROR
                TIMEOUT 900
            )

            if(NOT OPENBLAS_MAKE_RESULT EQUAL 0)
                message(FATAL_ERROR "OpenBLAS 编译失败: ${OPENBLAS_MAKE_ERROR}")
            endif()

            # 安装到构建目录
            execute_process(
                COMMAND make ${OPENBLAS_MAKE_ARGS} install
                WORKING_DIRECTORY ${OPENBLAS_ROOT}
                RESULT_VARIABLE OPENBLAS_INSTALL_RESULT
                OUTPUT_VARIABLE OPENBLAS_INSTALL_OUTPUT
                ERROR_VARIABLE OPENBLAS_INSTALL_ERROR
                TIMEOUT 300
            )

            if(NOT OPENBLAS_INSTALL_RESULT EQUAL 0)
                message(FATAL_ERROR "OpenBLAS 安装失败: ${OPENBLAS_INSTALL_ERROR}")
            endif()

            message(STATUS "  [macOS] OpenBLAS 编译完成")
        else()
            message(STATUS "  [macOS] 使用已编译的 OpenBLAS 库")
        endif()

        # 验证库文件存在
        if(EXISTS ${OPENBLAS_STATIC_LIB})
            message(STATUS "  [macOS] OpenBLAS 静态库: ${OPENBLAS_STATIC_LIB}")

            # 导入预编译的静态库
            add_library(openblas STATIC IMPORTED GLOBAL)
            set_target_properties(openblas PROPERTIES
                IMPORTED_LOCATION "${OPENBLAS_STATIC_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${OPENBLAS_BUILD_DIR}/include"
            )

            set(OPENBLAS_AVAILABLE TRUE)
            message(STATUS "  [macOS] OpenBLAS 配置完成 (ARM64 优化静态库)")
        else()
            message(WARNING "  [macOS] OpenBLAS 编译后未找到库文件: ${OPENBLAS_STATIC_LIB}")
            set(OPENBLAS_AVAILABLE FALSE)
        endif()

    elseif(IS_LINUX OR IS_MACOS)
        # ============================================================
        # Linux 或 macOS x86_64：使用 CMake 构建或查找系统库
        # ============================================================
        message(STATUS "  [Unix] 尝试从源码编译 OpenBLAS...")

        option(BUILD_OPENBLAS_FROM_SOURCE "从源码编译 OpenBLAS" ON)

        if(BUILD_OPENBLAS_FROM_SOURCE)
            # 使用 CMake 构建 OpenBLAS
            set(OPENBLAS_CMAKE_ARGS
                -DCMAKE_INSTALL_PREFIX=${OPENBLAS_BUILD_DIR}
                -DCMAKE_BUILD_TYPE=Release
                -DBUILD_WITHOUT_LAPACK=OFF
                -DBUILD_WITHOUT_LAPACKE=OFF
                -DBUILD_STATIC_LIBS=ON
                -DBUILD_SHARED_LIBS=OFF
                -DNO_WARMUP=ON
                -DUSE_LOCKING=OFF
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            )

            if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
                list(APPEND OPENBLAS_CMAKE_ARGS
                    -DDYNAMIC_ARCH=ON
                    -DTARGET=CORE2
                )
            elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
                list(APPEND OPENBLAS_CMAKE_ARGS
                    -DDYNAMIC_ARCH=OFF
                    -DTARGET=ARMV8
                )
            endif()

            # 添加子目录编译
            add_subdirectory(${OPENBLAS_ROOT} ${OPENBLAS_BUILD_DIR} EXCLUDE_FROM_ALL)

            # OpenBLAS CMake 目标名称为 openblas（静态库）
            if(TARGET openblas_static)
                set_target_properties(openblas_static PROPERTIES GLOBAL TRUE)
                set(OPENBLAS_TARGET openblas_static)
            elseif(TARGET openblas)
                set_target_properties(openblas PROPERTIES GLOBAL TRUE)
                set(OPENBLAS_TARGET openblas)
            endif()

            set(OPENBLAS_AVAILABLE TRUE)
            message(STATUS "  [Unix] OpenBLAS 从源码编译完成")
        else()
            # 查找系统安装的 OpenBLAS
            find_package(BLAS QUIET)
            if(BLAS_FOUND)
                message(STATUS "  [Unix] 使用系统 BLAS: ${BLAS_LIBRARIES}")
                set(OPENBLAS_AVAILABLE TRUE)
            else()
                message(WARNING "  [Unix] 未找到系统 BLAS/OpenBLAS")
                set(OPENBLAS_AVAILABLE FALSE)
            endif()
        endif()

    else()
        # 其他平台：尝试查找预编译库
        set(OPENBLAS_PREBUILT_LIB "${OPENBLAS_ROOT}/lib/libopenblas.a")
        if(EXISTS ${OPENBLAS_PREBUILT_LIB})
            add_library(openblas STATIC IMPORTED GLOBAL)
            set_target_properties(openblas PROPERTIES
                IMPORTED_LOCATION "${OPENBLAS_PREBUILT_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${OPENBLAS_ROOT}/include"
            )
            set(OPENBLAS_AVAILABLE TRUE)
            message(STATUS "  [Other] 使用预编译 OpenBLAS 库")
        else()
            set(OPENBLAS_AVAILABLE FALSE)
            message(STATUS "  OpenBLAS 在此平台不可用")
        endif()
    endif()

    # 全局状态变量供其他组件使用
    if(OPENBLAS_AVAILABLE)
        set(OPENBLAS_ROOT_DIR ${OPENBLAS_BUILD_DIR} CACHE PATH "OpenBLAS installation root" FORCE)
        message(STATUS "  ✓ OpenBLAS 可用")
    else()
        message(STATUS "  ✗ OpenBLAS 不可用")
    endif()

else()
    set(OPENBLAS_AVAILABLE FALSE)
    message(STATUS "OpenBLAS support disabled (use -DUSE_OPENBLAS=ON to enable)")
endif()

# ============================================================
# MUMPS 配置（可选，需预编译 MUMPS + OpenBLAS）
# 跨平台支持：
# - Windows: 使用 build_windows.ps1 预编译的静态库 (.lib)
# - Linux:   从源码编译或使用系统包管理器
# - macOS:   从源码编译（依赖 OpenBLAS）
# ============================================================
option(USE_MUMPS "启用 MUMPS 并行直接求解器后端" OFF)
option(BUILD_MUMPS_FROM_SOURCE "从源码编译 MUMPS (需要 OpenBLAS)" ON)

if(USE_MUMPS)
    message(STATUS "==========================================")
    message(STATUS "  MUMPS 5.8.2 配置 (跨平台)")
    message(STATUS "==========================================")

    # 设置MUMPS根目录
    set(MUMPS_ROOT "${CMAKE_SOURCE_DIR}/lib/MUMPS")

    if(IS_WINDOWS)
        # Windows平台：使用预编译的静态库
        # 原因说明：
        # - MUMPS依赖Fortran编译器，Windows上通常使用Intel Fortran或MinGW-gfortran预编译
        # - OpenBLAS同样需要预编译，避免在Windows上从源码构建的复杂性
        # - .lib格式是MSVC标准静态库格式

        set(OPENBLAS_ROOT "${CMAKE_SOURCE_DIR}/lib/OpenBLAS")

        # 检查Windows预编译库是否存在
        set(MUMPS_LIBS
            "${MUMPS_ROOT}/lib/mumps.lib"
            "${MUMPS_ROOT}/lib/pord.lib"
            "${MUMPS_ROOT}/lib/mpiseq.lib"
        )
        set(MUMPS_LIBS_FOUND TRUE)
        foreach(_lib ${MUMPS_LIBS})
            if(NOT EXISTS ${_lib})
                message(WARNING "  [Windows] MUMPS预编译库未找到: ${_lib}")
                set(MUMPS_LIBS_FOUND FALSE)
            endif()
        endforeach()

        # 检查OpenBLAS Windows版本
        if(EXISTS "${OPENBLAS_ROOT}/lib/libopenblas.dll.a")
            set(OPENBLAS_LIB "${OPENBLAS_ROOT}/lib/libopenblas.dll.a")
            message(STATUS "  [Windows] OpenBLAS找到 (导入库): ${OPENBLAS_LIB}")
        elseif(EXISTS "${OPENBLAS_ROOT}/lib/libopenblas.a")
            set(OPENBLAS_LIB "${OPENBLAS_ROOT}/lib/libopenblas.a")
            message(STATUS "  [Windows] OpenBLAS找到 (静态库): ${OPENBLAS_LIB}")
        else()
            message(WARNING "  [Windows] OpenBLAS库未找到 in ${OPENBLAS_ROOT}/lib")
            set(MUMPS_LIBS_FOUND FALSE)
        endif()

        if(MUMPS_LIBS_FOUND)
            # 导入预编译的MUMPS静态库（IMPORTED目标）
            add_library(mumps STATIC IMPORTED)
            set_target_properties(mumps PROPERTIES
                IMPORTED_LOCATION "${MUMPS_ROOT}/lib/mumps.lib"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/include;${MUMPS_ROOT}/libseq"
            )

            add_library(pord STATIC IMPORTED)
            set_target_properties(pord PROPERTIES
                IMPORTED_LOCATION "${MUMPS_ROOT}/lib/pord.lib"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/PORD/include"
            )

            add_library(mpiseq STATIC IMPORTED)
            set_target_properties(mpiseq PROPERTIES
                IMPORTED_LOCATION "${MUMPS_ROOT}/lib/mpiseq.lib"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/libseq"
            )

            add_library(openblas STATIC IMPORTED)
            set_target_properties(openblas PROPERTIES
                IMPORTED_LOCATION "${OPENBLAS_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${OPENBLAS_ROOT}/include"
            )

            # Intel Fortran运行时库（MUMPS Fortran组件需要）
            # 原因说明：MUMPS部分核心用Fortran编写，需要Intel Fortran运行时支持
            set(INTEL_FORTRAN_LIB_DIR "C:/Program Files (x86)/Intel/oneAPI/compiler/2023.2.0/windows/compiler/lib/intel64")
            if(EXISTS ${INTEL_FORTRAN_LIB_DIR})
                message(STATUS "  [Windows] Intel Fortran运行时目录: ${INTEL_FORTRAN_LIB_DIR}")
            else()
                message(WARNING "  [Windows] Intel Fortran运行时未找到: ${INTEL_FORTRAN_LIB_DIR}（如果链接失败请安装）")
            endif()

            set(MUMPS_AVAILABLE TRUE)
            include_directories(${MUMPS_ROOT}/include)
            include_directories(${MUMPS_ROOT}/libseq)
            message(STATUS "  [Windows] MUMPS配置完成 (预编译静态库)")
        else()
            set(MUMPS_AVAILABLE FALSE)
            message(WARNING "  [Windows] MUMPS预编译库未找到。请先运行 build_windows.ps1 构建依赖库。")
        endif()

    elseif(IS_LINUX OR IS_MACOS)
        # Unix平台：支持从源码编译或使用系统库
        # 策略优先级：
        # 1. 从源码编译 MUMPS (BUILD_MUMPS_FROM_SOURCE=ON, 推荐)
        # 2. 查找系统安装的 MUMPS 库

        if(BUILD_MUMPS_FROM_SOURCE AND OPENBLAS_AVAILABLE)
            # ============================================================
            # 从源码编译 MUMPS（推荐方式）
            # 需要依赖：OpenBLAS (已通过 USE_OPENBLAS=ON 编译)
            #          gfortran (Fortran 编译器)
            # ============================================================
            message(STATUS "  [Unix] 从源码编译 MUMPS...")

            # 检查 Fortran 编译器
            find_program(GFORTRAN_EXECUTABLE gfortran)
            if(NOT GFORTRAN_EXECUTABLE)
                message(FATAL_ERROR "MUMPS 需要 Fortran 编译器 (gfortran)。请安装: brew install gcc")
            endif()
            message(STATUS "    [Unix] Fortran 编译器: ${GFORTRAN_EXECUTABLE}")

            # 设置 OpenBLAS 路径供 MUMPS 使用
            set(MUMPS_OPENBLAS_DIR "${OPENBLAS_ROOT_DIR}")
            if(NOT MUMPS_OPENBLAS_DIR OR NOT EXISTS "${MUMPS_OPENBLAS_DIR}/lib/libopenblas.a")
                # 回退到预编译位置
                set(MUMPS_OPENBLAS_DIR "${CMAKE_BINARY_DIR}/openblas_build")
                if(NOT EXISTS "${MUMPS_OPENBLAS_DIR}/lib/libopenblas.a")
                    message(FATAL_ERROR "MUMPS 需要 OpenBLAS，但未找到编译好的库文件")
                endif()
            endif()

            message(STATUS "    [Unix] 使用 OpenBLAS: ${MUMPS_OPENBLAS_DIR}")

            # 配置 MUMPS CMake 参数
            set(MUMPS_CMAKE_ARGS
                -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/mumps_build
                -DCMAKE_BUILD_TYPE=Release
                -DMUMPS_USE_SEQMPI=ON
                -DMUMPS_BUILD_PORD=ON
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_Fortran_COMPILER=${GFORTRAN_EXECUTABLE}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DOpenBLAS_ROOT=${MUMPS_OPENBLAS_DIR}
            )

            # 添加子目录编译 MUMPS
            set(MUMPS_BUILD_DIR "${CMAKE_BINARY_DIR}/mumps_src_build")

            # 获取 CPU 核心数用于并行编译
            execute_process(
                COMMAND sysctl -n hw.ncpu
                OUTPUT_VARIABLE MUMPS_NUM_CORES
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            # 创建自定义目标来跟踪 MUMPS 构建
            add_custom_target(mumps_build ALL
                COMMAND ${CMAKE_COMMAND} ${MUMPS_CMAKE_ARGS} -S ${MUMPS_ROOT} -B ${MUMPS_BUILD_DIR}
                COMMAND ${CMAKE_COMMAND} --build ${MUMPS_BUILD_DIR} --config Release -j${MUMPS_NUM_CORES}
                COMMENT "Building MUMPS from source..."
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                BYPRODUCTS
                    ${MUMPS_BUILD_DIR}/lib/libmumps_common.a
                    ${MUMPS_BUILD_DIR}/lib/libdmumps.a
                    ${MUMPS_BUILD_DIR}/lib/libzmumps.a
                    ${MUMPS_BUILD_DIR}/lib/libpord.a
                    ${MUMPS_BUILD_DIR}/lib/libmpiseq.a
            )

            # 导入 MUMPS 三层库架构
            # mumps_common: 公共 Fortran/C 模块（只链接一次，避免重复符号）
            add_library(mumps_common STATIC IMPORTED GLOBAL)
            set_target_properties(mumps_common PROPERTIES
                IMPORTED_LOCATION "${MUMPS_BUILD_DIR}/lib/libmumps_common.a"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/include;${MUMPS_ROOT}/libseq"
            )
            add_dependencies(mumps_common mumps_build)

            # dmups: 实数 MUMPS（依赖 mumps_common）
            add_library(dmumps STATIC IMPORTED GLOBAL)
            set_target_properties(dmumps PROPERTIES
                IMPORTED_LOCATION "${MUMPS_BUILD_DIR}/lib/libdmumps.a"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/include;${MUMPS_ROOT}/libseq"
            )
            add_dependencies(dmumps mumps_build)

            # zmumps: 复数 MUMPS（依赖 mumps_common）
            add_library(zmumps STATIC IMPORTED GLOBAL)
            set_target_properties(zmumps PROPERTIES
                IMPORTED_LOCATION "${MUMPS_BUILD_DIR}/lib/libzmumps.a"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/include;${MUMPS_ROOT}/libseq"
            )
            add_dependencies(zmumps mumps_build)

            add_library(pord_mumps STATIC IMPORTED GLOBAL)
            set_target_properties(pord_mumps PROPERTIES
                IMPORTED_LOCATION "${MUMPS_BUILD_DIR}/lib/libpord.a"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/PORD/include"
            )
            add_dependencies(pord_mumps mumps_build)

            add_library(mpiseq STATIC IMPORTED GLOBAL)
            set_target_properties(mpiseq PROPERTIES
                IMPORTED_LOCATION "${MUMPS_BUILD_DIR}/lib/libmpiseq.a"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/libseq"
            )
            add_dependencies(mpiseq mumps_build)

            # MUMPS 组件集合（三层架构：common + dmumps + zmumps + 依赖）
            add_library(mumps_components INTERFACE IMPORTED GLOBAL)
            target_link_libraries(mumps_components INTERFACE mumps_common dmumps zmumps pord_mumps mpiseq openblas)

            set(MUMPS_AVAILABLE TRUE)
            include_directories(${MUMPS_ROOT}/include)
            include_directories(${MUMPS_ROOT}/libseq)
            message(STATUS "  [Unix] MUMPS 从源码编译完成 (三层库架构)")
            message(STATUS "    [Unix] MUMPS 公共库: ${MUMPS_BUILD_DIR}/lib/libmumps_common.a")
            message(STATUS "    [Unix] MUMPS 实数库: ${MUMPS_BUILD_DIR}/lib/libdmumps.a")
            message(STATUS "    [Unix] MUMPS 复数库: ${MUMPS_BUILD_DIR}/lib/libzmumps.a")

        else()
            # ============================================================
            # 查找系统安装的 MUMPS（回退方案）
            # ============================================================
            message(STATUS "  [Unix] 尝试查找系统 MUMPS 库...")

            find_library(MUMPS_SYSTEM_LIB
                NAMES mumps mumps_common dmumps
                PATHS
                    /usr/local/lib
                    /usr/lib
                    /opt/homebrew/lib
                    ${MUMPS_ROOT}/lib
                    ${CMAKE_BINARY_DIR}/lib
                NO_DEFAULT_PATH
            )

            find_library(PORD_SYSTEM_LIB
                NAMES pord
                PATHS
                    /usr/local/lib
                    /usr/lib
                    /opt/homebrew/lib
                    ${MUMPS_ROOT}/lib
                    ${CMAKE_BINARY_DIR}/lib
                NO_DEFAULT_PATH
            )

            find_library(MPISEQ_SYSTEM_LIB
                NAMES mpiseq
                PATHS
                    /usr/local/lib
                    /usr/lib
                    /opt/homebrew/lib
                    ${MUMPS_ROOT}/lib
                    ${CMAKE_BINARY_DIR}/lib
                NO_DEFAULT_PATH
            )

            # 查找OpenBLAS或其他BLAS实现
            if(OPENBLAS_AVAILABLE AND TARGET openblas)
                set(OPENBLAS_SYSTEM_LIB openblas)
                message(STATUS "    [优先] 使用从源码编译的 OpenBLAS (目标: openblas)")
            else()
                find_library(OPENBLAS_SYSTEM_LIB
                    NAMES openblas blas cblas
                    PATHS
                        /usr/local/lib
                        /usr/lib
                        /opt/homebrew/lib
                        /opt/homebrew/opt/openblas/lib
                        ${CMAKE_SOURCE_DIR}/lib/OpenBLAS/lib
                        ${CMAKE_BINARY_DIR}/openblas_build/lib
                    NO_DEFAULT_PATH
                )
            endif()

            # 检查是否所有必需库都找到了
            if(MUMPS_SYSTEM_LIB AND PORD_SYSTEM_LIB AND MPISEQ_SYSTEM_LIB AND OPENBLAS_SYSTEM_LIB)
                message(STATUS "  [Unix] 找到MUMPS库:")
                message(STATUS "    MUMPS:  ${MUMPS_SYSTEM_LIB}")
                message(STATUS "    PORD:   ${PORD_SYSTEM_LIB}")
                message(STATUS "    MPISEQ: ${MPISEQ_SYSTEM_LIB}")
                message(STATUS "    BLAS:   ${OPENBLAS_SYSTEM_LIB}")

                # 创建导入目标
                add_library(mumps SHARED IMPORTED)
                set_target_properties(mumps PROPERTIES
                    IMPORTED_LOCATION "${MUMPS_SYSTEM_LIB}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/include;${MUMPS_ROOT}/libseq"
                )

                add_library(pord SHARED IMPORTED)
                set_target_properties(pord PROPERTIES
                    IMPORTED_LOCATION "${PORD_SYSTEM_LIB}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/PORD/include"
                )

                add_library(mpiseq SHARED IMPORTED)
                set_target_properties(mpiseq PROPERTIES
                    IMPORTED_LOCATION "${MPISEQ_SYSTEM_LIB}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_ROOT}/libseq"
                )

                if(NOT TARGET openblas)
                    add_library(openblas SHARED IMPORTED)
                    set_target_properties(openblas PROPERTIES
                        IMPORTED_LOCATION "${OPENBLAS_SYSTEM_LIB}"
                    )
                endif()

                set(MUMPS_AVAILABLE TRUE)
                include_directories(${MUMPS_ROOT}/include)
                include_directories(${MUMPS_ROOT}/libseq)
                message(STATUS "  [Unix] MUMPS配置完成 (系统共享库)")

                if(IS_MACOS)
                    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
                    message(STATUS "  [macOS] 已启用RPATH自动设置")
                endif()

            else()
                # 未找到系统库时的提示信息
                set(MUMPS_AVAILABLE FALSE)
                if(IS_LINUX)
                    message(WARNING "  [Linux] 未找到完整的MUMPS库。请执行以下操作之一：")
                    message(WARNING "    1. 启用源码编译: cmake .. -DUSE_MUMPS=ON -DUSE_OPENBLAS=ON")
                    message(WARNING "    2. Ubuntu/Debian: sudo apt-get install libmumps-dev libopenblas-dev")
                    message(WARNING "    3. CentOS/RHEL: sudo yum install mumps-openblas-devel")
                elseif(IS_MACOS)
                    message(WARNING "  [macOS] 未找到完整的MUMPS库。请执行以下操作：")
                    message(WARNING "    1. 启用源码编译: cmake .. -DUSE_MUMPS=ON -DUSE_OPENBLAS=ON -DBUILD_MUMPS_FROM_SOURCE=ON")
                    message(WARNING "    2. 确保 gfortran 已安装: brew install gcc")
                endif()
            endif()
        endif()

    else()
        set(MUMPS_AVAILABLE FALSE)
        message(STATUS "  MUMPS不支持当前平台")
    endif()

else()
    set(MUMPS_AVAILABLE FALSE)
    message(STATUS "==========================================")
    message(STATUS "  MUMPS 支持已禁用（用户选项）")
    message(STATUS "==========================================")
endif()

message(STATUS "")

# ============================================================
# 包含项目子目录（在第三方库配置之后，确保编译定义可用）
# ============================================================
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(test)

# googletest配置（测试时启用）
add_subdirectory(lib/googletest)

# 创建可执行文件
add_executable(fetidp_solver src/main.cpp)

# 设置包含路径
target_include_directories(fetidp_solver
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持 shape_function_base.hpp
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen            # 新增：Eigen 库依赖
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

# 链接依赖库
target_link_libraries(fetidp_solver
    PRIVATE
        app_lib
        fe_em_lib
        fetidp_lib
        numeric_lib
        tool_lib
        spdlog_lib
        pugixml_lib
)

# ============================================================
# 平台特定的预处理定义（供C++源代码使用）
# 原因说明：
# - 这些宏在源码中通过 #ifdef 进行条件编译
# - 统一使用CMake管理，避免硬编码在各源文件中
# - 支持file_utils.cpp等文件中的Win32 API vs POSIX接口切换
# ============================================================

message(STATUS "==========================================")
message(STATUS "  平台预处理定义")
message(STATUS "==========================================")

if(IS_WINDOWS)
    # Windows平台宏（与config.h配合使用）
    add_compile_definitions(
        PLATFORM_WINDOWS=1
        PLATFORM_LINUX=0
        PLATFORM_MACOS=0
    )
    message(STATUS "  [Windows] 定义: PLATFORM_WINDOWS=1")

    # Windows特定的链接库
    # - ws2_32: Windows Socket API（网络功能）
    # - psapi: Process Status API（进程信息查询）
    # - shlwapi: Shell Lightweight Utility APIs（路径操作）
    link_libraries(ws2_32 psapi shlwapi)
    message(STATUS "  [Windows] 链接系统库: ws2_32, psapi, shlwapi")

elseif(IS_LINUX)
    # Linux平台宏
    add_compile_definitions(
        PLATFORM_WINDOWS=0
        PLATFORM_LINUX=1
        PLATFORM_MACOS=0
        _GNU_SOURCE  # 启用GNU扩展（getline、strdup等）
    )
    message(STATUS "  [Linux] 定义: PLATFORM_LINUX=1, _GNU_SOURCE")

    # Linux需要显式链接pthread（POSIX线程）
    if(TARGET Threads::Threads)
        link_libraries(Threads::Threads)
        message(STATUS "  [Linux] 链接线程库: pthread (via CMake)")
    else()
        # fallback：直接链接
        find_package(Threads REQUIRED)
        link_libraries(${CMAKE_THREAD_LIBS_INIT})
        message(STATUS "  [Linux] 链接线程库: ${CMAKE_THREAD_LIBS_INIT}")
    endif()

    # Linux常用系统库
    # - dl: 动态链接加载器（dlopen/dlsym）
    # - rt: 实时扩展（clock_gettime等）
    link_libraries(dl rt)
    message(STATUS "  [Linux] 链接系统库: dl, rt")

elseif(IS_MACOS)
    # macOS平台宏
    add_compile_definitions(
        PLATFORM_WINDOWS=0
        PLATFORM_LINUX=0
        PLATFORM_MACOS=1
    )
    message(STATUS "  [macOS] 定义: PLATFORM_MACOS=1")

    # macOS也需要pthread（基于BSD实现）
    if(TARGET Threads::Threads)
        link_libraries(Threads::Threads)
        message(STATUS "  [macOS] 链接线程库: pthread (via CMake)")
    else()
        find_package(Threads REQUIRED)
        link_libraries(${CMAKE_THREAD_LIBS_INIT})
        message(STATUS "  [macOS] 链接线程库: ${CMAKE_THREAD_LIBS_INIT}")
    endif()

    # macOS框架（Foundation用于基础系统服务）
    # 原因说明：某些高级功能可能需要Objective-C桥接
    # link_libraries("-framework Foundation")  # 暂不启用，按需添加

else()
    # 未知平台的保守定义
    add_compile_definitions(
        PLATFORM_WINDOWS=0
        PLATFORM_LINUX=1  # 默认假设为类Unix系统
        PLATFORM_MACOS=0
    )
    message(WARNING "  使用默认的Unix平台定义")
endif()

message(STATUS "")

# ============================================================
# 平台特定的安装路径配置
# 原因说明：
# - Windows: 通常安装在Program Files或用户自定义目录
# - Linux:   遵循FHS标准（/usr/local用于本地编译软件）
# - macOS:   使用/usr/local或/opt/homebrew（Apple Silicon）
# - CMAKE_INSTALL_PREFIX可通过-DCMAKE_INSTALL_PREFIX=...自定义
# ============================================================

message(STATUS "==========================================")
message(STATUS "  安装路径配置")
message(STATUS "==========================================")

if(IS_WINDOWS)
    # Windows默认安装路径（可在cmake命令行覆盖）
    if(NOT CMAKE_INSTALL_PREFIX)
        set(CMAKE_INSTALL_PREFIX "C:/Program Files/FETIDP_EM" CACHE PATH "Installation prefix" FORCE)
    endif()
    message(STATUS "  [Windows] 安装前缀: ${CMAKE_INSTALL_PREFIX}")

    # Windows下安装规则（保持Visual Studio惯例）
    install(TARGETS fetidp_solver
        RUNTIME DESTINATION bin
    )

    # 安装依赖的DLL（如果使用动态链接）
    # 注意：静态链接时不需要此步骤

elseif(IS_LINUX)
    # Linux遵循FHS (Filesystem Hierarchy Standard)
    if(NOT CMAKE_INSTALL_PREFIX)
        set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix" FORCE)
    endif()
    message(STATUS "  [Linux] 安装前缀: ${CMAKE_INSTALL_PREFIX}")

    # Linux标准目录结构
    install(TARGETS fetidp_solver
        RUNTIME DESTINATION bin
    )

    # 可选：安装头文件和库文件到系统目录
    # install(DIRECTORY include/ DESTINATION include)
    # install(TARGETS tool_lib numeric_lib ARCHIVE DESTINATION lib)

elseif(IS_MACOS)
    # macOS安装路径（区分Intel和Apple Silicon）
    if(NOT CMAKE_INSTALL_PREFIX)
        if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
            # Apple Silicon Mac: 使用Homebrew默认路径
            set(CMAKE_INSTALL_PREFIX "/opt/homebrew" CACHE PATH "Installation prefix" FORCE)
        else()
            # Intel Mac: 使用传统Unix路径
            set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix" FORCE)
        endif()
    endif()
    message(STATUS "  [macOS] 安装前缀: ${CMAKE_INSTALL_PREFIX} (架构: ${CMAKE_SYSTEM_PROCESSOR})")

    # macOS安装规则
    install(TARGETS fetidp_solver
        RUNTIME DESTINATION bin
    )

else()
    # 其他平台：使用通用默认值
    if(NOT CMAKE_INSTALL_PREFIX)
        set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation prefix" FORCE)
    endif()
    message(STATUS "  [Generic] 安装前缀: ${CMAKE_INSTALL_PREFIX}")
endif()

message(STATUS "")

# ============================================================
# 跨平台编译配置总结报告
# 输出所有关键配置信息，便于调试和问题排查
# ============================================================

message(STATUS "╔══════════════════════════════════════════════╗")
message(STATUS "║     FETI-DP 跨平台编译配置完成               ║")
message(STATUS "╠══════════════════════════════════════════════╣")
message(STATUS "║  项目版本: ${PROJECT_VERSION}                        ║")
message(STATUS "║  操作系统: ${CMAKE_SYSTEM_NAME}                    ║")
message(STATUS "║  编译器:   ${COMPILER_FAMILY} (${CMAKE_CXX_COMPILER_VERSION})           ║")
message(STATUS "║  架构:     ${CMAKE_SYSTEM_PROCESSOR}                       ║")
message(STATUS "║  构建类型: ${CMAKE_BUILD_TYPE}                      ║")
message(STATUS "╠══════════════════════════════════════════════╣")
message(STATUS "║  功能开关:                                   ║")
message(STATUS "║    MPI:      ${USE_MPI}                            ║")
message(STATUS "║    OpenMP:   ${USE_OPENMP}                           ║")
message(STATUS "║    SuperLU:  ${USE_SUPERLU_MT}                         ║")
message(STATUS "║    MUMPS:    ${USE_MUMPS}                          ║")
message(STATUS "╠══════════════════════════════════════════════╣")
message(STATUS "║  输出目录:                                   ║")
message(STATUS "║    可执行:  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}  ║")
message(STATUS "║    库文件:  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}  ║")
message(STATUS "║    安装至:  ${CMAKE_INSTALL_PREFIX}            ║")
message(STATUS "╚══════════════════════════════════════════════╝")
message(STATUS "")

# 输出下一步操作提示
if(IS_WINDOWS)
    message(STATUS "构建命令示例:")
    message(STATUS "  mkdir build && cd build")
    message(STATUS "  cmake .. -DUSE_MPI=ON -DUSE_MUMPS=ON")
    message(STATUS "  cmake --build . --config Release")
    message(STATUS "")
elseif(IS_LINUX OR IS_MACOS)
    message(STATUS "构建命令示例:")
    message(STATUS "  mkdir build && cd build")
    message(STATUS "  cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON")
    message(STATUS "  make -j$(nproc)")  # Linux
    # message(STATUS "  make -j$(sysctl -n hw.ncpu)")  # macOS
    message(STATUS "")
endif()
