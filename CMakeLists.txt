# FETI-DP电磁有限元求解器 - 根CMake配置文件
# 负责全局编译选项、第三方库依赖、子目录包含

cmake_minimum_required(VERSION 3.20)
project(FE_EM_FETIDP VERSION 1.0.0 LANGUAGES CXX C)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI/OpenMP 编译选项
option(USE_MPI "Enable MPI support" OFF)
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_SUPERLU_MT "Enable SuperLU_MT support" ON)

# 全局编译选项
if(MSVC)
    # MSVC编译器选项
    add_compile_options(/W4 /utf-8)
    # 禁用特定警告：C4100（未引用参数），C4101（未引用局部变量）
    add_compile_options(/wd4100 /wd4101)
    
    # OpenMP支持（MSVC）
    if(USE_OPENMP)
        add_compile_options(/openmp)
        message(STATUS "OpenMP enabled for MSVC")
    endif()
else()
    # GCC/Clang编译器选项
    add_compile_options(-Wall -Wextra -Werror -pedantic)
    
    # OpenMP支持（GCC/Clang）
    if(USE_OPENMP)
        find_package(OpenMP REQUIRED)
        if(OpenMP_FOUND)
            add_compile_options(${OpenMP_CXX_FLAGS})
            message(STATUS "OpenMP enabled: ${OpenMP_CXX_VERSION}")
        endif()
    endif()
endif()

# MPI支持
if(USE_MPI)
    # 使用本地MPI库路径
    set(MPI_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/Microsoft SDKs/MPI/Include")
    if(EXISTS "${MPI_INCLUDE_PATH}/mpi.h")
        message(STATUS "Found local MPI library at: ${MPI_INCLUDE_PATH}")
        add_definitions(-DUSE_MPI=1)
        include_directories(${MPI_INCLUDE_PATH})
        # 设置MPI库路径（Windows下可能需要链接msmpi.lib）
        if(WIN32)
            set(MPI_LIBRARIES "msmpi.lib")
            link_libraries(${MPI_LIBRARIES})
        endif()
    else()
        message(FATAL_ERROR "MPI header not found at: ${MPI_INCLUDE_PATH}")
    endif()
else()
    add_definitions(-DUSE_MPI=0)
    add_definitions(-DMPI_ROOT="")
endif()

# OpenMP编译定义
if(USE_OPENMP)
    add_definitions(-DUSE_OPENMP=1)
else()
    add_definitions(-DUSE_OPENMP=0)
endif()

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

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

# spdlog库配置（需要编译部分源文件）
add_library(spdlog_lib
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/spdlog.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/cfg.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/color_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/file_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/stdout_sinks.cpp
    ${CMAKE_SOURCE_DIR}/lib/spdlog/src/async.cpp
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

# 设置编译选项
if(MSVC)
    target_compile_options(spdlog_lib PRIVATE /W4)
else()
    target_compile_options(spdlog_lib PRIVATE -Wall -Wextra)
endif()

# ============================================================
# SuperLU_MT 库配置（多线程稀疏 LU 分解直接求解器）
# 使用 OpenMP 模式，避免 Windows pthread 兼容性问题
# ============================================================
if(USE_SUPERLU_MT)
    message(STATUS "Configuring SuperLU_MT with OpenMP...")
    
    # SuperLU_MT 编译选项设置
    set(PLAT "_OPENMP" CACHE STRING "SuperLU_MT threading flavor (OpenMP)" FORCE)
    set(enable_internal_blaslib ON CACHE BOOL "Build internal CBLAS for SuperLU_MT" FORCE)
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
# MUMPS 配置（可选，需预编译 MUMPS + OpenBLAS）
# 使用 build_windows.ps1 预编译的静态库
# ============================================================
option(USE_MUMPS "启用 MUMPS 并行直接求解器后端" OFF)

if(USE_MUMPS)
    message(STATUS "Configuring MUMPS 5.8.2 (double precision, sequential MPI)...")

    set(MUMPS_ROOT "${CMAKE_SOURCE_DIR}/lib/MUMPS")
    set(OPENBLAS_ROOT "${CMAKE_SOURCE_DIR}/lib/OpenBLAS")

    # 检查预编译库是否存在
    set(MUMPS_LIBS
        "${MUMPS_ROOT}/lib/mumps.lib"
        "${MUMPS_ROOT}/lib/pord.lib"
        "${MUMPS_ROOT}/lib/mpiseq.lib"
    )
    set(MUMPS_LIBS_FOUND TRUE)
    foreach(_lib ${MUMPS_LIBS})
        if(NOT EXISTS ${_lib})
            message(WARNING "  MUMPS pre-built library not found: ${_lib}")
            set(MUMPS_LIBS_FOUND FALSE)
        endif()
    endforeach()

    # 检查 OpenBLAS
    if(EXISTS "${OPENBLAS_ROOT}/lib/libopenblas.dll.a")
        set(OPENBLAS_LIB "${OPENBLAS_ROOT}/lib/libopenblas.dll.a")
        message(STATUS "  OpenBLAS found: ${OPENBLAS_LIB}")
    elseif(EXISTS "${OPENBLAS_ROOT}/lib/libopenblas.a")
        set(OPENBLAS_LIB "${OPENBLAS_ROOT}/lib/libopenblas.a")
        message(STATUS "  OpenBLAS found: ${OPENBLAS_LIB}")
    else()
        message(WARNING "  OpenBLAS library not found in ${OPENBLAS_ROOT}/lib")
        set(MUMPS_LIBS_FOUND FALSE)
    endif()

    if(MUMPS_LIBS_FOUND)
        # 导入预编译的 MUMPS 静态库
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

        set(INTEL_FORTRAN_LIB_DIR "C:/Program Files (x86)/Intel/oneAPI/compiler/2023.2.0/windows/compiler/lib/intel64")
        if(EXISTS ${INTEL_FORTRAN_LIB_DIR})
            message(STATUS "  Intel Fortran runtime lib dir: ${INTEL_FORTRAN_LIB_DIR}")
        else()
            message(WARNING "  Intel Fortran runtime lib dir not found: ${INTEL_FORTRAN_LIB_DIR}")
        endif()

        set(MUMPS_AVAILABLE TRUE)
        include_directories(${MUMPS_ROOT}/include)
        include_directories(${MUMPS_ROOT}/libseq)
        message(STATUS "  MUMPS configured successfully (pre-built libraries)")
    else()
        set(MUMPS_AVAILABLE FALSE)
        message(WARNING "  MUMPS pre-built libraries not found. Run build_windows.ps1 first.")
    endif()
else()
    set(MUMPS_AVAILABLE FALSE)
    message(STATUS "MUMPS support disabled by user option")
endif()

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