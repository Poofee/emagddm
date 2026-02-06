# FETI-DP电磁有限元求解器 - 根CMake配置文件
# 负责全局编译选项、第三方库依赖、子目录包含

cmake_minimum_required(VERSION 3.20)
project(FE_EM_FETIDP VERSION 1.0.0 LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI/OpenMP编译选项
option(USE_MPI "Enable MPI support" OFF)
option(USE_OPENMP "Enable OpenMP support" OFF)

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

# 包含子目录
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(test)

# 第三方库配置
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

# googletest配置（测试时启用）
add_subdirectory(lib/googletest)

# 创建可执行文件
add_executable(fetidp_solver src/main.cpp)

# 设置包含路径
target_include_directories(fetidp_solver
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
)

# 链接依赖库
target_link_libraries(fetidp_solver
    PRIVATE
        app_lib
        fe_em_lib
        fetidp_lib
        numeric_lib
        tool_lib
)