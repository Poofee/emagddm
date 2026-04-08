# FindSuperLUDIST.cmake
# 查找 SuperLU_DIST 库的 CMake 模块
# 
# 此模块尝试查找 SuperLU_DIST 库，支持两种模式：
# 1. 使用项目内建的 lib/superlu_dist 源码编译
# 2. 查找系统已安装的 SuperLU_DIST
#
# 输出变量:
#   SUPERLUDIST_FOUND        - 是否找到 SuperLU_DIST
#   SUPERLUDIST_INCLUDE_DIRS - 包含目录路径
#   SUPERLUDIST_LIBRARIES    - 需要链接的库
#   SUPERLUDIST_VERSION      - 版本号
#
# 导入目标:
#   SuperLU_DIST::superlu_dist - 可链接的导入目标

cmake_minimum_required(VERSION 3.18)

# 默认使用项目内建的 superlu_dist 源码
set(SUPERLUDIST_USE_BUILTIN TRUE CACHE BOOL "Use built-in superlu_dist source")

if(SUPERLUDIST_USE_BUILTIN)
    # 检查内建源码是否存在
    if(EXISTS "${CMAKE_SOURCE_DIR}/lib/superlu_dist/CMakeLists.txt")
        message(STATUS "Found built-in SuperLU_DIST source at lib/superlu_dist")
        
        # 配置 SuperLU_DIST 编译选项
        set(enable_doc OFF CACHE BOOL "Disable documentation")
        set(enable_tests OFF CACHE BOOL "Disable SuperLU_DIST tests")
        set(enable_examples OFF CACHE BOOL "Disable examples")
        set(enable_python OFF CACHE BOOL "Disable Python")
        set(XSDK_ENABLE_Fortran OFF CACHE BOOL "Disable Fortran")
        
        # 精度选项（仅启用双精度）
        set(enable_double ON CACHE BOOL "Enable double precision")
        set(enable_single OFF CACHE BOOL "Disable single precision")
        set(enable_complex16 OFF CACHE BOOL "Disable complex16")
        
        # 依赖选项
        set(TPL_ENABLE_INTERNAL_BLASLIB ON CACHE BOOL "Build internal CBLAS")
        set(TPL_ENABLE_PARMETISLIB OFF CACHE BOOL "Disable ParMETIS (not needed for shared-memory)")
        set(TPL_ENABLE_COLAMDLIB OFF CACHE BOOL "Use internal COLAMD")
        set(TPL_ENABLE_LAPACKLIB OFF CACHE BOOL "Not needed")
        set(TPL_ENABLE_CUDALIB OFF CACHE BOOL "Disable CUDA")
        set(TPL_ENABLE_HIPLIB OFF CACHE BOOL "Disable HIP")
        set(TPL_ENABLE_MAGMALIB OFF CACHE BOOL "Disable MAGMA")
        
        # 添加 SuperLU_DIST 子目录
        add_subdirectory(${CMAKE_SOURCE_DIR}/lib/superlu_dist ${CMAKE_BINARY_DIR}/superlu_dist)
        
        # 设置输出变量
        set(SUPERLUDIST_FOUND TRUE)
        set(SUPERLUDIST_INCLUDE_DIRS 
            ${CMAKE_SOURCE_DIR}/lib/superlu_dist/SRC/include
            ${CMAKE_SOURCE_DIR}/lib/superlu_dist/CBLAS
        )
        
        # 库目标名称（根据配置确定）
        set(SUPERLUDIST_LIBRARIES superlu_dist)
        set(SUPERLUDIST_VERSION "9.2.1")
        
        message(STATUS "SuperLU_DIST configured: ${SUPERLUDIST_VERSION} (double precision only)")
    else()
        message(WARNING "Built-in SuperLU_DIST source not found at lib/superlu_dist")
        set(SUPERLUDIST_FOUND FALSE)
    endif()
else()
    # 尝试查找系统安装的 SuperLU_DIST
    find_path(SUPERLUDIST_INCLUDE_DIR 
        NAMES superlu_defs.h
        PATHS 
            /usr/include
            /usr/local/include
            $ENV{SUPERLUDIST_DIR}/include
        PATH_SUFFIXES superlu_dist
    )
    
    find_library(SUPERLUDIST_LIBRARY 
        NAMES superlu_dist superludist
        PATHS 
            /usr/lib
            /usr/local/lib
            $ENV{SUPERLUDIST_DIR}/lib
        PATH_SUFFIXES superlu_dist
    )
    
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SuperLUDIST
        REQUIRED_VARS SUPERLUDIST_LIBRARY SUPERLUDIST_INCLUDE_DIR
        VERSION_VAR SUPERLUDIST_VERSION
    )
    
    if(SUPERLUDIST_FOUND)
        set(SUPERLUDIST_INCLUDE_DIRS ${SUPERLUDIST_INCLUDE_DIR})
        set(SUPERLUDIST_LIBRARIES ${SUPERLUDIST_LIBRARY})
        
        # 创建导入目标
        if(NOT TARGET SuperLU_DIST::superlu_dist)
            add_library(SuperLU_DIST::superlu_dist UNKNOWN IMPORTED)
            set_target_properties(SuperLU_DIST::superlu_dist PROPERTIES
                IMPORTED_LOCATION "${SUPERLUDIST_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUDIST_INCLUDE_DIR}"
            )
        endif()
    endif()
endif()

# 标记为高级变量
mark_as_advanced(SUPERLUDIST_INCLUDE_DIRS SUPERLUDIST_LIBRARIES SUPERLUDIST_VERSION)
