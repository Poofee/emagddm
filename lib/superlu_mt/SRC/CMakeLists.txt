
# 获取所有 C 源文件
file (GLOB ALL_SOURCES "*.c")

# 当使用外部 BLAS 时，排除与 BLAS 重复的辅助函数文件
# 避免与 OpenBLAS 产生重复符号（icmax1, dlamch, lsame 等）
# 注意：只排除标准 BLAS 辅助函数，保留 SuperLU_MT 自定义函数（如 dcomplex.c 的 z_abs/z_div）
set(_BLAS_HELPER_FILES
    icmax1.c         # BLAS: 找最大绝对值索引 (整数)
    izmax1.c         # BLAS: 复数版 icmax1
    dlaran.c         # LAPACK: 随机数生成器
    dlamch.c         # BLAS: 双精度机器参数
    slamch.c         # BLAS: 单精度机器参数
    dlamc3.c         # dlamch 内部辅助
    slamc3.c         # slamch 内部辅助
    lsame.c          # BLAS: 字符串比较
    xerbla.c         # BLAS/LAPACK: 错误处理
    dlacon.c         # LAPACK: 1范数估计器 (双精度)
    slacon.c         # LAPACK: 1范数估计器 (单精度)
    clacon.c         # LAPACK: 1范数估计器 (复单精度)
    zlacon.c         # LAPACK: 1范数估计器 (复双精度)
    scsum1.c         # LAPACK: 向量1范数 (单精度复数)
    dzsum1.c         # LAPACK: 向量1范数 (双精度复数)
)

if(NOT enable_internal_blaslib AND (BLAS_FOUND OR TPL_BLAS_LIBRARIES))
    set(sources "")
    foreach(_SRC ${ALL_SOURCES})
        get_filename_component(_SRC_NAME ${_SRC} NAME)
        if(NOT _SRC_NAME IN_LIST _BLAS_HELPER_FILES)
            list(APPEND sources ${_SRC})
        endif()
    endforeach()
    list(LENGTH sources _SRC_COUNT)
    message(STATUS "  SuperLU_MT: 排除了 ${_BLAS_HELPER_FILES} 个 BLAS 辅助文件，剩余 ${_SRC_COUNT} 个源文件 (使用外部 BLAS)")
else()
    set(sources ${ALL_SOURCES})
endif()

add_library (superlu_mt${PLAT} ${sources})
set_target_properties (superlu_mt${PLAT} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_property(TARGET superlu_mt${PLAT} PROPERTY C_STANDARD 99)

# 当使用外部 BLAS 时，重命名与 OpenBLAS 冲突的符号
# _c_abs (来自 scomplex.c) → _superlu_c_abs
# _z_abs (来自 dcomplex.c) → _superlu_z_abs
# 原因：这 2 个函数同时存在于 SuperLU_MT 和 OpenBLAS 中
if(NOT enable_internal_blaslib AND (BLAS_FOUND OR TPL_BLAS_LIBRARIES))
    set_source_files_properties(dcomplex.c PROPERTIES
        COMPILE_DEFINITIONS "z_abs=superlu_z_abs"
    )
    set_source_files_properties(scomplex.c PROPERTIES
        COMPILE_DEFINITIONS "c_abs=superlu_c_abs"
    )
endif()

if (HAVE_LIB_M)
  target_link_libraries (superlu_mt${PLAT} PRIVATE m)
endif ()

target_link_libraries (superlu_mt${PLAT} PRIVATE ${BLAS_LIB})
target_compile_definitions (superlu_mt${PLAT} PRIVATE Add_)

# 禁用第三方代码的警告（避免 -Werror 导致编译失败）
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(superlu_mt${PLAT} PRIVATE
        -Wno-unused-variable
        -Wno-unused-but-set-variable
        -Wno-unused-parameter
        -Wno-unused-label
        -Wno-unused-function
        -Wno-maybe-uninitialized
        -Wno-comment
        -Wno-parentheses
        -Wno-unknown-pragmas
        -Wno-sequence-point
        -Wno-old-style-definition
        -Wno-implicit-fallthrough
    )
endif()

if (PLAT STREQUAL "_PTHREAD")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __PTHREAD)
  target_link_libraries (superlu_mt${PLAT} PRIVATE ${CMAKE_THREAD_LIBS_INIT})
elseif (PLAT STREQUAL "_OPENMP")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __OPENMP)
  target_link_libraries (superlu_mt${PLAT} PRIVATE OpenMP::OpenMP_C)
endif ()

if (LONGINT)
  target_compile_definitions (superlu_mt${PLAT} PUBLIC _LONGINT)
endif ()

target_include_directories (superlu_mt${PLAT} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (superlu_mt${PLAT} INTERFACE $<INSTALL_INTERFACE:include>)
install (TARGETS superlu_mt${PLAT} DESTINATION ${CMAKE_INSTALL_LIBDIR})

# 后处理：移除与 OpenBLAS 重复的符号（仅保留 SuperLU_MT 专有符号）
# dcomplex.c 包含 _c_abs/_z_abs（与 OpenBLAS 重复）和 _z_abs1/_z_div 等（SuperLU_MT 专有）
# 使用 objcopy --strip-symbol 移除重复符号
if(NOT enable_internal_blaslib AND (BLAS_FOUND OR TPL_BLAS_LIBRARIES))
    find_program(OBJCOPY objcopy)
    if(OBJCOPY)
        add_custom_command(TARGET superlu_mt${PLAT} POST_BUILD
            COMMAND ${OBJCOPY} --strip-symbol=_c_abs
                                --strip-symbol=_z_abs
                    $<TARGET_FILE:superlu_mt${PLAT}>
            COMMENT "Removing duplicate c_abs/z_abs symbols from SuperLU_MT"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
    endif()
endif()

file (GLOB headers "*.h")
install (FILES ${headers} DESTINATION ${SUPERLUMT_INSTALL_INCLUDEDIR})
