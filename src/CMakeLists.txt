# 源文件目录配置
# 使用 file(GLOB) 自动搜集各子目录的源文件，无需手动维护列表

# ============================================================
# 通用宏：自动搜集指定目录下的 .cpp 源文件
# ============================================================
macro(collect_sources VAR_NAME SRC_DIR)
    file(GLOB_RECURSE ${VAR_NAME}
        CONFIGURE_DEPENDS
        "${SRC_DIR}/*.cpp"
    )
endmacro()

# ============================================================
# 电磁物理层模块
# ============================================================
collect_sources(fe_em_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fe_em")
add_library(fe_em_lib ${fe_em_lib_SOURCES})

target_include_directories(fe_em_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持ElemType枚举
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen                # 新增：Eigen库（shape_function_base.hpp依赖）
)

if(USE_MPI)
    target_link_libraries(fe_em_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

target_link_libraries(fe_em_lib PRIVATE tool_lib spdlog_lib)

message(STATUS "  fe_em_lib  : ${fe_em_lib_SOURCES}")

# ============================================================
# FETI-DP核心层模块
# ============================================================
collect_sources(fetidp_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fetidp")
add_library(fetidp_lib ${fetidp_lib_SOURCES})

target_include_directories(fetidp_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fetidp
)

target_link_libraries(fetidp_lib
    PRIVATE
        include_lib
        fe_em_lib
)

if(USE_MPI)
    target_link_libraries(fetidp_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  fetidp_lib: ${fetidp_lib_SOURCES}")

# ============================================================
# 核心数值层模块
# ============================================================
collect_sources(numeric_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/numeric")

# MUMPS 源文件从 numeric_lib 中排除，放入独立的 mumps_lib
# SuperLU 源文件从 numeric_lib 中排除，放入独立的 superlu_lib
list(FILTER numeric_lib_SOURCES EXCLUDE REGEX "em_direct_solvers_mumps\\.cpp$")
list(FILTER numeric_lib_SOURCES EXCLUDE REGEX "em_direct_solvers_mumps_complex\\.cpp$")
list(FILTER numeric_lib_SOURCES EXCLUDE REGEX "em_direct_solvers_superlu\\.cpp$")

# 新架构的后端源文件分离
# MUMPS 后端放入 mumps_lib（隔离 Intel Fortran 运行时依赖）
list(FILTER numeric_lib_SOURCES EXCLUDE REGEX "mumps_solver_backend\\.cpp$")
# SuperLU 后端放入 superlu_lib
list(FILTER numeric_lib_SOURCES EXCLUDE REGEX "superlu_solver_backend\\.cpp$")

add_library(numeric_lib ${numeric_lib_SOURCES})

target_include_directories(numeric_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(numeric_lib
    PRIVATE
        include_lib
        spdlog_lib
        tool_lib
)

# SuperLU_MT 独立编译单元
if(SUPERLU_MT_AVAILABLE AND TARGET ${SUPERLU_MT_TARGET})
    add_library(superlu_lib STATIC
        "${CMAKE_CURRENT_SOURCE_DIR}/numeric/em_direct_solvers_superlu.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/numeric/superlu_solver_backend.cpp"  # 新架构：统一SuperLU后端
    )

    target_include_directories(superlu_lib
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/include/numeric
            ${CMAKE_SOURCE_DIR}/lib/eigen
            ${CMAKE_SOURCE_DIR}/lib/spdlog/include
    )

    target_compile_definitions(superlu_lib PRIVATE HAVE_SUPERLU)

    target_link_libraries(superlu_lib
        PUBLIC
            ${SUPERLU_MT_TARGET}
        PRIVATE
            include_lib
            spdlog_lib
    )

    # numeric_lib 需要 HAVE_SUPERLU 宏来编译 SuperLU 分支的 switch/case 代码
    # 同时需要链接 superlu_lib 来解析 SuperLU 符号
    target_compile_definitions(numeric_lib PRIVATE HAVE_SUPERLU)
    target_link_libraries(numeric_lib PUBLIC superlu_lib)

    message(STATUS "  superlu_lib: em_direct_solvers_superlu.cpp (HAVE_SUPERLU defined)")
    message(STATUS "  numeric_lib: HAVE_SUPERLU defined (SuperLU symbols resolved via superlu_lib)")
endif()

# MUMPS 独立编译单元（隔离 Intel Fortran 运行时依赖）
if(MUMPS_AVAILABLE)
    add_library(mumps_lib STATIC
        "${CMAKE_CURRENT_SOURCE_DIR}/numeric/em_direct_solvers_mumps.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/numeric/em_direct_solvers_mumps_complex.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/numeric/mumps_solver_backend.cpp"  # 新架构：统一MUMPS后端
    )

    target_include_directories(mumps_lib
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/include/numeric
            ${CMAKE_SOURCE_DIR}/lib/eigen
            ${CMAKE_SOURCE_DIR}/lib/spdlog/include
            ${MUMPS_ROOT}/include
            ${MUMPS_ROOT}/libseq
    )

    target_compile_definitions(mumps_lib PRIVATE HAVE_MUMPS)

    target_link_libraries(mumps_lib
        PUBLIC
            mumps pord mpiseq openblas
        PRIVATE
            include_lib
            spdlog_lib
    )

    # /nofor_main 仅对 Intel xilink.exe 有效，MSVC link.exe 不识别
    # 改用 MAIN__ stub 函数解决 Fortran 运行时 main 冲突

    if(DEFINED INTEL_FORTRAN_LIB_DIR AND EXISTS ${INTEL_FORTRAN_LIB_DIR})
        target_link_directories(mumps_lib PUBLIC ${INTEL_FORTRAN_LIB_DIR})
    endif()

    # numeric_lib 需要 HAVE_MUMPS 宏来编译 MUMPS 分支的 switch/case 代码
    # 同时需要链接 mumps_lib 来解析 MUMPS 符号
    # /nofor_main 通过 mumps_lib 的 PUBLIC link option 自动传播
    target_compile_definitions(numeric_lib PRIVATE HAVE_MUMPS)
    target_link_libraries(numeric_lib PUBLIC mumps_lib)

    message(STATUS "  mumps_lib: em_direct_solvers_mumps.cpp (HAVE_MUMPS defined)")
    message(STATUS "  numeric_lib: HAVE_MUMPS defined (MUMPS symbols resolved via mumps_lib)")
else()
    message(STATUS "  numeric_lib: MUMPS not available")
endif()

if(USE_MPI)
    target_link_libraries(numeric_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  numeric_lib: ${numeric_lib_SOURCES}")

# ============================================================
# 基础工具层模块
# ============================================================
collect_sources(tool_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tool")
add_library(tool_lib ${tool_lib_SOURCES})

target_include_directories(tool_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/fe_em          # 新增：支持EMMeshData（方案C架构）
        ${CMAKE_SOURCE_DIR}/include/numeric         # 新增：支持ElemType枚举（消除字符串）
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
        ${CMAKE_SOURCE_DIR}/lib/eigen                # 新增：Eigen库（shape_function_base.hpp依赖）
)

target_link_libraries(tool_lib
    PRIVATE
        spdlog_lib
        pugixml_lib
)

if(USE_MPI)
    target_link_libraries(tool_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  tool_lib   : ${tool_lib_SOURCES}")

# ============================================================
# 应用层模块
# ============================================================
collect_sources(app_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/app")
add_library(app_lib ${app_lib_SOURCES})

target_include_directories(app_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持 shape_function_base.hpp
        ${CMAKE_SOURCE_DIR}/lib/eigen            # 新增：Eigen 库依赖
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

target_link_libraries(app_lib
    PRIVATE
        include_lib
        fe_em_lib
        fetidp_lib
        numeric_lib
        tool_lib
)

if(MUMPS_AVAILABLE AND TARGET mumps_lib)
    target_link_libraries(app_lib PRIVATE mumps_lib)
endif()

if(USE_MPI)
    target_link_libraries(app_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  app_lib    : ${app_lib_SOURCES}")
