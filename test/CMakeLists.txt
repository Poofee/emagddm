# 测试目录配置
# 负责管理单元测试和集成测试

# 包含头文件和源文件目录
# 注意：include目录是头文件，不需要编译，直接包含路径即可

# 日志功能测试
add_executable(test_logger
    test_logger.cpp
)

target_include_directories(test_logger
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

target_link_libraries(test_logger
    PRIVATE
        tool_lib
        spdlog_lib
        pugixml_lib
)

# XML功能测试
add_executable(test_xml
    test_xml.cpp
)

target_include_directories(test_xml
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

target_link_libraries(test_xml
    PRIVATE
        tool_lib
        spdlog_lib
        pugixml_lib
)

# 稀疏矩阵功能测试
add_executable(test_sparse_matrix_stage1
    test_sparse_matrix_stage1.cpp
)

add_executable(test_sparse_matrix_stage2
    test_sparse_matrix_stage2.cpp
)

# 稀疏矩阵功能测试（Google Test版本）
add_executable(test_sparse_matrix_gtest_stage1
    test_sparse_matrix_gtest_stage1.cpp
)

add_executable(test_sparse_matrix_gtest_stage2
    test_sparse_matrix_gtest_stage2.cpp
)

# 暂时禁用test_sparse_matrix_gtest_stage3，因为API不兼容
# add_executable(test_sparse_matrix_gtest_stage3
#     test_sparse_matrix_gtest_stage3.cpp
# )
add_executable(test_sparse_matrix_gtest_stage3_simple
    test_sparse_matrix_gtest_stage3_simple.cpp
)

target_include_directories(test_sparse_matrix_stage1
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_include_directories(test_sparse_matrix_stage2
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_include_directories(test_sparse_matrix_gtest_stage1
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_include_directories(test_sparse_matrix_gtest_stage2
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

# 暂时禁用test_sparse_matrix_gtest_stage3配置
# target_include_directories(test_sparse_matrix_gtest_stage3
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
# )

target_include_directories(test_sparse_matrix_gtest_stage3_simple
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_sparse_matrix_stage1
    PRIVATE
        numeric_lib
)

target_link_libraries(test_sparse_matrix_stage2
    PRIVATE
        numeric_lib
)

target_link_libraries(test_sparse_matrix_gtest_stage1
    PRIVATE
        numeric_lib
        GTest::gtest
        GTest::gtest_main
)

target_link_libraries(test_sparse_matrix_gtest_stage2
    PRIVATE
        numeric_lib
        GTest::gtest
        GTest::gtest_main
)

# 暂时禁用test_sparse_matrix_gtest_stage3链接配置
# target_link_libraries(test_sparse_matrix_gtest_stage3
#     PRIVATE
#         numeric_lib
#         GTest::gtest
#         GTest::gtest_main
# )

target_link_libraries(test_sparse_matrix_gtest_stage3_simple
    PRIVATE
        numeric_lib
        GTest::gtest
        GTest::gtest_main
)

# 电磁模块工具类测试
add_executable(test_em_tools
    test_em_tools.cpp
)

target_include_directories(test_em_tools
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
)

target_link_libraries(test_em_tools
    PRIVATE
        tool_lib
        spdlog_lib
)

# 基础工具层模块测试（Google Test版本）
add_executable(test_tool_modules
    test_tool_modules.cpp
)

target_include_directories(test_tool_modules
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持 shape_function_base.hpp
        ${CMAKE_SOURCE_DIR}/lib/eigen            # 新增：Eigen 库依赖
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

target_link_libraries(test_tool_modules
    PRIVATE
        tool_lib
        spdlog_lib
        pugixml_lib
        GTest::gtest
        GTest::gtest_main
)

# Maxwell解析器测试
# OpenMP基础测试
add_executable(test_omp_basic
    test_omp_basic.cpp
)

target_include_directories(test_omp_basic
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_omp_basic
    PRIVATE
        tool_lib
)

# MPI基础测试
add_executable(test_mpi_basic
    test_mpi_basic.cpp
)

target_include_directories(test_mpi_basic
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_mpi_basic
    PRIVATE
        tool_lib
)

# Maxwell解析器测试
add_executable(test_maxwell_parser
    test_maxwell_parser.cpp
)

target_include_directories(test_maxwell_parser
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_maxwell_parser
    PRIVATE
        tool_lib
        spdlog_lib
)

# 材料数据提取测试
add_executable(test_material_extraction
    test_material_extraction.cpp
)

target_include_directories(test_material_extraction
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_material_extraction
    PRIVATE
        tool_lib
        spdlog_lib
)

# AEDT解析验证测试
add_executable(test_aedt_parsing
    test_aedt_parsing.cpp
)

target_include_directories(test_aedt_parsing
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

target_link_libraries(test_aedt_parsing
    PRIVATE
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

# 调试解析器
# add_executable(debug_parser
#     debug_parser.cpp
# )
# 
# target_include_directories(debug_parser
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
#         ${CMAKE_SOURCE_DIR}/lib/nlohmann
# )
# 
# target_link_libraries(debug_parser
#     PRIVATE
#         tool_lib
#         spdlog_lib
# )

# 启用测试
enable_testing()
add_test(NAME test_tool_modules COMMAND test_tool_modules)
add_test(NAME test_logger COMMAND test_logger)
add_test(NAME test_em_tools COMMAND test_em_tools)
add_test(NAME test_sparse_matrix_stage1 COMMAND test_sparse_matrix_stage1)
add_test(NAME test_sparse_matrix_stage2 COMMAND test_sparse_matrix_stage2)
add_test(NAME test_sparse_matrix_gtest_stage1 COMMAND test_sparse_matrix_gtest_stage1)
add_test(NAME test_sparse_matrix_gtest_stage2 COMMAND test_sparse_matrix_gtest_stage2)
add_test(NAME test_omp_basic COMMAND test_omp_basic)
add_test(NAME test_mpi_basic COMMAND test_mpi_basic)
add_test(NAME test_maxwell_parser COMMAND test_maxwell_parser)
add_test(NAME test_material_extraction COMMAND test_material_extraction)
add_test(NAME test_aedt_parsing COMMAND test_aedt_parsing)

# 图像提取测试
add_executable(test_image_extraction
    test_image_extraction.cpp
)

target_include_directories(test_image_extraction
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_image_extraction
    PRIVATE
        tool_lib
        spdlog_lib
)

add_test(NAME test_image_extraction COMMAND test_image_extraction)

# 形函数模块测试
add_executable(test_shape_function
    test_shape_function.cpp
)

target_include_directories(test_shape_function
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_shape_function
    PRIVATE
        numeric_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_shape_function COMMAND test_shape_function)

# 电磁场有限元单元矩阵计算模块测试
add_executable(test_em_element_integrator
    test_em_element_integrator.cpp
)

target_include_directories(test_em_element_integrator
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_em_element_integrator
    PRIVATE
        numeric_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_em_element_integrator COMMAND test_em_element_integrator)

# ==================== 电磁网格模块测试（方案C架构）====================

# 二维静磁标量位网格测试（TRI3单元）
add_executable(test_2d_mesh
    ../tests/test_2d_mesh.cpp
)

target_include_directories(test_2d_mesh
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持ElemType枚举
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen                # 新增：Eigen库（shape_function_base.hpp依赖）
)

target_link_libraries(test_2d_mesh
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
)

add_test(NAME test_2d_mesh COMMAND test_2d_mesh)

# 三维涡流场A-V混合格式测试（TET4单元）
add_executable(test_3d_mesh
    ../tests/test_3d_mesh.cpp
)

target_include_directories(test_3d_mesh
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric      # 新增：支持ElemType枚举
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen                # 新增：Eigen库（shape_function_base.hpp依赖）
)

target_link_libraries(test_3d_mesh
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
)

add_test(NAME test_3d_mesh COMMAND test_3d_mesh)

# 一阶金字塔纯棱边格式DOF管理测试（PYRAMID5_EDGE单元）
add_executable(test_pyramid5_edge
    test_pyramid5_edge.cpp
)

target_include_directories(test_pyramid5_edge
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_pyramid5_edge
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_pyramid5_edge COMMAND test_pyramid5_edge)

# PRISM15纯标量格式DOF管理测试
add_executable(test_prism15_scalar
    test_prism15_scalar.cpp
)

target_include_directories(test_prism15_scalar
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_prism15_scalar
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_prism15_scalar COMMAND test_prism15_scalar)

# PRISM6混合A-V格式DOF管理完整测试（GlobalEdgeIDGenerator + EMDOFManager）
add_executable(test_prism6_mixed_av
    test_prism6_mixed_av.cpp
)

target_include_directories(test_prism6_mixed_av
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_prism6_mixed_av
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_prism6_mixed_av COMMAND test_prism6_mixed_av)

# EMAssembly全局装配模块完整测试套件
add_executable(test_assembly
    test_assembly.cpp
)

target_include_directories(test_assembly
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_assembly
    PRIVATE
        numeric_lib
        fe_em_lib
        tool_lib
        spdlog_lib
)

add_test(NAME test_assembly COMMAND test_assembly)

# ==================== 线性求解器测试：A-V混合涡流场（PRISM6混合单元）====================
add_executable(test_linear_solver_mixed_av
    test_linear_solver_mixed_av.cpp
)

target_include_directories(test_linear_solver_mixed_av
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_linear_solver_mixed_av
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
)

add_test(NAME test_linear_solver_mixed_av COMMAND test_linear_solver_mixed_av)

# ==================== 线性求解器棱边单元测试（PYRAMID5_EDGE）====================

add_executable(test_linear_solver_edge
    test_linear_solver_edge.cpp
)

target_include_directories(test_linear_solver_edge
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_linear_solver_edge
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_linear_solver_edge COMMAND test_linear_solver_edge)

# ==================== 线性求解器标量场测试（PRISM15单元）====================
# Task 6: 测试用例1 - 标量场泊松方程求解（PRISM15单元）
# 验证 SymmetricDirectSolver 在 15×15 对称正定刚度矩阵上的正确性

add_executable(test_linear_solver_scalar
    test_linear_solver_scalar.cpp
)

target_include_directories(test_linear_solver_scalar
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_linear_solver_scalar
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_linear_solver_scalar COMMAND test_linear_solver_scalar)

# ==================== SuperLU_MT 集成测试 ====================
# 验证 SuperLU_MT 多线程稀疏LU求解器在项目中的正确集成
# 包括：后端可用性检测、SPD矩阵求解、多次solve一致性、与Eigen交叉验证等
add_executable(test_superlu_mt_integration
    test_superlu_mt_integration.cpp
)

target_include_directories(test_superlu_mt_integration
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_superlu_mt_integration
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

# 若SuperLU_MT可用，链接对应库并启用SuperLU后端
if(SUPERLU_MT_AVAILABLE AND TARGET superlu_lib)
    target_compile_definitions(test_superlu_mt_integration PRIVATE HAVE_SUPERLU)
endif()

add_test(NAME test_superlu_mt_integration COMMAND test_superlu_mt_integration)

# ==================== SuperLU_MT 诊断测试（用于定位崩溃）====================
add_executable(test_superlu_diagnostic
    test_superlu_diagnostic.cpp
)

target_include_directories(test_superlu_diagnostic
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_superlu_diagnostic
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

message(STATUS "  [TEST-CFG] SUPERLU_MT_AVAILABLE=${SUPERLU_MT_AVAILABLE}")
message(STATUS "  [TEST-CFG] SUPERLU_MT_TARGET=${SUPERLU_MT_TARGET}")
if(SUPERLU_MT_AVAILABLE AND TARGET superlu_lib)
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_SUPERLU for test_superlu_diagnostic")
    target_compile_definitions(test_superlu_diagnostic PRIVATE HAVE_SUPERLU)
else()
    message(STATUS "  [TEST-CFG] -> SuperLU NOT available for test_superlu_diagnostic")
endif()

add_test(NAME test_superlu_diagnostic COMMAND test_superlu_diagnostic)

# ==================== SuperLU_MT 全面功能与性能测试 ====================
add_executable(test_superlu_comprehensive
    test_superlu_comprehensive.cpp
)

target_include_directories(test_superlu_comprehensive
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_superlu_comprehensive
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

message(STATUS "  [TEST-CFG] Configuring test_superlu_comprehensive...")
if(SUPERLU_MT_AVAILABLE AND TARGET superlu_lib)
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_SUPERLU for test_superlu_comprehensive")
    target_compile_definitions(test_superlu_comprehensive PRIVATE HAVE_SUPERLU)
else()
    message(STATUS "  [TEST-CFG] -> SuperLU NOT available for test_superlu_comprehensive")
endif()

add_test(NAME test_superlu_comprehensive COMMAND test_superlu_comprehensive)

# ==================== SuperluContext RAII 封装全面测试 ====================
# 覆盖 SuperluContext 类的所有公开接口：构造/析构、initialize、factorize、solve、
# reset、移动语义、状态机转换、内存安全、不同规模矩阵
add_executable(test_superlu_context
    test_superlu_context.cpp
)

target_include_directories(test_superlu_context
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_superlu_context
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

message(STATUS "  [TEST-CFG] Configuring test_superlu_context...")
if(SUPERLU_MT_AVAILABLE AND TARGET superlu_lib)
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_SUPERLU for test_superlu_context")
    target_compile_definitions(test_superlu_context PRIVATE HAVE_SUPERLU)
else()
    message(STATUS "  [TEST-CFG] -> SuperLU NOT available for test_superlu_context")
endif()

add_test(NAME test_superlu_context COMMAND test_superlu_context)

# ==================== MumpsContext RAII 封装全面测试 ====================
# 覆盖 MumpsContext 类的所有公开接口：构造/析构、initialize、factorize、solve、
# reset、移动语义、状态机转换、内存安全、精度验证、与 SymmetricDirectSolver 集成
add_executable(test_mumps_context
    test_mumps_context.cpp
)

target_include_directories(test_mumps_context
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_mumps_context
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

message(STATUS "  [TEST-CFG] Configuring test_mumps_context...")
if(MUMPS_AVAILABLE)
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_MUMPS for test_mumps_context")
    target_link_libraries(test_mumps_context PRIVATE mumps_lib)
    target_compile_definitions(test_mumps_context PRIVATE HAVE_MUMPS)
    target_include_directories(test_mumps_context PRIVATE ${MUMPS_ROOT}/include ${MUMPS_ROOT}/libseq)
else()
    message(STATUS "  [TEST-CFG] -> MUMPS NOT available for test_mumps_context")
endif()

add_test(NAME test_mumps_context COMMAND test_mumps_context)

# ============================================================
# test_complex_mumps_context
# ============================================================
add_executable(test_complex_mumps_context test_complex_mumps_context.cpp)

target_include_directories(test_complex_mumps_context
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(test_complex_mumps_context
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

message(STATUS "  [TEST-CFG] Configuring test_complex_mumps_context...")
if(MUMPS_AVAILABLE)
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_MUMPS for test_complex_mumps_context")
    target_link_libraries(test_complex_mumps_context PRIVATE mumps_lib)
    target_compile_definitions(test_complex_mumps_context PRIVATE HAVE_MUMPS)
    target_include_directories(test_complex_mumps_context PRIVATE ${MUMPS_ROOT}/include ${MUMPS_ROOT}/libseq)
else()
    message(STATUS "  [TEST-CFG] -> MUMPS NOT available for test_complex_mumps_context")
endif()

add_test(NAME test_complex_mumps_context COMMAND test_complex_mumps_context)
