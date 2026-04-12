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
)

target_link_libraries(test_sparse_matrix_gtest_stage2
    PRIVATE
        numeric_lib
        GTest::gtest
)

target_link_libraries(test_sparse_matrix_gtest_stage3_simple
    PRIVATE
        numeric_lib
        GTest::gtest
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

if(USE_OPENMP AND TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(test_tool_modules PRIVATE OpenMP::OpenMP_CXX)
endif()

# Maxwell解析器测试
# OpenMP基础测试（仅在启用OpenMP时编译）
if(USE_OPENMP)
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
endif()

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
# 注意：此测试使用旧的求解器 API（GeneralDirectSolver），需要迁移到新的 UnifiedDirectSolver 架构。暂时禁用。
# add_executable(test_linear_solver_mixed_av
#     test_linear_solver_mixed_av.cpp
# )
#
# target_include_directories(test_linear_solver_mixed_av
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
#         ${CMAKE_SOURCE_DIR}/include/numeric
#         ${CMAKE_SOURCE_DIR}/include/tool
#         ${CMAKE_SOURCE_DIR}/lib/spdlog/include
#         ${CMAKE_SOURCE_DIR}/lib/eigen
# )
#
# target_link_libraries(test_linear_solver_mixed_av
#     PRIVATE
#         numeric_lib
#         tool_lib
#         spdlog_lib
# )
#
# add_test(NAME test_linear_solver_mixed_av COMMAND test_linear_solver_mixed_av)

# ==================== 线性求解器棱边单元测试（PYRAMID5_EDGE）====================
# 注意：此测试使用旧的求解器 API（SymmetricIndefiniteDirectSolver/GeneralDirectSolver），
#       需要迁移到新的 UnifiedDirectSolver 架构。暂时禁用。
# add_executable(test_linear_solver_edge
#     test_linear_solver_edge.cpp
# )
#
# target_include_directories(test_linear_solver_edge
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
#         ${CMAKE_SOURCE_DIR}/include/numeric
#         ${CMAKE_SOURCE_DIR}/include/tool
#         ${CMAKE_SOURCE_DIR}/lib/spdlog/include
#         ${CMAKE_SOURCE_DIR}/lib/eigen
# )
#
# target_link_libraries(test_linear_solver_edge
#     PRIVATE
#         numeric_lib
#         tool_lib
#         spdlog_lib
#         GTest::gtest
#         GTest::gtest_main
# )
#
# add_test(NAME test_linear_solver_edge COMMAND test_linear_solver_edge)

# ==================== 线性求解器标量场测试（PRISM15单元）====================
# 注意：此测试使用旧的求解器 API（SymmetricDirectSolver），需要迁移到新的 UnifiedDirectSolver 架构。暂时禁用。
# add_executable(test_linear_solver_scalar
#     test_linear_solver_scalar.cpp
# )
#
# target_include_directories(test_linear_solver_scalar
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
#         ${CMAKE_SOURCE_DIR}/include/numeric
#         ${CMAKE_SOURCE_DIR}/include/tool
#         ${CMAKE_SOURCE_DIR}/lib/spdlog/include
#         ${CMAKE_SOURCE_DIR}/lib/eigen
# )
#
# target_link_libraries(test_linear_solver_scalar
#     PRIVATE
#         numeric_lib
#         tool_lib
#         spdlog_lib
#         GTest::gtest
#         GTest::gtest_main
# )
#
# add_test(NAME test_linear_solver_scalar COMMAND test_linear_solver_scalar)

# ==================== MumpsContext RAII 封装全面测试 ====================
# 仅在 MUMPS 可用时编译这些测试
if(MUMPS_AVAILABLE)
    message(STATUS "  [TEST-CFG] Configuring test_mumps_context...")
    message(STATUS "  [TEST-CFG] -> Enabling HAVE_MUMPS for test_mumps_context")

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
            ${MUMPS_ROOT}/include
            ${MUMPS_ROOT}/libseq
    )

    target_link_libraries(test_mumps_context
        PRIVATE
            numeric_lib
            tool_lib
            spdlog_lib
            GTest::gtest
            mumps_lib
    )

    if(WIN32)
        target_link_libraries(test_mumps_context PRIVATE ifport)
    else()
        target_link_libraries(test_mumps_context PRIVATE gfortran)
    endif()

    target_compile_definitions(test_mumps_context PRIVATE HAVE_MUMPS)
    add_test(NAME test_mumps_context COMMAND test_mumps_context)

    # ============================================================
    # test_complex_mumps_context
    # ============================================================
    message(STATUS "  [TEST-CFG] Configuring test_complex_mumps_context...")

    add_executable(test_complex_mumps_context test_complex_mumps_context.cpp)

    target_include_directories(test_complex_mumps_context
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/include/numeric
            ${CMAKE_SOURCE_DIR}/include/tool
            ${CMAKE_SOURCE_DIR}/lib/eigen
            ${CMAKE_SOURCE_DIR}/lib/spdlog/include
            ${MUMPS_ROOT}/include
            ${MUMPS_ROOT}/libseq
    )

    target_link_libraries(test_complex_mumps_context
        PRIVATE
            numeric_lib
            tool_lib
            spdlog_lib
            GTest::gtest
            GTest::gtest_main
            mumps_lib
    )

    if(WIN32)
        target_link_libraries(test_complex_mumps_context PRIVATE ifport)
    else()
        target_link_libraries(test_complex_mumps_context PRIVATE gfortran)
    endif()

    target_compile_definitions(test_complex_mumps_context PRIVATE HAVE_MUMPS)
    add_test(NAME test_complex_mumps_context COMMAND test_complex_mumps_context)
else()
    message(STATUS "  [TEST-CFG] MUMPS NOT available - skipping MUMPS tests")
endif()

# ==================== Eigen 扩展接口测试 ====================
add_executable(test_eigen_extension
    test_eigen_extension_simple.cpp
)

target_include_directories(test_eigen_extension
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(test_eigen_extension
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
)

add_test(NAME test_eigen_extension COMMAND test_eigen_extension)

# ==================== 实数/复数统一架构线性求解器完整功能验证 ====================
# 注意：此测试使用了旧的 SolverType 枚举值（如 EIGEN_SYMMETRIC_DIRECT），
#       但新 API 使用不同的命名（如 SYMMETRIC_DIRECT）且 SolverType 是 EMSolverFactory 的内部枚举。
#       需要更新测试代码以匹配新的 API。暂时禁用。
# add_executable(test_unified_solver
#     test_unified_solver.cpp
# )
# target_include_directories(test_unified_solver
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/include
#         ${CMAKE_SOURCE_DIR}/include/numeric
#         ${CMAKE_SOURCE_DIR}/include/tool
#         ${CMAKE_SOURCE_DIR}/lib/eigen
#         ${CMAKE_SOURCE_DIR}/lib/spdlog/include
# )
#
# target_link_libraries(test_unified_solver
#     PRIVATE
#         numeric_lib
#         tool_lib
#         spdlog_lib
#         GTest::gtest
#         GTest::gtest_main
# )
#
# add_test(NAME test_unified_solver COMMAND test_unified_solver)

# ==================== CRTP 重构后 Eigen 求解器后端单元测试 ====================
add_executable(test_eigen_solver_backends_crtp
    test_eigen_solver_backends_crtp.cpp
)

target_include_directories(test_eigen_solver_backends_crtp
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(test_eigen_solver_backends_crtp
    PRIVATE
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
)

add_test(NAME test_eigen_solver_backends_crtp COMMAND test_eigen_solver_backends_crtp)

# ==================== 静电场求解器完整测试套件 ====================
add_executable(test_electrostatic_solver
    test_electrostatic_solver.cpp
)

target_include_directories(test_electrostatic_solver
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/solver
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(test_electrostatic_solver
    PRIVATE
        solver_lib
        fe_em_lib
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
)

add_test(NAME test_electrostatic_solver COMMAND test_electrostatic_solver)

# ==================== 直流电流场求解器完整测试套件 ====================
add_executable(test_dc_conduction_solver
    test_dc_conduction_solver.cpp
)

target_include_directories(test_dc_conduction_solver
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/solver
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/eigen
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
)

target_link_libraries(test_dc_conduction_solver
    PRIVATE
        solver_lib
        fe_em_lib
        numeric_lib
        tool_lib
        spdlog_lib
        GTest::gtest
)

add_test(NAME test_dc_conduction_solver COMMAND test_dc_conduction_solver)

# ==================== Tree Gauge模块：图结构提取器测试 ====================
add_executable(test_mesh_graph_extractor
    test_mesh_graph_extractor.cpp
)

target_include_directories(test_mesh_graph_extractor
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_mesh_graph_extractor
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_mesh_graph_extractor COMMAND test_mesh_graph_extractor)

# ==================== Tree Gauge模块：生成树构建器测试 ====================
add_executable(test_spanning_tree
    test_spanning_tree.cpp
)

target_include_directories(test_spanning_tree
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_spanning_tree
    PRIVATE
        tool_lib
        fe_em_lib
        spdlog_lib
        GTest::gtest
        GTest::gtest_main
)

add_test(NAME test_spanning_tree COMMAND test_spanning_tree)

# ==================== Tree Gauge模块：多连通域检测器测试 ====================
add_executable(test_multi_connected_domain
    test_multi_connected_domain.cpp
)

target_include_directories(test_multi_connected_domain
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_multi_connected_domain
    PRIVATE
    tool_lib
    fe_em_lib
    spdlog_lib
    GTest::gtest
    GTest::gtest_main
)

add_test(NAME test_multi_connected_domain COMMAND test_multi_connected_domain)

# ==================== Tree Gauge模块：核心应用类TreeGauge测试 ====================
add_executable(test_tree_gauge
    test_tree_gauge.cpp
)

target_include_directories(test_tree_gauge
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_tree_gauge
    PRIVATE
    tool_lib
    fe_em_lib
    spdlog_lib
    GTest::gtest
    GTest::gtest_main
)

add_test(NAME test_tree_gauge COMMAND test_tree_gauge)

# ==================== Tree Gauge模块：与EMDOFManager集成测试 ====================
add_executable(test_tree_gauge_integration
    test_tree_gauge_integration.cpp
)

target_include_directories(test_tree_gauge_integration
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/include/numeric
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/eigen
)

target_link_libraries(test_tree_gauge_integration
    PRIVATE
    tool_lib
    fe_em_lib
    spdlog_lib
    GTest::gtest
    GTest::gtest_main
)

add_test(NAME test_tree_gauge_integration COMMAND test_tree_gauge_integration)

# 为所有链接了 fe_em_lib 或 numeric_lib 的测试目标添加 OpenMP 链接
# （因为 Eigen 库在启用 OpenMP 时会使用 OpenMP 并行功能）
if(USE_OPENMP AND TARGET OpenMP::OpenMP_CXX)
    foreach(test_target IN ITEMS
        test_2d_mesh test_3d_mesh
        test_pyramid5_edge test_prism15_scalar test_prism6_mixed_av
        test_assembly test_shape_function test_em_element_integrator
        test_linear_solver_mixed_av test_linear_solver_edge test_linear_solver_scalar
        test_sparse_matrix_stage1 test_sparse_matrix_stage2
        test_sparse_matrix_gtest_stage1 test_sparse_matrix_gtest_stage2 test_sparse_matrix_gtest_stage3_simple
        test_unified_solver test_eigen_extension
        test_mumps_context test_complex_mumps_context
        test_omp_basic
        test_eigen_solver_backends_crtp  # 新增：CRTP 重构测试
        test_electrostatic_solver         # 新增：静电场求解器测试
        test_dc_conduction_solver          # 新增：直流电流场求解器测试
        test_mesh_graph_extractor          # 新增：Tree Gauge图结构提取器测试
        test_spanning_tree                 # 新增：Tree Gauge生成树构建器测试
        test_multi_connected_domain        # 新增：Tree Gauge多连通域检测器测试
        test_tree_gauge                    # 新增：Tree Gauge核心应用类测试
        test_tree_gauge_integration         # 新增：Tree Gauge与DOF管理器集成测试
    )
        if(TARGET ${test_target})
            target_link_libraries(${test_target} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endforeach()
endif()
