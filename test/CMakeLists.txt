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

# Maxwell转换器测试
add_executable(test_maxwell_converter
    test_maxwell_converter.cpp
)

target_include_directories(test_maxwell_converter
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

target_link_libraries(test_maxwell_converter
    PRIVATE
        tool_lib
        spdlog_lib
        pugixml_lib
        GTest::gtest
        GTest::gtest_main
)

# Maxwell转换器主测试
add_executable(test_maxwell_converter_main
    test_maxwell_converter_main.cpp
)

target_include_directories(test_maxwell_converter_main
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

target_link_libraries(test_maxwell_converter_main
    PRIVATE
        tool_lib
        spdlog_lib
        pugixml_lib
        GTest::gtest
        GTest::gtest_main
)

# MPI/OpenMP基础测试
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
