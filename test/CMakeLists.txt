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

target_include_directories(test_sparse_matrix_stage1
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_sparse_matrix_stage1
    PRIVATE
        numeric_lib
)

# 启用测试
enable_testing()
add_test(NAME test_logger COMMAND test_logger)
add_test(NAME test_sparse_matrix_stage1 COMMAND test_sparse_matrix_stage1)