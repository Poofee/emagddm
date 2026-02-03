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
)

target_link_libraries(test_logger
    PRIVATE
        tool_lib
        spdlog_lib
)

# 启用测试
enable_testing()
add_test(NAME test_logger COMMAND test_logger)