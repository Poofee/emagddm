# 源文件目录配置
# 负责管理所有模块的源文件编译

# 包含头文件目录
# 注意：include目录是头文件，不需要编译，直接包含路径即可

# 电磁物理层模块
add_library(fe_em_lib
    fe_em/fe_em_mesh.cpp
)

target_include_directories(fe_em_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
)

# FETI-DP核心层模块
add_library(fetidp_lib
    fetidp/fetidp_context.cpp
)

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

# 核心数值层模块
add_library(numeric_lib
    numeric/sparse_matrix.cpp
)

target_include_directories(numeric_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
)

target_link_libraries(numeric_lib
    PRIVATE
        include_lib
)

# 基础工具层模块
add_library(tool_lib
    tool/log_interface.cpp
    tool/spdlog_adapter.cpp
    tool/logger.cpp
    tool/xml_interface.cpp
    tool/em_enums.cpp
    tool/em_exception.cpp
    tool/id_generator.cpp
    tool/file_utils.cpp
)

target_include_directories(tool_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
)

target_link_libraries(tool_lib
    PRIVATE
        spdlog_lib
        pugixml_lib
)

target_include_directories(tool_lib PRIVATE
    ${CMAKE_SOURCE_DIR}/lib/spdlog/include
    ${CMAKE_SOURCE_DIR}/lib/pugixml/src
)

# 应用层模块
add_library(app_lib
    app/solver_app.cpp
)

target_include_directories(app_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
)

target_link_libraries(app_lib
    PRIVATE
        include_lib
        fe_em_lib
        fetidp_lib
        numeric_lib
        tool_lib
)