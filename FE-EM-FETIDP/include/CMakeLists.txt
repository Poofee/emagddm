# 头文件目录配置
# 负责管理所有模块的头文件路径

# 设置头文件包含路径
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/fe_em
    ${CMAKE_CURRENT_SOURCE_DIR}/fetidp
    ${CMAKE_CURRENT_SOURCE_DIR}/numeric
    ${CMAKE_CURRENT_SOURCE_DIR}/tool
    ${CMAKE_CURRENT_SOURCE_DIR}/app
)

# 创建头文件库（空库，仅用于依赖管理）
add_library(include_lib INTERFACE)

target_include_directories(include_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/fe_em
    ${CMAKE_CURRENT_SOURCE_DIR}/fetidp
    ${CMAKE_CURRENT_SOURCE_DIR}/numeric
    ${CMAKE_CURRENT_SOURCE_DIR}/tool
    ${CMAKE_CURRENT_SOURCE_DIR}/app
)