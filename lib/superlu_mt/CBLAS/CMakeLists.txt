cmake_minimum_required(VERSION 3.13)
project(CBLAS C)

# 仅编译双精度BLAS函数（项目仅需双精度求解）
file(GLOB CBLAS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/d*.c")

add_library(blas STATIC ${CBLAS_SOURCES})

target_include_directories(blas PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../SRC
)

# 禁用第三方代码的警告（避免 -Werror 导致编译失败）
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(blas PRIVATE
        -Wno-unused-variable
        -Wno-unused-but-set-variable
        -Wno-parentheses
        -Wno-old-style-definition
    )
endif()
