cmake_minimum_required(VERSION 3.13)
project(CBLAS C)

# 仅编译双精度BLAS函数（项目仅需双精度求解）
file(GLOB CBLAS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/d*.c")

add_library(blas STATIC ${CBLAS_SOURCES})

target_include_directories(blas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../SRC
)
