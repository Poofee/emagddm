cmake_minimum_required(VERSION 3.20)

project(mpiseq C Fortran)

add_library(mpiseq STATIC
    elapse.c
    mpic.c
    mpi.f
)

target_include_directories(mpiseq PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set_target_properties(mpiseq PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
