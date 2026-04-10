cmake_minimum_required(VERSION 3.20)

project(pord C)

add_library(pord STATIC
    lib/bucket.c
    lib/ddbisect.c
    lib/ddcreate.c
    lib/gbipart.c
    lib/gbisect.c
    lib/gelim.c
    lib/graph.c
    lib/interface.c
    lib/minpriority.c
    lib/multisector.c
    lib/nestdiss.c
    lib/sort.c
    lib/symbfac.c
    lib/tree.c
)

target_include_directories(pord PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties(pord PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
