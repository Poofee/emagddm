
file (GLOB sources "*.c")
add_library (superlu_mt${PLAT} ${sources})
set_target_properties (superlu_mt${PLAT} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_property(TARGET superlu_mt${PLAT} PROPERTY C_STANDARD 99)

if (HAVE_LIB_M)
  target_link_libraries (superlu_mt${PLAT} PRIVATE m)
endif ()

target_link_libraries (superlu_mt${PLAT} PRIVATE ${BLAS_LIB})
target_compile_definitions (superlu_mt${PLAT} PRIVATE Add_)

# 禁用第三方代码的警告（避免 -Werror 导致编译失败）
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(superlu_mt${PLAT} PRIVATE
        -Wno-unused-variable
        -Wno-unused-but-set-variable
        -Wno-unused-parameter
        -Wno-unused-label
        -Wno-unused-function
        -Wno-maybe-uninitialized
        -Wno-comment
        -Wno-parentheses
        -Wno-unknown-pragmas
        -Wno-sequence-point
        -Wno-old-style-definition
    )
endif()

if (PLAT STREQUAL "_PTHREAD")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __PTHREAD)
  target_link_libraries (superlu_mt${PLAT} PRIVATE ${CMAKE_THREAD_LIBS_INIT})
elseif (PLAT STREQUAL "_OPENMP")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __OPENMP)
  target_link_libraries (superlu_mt${PLAT} PRIVATE OpenMP::OpenMP_C)
endif ()

if (LONGINT)
  target_compile_definitions (superlu_mt${PLAT} PUBLIC _LONGINT)
endif ()

target_include_directories (superlu_mt${PLAT} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (superlu_mt${PLAT} INTERFACE $<INSTALL_INTERFACE:include>)
install (TARGETS superlu_mt${PLAT} DESTINATION ${CMAKE_INSTALL_LIBDIR})

file (GLOB headers "*.h")
install (FILES ${headers} DESTINATION ${SUPERLUMT_INSTALL_INCLUDEDIR})
