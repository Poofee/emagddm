# 源文件目录配置
# 使用 file(GLOB) 自动搜集各子目录的源文件，无需手动维护列表

# ============================================================
# 通用宏：自动搜集指定目录下的 .cpp 源文件
# ============================================================
macro(collect_sources VAR_NAME SRC_DIR)
    file(GLOB_RECURSE ${VAR_NAME}
        CONFIGURE_DEPENDS
        "${SRC_DIR}/*.cpp"
    )
endmacro()

# ============================================================
# 电磁物理层模块
# ============================================================
collect_sources(fe_em_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fe_em")
add_library(fe_em_lib ${fe_em_lib_SOURCES})

target_include_directories(fe_em_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
)

if(USE_MPI)
    target_link_libraries(fe_em_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  fe_em_lib  : ${fe_em_lib_SOURCES}")

# ============================================================
# FETI-DP核心层模块
# ============================================================
collect_sources(fetidp_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fetidp")
add_library(fetidp_lib ${fetidp_lib_SOURCES})

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

if(USE_MPI)
    target_link_libraries(fetidp_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  fetidp_lib: ${fetidp_lib_SOURCES}")

# ============================================================
# 核心数值层模块
# ============================================================
collect_sources(numeric_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/numeric")
add_library(numeric_lib ${numeric_lib_SOURCES})

target_include_directories(numeric_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/numeric
)

target_link_libraries(numeric_lib
    PRIVATE
        include_lib
)

if(USE_MPI)
    target_link_libraries(numeric_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  numeric_lib: ${numeric_lib_SOURCES}")

# ============================================================
# 基础工具层模块
# ============================================================
collect_sources(tool_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tool")
add_library(tool_lib ${tool_lib_SOURCES})

target_include_directories(tool_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/pugixml/src
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

target_link_libraries(tool_lib
    PRIVATE
        spdlog_lib
        pugixml_lib
)

if(USE_MPI)
    target_link_libraries(tool_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  tool_lib   : ${tool_lib_SOURCES}")

# ============================================================
# 应用层模块
# ============================================================
collect_sources(app_lib_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/app")
add_library(app_lib ${app_lib_SOURCES})

target_include_directories(app_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
        ${CMAKE_SOURCE_DIR}/lib/spdlog/include
        ${CMAKE_SOURCE_DIR}/lib/nlohmann
)

target_link_libraries(app_lib
    PRIVATE
        include_lib
        fe_em_lib
        fetidp_lib
        numeric_lib
        tool_lib
)

if(USE_MPI)
    target_link_libraries(app_lib PRIVATE ${MPI_CXX_LIBRARIES})
endif()

message(STATUS "  app_lib    : ${app_lib_SOURCES}")
