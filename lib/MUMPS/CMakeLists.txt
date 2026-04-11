cmake_minimum_required(VERSION 3.20)

project(MUMPS LANGUAGES C Fortran)

message(STATUS "========== 配置 MUMPS 5.8.2 (双精度顺序模式) ==========")

# ==================== 选项 ====================
option(MUMPS_USE_SEQMPI "使用顺序 MPI 桩（无需真实 MPI）" ON)
option(MUMPS_BUILD_PORD "构建 PORD 排序库" ON)

# ==================== Fortran 编译器检测 ====================
if(NOT CMAKE_Fortran_COMPILER AND WIN32)
    set(INTEL_FORTRAN_PATHS
        "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/bin/intel64/ifort.exe"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/2023.2.0/windows/bin/intel64/ifort.exe"
    )
    foreach(_ifort_path ${INTEL_FORTRAN_PATHS})
        if(EXISTS ${_ifort_path})
            set(CMAKE_Fortran_COMPILER ${_ifort_path} CACHE FILEPATH "Intel Fortran Compiler")
            message(STATUS "找到 Intel Fortran: ${_ifort_path}")
            break()
        endif()
    endforeach()
endif()

enable_language(Fortran)

# ==================== 编译选项 ====================
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /fpp /4Yportlib /assume:buffered_io /names:lowercase /assume:underscore")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} /O2 /MD")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od /MDd /debug:full")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch -fno-second-underscore")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")
endif()

# ==================== C 编译选项 ====================
if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /DADD_")
    if(MUMPS_USE_SEQMPI)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /DMUMPS_SEQ_MPI")
    endif()
else()
    add_compile_definitions(Add_)
    if(MUMPS_USE_SEQMPI)
        add_compile_definitions(MUMPS_SEQ_MPI)
    endif()
endif()

message(STATUS "Fortran 编译器: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

# ==================== OpenBLAS ====================
set(OPENBLAS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../OpenBLAS")
find_package(OpenBLAS QUIET PATHS ${OPENBLAS_ROOT} NO_DEFAULT_PATH)

if(OpenBLAS_FOUND)
    message(STATUS "MUMPS: 使用 OpenBLAS ${OpenBLAS_VERSION}")
    set(MUMPS_BLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
    set(MUMPS_BLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
else()
    # 回退：手动指定
    if(EXISTS ${OPENBLAS_ROOT}/lib/libopenblas.a)
        set(MUMPS_BLAS_LIBRARIES ${OPENBLAS_ROOT}/lib/libopenblas.a)
        set(MUMPS_BLAS_INCLUDE_DIRS ${OPENBLAS_ROOT}/include)
        message(STATUS "MUMPS: 手动配置 OpenBLAS")
    else()
        message(FATAL_ERROR "MUMPS: 未找到 OpenBLAS (${OPENBLAS_ROOT})")
    endif()
endif()

# ==================== 子目录 ====================
if(MUMPS_BUILD_PORD)
    add_subdirectory(PORD)
endif()

if(MUMPS_USE_SEQMPI)
    add_subdirectory(libseq)
endif()

add_subdirectory(src)

message(STATUS "========== MUMPS 配置完成 ==========")
