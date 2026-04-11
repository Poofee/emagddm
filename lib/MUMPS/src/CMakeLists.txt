cmake_minimum_required(VERSION 3.20)

project(mumps C Fortran)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/mumps_int_def32_h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/mumps_int_def.h
    COPYONLY
)

# ==================== MUMPS 源文件定义 ====================
# 架构：公共代码(mumps_common) + dmups专用 + zmumps专用
# 避免 libdmumps.a 和 libzmumps.a 之间的重复符号

set(MUMPS_COMMON_C_SOURCES
    mumps_addr.c
    mumps_common.c
    mumps_io.c
    mumps_io_basic.c
    mumps_io_err.c
    mumps_io_thread.c
    mumps_metis.c
    mumps_metis_int.c
    mumps_metis64.c
    mumps_numa.c
    mumps_pord.c
    mumps_register_thread.c
    mumps_thread.c
    mumps_thread_affinity.c
    mumps_flytes.c
    mumps_config_file_C.c
    mumps_save_restore_C.c
    mumps_scotch.c
    mumps_scotch_int.c
    mumps_scotch64.c
)

set(MUMPS_COMMON_F_SOURCES
    ana_AMDMF.F
    ana_blk.F
    ana_blk_m.F
    ana_omp_m.F
    ana_orderings.F
    ana_orderings_wrappers_m.F
    ana_set_ordering.F
    bcast_errors.F
    comp_tps_m.F
    double_linked_list.F
    estim_flops.F
    fac_asm_build_sort_index_m.F
    fac_asm_build_sort_index_ELT_m.F
    fac_descband_data_m.F
    fac_future_niv2_mod.F
    fac_maprow_data_m.F
    front_data_mgt_m.F
    lr_common.F
    lr_stats.F
    mumps_comm_buffer_common.F
    mumps_intr_types_common.F
    mumps_l0_omp_m.F
    mumps_load.F
    mumps_memory_mod.F
    mumps_ooc_common.F
    mumps_pivnul_mod.F
    mumps_print_defined.F
    mumps_static_mapping.F
    mumps_type2_blocking.F
    mumps_version.F
    omp_tps_common_m.F
    sol_common.F
    sol_ds_common_m.F
    sol_omp_common_m.F
    tools_common.F
    tools_common_m.F
)

set(DMUMPS_SOURCES
    dana_LDLT_preprocess.F
    dana_aux.F
    dana_aux_ELT.F
    dana_aux_par.F
    dana_dist_m.F
    dana_driver.F
    dana_lr.F
    dana_mtrans.F
    dana_reordertree.F
    darrowheads.F
    dbcast_int.F
    dend_driver.F
    dfac_asm.F
    dfac_asm_ELT.F
    dfac_asm_master_ELT_m.F
    dfac_asm_master_m.F
    dfac_b.F
    dfac_compact_factors_m.F
    dfac_determinant.F
    dfac_diag.F
    dfac_dist_arrowheads_omp.F
    dfac_distrib_ELT.F
    dfac_distrib_distentry.F
    dfac_driver.F
    dfac_front_LDLT_type1.F
    dfac_front_LDLT_type2.F
    dfac_front_LU_type1.F
    dfac_front_LU_type2.F
    dfac_front_aux.F
    dfac_front_type2_aux.F
    dfac_lastrtnelind.F
    dfac_lr.F
    dfac_mem_alloc_cb.F
    dfac_mem_compress_cb.F
    dfac_mem_dynamic.F
    dfac_mem_free_block_cb.F
    dfac_mem_stack.F
    dfac_mem_stack_aux.F
    dfac_omp_m.F
    dfac_par_m.F
    dfac_process_band.F
    dfac_process_bf.F
    dfac_process_blfac_slave.F
    dfac_process_blocfacto.F
    dfac_process_blocfacto_LDLT.F
    dfac_process_contrib_type1.F
    dfac_process_contrib_type2.F
    dfac_process_contrib_type3.F
    dfac_process_end_facto_slave.F
    dfac_process_maprow.F
    dfac_process_master2.F
    dfac_process_message.F
    dfac_process_root2slave.F
    dfac_process_root2son.F
    dfac_process_rtnelind.F
    dfac_root_parallel.F
    dfac_scalings.F
    dfac_scalings_simScaleAbs.F
    dfac_scalings_simScale_util.F
    dfac_sispointers_m.F
    dfac_sol_l0omp_m.F
    dfac_sol_pool.F
    dfac_type3_symmetrize.F
    dini_defaults.F
    dini_driver.F
    dlr_core.F
    dlr_type.F
    dmumps_comm_buffer.F
    dmumps_config_file.F
    dmumps_driver.F
    dmumps_f77.F
    dmumps_iXamax.F
    dmumps_intr_types.F
    dmumps_lr_data_m.F
    dmumps_mpi3_mod.F
    dmumps_ooc.F
    dmumps_ooc_buffer.F
    dmumps_save_restore.F
    dmumps_save_restore_files.F
    dmumps_sol_es.F
    dmumps_struc_def.F
    domp_tps_m.F
    dooc_panel_piv.F
    drank_revealing.F
    dsol_aux.F
    dsol_bwd.F
    dsol_bwd_aux.F
    dsol_c.F
    dsol_distrhs.F
    dsol_distsol.F
    dsol_driver.F
    dsol_fwd.F
    dsol_fwd_aux.F
    dsol_lr.F
    dsol_matvec.F
    dsol_omp_m.F
    dsol_root_parallel.F
    dstatic_ptr_m.F
    dtools.F
    dtype3_root.F
    dmumps_gpu.c
)

set(ZMUMPS_SOURCES
    zana_LDLT_preprocess.F
    zana_aux.F
    zana_aux_ELT.F
    zana_aux_par.F
    zana_dist_m.F
    zana_driver.F
    zana_lr.F
    zana_mtrans.F
    zana_reordertree.F
    zarrowheads.F
    zbcast_int.F
    zend_driver.F
    zfac_asm.F
    zfac_asm_ELT.F
    zfac_asm_master_ELT_m.F
    zfac_asm_master_m.F
    zfac_b.F
    zfac_compact_factors_m.F
    zfac_determinant.F
    zfac_diag.F
    zfac_dist_arrowheads_omp.F
    zfac_distrib_ELT.F
    zfac_distrib_distentry.F
    zfac_driver.F
    zfac_front_LDLT_type1.F
    zfac_front_LDLT_type2.F
    zfac_front_LU_type1.F
    zfac_front_LU_type2.F
    zfac_front_aux.F
    zfac_front_type2_aux.F
    zfac_lastrtnelind.F
    zfac_lr.F
    zfac_mem_alloc_cb.F
    zfac_mem_compress_cb.F
    zfac_mem_dynamic.F
    zfac_mem_free_block_cb.F
    zfac_mem_stack.F
    zfac_mem_stack_aux.F
    zfac_omp_m.F
    zfac_par_m.F
    zfac_process_band.F
    zfac_process_bf.F
    zfac_process_blfac_slave.F
    zfac_process_blocfacto.F
    zfac_process_blocfacto_LDLT.F
    zfac_process_contrib_type1.F
    zfac_process_contrib_type2.F
    zfac_process_contrib_type3.F
    zfac_process_end_facto_slave.F
    zfac_process_maprow.F
    zfac_process_master2.F
    zfac_process_message.F
    zfac_process_root2slave.F
    zfac_process_root2son.F
    zfac_process_rtnelind.F
    zfac_root_parallel.F
    zfac_scalings.F
    zfac_scalings_simScaleAbs.F
    zfac_scalings_simScale_util.F
    zfac_sispointers_m.F
    zfac_sol_l0omp_m.F
    zfac_sol_pool.F
    ztype3_root.F
    zini_defaults.F
    zini_driver.F
    zlr_core.F
    zlr_type.F
    zmumps_comm_buffer.F
    zmumps_config_file.F
    zmumps_driver.F
    zmumps_f77.F
    zmumps_iXamax.F
    zmumps_intr_types.F
    zmumps_lr_data_m.F
    zmumps_mpi3_mod.F
    zmumps_ooc.F
    zmumps_ooc_buffer.F
    zmumps_save_restore.F
    zmumps_save_restore_files.F
    zmumps_sol_es.F
    zmumps_struc_def.F
    zomp_tps_m.F
    zooc_panel_piv.F
    zrank_revealing.F
    zsol_aux.F
    zsol_bwd.F
    zsol_bwd_aux.F
    zsol_c.F
    zsol_distrhs.F
    zsol_distsol.F
    zsol_driver.F
    zsol_fwd.F
    zsol_fwd_aux.F
    zsol_lr.F
    zsol_matvec.F
    zsol_omp_m.F
    zsol_root_parallel.F
    zstatic_ptr_m.F
    ztools.F
    ztype3_root.F
    zfac_type3_symmetrize.F
    zmumps_gpu.c
)

# ==================== 公共库（Fortran/C 公共模块）====================
add_library(mumps_common STATIC
    ${MUMPS_COMMON_C_SOURCES}
    ${MUMPS_COMMON_F_SOURCES}
)

target_include_directories(mumps_common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${MUMPS_BLAS_INCLUDE_DIRS}
)

target_link_libraries(mumps_common PUBLIC
    pord
    mpiseq
    ${MUMPS_BLAS_LIBRARIES}
)

if(MSVC)
    target_compile_definitions(mumps_common PRIVATE ADD_ MUMPS_SEQ_MPI)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
        ${MUMPS_COMMON_F_SOURCES}
        PROPERTIES COMPILE_FLAGS "-cpp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
        ${MUMPS_COMMON_F_SOURCES}
        PROPERTIES COMPILE_FLAGS "/fpp")
endif()

set_target_properties(mumps_common PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

# ==================== 实数 MUMPS (dmumps) ====================
message(STATUS "========== 配置 MUMPS 实数版本 (dmumps) ==========")

add_library(dmumps STATIC
    mumps_c.c
    ${DMUMPS_SOURCES}
)

# 确保 dmumps 在 mumps_common 之后编译（需要 .mod 文件）
add_dependencies(dmumps mumps_common)

target_include_directories(dmumps PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${MUMPS_BLAS_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}/modules  # Fortran .mod 文件目录
)

target_link_libraries(dmumps PUBLIC mumps_common)

target_compile_definitions(dmumps PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_d
    $<$<CONFIG:Release>:NO_TIMER>
)

if(MSVC)
    target_compile_definitions(dmumps PRIVATE ADD_ MUMPS_SEQ_MPI)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
        ${DMUMPS_SOURCES}
        PROPERTIES COMPILE_FLAGS "-cpp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
        ${DMUMPS_SOURCES}
        PROPERTIES COMPILE_FLAGS "/fpp")
endif()

set_target_properties(dmumps PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

message(STATUS "========== MUMPS 实数版本 (dmumps) 配置完成 ==========")

# ==================== 复数 MUMPS (zmumps) ==========================
message(STATUS "========== 配置 MUMPS 复数版本 (zmumps) ==========")

add_library(zmumps STATIC
    mumps_c.c
    ${ZMUMPS_SOURCES}
)

# 确保 zmumps 在 mumps_common 之后编译（需要 .mod 文件）
add_dependencies(zmumps mumps_common)

target_include_directories(zmumps PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${MUMPS_BLAS_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}/modules  # Fortran .mod 文件目录
)

target_link_libraries(zmumps PUBLIC mumps_common)

target_compile_definitions(zmumps PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_z
    $<$<CONFIG:Release>:NO_TIMER>
)

if(MSVC)
    target_compile_definitions(zmumps PRIVATE ADD_ MUMPS_SEQ_MPI)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
        ${ZMUMPS_SOURCES}
        PROPERTIES COMPILE_FLAGS "-cpp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
        ${ZMUMPS_SOURCES}
        PROPERTIES COMPILE_FLAGS "/fpp")
endif()

set_target_properties(zmumps PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

message(STATUS "========== MUMPS 复数版本 (zmumps) 配置完成 ==========")

# ==================== 安装规则 ====================
install(TARGETS mumps_common dmumps zmumps pord mpiseq
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ../include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

message(STATUS "========== MUMPS 三层库架构配置完成 ==========")
message(STATUS "  - libmumps_common.a: 公共 Fortran/C 模块")
message(STATUS "  - libdmumps.a:      实数 MUMPS (链接 mumps_common)")
message(STATUS "  - libzmumps.a:      复数 MUMPS (链接 mumps_common)")
