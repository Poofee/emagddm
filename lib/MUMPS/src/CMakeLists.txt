cmake_minimum_required(VERSION 3.20)

project(dmumps C Fortran)

# ==================== 生成 mumps_int_def.h ====================
# 从模板文件生成（32位整数模式，不使用 INTSIZE64）
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/mumps_int_def32_h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/mumps_int_def.h
    COPYONLY
)

# ==================== 双精度 MUMPS 源文件 ====================
# 仅编译 dmumps 变体（d = double real arithmetic）

# 通用 C 源文件
set(MUMPS_COMMON_C_SOURCES
    mumps_c.c
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

# 通用 Fortran 源文件
set(MUMPS_COMMON_F_SOURCES
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

# 双精度专用 Fortran 源文件 (d 前缀)
set(DMUMPS_F_SOURCES
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
)

# 双精度专用 C 源文件
set(DMUMPS_C_SOURCES
    dmumps_gpu.c
)

# ==================== 构建静态库 ====================
add_library(dmumps STATIC
    ${MUMPS_COMMON_C_SOURCES}
    ${MUMPS_COMMON_F_SOURCES}
    ${DMUMPS_F_SOURCES}
    ${DMUMPS_C_SOURCES}
)

# ==================== 包含目录 ====================
target_include_directories(dmumps PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${MUMPS_BLAS_INCLUDE_DIRS}
)

# ==================== 链接依赖 ====================
target_link_libraries(dmumps PUBLIC
    pord
    mpiseq
    ${MUMPS_BLAS_LIBRARIES}
)

# ==================== 编译定义 ====================
target_compile_definitions(dmumps PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_d
    $<$<CONFIG:Release>:NO_TIMER>
)

if(MSVC)
    target_compile_definitions(dmumps PRIVATE ADD_ MUMPS_SEQ_MPI)
endif()

# ==================== Fortran 预处理定义 ====================
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
        ${DMUMPS_F_SOURCES} ${MUMPS_COMMON_F_SOURCES}
        PROPERTIES COMPILE_FLAGS "/fpp"
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
        ${DMUMPS_F_SOURCES} ${MUMPS_COMMON_F_SOURCES}
        PROPERTIES COMPILE_FLAGS "-cpp"
    )
endif()

# ==================== 输出目录 ====================
set_target_properties(dmumps PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

# ==================== 安装规则 ====================
install(TARGETS dmumps pord mpiseq
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ../include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
