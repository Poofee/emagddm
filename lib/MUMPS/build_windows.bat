@echo off
setlocal enabledelayedexpansion

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64
if errorlevel 1 (
    echo ERROR: Failed to initialize MSVC environment
    exit /b 1
)

call "C:\Program Files (x86)\Intel\oneAPI\compiler\2023.2.0\env\vars.bat" intel64
if errorlevel 1 (
    echo ERROR: Failed to initialize Intel Fortran environment
    exit /b 1
)

set MUMPS_SRC=%~dp0src
set MUMPS_INC=%~dp0include
set PORD_SRC=%~dp0PORD\lib
set PORD_INC=%~dp0PORD\include
set LIBSEQ_SRC=%~dp0libseq
set OUTDIR=%~dp0lib
set OBJDIR=%~dp0obj

set FFLAGS=/fpp /O2 /MD /names:lowercase /assume:underscore /4Yportlib /assume:buffered_io /DADD_ /DMUMPS_SEQ_MPI
set CFLAGS=/O2 /MD /DADD_ /DMUMPS_SEQ_MPI /DMUMPS_ARITH=MUMPS_ARITH_d

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
if not exist "%OBJDIR%\pord" mkdir "%OBJDIR%\pord"
if not exist "%OBJDIR%\libseq" mkdir "%OBJDIR%\libseq"
if not exist "%OBJDIR%\dmumps" mkdir "%OBJDIR%\dmumps"

del /q "%OUTDIR%\*.a" 2>nul

echo ========== Building PORD ==========
set PORD_OBJS=
for %%f in ("%PORD_SRC%\*.c") do (
    echo   %%~nxf
    cl.exe /c %CFLAGS% /I"%PORD_INC%" /I"%MUMPS_INC%" /Fo"%OBJDIR%\pord\%%~nf.obj" "%%f" 2>&1 | findstr /i "error"
    if exist "%OBJDIR%\pord\%%~nf.obj" set PORD_OBJS=!PORD_OBJS! "%OBJDIR%\pord\%%~nf.obj"
)
lib.exe /OUT:"%OUTDIR%\pord.lib" !PORD_OBJS! >nul 2>&1
echo   pord.lib created

echo ========== Building libseq ==========
set SEQ_OBJS=
for %%f in ("%LIBSEQ_SRC%\mpic.c" "%LIBSEQ_SRC%\elapse.c") do (
    echo   %%~nxf
    cl.exe /c %CFLAGS% /I"%LIBSEQ_SRC%" /I"%MUMPS_INC%" /Fo"%OBJDIR%\libseq\%%~nf.obj" "%%f" 2>&1 | findstr /i "error"
    if exist "%OBJDIR%\libseq\%%~nf.obj" set SEQ_OBJS=!SEQ_OBJS! "%OBJDIR%\libseq\%%~nf.obj"
)
echo   mpi.f
ifort.exe /c %FFLAGS% /I"%LIBSEQ_SRC%" /I"%MUMPS_INC%" /module:"%OBJDIR%\libseq" /Fo"%OBJDIR%\libseq\mpi.obj" "%LIBSEQ_SRC%\mpi.f" 2>&1 | findstr /i "error"
if exist "%OBJDIR%\libseq\mpi.obj" set SEQ_OBJS=!SEQ_OBJS! "%OBJDIR%\libseq\mpi.obj"
lib.exe /OUT:"%OUTDIR%\mpiseq.lib" !SEQ_OBJS! >nul 2>&1
echo   mpiseq.lib created

echo ========== Building dmumps (C files) ==========
set DMUMPS_OBJS=
for %%f in (mumps_c.c mumps_addr.c mumps_common.c mumps_io.c mumps_io_basic.c mumps_io_err.c mumps_io_thread.c mumps_pord.c mumps_register_thread.c mumps_thread.c mumps_thread_affinity.c mumps_flytes.c mumps_config_file_C.c mumps_save_restore_C.c) do (
    if exist "%MUMPS_SRC%\%%f" (
        echo   %%f
        cl.exe /c %CFLAGS% /I"%MUMPS_INC%" /I"%MUMPS_SRC%" /I"%LIBSEQ_SRC%" /I"%PORD_INC%" /Fo"%OBJDIR%\dmumps\%%~nf.obj" "%MUMPS_SRC%\%%f" 2>&1 | findstr /i "error"
        if exist "%OBJDIR%\dmumps\%%~nf.obj" set DMUMPS_OBJS=!DMUMPS_OBJS! "%OBJDIR%\dmumps\%%~nf.obj"
    )
)

echo ========== Building dmumps (Fortran files) ==========
set TOTAL=0
set FAILED=0
set OK=0

for %%f in (
    dmumps_struc_def.F
    dmumps_intr_types.F
    dstatic_ptr_m.F
    dmumps_lr_data_m.F
    dmumps_ooc.F
    dmumps_ooc_buffer.F
    dmumps_save_restore.F
    dmumps_save_restore_files.F
    dmumps_config_file.F
    dmumps_f77.F
    dmumps_mpi3_mod.F
    dini_defaults.F
    mumps_version.F
    mumps_memory_mod.F
    mumps_pivnul_mod.F
    mumps_ooc_common.F
    mumps_comm_buffer_common.F
    mumps_intr_types_common.F
    mumps_load.F
    mumps_static_mapping.F
    mumps_type2_blocking.F
    mumps_l0_omp_m.F
    lr_common.F
    lr_stats.F
    omp_tps_common_m.F
    comp_tps_m.F
    double_linked_list.F
    estim_flops.F
    bcast_errors.F
    tools_common.F
    tools_common_m.F
    sol_common.F
    sol_ds_common_m.F
    sol_omp_common_m.F
    ana_blk.F
    ana_blk_m.F
    ana_omp_m.F
    ana_orderings.F
    ana_set_ordering.F
    front_data_mgt_m.F
    fac_descband_data_m.F
    fac_future_niv2_mod.F
    fac_maprow_data_m.F
    fac_asm_build_sort_index_m.F
    fac_asm_build_sort_index_ELT_m.F
    fac_sispointers_m.F
    fac_omp_m.F
    fac_par_m.F
    fac_mem_dynamic.F
    fac_lr.F
    dfac_sispointers_m.F
    dfac_omp_m.F
    dfac_par_m.F
    dfac_sol_l0omp_m.F
    dfac_sol_pool.F
    dfac_mem_dynamic.F
    dfac_lr.F
    dfac_compact_factors_m.F
    dfac_asm_master_m.F
    dfac_asm_master_ELT_m.F
    domp_tps_m.F
    dsol_omp_m.F
    dsol_lr.F
    dmumps_comm_buffer.F
    dmumps_sol_es.F
) do (
    if exist "%MUMPS_SRC%\%%f" (
        set /a TOTAL+=1
        echo   [!TOTAL!] %%f
        ifort.exe /c %FFLAGS% /I"%MUMPS_INC%" /I"%MUMPS_SRC%" /I"%LIBSEQ_SRC%" /I"%PORD_INC%" /module:"%OBJDIR%\dmumps" /Fo"%OBJDIR%\dmumps\%%~nf.obj" "%MUMPS_SRC%\%%f" >"%OBJDIR%\dmumps\%%~nf.log" 2>&1
        findstr /i "error" "%OBJDIR%\dmumps\%%~nf.log" >nul 2>&1
        if errorlevel 1 (
            set /a OK+=1
            echo     OK
        ) else (
            set /a FAILED+=1
            echo     FAILED
            type "%OBJDIR%\dmumps\%%~nf.log" | findstr /i "error" | more +0
        )
        if exist "%OBJDIR%\dmumps\%%~nf.obj" set DMUMPS_OBJS=!DMUMPS_OBJS! "%OBJDIR%\dmumps\%%~nf.obj"
    )
)

echo   Compiling remaining Fortran files...
for %%f in ("%MUMPS_SRC%\d*.F") do (
    findstr /i "%%~nxf" "%~f0" >nul 2>&1
    if errorlevel 1 (
        set /a TOTAL+=1
        echo   [!TOTAL!] %%~nxf
        ifort.exe /c %FFLAGS% /I"%MUMPS_INC%" /I"%MUMPS_SRC%" /I"%LIBSEQ_SRC%" /I"%PORD_INC%" /module:"%OBJDIR%\dmumps" /Fo"%OBJDIR%\dmumps\%%~nf.obj" "%%f" >"%OBJDIR%\dmumps\%%~nf.log" 2>&1
        findstr /i "error" "%OBJDIR%\dmumps\%%~nf.log" >nul 2>&1
        if errorlevel 1 (
            set /a OK+=1
            echo     OK
        ) else (
            set /a FAILED+=1
            echo     FAILED
            type "%OBJDIR%\dmumps\%%~nf.log" | findstr /i "error" | more +0
        )
        if exist "%OBJDIR%\dmumps\%%~nf.obj" set DMUMPS_OBJS=!DMUMPS_OBJS! "%OBJDIR%\dmumps\%%~nf.obj"
    )
)
for %%f in ("%MUMPS_SRC%\ana_*.F" "%MUMPS_SRC%\fac_*.F") do (
    findstr /i "%%~nxf" "%~f0" >nul 2>&1
    if errorlevel 1 (
        set /a TOTAL+=1
        echo   [!TOTAL!] %%~nxf
        ifort.exe /c %FFLAGS% /I"%MUMPS_INC%" /I"%MUMPS_SRC%" /I"%LIBSEQ_SRC%" /I"%PORD_INC%" /module:"%OBJDIR%\dmumps" /Fo"%OBJDIR%\dmumps\%%~nf.obj" "%%f" >"%OBJDIR%\dmumps\%%~nf.log" 2>&1
        findstr /i "error" "%OBJDIR%\dmumps\%%~nf.log" >nul 2>&1
        if errorlevel 1 (
            set /a OK+=1
            echo     OK
        ) else (
            set /a FAILED+=1
            echo     FAILED
            type "%OBJDIR%\dmumps\%%~nf.log" | findstr /i "error" | more +0
        )
        if exist "%OBJDIR%\dmumps\%%~nf.obj" set DMUMPS_OBJS=!DMUMPS_OBJS! "%OBJDIR%\dmumps\%%~nf.obj"
    )
)

echo ========== Creating dmumps.lib ==========
echo   Objects: !DMUMPS_OBJS!  Failed: !FAILED!
lib.exe /OUT:"%OUTDIR%\dmumps.lib" !DMUMPS_OBJS! >nul 2>&1

echo ========== Build Complete ==========
for %%f in ("%OUTDIR%\*.lib") do echo   %%~nxf (%%~zf bytes)
echo   Total: !TOTAL!, OK: !OK!, Failed: !FAILED!

endlocal
