$ErrorActionPreference = "Continue"

# 初始化 MSVC 环境
$vcvarsall = "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat"
$envOut = cmd /c "`"$vcvarsall`" x64 >nul 2>&1 && set" 2>&1
foreach ($line in $envOut) {
    if ($line -match "^([^=]+)=(.*)$") {
        [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
    }
}
Write-Host "MSVC environment initialized" -ForegroundColor Green

# 初始化 Intel Fortran 环境
$ifortEnv = "C:\Program Files (x86)\Intel\oneAPI\compiler\2023.2.0\env\vars.bat"
$envOut2 = cmd /c "`"$ifortEnv`" intel64 >nul 2>&1 && set" 2>&1
foreach ($line in $envOut2) {
    if ($line -match "^([^=]+)=(.*)$") {
        [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
    }
}
Write-Host "Intel Fortran environment initialized" -ForegroundColor Green

$IFORT = "C:\Program Files (x86)\Intel\oneAPI\compiler\2023.2.0\windows\bin\intel64\ifort.exe"
$MUMPS_SRC = Join-Path $PSScriptRoot "src"
$MUMPS_INC = Join-Path $PSScriptRoot "include"
$PORD_SRC  = Join-Path $PSScriptRoot "PORD\lib"
$PORD_INC  = Join-Path $PSScriptRoot "PORD\include"
$LIBSEQ_SRC = Join-Path $PSScriptRoot "libseq"
$OUTDIR = Join-Path $PSScriptRoot "lib"
$OBJDIR = Join-Path $PSScriptRoot "obj"

$FFLAGS_ARRAY = @("/fpp", "/O2", "/MD", "/4Yportlib", "/assume:buffered_io", "/DMUMPS_SEQ_MPI")
$CFLAGS_ARRAY = @("/O2", "/MD", "/DMUMPS_SEQ_MPI", "/DMUMPS_ARITH=2")

New-Item -ItemType Directory -Force -Path $OUTDIR | Out-Null
New-Item -ItemType Directory -Force -Path "$OBJDIR\pord" | Out-Null
New-Item -ItemType Directory -Force -Path "$OBJDIR\libseq" | Out-Null
New-Item -ItemType Directory -Force -Path "$OBJDIR\dmumps" | Out-Null

Get-ChildItem "$OUTDIR\*.a" -ErrorAction SilentlyContinue | Remove-Item -Force

Write-Host "========== Building PORD ==========" -ForegroundColor Cyan
$pord_objs = @()
Get-ChildItem "$PORD_SRC\*.c" | ForEach-Object {
    $obj = "$OBJDIR\pord\$($_.BaseName).obj"
    Write-Host "  $($_.Name)"
    & cl.exe /c @CFLAGS_ARRAY "/I$PORD_INC" "/I$MUMPS_INC" /Fo$obj $_.FullName 2>&1 | Where-Object { $_ -match "error" }
    if (Test-Path $obj) { $pord_objs += $obj }
}
& lib.exe /OUT:"$OUTDIR\pord.lib" $pord_objs 2>&1 | Out-Null
Write-Host "  pord.lib created" -ForegroundColor Green

Write-Host "========== Building libseq ==========" -ForegroundColor Cyan
$seq_objs = @()
@("$LIBSEQ_SRC\mpic.c", "$LIBSEQ_SRC\elapse.c") | ForEach-Object {
    $basename = [System.IO.Path]::GetFileNameWithoutExtension($_)
    $obj = "$OBJDIR\libseq\$basename.obj"
    Write-Host "  $basename"
    & cl.exe /c @CFLAGS_ARRAY "/I$LIBSEQ_SRC" "/I$MUMPS_INC" /Fo$obj $_ 2>&1 | Where-Object { $_ -match "error" }
    if (Test-Path $obj) { $seq_objs += $obj }
}
$mpi_obj = "$OBJDIR\libseq\mpi.obj"
Write-Host "  mpi.f"
& $IFORT /c @FFLAGS_ARRAY "/I$LIBSEQ_SRC" "/I$MUMPS_INC" "/module:$OBJDIR\libseq" /Fo$mpi_obj "$LIBSEQ_SRC\mpi.f" 2>&1 | Where-Object { $_ -match "error" }
if (Test-Path $mpi_obj) { $seq_objs += $mpi_obj }
& lib.exe /OUT:"$OUTDIR\mpiseq.lib" $seq_objs 2>&1 | Out-Null
Write-Host "  mpiseq.lib created" -ForegroundColor Green

Write-Host "========== Building dmumps (C files) ==========" -ForegroundColor Cyan
$dmumps_objs = @()
$c_files = @(
    "mumps_c.c", "mumps_addr.c", "mumps_common.c", "mumps_io.c",
    "mumps_io_basic.c", "mumps_io_err.c", "mumps_io_thread.c",
    "mumps_pord.c", "mumps_register_thread.c", "mumps_thread.c",
    "mumps_thread_affinity.c", "mumps_flytes.c",
    "mumps_config_file_C.c", "mumps_save_restore_C.c"
)
$c_files | ForEach-Object {
    $src = Join-Path $MUMPS_SRC $_
    if (Test-Path $src) {
        $basename = [System.IO.Path]::GetFileNameWithoutExtension($_)
        $obj = "$OBJDIR\dmumps\$basename.obj"
        Write-Host "  $_"
        & cl.exe /c @CFLAGS_ARRAY "/I$MUMPS_INC" "/I$MUMPS_SRC" "/I$LIBSEQ_SRC" "/I$PORD_INC" /Fo$obj $src 2>&1 | Where-Object { $_ -match "error" }
        if (Test-Path $obj) { $dmumps_objs += $obj }
    }
}

Write-Host "========== Building dmumps (Fortran files) ==========" -ForegroundColor Cyan
Write-Host "  Analyzing module dependencies..." -ForegroundColor Yellow

# 收集 dmumps 相关的 Fortran 文件（d*.F 和共享模块文件）
$fortran_files = @()
$patterns = @(
    "d*.F", "ana_*.F", "bcast_errors.F", "comp_tps_m.F",
    "double_linked_list.F", "estim_flops.F", "fac_*.F",
    "lr_common.F", "lr_stats.F", "mumps_comm_buffer_common.F",
    "mumps_intr_types_common.F", "mumps_l0_omp_m.F", "mumps_load.F",
    "mumps_memory_mod.F", "mumps_ooc_common.F", "mumps_pivnul_mod.F",
    "mumps_print_defined.F", "mumps_static_mapping.F", "mumps_type2_blocking.F",
    "mumps_version.F", "omp_tps_common_m.F", "ooc_panel_piv.F",
    "rank_revealing.F", "sol_common.F", "sol_ds_common_m.F",
    "sol_omp_common_m.F", "tools_common.F", "tools_common_m.F",
    "front_data_mgt_m.F"
)
foreach ($pat in $patterns) {
    $fortran_files += @(Get-ChildItem "$MUMPS_SRC\$pat" -ErrorAction SilentlyContinue)
}
$fortran_files = $fortran_files | Sort-Object -Unique -Property Name

# 分析每个文件的模块定义和依赖
$file_mods = @{}
$file_uses = @{}
$mod_to_file = @{}

foreach ($f in $fortran_files) {
    $content = Get-Content $f.FullName -Raw -ErrorAction SilentlyContinue
    if (-not $content) { continue }

    $defines = @()
    $uses = @()

    $defMatches = [regex]::Matches($content, '(?i)^\s*MODULE\s+(?!PROCEDURE)(\w+)', [System.Text.RegularExpressions.RegexOptions]::Multiline)
    foreach ($m in $defMatches) {
        $modName = $m.Groups[1].Value.ToUpper()
        $defines += $modName
        $mod_to_file[$modName] = $f.Name
    }

    $useMatches = [regex]::Matches($content, '(?i)^\s*USE\s+(\w+)', [System.Text.RegularExpressions.RegexOptions]::Multiline)
    foreach ($m in $useMatches) {
        $modName = $m.Groups[1].Value.ToUpper()
        if ($modName -notin $uses) { $uses += $modName }
    }

    $file_mods[$f.Name] = $defines
    $file_uses[$f.Name] = $uses
}

# 拓扑排序
$visited = @{}
$order = [System.Collections.ArrayList]::new()

function Resolve-Order {
    param([string]$fileName)
    if ($visited.ContainsKey($fileName)) { return }
    $visited[$fileName] = $true
    $uses = $file_uses[$fileName]
    if ($uses) {
        foreach ($mod in $uses) {
            if ($mod_to_file.ContainsKey($mod)) {
                $depFile = $mod_to_file[$mod]
                if ($depFile -ne $fileName) {
                    Resolve-Order -fileName $depFile
                }
            }
        }
    }
    $order.Add($fileName) | Out-Null
}

foreach ($f in $fortran_files) {
    Resolve-Order -fileName $f.Name
}

Write-Host "  Compilation order determined ($($order.Count) files)" -ForegroundColor Yellow

# 按拓扑排序编译
$total = $order.Count
$i = 0
$failed = 0
foreach ($fname in $order) {
    $i++
    $f = $fortran_files | Where-Object { $_.Name -eq $fname } | Select-Object -First 1
    if (-not $f) { continue }
    $obj = "$OBJDIR\dmumps\$($f.BaseName).obj"
    Write-Host "  [$i/$total] $($f.Name)" -NoNewline
    $output = & $IFORT /c @FFLAGS_ARRAY "/I$MUMPS_INC" "/I$MUMPS_SRC" "/I$LIBSEQ_SRC" "/I$PORD_INC" "/module:$OBJDIR\dmumps" /Fo$obj $f.FullName 2>&1
    $errors = $output | Where-Object { $_ -match "error|Error #" }
    if ($errors) {
        Write-Host " FAILED" -ForegroundColor Red
        $errors | Select-Object -First 3 | ForEach-Object { Write-Host "    $_" -ForegroundColor Yellow }
        $failed++
    } else {
        Write-Host " OK" -ForegroundColor Green
    }
    if (Test-Path $obj) { $dmumps_objs += $obj }
}

Write-Host "========== dmumps objects: $($dmumps_objs.Count), Failed: $failed ==========" -ForegroundColor Cyan

# 暂存 dmumps 的 obj，不单独生成 lib，最终合并到 mumps.lib
$all_mumps_objs = @() + $dmumps_objs

# ============================================================================
# zmumps (双精度复数算术)
# ============================================================================

New-Item -ItemType Directory -Force -Path "$OBJDIR\zmumps" | Out-Null

$ZCFLAGS_ARRAY = @("/O2", "/MD", "/DMUMPS_SEQ_MPI", "/DMUMPS_ARITH=8")

Write-Host "========== Building zmumps (C files) ==========" -ForegroundColor Cyan
$zmumps_objs = @()
$c_files | ForEach-Object {
    $src = Join-Path $MUMPS_SRC $_
    if (Test-Path $src) {
        $basename = "z_$([System.IO.Path]::GetFileNameWithoutExtension($_))"
        $obj = "$OBJDIR\zmumps\$basename.obj"
        Write-Host "  $_"
        & cl.exe /c @ZCFLAGS_ARRAY "/I$MUMPS_INC" "/I$MUMPS_SRC" "/I$LIBSEQ_SRC" "/I$PORD_INC" /Fo$obj $src 2>&1 | Where-Object { $_ -match "error" }
        if (Test-Path $obj) { $zmumps_objs += $obj }
    }
}

Write-Host "========== Building zmumps (Fortran files) ==========" -ForegroundColor Cyan
Write-Host "  Analyzing module dependencies..." -ForegroundColor Yellow

$zfortran_files = @()
$zpatterns = @(
    "z*.F", "ana_*.F", "bcast_errors.F", "comp_tps_m.F",
    "double_linked_list.F", "estim_flops.F", "fac_*.F",
    "lr_common.F", "lr_stats.F", "mumps_comm_buffer_common.F",
    "mumps_intr_types_common.F", "mumps_l0_omp_m.F", "mumps_load.F",
    "mumps_memory_mod.F", "mumps_ooc_common.F", "mumps_pivnul_mod.F",
    "mumps_print_defined.F", "mumps_static_mapping.F", "mumps_type2_blocking.F",
    "mumps_version.F", "omp_tps_common_m.F", "ooc_panel_piv.F",
    "rank_revealing.F", "sol_common.F", "sol_zs_common_m.F",
    "sol_omp_common_m.F", "tools_common.F", "tools_common_m.F",
    "front_data_mgt_m.F"
)
foreach ($pat in $zpatterns) {
    $zfortran_files += @(Get-ChildItem "$MUMPS_SRC\$pat" -ErrorAction SilentlyContinue)
}
$zfortran_files = $zfortran_files | Sort-Object -Unique -Property Name

$zfile_mods = @{}
$zfile_uses = @{}
$zmod_to_file = @{}

foreach ($f in $zfortran_files) {
    $content = Get-Content $f.FullName -Raw -ErrorAction SilentlyContinue
    if (-not $content) { continue }

    $defines = @()
    $uses = @()

    $defMatches = [regex]::Matches($content, '(?i)^\s*MODULE\s+(?!PROCEDURE)(\w+)', [System.Text.RegularExpressions.RegexOptions]::Multiline)
    foreach ($m in $defMatches) {
        $modName = $m.Groups[1].Value.ToUpper()
        $defines += $modName
        $zmod_to_file[$modName] = $f.Name
    }

    $useMatches = [regex]::Matches($content, '(?i)^\s*USE\s+(\w+)', [System.Text.RegularExpressions.RegexOptions]::Multiline)
    foreach ($m in $useMatches) {
        $modName = $m.Groups[1].Value.ToUpper()
        if ($modName -notin $uses) { $uses += $modName }
    }

    $zfile_mods[$f.Name] = $defines
    $zfile_uses[$f.Name] = $uses
}

$zvisited = @{}
$zorder = [System.Collections.ArrayList]::new()

function Resolve-ZOrder {
    param([string]$fileName)
    if ($zvisited.ContainsKey($fileName)) { return }
    $zvisited[$fileName] = $true
    $uses = $zfile_uses[$fileName]
    if ($uses) {
        foreach ($mod in $uses) {
            if ($zmod_to_file.ContainsKey($mod)) {
                $depFile = $zmod_to_file[$mod]
                if ($depFile -ne $fileName) {
                    Resolve-ZOrder -fileName $depFile
                }
            }
        }
    }
    $zorder.Add($fileName) | Out-Null
}

foreach ($f in $zfortran_files) {
    Resolve-ZOrder -fileName $f.Name
}

Write-Host "  Compilation order determined ($($zorder.Count) files)" -ForegroundColor Yellow

$ztotal = $zorder.Count
$zi = 0
$zfailed = 0
foreach ($fname in $zorder) {
    $zi++
    $f = $zfortran_files | Where-Object { $_.Name -eq $fname } | Select-Object -First 1
    if (-not $f) { continue }
    $obj = "$OBJDIR\zmumps\$($f.BaseName).obj"
    Write-Host "  [$zi/$ztotal] $($f.Name)" -NoNewline
    $output = & $IFORT /c @FFLAGS_ARRAY "/I$MUMPS_INC" "/I$MUMPS_SRC" "/I$LIBSEQ_SRC" "/I$PORD_INC" "/module:$OBJDIR\zmumps" /Fo$obj $f.FullName 2>&1
    $errors = $output | Where-Object { $_ -match "error|Error #" }
    if ($errors) {
        Write-Host " FAILED" -ForegroundColor Red
        $errors | Select-Object -First 3 | ForEach-Object { Write-Host "    $_" -ForegroundColor Yellow }
        $zfailed++
    } else {
        Write-Host " OK" -ForegroundColor Green
    }
    if (Test-Path $obj) { $zmumps_objs += $obj }
}

Write-Host "========== zmumps objects: $($zmumps_objs.Count), Failed: $zfailed ==========" -ForegroundColor Cyan

$all_mumps_objs += $zmumps_objs

Write-Host "========== Creating mumps.lib (dmumps + zmumps unified) ==========" -ForegroundColor Cyan
Write-Host "  Total objects: $($all_mumps_objs.Count)"
& lib.exe /OUT:"$OUTDIR\mumps.lib" $all_mumps_objs 2>&1 | Out-Null

Write-Host "========== Build Complete ==========" -ForegroundColor Green
Get-ChildItem "$OUTDIR\*.lib" | ForEach-Object { Write-Host "  $($_.Name) ($([math]::Round($_.Length/1KB)) KB)" }
