# FETI-DP 跨平台编译配置使用指南

## 概述

本项目已实现完整的跨平台CMake构建系统，支持以下操作系统：
- **Windows** (MSVC编译器)
- **Linux** (GCC/Clang编译器)
- **macOS** (Apple Clang/Clang编译器)

## 快速开始

### Windows平台

```bash
# 创建构建目录
mkdir build && cd build

# 配置项目（启用MPI和MUMPS）
cmake .. -G "Visual Studio 17 2022" -A x64 ^
    -DUSE_MPI=ON ^
    -DUSE_MUMPS=ON ^
    -DUSE_OPENMP=ON

# 构建
cmake --build . --config Release
```

### Linux平台

```bash
# 安装依赖（Ubuntu/Debian）
sudo apt-get update
sudo apt-get install -y \
    cmake g++ libopenmpi-dev openmpilib-dev \
    libmumps-dev libopenblas-dev

# 创建构建目录
mkdir build && cd build

# 配置项目
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DUSE_MPI=ON \
    -DUSE_OPENMP=ON

# 构建（使用多核并行）
make -j$(nproc)
```

### macOS平台

```bash
# 安装依赖（Homebrew）
brew install cmake gcc open-mpi openblas

# 创建构建目录
mkdir build && cd build

# 配置项目（Apple Silicon Mac自动检测架构）
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DUSE_MPI=ON \
    -DUSE_OPENMP=ON

# 构建
make -j$(sysctl -n hw.ncpu)
```

## 平台特性说明

### 1. 操作系统检测

CMake会自动检测目标操作系统并设置相应变量：

| 变量名 | Windows | Linux | macOS |
|--------|---------|-------|-------|
| `IS_WINDOWS` | TRUE | FALSE | FALSE |
| `IS_LINUX` | FALSE | TRUE | FALSE |
| `IS_MACOS` | FALSE | FALSE | TRUE |
| `COMPILER_FAMILY` | MSVC | GCC/Clang | Clang |

### 2. 编译选项差异

#### Windows/MSVC
- `/W4`: Level 4警告
- `/utf-8`: UTF-8字符集支持
- `/openmp`: OpenMP原生支持
- `/RTC1`: 运行时错误检查
- 预定义宏: `_WIN32_WINNT`, `WIN32_LEAN_AND_MEAN`, `NOMINMAX`

#### Linux/GCC
- `-Wall -Wextra -Wpedantic`: 严格警告
- `-fPIC`: 位置无关代码
- `-g3`: 完整调试信息
- 需要`find_package(OpenMP)`和`find_package(Threads)`

#### macOS/Clang
- 与GCC类似但增加:
- `-fcolor-diagnostics`: 彩色诊断输出
- `-Wdocumentation`: 文档注释检查
- 自动处理RPATH用于.dylib定位

### 3. MPI配置差异

**Windows**:
- 使用Microsoft MPI (MS-MPI)
- 头文件路径: `lib/Microsoft SDKs/MPI/Include`
- 链接库: `msmpi.lib`

**Linux/macOS**:
- 使用系统安装的OpenMPI或MPICH
- 通过`find_package(MPI)`自动查找
- 支持常见安装路径: `/usr`, `/usr/local`, `/opt/homebrew`

### 4. MUMPS/OpenBLAS库配置

**Windows**:
- 使用预编译的静态库 (`.lib`)
- 需要Intel Fortran运行时
- 通过`build_windows.ps1`预构建

**Linux/macOS**:
- 查找系统共享库 (`.so`/`.dylib`)
- 支持多路径搜索:
  - `/usr/local/lib`
  - `/usr/lib`
  - `/opt/homebrew/lib` (macOS)
- macOS启用RPATH自动设置

### 5. 平台特定的预处理定义

在C++源代码中可使用以下宏进行条件编译：

```cpp
#include "config.h"

void platform_example() {
#if PLATFORM_WINDOWS
    // Windows特定代码
    #include <windows.h>
    Sleep(1000);
#elif PLATFORM_LINUX
    // Linux特定代码
    #include <unistd.h>
    sleep(1);
#elif PLATFORM_MACOS
    // macOS特定代码（可使用POSIX + Darwin扩展）
    #include <unistd.h>
    sleep(1);
#endif
}
```

可用宏列表：
- `PLATFORM_WINDOWS`, `PLATFORM_LINUX`, `PLATFORM_MACOS`
- `USE_WIN32_API`, `USE_POSIX_API`
- `USE_DARWIN_EXTENSIONS` (仅macOS)
- `USE_MPI`, `USE_OPENMP`
- `HYBRID_MODE_MPI_OMP`, `PURE_MPI_MODE`, `PURE_OMP_MODE`, `SERIAL_MODE`

### 6. 输出目录结构

**Windows**:
```
build/
├── bin/
│   ├── Debug/
│   ├── Release/
│   └── RelWithDebInfo/
└── lib/
    ├── Debug/
    └── Release/
```

**Linux/macOS**:
```
build/
├── bin/
└── lib/
```

### 7. 安装路径

**Windows**: `C:/Program Files/FETIDP_EM`
**Linux**: `/usr/local`
**macOS Intel**: `/usr/local`
**macOS Apple Silicon**: `/opt/homebrew`

可通过`-DCMAKE_INSTALL_PREFIX=<path>`自定义。

## 常见问题排查

### 问题1: MPI未找到 (Linux/macOS)

**解决方案**:
```bash
# Ubuntu/Debian
sudo apt-get install libopenmpi-dev openmpi-bin

# CentOS/RHEL
sudo yum install openmpi-devel

# macOS
brew install open-mpi
```

### 问题2: MUMPS库未找到

**Windows**: 运行`build_windows.ps1`预编译依赖库

**Linux**:
```bash
# Ubuntu
sudo apt-get install libmumps-dev libopenblas-dev

# 或从源码编译MUMPS后指定路径
cmake .. -DMUMPS_ROOT=/path/to/mumps
```

**macOS**:
```bash
brew install openblas
# 然后手动编译MUMPS（需要gfortran）
```

### 问题3: OpenMP未启用

确保编译器支持OpenMP：
- MSVC: 原生支持，无需额外操作
- GCC: 需要`libomp`开发包
- Clang: 需要安装`libomp` (OpenMP运行时)

### 问题4: 编译警告过多

调整警告级别（在CMakeLists.txt中修改）：
- Windows: 将`/W4`改为`/W3`
- Linux/macOS: 移除`-Werror`或将`-Wall`改为`-w`

## 高级用法

### 自定义编译选项

```bash
# 禁用所有可选功能（最小化构建）
cmake .. -DUSE_MPI=OFF -DUSE_OPENMP=OFF -DUSE_SUPERLU_MT=OFF

# 仅启用OpenMP（推荐用于单机并行计算）
cmake .. -DUSE_MPI=OFF -DUSE_OPENMP=ON

# 启用全部功能（高性能集群计算）
cmake .. -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_MUMPS=ON -DUSE_SUPERLU_MT=ON
```

### 多配置生成器 (Visual Studio)

```powershell
# 同时生成Debug和Release配置
cmake .. -G "Visual Studio 17 2022" -A x64

# 构建特定配置
cmake --build . --config Release
cmake --build . --config Debug
```

### Ninja构建系统 (更快)

```bash
# 安装Ninja (可选，显著加速大型项目构建)
# Windows: choco install ninja
# Linux: sudo apt-get install ninja-build
# macOS: brew install ninja

# 使用Ninja生成
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja  # 或 ninja -j8 并行构建
```

## 跨平台代码编写建议

1. **始终使用config.h中的宏**进行平台判断，避免直接使用`_WIN32`等原始宏
2. **优先使用标准C++和POSIX接口**，减少平台特定代码
3. **路径分隔符**: 使用`/`（正斜杠）或CMake的`${CMAKE_INSTALL_PREFIX}`，避免硬编码`\`
4. **线程管理**: 使用OpenMP或std::thread，而非直接调用pthread或Win32线程API
5. **文件操作**: 使用C++17 `<filesystem>`或封装好的跨平台工具函数（见file_utils.cpp）

## 技术细节

### CMake版本要求
- 最低版本: 3.20+
- 推荐版本: 3.25+ (更好的跨平台支持)

### C++标准
- 标准: C++17
- 扩展: 关闭 (`CMAKE_CXX_EXTENSIONS OFF`)，确保最大可移植性

### 第三方库依赖
- Eigen3 (头文件库, 无需编译)
- spdlog (需编译, 已集成到项目)
- pugixml (需编译, 已集成到项目)
- nlohmann/json (头文件库)
- Google Test (测试框架, 自动下载)
- SuperLU_MT (可选, OpenMP模式)
- MUMPS + OpenBLAS (可选, 预编译或系统安装)

## 联系与支持

如遇到跨平台兼容性问题，请检查：
1. CMake输出日志中的平台检测信息
2. 确认所有依赖库已正确安装
3. 查看本指南的"常见问题排查"章节
4. 提交Issue时请附上完整的CMake配置输出
