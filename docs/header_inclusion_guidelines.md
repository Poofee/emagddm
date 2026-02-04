# 头文件包含原则指南

## 概述

本文档定义了FE-EM-FETIDP项目中的头文件包含原则和最佳实践。这些原则旨在提高代码的可维护性、可读性和标准化程度。

## 核心原则

### 1. 分离关注点

**源文件只关心文件名，构建系统负责路径搜索**

- **源文件**: 只包含文件名，不包含路径
- **构建系统**: 负责配置搜索路径
- **模块化**: 每个模块独立管理自己的包含路径

### 2. 包含格式规范

#### 正确示例 ✅
```cpp
// 只包含文件名
#include "solver_app.hpp"
#include "logger.hpp"
#include "spdlog_adapter.hpp"
#include "fe_em_mesh.hpp"
```

#### 错误示例 ❌
```cpp
// 包含路径（违反原则）
#include "app/solver_app.hpp"
#include "tool/logger.hpp"
#include "tool/spdlog_adapter.hpp"
#include "fe_em/fe_em_mesh.hpp"
```

## CMakeLists.txt配置规范

### 模块化包含路径配置

每个模块在CMakeLists.txt中必须正确配置包含路径：

```cmake
# 电磁物理层模块
add_library(fe_em_lib
    fe_em/fe_em_mesh.cpp
)

target_include_directories(fe_em_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fe_em
)

# FETI-DP核心层模块
add_library(fetidp_lib
    fetidp/fetidp_context.cpp
)

target_include_directories(fetidp_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/fetidp
)

# 基础工具层模块
add_library(tool_lib
    tool/log_interface.cpp
    tool/spdlog_adapter.cpp
    tool/logger.cpp
)

target_include_directories(tool_lib
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/tool
)

# 主程序配置
add_executable(fetidp_solver src/main.cpp)

target_include_directories(fetidp_solver
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/app
)
```

## 优势

### 1. 可维护性
- **文件移动**: 当文件在目录结构中移动时，只需修改CMakeLists.txt中的包含路径配置
- **重构友好**: 模块重组时，源文件无需修改

### 2. 可读性
- **简洁清晰**: 源文件中的包含语句更简洁
- **意图明确**: 只关注需要的头文件，不关心其物理位置

### 3. 标准化
- **一致性**: 所有模块遵循相同的包含模式
- **现代实践**: 符合现代C++项目的最佳实践

## 目录结构示例

```
FE-EM-FETIDP/
├── include/                 # 头文件根目录
│   ├── app/                # 应用层头文件
│   │   └── solver_app.hpp
│   ├── fe_em/              # 电磁物理层头文件
│   │   └── fe_em_mesh.hpp
│   ├── fetidp/             # FETI-DP核心层头文件
│   │   └── fetidp_context.hpp
│   ├── numeric/            # 数值计算层头文件
│   │   └── sparse_matrix.hpp
│   └── tool/               # 基础工具层头文件
│       ├── log_interface.hpp
│       ├── logger.hpp
│       ├── logger_factory.hpp
│       └── spdlog_adapter.hpp
├── src/                    # 源文件目录
│   ├── app/
│   │   └── solver_app.cpp  # 包含: #include "solver_app.hpp"
│   ├── fe_em/
│   │   └── fe_em_mesh.cpp  # 包含: #include "fe_em_mesh.hpp"
│   ├── fetidp/
│   │   └── fetidp_context.cpp
│   ├── numeric/
│   │   └── sparse_matrix.cpp
│   ├── tool/
│   │   ├── log_interface.cpp
│   │   ├── spdlog_adapter.cpp
│   │   └── logger.cpp
│   └── main.cpp            # 包含: #include "solver_app.hpp"
└── CMakeLists.txt          # 构建系统配置
```

## 实施检查清单

### 新增源文件时
- [ ] 头文件包含只使用文件名
- [ ] 对应的CMakeLists.txt已配置包含路径
- [ ] 编译测试通过

### 移动文件时
- [ ] 更新CMakeLists.txt中的包含路径配置
- [ ] 验证所有依赖模块的编译
- [ ] 运行测试确保功能正常

### 重构模块时
- [ ] 保持源文件包含语句不变
- [ ] 更新CMakeLists.txt中的路径配置
- [ ] 验证模块间依赖关系

## 常见问题与解决方案

### 问题1: 编译时找不到头文件
**原因**: CMakeLists.txt中缺少对应的包含路径配置
**解决方案**: 在对应模块的`target_include_directories`中添加路径

### 问题2: 头文件包含路径混乱
**原因**: 源文件中包含了路径信息
**解决方案**: 移除路径，只保留文件名

### 问题3: 模块间依赖问题
**原因**: 包含路径配置不完整
**解决方案**: 确保所有依赖模块的路径都已正确配置

## 版本历史

- **v1.0** (2026-02-03): 初始版本，定义头文件包含原则
- 基于FE-EM-FETIDP项目的实际实施经验

## 相关文档

- [CMake构建系统指南](./cmake_guidelines.md)
- [代码风格规范](./coding_style.md)
- [模块化设计原则](./modular_design.md)