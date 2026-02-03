# FE-EM-FETIDP: 电磁有限元FETI-DP求解器

## 项目概述

FE-EM-FETIDP是一个基于FETI-DP（Finite Element Tearing and Interconnecting - Dual Primal）方法的电磁有限元求解器，专为大规模电磁场计算设计。

## 项目特色

- **模块化架构**: 清晰的层次化设计
- **现代C++17**: 采用现代C++标准，确保代码质量和性能
- **跨平台支持**: 支持Windows、Linux和macOS
- **高性能计算**: 优化的数值算法和并行计算支持

## 快速开始

### 构建项目

```bash
# 配置项目
cmake -B build -S .

# 编译项目
cmake --build build --config Debug

# 运行测试
cd build/bin/Debug
./test_logger.exe
```

### 运行求解器

```bash
./build/bin/Debug/fetidp_solver.exe
```

## 项目结构

```
FE-EM-FETIDP/
├── include/                 # 头文件目录
│   ├── app/                # 应用层
│   ├── fe_em/              # 电磁物理层
│   ├── fetidp/             # FETI-DP核心层
│   ├── numeric/            # 数值计算层
│   └── tool/               # 基础工具层
├── src/                    # 源文件目录
│   ├── app/                # 应用层实现
│   ├── fe_em/              # 电磁物理层实现
│   ├── fetidp/             # FETI-DP核心层实现
│   ├── numeric/            # 数值计算层实现
│   ├── tool/               # 基础工具层实现
│   └── main.cpp            # 主程序入口
├── test/                   # 测试代码
├── docs/                   # 项目文档
│   └── header_inclusion_guidelines.md  # 头文件包含原则
├── config/                 # 配置文件
├── mesh/                   # 网格文件
└── lib/                    # 第三方库
```

## 核心设计原则

### 头文件包含原则

本项目采用**分离关注点**的头文件包含原则：

#### 基本原则
- **源文件**: 只包含文件名，不包含路径
- **构建系统**: 负责配置搜索路径
- **模块化**: 每个模块独立管理包含路径

#### 示例

**正确示例 ✅**
```cpp
#include "solver_app.hpp"
#include "logger.hpp"
#include "fe_em_mesh.hpp"
```

**错误示例 ❌**
```cpp
#include "app/solver_app.hpp"
#include "tool/logger.hpp"
#include "fe_em/fe_em_mesh.hpp"
```

#### 优势
- **可维护性**: 文件移动时只需修改CMakeLists.txt
- **可读性**: 源文件更简洁清晰
- **标准化**: 符合现代C++项目最佳实践

详细指南请参考：[头文件包含原则指南](docs/header_inclusion_guidelines.md)

## 模块架构

### 1. 基础工具层 (Tool Layer)
- **日志系统**: 支持多级别日志输出
- **配置文件**: 统一的配置管理
- **工具函数**: 通用工具和辅助函数

### 2. 数值计算层 (Numeric Layer)
- **稀疏矩阵**: 高效的稀疏矩阵运算
- **线性求解器**: 各种线性系统求解方法
- **数值积分**: 有限元数值积分

### 3. 电磁物理层 (FE-EM Layer)
- **网格管理**: 电磁网格的加载和管理
- **物理模型**: 电磁场基本方程和边界条件
- **单元计算**: 有限元单元矩阵计算

### 4. FETI-DP核心层 (FETI-DP Layer)
- **区域分解**: 计算域的区域划分
- **界面条件**: 区域间界面条件处理
- **并行求解**: 分布式并行求解算法

### 5. 应用层 (Application Layer)
- **求解器应用**: 主求解器应用程序
- **用户接口**: 用户交互和配置接口
- **结果输出**: 计算结果的后处理

## 开发指南

### 代码风格
- 遵循现代C++17标准
- 使用Doxygen格式注释
- 保持一致的命名约定

### 构建系统
- 使用CMake作为构建系统
- 支持跨平台编译
- 模块化的依赖管理

### 测试
- 单元测试覆盖核心功能
- 集成测试验证模块间协作
- 性能测试确保算法效率

## 依赖项

### 核心依赖
- **Eigen3**: 线性代数库
- **spdlog**: 高性能日志库

### 可选依赖
- **MPI**: 消息传递接口（并行计算）
- **HDF5**: 科学数据格式（结果输出）

## 贡献指南

1. Fork项目仓库
2. 创建功能分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 创建Pull Request

## 许可证

本项目采用MIT许可证。详见 [LICENSE](LICENSE) 文件。

## 联系方式

- 项目维护者: Poofee
- 邮箱: [请在此处添加联系方式]
- 项目主页: [请在此处添加项目主页]

## 致谢

感谢所有为项目做出贡献的开发者！

---

**注意**: 本项目仍在积极开发中，API可能会发生变化。