# 工业级电磁场有限元线性求解器模块 Spec

## Why
当前电磁场有限元系统缺少统一的线性求解器模块，无法支撑从标量场到矢量场、从静态场到瞬态场的全场景求解需求。需要构建工业级线性求解器系统，与现有装配模块、DOF管理模块100%无缝对接。

## What Changes
- **新增** 线性求解器抽象基类（EMLinearSolverBase），定义统一接口规范
- **新增** CSR矩阵与Eigen稀疏矩阵转换工具类（SparseConverter）
- **新增** 3个直接求解器实现（**支持多后端策略模式**：Eigen/SuperLU/MUMPS可选）
- **新增** DirectBackendType 枚举和 DirectBackendManager 后端管理器
- **新增** 2个迭代求解器实现：CGSolver、BiCGSTABSolver
- **新增** 代数多重网格（AMG）预条件子框架（含标量AMG实现和H(curl) AMG接口预留）
- **新增** 求解器工厂类（EMSolverFactory），支持场景自动适配
- **新增** 3个完整测试用例覆盖全电磁场景
- **零修改** 所有已有代码保持不变，仅扩展numeric命名空间

## Impact
- Affected specs: 数值计算层核心能力扩展
- Affected code: include/numeric/、src/numeric/ 目录新增文件
- External dependencies: 
  - **必需**: Eigen 3.4+（已集成）
  - **可选**: SuperLU（条件编译 #ifdef HAVE_SUPERLU）
  - **可选**: MUMPS（条件编译 #ifdef HAVE_MUMPS，需 MPI）

## ADDED Requirements

### Requirement: 求解器抽象基类接口
系统 SHALL 提供统一线性求解器抽象基类 EMLinearSolverBase，包含以下核心接口：
1. set_matrix() - 设置系数矩阵（支持瞬态多步复用）
2. solve() - 求解方程组 Ax = b
3. get_solver_name() - 获取求解器标识
4. clear() - 释放矩阵分解资源

#### Scenario: 基础接口调用流程
- **WHEN** 用户创建求解器实例并调用 set_matrix(A) 后多次调用 solve(b)
- **THEN** 首次solve执行矩阵分解并缓存，后续solve复用分解结果提升性能
- **THEN** 返回 SolverResult 包含解向量、状态码、迭代次数、残差范数、耗时统计

### Requirement: CSR-Eigen矩阵转换工具
系统 SHALL 提供 SparseConverter 工具类，实现 numeric::CsrMatrix<double> 与 Eigen::SparseMatrix<double> 的双向转换：
1. to_eigen() - 将自定义CSR转为Eigen格式
2. from_eigen() - 将Eigen格式转回自定义CSR
3. 转换过程保证数据一致性，支持维度校验和性能优化

#### Scenario: 矩阵转换正确性验证
- **WHEN** 用户将 CsrMatrix<double> 通过 to_eigen 转换后再通过 from_eigen 转回
- **THEN** 原始矩阵与往返转换后矩阵的行数、列数、非零元数量完全一致
- **THEN** 对应位置的元素值误差小于机器精度（1e-15）

### Requirement: 直接求解器实现（工业级鲁棒性 + 多后端支持）
系统 SHALL 实现3个直接求解器派生类，**每个求解器支持可插拔的底层算法后端选择**：

#### 架构设计：策略模式（Strategy Pattern）
- **核心思想**: 每个直接求解器类内部使用策略模式，通过枚举参数选择底层算法实现
- **后端类型枚举 DirectBackendType**:
  ```cpp
  enum class DirectBackendType {
      EIGEN,       // Eigen内置稀疏求解器（默认，无需额外依赖）
      SUPERLU,     // SuperLU（开源高性能直接求解器）
      MUMPS        // MUMPS（并行分布式直接求解器）
  };
  ```
- **接口统一**: 无论选择哪种后端，对外接口完全一致（set_matrix/solve/clear）
- **运行时切换**: 支持构造时或通过 set_backend() 方法动态切换后端
- **优雅降级**: 如果请求的后端库未编译/未安装，自动回退到 Eigen 后端并输出 WARNING 日志

#### 3.1 SymmetricDirectSolver（对称正定矩阵求解器）
- **适配场景**: 二维/三维静电场/静磁场（标量位）、瞬态场有效刚度矩阵
- **支持的后端算法**:
  - **EIGEN 后端**: Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> （默认）
    - 符号分解 + LDL^T 分解（Cholesky变体）
    - 优先支持 Intel MKL Pardiso（条件编译 #ifdef EIGEN_USE_MKL_ALL）
  - **SUPERLU 后端**: SuperLU_DIST 或 SuperLU_MT
    - 对称矩阵优化存储（仅存上三角）
    - 支持多线程并行分解
  - **MUMPS 后端**: MUMPS (MUltifrontal Massively Parallel Sparse direct Solver)
    - 分布式内存并行求解（需 MPI 支持）
    - 对称正定矩阵专用对称模式 (ICNTL(28)=1, ICNTL(29)=2)
- **核心功能（与后端无关）**:
  - 符号分解仅执行一次，数值分解可复用（瞬态多步优化）
  - 自动检测矩阵正定性，失败时返回 NUMERICAL_ERROR 状态
  - 支持矩阵对称性校验（容差可配置，默认 1e-10）
  - 后端性能对比日志输出（各后端的分解时间、求解时间）

#### 3.2 SymmetricIndefiniteDirectSolver（对称半正定/不定矩阵求解器）
- **适配场景**: 三维棱边单元静磁场（H(curl)半正定）、A-V混合涡流场（不定）
- **支持的后端算法**:
  - **EIGEN 后端**: Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> （默认）
    - 稳定的 LDL^T 分解，处理不定对角元
    - 自动处理零主元（添加微小扰动保证数值稳定）
  - **SUPERLU 后端**: SuperLU 不定对称模式
    - 完全选主元的 LU 分解（保证稳定性）
    - 支持对称存储优化
  - **MUMPS 后端**: MUMPS 非对称通用模式 + 对称预处理
    - 处理奇异矩阵的鲁棒性强
    - 支持秩检测和零空间计算
- **核心功能（与后端无关）**:
  - 自动处理矩阵奇异性（棱边单元零梯度空间）
  - 奇异矩阵检测与零空间维度报告
  - 正则化参数配置（regularization_epsilon，默认 1e-14）
  - 符号分解复用机制

#### 3.3 GeneralDirectSolver（通用非对称矩阵求解器）
- **适配场景**: 非对称电磁问题、谐波场、复杂耦合场
- **支持的后端算法**:
  - **EIGEN 后端**: Eigen::SparseLU<Eigen::SparseMatrix<double>> （默认）
    - 列主序 LU 分解，部分选主元策略（阈值选主元）
    - 支持任意非对称方阵的鲁棒求解
  - **SUPERLU 后端**: SuperLU 完全选主元模式
    - 全选主元保证数值稳定性
    - 适合病态矩阵（条件数 < 1e15）
  - **MUMPS 后端**: MUMPS 非对称并行模式
    - 分布式大规模非对称求解
    - 支持复杂非对称结构（如涡流场的复数矩阵）
- **核心功能（与后端无关）**:
  - 支持任意非对称方阵的LU分解
  - 鲁棒性优先，适配病态矩阵（条件数<1e12）
  - 病态矩阵警告机制（条件数估算 > 1e10 时输出 WARNING 日志）

#### 多后端管理辅助类
系统 SHALL 提供 DirectBackendManager 辅助类：
```cpp
class DirectBackendManager {
public:
    static bool isBackendAvailable(DirectBackendType type);
    static std::string getBackendName(DirectBackendType type);
    static std::vector<DirectBackendType> getAvailableBackends();
};
```
- 运行时检测哪些后端可用（基于编译期宏定义 #ifdef HAVE_SUPERLU / #ifdef HAVE_MUMPS）
- 提供友好的后端名称字符串用于日志输出
- 查询当前系统所有已编译启用的后端列表

#### Scenario: 直接求解器多后端工作流
- **WHEN** 用户创建 SymmetricDirectSolver 并设置 backend_type = MUMPS
- **THEN** 如果 MUMPS 已编译启用，使用 MUMPS 并行后端进行分布式求解
- **THEN** 如果 MUMPS 未启用，自动回退到 EIGEN 后端并输出 WARNING: "MUMPS backend not available, fallback to Eigen"
- **WHEN** 用户在运行时调用 set_backend(SUPERLU) 切换后端
- **THEN** 清空已有分解结果，使用新后端重新执行符号+数值分解
- **WHEN** 用户使用任一后端求解15×15对称正定刚度矩阵
- **THEN** 所有后端均返回 SUCCESS 状态，残差范数 < 1e-10，结果误差 < 1e-6

### Requirement: 迭代求解器实现（高性能大规模问题）
系统 SHALL 实现2个迭代求解器派生类及预条件子系统：

#### 4.1 CGSolver（共轭梯度求解器）
- 底层算法：Eigen::ConjugateGradient
- 适配场景：对称正定/半正定矩阵（标量场、棱边单元静磁场）
- 核心功能：
  - 可配置收敛阈值（tolerance）、最大迭代次数（max_iterations）
  - 支持对角预条件子（Jacobi）、ILU0预条件子
  - 迭代过程实时残差监控，返回收敛曲线数据
  - 自动检测发散情况，返回 DIVERGED 或 MAX_ITER_REACHED 状态

#### 4.2 BiCGSTABSolver（稳定双共轭梯度求解器）
- 底层算法：Eigen::BiCGSTAB
- 适配场景：对称不定矩阵、非对称矩阵（A-V混合涡流场、谐波场）
- 核心功能：
  - 可配置收敛阈值、最大迭代次数
  - 支持ILU0预条件子
  - 适配非对称矩阵的稳定迭代，避免标准BiCG的不稳定性

#### 4.3 AMG预条件子框架（代数多重网格）
- 标量AMG预条件子（ScalarAMG）：
  - 适配标量场SPD矩阵的大规模求解
  - 实现经典Ruge-Stüben粗化策略或聚合式AMG
  - 完全兼容Eigen迭代求解器的预条件子接口规范
- H(curl)专用AMG接口预留（HiptmairAMG）：
  - 适配三维棱边单元H(curl)空间大规模问题
  - 预留接口签名，等待后续完整实现
  - 支持高对比度材料电磁问题的收敛加速

#### Scenario: 迭代求解器收敛性测试
- **WHEN** 用户使用 CGSolver 配合Jacobi预条件子求解1000×1000 SPD矩阵
- **THEN** 在合理迭代次数内（<500次）收敛至 tolerance=1e-8
- **THEN** 返回正确的迭代次数和最终残差范数

### Requirement: 求解器工厂模式
系统 SHALL 提供 EMSolverFactory 工厂类，根据物理场景枚举自动创建对应求解器实例：
- SYMMETRIC_DIRECT → SymmetricDirectSolver
- SYMMETRIC_INDEFINITE → SymmetricIndefiniteDirectSolver
- GENERAL_DIRECT → GeneralDirectSolver
- CG → CGSolver
- BICGSTAB → BiCGSTABSolver

#### Scenario: 工厂模式快速创建
- **WHEN** 用户调用 EMSolverFactory::create_solver(EMSolverFactory::SolverType::CG)
- **THEN** 返回 std::unique_ptr<EMLinearSolverBase> 指向 CGSolver 实例
- **THEN** 返回的指针可直接调用 set_matrix() 和 solve() 接口

## MODIFIED Requirements
无（本模块为全新增量开发，不修改任何已有代码）

## REMOVED Requirements
无

## 技术约束
1. **语言标准**: C++17及以上，禁止使用C++20特性
2. **核心依赖**: STL + Eigen 3.4+（必需）
3. **可选依赖**:
   - SuperLU（条件编译 #ifdef HAVE_SUPERLU，高性能直接求解器）
   - MUMPS（条件编译 #ifdef HAVE_MUMPS，并行分布式求解器，需 MPI）
   - 所有可选依赖未安装时自动回退到 Eigen 后端
4. **内存管理**: 禁止裸指针管理动态资源，优先使用智能指针
5. **线程安全**: 无全局变量，所有状态封装在实例对象内
6. **日志规范**: 使用 FEEM_INFO/FEEM_DEBUG/FEEM_ERROR 等宏输出日志
7. **头文件规范**: 仅包含文件名，不包含路径（由CMakeLists.txt管理路径）
8. **命名规范**: 类名PascalCase、函数camelCase、变量snake_case、常量UPPER_SNAKE_CASE
9. **注释要求**: 关键逻辑必须注释，使用Doxygen风格文档注释
10. **输入校验**: 所有公开接口必须进行参数合法性检查
11. **错误处理**: 使用 SolverStatus 枚举明确反馈错误类型，禁止静默失败
12. **后端可插拔性**: 直接求解器必须支持运行时切换后端（策略模式），接口完全统一
