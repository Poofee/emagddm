# 工业级实数/复数统一架构线性求解器模块 Spec

## Why
当前电磁场有限元求解器已有完整的求解器体系（直接求解器、迭代求解器、工厂类），但**仅支持实数矩阵**（`CsrMatrix<double>`、`Eigen::VectorXd`）。需要扩展现有代码以同时支持复数矩阵（时谐场场景），实现上层接口统一、底层模板特化的高性能架构。

## What Changes
- 扩展 `SparseMatrixBase` 抽象基类，添加 `merge_duplicates()`、`get_eigen_real()`、`get_eigen_complex()` 接口
- 扩展 `CooMatrix<T>` 和 `CsrMatrix<T>` 模板类，实现上述新接口
- 扩展 `Vector<T>` 模板类，添加 `get_eigen_real()` 和 `get_eigen_complex()` 接口
- 扩展 `EMLinearSolverBase` 抽象基类，添加复数版本的纯虚函数接口
- 扩展 `SolverResult` 结构体，添加复数解向量字段
- 扩展 `SparseConverter` 转换工具类，添加复数矩阵转换方法
- **扩展三个直接求解器**（`em_direct_solvers.h/.cpp`），添加复数版本的重载方法
- **扩展两个迭代求解器**（`em_iterative_solvers.h/.cpp`），添加复数版本的重载方法
- **扩展工厂类**（`em_solver_factory.h`），支持创建复数场景求解器
- 新增测试用例

## Impact
- Affected specs: 实数/复数统一求解能力
- Affected code:
  - `include/numeric/sparse_base.hpp` - 扩展抽象基类
  - `include/numeric/coo_matrix.hpp` + `coo_matrix_impl.hpp` - 扩展COO矩阵
  - `include/numeric/csr_matrix.hpp` + `csr_matrix_impl.hpp` - 扩展CSR矩阵
  - `include/numeric/vector.hpp` - 扩展向量类
  - `include/numeric/em_linear_solver.h` - 扩展求解器基类
  - `include/numeric/em_sparse_converter.h` + `.cpp` - 扩展转换工具
  - `include/numeric/em_direct_solvers.h` + `.cpp` - **扩展直接求解器（已有）**
  - `include/numeric/em_iterative_solvers.h` + `.cpp` - **扩展迭代求解器（已有）**
  - `include/numeric/em_solver_factory.h` - **扩展工厂类（已有）**

## ADDED Requirements

### Requirement: 扩展稀疏矩阵和向量接口
系统 SHALL 在现有数据结构基础上扩展求解器所需接口。

#### Scenario: SparseMatrixBase 扩展
- **WHEN** 扩展 `SparseMatrixBase`
- **THEN** 添加 `merge_duplicates()`、`get_eigen_real()`、`get_eigen_complex()` 纯虚函数

#### Scenario: CooMatrix/CsrMatrix/Vector 扩展
- **WHEN** 扩展模板类
- **THEN** 实现上述新接口，保持原有接口不变

### Requirement: 扩展求解器基类支持复数
系统 SHALL 扩展 `EMLinearSolverBase`，添加复数版本的纯虚函数接口。

#### Scenario: 复数 set_matrix 和 solve 接口
- **WHEN** 扩展 `EMLinearSolverBase`
- **THEN** 添加 `set_matrix(const CsrMatrix<std::complex<double>>& A)` 纯虚函数
- **AND** 添加 `solve(const Eigen::VectorXcd& b)` 纯虚函数

#### Scenario: SolverResult 扩展
- **WHEN** 扩展 `SolverResult`
- **THEN** 添加 `Eigen::VectorXcd x_complex` 字段

### Requirement: 扩展 SparseConverter 支持复数
系统 SHALL 扩展 `SparseConverter`，添加复数矩阵转换方法。

#### Scenario: 复数转换方法
- **WHEN** 扩展 `SparseConverter`
- **THEN** 添加 `to_eigen_complex()` 和 `from_eigen_complex()` 静态方法

### Requirement: 扩展直接求解器支持复数
系统 SHALL 扩展现有三个直接求解器，添加复数版本重载。

#### Scenario: 三个直接求解器复数支持
- **WHEN** 扩展 `SymmetricDirectSolver`、`SymmetricIndefiniteDirectSolver`、`GeneralDirectSolver`
- **THEN** 各自添加复数版本的 `set_matrix` 和 `solve` 重载方法
- **AND** 内部维护独立的复数 Eigen 求解器实例

### Requirement: 扩展迭代求解器支持复数
系统 SHALL 扩展现有两个迭代求解器，添加复数版本重载。

#### Scenario: CGSolver 和 BiCGSTABSolver 复数支持
- **WHEN** 扩展 `CGSolver` 和 `BiCGSTABSolver`
- **THEN** 各自添加复数版本的 `set_matrix` 和 `solve` 重载方法
- **AND** 内部维护独立的复数 Eigen 迭代求解器实例

### Requirement: 扩展工厂类支持复数场景
系统 SHALL 扩展 `EMSolverFactory`，支持创建复数场景的求解器。

#### Scenario: 工厂类复数支持
- **WHEN** 扩展 `EMSolverFactory`
- **THEN** 可选地支持指定是否为复数矩阵场景
- **AND** 创建的求解器自动适配实数或复数接口

### Requirement: 测试用例
系统 SHALL 提供完整的测试用例。

## MODIFIED Requirements

### Requirement: 所有现有求解器类需添加复数重载
详见 ADDED Requirements 中的具体描述。

## REMOVED Requirements
无移除需求。
