# Tasks

- [ ] Task 1: 扩展现有数据结构接口
  - [ ] SubTask 1.1: 扩展 `SparseMatrixBase`（`sparse_base.hpp`），添加 `merge_duplicates()`、`get_eigen_real()`、`get_eigen_complex()` 接口
  - [ ] SubTask 1.2: 扩展 `CooMatrix<T>`（`coo_matrix.hpp` + `coo_matrix_impl.hpp`），实现新增接口
  - [ ] SubTask 1.3: 扩展 `CsrMatrix<T>`（`csr_matrix.hpp` + `csr_matrix_impl.hpp`），实现新增接口
  - [ ] SubTask 1.4: 扩展 `Vector<T>`（`vector.hpp`），添加 `get_eigen_real()` 和 `get_eigen_complex()` 接口

- [ ] Task 2: 扩展线性求解器基类和结果结构体
  - [ ] SubTask 2.1: 扩展 `EMLinearSolverBase`（`em_linear_solver.h`），添加复数版本纯虚函数
  - [ ] SubTask 2.2: 扩展 `SolverResult`（`em_linear_solver.h`），添加复数字段

- [ ] Task 3: 扩展 SparseConverter 转换工具
  - [ ] SubTask 3.1: 扩展 `SparseConverter`（`em_sparse_converter.h` + `.cpp`），添加复数转换方法

- [ ] Task 4: 扩展直接求解器支持复数
  - [ ] SubTask 4.1: 扩展 `SymmetricDirectSolver`（`em_direct_solvers.h` + `.cpp`），添加复数重载
  - [ ] SubTask 4.2: 扩展 `SymmetricIndefiniteDirectSolver`（`em_direct_solvers.h` + `.cpp`），添加复数重载
  - [ ] SubTask 4.3: 扩展 `GeneralDirectSolver`（`em_direct_solvers.h` + `.cpp`），添加复数重载

- [ ] Task 5: 扩展迭代求解器支持复数
  - [ ] SubTask 5.1: 扩展 `CGSolver`（`em_iterative_solvers.h` + `.cpp`），添加复数重载
  - [ ] SubTask 5.2: 扩展 `BiCGSTABSolver`（`em_iterative_solvers.h` + `.cpp`），添加复数重载

- [ ] Task 6: 扩展工厂类支持复数场景
  - [ ] SubTask 6.1: 扩展 `EMSolverFactory`（`em_solver_factory.h`），支持复数场景

- [ ] Task 7: 创建测试用例
  - [ ] SubTask 7.1: 实数矩阵求解测试（PRISM15 二阶三棱柱标量单元）
  - [ ] SubTask 7.2: 复数矩阵求解测试（PRISM6 混合单元，时谐场 50Hz）
  - [ ] SubTask 7.3: 迭代求解器测试（PYRAMID5_EDGE 棱边单元）

- [ ] Task 8: 更新 CMake 构建配置
  - [ ] SubTask 8.1: 更新 CMakeLists.txt（如需要）

# Task Dependencies
- [Task 2] depends on [Task 1]
- [Task 3] depends on [Task 1]
- [Task 4] depends on [Task 2, Task 3]
- [Task 5] depends on [Task 2, Task 3]
- [Task 6] depends on [Task 4, Task 5]
- [Task 7] depends on [Task 4, Task 5, Task 6]
- [Task 8] depends on [Task 1..7]
