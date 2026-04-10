# Checklist

- [ ] SparseMatrixBase 扩展正确
  - [ ] 添加了 `merge_duplicates()` 纯虚函数
  - [ ] 添加了 `get_eigen_real()` 纯虚函数
  - [ ] 添加了 `get_eigen_complex()` 纯虚函数
  - [ ] 原有接口保持不变

- [ ] CooMatrix 扩展正确
  - [ ] 实现了 `merge_duplicates()` 方法
  - [ ] 实现了 `get_eigen_real()` 方法
  - [ ] 实现了 `get_eigen_complex()` 方法
  - [ ] 原有接口保持不变

- [ ] CsrMatrix 扩展正确
  - [ ] 实现了 `merge_duplicates()` 方法
  - [ ] 实现了 `get_eigen_real()` 方法
  - [ ] 实现了 `get_eigen_complex()` 方法
  - [ ] 原有接口保持不变

- [ ] Vector 扩展正确
  - [ ] 实现了 `get_eigen_real()` 方法
  - [ ] 实现了 `get_eigen_complex()` 方法
  - [ ] 原有接口保持不变

- [ ] EMLinearSolverBase 扩展正确
  - [ ] 添加了复数版本 `set_matrix(const CsrMatrix<std::complex<double>>& A)` 纯虚函数
  - [ ] 添加了复数版本 `solve(const Eigen::VectorXcd& b)` 纯虚函数
  - [ ] 原有实数接口保持不变

- [ ] SolverResult 扩展正确
  - [ ] 添加了 `Eigen::VectorXcd x_complex` 字段
  - [ ] 原有字段保持不变

- [ ] SparseConverter 扩展正确
  - [ ] 添加了 `to_eigen_complex()` 静态方法
  - [ ] 添加了 `from_eigen_complex()` 静态方法
  - [ ] 原有方法保持不变

- [ ] SymmetricDirectSolver 复数支持正确
  - [ ] 添加了复数 `set_matrix` 重载
  - [ ] 添加了复数 `solve` 重载
  - [ ] 内部维护复数 Eigen 求解器实例
  - [ ] 原有实数版本保持不变

- [ ] SymmetricIndefiniteDirectSolver 复数支持正确
  - [ ] 添加了复数 `set_matrix` 和 `solve` 重载
  - [ ] 内部维护复数 Eigen 求解器实例
  - [ ] 原有实数版本保持不变

- [ ] GeneralDirectSolver 复数支持正确
  - [ ] 添加了复数 `set_matrix` 和 `solve` 重载
  - [ ] 内部维护复数 Eigen 求解器实例
  - [ ] 原有实数版本保持不变

- [ ] CGSolver 复数支持正确
  - [ ] 添加了复数 `set_matrix` 和 `solve` 重载
  - [ ] 内部维护复数 Eigen 迭代求解器实例
  - [ ] 原有实数版本保持不变

- [ ] BiCGSTABSolver 复数支持正确
  - [ ] 添加了复数 `set_matrix` 和 `solve` 重载
  - [ ] 内部维护复数 Eigen 迭代求解器实例
  - [ ] 原有实数版本保持不变

- [ ] EMSolverFactory 扩展正确
  - [ ] 支持创建复数场景的求解器（可选）

- [ ] 测试用例通过
  - [ ] 实数矩阵求解测试通过
  - [ ] 复数矩阵求解测试通过
  - [ ] 迭代求解器测试通过

- [ ] 向后兼容性保证
  - [ ] 所有原有接口保持不变
  - [ ] 现有代码无需修改即可使用
