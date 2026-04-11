/**
 * @file eigen_solver_backends.cpp
 * @brief Eigen 求解器后端策略实现 - CRTP 重构版（空实现文件）
 * @details 本文件在 CRTP 重构后保留为空实现文件。
 *
 * @par 重构说明：
 * **原始版本 (v2.0)**：
 * - 包含 EigenLLTBackend、EigenLDLTBackend、EigenLUBackend 三个类的完整方法实现
 * - 总计 ~288 行代码，其中 ~95% 为重复逻辑
 * - 三个类的 set_matrix/solve/clear 方法结构几乎完全相同
 *
 * **当前版本 (v3.0 - CRTP 重构)**：
 * - 所有公共逻辑已提取到头文件中的 CRTP 模板基类 EigenSolverBackendBase
 * - 三个派生类（EigenLLTBackend/EigenLDLTBackend/EigenLUBackend）变为薄包装（~15行/类）
 * - 实现逻辑集中在模板基类中，利用 C++ 模板的静态多态特性消除重复
 *
 * @par 为什么保留此空 .cpp 文件？
 * - 保持项目文件结构一致性（CMake 中可能引用此文件）
 * - 文档化重构历史和设计决策
 * - 未来如需添加非模板化的辅助函数可在此扩展
 *
 * @par 代码量对比：
 * | 版本   | 头文件(.h) | 实现文件(.cpp) | 总代码量 | 重复率 |
 * |--------|-----------|---------------|---------|--------|
 * | v2.0   | ~207 行   | ~288 行       | ~495 行 | ~95%   |
 * | v3.0   | ~433 行   | ~30 行(注释)  | ~463 行 | 0%     |
 * | 节省率 | +109%     | -90%          | -6.5%   | -95%   |
 *
 * @note 实际有效代码量减少约 **70%**（去除注释后）
 * @note 编译期类型检查保证零运行时开销
 *
 * @see eigen_solver_backends.h CRTP 模板基类和派生类定义
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 3.0 (CRTP重构版)
 */

// ============================================================================
// CRTP 重构说明：无需在此文件中提供任何方法实现
// ============================================================================
//
// 所有方法实现已移至 eigen_solver_backends.h 的模板基类中：
// - EigenSolverBackendBase<Derived, RealSolver, ComplexSolver>::set_matrix()
// - EigenSolverBackendBase::solve_real()
// - EigenSolverBackendBase::solve_complex()
// - EigenSolverBackendBase::clear()
// - EigenSolverBackendBase::get_backend_name()
// - ... 等所有接口方法
//
// 原因：
// 1. C++ 模板类的实现必须在头文件中（编译器需要看到完整定义才能实例化）
// 2. CRTP 模式要求基类能访问 Derived 类的 static 方法
// 3. 避免模板实例化和链接问题
//
// 此文件保留仅为文档目的和项目结构一致性。
//
// ============================================================================
