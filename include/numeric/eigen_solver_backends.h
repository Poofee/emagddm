/**
 * @file eigen_solver_backends.h
 * @brief Eigen 求解器后端策略实现 - LLT/LDLT/LU 三种分解算法
 * @details 提供3个基于 Eigen 库的求解器后端实现，封装为 SolverBackend 接口：
 *          - EigenLLTBackend：对称正定矩阵 Cholesky 分解 (A = LL^T)
 *          - EigenLDLTBackend：对称不定矩阵 LDL^T 分解 (A = LDL^T)
 *          - EigenLUBackend：通用非对称矩阵 LU 分解 (PA = LU)
 *
 * @par 设计目标：
 * - **零代码重复**：公共逻辑由 UnifiedDirectSolver 处理，本文件仅包含 Eigen 特定逻辑
 * - **性能优化**：利用 Eigen 的稀疏矩阵分解器，支持瞬态复用（一次分解多次求解）
 * - **完整功能**：同时支持实数和复数线性系统
 *
 * @par 适用场景：
 * - **EigenLLTBackend**：静电场、静磁场SPD系统、热传导方程（对称正定矩阵）
 * - **EigenLDLTBackend**：Saddle Point问题、混合有限元、半正定奇异矩阵（对称不定）
 * - **EigenLUBackend**：非对称系统、复数系统、通用稀疏系统（最通用的选择）
 *
 * @par 性能特征：
 * | 后端类型     | 时间复杂度    | 数值稳定性  | 内存占用   | 对称性要求 |
 * |------------|--------------|-----------|----------|----------|
 * | LLT        | O(n³/3)      | 最优       | ~1.5n²    | 必须SPD   |
 * | LDLT       | O(n³/3)      | 良好       | ~1.5n²    | 对称即可  |
 * | LU         | O(n³/3)      | 良好(选主元)| ~2n²     | 无要求    |
 *
 * @see SolverBackend 策略接口基类
 * @see UnifiedDirectSolver 统一调度类
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版)
 */

#pragma once

#include "solver_backend.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <memory>

namespace numeric {

/**
 * @class EigenLLTBackend
 * @brief Eigen Cholesky 分解后端（对称正定矩阵）
 * @details 使用 Eigen::SimplicialLLT 执行稀疏 Cholesky 分解 A = LL^T，
 *          专为对称正定(SPD)矩阵设计，具有最优的数值稳定性和效率。
 *
 * @par 算法特点：
 * - 仅使用矩阵的下三角部分（自动利用对称性）
 * - 无需选主元（SPD矩阵保证数值稳定性）
 * - 使用近似最小度排序(AMD)减少填入
 * - 时间复杂度：分解O(n³/3)，求解O(n²)
 *
 * @par 内部状态：
 * - eigen_matrix_：缓存的实数 Eigen 稀疏矩阵
 * - eigen_matrix_complex_：缓存的复数 Eigen 稀疏矩阵
 * - solver_：实数 SimplicialLLT 分解器实例（缓存L因子）
 * - solver_complex_：复数 SimplicialLLT 分解器实例（缓存L因子）
 *
 * @par 典型应用：
 * @code
 * auto backend = std::make_unique<EigenLLTBackend>();
 * backend->set_matrix(spd_stiffness_matrix);  // 触发Cholesky分解
 * Eigen::VectorXd x = backend->solve_real(rhs_vector);  // 前代+回代
 * @endcode
 */
class EigenLLTBackend : public SolverBackend {
public:
    /**
     * @brief 默认构造函数，初始化内部状态
     */
    EigenLLTBackend() = default;

    /**
     * @brief 析构函数，自动清理资源（RAII）
     */
    ~EigenLLTBackend() override = default;

    void set_matrix(const CsrMatrix<double>& A) override;
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override;
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override;
    void clear() override;

    std::string get_backend_name() const override { return "EigenLLT"; }
    bool supports_complex() const override { return true; }
    bool is_symmetric_only() const override { return true; }
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override { return &eigen_matrix_; }
    const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const override { return &eigen_matrix_complex_; }

private:
    Eigen::SparseMatrix<double> eigen_matrix_;  ///< 缓存的实数系数矩阵
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;  ///< 缓存的复数系数矩阵

    std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> solver_;  ///< 实数Cholesky分解器
    std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>>> solver_complex_;  ///< 复数Cholesky分解器

    bool matrix_set_ = false;       ///< 实数矩阵是否已设置并分解
    bool matrix_set_complex_ = false;  ///< 复数矩阵是否已设置并分解
};

/**
 * @class EigenLDLTBackend
 * @brief Eigen LDL^T 分解后端（对称不定矩阵）
 * @details 使用 Eigen::SimplicialLDLT 执行稀疏 LDL^T 分解 A = LDL^T，
 *          专为对称但未必正定的矩阵设计，可处理接近奇异的情况。
 *
 * @par 算法特点：
 * - 支持对角元出现负值或零（处理不定矩阵）
 * - 自动处理 2x2 主元旋转（提升数值稳定性）
 * - 可配合正则化技术处理奇异矩阵
 * - 时间复杂度：分解O(n³/3)，求解O(n²)
 *
 * @par 与 EigenLLTBackend 的区别：
 * - LLT 要求输入必须 SPD（否则分解失败或结果无意义）
 * - LDLT 允许输入对称但不定（D 的对角元可为负或零）
 * - LDLT 的数值稳定性略低于 LLT（对于 SPD 矩阵）
 *
 * @par 典型应用：
 * @code
 * auto backend = std::make_unique<EigenLDLTBackend>();
 * backend->set_matrix(saddle_point_matrix);  // 对称不定矩阵LDL^T分解
 * Eigen::VectorXd x = backend->solve_real(rhs_vector);
 * @endcode
 */
class EigenLDLTBackend : public SolverBackend {
public:
    EigenLDLTBackend() = default;
    ~EigenLDLTBackend() override = default;

    void set_matrix(const CsrMatrix<double>& A) override;
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override;
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override;
    void clear() override;

    std::string get_backend_name() const override { return "EigenLDLT"; }
    bool supports_complex() const override { return true; }
    bool is_symmetric_only() const override { return true; }
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override { return &eigen_matrix_; }
    const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const override { return &eigen_matrix_complex_; }

private:
    Eigen::SparseMatrix<double> eigen_matrix_;
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;

    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> solver_;  ///< 实数LDLT分解器
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>> solver_complex_;  ///< 复数LDLT分解器

    bool matrix_set_ = false;
    bool matrix_set_complex_ = false;
};

/**
 * @class EigenLUBackend
 * @brief Eigen LU 分解后端（通用非对称矩阵）
 * @details 使用 Eigen::SparseLU 执行带部分选主元的 LU 分解 PA = LU，
 *          是最通用的直接求解器后端，适用于所有类型的方阵。
 *
 * @par 算法特点：
 * - 带列选主元（partial pivoting）保证数值稳定性
 * - 支持完全非对称矩阵（无对称性假设）
 * - 使用 AMD + COLAMD 双重排序优化 fill-in
 * - 时间复杂度：分解O(n³/3)，求解O(n²)
 *
 * @par 内存开销：
 * - LU 分解的 fill-in 通常高于 Cholesky/LDL^T（不对称性导致更多非零元）
 * - 对于大型非对称系统，内存消耗可能是问题
 * - 建议大规模问题考虑迭代法（GMRES/BiCGSTAB + 预条件子）
 *
 * @par 典型应用：
 * @code
 * auto backend = std::make_unique<EigenLUBackend>();
 * backend->set_matrix(general_matrix);  // 通用方阵LU分解
 * Eigen::VectorXd x = backend->solve_real(rhs_vector);
 * @endcode
 */
class EigenLUBackend : public SolverBackend {
public:
    EigenLUBackend() = default;
    ~EigenLUBackend() override = default;

    void set_matrix(const CsrMatrix<double>& A) override;
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override;
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override;
    void clear() override;

    std::string get_backend_name() const override { return "EigenLU"; }
    bool supports_complex() const override { return true; }
    bool is_symmetric_only() const override { return false; }
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override { return &eigen_matrix_; }
    const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const override { return &eigen_matrix_complex_; }

private:
    Eigen::SparseMatrix<double> eigen_matrix_;
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;

    std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> solver_;  ///< 实数LU分解器
    std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>> solver_complex_;  ///< 复数LU分解器

    bool matrix_set_ = false;
    bool matrix_set_complex_ = false;
};

} // namespace numeric
