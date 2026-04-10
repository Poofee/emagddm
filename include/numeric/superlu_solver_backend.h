/**
 * @file superlu_solver_backend.h
 * @brief 统一 SuperLU_MT 求解器后端 - 高性能并行 LU 分解
 * @details 提供 SuperLUBackend 类，封装 SuperLU_MT 库的并行 LU 分解功能。
 *          与 MUMPSBackend 类似的设计理念，通过统一的接口消除代码重复。
 *
 * @par 适用场景：
 * - 中大规模稀疏线性系统（1万~50万 DOF）
 * - 需要多线程并行加速的 LU 分解
 * - SuperLU_MT 的列选主元 LU 分解比 Eigen::SparseLU 更稳定（对于病态矩阵）
 *
 * @par 性能特征：
 * - 支持多线程并行（PORD 排序 + 并行消去）
 * - 列选主元保证数值稳定性
 * - 适合共享内存多核系统
 *
 * @note 此类仅在 HAVE_SUPERLU 编译选项启用时可用
 *
 * @see SolverBackend 策略接口基类
 * @see SuperluContext SuperLU RAII 封装
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.1 (完整集成版)
 */

#pragma once

#include "solver_backend.h"
#include <memory>
#include <string>

namespace numeric {

#ifdef HAVE_SUPERLU

// 前向声明（完整定义在 .cpp 文件中通过 #include "em_direct_solvers.h" 获取）
class SuperluContext;

/**
 * @class SuperLUBackend
 * @brief 统一 SuperLU_MT 求解器后端（并行 LU 分解）
 * @details 封装 SuperLU_MT 库的并行 LU 分解功能，提供与 MUMPSBackend 类似的统一接口。
 *          目前仅实现实数版本（SuperLU_MT 的复数支持较复杂，可后续扩展）。
 *
 * @par 设计特点：
 * - 复用已有的 SuperluContext RAII 封装
 * - 统一的错误处理和日志输出
 * - 支持瞬态复用（一次分解多次求解）
 *
 * @par 典型使用流程：
 * @code
 * auto backend = std::make_unique<SuperLUBackend>();
 * backend->set_matrix(general_matrix);  // 触发并行 LU 分解
 * Eigen::VectorXd x = backend->solve_real(rhs);
 * @endcode
 */
class SuperLUBackend : public SolverBackend {
public:
    /**
     * @brief 默认构造函数
     */
    SuperLUBackend() = default;

    /**
     * @brief 虚析构函数（RAII 自动清理）
     */
    ~SuperLUBackend() override = default;

    // ===== 核心操作接口 =====

    void set_matrix(const CsrMatrix<double>& A) override;
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override;
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override;
    void clear() override;

    // ===== 属性查询接口 =====

    std::string get_backend_name() const override { return "SuperLU_MT"; }
    bool supports_complex() const override { return false; }  ///< 当前仅支持实数
    bool is_symmetric_only() const override { return false; }
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override { return &eigen_matrix_; }

private:
    std::unique_ptr<SuperluContext> superlu_ctx_;  ///< SuperLU MT RAII 封装

    bool matrix_set_ = false;  ///< 实数矩阵是否已设置并分解
    Eigen::SparseMatrix<double> eigen_matrix_;  ///< 缓存的实数 Eigen 矩阵
};

#endif // HAVE_SUPERLU

} // namespace numeric
