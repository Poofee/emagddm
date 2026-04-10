/**
 * @file mumps_solver_backend.h
 * @brief 统一 MUMPS 求解器后端 - 支持所有矩阵类型的单一实现
 * @details 提供 MUMPSBackend 类，统一封装 MUMPS 直接求解器的实数和复数功能。
 *          **核心创新**：通过构造函数的 MatrixType 参数区分矩阵类型（SPD/SID/非对称），
 *          将原来3个求解器类中的 MUMPS 代码（70%+ 重复）合并为1个类！
 *
 * @par 重构前后对比：
 * **重构前**（代码重复率 70%+）：
 * ```
 * SymmetricDirectSolver::decompose_with_mumps()  // sym=1
 * SymmetricIndefiniteDirectSolver::decompose_with_mumps()  // sym=2
 * GeneralDirectSolver::decompose_with_mumps()  // sym=0
 * + 各自的 solve_with_mumps() 完全相同
 * + 复数版本又重复一遍
 * ```
 *
 * **重构后**（零重复）：
 * ```
 * MUMPSBackend(MatrixType::SYMMETRIC_POSITIVE_DEFINITE)  // sym=1
 * MUMPSBackend(MatrixType::SYMMETRIC_INDEFINITE)        // sym=2
 * MUMPSBackend(MatrixType::UNSYMMETRIC)                 // sym=0
 * 所有公共逻辑只写一次！
 * ```
 *
 * @par 矩阵类型映射：
 * | MatrixType 枚举值          | MUMPS sym 参数 | 对应原求解器           | 适用场景               |
 * |--------------------------|---------------|---------------------|----------------------|
 * | SYMMETRIC_POSITIVE_DEFINITE | 1             | SymmetricDirectSolver | 静电场、静磁场SPD系统 |
 * | SYMMETRIC_INDEFINITE       | 2             | SymmetricIndefiniteDirectSolver | Saddle Point问题 |
 * | UNSYMMETRIC                | 0             | GeneralDirectSolver   | 非对称通用系统       |
 *
 * @par 内部架构：
 * - 复用已有的 MumsContext（实数MUMPS RAII封装）
 * - 复用已有的 ComplexMumpsContext（复数MUMPS RAII封装）
 * - 统一管理 CSR→COO 转换和分解流程
 * - 同时支持实数和复数模式（互不干扰）
 *
 * @par 性能优势：
 * - MPI并行分布式计算能力（适合大规模问题 >10万DOF）
 * - 优秀的内存管理和负载均衡
 * - 瞬态复用：set_matrix() 分解一次，solve() 多次复用
 *
 * @see SolverBackend 策略接口基类
 * @see MumpsContext 实数MUMPS RAII封装
 * @see ComplexMumpsContext 复数MUMPS RAII封装
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

#ifdef HAVE_MUMPS

// 前向声明（完整定义在 .cpp 文件中通过 #include "em_direct_solvers.h" 获取）
class MumpsContext;
class ComplexMumpsContext;

/**
 * @class MUMPSBackend
 * @brief 统一 MUMPS 直接求解器后端（支持所有矩阵类型）
 * @details 封装 MUMPS 库的全部功能（实数+复数），通过 MatrixType 枚举区分矩阵对称性，
 *          完全消除原三个求解器类中的 MUMPS 代码重复。
 *
 * @par 设计原则：
 * 1. **零代码重复**：所有 MUMPS 接口调用只写一次
 * 2. **类型安全**：通过枚举参数而非魔术数字指定矩阵类型
 * 3. **资源安全**：RAII 管理 MUMPS 生命周期（依赖 MumpsContext/ComplexMumpsContext）
 * 4. **双模式支持**：同时支持实数和复数线性系统，状态独立管理
 *
 * @par 典型使用流程：
 * @code
 * // 创建 MUMPS 后端（指定矩阵类型为对称正定）
 * auto backend = std::make_unique<MUMPSBackend>(MUMPSBackend::MatrixType::SYMMETRIC_POSITIVE_DEFINITE);
 *
 * // 设置系数矩阵（自动完成初始化、CSR→COO转换、分析+分解）
 * backend->set_matrix(spd_matrix);
 *
 * // 多次求解（复用缓存的分解结果）
 * for (auto& rhs : rhs_list) {
 *     Eigen::VectorXd x = backend->solve_real(rhs);
 * }
 *
 * // 清理资源（或依赖析构函数自动清理）
 * backend->clear();
 * @endcode
 */
class MUMPSBackend : public SolverBackend {
public:
    /**
     * @enum MatrixType
     * @brief MUMPS 矩阵对称性类型枚举
     * @details 对应原三个直接求解器的矩阵假设，映射到 MUMPS 的 sym 参数
     */
    enum class MatrixType {
        SYMMETRIC_POSITIVE_DEFINITE,  ///< 对称正定(SPD)矩阵，对应 SymmetricDirectSolver (sym=1)
        SYMMETRIC_INDEFINITE,         ///< 对称不定(SID)矩阵，对应 SymmetricIndefiniteDirectSolver (sym=2)
        UNSYMMETRIC                   ///< 非对称(一般)矩阵，对应 GeneralDirectSolver (sym=0)
    };

    /**
     * @brief 构造函数：指定矩阵对称性类型
     * @param matrix_type 矩阵类型枚举值，决定 MUMPS 的 sym 参数
     */
    explicit MUMPSBackend(MatrixType matrix_type);

    ~MUMPSBackend() override = default;

    void set_matrix(const CsrMatrix<double>& A) override;
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override;
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override;
    void clear() override;

    std::string get_backend_name() const override;
    bool supports_complex() const override { return true; }
    bool is_symmetric_only() const override;
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override { return &eigen_matrix_; }
    const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const override { return &eigen_matrix_complex_; }

private:
    MatrixType matrix_type_;
    std::unique_ptr<MumpsContext> mumps_ctx_;
    std::unique_ptr<ComplexMumpsContext> mumps_ctx_complex_;
    bool matrix_set_ = false;
    bool matrix_set_complex_ = false;
    Eigen::SparseMatrix<double> eigen_matrix_;
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;

    int get_mumps_sym() const;
};

#endif // HAVE_MUMPS

} // namespace numeric
