/**
 * @file superlu_solver_backend.cpp
 * @brief 统一 SuperLU_MT 求解器后端完整实现
 * @details 实现 SuperLUBackend 类的全部方法，封装 SuperLU_MT 库的并行 LU 分解功能。
 *
 * @par 实现策略：
 * - 复用已有的 SuperluContext RAII 封装类
 * - 实现完整的 Eigen (CSC) → SuperLU (CSC) 格式转换
 * - 支持多线程并行 LU 分解（默认使用 1 线程，可配置）
 * - 完整的错误处理和日志输出
 *
 * @par 代码迁移来源：
 * - 原代码位置：SymmetricDirectSolver::decompose_with_superlu() (direct_solvers_superlu.cpp:395)
 *
 * @par 性能特征：
 * - 支持多线程并行消元（PORD 排序 + 并行 pdgstrf）
 * - 列选主元保证数值稳定性
 * - 适合中大规模稀疏系统（1万~50万 DOF）
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.1 (完整集成版)
 */

#include "superlu_solver_backend.h"
#include "direct_solvers.h"     // 使用 SuperluContext
#include "sparse_converter.h"   // SparseConverter 工具类
#include "logger_factory.hpp"
#include <stdexcept>
#include <cstring>

#ifdef HAVE_SUPERLU

namespace numeric {

// ============================================================================
// 实数矩阵设置与分解
// ============================================================================

/**
 * @brief 设置实数稀疏系数矩阵并执行 SuperLU_MT LU 分解
 * @param A CSR格式的实数稀疏系数矩阵
 *
 * @details 统一的实数矩阵分解流程：
 * 1. 将 CSR 矩阵转换为 Eigen CSC 格式并缓存（用于残差计算）
 * 2. 从 Eigen 矩阵提取 CSC 数据（nzval/rowind/colptr）
 * 3. 创建 SuperluContext 实例并初始化
 * 4. 执行并行 LU 分解（pdgstrf_init → pdgstrf → pxgstrf_finalize）
 * 5. 缓存 L、U 因子供后续多次求解复用
 *
 * @note 此方法是瞬态优化的核心：一次分解，多次求解复用
 */
void SuperLUBackend::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("SuperLUBackend::set_matrix - 开始设置矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        const int_t n = static_cast<int_t>(A.rows());

        // 步骤1：将 CSR 矩阵转换为 Eigen CSC 格式并缓存
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();

        const int_t nnz = static_cast<int_t>(eigen_matrix_.nonZeros());

        // 步骤2：从 Eigen CSC 矩阵提取数据数组（SuperLU_MT 使用 CSC 格式输入）
        // 注意：Eigen 的 SparseMatrix 内部就是 CSC 格式，可直接访问内部指针
        double* nzval = static_cast<double*>(malloc(nnz * sizeof(double)));
        int_t* rowind = static_cast<int_t*>(malloc(nnz * sizeof(int_t)));
        int_t* colptr = static_cast<int_t*>((malloc((n + 1) * sizeof(int_t))));

        if (!nzval || !rowind || !colptr) {
            free(nzval); free(rowind); free(colptr);
            throw std::bad_alloc();
        }

        // 复制 Eigen CSC 数据到 SuperLU 数组（Eigen 是 CSC 格式，直接映射）
        memcpy(nzval, eigen_matrix_.valuePtr(), nnz * sizeof(double));
        memcpy(rowind, eigen_matrix_.innerIndexPtr(), nnz * sizeof(int_t));
        memcpy(colptr, eigen_matrix_.outerIndexPtr(), (n + 1) * sizeof(int_t));

        // 步骤3：创建 SuperLU 上下文实例
        superlu_ctx_ = std::make_unique<SuperluContext>();

        // 步骤4：初始化 SuperLU（创建系数矩阵、计算列置换、分配辅助数组）
        if (!superlu_ctx_->initialize(n, nnz, nzval, rowind, colptr)) {
            FEEM_ERROR("SuperLUBackend::set_matrix - SuperLU 初始化失败");
            free(nzval); free(rowind); free(colptr);
            superlu_ctx_.reset();
            eigen_matrix_ = Eigen::SparseMatrix<double>();
            matrix_set_ = false;
            return;
        }

        // 步骤5：执行并行 LU 分解（默认使用 1 线程，可根据需要调整）
        int_t info = superlu_ctx_->factorize(1);

        if (info != 0) {
            std::string error_msg;
            if (info > 0 && info <= n) {
                error_msg = "矩阵可能奇异 (U(" + std::to_string(info) + "," +
                            std::to_string(info) + ")=0)";
            } else if (info < 0) {
                error_msg = "参数错误或内存不足 (info=" + std::to_string(info) + ")";
            } else {
                error_msg = "内存分配失败 (需 " + std::to_string(info - n) + " 字节)";
            }
            FEEM_ERROR("SuperLUBackend::set_matrix - LU分解失败: {}", error_msg);
            superlu_ctx_.reset();
            eigen_matrix_ = Eigen::SparseMatrix<double>();
            matrix_set_ = false;
            return;
        }

        matrix_set_ = true;

        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_INFO("SuperLUBackend::set_matrix 完成, 维度: {}x{}, NNZ={}, 耗时: {:.3f} ms",
                  n, n, nnz, duration_ms);

    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("SuperLUBackend::set_matrix - 内存分配失败: {}", e.what());
        superlu_ctx_.reset();
        eigen_matrix_ = Eigen::SparseMatrix<double>();
        matrix_set_ = false;

    } catch (const std::exception& e) {
        FEEM_ERROR("SuperLUBackend::set_matrix - 异常: {}", e.what());
        superlu_ctx_.reset();
        eigen_matrix_ = Eigen::SparseMatrix<double>();
        matrix_set_ = false;
    }
}

/**
 * @brief 使用缓存的 L、U 因子执行三角回代求解
 * @param b 右端项向量
 * @return Eigen::VectorXd 解向量 x
 *
 * @details 委托给 SuperluContext::solve() 执行实际的三角回代（dgstrs）
 */
Eigen::VectorXd SuperLUBackend::solve_real(const Eigen::VectorXd& b) {
    if (!matrix_set_ || !superlu_ctx_ || !superlu_ctx_->is_factored()) {
        throw std::runtime_error("SuperLUBackend: 矩阵未设置或分解失败，请先调用 set_matrix()");
    }
    if (b.size() != eigen_matrix_.rows()) {
        throw std::invalid_argument("SuperLUBackend: 右端项向量维度不匹配");
    }

    FEEM_DEBUG("SuperLUBackend::solve_real - 开始SuperLU三角回代");

    SolverResult result = superlu_ctx_->solve(b);

    if (result.status != SolverStatus::SUCCESS) {
        throw std::runtime_error("SuperLUBackend: 求解失败 - " + result.error_msg);
    }

    return result.x;
}

// ============================================================================
// 复数支持（当前版本暂不支持）
// ============================================================================

void SuperLUBackend::set_matrix(const CsrMatrix<std::complex<double>>& /*A*/) {
    FEEM_ERROR("SuperLUBackend::set_matrix(复数) - 当前版本不支持复数运算");
    throw std::runtime_error("SuperLUBackend: 不支持复数矩阵（当前仅实现实数版本）");
}

Eigen::VectorXcd SuperLUBackend::solve_complex(const Eigen::VectorXcd& /*b*/) {
    throw std::runtime_error("SuperLUBackend: 不支持复数求解（当前仅实现实数版本）");
}

// ============================================================================
// 资源清理
// ============================================================================

void SuperLUBackend::clear() {
    FEEM_DEBUG("SuperLUBackend::clear - 释放所有资源");

    if (superlu_ctx_) {
        superlu_ctx_->reset();
        superlu_ctx_.reset();
    }

    eigen_matrix_ = Eigen::SparseMatrix<double>();
    matrix_set_ = false;
}

} // namespace numeric

#endif // HAVE_SUPERLU
