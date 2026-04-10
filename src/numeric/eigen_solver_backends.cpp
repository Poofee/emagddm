/**
 * @file eigen_solver_backends.cpp
 * @brief Eigen 求解器后端策略实现（LLT/LDLT/LU 分解）
 * @details 实现3个 Eigen 求解器后端类的全部方法，封装 Eigen 库的稀疏直接求解器。
 *
 * @par 实现要点：
 * - 复用 SparseConverter 工具类进行 CSR ↔ Eigen 格式转换
 * - 利用 Eigen 分解器的瞬态复用机制（一次分解多次求解）
 * - 完整的错误处理和日志输出
 * - 实数/复数双模式支持
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版)
 */

#include "eigen_solver_backends.h"
#include "em_sparse_converter.h"
#include "logger_factory.hpp"
#include <stdexcept>

namespace numeric {

// ============================================================================
// EigenLLTBackend 实现（Cholesky 分解：A = LL^T）
// ============================================================================

void EigenLLTBackend::set_matrix(const CsrMatrix<double>& A) {
    FEEM_DEBUG("EigenLLTBackend::set_matrix - 开始设置实数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();

        solver_ = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>();
        solver_->compute(eigen_matrix_);

        if (solver_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLLTBackend::set_matrix - Cholesky分解失败（矩阵可能非正定或奇异）");
            solver_.reset();
            matrix_set_ = false;
            return;
        }

        matrix_set_ = true;
        FEEM_DEBUG("EigenLLTBackend::set_matrix - Cholesky分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLLTBackend::set_matrix - 异常: {}", e.what());
        solver_.reset();
        matrix_set_ = false;
    }
}

Eigen::VectorXd EigenLLTBackend::solve_real(const Eigen::VectorXd& b) {
    if (!matrix_set_ || !solver_) {
        throw std::runtime_error("EigenLLTBackend: 矩阵未设置或分解失败，请先调用 set_matrix()");
    }
    if (b.size() != eigen_matrix_.rows()) {
        throw std::invalid_argument("EigenLLTBackend: 右端项向量维度不匹配");
    }

    return solver_->solve(b);
}

void EigenLLTBackend::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    FEEM_DEBUG("EigenLLTBackend::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();

        solver_complex_ = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>>>();
        solver_complex_->compute(eigen_matrix_complex_);

        if (solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLLTBackend::set_matrix(复数) - 复数Cholesky分解失败");
            solver_complex_.reset();
            matrix_set_complex_ = false;
            return;
        }

        matrix_set_complex_ = true;
        FEEM_DEBUG("EigenLLTBackend::set_matrix(复数) - 复数Cholesky分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLLTBackend::set_matrix(复数) - 异常: {}", e.what());
        solver_complex_.reset();
        matrix_set_complex_ = false;
    }
}

Eigen::VectorXcd EigenLLTBackend::solve_complex(const Eigen::VectorXcd& b) {
    if (!matrix_set_complex_ || !solver_complex_) {
        throw std::runtime_error("EigenLLTBackend: 复数矩阵未设置或分解失败，请先调用 set_matrix(复数)");
    }
    if (b.size() != eigen_matrix_complex_.rows()) {
        throw std::invalid_argument("EigenLLTBackend: 复数右端项向量维度不匹配");
    }

    return solver_complex_->solve(b);
}

void EigenLLTBackend::clear() {
    FEEM_DEBUG("EigenLLTBackend::clear - 释放所有资源");
    solver_.reset();
    solver_complex_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_ = false;
    matrix_set_complex_ = false;
}

// ============================================================================
// EigenLDLTBackend 实现（LDL^T 分解：A = LDL^T）
// ============================================================================

void EigenLDLTBackend::set_matrix(const CsrMatrix<double>& A) {
    FEEM_DEBUG("EigenLDLTBackend::set_matrix - 开始设置实数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();

        solver_ = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>();
        solver_->compute(eigen_matrix_);

        if (solver_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLDLTBackend::set_matrix - LDL^T分解失败");
            solver_.reset();
            matrix_set_ = false;
            return;
        }

        matrix_set_ = true;
        FEEM_DEBUG("EigenLDLTBackend::set_matrix - LDL^T分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLDLTBackend::set_matrix - 异常: {}", e.what());
        solver_.reset();
        matrix_set_ = false;
    }
}

Eigen::VectorXd EigenLDLTBackend::solve_real(const Eigen::VectorXd& b) {
    if (!matrix_set_ || !solver_) {
        throw std::runtime_error("EigenLDLTBackend: 矩阵未设置或分解失败，请先调用 set_matrix()");
    }
    if (b.size() != eigen_matrix_.rows()) {
        throw std::invalid_argument("EigenLDLTBackend: 右端项向量维度不匹配");
    }

    return solver_->solve(b);
}

void EigenLDLTBackend::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    FEEM_DEBUG("EigenLDLTBackend::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();

        solver_complex_ = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>>();
        solver_complex_->compute(eigen_matrix_complex_);

        if (solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLDLTBackend::set_matrix(复数) - 复数LDL^T分解失败");
            solver_complex_.reset();
            matrix_set_complex_ = false;
            return;
        }

        matrix_set_complex_ = true;
        FEEM_DEBUG("EigenLDLTBackend::set_matrix(复数) - 复数LDL^H分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLDLTBackend::set_matrix(复数) - 异常: {}", e.what());
        solver_complex_.reset();
        matrix_set_complex_ = false;
    }
}

Eigen::VectorXcd EigenLDLTBackend::solve_complex(const Eigen::VectorXcd& b) {
    if (!matrix_set_complex_ || !solver_complex_) {
        throw std::runtime_error("EigenLDLTBackend: 复数矩阵未设置或分解失败，请先调用 set_matrix(复数)");
    }
    if (b.size() != eigen_matrix_complex_.rows()) {
        throw std::invalid_argument("EigenLDLTBackend: 复数右端项向量维度不匹配");
    }

    return solver_complex_->solve(b);
}

void EigenLDLTBackend::clear() {
    FEEM_DEBUG("EigenLDLTBackend::clear - 释放所有资源");
    solver_.reset();
    solver_complex_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_ = false;
    matrix_set_complex_ = false;
}

// ============================================================================
// EigenLUBackend 实现（LU 分解：PA = LU）
// ============================================================================

void EigenLUBackend::set_matrix(const CsrMatrix<double>& A) {
    FEEM_DEBUG("EigenLUBackend::set_matrix - 开始设置实数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();

        solver_ = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
        solver_->compute(eigen_matrix_);

        if (solver_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLUBackend::set_matrix - LU分解失败（矩阵可能奇异）");
            solver_.reset();
            matrix_set_ = false;
            return;
        }

        matrix_set_ = true;
        FEEM_DEBUG("EigenLUBackend::set_matrix - LU分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLUBackend::set_matrix - 异常: {}", e.what());
        solver_.reset();
        matrix_set_ = false;
    }
}

Eigen::VectorXd EigenLUBackend::solve_real(const Eigen::VectorXd& b) {
    if (!matrix_set_ || !solver_) {
        throw std::runtime_error("EigenLUBackend: 矩阵未设置或分解失败，请先调用 set_matrix()");
    }
    if (b.size() != eigen_matrix_.rows()) {
        throw std::invalid_argument("EigenLUBackend: 右端项向量维度不匹配");
    }

    return solver_->solve(b);
}

void EigenLUBackend::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    FEEM_DEBUG("EigenLUBackend::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}", A.rows(), A.cols());

    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();

        solver_complex_ = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>();
        solver_complex_->compute(eigen_matrix_complex_);

        if (solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("EigenLUBackend::set_matrix(复数) - 复数LU分解失败");
            solver_complex_.reset();
            matrix_set_complex_ = false;
            return;
        }

        matrix_set_complex_ = true;
        FEEM_DEBUG("EigenLUBackend::set_matrix(复数) - 复数LU分解成功");
    } catch (const std::exception& e) {
        FEEM_ERROR("EigenLUBackend::set_matrix(复数) - 异常: {}", e.what());
        solver_complex_.reset();
        matrix_set_complex_ = false;
    }
}

Eigen::VectorXcd EigenLUBackend::solve_complex(const Eigen::VectorXcd& b) {
    if (!matrix_set_complex_ || !solver_complex_) {
        throw std::runtime_error("EigenLUBackend: 复数矩阵未设置或分解失败，请先调用 set_matrix(复数)");
    }
    if (b.size() != eigen_matrix_complex_.rows()) {
        throw std::invalid_argument("EigenLUBackend: 复数右端项向量维度不匹配");
    }

    return solver_complex_->solve(b);
}

void EigenLUBackend::clear() {
    FEEM_DEBUG("EigenLUBackend::clear - 释放所有资源");
    solver_.reset();
    solver_complex_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_ = false;
    matrix_set_complex_ = false;
}

} // namespace numeric
