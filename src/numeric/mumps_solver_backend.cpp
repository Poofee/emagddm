/**
 * @file mumps_solver_backend.cpp
 * @brief 统一 MUMPS 求解器后端完整实现（零重复版本）
 * @details 实现 MUMPSBackend 类的全部方法，将原3个求解器类中的 MUMPS 代码
 *          （约 70% 重复率）合并为单一实现。
 *
 * @par 实现策略：
 * - 复用已有的 MumpsContext 和 ComplexMumpsContext RAII 封装
 * - 通过 matrix_type_ 成员变量区分矩阵类型（唯一区别点：sym 参数）
 * - 所有公共逻辑（输入校验、CSR→COO转换、错误处理）只写一次
 * - 完整支持实数和复数双模式
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.1 (完整集成版)
 */

#include "mumps_solver_backend.h"
#include "em_solver_backends.hpp"  // 必须：定义 HAVE_MUMPS 宏 + MUMPS 头文件
#include "direct_solvers.h"     // 使用 MumpsContext / ComplexMumpsContext / extract_complex
#include "sparse_converter.h"   // SparseConverter 工具类
#include "logger_factory.hpp"
#include <stdexcept>

namespace numeric {
#ifdef HAVE_MUMPS

// ============================================================================
// 构造函数与属性查询
// ============================================================================

/**
 * @brief 构造函数：指定 MUMPS 矩阵对称性类型
 */
MUMPSBackend::MUMPSBackend(MatrixType matrix_type)
    : matrix_type_(matrix_type) {

    FEEM_INFO("MUMPSBackend 构造完成, 矩阵类型: {}",
              matrix_type_ == MatrixType::SYMMETRIC_POSITIVE_DEFINITE ? "对称正定(SPD)" :
              matrix_type_ == MatrixType::SYMMETRIC_INDEFINITE ? "对称不定(SID)" : "非对称(General)");
}

std::string MUMPSBackend::get_backend_name() const {
    switch (matrix_type_) {
        case MatrixType::SYMMETRIC_POSITIVE_DEFINITE:
            return "MUMPS_SPD";
        case MatrixType::SYMMETRIC_INDEFINITE:
            return "MUMPS_SID";
        case MatrixType::UNSYMMETRIC:
            return "MUMPS_LU";
        default:
            return "MUMPS_Unknown";
    }
}

bool MUMPSBackend::is_symmetric_only() const {
    return matrix_type_ != MatrixType::UNSYMMETRIC;
}

int MUMPSBackend::get_mumps_sym() const {
    switch (matrix_type_) {
        case MatrixType::SYMMETRIC_POSITIVE_DEFINITE:
            return 1;
        case MatrixType::SYMMETRIC_INDEFINITE:
            return 2;
        case MatrixType::UNSYMMETRIC:
            return 0;
        default:
            return 0;
    }
}

// ============================================================================
// 实数矩阵设置与分解
// ============================================================================

void MUMPSBackend::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("MUMPSBackend::set_matrix(实数) - 开始设置矩阵, 维度: {}x{}, 类型: {}",
               A.rows(), A.cols(), get_backend_name());

    try {
        const int n = A.rows();

        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();

        mumps_ctx_ = std::make_unique<MumpsContext>();

        int sym = get_mumps_sym();
        FEEM_DEBUG("MUMPSBackend::set_matrix(实数) - 调用 factorize_from_csr, sym={}", sym);

        int info = mumps_ctx_->factorize_from_csr(A, sym);

        if (info != 0) {
            FEEM_ERROR("MUMPSBackend::set_matrix(实数) - MUMPS分解失败 (info={}, sym={})", info, sym);
            mumps_ctx_.reset();
            eigen_matrix_ = Eigen::SparseMatrix<double>();
            matrix_set_ = false;
            return;
        }

        matrix_set_ = true;

        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_INFO("MUMPSBackend::set_matrix(实数) 完成, 维度: {}x{}, 类型: {}, 耗时: {:.3f} ms",
                  n, n, get_backend_name(), duration_ms);

    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("MUMPSBackend::set_matrix(实数) - 内存分配失败: {}", e.what());
        mumps_ctx_.reset();
        eigen_matrix_ = Eigen::SparseMatrix<double>();
        matrix_set_ = false;

    } catch (const std::exception& e) {
        FEEM_ERROR("MUMPSBackend::set_matrix(实数) - 异常: {}", e.what());
        mumps_ctx_.reset();
        eigen_matrix_ = Eigen::SparseMatrix<double>();
        matrix_set_ = false;
    }
}

Eigen::VectorXd MUMPSBackend::solve_real(const Eigen::VectorXd& b) {
    if (!matrix_set_ || !mumps_ctx_ || !mumps_ctx_->is_factored()) {
        throw std::runtime_error("MUMPSBackend: 矩阵未设置或分解失败，请先调用 set_matrix()");
    }
    if (b.size() != eigen_matrix_.rows()) {
        throw std::invalid_argument("MUMPSBackend: 右端项向量维度不匹配");
    }

    FEEM_DEBUG("MUMPSBackend::solve_real - 开始MUMPS求解");

    SolverResult result = mumps_ctx_->solve(b);

    if (result.status != SolverStatus::SUCCESS) {
        throw std::runtime_error("MUMPSBackend: 求解失败 - " + result.error_msg);
    }

    return result.x;
}

// ============================================================================
// 复数矩阵设置与分解
// ============================================================================

void MUMPSBackend::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("MUMPSBackend::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}, 类型: {}",
               A.rows(), A.cols(), get_backend_name());

    try {
        const int n = A.rows();

        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();

        mumps_ctx_complex_ = std::make_unique<ComplexMumpsContext>();

        int sym = get_mumps_sym();
        FEEM_DEBUG("MUMPSBackend::set_matrix(复数) - 调用 factorize_from_csr, sym={}", sym);

        int info = mumps_ctx_complex_->factorize_from_csr(A, sym);

        if (info != 0) {
            FEEM_ERROR("MUMPSBackend::set_matrix(复数) - MUMPS复数分解失败 (info={}, sym={})", info, sym);
            mumps_ctx_complex_.reset();
            eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
            matrix_set_complex_ = false;
            return;
        }

        matrix_set_complex_ = true;

        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_INFO("MUMPSBackend::set_matrix(复数) 完成, 维度: {}x{}, 类型: {}, 耗时: {:.3f} ms",
                  n, n, get_backend_name(), duration_ms);

    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("MUMPSBackend::set_matrix(复数) - 内存分配失败: {}", e.what());
        mumps_ctx_complex_.reset();
        eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
        matrix_set_complex_ = false;

    } catch (const std::exception& e) {
        FEEM_ERROR("MUMPSBackend::set_matrix(复数) - 异常: {}", e.what());
        mumps_ctx_complex_.reset();
        eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
        matrix_set_complex_ = false;
    }
}

/**
 * @brief 使用缓存的 MUMPS 复数分解结果求解复数线性系统
 * @param b 复数右端项向量
 * @return Eigen::VectorXcd 复数解向量 x
 *
 * @details 委托给 ComplexMumpsContext::solve() 执行复数 MUMPS 求解，
 *          并将返回的交错实数向量转换为 VectorXcd。
 *
 * @note ComplexMumpsContext::solve() 返回的 raw_result.x 为实数向量（长度 2n），按 [Re,Im,...] 排列
 * @note 使用 extract_complex 辅助函数进行格式转换
 *
 * @note 原代码位置：solve_mumps_complex_helper() (direct_solvers_mumps_complex.cpp:121)
 */
Eigen::VectorXcd MUMPSBackend::solve_complex(const Eigen::VectorXcd& b) {
    if (!matrix_set_complex_ || !mumps_ctx_complex_ || !mumps_ctx_complex_->is_factored()) {
        throw std::runtime_error("MUMPSBackend: 复数矩阵未设置或分解失败，请先调用 set_matrix(复数)");
    }
    if (b.size() != eigen_matrix_complex_.rows()) {
        throw std::invalid_argument("MUMPSBackend: 复数右端项向量维度不匹配");
    }

    FEEM_DEBUG("MUMPSBackend::solve_complex - 开始MUMPS复数求解");

    SolverResult raw_result = mumps_ctx_complex_->solve(b);

    if (raw_result.status != SolverStatus::SUCCESS) {
        throw std::runtime_error("MUMPSBackend: 复数求解失败 - " + raw_result.error_msg);
    }

    // ComplexMumpsContext::solve() 返回的 raw_result.x 为实数向量（长度 2n），按 [Re,Im,...] 排列
    // 使用 extract_complex 辅助函数转换为 VectorXcd
    return extract_complex(raw_result);
}

// ============================================================================
// 资源清理
// ============================================================================

/**
 * @brief 清理所有 MUMPS 资源并重置状态
 * @details 释放实数和复数 MUMPS 上下文，重置所有内部标志
 */
void MUMPSBackend::clear() {
    FEEM_DEBUG("MUMPSBackend::clear - 释放所有资源");

    if (mumps_ctx_) {
        mumps_ctx_->reset();
        mumps_ctx_.reset();
    }

    if (mumps_ctx_complex_) {
        mumps_ctx_complex_->reset();
        mumps_ctx_complex_.reset();
    }

    eigen_matrix_ = Eigen::SparseMatrix<double>();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_ = false;
    matrix_set_complex_ = false;
}
#endif

} // namespace numeric
