/**
 * @file unified_direct_solver.cpp
 * @brief 统一直接求解器调度类完整实现（模板方法模式）
 * @details 实现 UnifiedDirectSolver 类的全部方法，包含所有直接求解器的公共逻辑。
 *
 * @par 实现要点：
 * - **输入校验**：统一的矩阵/向量合法性检查
 * - **后端委托**：通过 SolverBackend 接口调用具体后端
 * - **性能监控**：高精度耗时统计、残差计算
 * - **错误处理**：统一异常捕获和状态码转换
 * - **日志输出**：详细的调试信息和性能指标
 *
 * @par 代码重复消除统计：
 * **重构前**：3个求解器类 × 4个公共方法 = 12份重复代码
 * **重构后**：1个类 × 4个方法 = 4份代码（减少 67%！）
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版)
 */

#include "unified_direct_solver.h"
#include "logger_factory.hpp"
#include <stdexcept>
#include <string>

namespace numeric {

// ============================================================================
// 构造函数与属性查询
// ============================================================================

UnifiedDirectSolver::UnifiedDirectSolver(std::unique_ptr<SolverBackend> backend)
    : backend_(std::move(backend)) {

    if (!backend_) {
        throw std::invalid_argument("UnifiedDirectSolver: 后端对象不能为 nullptr");
    }

    FEEM_INFO("UnifiedDirectSolver 构造完成, 后端类型: {}", backend_->get_backend_name());
}

std::string UnifiedDirectSolver::get_solver_name() const {
    return std::string("UnifiedDirect_") + backend_->get_backend_name();
}

// ============================================================================
// 实数矩阵设置与求解（统一实现）
// ============================================================================

void UnifiedDirectSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("UnifiedDirectSolver::set_matrix(实数) - 开始设置矩阵, 维度: {}x{}, 后端: {}",
               A.rows(), A.cols(), backend_->get_backend_name());

    // 步骤1：输入合法性校验（统一实现）
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(实数) - 矩阵维度非法: rows={}, cols={}",
                   A.rows(), A.cols());
        matrix_set_ = false;
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(实数) - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        matrix_set_ = false;
        return;
    }

    // 步骤2：委托给后端执行实际的分解操作
    try {
        backend_->set_matrix(A);
        matrix_set_ = true;

        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_INFO("UnifiedDirectSolver::set_matrix(实数) 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
                  A.rows(), A.cols(), backend_->get_backend_name(), duration_ms);

    } catch (const std::exception& e) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(实数) - 后端异常: {}", e.what());
        matrix_set_ = false;
    }
}

SolverResult UnifiedDirectSolver::solve(const Eigen::VectorXd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("UnifiedDirectSolver::solve(实数) - 开始求解, 后端: {}", backend_->get_backend_name());

    SolverResult result;

    // 步骤1：输入校验（统一实现）
    if (!matrix_set_) {
        FEEM_ERROR("UnifiedDirectSolver::solve(实数) - 矩阵未设置，请先调用 set_matrix()");
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "矩阵未设置，请先调用 set_matrix()";
        return result;
    }

    try {
        // 步骤2：维度匹配检查（需要从后端获取矩阵信息或依赖异常）
        // 注意：这里不直接访问后端的内部状态，依赖后端的异常机制

        // 步骤3：委托给后端执行实际求解
        Eigen::VectorXd x = backend_->solve_real(b);

        // 步骤4：计算残差范数用于质量评估
        result.x = std::move(x);
        result.residual_norm = compute_residual_norm_real(result.x, b);
        result.status = SolverStatus::SUCCESS;

        auto end_time = std::chrono::high_resolution_clock::now();
        result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_DEBUG("UnifiedDirectSolver::solve(实数) 成功, 残差: {:.6e}, 耗时: {:.3f} ms",
                   result.residual_norm, result.solve_time_ms);

    } catch (const std::invalid_argument& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(实数) - 参数错误: {}", e.what());
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = e.what();

    } catch (const std::runtime_error& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(实数) - 求解失败: {}", e.what());
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = e.what();

    } catch (const std::exception& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(实数) - 异常: {}", e.what());
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = std::string("未知异常: ") + e.what();
    }

    return result;
}

// ============================================================================
// 复数矩阵设置与求解（统一实现）
// ============================================================================

void UnifiedDirectSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("UnifiedDirectSolver::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}, 后端: {}",
               A.rows(), A.cols(), backend_->get_backend_name());

    // 步骤1：检查后端是否支持复数运算
    if (!backend_->supports_complex()) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(复数) - 后端 {} 不支持复数运算",
                   backend_->get_backend_name());
        matrix_set_complex_ = false;
        return;
    }

    // 步骤2：输入合法性校验
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(复数) - 矩阵维度非法: rows={}, cols={}",
                   A.rows(), A.cols());
        matrix_set_complex_ = false;
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(复数) - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        matrix_set_complex_ = false;
        return;
    }

    // 步骤3：委托给后端执行实际的分解操作
    try {
        backend_->set_matrix(A);
        matrix_set_complex_ = true;

        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_INFO("UnifiedDirectSolver::set_matrix(复数) 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
                  A.rows(), A.cols(), backend_->get_backend_name(), duration_ms);

    } catch (const std::exception& e) {
        FEEM_ERROR("UnifiedDirectSolver::set_matrix(复数) - 后端异常: {}", e.what());
        matrix_set_complex_ = false;
    }
}

SolverResult UnifiedDirectSolver::solve(const Eigen::VectorXcd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    FEEM_DEBUG("UnifiedDirectSolver::solve(复数) - 开始复数求解, 后端: {}", backend_->get_backend_name());

    SolverResult result;

    // 步骤1：输入校验
    if (!matrix_set_complex_) {
        FEEM_ERROR("UnifiedDirectSolver::solve(复数) - 复数矩阵未设置，请先调用 set_matrix(复数)");
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "复数矩阵未设置，请先调用 set_matrix(复数版本)";
        return result;
    }

    if (!backend_->supports_complex()) {
        FEEM_ERROR("UnifiedDirectSolver::solve(复数) - 后端不支持复数运算");
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "当前后端不支持复数运算";
        return result;
    }

    try {
        // 步骤2：委托给后端执行实际求解
        Eigen::VectorXcd x = backend_->solve_complex(b);

        // 步骤3：计算残差范数
        result.x_complex = std::move(x);
        result.residual_norm = compute_residual_norm_complex(result.x_complex, b);
        result.status = SolverStatus::SUCCESS;

        auto end_time = std::chrono::high_resolution_clock::now();
        result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        FEEM_DEBUG("UnifiedDirectSolver::solve(复数) 成功, 残差: {:.6e}, 耗时: {:.3f} ms",
                   result.residual_norm, result.solve_time_ms);

    } catch (const std::invalid_argument& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(复数) - 参数错误: {}", e.what());
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = e.what();

    } catch (const std::runtime_error& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(复数) - 求解失败: {}", e.what());
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = e.what();

    } catch (const std::exception& e) {
        FEEM_ERROR("UnifiedDirectSolver::solve(复数) - 异常: {}", e.what());
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = std::string("未知异常: ") + e.what();
    }

    return result;
}

// ============================================================================
// 资源清理
// ============================================================================

void UnifiedDirectSolver::clear() {
    FEEM_DEBUG("UnifiedDirectSolver::clear - 释放所有资源, 后端: {}", backend_->get_backend_name());

    if (backend_) {
        backend_->clear();
    }

    matrix_set_ = false;
    matrix_set_complex_ = false;
}

// ============================================================================
// 私有辅助方法：精确残差计算
// ============================================================================

/**
 * @brief 计算实数解的精确残差范数 ||Ax - b||_2
 * @param x 解向量
 * @param b 右端项向量
 * @return double 残差的2-范数
 *
 * @details 通过后端的 get_eigen_matrix_real() 接口获取缓存的系数矩阵，
 *          计算精确的相对残差范数，用于评估求解质量。
 *
 * @par 计算公式：
 * residual = ||Ax - b||_2 / ||b||_2  （相对残差）
 */
double UnifiedDirectSolver::compute_residual_norm_real(
    const Eigen::VectorXd& x, const Eigen::VectorXd& b) {

    // 从后端获取缓存的系数矩阵
    const Eigen::SparseMatrix<double>* A = backend_->get_eigen_matrix_real();

    if (!A || A->rows() == 0) {
        FEEM_WARN("UnifiedDirectSolver::compute_residual_norm_real - 矩阵未缓存，无法计算残差");
        return 0.0;
    }

    // 计算残差 r = Ax - b
    Eigen::VectorXd residual = (*A) * x - b;

    // 计算2-范数
    double norm_residual = residual.norm();
    double norm_b = b.norm();

    // 返回相对残差（避免除以零）
    if (norm_b > 1e-15) {
        return norm_residual / norm_b;
    } else {
        return norm_residual;
    }
}

/**
 * @brief 计算复数解的精确残差范数 ||Ax - b||_2
 * @param x 复数解向量
 * @param b 复数右端项向量
 * @return double 残差的2-范数
 *
 * @details 复数版本的精确残差计算。
 */
double UnifiedDirectSolver::compute_residual_norm_complex(
    const Eigen::VectorXcd& x, const Eigen::VectorXcd& b) {

    // 从后端获取缓存的复数系数矩阵
    const Eigen::SparseMatrix<std::complex<double>>* A = backend_->get_eigen_matrix_complex();

    if (!A || A->rows() == 0) {
        FEEM_WARN("UnifiedDirectSolver::compute_residual_norm_complex - 矩阵未缓存，无法计算残差");
        return 0.0;
    }

    // 计算复数残差 r = Ax - b
    Eigen::VectorXcd residual = (*A) * x - b;

    // 计算2-范数
    double norm_residual = residual.norm();
    double norm_b = b.norm();

    // 返回相对残差（避免除以零）
    if (norm_b > 1e-15) {
        return norm_residual / norm_b;
    } else {
        return norm_residual;
    }
}

} // namespace numeric
