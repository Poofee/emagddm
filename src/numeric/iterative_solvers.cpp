/**
 * @file iterative_solvers.cpp
 * @brief 核心数值层 - 高性能迭代求解器与AMG预条件子框架实现
 * @details 实现CGSolver、BiCGSTABSolver、ScalarAMG等组件的完整功能，
 *          包括预条件子设置、迭代求解核心逻辑、AMG多层次构建和V/W循环。
 *
 * @par 实现策略：
 * - CG/BiCGSTAB: 基于Eigen迭代求解器框架，手动实现预条件子集成以获得更细粒度的控制
 * - ScalarAMG: 简化但功能完整的AMG实现，采用交替粗化和直接插值
 * - 所有求解器支持收敛历史记录、发散检测、性能统计
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "iterative_solvers.h"
#include "logger_factory.hpp"

namespace numeric {

// ==================== ConvergenceHistory 实现 ====================

void ConvergenceHistory::record_iteration(int iter, double residual_norm) {
    iterations_.push_back(iter);
    residuals_.push_back(residual_norm);
}

std::vector<int> ConvergenceHistory::get_iterations() const {
    return iterations_;
}

std::vector<double> ConvergenceHistory::get_residuals() const {
    return residuals_;
}

bool ConvergenceHistory::has_converged(double tolerance) const {
    if (residuals_.empty()) {
        return false;
    }
    return residuals_.back() < tolerance;
}

int ConvergenceHistory::total_iterations() const {
    return static_cast<int>(iterations_.size());
}

double ConvergenceHistory::final_residual() const {
    if (residuals_.empty()) {
        return -1.0;
    }
    return residuals_.back();
}

void ConvergenceHistory::clear() {
    iterations_.clear();
    residuals_.clear();
}

// ==================== CGSolver 实现 ====================

CGSolver::CGSolver(const CGConfig& config)
    : config_(config) {
    FEEM_DEBUG("CGSolver 构造完成, tolerance={}, max_iter={}, precond={}",
               config_.tolerance, config_.max_iterations,
               static_cast<int>(config_.preconditioner));
}

void CGSolver::set_config(const CGConfig& config) {
    config_ = config;
    matrix_set_ = false;  // 配置变更后需重新设置矩阵
    FEEM_DEBUG("CGSolver 配置已更新");
}

CGConfig CGSolver::get_config() const {
    return config_;
}

ConvergenceHistory CGSolver::get_convergence_history() const {
    return history_;
}

void CGSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    eigen_matrix_ = SparseConverter::to_eigen(A);
    matrix_set_ = true;

    setup_preconditioner();
    history_.clear();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);

    FEEM_INFO("CGSolver::set_matrix 完成, 矩阵尺寸: {}x{}, 非零元数: {}, 耗时: {:.3f}ms",
              A.rows(), A.cols(), A.nnz(), duration.count());
}

std::string CGSolver::get_solver_name() const {
    std::string name = "CGSolver";
    switch (config_.preconditioner) {
        case CGConfig::PreconditionerType::NONE:
            name += "_None";
            break;
        case CGConfig::PreconditionerType::JACOBI:
            name += "_Jacobi";
            break;
        case CGConfig::PreconditionerType::ILU0:
            name += "_ILU0";
            break;
    }
    return name;
}

void CGSolver::clear() {
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    jacobi_diag_ = Eigen::VectorXd();
    // 注意: IncompleteLUT不支持拷贝/移动赋值，通过matrix_set_=false确保下次set_matrix时重新初始化
    matrix_set_ = false;
    history_.clear();

    // 清理复数版本状态
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_complex_ = false;
    eigen_solver_complex_.reset();

    FEEM_DEBUG("CGSolver 资源已清理（含复数版本）");
}

/**
 * @brief 设置复数系数矩阵并配置预条件子（CGSolver复数版本）
 * @param A CSR格式复数稀疏系数矩阵
 *
 * @details 实现步骤：
 * 1. 输入校验（矩阵非空、方阵）
 * 2. 使用SparseConverter将CsrMatrix转换为Eigen格式
 * 3. 创建Eigen ConjugateGradient求解器实例（带IncompleteCholesky预条件子）
 * 4. 配置收敛容差和最大迭代次数
 * 5. 执行compute()预计算预条件子
 * 6. 设置matrix_set_complex_标志
 */
void CGSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (A.rows() == 0 || A.cols() == 0) {
        FEEM_ERROR("CGSolver::set_matrix(复数) 错误: 矩阵为空");
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("CGSolver::set_matrix(复数) 错误: 矩阵不是方阵, {}x{}", A.rows(), A.cols());
        return;
    }

    // 转换为Eigen格式
    eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);

    // 创建或重置复数Eigen求解器实例
    eigen_solver_complex_ = std::make_unique<
        Eigen::ConjugateGradient<
            Eigen::SparseMatrix<std::complex<double>>,
            Eigen::Lower|Eigen::Upper,
            Eigen::IncompleteCholesky<std::complex<double>>>>();

    // 配置求解器参数
    eigen_solver_complex_->setTolerance(config_.tolerance);
    eigen_solver_complex_->setMaxIterations(config_.max_iterations);

    // 预计算：分析矩阵并构建预条件子
    eigen_solver_complex_->compute(eigen_matrix_complex_);

    // 检查预计算是否成功
    if (eigen_solver_complex_->info() != Eigen::Success) {
        FEEM_ERROR("CGSolver::set_matrix(复数) 错误: Eigen ConjugateGradient compute()失败");
        matrix_set_complex_ = false;
        return;
    }

    matrix_set_complex_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);

    FEEM_INFO("CGSolver::set_matrix(复数) 完成, 矩阵尺寸: {}x{}, 非零元数: {}, 耗时: {:.3f}ms",
              A.rows(), A.cols(), A.nnz(), duration.count());
}

/**
 * @brief 执行CG迭代求解复数线性系统（CGSolver复数版本）
 * @param b 复数右端项向量
 * @return SolverResult 求解结果，复数解存储在x_complex字段
 *
 * @details 求解流程：
 * 1. 检查复数矩阵是否已设置（matrix_set_complex_）
 * 2. 验证向量维度与矩阵匹配
 * 3. 调用Eigen ConjugateGradient求解器的solve()方法
 * 4. 提取解向量、迭代次数、误差估计值
 * 5. 判断收敛状态并组装返回结果
 */
SolverResult CGSolver::solve(const Eigen::VectorXcd& b) {
    SolverResult result;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验：检查复数矩阵是否已设置
    if (!matrix_set_complex_) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "复数矩阵未设置，请先调用set_matrix(CsrMatrix<complex>)";
        FEEM_ERROR("CGSolver::solve(复数) 错误: {}", result.error_msg);
        return result;
    }

    // 输入校验：维度匹配检查
    if (b.size() != eigen_matrix_complex_.rows()) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "复数右端项向量维度与矩阵不匹配";
        FEEM_ERROR("CGSolver::solve(复数) 错误: {}, 矩阵维度{}, 向量维度{}",
                   result.error_msg, eigen_matrix_complex_.rows(), b.size());
        return result;
    }

    // 执行Eigen ConjugateGradient求解
    Eigen::VectorXcd x = eigen_solver_complex_->solve(b);

    // 提取求解结果信息
    result.x_complex = x;
    result.iterations = eigen_solver_complex_->iterations();
    result.residual_norm = eigen_solver_complex_->error();

    // 判断收敛状态
    if (eigen_solver_complex_->info() == Eigen::Success && result.residual_norm < config_.tolerance) {
        result.status = SolverStatus::SUCCESS;
    } else if (result.iterations >= config_.max_iterations) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        result.error_msg = "达到最大迭代次数但未收敛";
        FEEM_WARN("CGSolver::solve(复数): 达到最大迭代次数{}, 最终误差={:.6e}",
                  config_.max_iterations, result.residual_norm);
    } else {
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = "Eigen求解器报告数值错误";
        FEEM_ERROR("CGSolver::solve(复数): 数值错误, error={:.6e}", result.residual_norm);
    }

    // 统计耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    // 记录收敛历史（复数版本使用独立的history或可扩展）
    if (result.iterations > 0) {
        history_.record_iteration(result.iterations, result.residual_norm);
    }

    FEEM_INFO("CGSolver::solve(复数) 完成, 状态: {}, 迭代次数: {}, 误差: {:.6e}, 耗时: {:.3f}ms",
              static_cast<int>(result.status), result.iterations,
              result.residual_norm, result.solve_time_ms);

    return result;
}


void CGSolver::setup_preconditioner() {
    switch (config_.preconditioner) {
        case CGConfig::PreconditionerType::JACOBI: {
            int n = eigen_matrix_.rows();
            jacobi_diag_ = Eigen::VectorXd(n);
            for (int i = 0; i < n; ++i) {
                double diag = eigen_matrix_.coeff(i, i);
                if (std::abs(diag) > 1e-15) {
                    jacobi_diag_(i) = 1.0 / diag;
                } else {
                    jacobi_diag_(i) = 1.0;
                    FEEM_WARN("CGSolver Jacobi预条件子: 对角线元素{}接近零, 使用默认值", i);
                }
            }
            FEEM_DEBUG("CGSolver Jacobi预条件子设置完成");
            break;
        }
        case CGConfig::PreconditionerType::ILU0: {
            ilu0_precond_.compute(eigen_matrix_);
            if (ilu0_precond_.info() != Eigen::Success) {
                FEEM_ERROR("CGSolver ILU(0)分解失败, 回退到无预条件子模式");
                config_.preconditioner = CGConfig::PreconditionerType::NONE;
            } else {
                FEEM_DEBUG("CGSolver ILU(0)预条件子设置完成");
            }
            break;
        }
        case CGConfig::PreconditionerType::NONE:
        default:
            break;
    }
}

SolverResult CGSolver::solve(const Eigen::VectorXd& b) {
    SolverResult result;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (!matrix_set_) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "矩阵未设置，请先调用set_matrix()";
        FEEM_ERROR("CGSolver::solve 错误: {}", result.error_msg);
        return result;
    }

    if (b.size() != eigen_matrix_.rows()) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "右端项向量维度与矩阵不匹配";
        FEEM_ERROR("CGSolver::solve 错误: {}, 矩阵维度{}, 向量维度{}",
                   result.error_msg, eigen_matrix_.rows(), b.size());
        return result;
    }

    // 选择预条件子分支
    switch (config_.preconditioner) {
        case CGConfig::PreconditionerType::NONE:
            result = solve_no_preconditioner(b);
            break;
        case CGConfig::PreconditionerType::JACOBI:
            result = solve_with_jacobi(b);
            break;
        case CGConfig::PreconditionerType::ILU0:
            result = solve_with_ilu0(b);
            break;
    }

    // 记录收敛历史
    if (result.iterations > 0) {
        history_.record_iteration(result.iterations, result.residual_norm);
    }

    // 统计耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("CGSolver::solve 完成, 状态: {}, 迭代次数: {}, 残差: {:.6e}, 耗时: {:.3f}ms",
              static_cast<int>(result.status), result.iterations,
              result.residual_norm, result.solve_time_ms);

    return result;
}

SolverResult CGSolver::solve_no_preconditioner(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd p = r;
    double rs_old = r.squaredNorm();
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = std::sqrt(rs_old) / b_norm;
    int iter;
    double prev_residual = residual_norm;

    for (iter = 1; iter <= max_iter; ++iter) {
        Eigen::VectorXd Ap = eigen_matrix_ * p;
        double pAp = p.dot(Ap);

        if (std::abs(pAp) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "矩阵可能奇异或数值不稳定（p^T Ap ≈ 0）";
            FEEM_WARN("CGSolver: 第{}次迭代检测到pAp≈0, 可能矩阵奇异", iter);
            break;
        }

        double alpha = rs_old / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        double rs_new = r.squaredNorm();
        residual_norm = std::sqrt(rs_new) / b_norm;

        // 记录残差（按配置间隔）
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("CGSolver 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测：残差异常增长
        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长，可能发散";
            FEEM_WARN("CGSolver: 可能发散于第{}次迭代, 残差从{:.2e}增至{:.2e}",
                       iter, prev_residual, residual_norm);
            break;
        }

        double beta = rs_new / rs_old;
        p = r + beta * p;
        rs_old = rs_new;
        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("CGSolver: 达到最大迭代次数{}, 最终残差={:.6e}", max_iter, residual_norm);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

SolverResult CGSolver::solve_with_jacobi(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd z = jacobi_diag_.cwiseProduct(r);  // Jacobi预条件: z = D^{-1}*r
    Eigen::VectorXd p = z;
    double rz_old = r.dot(z);
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = r.norm() / b_norm;
    int iter;
    double prev_residual = residual_norm;

    for (iter = 1; iter <= max_iter; ++iter) {
        Eigen::VectorXd Ap = eigen_matrix_ * p;
        double pAp = p.dot(Ap);

        if (std::abs(pAp) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "矩阵可能奇异或数值不稳定（p^T Ap ≈ 0）";
            break;
        }

        double alpha = rz_old / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        z = jacobi_diag_.cwiseProduct(r);  // 应用Jacobi预条件子
        double rz_new = r.dot(z);
        residual_norm = r.norm() / b_norm;

        // 残差记录
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("CGSolver+Jacobi 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测
        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长，可能发散";
            FEEM_WARN("CGSolver+Jacobi: 可能发散于第{}次迭代", iter);
            break;
        }

        double beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("CGSolver+Jacobi: 达到最大迭代次数{}", max_iter);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

SolverResult CGSolver::solve_with_ilu0(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd z = ilu0_precond_.solve(r);  // ILU0预条件
    Eigen::VectorXd p = z;
    double rz_old = r.dot(z);
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = r.norm() / b_norm;
    int iter;
    double prev_residual = residual_norm;

    for (iter = 1; iter <= max_iter; ++iter) {
        Eigen::VectorXd Ap = eigen_matrix_ * p;
        double pAp = p.dot(Ap);

        if (std::abs(pAp) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "矩阵可能奇异或数值不稳定（p^T Ap ≈ 0）";
            break;
        }

        double alpha = rz_old / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        z = ilu0_precond_.solve(r);  // 应用ILU0预条件子
        double rz_new = r.dot(z);
        residual_norm = r.norm() / b_norm;

        // 残差记录
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("CGSolver+ILU0 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测
        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长，可能发散";
            FEEM_WARN("CGSolver+ILU0: 可能发散于第{}次迭代", iter);
            break;
        }

        double beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("CGSolver+ILU0: 达到最大迭代次数{}", max_iter);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

// ==================== BiCGSTABSolver 实现 ====================

BiCGSTABSolver::BiCGSTABSolver(const BiCGSTABConfig& config)
    : config_(config) {
    FEEM_DEBUG("BiCGSTABSolver 构造完成, tolerance={}, max_iter={}",
               config_.tolerance, config_.max_iterations);
}

void BiCGSTABSolver::set_config(const BiCGSTABConfig& config) {
    config_ = config;
    matrix_set_ = false;
    FEEM_DEBUG("BiCGSTABSolver 配置已更新");
}

BiCGSTABConfig BiCGSTABSolver::get_config() const {
    return config_;
}

ConvergenceHistory BiCGSTABSolver::get_convergence_history() const {
    return history_;
}

void BiCGSTABSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    eigen_matrix_ = SparseConverter::to_eigen(A);
    matrix_set_ = true;

    setup_preconditioner();
    history_.clear();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);

    FEEM_INFO("BiCGSTABSolver::set_matrix 完成, 矩阵尺寸: {}x{}, 非零元数: {}, 耗时: {:.3f}ms",
              A.rows(), A.cols(), A.nnz(), duration.count());
}

std::string BiCGSTABSolver::get_solver_name() const {
    std::string name = "BiCGSTABSolver";
    switch (config_.preconditioner) {
        case BiCGSTABConfig::PreconditionerType::NONE:
            name += "_None";
            break;
        case BiCGSTABConfig::PreconditionerType::JACOBI:
            name += "_Jacobi";
            break;
        case BiCGSTABConfig::PreconditionerType::ILU0:
            name += "_ILU0";
            break;
    }
    return name;
}

void BiCGSTABSolver::clear() {
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    jacobi_diag_ = Eigen::VectorXd();
    // 注意: IncompleteLUT不支持拷贝/移动赋值，通过matrix_set_=false确保下次set_matrix时重新初始化
    matrix_set_ = false;
    history_.clear();

    // 清理复数版本状态
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_complex_ = false;
    eigen_solver_complex_.reset();

    FEEM_DEBUG("BiCGSTABSolver 资源已清理（含复数版本）");
}

/**
 * @brief 设置复数系数矩阵并配置预条件子（BiCGSTABSolver复数版本）
 * @param A CSR格式复数稀疏系数矩阵
 *
 * @details 实现步骤：
 * 1. 输入校验（矩阵非空、方阵）
 * 2. 使用SparseConverter将CsrMatrix转换为Eigen格式
 * 3. 创建Eigen BiCGSTAB求解器实例（带IncompleteLUT预条件子）
 * 4. 配置收敛容差和最大迭代次数
 * 5. 执行compute()预计算预条件子
 * 6. 设置matrix_set_complex_标志
 */
void BiCGSTABSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (A.rows() == 0 || A.cols() == 0) {
        FEEM_ERROR("BiCGSTABSolver::set_matrix(复数) 错误: 矩阵为空");
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("BiCGSTABSolver::set_matrix(复数) 错误: 矩阵不是方阵, {}x{}", A.rows(), A.cols());
        return;
    }

    // 转换为Eigen格式
    eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);

    // 创建或重置复数Eigen求解器实例
    eigen_solver_complex_ = std::make_unique<
        Eigen::BiCGSTAB<
            Eigen::SparseMatrix<std::complex<double>>,
            Eigen::IncompleteLUT<std::complex<double>>>>();

    // 配置求解器参数
    eigen_solver_complex_->setTolerance(config_.tolerance);
    eigen_solver_complex_->setMaxIterations(config_.max_iterations);

    // 预计算：分析矩阵并构建预条件子
    eigen_solver_complex_->compute(eigen_matrix_complex_);

    // 检查预计算是否成功
    if (eigen_solver_complex_->info() != Eigen::Success) {
        FEEM_ERROR("BiCGSTABSolver::set_matrix(复数) 错误: Eigen BiCGSTAB compute()失败");
        matrix_set_complex_ = false;
        return;
    }

    matrix_set_complex_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);

    FEEM_INFO("BiCGSTABSolver::set_matrix(复数) 完成, 矩阵尺寸: {}x{}, 非零元数: {}, 耗时: {:.3f}ms",
              A.rows(), A.cols(), A.nnz(), duration.count());
}

/**
 * @brief 执行BiCGSTAB迭代求解复数线性系统（BiCGSTABSolver复数版本）
 * @param b 复数右端项向量
 * @return SolverResult 求解结果，复数解存储在x_complex字段
 *
 * @details 求解流程：
 * 1. 检查复数矩阵是否已设置（matrix_set_complex_）
 * 2. 验证向量维度与矩阵匹配
 * 3. 调用Eigen BiCGSTAB求解器的solve()方法
 * 4. 提取解向量、迭代次数、误差估计值
 * 5. 判断收敛状态并组装返回结果
 *
 * @note BiCGSTAB适用于非对称/非正定复数矩阵（如涡流场、谐波电磁场问题）
 */
SolverResult BiCGSTABSolver::solve(const Eigen::VectorXcd& b) {
    SolverResult result;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验：检查复数矩阵是否已设置
    if (!matrix_set_complex_) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "复数矩阵未设置，请先调用set_matrix(CsrMatrix<complex>)";
        FEEM_ERROR("BiCGSTABSolver::solve(复数) 错误: {}", result.error_msg);
        return result;
    }

    // 输入校验：维度匹配检查
    if (b.size() != eigen_matrix_complex_.rows()) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "复数右端项向量维度与矩阵不匹配";
        FEEM_ERROR("BiCGSTABSolver::solve(复数) 错误: {}, 矩阵维度{}, 向量维度{}",
                   result.error_msg, eigen_matrix_complex_.rows(), b.size());
        return result;
    }

    // 执行Eigen BiCGSTAB求解
    Eigen::VectorXcd x = eigen_solver_complex_->solve(b);

    // 提取求解结果信息
    result.x_complex = x;
    result.iterations = eigen_solver_complex_->iterations();
    result.residual_norm = eigen_solver_complex_->error();

    // 判断收敛状态
    if (eigen_solver_complex_->info() == Eigen::Success && result.residual_norm < config_.tolerance) {
        result.status = SolverStatus::SUCCESS;
    } else if (result.iterations >= config_.max_iterations) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        result.error_msg = "达到最大迭代次数但未收敛";
        FEEM_WARN("BiCGSTABSolver::solve(复数): 达到最大迭代次数{}, 最终误差={:.6e}",
                  config_.max_iterations, result.residual_norm);
    } else {
        result.status = SolverStatus::NUMERICAL_ERROR;
        result.error_msg = "Eigen求解器报告数值错误";
        FEEM_ERROR("BiCGSTABSolver::solve(复数): 数值错误, error={:.6e}", result.residual_norm);
    }

    // 统计耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    // 记录收敛历史（复数版本使用独立的history或可扩展）
    if (result.iterations > 0) {
        history_.record_iteration(result.iterations, result.residual_norm);
    }

    FEEM_INFO("BiCGSTABSolver::solve(复数) 完成, 状态: {}, 迭代次数: {}, 误差: {:.6e}, 耗时: {:.3f}ms",
              static_cast<int>(result.status), result.iterations,
              result.residual_norm, result.solve_time_ms);

    return result;
}


void BiCGSTABSolver::setup_preconditioner() {
    switch (config_.preconditioner) {
        case BiCGSTABConfig::PreconditionerType::JACOBI: {
            int n = eigen_matrix_.rows();
            jacobi_diag_ = Eigen::VectorXd(n);
            for (int i = 0; i < n; ++i) {
                double diag = eigen_matrix_.coeff(i, i);
                if (std::abs(diag) > 1e-15) {
                    jacobi_diag_(i) = 1.0 / diag;
                } else {
                    jacobi_diag_(i) = 1.0;
                    FEEM_WARN("BiCGSTABSolver Jacobi预条件子: 对角线元素{}接近零", i);
                }
            }
            FEEM_DEBUG("BiCGSTABSolver Jacobi预条件子设置完成");
            break;
        }
        case BiCGSTABConfig::PreconditionerType::ILU0: {
            ilu0_precond_.compute(eigen_matrix_);
            if (ilu0_precond_.info() != Eigen::Success) {
                FEEM_ERROR("BiCGSTABSolver ILU(0)分解失败, 回退到无预条件子模式");
                config_.preconditioner = BiCGSTABConfig::PreconditionerType::NONE;
            } else {
                FEEM_DEBUG("BiCGSTABSolver ILU(0)预条件子设置完成");
            }
            break;
        }
        case BiCGSTABConfig::PreconditionerType::NONE:
        default:
            break;
    }
}

SolverResult BiCGSTABSolver::solve(const Eigen::VectorXd& b) {
    SolverResult result;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (!matrix_set_) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "矩阵未设置，请先调用set_matrix()";
        FEEM_ERROR("BiCGSTABSolver::solve 错误: {}", result.error_msg);
        return result;
    }

    if (b.size() != eigen_matrix_.rows()) {
        result.status = SolverStatus::INVALID_INPUT;
        result.error_msg = "右端项向量维度与矩阵不匹配";
        FEEM_ERROR("BiCGSTABSolver::solve 错误: {}", result.error_msg);
        return result;
    }

    // 选择预条件子分支
    switch (config_.preconditioner) {
        case BiCGSTABConfig::PreconditionerType::NONE:
            result = solve_no_preconditioner(b);
            break;
        case BiCGSTABConfig::PreconditionerType::JACOBI:
            result = solve_with_jacobi(b);
            break;
        case BiCGSTABConfig::PreconditionerType::ILU0:
            result = solve_with_ilu0(b);
            break;
    }

    // 记录收敛历史
    if (result.iterations > 0) {
        history_.record_iteration(result.iterations, result.residual_norm);
    }

    // 统计耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("BiCGSTABSolver::solve 完成, 状态: {}, 迭代次数: {}, 残差: {:.6e}, 耗时: {:.3f}ms",
              static_cast<int>(result.status), result.iterations,
              result.residual_norm, result.solve_time_ms);

    return result;
}

SolverResult BiCGSTABSolver::solve_no_preconditioner(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd r_hat = r;  // 任意选择，通常 r_hat = r
    Eigen::VectorXd p = r;
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    double rho = r.dot(r_hat);
    double alpha = 1.0;
    double omega = 1.0;
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = r.norm() / b_norm;
    int iter;
    double prev_residual = residual_norm;
    int divergence_count = 0;

    for (iter = 1; iter <= max_iter; ++iter) {
        double rho_old = rho;
        rho = r.dot(r_hat);

        if (std::abs(rho) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "BiCGSTAB算法中断（rho ≈ 0）";
            FEEM_WARN("BiCGSTABSolver: 第{}次迭代rho≈0", iter);
            break;
        }

        double beta = (rho / rho_old) * (alpha / omega);
        p = r + beta * (p - omega * v);
        v = eigen_matrix_ * p;

        double rv = r_hat.dot(v);
        if (std::abs(rv) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "算法中断（r_hat^T v ≈ 0）";
            break;
        }
        alpha = rho / rv;

        Eigen::VectorXd s = r - alpha * v;
        residual_norm = s.norm() / b_norm;

        // 检查中间残差是否足够小（提前退出优化）
        if (residual_norm < tol) {
            x += alpha * p;
            result.x = x;
            result.status = SolverStatus::SUCCESS;
            result.iterations = iter;
            result.residual_norm = residual_norm;
            history_.record_iteration(iter, residual_norm);
            return result;
        }

        // 应用预条件子（此处为恒等，z = s）
        Eigen::VectorXd z = s;  // 无预条件子
        Eigen::VectorXd t = eigen_matrix_ * z;

        double t_dot_t = t.dot(t);
        if (std::abs(t_dot_t) < 1e-30) {
            omega = 0.0;
        } else {
            omega = z.dot(t) / t_dot_t;
        }

        x += alpha * p + omega * z;
        r = s - omega * t;
        residual_norm = r.norm() / b_norm;

        // 残差记录
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("BiCGSTABSolver 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测（比CG更严格，因为BiCGSTAB更容易出现振荡）
        if (residual_norm > prev_residual * 1.5 && iter > 20) {
            divergence_count++;
            if (divergence_count >= 3) {
                result.status = SolverStatus::DIVERGED;
                result.error_msg = "检测到持续残差增长，可能发散";
                FEEM_WARN("BiCGSTABSolver: 可能发散于第{}次迭代, 连续{}次残差增长",
                           iter, divergence_count);
                break;
            }
        } else {
            divergence_count = 0;
        }

        // 绝对发散检测
        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长，可能发散";
            FEEM_WARN("BiCGSTABSolver: 可能发散于第{}次迭代, 残差从{:.2e}增至{:.2e}",
                       iter, prev_residual, residual_norm);
            break;
        }

        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("BiCGSTABSolver: 达到最大迭代次数{}, 最终残差={:.6e}", max_iter, residual_norm);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

SolverResult BiCGSTABSolver::solve_with_jacobi(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd r_hat = r;
    Eigen::VectorXd p = r;
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    double rho = r.dot(r_hat);
    double alpha = 1.0;
    double omega = 1.0;
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = r.norm() / b_norm;
    int iter;
    double prev_residual = residual_norm;
    int divergence_count = 0;

    for (iter = 1; iter <= max_iter; ++iter) {
        double rho_old = rho;
        rho = r.dot(r_hat);

        if (std::abs(rho) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "BiCGSTAB算法中断（rho ≈ 0）";
            break;
        }

        double beta = (rho / rho_old) * (alpha / omega);
        p = r + beta * (p - omega * v);
        v = eigen_matrix_ * p;

        double rv = r_hat.dot(v);
        if (std::abs(rv) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "算法中断（r_hat^T v ≈ 0）";
            break;
        }
        alpha = rho / rv;

        Eigen::VectorXd s = r - alpha * v;
        residual_norm = s.norm() / b_norm;

        if (residual_norm < tol) {
            x += alpha * p;
            result.x = x;
            result.status = SolverStatus::SUCCESS;
            result.iterations = iter;
            result.residual_norm = residual_norm;
            history_.record_iteration(iter, residual_norm);
            return result;
        }

        // 应用Jacobi预条件子
        Eigen::VectorXd z = jacobi_diag_.cwiseProduct(s);
        Eigen::VectorXd t = eigen_matrix_ * z;

        double t_dot_t = t.dot(t);
        omega = (std::abs(t_dot_t) < 1e-30) ? 0.0 : z.dot(t) / t_dot_t;

        x += alpha * p + omega * z;
        r = s - omega * t;
        residual_norm = r.norm() / b_norm;

        // 残差记录
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("BiCGSTABSolver+Jacobi 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测
        if (residual_norm > prev_residual * 1.5 && iter > 20) {
            divergence_count++;
            if (divergence_count >= 3) {
                result.status = SolverStatus::DIVERGED;
                result.error_msg = "检测到持续残差增长";
                FEEM_WARN("BiCGSTABSolver+Jacobi: 可能发散于第{}次迭代", iter);
                break;
            }
        } else {
            divergence_count = 0;
        }

        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长";
            FEEM_WARN("BiCGSTABSolver+Jacobi: 可能发散于第{}次迭代", iter);
            break;
        }

        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("BiCGSTABSolver+Jacobi: 达到最大迭代次数{}", max_iter);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

SolverResult BiCGSTABSolver::solve_with_ilu0(const Eigen::VectorXd& b) {
    SolverResult result;
    int n = eigen_matrix_.rows();
    int max_iter = config_.max_iterations;
    double tol = config_.tolerance;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r = b - eigen_matrix_ * x;
    Eigen::VectorXd r_hat = r;
    Eigen::VectorXd p = r;
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    double rho = r.dot(r_hat);
    double alpha = 1.0;
    double omega = 1.0;
    double b_norm = b.norm();

    if (b_norm < 1e-15) {
        result.x = x;
        result.status = SolverStatus::SUCCESS;
        result.iterations = 0;
        result.residual_norm = 0.0;
        return result;
    }

    double residual_norm = r.norm() / b_norm;
    int iter;
    double prev_residual = residual_norm;
    int divergence_count = 0;

    for (iter = 1; iter <= max_iter; ++iter) {
        double rho_old = rho;
        rho = r.dot(r_hat);

        if (std::abs(rho) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "BiCGSTAB算法中断（rho ≈ 0）";
            break;
        }

        double beta = (rho / rho_old) * (alpha / omega);
        p = r + beta * (p - omega * v);
        v = eigen_matrix_ * p;

        double rv = r_hat.dot(v);
        if (std::abs(rv) < 1e-30) {
            result.status = SolverStatus::NUMERICAL_ERROR;
            result.error_msg = "算法中断（r_hat^T v ≈ 0）";
            break;
        }
        alpha = rho / rv;

        Eigen::VectorXd s = r - alpha * v;
        residual_norm = s.norm() / b_norm;

        if (residual_norm < tol) {
            x += alpha * p;
            result.x = x;
            result.status = SolverStatus::SUCCESS;
            result.iterations = iter;
            result.residual_norm = residual_norm;
            history_.record_iteration(iter, residual_norm);
            return result;
        }

        // 应用ILU0预条件子
        Eigen::VectorXd z = ilu0_precond_.solve(s);
        Eigen::VectorXd t = eigen_matrix_ * z;

        double t_dot_t = t.dot(t);
        omega = (std::abs(t_dot_t) < 1e-30) ? 0.0 : z.dot(t) / t_dot_t;

        x += alpha * p + omega * z;
        r = s - omega * t;
        residual_norm = r.norm() / b_norm;

        // 残差记录
        if (iter % config_.residual_monitor_interval == 0 || iter == 1) {
            history_.record_iteration(iter, residual_norm);
            FEEM_DEBUG("BiCGSTABSolver+ILU0 迭代 {}: 相对残差={:.6e}", iter, residual_norm);
        }

        // 收敛检查
        if (residual_norm < tol) {
            result.status = SolverStatus::SUCCESS;
            break;
        }

        // 发散检测
        if (residual_norm > prev_residual * 1.5 && iter > 20) {
            divergence_count++;
            if (divergence_count >= 3) {
                result.status = SolverStatus::DIVERGED;
                result.error_msg = "检测到持续残差增长";
                FEEM_WARN("BiCGSTABSolver+ILU0: 可能发散于第{}次迭代", iter);
                break;
            }
        } else {
            divergence_count = 0;
        }

        if (iter > 10 && residual_norm > prev_residual * 1000.0) {
            result.status = SolverStatus::DIVERGED;
            result.error_msg = "检测到残差异常增长";
            FEEM_WARN("BiCGSTABSolver+ILU0: 可能发散于第{}次迭代", iter);
            break;
        }

        prev_residual = residual_norm;
    }

    if (iter > max_iter) {
        result.status = SolverStatus::MAX_ITER_REACHED;
        FEEM_WARN("BiCGSTABSolver+ILU0: 达到最大迭代次数{}", max_iter);
    }

    result.x = x;
    result.iterations = (iter <= max_iter) ? iter : max_iter;
    result.residual_norm = residual_norm;

    return result;
}

// ==================== ScalarAMG 实现 ====================

ScalarAMG::ScalarAMG(int max_levels, int min_coarse_size)
    : max_levels_(max_levels)
    , min_coarse_size_(min_coarse_size)
    , strong_threshold_(0.25)
    , smoother_iterations_(2)
    , cycle_type_(CycleType::V_CYCLE)
    , smoother_type_(SmootherType::JACOBI)
    , jacobi_damping_(0.8) {
    FEEM_DEBUG("ScalarAMG 构造完成, max_levels={}, min_coarse_size={}",
               max_levels_, min_coarse_size_);
}

void ScalarAMG::compute(const Eigen::SparseMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 清除旧数据
    A_.clear();
    P_.clear();
    R_.clear();
    initialized_ = false;

    // 存储最细层矩阵（确保压缩存储）
    Eigen::SparseMatrix<double> A_compressed = A;
    A_compressed.makeCompressed();
    A_.push_back(A_compressed);

    Eigen::SparseMatrix<double> current_A = A;
    FEEM_INFO("ScalarAMG::compute 开始构建层次结构, 最细网格尺寸: {}x{}", A.rows(), A.cols());

    // 逐层粗化
    for (int level = 0; level < max_levels_ - 1; ++level) {
        // 检查是否达到最小粗网格尺寸
        if (current_A.rows() <= min_coarse_size_) {
            FEEM_DEBUG("AMG 第{}层达到最小尺寸{}, 停止粗化", level, current_A.rows());
            break;
        }

        // 构建插值算子和C/F点分离
        auto [P, coarse_indices] = build_interpolation(current_A);

        int n_fine = current_A.rows();
        int n_coarse = P.cols();

        if (n_coarse <= 0 || n_coarse >= n_fine) {
            FEEM_WARN("AMG 第{}层粗化异常: 细网格{}, 粗网格{}, 停止粗化",
                      level, n_fine, n_coarse);
            break;
        }

        // 存储插值算子
        Eigen::SparseMatrix<double> P_compressed = P;
        P_compressed.makeCompressed();
        P_.push_back(P_compressed);

        // 限制算子 R = P^T
        Eigen::SparseMatrix<double> R = P.transpose().eval();
        R.makeCompressed();
        R_.push_back(R);

        // Galerkin粗网格矩阵: A_coarse = R * A_fine * P
        Eigen::SparseMatrix<double> A_coarse = (R * current_A * P).pruned().eval();
        A_coarse.makeCompressed();
        A_.push_back(A_coarse);

        current_A = A_coarse;

        FEEM_DEBUG("AMG 层 {}: 细网格{} → 粗网格{}, 插值非零元: {}, 粗网格非零元: {}",
                   level, n_fine, n_coarse, P.nonZeros(), A_coarse.nonZeros());
    }

    initialized_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);

    FEEM_INFO("ScalarAMG::compute 完成, 总层数: {}, 耗时: {:.3f}ms",
              A_.size(), duration.count());

    // 打印层次统计信息
    for (size_t i = 0; i < A_.size(); ++i) {
        FEEM_DEBUG("  层 {}: 尺寸={}x{}, 非零元={}", i, A_[i].rows(), A_[i].cols(), A_[i].nonZeros());
    }
}

Eigen::VectorXd ScalarAMG::apply(const Eigen::VectorXd& v) const {
    if (!initialized_) {
        FEEM_ERROR("ScalarAMG::apply 错误: AMG未初始化，请先调用compute()");
        return v;
    }

    if (v.size() != A_[0].rows()) {
        FEEM_ERROR("ScalarAMG::apply 错误: 向量维度{}与矩阵维度{}不匹配",
                   v.size(), A_[0].rows());
        return v;
    }

    Eigen::VectorXd x = Eigen::VectorXd::Zero(v.size());
    apply_v_cycle(0, x, v, 1);

    return x;
}

void ScalarAMG::set_strong_threshold(double theta) {
    strong_threshold_ = theta;
    FEEM_DEBUG("ScalarAMG 强连接阈值设置为: {}", theta);
}

void ScalarAMG::set_smoother_iterations(int iter) {
    smoother_iterations_ = iter;
    FEEM_DEBUG("ScalarAMG 光滑器迭代次数设置为: {}", iter);
}

void ScalarAMG::set_cycle_type(CycleType cycle) {
    cycle_type_ = cycle;
    FEEM_DEBUG("ScalarAMG 循环类型设置为: {}", (cycle == CycleType::V_CYCLE) ? "V-CYCLE" : "W-CYCLE");
}

void ScalarAMG::set_smoother_type(SmootherType smoother) {
    smoother_type_ = smoother;
}

void ScalarAMG::set_jacobi_damping(double omega) {
    jacobi_damping_ = omega;
    FEEM_DEBUG("ScalarAMG Jacobi阻尼因子设置为: {}", omega);
}

int ScalarAMG::get_num_levels() const {
    return static_cast<int>(A_.size());
}

std::vector<int> ScalarAMG::get_level_sizes() const {
    std::vector<int> sizes;
    for (const auto& mat : A_) {
        sizes.push_back(mat.rows());
    }
    return sizes;
}

bool ScalarAMG::is_initialized() const {
    return initialized_;
}

std::pair<Eigen::SparseMatrix<double>, std::vector<int>>
ScalarAMG::build_interpolation(const Eigen::SparseMatrix<double>& A) {
    int n = A.rows();
    std::vector<bool> is_coarse(n, false);
    std::vector<int> coarse_indices;

    // 简化的交替粗化策略：偶数索引选为C点
    for (int i = 0; i < n; ++i) {
        if (i % 2 == 0) {
            is_coarse[i] = true;
            coarse_indices.push_back(i);
        }
    }

    int n_coarse = static_cast<int>(coarse_indices.size());
    std::vector<Eigen::Triplet<double>> triplets;

    // 构建插值算子P（n × n_coarse）
    for (int i = 0; i < n; ++i) {
        if (is_coarse[i]) {
            // 粗网点：单位插值（保持原值）
            auto it = std::find(coarse_indices.begin(), coarse_indices.end(), i);
            int coarse_idx = static_cast<int>(std::distance(coarse_indices.begin(), it));
            triplets.emplace_back(i, coarse_idx, 1.0);
        } else {
            // 细网点：基于邻居粗网点的加权插值
            double diag = A.coeff(i, i);
            if (std::abs(diag) < 1e-15) {
                diag = 1.0;  // 避免除零
            }

            std::vector<std::pair<int, double>> neighbors;

            // 收集强连接的粗网点邻居
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                int j = it.col();
                if (is_coarse[j]) {
                    double weight = std::abs(it.value());
                    neighbors.emplace_back(j, weight);
                }
            }

            // 归一化权重并插入三元组
            for (auto& [j, weight] : neighbors) {
                auto it = std::find(coarse_indices.begin(), coarse_indices.end(), j);
                int coarse_idx = static_cast<int>(std::distance(coarse_indices.begin(), it));
                // 插值权重：-a_ij / a_ii（经典直接插值公式）
                double interp_weight = -A.coeff(i, j) / diag;
                triplets.emplace_back(i, coarse_idx, interp_weight);
            }

            // 如果没有粗网点邻居，使用均匀分配到最近的几个C点
            if (neighbors.empty()) {
                FEEM_WARN("AMG插值: 细网点{}没有粗网点邻居, 使用默认插值", i);
                // 找到最近的粗网点
                int nearest_coarse = coarse_indices[0];
                double min_dist = std::abs(i - nearest_coarse);
                for (int ci : coarse_indices) {
                    double dist = std::abs(i - ci);
                    if (dist < min_dist) {
                        min_dist = dist;
                        nearest_coarse = ci;
                    }
                }
                auto it = std::find(coarse_indices.begin(), coarse_indices.end(), nearest_coarse);
                int coarse_idx = static_cast<int>(std::distance(coarse_indices.begin(), it));
                triplets.emplace_back(i, coarse_idx, 1.0);
            }
        }
    }

    Eigen::SparseMatrix<double> P(n, n_coarse);
    P.setFromTriplets(triplets.begin(), triplets.end());
    P.makeCompressed();

    return {P, coarse_indices};
}

void ScalarAMG::apply_v_cycle(int level, Eigen::VectorXd& x,
                              const Eigen::VectorXd& b, int depth) const {
    if (level >= static_cast<int>(A_.size()) - 1) {
        // 最粗层：使用直接求解器精确求解（转换为稠密矩阵）
        Eigen::MatrixXd A_dense = A_[level];
        x = A_dense.colPivHouseholderQr().solve(b);
        return;
    }

    int half_smoother_iters = std::max(1, smoother_iterations_ / 2);

    // 前光滑：消除当前层的高频误差分量
    for (int i = 0; i < half_smoother_iters; ++i) {
        if (smoother_type_ == SmootherType::JACOBI) {
            x = jacobi_smooth(A_[level], x, b);
        } else {
            gauss_seidel_smooth_forward(A_[level], x, b);
        }
    }

    // 计算残差并限制到粗网格
    Eigen::VectorXd r = b - A_[level] * x;
    Eigen::VectorXd r_coarse = R_[level] * r;

    // 在粗网格上递归求解残差方程 A_c * e_c = r_c
    Eigen::VectorXd e_coarse = Eigen::VectorXd::Zero(r_coarse.size());

    if (cycle_type_ == CycleType::W_CYCLE && depth < 2) {
        // W循环：在中间层访问两次粗网格
        apply_v_cycle(level + 1, e_coarse, r_coarse, depth + 1);
        Eigen::VectorXd temp = e_coarse;
        apply_v_cycle(level + 1, e_coarse, r_coarse - A_[level + 1] * temp, depth + 1);
    } else {
        // V循环：标准单次递归
        apply_v_cycle(level + 1, e_coarse, r_coarse, depth + 1);
    }

    // 插值校正：将粗网格解加回细网格
    x += P_[level] * e_coarse;

    // 后光滑：消除插值引入的高频误差
    for (int i = 0; i < half_smoother_iters; ++i) {
        if (smoother_type_ == SmootherType::JACOBI) {
            x = jacobi_smooth(A_[level], x, b);
        } else {
            gauss_seidel_smooth_backward(A_[level], x, b);
        }
    }
}

Eigen::VectorXd ScalarAMG::jacobi_smooth(const Eigen::SparseMatrix<double>& A,
                                         const Eigen::VectorXd& x,
                                         const Eigen::VectorXd& b) const {
    Eigen::VectorXd x_new = x;
    Eigen::VectorXd residual = b - A * x;

    // Jacobi更新: x_new = x + ω * D^{-1} * residual
    for (int i = 0; i < A.rows(); ++i) {
        double diag = A.coeff(i, i);
        if (std::abs(diag) > 1e-15) {
            x_new(i) += jacobi_damping_ * residual(i) / diag;
        }
    }

    return x_new;
}

void ScalarAMG::gauss_seidel_smooth_forward(const Eigen::SparseMatrix<double>& A,
                                            Eigen::VectorXd& x,
                                            const Eigen::VectorXd& b) const {
    // 前向GS光滑：逐行更新，使用最新值
    for (int i = 0; i < A.rows(); ++i) {
        double diag = A.coeff(i, i);
        if (std::abs(diag) < 1e-15) {
            continue;
        }

        double sum = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            int j = it.col();
            if (j != i) {
                sum += it.value() * x(j);
            }
        }
        x(i) = (b(i) - sum) / diag;
    }
}

void ScalarAMG::gauss_seidel_smooth_backward(const Eigen::SparseMatrix<double>& A,
                                             Eigen::VectorXd& x,
                                             const Eigen::VectorXd& b) const {
    // 后向GS光滑：从最后一行向前逐行更新
    for (int i = A.rows() - 1; i >= 0; --i) {
        double diag = A.coeff(i, i);
        if (std::abs(diag) < 1e-15) {
            continue;
        }

        double sum = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            int j = it.col();
            if (j != i) {
                sum += it.value() * x(j);
            }
        }
        x(i) = (b(i) - sum) / diag;
    }
}

// ==================== 工厂函数实现 ====================

std::unique_ptr<int> create_preconditioner(
    IterativeSolverConfig::PreconditionerType type,
    const Eigen::SparseMatrix<double>& /*A*/) {

    // 当前工厂函数主要用于接口预留
    // 实际的预条件子在CGSolver/BiCGSTABSolver内部根据类型直接创建
    // 未来可扩展为返回统一的预条件子基类指针

    switch (type) {
        case IterativeSolverConfig::PreconditionerType::NONE:
            return nullptr;
        case IterativeSolverConfig::PreconditionerType::JACOBI:
        case IterativeSolverConfig::PreconditionerType::ILU0:
            // 预条件子在求解器内部创建，这里返回空指针作为占位符
            return nullptr;
        default:
            return nullptr;
    }
}

} // namespace numeric
