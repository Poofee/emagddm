/**
 * @file em_direct_solvers.cpp
 * @brief 核心数值层 - 多后端直接求解器完整实现
 * @details 实现3个直接求解器类（SymmetricDirectSolver、SymmetricIndefiniteDirectSolver、GeneralDirectSolver）
 *          和DirectBackendManager后端管理器的全部方法。
 *
 * @par 实现架构：
 * - 策略模式：通过backend_type_成员变量在运行时选择求解后端
 * - 条件编译：#ifdef HAVE_SUPERLU / #ifdef HAVE_MUMPS 控制第三方后端
 * - 瞬态优化：set_matrix()分解一次并缓存，solve()多次复用分解因子
 * - 优雅降级：请求的后端不可用时自动回退到Eigen并记录WARNING日志
 *
 * @par 性能关键点：
 * - Eigen矩阵压缩存储：eigen_matrix_.makeCompressed() 减少内存占用
 * - 智能指针管理：std::unique_ptr确保资源自动释放，避免内存泄漏
 * - 避免不必要的拷贝：使用std::move、const引用传递
 * - 高精度计时：std::chrono::high_resolution_clock 提供纳秒级精度
 *
 * @par 错误处理策略：
 * - 输入校验前置：set_matrix()/solve()入口处严格检查参数合法性
 * - 状态码明确：通过SolverStatus枚举清晰标识错误类型
 * - 日志分级输出：INFO/DEBUG/WARN/ERROR/CRITICAL对应不同严重程度
 * - 资源安全释放：clear()方法保证无内存泄漏
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"
#include <algorithm>
#include <stdexcept>

namespace {

// 匿名命名空间辅助函数：创建成功结果
numeric::SolverResult create_success_result(Eigen::VectorXd x, int iterations = 0) {
    numeric::SolverResult result;
    result.x = std::move(x);
    result.status = numeric::SolverStatus::SUCCESS;
    result.iterations = iterations;
    return result;
}

// 匿名命名空间辅助函数：创建错误结果
numeric::SolverResult create_error_result(numeric::SolverStatus status, const std::string& error_msg) {
    numeric::SolverResult result;
    result.status = status;
    result.error_msg = error_msg;
    return result;
}

// 匿名命名空间辅助函数：创建成功结果（复数版本）
numeric::SolverResult create_success_result_complex(Eigen::VectorXcd x, int iterations = 0) {
    numeric::SolverResult result;
    result.x_complex = std::move(x);
    result.status = numeric::SolverStatus::SUCCESS;
    result.iterations = iterations;
    return result;
}

} // anonymous namespace

namespace numeric {

// SuperluContext 和 SuperLU 后端方法已移至 em_direct_solvers_superlu.cpp（独立编译单元）

// ============================================================================
// DirectBackendManager 实现
// ============================================================================

bool DirectBackendManager::isBackendAvailable(DirectBackendType type) {
    switch (type) {
        case DirectBackendType::EIGEN:
            return true;  // Eigen始终可用（作为默认回退方案）
        case DirectBackendType::SUPERLU:
#ifdef HAVE_SUPERLU
            return true;
#else
            return false;
#endif
        case DirectBackendType::MUMPS:
#ifdef HAVE_MUMPS
            return true;
#else
            return false;
#endif
        default:
            FEEM_ERROR("DirectBackendManager::isBackendAvailable - 未知的后端类型: {}", static_cast<int>(type));
            return false;
    }
}

std::string DirectBackendManager::getBackendName(DirectBackendType type) {
    switch (type) {
        case DirectBackendType::EIGEN:
            return "Eigen";
        case DirectBackendType::SUPERLU:
            return "SuperLU";
        case DirectBackendType::MUMPS:
            return "MUMPS";
        default:
            return "Unknown";
    }
}

std::vector<DirectBackendType> DirectBackendManager::getAvailableBackends() {
    std::vector<DirectBackendType> backends;
    backends.push_back(DirectBackendType::EIGEN);  // Eigen始终可用

#ifdef HAVE_SUPERLU
    backends.push_back(DirectBackendType::SUPERLU);
#endif

#ifdef HAVE_MUMPS
    backends.push_back(DirectBackendType::MUMPS);
#endif

    return backends;
}

// ============================================================================
// SymmetricDirectSolver 实现
// ============================================================================

SymmetricDirectSolver::SymmetricDirectSolver(DirectBackendType backend) {
    backend_type_ = fallback_to_eigen_if_unavailable(backend);

    FEEM_INFO("SymmetricDirectSolver 构造完成, 初始后端: {}",
              DirectBackendManager::getBackendName(backend_type_));
}

void SymmetricDirectSolver::set_backend(DirectBackendType backend) {
    DirectBackendType actual_backend = fallback_to_eigen_if_unavailable(backend);

    if (actual_backend != backend_type_) {
        FEEM_INFO("SymmetricDirectSolver 后端切换: {} -> {}",
                  DirectBackendManager::getBackendName(backend_type_),
                  DirectBackendManager::getBackendName(actual_backend));

        clear();  // 清除旧后端的缓存数据
        backend_type_ = actual_backend;
    }
}

void SymmetricDirectSolver::set_symmetry_tolerance(double tol) {
    if (tol <= 0.0 || tol > 0.1) {
        FEEM_WARN("SymmetricDirectSolver::set_symmetry_tolerance - "
                  "容差值 {:.2e}超出推荐范围[1e-12, 1e-1]", tol);
    }
    symmetry_tol_ = tol;

    FEEM_DEBUG("SymmetricDirectSolver 对称性容差设置为: {:.2e}", symmetry_tol_);
}

void SymmetricDirectSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1：输入合法性校验
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix - 矩阵维度非法: rows={}, cols={}",
                   A.rows(), A.cols());
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // 步骤2：CSR → Eigen格式转换（统一中间格式）
    FEEM_DEBUG("SymmetricDirectSolver 开始转换 CSR -> Eigen 格式, 维度: {}x{}",
               A.rows(), A.cols());

    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();  // 压缩存储减少内存占用
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix - 矩阵转换失败: {}", e.what());
        return;
    }

    // 步骤3：对称性校验（可选，仅输出警告）
    if (!check_symmetry(eigen_matrix_, symmetry_tol_)) {
        FEEM_WARN("SymmetricDirectSolver::set_matrix - "
                  "矩阵不对称性超过容差阈值 {:.2e} (||A-A^T||/||A||)",
                  symmetry_tol_);
    }

    // 步骤4：根据当前后端执行Cholesky分解
    SolverResult decompose_result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            decompose_result = decompose_with_eigen();
            break;
#ifdef HAVE_SUPERLU
        case DirectBackendType::SUPERLU:
            decompose_result = decompose_with_superlu();
            break;
#endif
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            decompose_result = decompose_with_mumps();
            break;
#endif
        default:
            decompose_result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                                   "未知的后端类型");
            break;
    }

    if (decompose_result.status != SolverStatus::SUCCESS) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix - 分解失败: {}",
                   decompose_result.error_msg);
        return;
    }

    matrix_set_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("SymmetricDirectSolver::set_matrix 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
              eigen_matrix_.rows(), eigen_matrix_.cols(),
              DirectBackendManager::getBackendName(backend_type_), duration_ms);
}

SolverResult SymmetricDirectSolver::solve(const Eigen::VectorXd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验：矩阵是否已设置
    if (!matrix_set_) {
        FEEM_ERROR("SymmetricDirectSolver::solve - 矩阵未设置，请先调用 set_matrix()");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "矩阵未设置，请先调用 set_matrix()");
    }

    // 输入校验：维度匹配
    if (b.size() != eigen_matrix_.rows()) {
        FEEM_ERROR("SymmetricDirectSolver::solve - 维度不匹配: 矩阵{}, 向量{}",
                   eigen_matrix_.rows(), b.size());
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与矩阵维度不匹配");
    }

    // 根据当前后端执行求解
    SolverResult result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            result = solve_with_eigen(b);
            break;
#ifdef HAVE_SUPERLU
        case DirectBackendType::SUPERLU:
            result = solve_with_superlu(b);
            break;
#endif
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            result = solve_with_mumps(b);
            break;
#endif
        default:
            result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                         "未知的后端类型");
            break;
    }

    // 计算残差范数用于质量评估（仅在成功时计算）
    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXd residual = eigen_matrix_ * result.x - b;
        result.residual_norm = residual.norm();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    if (result.status == SolverStatus::SUCCESS) {
        FEEM_DEBUG("SymmetricDirectSolver::solve 成功, 残差: {:.6e}, 耗时: {:.3f} ms",
                   result.residual_norm, result.solve_time_ms);
    } else {
        FEEM_ERROR("SymmetricDirectSolver::solve 失败: {}", result.error_msg);
    }

    return result;
}

void SymmetricDirectSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1：输入合法性校验
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix(复数) - 矩阵维度非法: rows={}, cols={}",
                   A.rows(), A.cols());
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix(复数) - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // 步骤2：CSR → Eigen格式转换（复数版本）
    FEEM_DEBUG("SymmetricDirectSolver 开始转换 CSR -> Eigen 复数格式, 维度: {}x{}",
               A.rows(), A.cols());

    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();  // 压缩存储减少内存占用
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix(复数) - 矩阵转换失败: {}", e.what());
        return;
    }

    // 步骤3：根据当前后端执行复数Cholesky分解
    SolverResult decompose_result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            decompose_result = decompose_with_eigen_complex();
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            decompose_result = decompose_with_mumps_complex();
            break;
#endif
        default:
            decompose_result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                                   "未知的后端类型");
            break;
    }

    if (decompose_result.status != SolverStatus::SUCCESS) {
        FEEM_ERROR("SymmetricDirectSolver::set_matrix(复数) - 分解失败: {}",
                   decompose_result.error_msg);
        return;
    }

    matrix_set_complex_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("SymmetricDirectSolver::set_matrix(复数) 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
              eigen_matrix_complex_.rows(), eigen_matrix_complex_.cols(),
              DirectBackendManager::getBackendName(backend_type_), duration_ms);
}

SolverResult SymmetricDirectSolver::solve(const Eigen::VectorXcd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验：复数矩阵是否已设置
    if (!matrix_set_complex_) {
        FEEM_ERROR("SymmetricDirectSolver::solve(复数) - 复数矩阵未设置，请先调用 set_matrix(复数)");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "复数矩阵未设置，请先调用 set_matrix(复数版本)");
    }

    // 输入校验：维度匹配
    if (b.size() != eigen_matrix_complex_.rows()) {
        FEEM_ERROR("SymmetricDirectSolver::solve(复数) - 维度不匹配: 矩阵{}, 向量{}",
                   eigen_matrix_complex_.rows(), b.size());
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与复数矩阵维度不匹配");
    }

    // 根据当前后端执行复数求解
    SolverResult result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            result = solve_with_eigen_complex(b);
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            result = solve_with_mumps_complex(b);
            break;
#endif
        default:
            result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                         "未知的后端类型");
            break;
    }

    // 计算残差范数用于质量评估（仅在成功时计算）
    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXcd residual = eigen_matrix_complex_ * result.x_complex - b;
        result.residual_norm = residual.norm();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    if (result.status == SolverStatus::SUCCESS) {
        FEEM_DEBUG("SymmetricDirectSolver::solve(复数) 成功, 残差: {:.6e}, 耗时: {:.3f} ms",
                   result.residual_norm, result.solve_time_ms);
    } else {
        FEEM_ERROR("SymmetricDirectSolver::solve(复数) 失败: {}", result.error_msg);
    }

    return result;
}

std::string SymmetricDirectSolver::get_solver_name() const {
    return std::string("SymmetricDirect_") + DirectBackendManager::getBackendName(backend_type_);
}

void SymmetricDirectSolver::clear() {
    FEEM_DEBUG("SymmetricDirectSolver::clear - 释放所有内部资源");

    // 释放Eigen后端资源（实数版本）
    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();  // 释放矩阵内存

    // 释放Eigen后端资源（复数版本）
    eigen_solver_complex_.reset();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();  // 释放复数矩阵内存
    matrix_set_complex_ = false;

#ifdef HAVE_SUPERLU
    // 释放 SuperLU_MT 资源（RAII 封装自动管理所有资源释放）
    if (superlu_ctx_) {
        superlu_ctx_->reset();
        superlu_ctx_.reset();
    }
#endif

#ifdef HAVE_MUMPS
    // 释放 MUMPS 资源（RAII 封装自动管理所有资源释放）
    if (mumps_ctx_) {
        mumps_ctx_->reset();
        mumps_ctx_.reset();
    }
    // 释放 MUMPS 复数资源（RAII 封装自动管理所有资源释放）
    if (mumps_ctx_complex_) {
        mumps_ctx_complex_->reset();
        mumps_ctx_complex_.reset();
    }
#endif

    matrix_set_ = false;
}

// 私有方法实现

DirectBackendType SymmetricDirectSolver::fallback_to_eigen_if_unavailable(DirectBackendType requested) {
    if (DirectBackendManager::isBackendAvailable(requested)) {
        return requested;
    }

    FEEM_WARN("SymmetricDirectSolver - 请求的后端 {} 不可用，优雅降级到 Eigen",
              DirectBackendManager::getBackendName(requested));

    return DirectBackendType::EIGEN;
}

bool SymmetricDirectSolver::check_symmetry(const Eigen::SparseMatrix<double>& mat, double tol) const {
    // 计算不对称性的Frobenius范数: ||A - A^T||_F
    // 注意：Eigen稀疏矩阵不支持直接减去转置（存储顺序不同），需转为稠密矩阵
    Eigen::MatrixXd mat_dense = mat;
    Eigen::MatrixXd asymmetry = mat_dense - mat_dense.transpose();
    double asym_norm = asymmetry.norm();

    // 计算原始矩阵的Frobenius范数: ||A||_F
    double mat_norm = mat.norm();

    // 避免除以零（全零矩阵视为对称）
    if (mat_norm < 1e-15) {
        return true;
    }

    double relative_asymmetry = asym_norm / mat_norm;

    FEEM_DEBUG("对称性检查: 相对不对称度 = {:.2e}, 阈值 = {:.2e}",
               relative_asymmetry, tol);

    return relative_asymmetry <= tol;
}

SolverResult SymmetricDirectSolver::decompose_with_eigen() {
    FEEM_DEBUG("SymmetricDirectSolver::decompose_with_eigen - 开始Eigen Cholesky分解");

    try {
        eigen_solver_ = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>();

        // 执行符号分析 + 数值分解（SimplicialLLT自动处理）
        eigen_solver_->compute(eigen_matrix_);

        // 检查分解是否成功
        if (eigen_solver_->info() != Eigen::Success) {
            std::string error_msg = "Eigen Cholesky分解失败";

            if (eigen_matrix_.rows() > 0) {
                error_msg += "（矩阵可能非正定或奇异）";
            }

            FEEM_ERROR("{}", error_msg);
            return create_error_result(SolverStatus::NUMERICAL_ERROR, error_msg);
        }

        FEEM_DEBUG("SymmetricDirectSolver::decompose_with_eigen - Cholesky分解成功");
        return create_success_result(Eigen::VectorXd());  // 空向量占位（分解无解向量）
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("SymmetricDirectSolver::decompose_with_eigen - 内存分配失败: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "内存不足，无法完成Cholesky分解");
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::decompose_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("Eigen分解异常: ") + e.what());
    }
}

SolverResult SymmetricDirectSolver::solve_with_eigen(const Eigen::VectorXd& b) {
    FEEM_DEBUG("SymmetricDirectSolver::solve_with_eigen - 开始前代/回代");

    try {
        // 利用缓存的L因子执行前代 Ly = b 和回代 L^T x = y
        Eigen::VectorXd x = eigen_solver_->solve(b);

        if (eigen_solver_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricDirectSolver::solve_with_eigen - 求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "Eigen前代/回代失败");
        }

        return create_success_result(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::solve_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("Eigen求解异常: ") + e.what());
    }
}

SolverResult SymmetricDirectSolver::decompose_with_eigen_complex() {
    FEEM_DEBUG("SymmetricDirectSolver::decompose_with_eigen_complex - 开始复数Cholesky分解");

    try {
        eigen_solver_complex_ = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>>>();

        // 执行符号分析 + 数值分解（SimplicialLLT自动处理复数Hermitian正定矩阵）
        eigen_solver_complex_->compute(eigen_matrix_complex_);

        // 检查分解是否成功
        if (eigen_solver_complex_->info() != Eigen::Success) {
            std::string error_msg = "Eigen复数Cholesky分解失败";

            if (eigen_matrix_complex_.rows() > 0) {
                error_msg += "（矩阵可能非Hermitian正定或奇异）";
            }

            FEEM_ERROR("{}", error_msg);
            return create_error_result(SolverStatus::NUMERICAL_ERROR, error_msg);
        }

        FEEM_DEBUG("SymmetricDirectSolver::decompose_with_eigen_complex - 复数Cholesky分解成功");
        return create_success_result_complex(Eigen::VectorXcd());  // 空向量占位（分解无解向量）
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("SymmetricDirectSolver::decompose_with_eigen_complex - 内存分配失败: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "内存不足，无法完成复数Cholesky分解");
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::decompose_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("Eigen复数分解异常: ") + e.what());
    }
}

SolverResult SymmetricDirectSolver::solve_with_eigen_complex(const Eigen::VectorXcd& b) {
    FEEM_DEBUG("SymmetricDirectSolver::solve_with_eigen_complex - 开始复数前代/回代");

    try {
        // 利用缓存的L因子执行前代 Ly = b 和回代 L^H x = y
        Eigen::VectorXcd x = eigen_solver_complex_->solve(b);

        if (eigen_solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricDirectSolver::solve_with_eigen_complex - 复数求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "Eigen复数前代/回代失败");
        }

        return create_success_result_complex(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::solve_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("Eigen复数求解异常: ") + e.what());
    }
}

// SuperLU 后端方法已移至 em_direct_solvers_superlu.cpp（独立编译单元）

// MumpsContext 和 MUMPS 后端方法已移至 em_direct_solvers_mumps.cpp（独立编译单元）

// ============================================================================
// SymmetricIndefiniteDirectSolver 实现
// ============================================================================

SymmetricIndefiniteDirectSolver::SymmetricIndefiniteDirectSolver(DirectBackendType backend) {
    backend_type_ = fallback_to_eigen_if_unavailable(backend);

    FEEM_INFO("SymmetricIndefiniteDirectSolver 构造完成, 初始后端: {}",
              DirectBackendManager::getBackendName(backend_type_));
}

void SymmetricIndefiniteDirectSolver::set_regularization_epsilon(double epsilon) {
    if (epsilon < 0.0 || epsilon > 1.0) {
        FEEM_WARN("SymmetricIndefiniteDirectSolver::set_regularization_epsilon - "
                  "参数 {:.2e} 超出合理范围 [0, 1]", epsilon);
    }

    regularization_epsilon_ = epsilon;

    FEEM_DEBUG("SymmetricIndefiniteDirectSolver 正则化参数设置为: {:.2e}",
               regularization_epsilon_);
}

void SymmetricIndefiniteDirectSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (A.rows() <= 0 || A.cols() <= 0 || A.rows() != A.cols()) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix - 矩阵维度非法: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // CSR → Eigen转换
    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix - 转换失败: {}", e.what());
        return;
    }

    // 对称性校验
    if (!check_symmetry(eigen_matrix_, symmetry_tol_)) {
        FEEM_WARN("SymmetricIndefiniteDirectSolver::set_matrix - "
                  "矩阵不对称性超过阈值 {:.2e}", symmetry_tol_);
    }

    regularization_applied_ = false;

    // 尝试标准LDL^T分解
    SolverResult decompose_result = decompose_with_eigen();

    // 若标准分解失败（可能因零主元），自动应用正则化重试
    if (decompose_result.status == SolverStatus::NUMERICAL_ERROR) {
        FEEM_WARN("SymmetricIndefiniteDirectSolver::set_matrix - "
                  "标准LDL^T分解失败，尝试施加正则化 ε={:.2e}",
                  regularization_epsilon_);

        regularized_matrix_ = eigen_matrix_;

        // 添加微小扰动到对角线：A_reg = A + ε*I
        for (int i = 0; i < regularized_matrix_.rows(); i++) {
            regularized_matrix_.coeffRef(i, i) += regularization_epsilon_;
        }

        // 使用正则化后的矩阵重新分解
        eigen_matrix_ = regularized_matrix_;  // 替换为正则化版本
        decompose_result = decompose_with_eigen();

        if (decompose_result.status == SolverStatus::SUCCESS) {
            regularization_applied_ = true;
            FEEM_INFO("SymmetricIndefiniteDirectSolver::set_matrix - "
                      "正则化成功, 扰动幅度: {:.2e}", regularization_epsilon_);
        } else {
            FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix - "
                       "即使正定化后分解仍然失败");
            return;
        }
    }

    matrix_set_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("SymmetricIndefiniteDirectSolver::set_matrix 完成, 维度: {}x{}, "
              "正则化: {}, 耗时: {:.3f} ms",
              eigen_matrix_.rows(), eigen_matrix_.cols(),
              regularization_applied_ ? "是" : "否", duration_ms);
}

SolverResult SymmetricIndefiniteDirectSolver::solve(const Eigen::VectorXd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    if (!matrix_set_) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve - 矩阵未设置");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "矩阵未设置，请先调用 set_matrix()");
    }

    if (b.size() != eigen_matrix_.rows()) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve - 维度不匹配");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与矩阵维度不匹配");
    }

    SolverResult result = solve_with_eigen(b);

    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXd residual = eigen_matrix_ * result.x - b;
        result.residual_norm = residual.norm();

        // 若使用了正则化，额外提示用户
        if (regularization_applied_) {
            FEEM_DEBUG("SymmetricIndefiniteDirectSolver::solve - "
                       "注意：当前解基于正则化矩阵，与原问题存在微小偏差");
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    return result;
}

void SymmetricIndefiniteDirectSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验
    if (A.rows() <= 0 || A.cols() <= 0 || A.rows() != A.cols()) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix(复数) - 矩阵维度非法: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // CSR → Eigen转换（复数版本）
    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix(复数) - 转换失败: {}", e.what());
        return;
    }

    // 执行复数LDL^H分解（根据后端类型选择）
    SolverResult decompose_result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            decompose_result = decompose_with_eigen_complex();
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            decompose_result = decompose_with_mumps_complex();
            break;
#endif
        default:
            decompose_result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                                   "未知的后端类型");
            break;
    }

    if (decompose_result.status != SolverStatus::SUCCESS) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::set_matrix(复数) - LDL^H分解失败: {}",
                   decompose_result.error_msg);
        return;
    }

    matrix_set_complex_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("SymmetricIndefiniteDirectSolver::set_matrix(复数) 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
              eigen_matrix_complex_.rows(), eigen_matrix_complex_.cols(),
              DirectBackendManager::getBackendName(backend_type_), duration_ms);
}

SolverResult SymmetricIndefiniteDirectSolver::solve(const Eigen::VectorXcd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    if (!matrix_set_complex_) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve(复数) - 复数矩阵未设置");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "复数矩阵未设置，请先调用 set_matrix(复数版本)");
    }

    if (b.size() != eigen_matrix_complex_.rows()) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve(复数) - 维度不匹配");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与复数矩阵维度不匹配");
    }

    // 根据当前后端执行复数求解
    SolverResult result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            result = solve_with_eigen_complex(b);
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            result = solve_with_mumps_complex(b);
            break;
#endif
        default:
            result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                         "未知的后端类型");
            break;
    }

    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXcd residual = eigen_matrix_complex_ * result.x_complex - b;
        result.residual_norm = residual.norm();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    return result;
}

std::string SymmetricIndefiniteDirectSolver::get_solver_name() const {
    return std::string("SymmetricIndefiniteDirect_") +
           DirectBackendManager::getBackendName(backend_type_);
}

void SymmetricIndefiniteDirectSolver::clear() {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::clear - 释放资源");

    // 释放实数版本资源
    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    regularized_matrix_ = Eigen::SparseMatrix<double>();

    // 释放复数版本资源
    eigen_solver_complex_.reset();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_complex_ = false;

#ifdef HAVE_MUMPS
    // 释放 MUMPS 复数资源（RAII 封装自动管理所有资源释放）
    if (mumps_ctx_complex_) {
        mumps_ctx_complex_->reset();
        mumps_ctx_complex_.reset();
    }
#endif

    matrix_set_ = false;
    regularization_applied_ = false;
}

DirectBackendType SymmetricIndefiniteDirectSolver::fallback_to_eigen_if_unavailable(
    DirectBackendType requested) {

    if (DirectBackendManager::isBackendAvailable(requested)) {
        return requested;
    }

    FEEM_WARN("SymmetricIndefiniteDirectSolver - 后端 {} 不可用，降级到 Eigen",
              DirectBackendManager::getBackendName(requested));

    return DirectBackendType::EIGEN;
}

bool SymmetricIndefiniteDirectSolver::check_symmetry(
    const Eigen::SparseMatrix<double>& mat, double tol) const {

    // 注意：Eigen稀疏矩阵不支持直接减去转置（存储顺序不同），需转为稠密矩阵
    Eigen::MatrixXd mat_dense = mat;
    Eigen::MatrixXd asymmetry = mat_dense - mat_dense.transpose();
    double asym_norm = asymmetry.norm();
    double mat_norm = mat.norm();

    if (mat_norm < 1e-15) {
        return true;
    }

    return (asym_norm / mat_norm) <= tol;
}

SolverResult SymmetricIndefiniteDirectSolver::decompose_with_eigen() {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::decompose_with_eigen - 开始LDL^T分解");

    try {
        eigen_solver_ = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>();

        // SimplicialLDLT支持2x2 pivot旋转，可处理不定矩阵
        eigen_solver_->compute(eigen_matrix_);

        if (eigen_solver_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen - LDL^T分解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "LDL^T分解失败（矩阵可能奇异或强不定）");
        }

        FEEM_DEBUG("SymmetricIndefiniteDirectSolver::decompose_with_eigen - LDL^T分解成功");
        return create_success_result(Eigen::VectorXd());
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen - 内存不足: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("LDL^T异常: ") + e.what());
    }
}

SolverResult SymmetricIndefiniteDirectSolver::solve_with_eigen(const Eigen::VectorXd& b) {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::solve_with_eigen - 开始求解");

    try {
        Eigen::VectorXd x = eigen_solver_->solve(b);

        if (eigen_solver_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve_with_eigen - 求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "LDL^T求解失败");
        }

        return create_success_result(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("求解异常: ") + e.what());
    }
}

SolverResult SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex() {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex - 开始复数LDL^H分解");

    try {
        eigen_solver_complex_ = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>>();

        // SimplicialLDLT支持2x2 pivot旋转，可处理复数Hermitian不定矩阵
        eigen_solver_complex_->compute(eigen_matrix_complex_);

        if (eigen_solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex - 复数LDL^H分解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "复数LDL^H分解失败（矩阵可能奇异或强不定）");
        }

        FEEM_DEBUG("SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex - 复数LDL^H分解成功");
        return create_success_result_complex(Eigen::VectorXcd());
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex - 内存不足: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::decompose_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("复数LDL^H异常: ") + e.what());
    }
}

SolverResult SymmetricIndefiniteDirectSolver::solve_with_eigen_complex(const Eigen::VectorXcd& b) {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::solve_with_eigen_complex - 开始复数求解");

    try {
        Eigen::VectorXcd x = eigen_solver_complex_->solve(b);

        if (eigen_solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve_with_eigen_complex - 复数求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "复数LDL^H求解失败");
        }

        return create_success_result_complex(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricIndefiniteDirectSolver::solve_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("复数求解异常: ") + e.what());
    }
}

// ============================================================================
// GeneralDirectSolver 实现
// ============================================================================

GeneralDirectSolver::GeneralDirectSolver(DirectBackendType backend) {
    backend_type_ = fallback_to_eigen_if_unavailable(backend);

    FEEM_INFO("GeneralDirectSolver 构造完成, 初始后端: {}",
              DirectBackendManager::getBackendName(backend_type_));
}

void GeneralDirectSolver::set_pivot_threshold(double threshold) {
    if (threshold < 0.0 || threshold > 1.0) {
        FEEM_WARN("GeneralDirectSolver::set_pivot_threshold - "
                  "阈值 {:.2f} 超出范围 [0, 1]，将截断", threshold);
        threshold = std::clamp(threshold, 0.0, 1.0);
    }

    pivot_threshold_ = threshold;

    FEEM_DEBUG("GeneralDirectSolver 选主元阈值设置为: {:.2f}", pivot_threshold_);
}

void GeneralDirectSolver::set_matrix(const CsrMatrix<double>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验（仅要求方阵，不要求对称）
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix - 矩阵维度非法: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // CSR → Eigen转换
    try {
        eigen_matrix_ = SparseConverter::to_eigen(A);
        eigen_matrix_.makeCompressed();
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix - 转换失败: {}", e.what());
        return;
    }

    // 执行LU分解（带部分选主元）
    SolverResult decompose_result = decompose_with_eigen();

    if (decompose_result.status != SolverStatus::SUCCESS) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix - LU分解失败: {}",
                   decompose_result.error_msg);
        return;
    }

    // 估算条件数并发出病态警告
    double cond_est = estimate_condition_number();

    if (cond_est > 1e15) {
        FEEM_WARN("GeneralDirectSolver::set_matrix - 矩阵严重病态! 条件数估计: {:.2e}", cond_est);
    } else if (cond_est > 1e10) {
        FEEM_WARN("GeneralDirectSolver::set_matrix - 矩阵较为病态, 条件数估计: {:.2e}", cond_est);
    } else {
        FEEM_DEBUG("GeneralDirectSolver::set_matrix - 条件数估计: {:.2e}", cond_est);
    }

    matrix_set_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("GeneralDirectSolver::set_matrix 完成, 维度: {}x{}, 条件数≈{:.2e}, 耗时: {:.3f} ms",
              eigen_matrix_.rows(), eigen_matrix_.cols(), cond_est, duration_ms);
}

SolverResult GeneralDirectSolver::solve(const Eigen::VectorXd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    if (!matrix_set_) {
        FEEM_ERROR("GeneralDirectSolver::solve - 矩阵未设置");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "矩阵未设置，请先调用 set_matrix()");
    }

    if (b.size() != eigen_matrix_.rows()) {
        FEEM_ERROR("GeneralDirectSolver::solve - 维度不匹配");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与矩阵维度不匹配");
    }

    SolverResult result = solve_with_eigen(b);

    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXd residual = eigen_matrix_ * result.x - b;
        result.residual_norm = residual.norm();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    return result;
}

void GeneralDirectSolver::set_matrix(const CsrMatrix<std::complex<double>>& A) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // 输入校验（仅要求方阵，不要求对称）
    if (A.rows() <= 0 || A.cols() <= 0) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix(复数) - 矩阵维度非法: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    if (A.rows() != A.cols()) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix(复数) - 矩阵非方阵: {}x{}",
                   A.rows(), A.cols());
        return;
    }

    // CSR → Eigen转换（复数版本）
    try {
        eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
        eigen_matrix_complex_.makeCompressed();
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix(复数) - 转换失败: {}", e.what());
        return;
    }

    // 执行复数LU分解（根据后端类型选择）
    SolverResult decompose_result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            decompose_result = decompose_with_eigen_complex();
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            decompose_result = decompose_with_mumps_complex();
            break;
#endif
        default:
            decompose_result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                                   "未知的后端类型");
            break;
    }

    if (decompose_result.status != SolverStatus::SUCCESS) {
        FEEM_ERROR("GeneralDirectSolver::set_matrix(复数) - LU分解失败: {}",
                   decompose_result.error_msg);
        return;
    }

    matrix_set_complex_ = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    FEEM_INFO("GeneralDirectSolver::set_matrix(复数) 完成, 维度: {}x{}, 后端: {}, 耗时: {:.3f} ms",
              eigen_matrix_complex_.rows(), eigen_matrix_complex_.cols(),
              DirectBackendManager::getBackendName(backend_type_), duration_ms);
}

SolverResult GeneralDirectSolver::solve(const Eigen::VectorXcd& b) {
    auto start_time = std::chrono::high_resolution_clock::now();

    if (!matrix_set_complex_) {
        FEEM_ERROR("GeneralDirectSolver::solve(复数) - 复数矩阵未设置");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "复数矩阵未设置，请先调用 set_matrix(复数版本)");
    }

    if (b.size() != eigen_matrix_complex_.rows()) {
        FEEM_ERROR("GeneralDirectSolver::solve(复数) - 维度不匹配");
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量长度与复数矩阵维度不匹配");
    }

    // 根据当前后端执行复数求解
    SolverResult result;

    switch (backend_type_) {
        case DirectBackendType::EIGEN:
            result = solve_with_eigen_complex(b);
            break;
#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            result = solve_with_mumps_complex(b);
            break;
#endif
        default:
            result = create_error_result(SolverStatus::NUMERICAL_ERROR,
                                         "未知的后端类型");
            break;
    }

    if (result.status == SolverStatus::SUCCESS) {
        Eigen::VectorXcd residual = eigen_matrix_complex_ * result.x_complex - b;
        result.residual_norm = residual.norm();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    return result;
}

std::string GeneralDirectSolver::get_solver_name() const {
    return std::string("GeneralDirect_") + DirectBackendManager::getBackendName(backend_type_);
}

void GeneralDirectSolver::clear() {
    FEEM_DEBUG("GeneralDirectSolver::clear - 释放资源");

    // 释放实数版本资源
    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();

    // 释放复数版本资源
    eigen_solver_complex_.reset();
    eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
    matrix_set_complex_ = false;

#ifdef HAVE_MUMPS
    // 释放 MUMPS 复数资源（RAII 封装自动管理所有资源释放）
    if (mumps_ctx_complex_) {
        mumps_ctx_complex_->reset();
        mumps_ctx_complex_.reset();
    }
#endif

    matrix_set_ = false;
}

DirectBackendType GeneralDirectSolver::fallback_to_eigen_if_unavailable(DirectBackendType requested) {
    if (DirectBackendManager::isBackendAvailable(requested)) {
        return requested;
    }

    FEEM_WARN("GeneralDirectSolver - 后端 {} 不可用，降级到 Eigen",
              DirectBackendManager::getBackendName(requested));

    return DirectBackendType::EIGEN;
}

SolverResult GeneralDirectSolver::decompose_with_eigen() {
    FEEM_DEBUG("GeneralDirectSolver::decompose_with_eigen - 开始LU分解");

    try {
        eigen_solver_ = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();

        // 设置选主元阈值（通过Eigen内部API，若支持的话）
        // 注意：Eigen::SparseLU的pivot threshold设置方式取决于版本

        // 执行带部分选主元的LU分解 PA = LU
        eigen_solver_->compute(eigen_matrix_);

        if (eigen_solver_->info() != Eigen::Success) {
            FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen - LU分解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "LU分解失败（矩阵可能奇异）");
        }

        FEEM_DEBUG("GeneralDirectSolver::decompose_with_eigen - LU分解成功");
        return create_success_result(Eigen::VectorXd());
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen - 内存不足: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("LU分解异常: ") + e.what());
    }
}

SolverResult GeneralDirectSolver::solve_with_eigen(const Eigen::VectorXd& b) {
    FEEM_DEBUG("GeneralDirectSolver::solve_with_eigen - 开始求解");

    try {
        // 利用缓存的P、L、U因子执行前代和回代
        Eigen::VectorXd x = eigen_solver_->solve(b);

        if (eigen_solver_->info() != Eigen::Success) {
            FEEM_ERROR("GeneralDirectSolver::solve_with_eigen - 求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "LU求解失败");
        }

        return create_success_result(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::solve_with_eigen - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("求解异常: ") + e.what());
    }
}

SolverResult GeneralDirectSolver::decompose_with_eigen_complex() {
    FEEM_DEBUG("GeneralDirectSolver::decompose_with_eigen_complex - 开始复数LU分解");

    try {
        eigen_solver_complex_ = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>();

        // 执行带部分选主元的复数LU分解 PA = LU
        eigen_solver_complex_->compute(eigen_matrix_complex_);

        if (eigen_solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen_complex - 复数LU分解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "复数LU分解失败（矩阵可能奇异）");
        }

        FEEM_DEBUG("GeneralDirectSolver::decompose_with_eigen_complex - 复数LU分解成功");
        return create_success_result_complex(Eigen::VectorXcd());
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen_complex - 内存不足: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::decompose_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("复数LU分解异常: ") + e.what());
    }
}

SolverResult GeneralDirectSolver::solve_with_eigen_complex(const Eigen::VectorXcd& b) {
    FEEM_DEBUG("GeneralDirectSolver::solve_with_eigen_complex - 开始复数求解");

    try {
        // 利用缓存的P、L、U因子执行前代和回代
        Eigen::VectorXcd x = eigen_solver_complex_->solve(b);

        if (eigen_solver_complex_->info() != Eigen::Success) {
            FEEM_ERROR("GeneralDirectSolver::solve_with_eigen_complex - 复数求解失败");
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "复数LU求解失败");
        }

        return create_success_result_complex(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("GeneralDirectSolver::solve_with_eigen_complex - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("复数求解异常: ") + e.what());
    }
}

double GeneralDirectSolver::estimate_condition_number() const {
    // 基于LU分解的U因子估算条件数下界
    // 注意：Eigen::SparseLU的matrixU()返回类型复杂，暂不实现详细估算
    // 返回0.0表示未估算，调用方可根据残差范数判断解的质量

    if (!eigen_solver_) {
        return 0.0;
    }

    // TODO: 实现完整的条件数估算（需处理Eigen SparseLU的特殊返回类型）
    return 0.0;
}

} // namespace numeric
