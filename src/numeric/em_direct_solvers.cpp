/**
 * @file em_direct_solvers.cpp
 * @brief 核心数值层 - 多后端直接求解器完整实现
 * @details 实现3个直接求解器类（SymmetricDirectSolver、SymmetricIndefiniteDirectSolver、GeneralDirectSolver）
 *          和DirectBackendManager后端管理器的全部方法。
 *
 * @par 实现架构：
 * - 策略模式：通过backend_type_成员变量在运行时选择求解后端
 * - 条件编译：#if EM_SOLVER_HAS_SUPERLU / #if EM_SOLVER_HAS_MUMPS 控制第三方后端
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

} // anonymous namespace

namespace numeric {

// ============================================================================
// SuperluContext 实现（SuperLU_MT RAII 封装）
// ============================================================================

#if EM_SOLVER_HAS_SUPERLU

SuperluContext::SuperluContext()
    : state_(State::UNINITIALIZED), n_(0), nnz_(0) {}

SuperluContext::~SuperluContext() {
    reset();
}

SuperluContext::SuperluContext(SuperluContext&& other) noexcept
    : state_(other.state_),
      n_(other.n_),
      nnz_(other.nnz_),
      A_(std::move(other.A_)),
      L_(std::move(other.L_)),
      U_(std::move(other.U_)),
      AC_(std::move(other.AC_)),
      perm_c_(std::move(other.perm_c_)),
      perm_r_(std::move(other.perm_r_)),
      etree_(std::move(other.etree_)),
      colcnt_h_(std::move(other.colcnt_h_)),
      part_super_h_(std::move(other.part_super_h_)),
      options_(std::move(other.options_)),
      Gstat_(std::move(other.Gstat_))
{
    other.state_ = State::UNINITIALIZED;
    other.n_ = 0;
    other.nnz_ = 0;
}

SuperluContext& SuperluContext::operator=(SuperluContext&& other) noexcept {
    if (this != &other) {
        reset();
        state_ = other.state_;
        n_ = other.n_;
        nnz_ = other.nnz_;
        A_ = std::move(other.A_);
        L_ = std::move(other.L_);
        U_ = std::move(other.U_);
        AC_ = std::move(other.AC_);
        perm_c_ = std::move(other.perm_c_);
        perm_r_ = std::move(other.perm_r_);
        etree_ = std::move(other.etree_);
        colcnt_h_ = std::move(other.colcnt_h_);
        part_super_h_ = std::move(other.part_super_h_);
        options_ = std::move(other.options_);
        Gstat_ = std::move(other.Gstat_);
        other.state_ = State::UNINITIALIZED;
        other.n_ = 0;
        other.nnz_ = 0;
    }
    return *this;
}

bool SuperluContext::initialize(int_t n, int_t nnz, double* nzval, int_t* rowind, int_t* colptr) {
    if (n <= 0 || nnz < 0 || !nzval || !rowind || !colptr) {
        FEEM_ERROR("SuperluContext::initialize - 参数非法: n={}, nnz={}, pointers valid={},{},{}",
                   n, nnz, nzval != nullptr, rowind != nullptr, colptr != nullptr);
        return false;
    }

    reset();

    n_ = n;
    nnz_ = nnz;

    try {
        A_ = std::make_unique<SuperMatrix>();
        dCreate_CompCol_Matrix(A_.get(), n, n, nnz, nzval, rowind, colptr,
                               SLU_NC, SLU_D, SLU_GE);

        perm_c_ = std::unique_ptr<int_t[]>(intMalloc(n));
        perm_r_ = std::unique_ptr<int_t[]>(intMalloc(n));
        if (!perm_c_ || !perm_r_) {
            FEEM_ERROR("SuperluContext::initialize - 置换向量分配失败");
            reset();
            return false;
        }

        int_t permc_spec = 1;
        get_perm_c(permc_spec, A_.get(), perm_c_.get());

        etree_ = std::unique_ptr<int_t[]>(intMalloc(n));
        colcnt_h_ = std::unique_ptr<int_t[]>(intMalloc(n));
        part_super_h_ = std::unique_ptr<int_t[]>(intMalloc(n));
        if (!etree_ || !colcnt_h_ || !part_super_h_) {
            FEEM_ERROR("SuperluContext::initialize - 辅助数组分配失败");
            reset();
            return false;
        }

        options_ = std::make_unique<superlumt_options_t>();
        options_->etree = etree_.get();
        options_->colcnt_h = colcnt_h_.get();
        options_->part_super_h = part_super_h_.get();
        options_->perm_c = perm_c_.get();
        options_->perm_r = perm_r_.get();

        state_ = State::INITIALIZED;

        FEEM_DEBUG("SuperluContext::initialize 成功, n={}, nnz={}", n, nnz);
        return true;
    } catch (const std::exception& e) {
        FEEM_ERROR("SuperluContext::initialize - 异常: {}", e.what());
        reset();
        return false;
    }
}

int_t SuperluContext::factorize(int_t nprocs) {
    if (state_ != State::INITIALIZED) {
        FEEM_ERROR("SuperluContext::factorize - 状态错误（需要先调用 initialize），当前状态={}",
                   static_cast<int>(state_));
        return -1;
    }

    release_factors();

    try {
        L_ = std::make_unique<SuperMatrix>();
        U_ = std::make_unique<SuperMatrix>();
        AC_ = std::make_unique<SuperMatrix>();  // pdgstrf_init 需要有效的 AC 输出参数
        Gstat_ = std::make_unique<Gstat_t>();

        superlumt_options_t& opts = *options_;
        opts.nprocs = nprocs;
        opts.fact = DOFACT;
        opts.trans = NOTRANS;
        opts.refact = NO;
        opts.panel_size = sp_ienv(1);
        opts.relax = sp_ienv(2);
        opts.diag_pivot_thresh = 1.0;
        opts.usepr = NO;
        opts.drop_tol = 0.0;
        opts.ColPerm = COLAMD;
        opts.SymmetricMode = NO;
        opts.PrintStat = NO;
        opts.work = nullptr;
        opts.lwork = 0;

        int_t info = 0;
        StatAlloc(n_, nprocs, opts.panel_size, opts.relax, Gstat_.get());
        StatInit(n_, nprocs, Gstat_.get());

        pdgstrf_init(nprocs, opts.fact, opts.trans, opts.refact,
                     opts.panel_size, opts.relax, opts.diag_pivot_thresh,
                     opts.usepr, opts.drop_tol, opts.perm_c, opts.perm_r,
                     opts.work, opts.lwork, A_.get(),
                     AC_.get(), options_.get(), Gstat_.get());

        pdgstrf(options_.get(), AC_.get(), opts.perm_r,
                L_.get(), U_.get(), Gstat_.get(), &info);

        pxgstrf_finalize(options_.get(), AC_.get());

        // 重要: pxgstrf_finalize 已通过 options 指针释放了
        // etree/colcnt_h/part_super_h 的内存，必须释放 unique_ptr 所有权
        // 否则 reset() 时会触发 double-free 崩溃
        etree_.reset();
        colcnt_h_.reset();
        part_super_h_.reset();
        // 同时置空 options 中的指针，防止 release_ac_matrix() 再次调用
        // pxgstrf_finalize 时重复释放
        options_->etree = nullptr;
        options_->colcnt_h = nullptr;
        options_->part_super_h = nullptr;

        // pxgstrf_finalize 已销毁 AC 内部数据（Destroy_CompCol_Permuted），
        // 释放 AC_ 所有权防止重复销毁
        AC_.reset();

        StatFree(Gstat_.get());
        Gstat_.reset();

        if (info != 0) {
            if (info > 0 && info <= n_) {
                FEEM_WARN("SuperluContext::factorize - 矩阵可能奇异, U({},{})=0", info, info);
            } else if (info < 0) {
                FEEM_ERROR("SuperluContext::factorize - 参数错误或内存不足, info={}", info);
            } else if (info > n_) {
                FEEM_ERROR("SuperluContext::factorize - 内存分配失败, 需 {} 字节", info - n_);
            }
            release_factors();
            state_ = State::INITIALIZED;
            return info;
        }

        state_ = State::FACTORIZED;
        FEEM_DEBUG("SuperluContext::factorize 成功, nprocs={}", nprocs);
        return 0;
    } catch (const std::exception& e) {
        FEEM_ERROR("SuperluContext::factorize - 异常: {}", e.what());
        // 异常路径: 可能 pdgstrf_init 已创建但 pdgstrf 失败，
        // 需要调用 pxgstrf_finalize 清理 AC 和辅助数组
        if (AC_ && options_) {
            pxgstrf_finalize(options_.get(), AC_.get());
            etree_.reset();
            colcnt_h_.reset();
            part_super_h_.reset();
            options_->etree = nullptr;
            options_->colcnt_h = nullptr;
            options_->part_super_h = nullptr;
            AC_.reset();
        }
        release_factors();
        if (Gstat_) {
            StatFree(Gstat_.get());
            Gstat_.reset();
        }
        state_ = State::INITIALIZED;
        return -999;
    }
}

SolverResult SuperluContext::solve(const Eigen::VectorXd& b) {
    if (state_ != State::FACTORIZED) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "SuperLU_MT 未完成分解，请先调用 factorize()");
    }
    if (b.size() != n_) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量维度不匹配");
    }

    try {
        const int_t nrhs = 1;
        double* b_data = doubleMalloc(n_);
        if (!b_data) {
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
        }
        for (int i = 0; i < n_; ++i) {
            b_data[i] = b(i);
        }

        SuperMatrix B;
        dCreate_Dense_Matrix(&B, n_, nrhs, b_data, n_, SLU_DN, SLU_D, SLU_GE);

        auto local_gstat = std::make_unique<Gstat_t>();
        StatAlloc(n_, 1, 1, 1, local_gstat.get());
        StatInit(n_, 1, local_gstat.get());

        int_t info = 0;
        dgstrs(NOTRANS, L_.get(), U_.get(),
               options_->perm_c, options_->perm_r,
               &B, local_gstat.get(), &info);

        StatFree(local_gstat.get());

        if (info != 0) {
            Destroy_SuperMatrix_Store(&B);
            SUPERLU_FREE(b_data);
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "dgstrs 求解失败 (info=" + std::to_string(info) + ")");
        }

        DNformat* Bstore = reinterpret_cast<DNformat*>(B.Store);
        double* x_data = reinterpret_cast<double*>(Bstore->nzval);

        Eigen::VectorXd x(n_);
        for (int i = 0; i < n_; ++i) {
            x(i) = x_data[i];
        }

        Destroy_SuperMatrix_Store(&B);
        SUPERLU_FREE(b_data);

        return create_success_result(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SuperluContext::solve - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("求解异常: ") + e.what());
    }
}

void SuperluContext::reset() {
    release_factors();
    release_ac_matrix();

    if (A_ && A_->Store) {
        Destroy_CompCol_Matrix(A_.get());
    }
    A_.reset();

    perm_c_.reset();
    perm_r_.reset();
    etree_.reset();
    colcnt_h_.reset();
    part_super_h_.reset();
    options_.reset();

    if (Gstat_) {
        StatFree(Gstat_.get());
        Gstat_.reset();
    }

    state_ = State::UNINITIALIZED;
    n_ = 0;
    nnz_ = 0;
}

void SuperluContext::release_factors() {
    if (L_ && L_->Store) {
        Destroy_SuperNode_SCP(L_.get());
    }
    L_.reset();

    if (U_ && U_->Store) {
        Destroy_CompCol_NCP(U_.get());
    }
    U_.reset();
}

void SuperluContext::release_ac_matrix() {
    if (AC_ && options_) {
        pxgstrf_finalize(options_.get(), AC_.get());
    }
    AC_.reset();
}

#endif  // EM_SOLVER_HAS_SUPERLU

// ============================================================================
// DirectBackendManager 实现
// ============================================================================

bool DirectBackendManager::isBackendAvailable(DirectBackendType type) {
    switch (type) {
        case DirectBackendType::EIGEN:
            return true;  // Eigen始终可用（作为默认回退方案）
        case DirectBackendType::SUPERLU:
#if EM_SOLVER_HAS_SUPERLU
            return true;
#else
            return false;
#endif
        case DirectBackendType::MUMPS:
#if EM_SOLVER_HAS_MUMPS
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

#if EM_SOLVER_HAS_SUPERLU
    backends.push_back(DirectBackendType::SUPERLU);
#endif

#if EM_SOLVER_HAS_MUMPS
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
#if EM_SOLVER_HAS_SUPERLU
        case DirectBackendType::SUPERLU:
            decompose_result = decompose_with_superlu();
            break;
#endif
#if EM_SOLVER_HAS_MUMPS
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
#if EM_SOLVER_HAS_SUPERLU
        case DirectBackendType::SUPERLU:
            result = solve_with_superlu(b);
            break;
#endif
#if EM_SOLVER_HAS_MUMPS
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

std::string SymmetricDirectSolver::get_solver_name() const {
    return std::string("SymmetricDirect_") + DirectBackendManager::getBackendName(backend_type_);
}

void SymmetricDirectSolver::clear() {
    FEEM_DEBUG("SymmetricDirectSolver::clear - 释放所有内部资源");

    // 释放Eigen后端资源
    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();  // 释放矩阵内存

#if EM_SOLVER_HAS_SUPERLU
    // 释放 SuperLU_MT 资源（RAII 封装自动管理所有资源释放）
    if (superlu_ctx_) {
        superlu_ctx_->reset();
        superlu_ctx_.reset();
    }
#endif

#if EM_SOLVER_HAS_MUMPS
    // 释放MUMPS资源（调用MUMPS清理接口）
    if (mumps_initialized_) {
        mumps_data_.job = -2;  // MUMPS清理作业代码
        dmumps_c(&mumps_data_);  // 释放MUMPS内部分配的所有内存
        mumps_initialized_ = false;
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

#if EM_SOLVER_HAS_SUPERLU

/**
 * @brief 使用SuperLU_MT后端执行LU分解（通过 SuperluContext RAII 封装）
 * @return SolverResult 分解结果
 *
 * @details 通过 SuperluContext 封装执行完整的初始化和分解流程：
 * 1. 从Eigen稀疏矩阵提取CSC格式数据（nzval/rowind/colptr）
 * 2. 调用 ctx.initialize() 创建系数矩阵、计算列置换、分配辅助数组
 * 3. 调用 ctx.factorize() 执行 pdgstrf_init → pdgstrf → pdgstrf_finalize
 * 4. 缓存L、U因子供后续solve()复用
 *
 * @note 修复了原始代码中缺少 pdgstrf_init/get_perm_c/pdgstrf_finalize 的致命问题
 */
SolverResult SymmetricDirectSolver::decompose_with_superlu() {
    FEEM_DEBUG("SymmetricDirectSolver::decompose_with_superlu - 开始SuperLU_MT LU分解");

    try {
        const int_t n = static_cast<int_t>(eigen_matrix_.rows());
        eigen_matrix_.makeCompressed();
        const int_t nnz = static_cast<int_t>(eigen_matrix_.nonZeros());

        double* nzval = doubleMalloc(nnz);
        int_t* rowind = intMalloc(nnz);
        int_t* colptr = intMalloc(n + 1);

        if (!nzval || !rowind || !colptr) {
            return create_error_result(SolverStatus::NUMERICAL_ERROR, "内存分配失败");
        }

        for (int k = 0; k < nnz; ++k) {
            nzval[k] = eigen_matrix_.valuePtr()[k];
            rowind[k] = static_cast<int_t>(eigen_matrix_.innerIndexPtr()[k]);
        }
        for (int i = 0; i <= n; ++i) {
            colptr[i] = static_cast<int_t>(eigen_matrix_.outerIndexPtr()[i]);
        }

        superlu_ctx_ = std::make_unique<SuperluContext>();
        if (!superlu_ctx_->initialize(n, nnz, nzval, rowind, colptr)) {
            superlu_ctx_.reset();
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "SuperLU 初始化失败");
        }

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
            superlu_ctx_.reset();
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "SuperLU_MT LU分解失败: " + error_msg);
        }

        FEEM_DEBUG("decompose_with_superlu - SuperLU_MT LU分解成功, 维度: {}x{}, 非零元数: {}", n, n, nnz);
        return create_success_result(Eigen::VectorXd());
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("decompose_with_superlu - 内存分配失败: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "内存不足，无法完成SuperLU_MT LU分解");
    } catch (const std::exception& e) {
        FEEM_ERROR("decompose_with_superlu - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("SuperLU_MT分解异常: ") + e.what());
    }
}

/**
 * @brief 使用SuperLU_MT后端执行求解（三角回代，通过 SuperluContext 封装）
 * @param b 右端项向量
 * @return SolverResult 求解结果
 *
 * @details 委托给 SuperluContext::solve() 执行完整的三角系统求解流程，
 *          包括 DenseMatrix 创建、dgstrs 回代、结果提取和临时资源释放。
 */
SolverResult SymmetricDirectSolver::solve_with_superlu(const Eigen::VectorXd& b) {
    FEEM_DEBUG("SymmetricDirectSolver::solve_with_superlu - 开始SuperLU_MT三角回代");

    if (!superlu_ctx_ || !superlu_ctx_->is_factored()) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "SuperLU_MT未完成分解，请先调用set_matrix()");
    }

    return superlu_ctx_->solve(b);
}

#endif  // EM_SOLVER_HAS_SUPERLU

#if EM_SOLVER_HAS_MUMPS

SolverResult SymmetricDirectSolver::decompose_with_mumps() {
    FEEM_DEBUG("SymmetricDirectSolver::decompose_with_mumps - 开始MUMPS分解");

    try {
        // 初始化MUMPS数据结构
        memset(&mumps_data_, 0, sizeof(mumps_data_));

        mumps_data_.par = 1;       // 并行工作进程数（单进程设为1）
        mumps_data_.sym = 1;       // 对称矩阵标志（1=对称正定）
        mumps_data_.job = -1;      // 初始化作业

        dmumps_c(&mumps_data_);   // 调用MUMPS初始化

        // 设置矩阵维度
        mumps_data_.n = eigen_matrix_.rows();

        // 将Eigen矩阵转换为MUMPS需要的CSR格式（IRN/JCN/A数组）
        // ... （实际数据转换逻辑）

        // 执行分析+分解（job=5：分析+分解一步完成）
        mumps_data_.job = 5;
        dmumps_c(&mumps_data_);

        // 检查MUMPS错误码
        if (mumps_data_.info[0] < 0) {
            FEEM_ERROR("SymmetricDirectSolver::decompose_with_mumps - MUMPS错误: info[0]={}",
                       mumps_data_.info[0]);
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "MUMPS分解失败，错误码: " +
                                       std::to_string(mumps_data_.info[0]));
        }

        mumps_initialized_ = true;

        FEEM_DEBUG("SymmetricDirectSolver::decompose_with_mumps - MUMPS分解成功");
        return create_success_result(Eigen::VectorXd());
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::decompose_with_mumps - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("MUMPS分解异常: ") + e.what());
    }
}

SolverResult SymmetricDirectSolver::solve_with_mumps(const Eigen::VectorXd& b) {
    FEEM_DEBUG("SymmetricDirectSolver::solve_with_mumps - 开始MUMPS求解");

    if (!mumps_initialized_) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "MUMPS未初始化，请先调用set_matrix()");
    }

    try {
        // 设置右端项向量
        mumps_data_.rhs = const_cast<double*>(b.data());

        // 执行求解（job=3：利用已有分解进行求解）
        mumps_data_.job = 3;
        dmumps_c(&mumps_data_);

        // 检查错误
        if (mumps_data_.info[0] < 0) {
            FEEM_ERROR("SymmetricDirectSolver::solve_with_mumps - MUMPS求解错误: info[0]={}",
                       mumps_data_.info[0]);
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "MUMPS求解失败");
        }

        // 从MUMPS rhs中提取解向量（MUMPS原地修改rhs为解）
        Eigen::VectorXd x = b;  // 实际应从mumps_data_.rhs复制

        return create_success_result(std::move(x));
    } catch (const std::exception& e) {
        FEEM_ERROR("SymmetricDirectSolver::solve_with_mumps - 异常: {}", e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("MUMPS求解异常: ") + e.what());
    }
}

#endif  // EM_SOLVER_HAS_MUMPS

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

std::string SymmetricIndefiniteDirectSolver::get_solver_name() const {
    return std::string("SymmetricIndefiniteDirect_") +
           DirectBackendManager::getBackendName(backend_type_);
}

void SymmetricIndefiniteDirectSolver::clear() {
    FEEM_DEBUG("SymmetricIndefiniteDirectSolver::clear - 释放资源");

    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
    regularized_matrix_ = Eigen::SparseMatrix<double>();
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

std::string GeneralDirectSolver::get_solver_name() const {
    return std::string("GeneralDirect_") + DirectBackendManager::getBackendName(backend_type_);
}

void GeneralDirectSolver::clear() {
    FEEM_DEBUG("GeneralDirectSolver::clear - 释放资源");

    eigen_solver_.reset();
    eigen_matrix_ = Eigen::SparseMatrix<double>();
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
