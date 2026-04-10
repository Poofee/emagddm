/**
 * @file em_direct_solvers_superlu.cpp
 * @brief SuperLU_MT 后端直接求解器实现（独立编译单元）
 * @details 将 SuperLU_MT 相关代码从 em_direct_solvers.cpp 拆分出来，
 *          与 MUMPS 拆分模式保持一致，便于独立管理第三方库依赖。
 *
 * 包含：
 * - SuperluContext RAII 封装类的完整实现
 * - SymmetricDirectSolver 的 SuperLU 后端方法（decompose_with_superlu / solve_with_superlu）
 *
 * @author Poofee
 * @date 2026-04-09
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"

#ifdef HAVE_SUPERLU

namespace {

numeric::SolverResult create_success_result(Eigen::VectorXd x, int iterations = 0) {
    numeric::SolverResult result;
    result.x = std::move(x);
    result.status = numeric::SolverStatus::SUCCESS;
    result.iterations = iterations;
    return result;
}

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

/**
 * @brief 初始化 SuperLU_MT 求解器
 * @param n 矩阵维度
 * @param nnz 非零元个数
 * @param nzval 非零值数组（CSC 格式）
 * @param rowind 行索引数组（CSC 格式）
 * @param colptr 列指针数组（CSC 格式）
 * @return true 初始化成功，false 初始化失败
 */
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

/**
 * @brief 执行 LU 分解
 * @param nprocs 线程数
 * @return 0 成功，非零 SuperLU 错误码
 */
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
        AC_ = std::make_unique<SuperMatrix>();
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

        etree_.reset();
        colcnt_h_.reset();
        part_super_h_.reset();
        options_->etree = nullptr;
        options_->colcnt_h = nullptr;
        options_->part_super_h = nullptr;

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

/**
 * @brief 使用 SuperLU_MT 执行三角回代求解
 * @param b 右端项向量
 * @return SolverResult 求解结果
 */
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

/**
 * @brief 重置 SuperLU_MT 实例，释放所有资源
 */
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

/**
 * @brief 释放 L、U 因子
 */
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

/**
 * @brief 释放 AC 矩阵（列置换后的矩阵）
 */
void SuperluContext::release_ac_matrix() {
    if (AC_ && options_) {
        pxgstrf_finalize(options_.get(), AC_.get());
    }
    AC_.reset();
}

// ============================================================================
// SymmetricDirectSolver SuperLU 后端
// ============================================================================

/**
 * @brief 使用SuperLU_MT后端执行LU分解（通过 SuperluContext RAII 封装）
 * @return SolverResult 分解结果
 *
 * @details 通过 SuperluContext 封装执行完整的初始化和分解流程：
 * 1. 从Eigen稀疏矩阵提取CSC格式数据（nzval/rowind/colptr）
 * 2. 调用 ctx.initialize() 创建系数矩阵、计算列置换、分配辅助数组
 * 3. 调用 ctx.factorize() 执行 pdgstrf_init → pdgstrf → pdgstrf_finalize
 * 4. 缓存L、U因子供后续solve()复用
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
 */
SolverResult SymmetricDirectSolver::solve_with_superlu(const Eigen::VectorXd& b) {
    FEEM_DEBUG("SymmetricDirectSolver::solve_with_superlu - 开始SuperLU_MT三角回代");

    if (!superlu_ctx_ || !superlu_ctx_->is_factored()) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "SuperLU_MT未完成分解，请先调用set_matrix()");
    }

    return superlu_ctx_->solve(b);
}

} // namespace numeric

#endif // HAVE_SUPERLU
