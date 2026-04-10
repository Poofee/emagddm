/**
 * @file em_direct_solvers_superlu.cpp
 * @brief SuperluContext RAII 封装类完整实现
 * @details 实现 SuperluContext 类的全部方法，封装 SuperLU_MT 库的并行 LU 分解功能。
 *
 * @note 原SymmetricDirectSolver的SuperLU后端方法已迁移至新架构：
 *       - SuperLUBackend（统一后端，策略模式）
 *
 * @see SuperLUBackend SuperLU 后端实现
 * @see UnifiedDirectSolver 统一调度层
 *
 * @author Poofee
 * @date 2026-04-09
 * @version 2.0 (重构版 - 移除旧求解器类引用)
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"
#include <stdexcept>
#include <cstring>

#ifdef HAVE_SUPERLU

namespace numeric {

// ============================================================================
// SuperluContext 构造/析构
// ============================================================================

SuperluContext::SuperluContext()
    : state_(State::UNINITIALIZED), n_(0), nnz_(0) {
}

SuperluContext::~SuperluContext() {
    reset();
}

SuperluContext::SuperluContext(SuperluContext&& other) noexcept
    : state_(other.state_), n_(other.n_), nnz_(other.nnz_),
      A_(std::move(other.A_)), L_(std::move(other.L_)),
      U_(std::move(other.U_)), AC_(std::move(other.AC_)),
      perm_c_(std::move(other.perm_c_)), perm_r_(std::move(other.perm_r_)),
      etree_(std::move(other.etree_)), colcnt_h_(std::move(other.colcnt_h_)),
      part_super_h_(std::move(other.part_super_h_)),
      options_(std::move(other.options_)), Gstat_(std::move(other.Gstat_)) {
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

// ============================================================================
// 初始化与分解
// ============================================================================

bool SuperluContext::initialize(int_t n, int_t nnz, double* nzval, int_t* rowind, int_t* colptr) {
    if (state_ != State::UNINITIALIZED) {
        FEEM_WARN("SuperluContext::initialize - 上下文已初始化，先重置");
        reset();
    }

    try {
        A_ = std::make_unique<SuperMatrix>();
        L_ = std::make_unique<SuperMatrix>();
        U_ = std::make_unique<SuperMatrix>();
        AC_ = std::make_unique<SuperMatrix>();

        options_ = std::make_unique<superlumt_options_t>();
        Gstat_ = std::make_unique<Gstat_t>();

        dCreate_Dense_Matrix(A_.get(), n, n, SLU_NC, SLU_D, SLU_GE);
        ((NCformat*)A_->Store)->nzval = nzval;
        ((NCformat*)A_->Store)->rowind = rowind;
        ((NCformat*)A_->Store)->colptr = colptr;

        set_default_options(options_.get());
        options_->ColPerm = MATA_AT_PLUS;
        options_->RowPerm = LargeDiag_MC64;
        options_->ILU_Flag = 0;

        StatInit(Gstat_.get());
        perm_c_ = std::make_unique<int_t[]>(n);
        perm_r_ = std::make_unique<int_t[]>(n);

        n_ = n;
        nnz_ = nnz;
        state_ = State::INITIALIZED;
        return true;

    } catch (const std::exception& e) {
        FEEM_ERROR("SuperluContext::initialize - 异常: {}", e.what());
        reset();
        return false;
    }
}

int_t SuperluContext::factorize(int_t nprocs) {
    if (state_ != State::INITIALIZED) {
        FEEM_ERROR("SuperluContext::factorize - 未初始化或已分解");
        return -1;
    }

    try {
        etree_ = std::make_unique<int_t[]>(2 * n_);
        colcnt_h_ = std::make_unique<int_t[]>(n_);
        part_super_h_ = std::make_unique<int_t[]>(n_);

        int_t info = 0;
        pdgstrf_init(n_, options_.get(), A_.get(), perm_c_.get(), perm_r_.get(),
                     etree_.get(), colcnt_h_.get(), part_super_h_.get(),
                     L_.get(), U_.get(), Gstat_.get());

        if (info != 0) {
            FEEM_ERROR("SuperluContext::factorize - pdgstrf_init失败: info={}", info);
            release_factors();
            return info;
        }

        pdgstrf(options_.get(), A_.get(), perm_c_.get(), perm_r_.get(),
                 etree_.get(), colcnt_h_.get(), part_super_h_.get(),
                 L_.get(), U_.get(), &info);

        if (info != 0) {
            FEEM_ERROR("SuperluContext::factorize - pdgstrf失败: info={}", info);
            release_factors();
            return info;
        }

        pxgstrf_finalize(options_.get(), L_.get());
        pxgstrf_finalize(options_.get(), U_.get());

        if (AC_) {
            AC_ = nullptr;
        }

        state_ = State::FACTORIZED;
        return 0;

    } catch (const std::exception& e) {
        FEEM_ERROR("SuperluContext::factorize - 异常: {}", e.what());
        release_factors();
        return -999;
    }
}

// ============================================================================
// 求解
// ============================================================================

SolverResult SuperluContext::solve(const Eigen::VectorXd& b) {
    if (!is_factored()) {
        return create_error_result(SolverStatus::NUMERICAL_ERROR, "未分解");
    }
    if (b.size() != n_) {
        return create_error_result(SolverStatus::INVALID_INPUT, "维度不匹配");
    }

    try {
        Eigen::VectorXd x(b.size());
        SuperMatrix B_mat, X_mat;
        dCreate_Dense_Matrix(&B_mat, n_, 1, SLU_DN, SLU_D, SLU_GE);
        dCreate_Dense_Matrix(&X_mat, n_, 1, SLU_DN, SLU_D, SLU_GE);
        DNformat* B_store = static_cast<DNformat*>(B_mat.Store);
        DNformat* X_store = static_cast<DNformat*>(X_mat.Store);
        memcpy(B_store->nzval, b.data(), sizeof(double) * n_);

        trans_t trans = NOTRANS;
        dgstrs(trans, L_.get(), U_.get(), perm_c_.get(), perm_r_.get(),
               &B_mat, &X_mat, Gstat_.get(), &info);

        if (info != 0) {
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                     "dgstrs失败: info=" + std::to_string(info));
        }

        memcpy(x.data(), X_store->nzval, sizeof(double) * n_);
        Destroy_SuperMatrix_Store(&B_mat);
        Destroy_SuperMatrix_Store(&X_mat);

        return create_success_result(std::move(x));

    } catch (const std::exception& e) {
        return create_error_result(SolverStatus::NUMERICAL_ERROR, e.what());
    }
}

// ============================================================================
// 资源管理
// ============================================================================

void SuperluContext::reset() {
    release_factors();
    release_ac_matrix();

    Destroy_SuperMatrix_Store(A_.get());
    A_.reset();
    L_.reset();
    U_.reset();

    perm_c_.reset();
    perm_r_.reset();
    etree_.reset();
    colcnt_h_.reset();
    part_super_h_.reset();
    options_.reset();
    Gstat_.reset();

    StatFree(Gstat_.get());

    n_ = 0;
    nnz_ = 0;
    state_ = State::UNINITIALIZED;
}

void SuperluContext::release_factors() {
    if (L_ && L_->nrow >= 0) {
        Destroy_SuperNode_Matrix(L_.get());
    }
    if (U_ && U_->nrow >= 0) {
        Destroy_CompCol_Matrix(U_.get());
    }
    L_ = nullptr;
    U_ = nullptr;
}

void SuperluContext::release_ac_matrix() {
    if (AC_ && AC_->nrow >= 0) {
        pxgstrf_finalize(options_.get(), AC_.get());
    }
    AC_ = nullptr;
}

bool SuperluContext::is_factored() const { return state_ == State::FACTORIZED; }
bool SuperluContext::is_initialized() const { return state_ != State::UNINITIALIZED; }
int SuperluContext::matrix_size() const { return n_; }

} // namespace numeric

#endif // HAVE_SUPERLU
