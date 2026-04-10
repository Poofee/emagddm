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
#include <cstdlib>

#ifdef HAVE_SUPERLU

namespace numeric {

// ============================================================================
// 辅助函数（与 MUMPS 实现保持一致）
// ============================================================================

static SolverResult create_success_result(Eigen::VectorXd x, int iterations = 0) {
    SolverResult result;
    result.status = SolverStatus::SUCCESS;
    result.x = std::move(x);
    result.iterations = iterations;
    return result;
}

static SolverResult create_error_result(SolverStatus status, const std::string& error_msg) {
    SolverResult result;
    result.status = status;
    result.error_msg = error_msg;
    return result;
}

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
      U_(std::move(other.U_)), B_(std::move(other.B_)),
      X_(std::move(other.X_)), perm_c_(std::move(other.perm_c_)),
      perm_r_(std::move(other.perm_r_)),
      options_(std::move(other.options_)) {
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
        B_ = std::move(other.B_);
        X_ = std::move(other.X_);
        perm_c_ = std::move(other.perm_c_);
        perm_r_ = std::move(other.perm_r_);
        options_ = std::move(other.options_);

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
        // 创建稀疏矩阵 A（列压缩格式）
        A_ = std::make_unique<SuperMatrix>();
        dCreate_CompCol_Matrix(A_.get(), n, n, nnz, nzval, rowind, colptr,
                               SLU_NC, SLU_D, SLU_GE);

        // 创建解和右端项矩阵
        B_ = std::make_unique<SuperMatrix>();
        X_ = std::make_unique<SuperMatrix>();

        // 配置 SuperLU_MT 选项
        options_ = std::make_unique<superlumt_options_t>();
        options_->nprocs = 1;  // 默认单线程，可通过 factorize 修改
        options_->fact = DOFACT;
        options_->trans = NOTRANS;
        options_->refact = NO;
        options_->panel_size = sp_ienv(1);
        options_->relax = sp_ienv(2);
        options_->usepr = NO;
        options_->drop_tol = 0.0;
        options_->diag_pivot_thresh = 1.0;
        options_->SymmetricMode = NO;
        options_->PrintStat = NO;

        // 分配排列向量和工作数组
        perm_c_ = std::unique_ptr<int_t[]>(static_cast<int_t*>(intMalloc(n)));
        perm_r_ = std::unique_ptr<int_t[]>(static_cast<int_t*>(intMalloc(n)));

        if (!perm_c_ || !perm_r_) {
            throw std::runtime_error("内存分配失败: perm_c/perm_r");
        }

        // 获取列排列（最小度排序）
        int_t permc_spec = 1;  // MMD on A'*A
        get_perm_c(permc_spec, A_.get(), perm_c_.get());

        // 分配选项中的工作数组
        options_->etree = static_cast<int_t*>(intMalloc(n));
        options_->colcnt_h = static_cast<int_t*>(intMalloc(n));
        options_->part_super_h = static_cast<int_t*>(intMalloc(n));

        if (!options_->etree || !options_->colcnt_h || !options_->part_super_h) {
            throw std::runtime_error("内存分配失败: 工作数组");
        }

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

    (void)nprocs;  // 通过 options_->nprocs 设置

    try {
        // 创建 L 和 U 矩阵
        L_ = std::make_unique<SuperMatrix>();
        U_ = std::make_unique<SuperMatrix>();

        // 分配均衡和误差界数组
        equed_t equed = NOEQUIL;
        auto R = std::unique_ptr<double[]>(static_cast<double*>(SUPERLU_MALLOC(n_ * sizeof(double))));
        auto C = std::unique_ptr<double[]>(static_cast<double*>(SUPERLU_MALLOC(n_ * sizeof(double))));
        auto ferr = std::unique_ptr<double[]>(static_cast<double*>(SUPERLU_MALLOC(sizeof(double))));
        auto berr = std::unique_ptr<double[]>(static_cast<double*>(SUPERLU_MALLOC(sizeof(double))));

        double rpg = 0.0, rcond = 0.0;
        superlu_memusage_t mem_usage{};
        int_t info = 0;

        // 执行 LU 分解（使用 pdgssvx 高级接口）
        pdgssvx(options_->nprocs, options_.get(), A_.get(),
                perm_c_.get(), perm_r_.get(),
                &equed, R.get(), C.get(),
                L_.get(), U_.get(),
                nullptr, nullptr,  // B 和 X 在 solve 中设置
                &rpg, &rcond,
                ferr.get(), berr.get(),
                &mem_usage, &info);

        if (info < 0 && info != -(n_ + 1)) {
            FEEM_ERROR("SuperluContext::factorize - pdgssvx失败: info={}", info);
            release_factors();
            return info;
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
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n_);

        // 创建右端项和解向量的密集矩阵
        auto rhsb = std::unique_ptr<double[]>(static_cast<double*>(doubleMalloc(n_)));
        auto rhsx = std::unique_ptr<double[]>(static_cast<double*>(doubleMalloc(n_)));

        if (!rhsb || !rhsx) {
            throw std::runtime_error("内存分配失败: rhs");
        }

        std::memcpy(rhsb.get(), b.data(), sizeof(double) * n_);
        std::memset(rhsx.get(), 0, sizeof(double) * n_);

        dCreate_Dense_Matrix(B_.get(), n_, 1, rhsb.get(), n_, SLU_DN, SLU_D, SLU_GE);
        dCreate_Dense_Matrix(X_.get(), n_, 1, rhsx.get(), n_, SLU_DN, SLU_D, SLU_GE);

        // 使用 dgstrs 求解（前代-回代）
        // 注意：dgstrs 直接在 B 矩阵中存储解
        trans_t trans = NOTRANS;
        Gstat_t Gstat{};
        StatInit(n_, 1, &Gstat);
        int_t info = 0;
        dgstrs(trans, L_.get(), U_.get(), perm_r_.get(), perm_c_.get(),
               B_.get(), &Gstat, &info);
        StatFree(&Gstat);

        if (info != 0) {
            Destroy_SuperMatrix_Store(B_.get());
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                     "dgstrs失败: info=" + std::to_string(info));
        }

        // 从 B 矩阵复制解向量
        DNformat* B_store = static_cast<DNformat*>(B_->Store);
        std::memcpy(x.data(), B_store->nzval, sizeof(double) * n_);

        // 清理密集矩阵
        Destroy_SuperMatrix_Store(B_.get());

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

    if (A_) {
        Destroy_CompCol_Matrix(A_.get());
        A_.reset();
    }
    if (B_) {
        Destroy_SuperMatrix_Store(B_.get());
        B_.reset();
    }
    if (X_) {
        Destroy_SuperMatrix_Store(X_.get());
        X_.reset();
    }

    perm_c_.reset();
    perm_r_.reset();

    if (options_) {
        SUPERLU_FREE(options_->etree);
        SUPERLU_FREE(options_->colcnt_h);
        SUPERLU_FREE(options_->part_super_h);
        options_.reset();
    }

    n_ = 0;
    nnz_ = 0;
    state_ = State::UNINITIALIZED;
}

void SuperluContext::release_factors() {
    if (L_ && L_->nrow >= 0) {
        Destroy_SuperNode_SCP(L_.get());
    }
    if (U_ && U_->nrow >= 0) {
        Destroy_CompCol_NCP(U_.get());
    }
    L_ = nullptr;
    U_ = nullptr;
}

bool SuperluContext::is_factored() const { return state_ == State::FACTORIZED; }
bool SuperluContext::is_initialized() const { return state_ != State::UNINITIALIZED; }
int SuperluContext::matrix_size() const { return n_; }

} // namespace numeric

#endif // HAVE_SUPERLU
