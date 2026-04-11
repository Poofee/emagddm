/**
 * @file direct_solvers.h
 * @brief 求解器后端辅助类定义（SuperLU/MUMPS RAII 封装）
 * @details 提供求解器后端所需的辅助数据结构封装：
 *          - SuperluContext：SuperLU_MT 数据结构 RAII 封装
 *          - MumpsContext：MUMPS 实数求解器 RAII 封装
 *          - ComplexMumpsContext：MUMPS 复数求解器 RAII 封装
 *
 * @note 原三个直接求解器类（SymmetricDirectSolver/SymmetricIndefiniteDirectSolver/GeneralDirectSolver）
 *       已迁移至新架构（UnifiedDirectSolver + SolverBackend 策略模式），本文件仅保留辅助类。
 *
 * @see SuperLUBackend SuperLU 后端实现
 * @see MUMPSBackend MUMPS 后端实现
 * @see UnifiedDirectSolver 统一调度层
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 2.0 (重构版 - 移除旧求解器类)
 */

#pragma once

#include "em_linear_solver.h"
#include "sparse_converter.h"
#include "em_solver_backends.hpp"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <memory>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>

#ifdef HAVE_SUPERLU
// SuperLU头文件已在em_solver_backends.hpp中包含
#endif

#ifdef HAVE_MUMPS
// MUMPS头文件已在em_solver_backends.hpp中包含
#endif

namespace numeric {

#ifdef HAVE_SUPERLU

/**
 * @class SuperluContext
 * @brief SuperLU_MT 数据结构 RAII 封装，管理完整的 LU 分解生命周期
 * @details 封装 SuperLU_MT 所有的 C 语言数据结构和资源，提供类型安全的 C++ 接口。
 */
class SuperluContext {
public:
    SuperluContext();
    ~SuperluContext();

    SuperluContext(const SuperluContext&) = delete;
    SuperluContext& operator=(const SuperluContext&) = delete;
    SuperluContext(SuperluContext&& other) noexcept;
    SuperluContext& operator=(SuperluContext&& other) noexcept;

    bool initialize(int_t n, int_t nnz, double* nzval, int_t* rowind, int_t* colptr);
    int_t factorize(int_t nprocs = 1);
    SolverResult solve(const Eigen::VectorXd& b);
    void reset();
    bool is_factored() const;
    bool is_initialized() const;
    int_t matrix_size() const;

private:
    enum class State { UNINITIALIZED, INITIALIZED, FACTORIZED };
    State state_ = State::UNINITIALIZED;
    int_t n_ = 0;
    int_t nnz_ = 0;

    std::unique_ptr<SuperMatrix> A_;
    std::unique_ptr<SuperMatrix> L_;
    std::unique_ptr<SuperMatrix> U_;
    std::unique_ptr<SuperMatrix> B_;  // 右端项矩阵
    std::unique_ptr<SuperMatrix> X_;  // 解向量矩阵

    std::unique_ptr<int_t[]> perm_c_;
    std::unique_ptr<int_t[]> perm_r_;

    std::unique_ptr<superlumt_options_t> options_;

    void release_factors();
};

#endif  // HAVE_SUPERLU

#ifdef HAVE_MUMPS

/**
 * @class MumpsContext
 * @brief MUMPS 双精度直接求解器 RAII 封装
 */
class MumpsContext {
public:
    MumpsContext();
    ~MumpsContext();

    MumpsContext(const MumpsContext&) = delete;
    MumpsContext& operator=(const MumpsContext&) = delete;
    MumpsContext(MumpsContext&& other) noexcept;
    MumpsContext& operator=(MumpsContext&& other) noexcept;

    bool initialize(int sym = 0);
    int factorize(int n, int nz, int* irn, int* jcn, double* a);
    SolverResult solve(const Eigen::VectorXd& b);
    void reset();
    bool is_factored() const;
    bool is_initialized() const;
    int matrix_size() const;

private:
    enum class State { UNINITIALIZED, INITIALIZED, FACTORIZED };
    State state_ = State::UNINITIALIZED;
    int n_ = 0;
    DMUMPS_STRUC_C mumps_data_;
    std::vector<double> rhs_buffer_;
    std::vector<int> irn_storage_;
    std::vector<int> jcn_storage_;
    std::vector<double> a_storage_;

    void configure_default_icntl();

public:
    int factorize_from_csr(const CsrMatrix<double>& csr, int sym = 0);
};

/**
 * @class ComplexMumpsContext
 * @brief MUMPS 双精度复数直接求解器 RAII 封装
 */
class ComplexMumpsContext {
public:
    ComplexMumpsContext();
    ~ComplexMumpsContext();

    ComplexMumpsContext(const ComplexMumpsContext&) = delete;
    ComplexMumpsContext& operator=(const ComplexMumpsContext&) = delete;
    ComplexMumpsContext(ComplexMumpsContext&& other) noexcept;
    ComplexMumpsContext& operator=(ComplexMumpsContext&& other) noexcept;

    bool initialize(int sym = 0);
    int factorize(int n, int nz, int* irn, int* jcn, std::complex<double>* a);
    SolverResult solve(const Eigen::VectorXcd& b);
    void reset();
    bool is_factored() const;
    bool is_initialized() const;
    int matrix_size() const;

private:
    enum class State { UNINITIALIZED, INITIALIZED, FACTORIZED };
    State state_ = State::UNINITIALIZED;
    int n_ = 0;
    ZMUMPS_STRUC_C zmumps_data_;
    std::vector<ZMUMPS_COMPLEX> rhs_buffer_;
    std::vector<int> irn_storage_;
    std::vector<int> jcn_storage_;
    std::vector<std::complex<double>> a_storage_;

    void configure_default_icntl();

public:
    int factorize_from_csr(const CsrMatrix<std::complex<double>>& csr, int sym = 0);
};

/**
 * @brief 从 SolverResult 的实数向量中提取复数解
 */
inline Eigen::VectorXcd extract_complex(const SolverResult& result) {
    const int n = result.x.size() / 2;
    Eigen::VectorXcd x(n);
    for (int i = 0; i < n; ++i) {
        x(i) = std::complex<double>(result.x(2 * i), result.x(2 * i + 1));
    }
    return x;
}

#endif  // HAVE_MUMPS


} // namespace numeric
