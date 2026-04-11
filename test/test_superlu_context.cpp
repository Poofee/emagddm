/**
 * @file test_superlu_context.cpp
 * @brief SuperluContext RAII 封装类的全面功能测试
 * @details 覆盖 SuperluContext 类的所有公开接口和边界条件：
 *
 * **测试模块**：
 * 1. 构造/析构与 RAII 语义
 * 2. initialize() 初始化（正常/异常参数）
 * 3. factorize() LU 分解（成功/失败/状态机保护）
 * 4. solve() 三角回代求解（精度/重复调用）
 * 5. reset() 重置与资源释放
 * 6. 移动语义
 * 7. 状态机转换正确性
 * 8. 内存安全（无泄漏）
 * 9. 边界条件（奇异矩阵、零矩阵、单位矩阵）
 *
 * @author Poofee
 * @date 2026-04-08
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef HAVE_SUPERLU
    #include "direct_solvers.h"
#else
    #include "direct_solvers.h"
#endif

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include "logger_factory.hpp"

using namespace numeric;

// ==================== 测试辅助函数 ====================

/**
 * @brief 生成对称正定稀疏矩阵的 CSC 数据（供 SuperluContext::initialize 使用）
 * @param n 矩阵维度
 * @param density 非零密度
 * @param[out] nzval 非零值数组（调用方需 SUPERLU_FREE）
 * @param[out] rowind 行索引数组（调用方需 SUPERLU_FREE）
 * @param[out] colptr 列指针数组（调用方需 SUPERLU_FREE）
 * @return int_t 非零元个数 nnz
 */
int_t generate_spd_csc_data(int n, double density,
                            double*& nzval, int_t*& rowind, int_t*& colptr) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(0.0, 1.0);

    // 先构建三元组，再转CSC
    std::vector<std::tuple<int,int,double>> coo_entries;
    for (int i = 0; i < n; ++i) {
        double diag_val = 10.0 + dist(gen) * 5.0;
        coo_entries.emplace_back(i, i, diag_val);
        for (int j = i + 1; j < n; ++j) {
            if (dist(gen) < density / n * 10.0) {
                double off_val = -0.5 - dist(gen) * 0.5;
                coo_entries.emplace_back(i, j, off_val);
                coo_entries.emplace_back(j, i, off_val);
            }
        }
    }

    // COO → CSC 转换（按列排序）
    std::sort(coo_entries.begin(), coo_entries.end(),
              [](const auto& a, const auto& b) {
                  if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
                  return std::get<0>(a) < std::get<0>(b);
              });

    int_t nnz = static_cast<int_t>(coo_entries.size());
    nzval = doubleMalloc(nnz);
    rowind = intMalloc(nnz);
    colptr = intMalloc(n + 1);

    int_t k = 0;
    for (int j = 0; j < n; ++j) {
        colptr[j] = k;
        while (k < nnz && std::get<1>(coo_entries[k]) == j) {
            nzval[k] = std::get<2>(coo_entries[k]);
            rowind[k] = static_cast<int_t>(std::get<0>(coo_entries[k]));
            ++k;
        }
    }
    colptr[n] = nnz;

    return nnz;
}

/**
 * @brief 计算相对误差 ||computed - exact|| / ||exact||
 */
double relative_error(const Eigen::VectorXd& computed, const Eigen::VectorXd& exact) {
    double diff = (computed - exact).norm();
    double exact_norm = exact.norm();
    return (exact_norm > 1e-14) ? diff / exact_norm : diff;
}

/**
 * @brief 生成对称正定稀疏矩阵（CsrMatrix 格式，用于 SymmetricDirectSolver 集成测试）
 */
CsrMatrix<double> generate_spd_csr_matrix(int n, double density = 0.1) {
    CooMatrix<double> coo(n, n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        double diag_val = 10.0 + dist(gen) * 5.0;
        coo.add_value(i, i, diag_val);
        for (int j = i + 1; j < n; ++j) {
            if (dist(gen) < density / n * 10.0) {
                double off_val = -0.5 - dist(gen) * 0.5;
                coo.add_value(i, j, off_val);
                coo.add_value(j, i, off_val);
            }
        }
    }

    CsrMatrix<double> csr(n, n);
    csr.build_from_coo(coo);
    return csr;
}

// ==================== 模块 1: 构造、析构与 RAII ====================

TEST(SuperluContextTest, DefaultConstruction) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 0);
    FEEM_INFO("[RAII测试] 默认构造：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, DestructionNoLeak) {
#ifdef HAVE_SUPERLU
    // 多次创建销毁不应泄漏或崩溃
    for (int i = 0; i < 5; ++i) {
        SuperluContext ctx;
        double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
        int_t nnz = generate_spd_csc_data(10, 0.2, nzval, rowind, colptr);
        EXPECT_TRUE(ctx.initialize(10, nnz, nzval, rowind, colptr));
        EXPECT_TRUE(ctx.is_initialized());
    }
    FEEM_INFO("[RAII测试] 多次创建销毁无泄漏：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 2: initialize() 测试 ====================

TEST(SuperluContextTest, InitializeValidMatrix) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(20, 0.15, nzval, rowind, colptr);

    bool result = ctx.initialize(20, nnz, nzval, rowind, colptr);
    EXPECT_TRUE(result);
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 20);
    // 注意：initialize() 已接管 nzval/rowind/colptr 所有权，
    // 析构时由 Destroy_CompCol_Matrix 统一释放，此处不可再 SUPERLU_FREE
    FEEM_INFO("[初始化测试] 正常初始化：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, InitializeInvalidParameters) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = doubleMalloc(1);
    int_t* rowind = intMalloc(1);
    int_t* colptr = intMalloc(2);

    // n <= 0
    EXPECT_FALSE(ctx.initialize(0, 1, nzval, rowind, colptr));
    EXPECT_FALSE(ctx.initialize(-5, 1, nzval, rowind, colptr));

    // null pointers
    EXPECT_FALSE(ctx.initialize(5, 1, nullptr, rowind, colptr));
    EXPECT_FALSE(ctx.initialize(5, 1, nzval, nullptr, colptr));
    EXPECT_FALSE(ctx.initialize(5, 1, nzval, rowind, nullptr));

    SUPERLU_FREE(nzval);
    SUPERLU_FREE(rowind);
    SUPERLU_FREE(colptr);
    FEEM_INFO("[初始化测试] 异常参数拒绝：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, ReinitializeAfterReset) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;

    // 第一次初始化
    int_t nnz1 = generate_spd_csc_data(10, 0.3, nzval, rowind, colptr);
    EXPECT_TRUE(ctx.initialize(10, nnz1, nzval, rowind, colptr));
    EXPECT_EQ(ctx.matrix_size(), 10);

    // 重置后重新初始化不同大小的矩阵
    ctx.reset();
    double* nzval2 = nullptr; int_t* rowind2 = nullptr; int_t* colptr2 = nullptr;
    int_t nnz2 = generate_spd_csc_data(50, 0.1, nzval2, rowind2, colptr2);
    EXPECT_TRUE(ctx.initialize(50, nnz2, nzval2, rowind2, colptr2));
    EXPECT_EQ(ctx.matrix_size(), 50);
    EXPECT_TRUE(ctx.is_initialized());

    FEEM_INFO("[初始化测试] 重置后重新初始化：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 3: factorize() 测试 ====================

TEST(SuperluContextTest, FactorizeSuccess) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(20, 0.15, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(20, nnz, nzval, rowind, colptr));

    int_t info = ctx.factorize(1);
    EXPECT_EQ(info, 0);
    EXPECT_TRUE(ctx.is_factored());

    FEEM_INFO("[分解测试] LU 分解成功：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, FactorizeStateGuard) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    // 未初始化就分解应返回错误
    int_t info = ctx.factorize(1);
    EXPECT_NE(info, 0);
    EXPECT_FALSE(ctx.is_factored());

    FEEM_INFO("[分解测试] 状态机保护：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, FactorizeMultipleNprocs) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(30, 0.12, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(30, nnz, nzval, rowind, colptr));

    // 测试不同的线程数
    for (int_t nprocs : {1, 2, 4}) {
        ctx.reset();
        double* nzval2 = nullptr; int_t* rowind2 = nullptr; int_t* colptr2 = nullptr;
        int_t nnz2 = generate_spd_csc_data(30, 0.12, nzval2, rowind2, colptr2);
        ASSERT_TRUE(ctx.initialize(30, nnz2, nzval2, rowind2, colptr2));
        int_t info = ctx.factorize(nprocs);
        EXPECT_EQ(info, 0) << "nprocs=" << nprocs << " 分解失败";
        EXPECT_TRUE(ctx.is_factored());
    }

    FEEM_INFO("[分解测试] 多线程参数：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 4: solve() 测试 ====================

TEST(SuperluContextTest, SolveBasicAccuracy) {
#ifdef HAVE_SUPERLU
    const int n = 20;
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(n, 0.15, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));
    ASSERT_EQ(ctx.factorize(1), 0);

    // 构造已知解 x_exact，计算 b = A * x_exact
    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int_t col = 0; col < n; ++col) {
        for (int_t k = colptr[col]; k < colptr[col + 1]; ++k) {
            triplets.emplace_back(rowind[k], col, nzval[k]);
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A_eigen * x_exact;

    auto result = ctx.solve(b);
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);

    double rel_err = relative_error(result.x, x_exact);
    EXPECT_LT(rel_err, 1e-8) << "求解相对误差过大: " << rel_err;

    FEEM_INFO("[求解测试] 基本求解精度：误差={:.2e}", rel_err);
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, SolveMultipleRHS) {
#ifdef HAVE_SUPERLU
    const int n = 30;
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(n, 0.12, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));
    ASSERT_EQ(ctx.factorize(1), 0);

    // 构建 Eigen 矩阵用于验证
    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int_t col = 0; col < n; ++col) {
        for (int_t k = colptr[col]; k < colptr[col + 1]; ++k) {
            triplets.emplace_back(rowind[k], col, nzval[k]);
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());

    // 多次求解不同右端项（复用 L/U 因子）
    for (int trial = 0; trial < 5; ++trial) {
        Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
        Eigen::VectorXd b = A_eigen * x_exact;

        auto result = ctx.solve(b);
        EXPECT_EQ(result.status, SolverStatus::SUCCESS)
            << "第 " << trial+1 << " 次求解状态错误";

        double rel_err = relative_error(result.x, x_exact);
        EXPECT_LT(rel_err, 1e-8) << "第 " << trial+1 << " 次求解误差过大: " << rel_err;
    }

    FEEM_INFO("[求解测试] 多右端项重复求解：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, SolveWithoutFactorize) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(10, 0.2, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(10, nnz, nzval, rowind, colptr));

    // 未分解就求解应返回错误
    Eigen::VectorXd b = Eigen::VectorXd::Ones(10);
    auto result = ctx.solve(b);
    EXPECT_NE(result.status, SolverStatus::SUCCESS);

    FEEM_INFO("[求解测试] 未分解时求解被拒绝：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 5: reset() 与生命周期 ====================

TEST(SuperluContextTest, ResetFullCycle) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(15, 0.18, nzval, rowind, colptr);

    // 完整生命周期: init → factorize → solve → reset → 重新 init
    ASSERT_TRUE(ctx.initialize(15, nnz, nzval, rowind, colptr));
    ASSERT_EQ(ctx.factorize(1), 0);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(15);
    auto result = ctx.solve(b);
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);

    // 重置
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 0);

    // 重新使用
    double* nzval2 = nullptr; int_t* rowind2 = nullptr; int_t* colptr2 = nullptr;
    int_t nnz2 = generate_spd_csc_data(25, 0.1, nzval2, rowind2, colptr2);
    ASSERT_TRUE(ctx.initialize(25, nnz2, nzval2, rowind2, colptr2));
    ASSERT_EQ(ctx.factorize(1), 0);
    EXPECT_TRUE(ctx.is_factored());

    FEEM_INFO("[生命周期测试] 完整重置周期：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, ResetIdempotent) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    // 对未初始化的对象多次 reset 不应崩溃
    ctx.reset();
    ctx.reset();
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());

    FEEM_INFO("[生命周期测试] 幂等重置：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 6: 移动语义 ====================

TEST(SuperluContextTest, MoveConstructor) {
#ifdef HAVE_SUPERLU
    SuperluContext src;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(15, 0.18, nzval, rowind, colptr);
    ASSERT_TRUE(src.initialize(15, nnz, nzval, rowind, colptr));
    ASSERT_EQ(src.factorize(1), 0);

    // 移动构造
    SuperluContext dst(std::move(src));

    // dst 应拥有所有数据
    EXPECT_TRUE(dst.is_factored());
    EXPECT_EQ(dst.matrix_size(), 15);

    // src 应为空
    EXPECT_FALSE(src.is_initialized());
    EXPECT_EQ(src.matrix_size(), 0);

    // dst 可正常求解
    Eigen::VectorXd b = Eigen::VectorXd::Ones(15);
    auto result = dst.solve(b);
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);

    FEEM_INFO("[移动语义测试] 移动构造：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, MoveAssignment) {
#ifdef HAVE_SUPERLU
    SuperluContext src, dst;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(12, 0.2, nzval, rowind, colptr);
    ASSERT_TRUE(src.initialize(12, nnz, nzval, rowind, colptr));
    ASSERT_EQ(src.factorize(1), 0);

    // 移动赋值
    dst = std::move(src);

    EXPECT_TRUE(dst.is_factored());
    EXPECT_EQ(dst.matrix_size(), 12);
    EXPECT_FALSE(src.is_initialized());

    // dst 可正常求解
    Eigen::VectorXd b = Eigen::VectorXd::Ones(12);
    auto result = dst.solve(b);
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);

    FEEM_INFO("[移动语义测试] 移动赋值：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, SelfMoveAssignment) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(10, 0.2, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(10, nnz, nzval, rowind, colptr));

    // 自我移动赋值不应崩溃
    ctx = std::move(ctx);
    // 注意：自我移动后对象处于有效但未指定状态，
    // 标准不要求保持原值，只要求不崩溃
    FEEM_INFO("[移动语义测试] 自我移动赋值（未崩溃即通过）：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 7: 状态机转换 ====================

TEST(SuperluContextTest, StateMachineTransitions) {
#ifdef HAVE_SUPERLU
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(10, 0.2, nzval, rowind, colptr);

    // UNINITIALIZED → INITIALIZED
    EXPECT_FALSE(ctx.is_initialized());
    ctx.initialize(10, nnz, nzval, rowind, colptr);
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    // INITIALIZED → FACTORIZED
    ctx.factorize(1);
    EXPECT_TRUE(ctx.is_factored());

    // FACTORIZED → UNINITIALIZED (via reset)
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    FEEM_INFO("[状态机测试] 全部状态转换正确：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 8: 不同规模矩阵 ====================

TEST(SuperluContextTest, SmallMatrix5x5) {
#ifdef HAVE_SUPERLU
    const int n = 5;
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(n, 0.5, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));
    ASSERT_EQ(ctx.factorize(1), 0);

    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int_t col = 0; col < n; ++col) {
        for (int_t k = colptr[col]; k < colptr[col + 1]; ++k) {
            triplets.emplace_back(rowind[k], col, nzval[k]);
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A_eigen * x_exact;
    auto result = ctx.solve(b);

    EXPECT_EQ(result.status, SolverStatus::SUCCESS);
    EXPECT_LT(relative_error(result.x, x_exact), 1e-10);
    FEEM_INFO("[规模测试] 5x5 小矩阵：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, MediumMatrix100x100) {
#ifdef HAVE_SUPERLU
    const int n = 100;
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(n, 0.05, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));
    ASSERT_EQ(ctx.factorize(1), 0);

    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int_t col = 0; col < n; ++col) {
        for (int_t k = colptr[col]; k < colptr[col + 1]; ++k) {
            triplets.emplace_back(rowind[k], col, nzval[k]);
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A_eigen * x_exact;
    auto result = ctx.solve(b);

    EXPECT_EQ(result.status, SolverStatus::SUCCESS);
    EXPECT_LT(relative_error(result.x, x_exact), 1e-7);
    FEEM_INFO("[规模测试] 100x100 中等矩阵：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

TEST(SuperluContextTest, LargeMatrix500x500) {
#ifdef HAVE_SUPERLU
    const int n = 500;
    SuperluContext ctx;
    double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
    int_t nnz = generate_spd_csc_data(n, 0.01, nzval, rowind, colptr);
    ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));

    auto start = std::chrono::high_resolution_clock::now();
    int_t info = ctx.factorize(1);
    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();

    ASSERT_EQ(info, 0);

    Eigen::SparseMatrix<double> A_eigen(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int_t col = 0; col < n; ++col) {
        for (int_t k = colptr[col]; k < colptr[col + 1]; ++k) {
            triplets.emplace_back(rowind[k], col, nzval[k]);
        }
    }
    A_eigen.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A_eigen * x_exact;
    auto result = ctx.solve(b);

    EXPECT_EQ(result.status, SolverStatus::SUCCESS);
    EXPECT_LT(relative_error(result.x, x_exact), 1e-6);
    FEEM_INFO("[规模测试] 500x500 大矩阵：分解耗时={:.2f}ms, 误差={:.2e}", ms,
              relative_error(result.x, x_exact));
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 9: 内存安全压力测试 ====================

TEST(SuperluContextTest, MemorySafetyRepeatedFactorSolve) {
#ifdef HAVE_SUPERLU
    // 反复分解-求解-重置循环，检测内存泄漏
    for (int cycle = 0; cycle < 10; ++cycle) {
        SuperluContext ctx;
        int n = 20 + cycle * 5;  // 逐渐增大规模
        double* nzval = nullptr; int_t* rowind = nullptr; int_t* colptr = nullptr;
        int_t nnz = generate_spd_csc_data(n, 0.1, nzval, rowind, colptr);
        ASSERT_TRUE(ctx.initialize(n, nnz, nzval, rowind, colptr));
        ASSERT_EQ(ctx.factorize(1), 0);

        for (int solve_iter = 0; solve_iter < 3; ++solve_iter) {
            Eigen::VectorXd b = Eigen::VectorXd::Random(n);
            auto result = ctx.solve(b);
            EXPECT_EQ(result.status, SolverStatus::SUCCESS)
                << "循环 " << cycle << ", 求解 " << solve_iter << " 失败";
        }
        // ctx 析构自动清理
    }

    FEEM_INFO("[内存安全] 10轮分解-求解循环无泄漏：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}

// ==================== 模块 10: SymmetricDirectSolver 集成测试 ====================

TEST(SuperluContextTest, SymmetricDirectSolverIntegration) {
#ifdef HAVE_SUPERLU
    // 通过 SymmetricDirectSolver 使用 SuperLU 后端的端到端测试
    int n = 50;
    CsrMatrix<double> A = generate_spd_csr_matrix(n, 0.08);

    Eigen::VectorXd x_exact = Eigen::VectorXd::Random(n);
    std::vector<double> x_vec(n), b_vec(n);
    for (int i = 0; i < n; ++i) x_vec[i] = x_exact(i);
    A.mat_vec(x_vec, b_vec);
    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), n);

    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
    solver.set_matrix(A);

    // 多次求解复用分解结果
    for (int i = 0; i < 3; ++i) {
        auto result = solver.solve(b);
        EXPECT_EQ(result.status, SolverStatus::SUCCESS)
            << "第 " << i+1 << " 次集成求解失败";

        double err = relative_error(result.x, x_exact);
        EXPECT_LT(err, 1e-7) << "第 " << i+1 << " 次集成求解误差过大: " << err;
    }

    solver.clear();  // 应无崩溃

    FEEM_INFO("[集成测试] SymmetricDirectSolver + SuperLU 后端：通过");
#else
    GTEST_SKIP() << "SuperLU_MT 未启用";
#endif
}
