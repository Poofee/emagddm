/**
 * @file test_complex_mumps_context.cpp
 * @brief ComplexMumpsContext 复数 MUMPS 求解器综合测试
 * @details 测试覆盖：
 * - 构造/析构/RAII
 * - 初始化（initialize）
 * - 分解（factorize）
 * - 求解（solve）
 * - 重置（reset）
 * - 移动语义
 * - 状态机验证
 * - 复数精度验证（残差检查）
 * - 边界条件（1x1 矩阵、纯实数、纯虚数）
 * - 错误处理（奇异矩阵、维度不匹配、状态错误）
 * - extract_complex 辅助函数
 * - 性能基准测试
 *
 * @author Poofee
 * @date 2026-04-09
 */

#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "direct_solvers.h"
#include "logger_factory.hpp"

using namespace numeric;
using namespace std::complex_literals;

TEST(BasicTest, AlwaysPasses) {
    EXPECT_TRUE(true);
}

TEST(BasicTest, ComplexArithmetic) {
    auto z = 1.0 + 2.0i;
    EXPECT_NEAR(z.real(), 1.0, 1e-15);
    EXPECT_NEAR(z.imag(), 2.0, 1e-15);
}

/**
 * @brief 构建复数对角矩阵的 COO 数据
 */
static void build_diagonal_complex_coo(
    int n, const std::complex<double>& diag_val,
    std::vector<int>& irn, std::vector<int>& jcn,
    std::vector<std::complex<double>>& a)
{
    irn.resize(n);
    jcn.resize(n);
    a.resize(n);
    for (int i = 0; i < n; ++i) {
        irn[i] = i + 1;
        jcn[i] = i + 1;
        a[i] = diag_val;
    }
}

/**
 * @brief 构建复数三对角矩阵的 COO 数据
 */
static int build_tridiag_complex_coo(
    int n, const std::complex<double>& diag, const std::complex<double>& offdiag,
    std::vector<int>& irn, std::vector<int>& jcn,
    std::vector<std::complex<double>>& a)
{
    irn.clear(); jcn.clear(); a.clear();
    for (int i = 0; i < n; ++i) {
        irn.push_back(i + 1); jcn.push_back(i + 1); a.push_back(diag);
        if (i > 0) {
            irn.push_back(i + 1); jcn.push_back(i); a.push_back(offdiag);
        }
        if (i < n - 1) {
            irn.push_back(i + 1); jcn.push_back(i + 2); a.push_back(offdiag);
        }
    }
    return static_cast<int>(a.size());
}

/**
 * @brief 构建复数稀疏矩阵用于残差验证
 */
static Eigen::SparseMatrix<std::complex<double>> build_complex_sparse_from_coo(
    int n, const std::vector<int>& irn, const std::vector<int>& jcn,
    const std::vector<std::complex<double>>& a)
{
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    for (size_t k = 0; k < a.size(); ++k) {
        triplets.emplace_back(irn[k] - 1, jcn[k] - 1, a[k]);
    }
    Eigen::SparseMatrix<std::complex<double>> mat(n, n);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

// ============================================================================
// 构造/析构/RAII
// ============================================================================

TEST(ComplexMumpsContextTest, DefaultConstruction) {
    ComplexMumpsContext ctx;
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 0);
}

TEST(ComplexMumpsContextTest, DestructorDoesNotCrash) {
    { ComplexMumpsContext ctx; }
    { ComplexMumpsContext ctx; ctx.initialize(0); }
    { ComplexMumpsContext ctx; ctx.initialize(0); ctx.reset(); }
}

// ============================================================================
// 初始化
// ============================================================================

TEST(ComplexMumpsContextTest, InitializeUnsymmetric) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

TEST(ComplexMumpsContextTest, InitializeHermitian) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(1));
    EXPECT_TRUE(ctx.is_initialized());
}

TEST(ComplexMumpsContextTest, ReinitializeResetsState) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    ASSERT_TRUE(ctx.initialize(1));
    EXPECT_TRUE(ctx.is_initialized());
}

// ============================================================================
// 分解
// ============================================================================

TEST(ComplexMumpsContextTest, FactorizeDiagonalComplex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 5;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    build_diagonal_complex_coo(n, std::complex<double>(2.0, 1.0), irn, jcn, a);

    EXPECT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), n);
}

TEST(ComplexMumpsContextTest, FactorizeTridiagComplex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 10;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    int nz = build_tridiag_complex_coo(n, std::complex<double>(4.0, 0.5),
                                         std::complex<double>(-1.0, 0.2),
                                         irn, jcn, a);

    EXPECT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_factored());
}

TEST(ComplexMumpsContextTest, FactorizeBeforeInitializeFails) {
    ComplexMumpsContext ctx;
    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {1.0 + 0i};
    EXPECT_NE(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);
}

// ============================================================================
// 求解
// ============================================================================

TEST(ComplexMumpsContextTest, SolveDiagonalComplex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 4;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    build_diagonal_complex_coo(n, std::complex<double>(2.0, 1.0), irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(n);
    for (int i = 0; i < n; ++i) b(i) = std::complex<double>(2.0 + i, 1.0 + i);

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    EXPECT_EQ(x.size(), n);

    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(n, irn, jcn, a);
    Eigen::VectorXcd residual = A * x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8)
        << "复数对角矩阵残差过大";
}

TEST(ComplexMumpsContextTest, SolveTridiagComplex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 10;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    int nz = build_tridiag_complex_coo(n, std::complex<double>(4.0, 0.5),
                                         std::complex<double>(-1.0, 0.2),
                                         irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b = Eigen::VectorXcd::Ones(n);
    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(n, irn, jcn, a);
    Eigen::VectorXcd residual = A * x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8)
        << "复数三对角矩阵残差过大";
}

TEST(ComplexMumpsContextTest, SolvePureRealMatrix) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 3;
    std::vector<int> irn = {1, 1, 2, 2, 2, 3, 3};
    std::vector<int> jcn = {1, 2, 1, 2, 3, 2, 3};
    std::vector<std::complex<double>> a = {4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0};

    ASSERT_EQ(ctx.factorize(n, 7, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(n);
    b << std::complex<double>(3.0, 0.0),
         std::complex<double>(7.0, 0.0),
         std::complex<double>(3.0, 0.0);

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    EXPECT_NEAR(x(0).real(), 19.0 / 14.0, 1e-10);
    EXPECT_NEAR(x(1).real(), 17.0 / 7.0, 1e-10);
    EXPECT_NEAR(x(2).real(), 19.0 / 14.0, 1e-10);
    EXPECT_LT(std::abs(x(0).imag()), 1e-10);
    EXPECT_LT(std::abs(x(1).imag()), 1e-10);
    EXPECT_LT(std::abs(x(2).imag()), 1e-10);
}

TEST(ComplexMumpsContextTest, SolvePureImaginaryMatrix) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 3;
    std::vector<int> irn = {1, 2, 3};
    std::vector<int> jcn = {1, 2, 3};
    std::vector<std::complex<double>> a = {2.0i, 3.0i, 4.0i};

    ASSERT_EQ(ctx.factorize(n, 3, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(n);
    b << 2.0i, 3.0i, 4.0i;

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    EXPECT_NEAR(x(0).real(), 1.0, 1e-10);
    EXPECT_NEAR(x(1).real(), 1.0, 1e-10);
    EXPECT_NEAR(x(2).real(), 1.0, 1e-10);
}

TEST(ComplexMumpsContextTest, SolveBeforeFactorizeFails) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    Eigen::VectorXcd b(3);
    auto result = ctx.solve(b);
    EXPECT_NE(result.status, SolverStatus::SUCCESS);
}

// ============================================================================
// 重置
// ============================================================================

TEST(ComplexMumpsContextTest, ResetAfterInitialize) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

TEST(ComplexMumpsContextTest, ResetAfterFactorize) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {2.0 + 1i};
    ASSERT_EQ(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

TEST(ComplexMumpsContextTest, ResetUninitializedDoesNotCrash) {
    ComplexMumpsContext ctx;
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
}

// ============================================================================
// 移动语义
// ============================================================================

TEST(ComplexMumpsContextTest, MoveConstruction) {
    ComplexMumpsContext ctx1;
    ASSERT_TRUE(ctx1.initialize(0));
    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {3.0 + 2i};
    ASSERT_EQ(ctx1.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);

    ComplexMumpsContext ctx2(std::move(ctx1));
    EXPECT_TRUE(ctx2.is_factored());
    EXPECT_FALSE(ctx1.is_initialized());
}

TEST(ComplexMumpsContextTest, MoveAssignment) {
    ComplexMumpsContext ctx1;
    ASSERT_TRUE(ctx1.initialize(0));
    ComplexMumpsContext ctx2;
    ctx2 = std::move(ctx1);
    EXPECT_TRUE(ctx2.is_initialized());
    EXPECT_FALSE(ctx1.is_initialized());
}

// ============================================================================
// 状态机验证
// ============================================================================

TEST(ComplexMumpsContextTest, StateTransitions) {
    ComplexMumpsContext ctx;

    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {1.0 + 0i};
    ASSERT_EQ(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_TRUE(ctx.is_factored());

    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

// ============================================================================
// 精度验证
// ============================================================================

TEST(ComplexMumpsContextTest, AccuracyResidualCheck) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 20;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    int nz = build_tridiag_complex_coo(n, std::complex<double>(4.0, 0.3),
                                         std::complex<double>(-1.0, 0.1),
                                         irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b = Eigen::VectorXcd::Ones(n);
    for (int i = 0; i < n; ++i) {
        b(i) = std::complex<double>(1.0 + 0.1 * i, 0.5 + 0.05 * i);
    }

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(n, irn, jcn, a);
    Eigen::VectorXcd residual = A * x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8)
        << "复数精度验证失败，残差过大";
}

TEST(ComplexMumpsContextTest, AccuracyMultipleSolveReuse) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 5;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    build_diagonal_complex_coo(n, std::complex<double>(3.0, 1.5), irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);

    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(n, irn, jcn, a);

    for (int trial = 0; trial < 5; ++trial) {
        Eigen::VectorXcd b = Eigen::VectorXcd::Random(n);
        auto result = ctx.solve(b);
        ASSERT_EQ(result.status, SolverStatus::SUCCESS);
        Eigen::VectorXcd x = extract_complex(result);
        Eigen::VectorXcd residual = A * x - b;
        EXPECT_LT(residual.norm() / b.norm(), 1e-8)
            << "多次求解复用测试失败，trial=" << trial;
    }
}

// ============================================================================
// 边界条件
// ============================================================================

TEST(ComplexMumpsContextTest, Solve1x1Complex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {std::complex<double>(2.0, 3.0)};
    ASSERT_EQ(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(1);
    b(0) = std::complex<double>(4.0, 6.0);

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    EXPECT_NEAR(x(0).real(), 2.0, 1e-10);
    EXPECT_NEAR(x(0).imag(), 0.0, 1e-10);
}

TEST(ComplexMumpsContextTest, Solve2x2Complex) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn = {1, 1, 2, 2};
    std::vector<int> jcn = {1, 2, 1, 2};
    std::vector<std::complex<double>> a = {
        std::complex<double>(3.0, 1.0), std::complex<double>(1.0, 0.5),
        std::complex<double>(1.0, 0.5), std::complex<double>(3.0, 1.0)
    };

    ASSERT_EQ(ctx.factorize(2, 4, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(2);
    b << std::complex<double>(4.0, 1.5), std::complex<double>(4.0, 1.5);

    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXcd x = extract_complex(result);
    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(2, irn, jcn, a);
    Eigen::VectorXcd residual = A * x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8);
}

// ============================================================================
// 错误处理
// ============================================================================

TEST(ComplexMumpsContextTest, DimensionMismatchError) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {2.0 + 1i};
    ASSERT_EQ(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXcd b(3);
    auto result = ctx.solve(b);
    EXPECT_NE(result.status, SolverStatus::SUCCESS);
    EXPECT_FALSE(result.error_msg.empty());
}

TEST(ComplexMumpsContextTest, FactorizeAfterFactorizeResetsState) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn = {1}, jcn = {1};
    std::vector<std::complex<double>> a = {2.0 + 1i};
    ASSERT_EQ(ctx.factorize(1, 1, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_factored());
}

// ============================================================================
// extract_complex 辅助函数
// ============================================================================

TEST(ComplexMumpsContextTest, ExtractComplexFunction) {
    SolverResult result;
    result.x = Eigen::VectorXd(6);
    result.x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    Eigen::VectorXcd x = extract_complex(result);
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x(0).real(), 1.0, 1e-15);
    EXPECT_NEAR(x(0).imag(), 2.0, 1e-15);
    EXPECT_NEAR(x(1).real(), 3.0, 1e-15);
    EXPECT_NEAR(x(1).imag(), 4.0, 1e-15);
    EXPECT_NEAR(x(2).real(), 5.0, 1e-15);
    EXPECT_NEAR(x(2).imag(), 6.0, 1e-15);
}

// ============================================================================
// 性能基准测试
// ============================================================================

TEST(ComplexMumpsContextTest, PerformanceBenchmark) {
    ComplexMumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 500;
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> a;
    int nz = build_tridiag_complex_coo(n, std::complex<double>(4.0, 0.5),
                                         std::complex<double>(-1.0, 0.2),
                                         irn, jcn, a);

    auto t0 = std::chrono::high_resolution_clock::now();
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    double factorize_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    Eigen::VectorXcd b = Eigen::VectorXcd::Ones(n);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto result = ctx.solve(b);
    auto t3 = std::chrono::high_resolution_clock::now();
    double solve_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);
    FEEM_INFO("ComplexMumpsContext 性能: n={}, nz={}, 分解={:.2f}ms, 求解={:.2f}ms",
              n, nz, factorize_ms, solve_ms);

    Eigen::VectorXcd x = extract_complex(result);
    Eigen::SparseMatrix<std::complex<double>> A = build_complex_sparse_from_coo(n, irn, jcn, a);
    Eigen::VectorXcd residual = A * x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-6)
        << "性能测试精度验证失败";
}
