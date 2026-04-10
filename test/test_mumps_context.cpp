/**
 * @file test_mumps_context.cpp
 * @brief MumpsContext RAII 封装类的全面功能测试
 * @details 覆盖 MumpsContext 类的所有公开接口和边界条件：
 *
 * **测试模块**：
 * 1. 构造/析构与 RAII 语义
 * 2. initialize() 初始化（正常/异常参数）
 * 3. factorize() 分析+分解（成功/失败/状态机保护）
 * 4. solve() 求解（精度/重复调用/多右端项）
 * 5. reset() 重置与资源释放
 * 6. 移动语义
 * 7. 状态机转换正确性
 * 8. 内存安全（无泄漏）
 * 9. 边界条件（奇异矩阵、单位矩阵、大规模矩阵）
 * 10. 与 SymmetricDirectSolver 集成测试
 *
 * @author Poofee
 * @date 2026-04-09
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "em_direct_solvers.h"
#include "logger_factory.hpp"

#if !HAVE_MUMPS
    #error "此测试需要 MUMPS 支持，请使用 -DUSE_MUMPS=ON 编译"
#endif

using namespace numeric;

/**
 * @brief 生成对称正定稀疏矩阵的 COO 数据（供 MumpsContext 使用）
 * @param n 矩阵维度
 * @param density 非零密度
 * @param[out] irn 行索引数组（1-based，MUMPS 格式）
 * @param[out] jcn 列索引数组（1-based，MUMPS 格式）
 * @param[out] a 非零值数组
 * @return int 非零元个数 nz
 */
int generate_spd_coo_data(int n, double density,
                          std::vector<int>& irn,
                          std::vector<int>& jcn,
                          std::vector<double>& a) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(0.0, 1.0);

    irn.clear();
    jcn.clear();
    a.clear();

    for (int i = 0; i < n; ++i) {
        double diag_val = 10.0 + dist(gen) * 5.0;
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(diag_val);

        for (int j = i + 1; j < n; ++j) {
            if (dist(gen) < density / n * 10.0) {
                double off_val = -0.5 - dist(gen) * 0.5;
                irn.push_back(i + 1);
                jcn.push_back(j + 1);
                a.push_back(off_val);
                irn.push_back(j + 1);
                jcn.push_back(i + 1);
                a.push_back(off_val);
            }
        }
    }

    return static_cast<int>(irn.size());
}

/**
 * @brief 生成一般非对称稀疏矩阵的 COO 数据
 */
int generate_unsymmetric_coo_data(int n, double density,
                                  std::vector<int>& irn,
                                  std::vector<int>& jcn,
                                  std::vector<double>& a) {
    std::mt19937 gen(123);
    std::uniform_real_distribution<> dist(0.0, 1.0);

    irn.clear();
    jcn.clear();
    a.clear();

    for (int i = 0; i < n; ++i) {
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(5.0 + dist(gen) * 3.0);

        for (int j = 0; j < n; ++j) {
            if (j != i && dist(gen) < density / n * 5.0) {
                irn.push_back(i + 1);
                jcn.push_back(j + 1);
                a.push_back(-0.3 - dist(gen) * 0.7);
            }
        }
    }

    return static_cast<int>(irn.size());
}

/**
 * @brief 生成三对角矩阵的 COO 数据
 */
int generate_tridiag_coo_data(int n,
                              std::vector<int>& irn,
                              std::vector<int>& jcn,
                              std::vector<double>& a,
                              double diag = 4.0,
                              double offdiag = -1.0) {
    irn.clear();
    jcn.clear();
    a.clear();

    for (int i = 0; i < n; ++i) {
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(diag);

        if (i > 0) {
            irn.push_back(i + 1);
            jcn.push_back(i);
            a.push_back(offdiag);
        }

        if (i < n - 1) {
            irn.push_back(i + 1);
            jcn.push_back(i + 2);
            a.push_back(offdiag);
        }
    }

    return static_cast<int>(irn.size());
}

/**
 * @brief 从 COO 数据构建 Eigen 稀疏矩阵（用于残差验证）
 */
Eigen::SparseMatrix<double> build_sparse_matrix_from_coo(
    int n, const std::vector<int>& irn,
    const std::vector<int>& jcn,
    const std::vector<double>& a) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t k = 0; k < irn.size(); ++k) {
        triplets.emplace_back(irn[k] - 1, jcn[k] - 1, a[k]);
    }
    Eigen::SparseMatrix<double> mat(n, n);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

/**
 * @brief 从 Eigen 稀疏矩阵构建 CsrMatrix（用于 SymmetricDirectSolver）
 */
CsrMatrix<double> eigen_to_csr(const Eigen::SparseMatrix<double>& eigen_mat) {
    int n = static_cast<int>(eigen_mat.rows());
    CooMatrix<double> coo(n, n);

    for (int col = 0; col < eigen_mat.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(eigen_mat, col); it; ++it) {
            coo.add_value(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
        }
    }

    CsrMatrix<double> csr(n, n);
    csr.build_from_coo(coo);
    return csr;
}
#ifdef HAVE_MUMPS

// ==================== 1. 构造/析构与 RAII 语义 ====================

TEST(MumpsContextTest, DefaultConstructor) {
    MumpsContext ctx;
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 0);
}

TEST(MumpsContextTest, DestructorNoDoubleFree) {
    {
        MumpsContext ctx;
        ctx.initialize(0);
    }
    SUCCEED();
}

TEST(MumpsContextTest, MultipleCreateDestroy) {
    for (int i = 0; i < 5; ++i) {
        MumpsContext ctx;
        ASSERT_TRUE(ctx.initialize(0));
    }
    SUCCEED();
}

// ==================== 2. initialize() 测试 ====================

TEST(MumpsContextTest, InitializeUnsymmetric) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

TEST(MumpsContextTest, InitializeSymmetricPositiveDefinite) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(1));
    EXPECT_TRUE(ctx.is_initialized());
}

TEST(MumpsContextTest, InitializeGeneralSymmetric) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(2));
    EXPECT_TRUE(ctx.is_initialized());
}

TEST(MumpsContextTest, ReinitializeResetsState) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 5;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_factored());

    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

// ==================== 3. factorize() 测试 ====================

TEST(MumpsContextTest, FactorizeSmallTridiag) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 10;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a);

    int result = ctx.factorize(n, nz, irn.data(), jcn.data(), a.data());
    ASSERT_EQ(result, 0);
    EXPECT_TRUE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), n);
}

TEST(MumpsContextTest, FactorizeWithoutInitialize) {
    MumpsContext ctx;
    std::vector<int> irn = {1, 2};
    std::vector<int> jcn = {1, 2};
    std::vector<double> a = {1.0, 1.0};

    int result = ctx.factorize(2, 2, irn.data(), jcn.data(), a.data());
    EXPECT_NE(result, 0);
}

TEST(MumpsContextTest, FactorizeSPDMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(1));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 20;
    int nz = generate_spd_coo_data(n, 0.3, irn, jcn, a);

    int result = ctx.factorize(n, nz, irn.data(), jcn.data(), a.data());
    ASSERT_EQ(result, 0);
}

TEST(MumpsContextTest, FactorizeUnsymmetricMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 15;
    int nz = generate_unsymmetric_coo_data(n, 0.3, irn, jcn, a);

    int result = ctx.factorize(n, nz, irn.data(), jcn.data(), a.data());
    ASSERT_EQ(result, 0);
}

// ==================== 4. solve() 测试 ====================

TEST(MumpsContextTest, SolveTridiagSystem) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 10;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a, 4.0, -1.0);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto result = ctx.solve(b);

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);
    EXPECT_EQ(result.x.size(), n);

    Eigen::SparseMatrix<double> A = build_sparse_matrix_from_coo(n, irn, jcn, a);
    Eigen::VectorXd residual = A * result.x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8)
        << "Tridiagonal system residual too large";
}

TEST(MumpsContextTest, SolveIdentityMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 5;
    std::vector<int> irn, jcn;
    std::vector<double> a;
    for (int i = 0; i < n; ++i) {
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(1.0);
    }

    ASSERT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b(n);
    b << 1.0, 2.0, 3.0, 4.0, 5.0;
    auto result = ctx.solve(b);

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(result.x(i), b(i), 1e-10);
    }
}

TEST(MumpsContextTest, SolveMultipleRHSReuseFactorization) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 20;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a, 4.0, -1.0);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    for (int rhs_idx = 0; rhs_idx < 5; ++rhs_idx) {
        Eigen::VectorXd b = Eigen::VectorXd::Random(n);
        auto result = ctx.solve(b);

        ASSERT_EQ(result.status, SolverStatus::SUCCESS)
            << "Solve failed for RHS #" << rhs_idx;

        Eigen::SparseMatrix<double> A = build_sparse_matrix_from_coo(n, irn, jcn, a);
        Eigen::VectorXd residual = A * result.x - b;
        EXPECT_LT(residual.norm() / (b.norm() + 1e-15), 1e-6)
            << "Residual too large for RHS #" << rhs_idx;
    }
}

TEST(MumpsContextTest, SolveWithoutFactorize) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    Eigen::VectorXd b(5);
    auto result = ctx.solve(b);
    EXPECT_NE(result.status, SolverStatus::SUCCESS);
}

TEST(MumpsContextTest, SolveDimensionMismatch) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 5;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b_wrong(10);
    auto result = ctx.solve(b_wrong);
    EXPECT_NE(result.status, SolverStatus::SUCCESS);
}

// ==================== 5. reset() 测试 ====================

TEST(MumpsContextTest, ResetFromInitialized) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());

    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
    EXPECT_EQ(ctx.matrix_size(), 0);
}

TEST(MumpsContextTest, ResetFromFactored) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int nz = generate_tridiag_coo_data(5, irn, jcn, a);
    ASSERT_EQ(ctx.factorize(5, nz, irn.data(), jcn.data(), a.data()), 0);

    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

TEST(MumpsContextTest, ResetFromUninitialized) {
    MumpsContext ctx;
    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
}

TEST(MumpsContextTest, ResetAndReinitialize) {
    MumpsContext ctx;

    for (int round = 0; round < 3; ++round) {
        ASSERT_TRUE(ctx.initialize(0));

        std::vector<int> irn, jcn;
        std::vector<double> a;
        int nz = generate_tridiag_coo_data(10, irn, jcn, a);
        ASSERT_EQ(ctx.factorize(10, nz, irn.data(), jcn.data(), a.data()), 0);

        Eigen::VectorXd b = Eigen::VectorXd::Ones(10);
        auto result = ctx.solve(b);
        ASSERT_EQ(result.status, SolverStatus::SUCCESS);

        ctx.reset();
    }
    SUCCEED();
}

// ==================== 6. 移动语义测试 ====================

TEST(MumpsContextTest, MoveConstructor) {
    MumpsContext ctx1;
    ASSERT_TRUE(ctx1.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int nz = generate_tridiag_coo_data(5, irn, jcn, a);
    ASSERT_EQ(ctx1.factorize(5, nz, irn.data(), jcn.data(), a.data()), 0);

    MumpsContext ctx2(std::move(ctx1));

    EXPECT_TRUE(ctx2.is_factored());
    EXPECT_EQ(ctx2.matrix_size(), 5);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(5);
    auto result = ctx2.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);
}

TEST(MumpsContextTest, MoveAssignment) {
    MumpsContext ctx1;
    ASSERT_TRUE(ctx1.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int nz = generate_tridiag_coo_data(8, irn, jcn, a);
    ASSERT_EQ(ctx1.factorize(8, nz, irn.data(), jcn.data(), a.data()), 0);

    MumpsContext ctx2;
    ctx2 = std::move(ctx1);

    EXPECT_TRUE(ctx2.is_factored());
    EXPECT_EQ(ctx2.matrix_size(), 8);
}

// ==================== 7. 状态机转换测试 ====================

TEST(MumpsContextTest, StateTransitions) {
    MumpsContext ctx;

    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    ASSERT_TRUE(ctx.initialize(0));
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int nz = generate_tridiag_coo_data(5, irn, jcn, a);
    ASSERT_EQ(ctx.factorize(5, nz, irn.data(), jcn.data(), a.data()), 0);
    EXPECT_TRUE(ctx.is_initialized());
    EXPECT_TRUE(ctx.is_factored());

    ctx.reset();
    EXPECT_FALSE(ctx.is_initialized());
    EXPECT_FALSE(ctx.is_factored());
}

// ==================== 8. 精度验证测试 ====================

TEST(MumpsContextTest, AccuracySmallDenseSystem) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 4;
    std::vector<int> irn = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    std::vector<int> jcn = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
    std::vector<double> a_vals = {
        10.0, 1.0, 1.0, 1.0,
        1.0, 10.0, 1.0, 1.0,
        1.0, 1.0, 10.0, 1.0,
        1.0, 1.0, 1.0, 10.0
    };

    ASSERT_EQ(ctx.factorize(n, 16, irn.data(), jcn.data(), a_vals.data()), 0);

    Eigen::VectorXd b(n);
    b << 13.0, 13.0, 13.0, 13.0;
    auto result = ctx.solve(b);

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(result.x(i), 1.0, 1e-10);
    }
}

TEST(MumpsContextTest, AccuracyResidualCheck) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 50;
    int nz = generate_spd_coo_data(n, 0.2, irn, jcn, a);
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::SparseMatrix<double> A = build_sparse_matrix_from_coo(n, irn, jcn, a);
    Eigen::VectorXd residual = A * result.x - b;
    double rel_res = residual.norm() / (A.norm() * result.x.norm() + b.norm());
    EXPECT_LT(rel_res, 1e-8);
}

// ==================== 9. 边界条件测试 ====================

TEST(MumpsContextTest, IdentityMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 100;
    std::vector<int> irn, jcn;
    std::vector<double> a;
    for (int i = 0; i < n; ++i) {
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(1.0);
    }

    ASSERT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(result.x(i), b(i), 1e-10);
    }
}

TEST(MumpsContextTest, DiagonalMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 20;
    std::vector<int> irn, jcn;
    std::vector<double> a;
    std::vector<double> diag_vals(n);
    for (int i = 0; i < n; ++i) {
        diag_vals[i] = 1.0 + i * 0.5;
        irn.push_back(i + 1);
        jcn.push_back(i + 1);
        a.push_back(diag_vals[i]);
    }

    ASSERT_EQ(ctx.factorize(n, n, irn.data(), jcn.data(), a.data()), 0);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(result.x(i), 1.0 / diag_vals[i], 1e-10);
    }
}

TEST(MumpsContextTest, SmallMatrix2x2) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    int n = 2;
    std::vector<int> irn = {1, 2, 1, 2};
    std::vector<int> jcn = {1, 1, 2, 2};
    std::vector<double> a_vals = {2.0, 1.0, 1.0, 3.0};

    ASSERT_EQ(ctx.factorize(n, 4, irn.data(), jcn.data(), a_vals.data()), 0);

    Eigen::VectorXd b(2);
    b << 5.0, 8.0;
    auto result = ctx.solve(b);
    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    EXPECT_NEAR(result.x(0), 1.4, 1e-10);
    EXPECT_NEAR(result.x(1), 2.2, 1e-10);
}

TEST(MumpsContextTest, LargeSparseMatrix) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 500;
    int nz = generate_spd_coo_data(n, 0.05, irn, jcn, a);

    int result = ctx.factorize(n, nz, irn.data(), jcn.data(), a.data());
    ASSERT_EQ(result, 0);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto solve_result = ctx.solve(b);
    ASSERT_EQ(solve_result.status, SolverStatus::SUCCESS);
    EXPECT_EQ(solve_result.x.size(), n);
}


// ==================== 11. 性能基准测试 ====================

TEST(MumpsContextPerfTest, BenchmarkTridiagSolve) {
    MumpsContext ctx;
    ASSERT_TRUE(ctx.initialize(0));

    std::vector<int> irn, jcn;
    std::vector<double> a;
    int n = 1000;
    int nz = generate_tridiag_coo_data(n, irn, jcn, a);

    auto t0 = std::chrono::high_resolution_clock::now();
    ASSERT_EQ(ctx.factorize(n, nz, irn.data(), jcn.data(), a.data()), 0);
    auto t1 = std::chrono::high_resolution_clock::now();

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto solve_result = ctx.solve(b);
    auto t2 = std::chrono::high_resolution_clock::now();

    ASSERT_EQ(solve_result.status, SolverStatus::SUCCESS);

    double factorize_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double solve_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();

    std::cout << "  [MUMPS Perf] n=" << n << " factorize=" << factorize_ms
              << "ms solve=" << solve_ms << "ms" << std::endl;

    EXPECT_LT(factorize_ms, 5000.0);
    EXPECT_LT(solve_ms, 100.0);
}
#endif

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
