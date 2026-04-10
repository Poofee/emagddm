/**
 * @file test_superlu_diagnostic.cpp
 * @brief SuperLU_MT 基本功能诊断测试
 * @details 验证 SuperLU_MT 后端的基本 set_matrix + solve 流程
 */

#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include "em_solver_backends.hpp"
#include "em_direct_solvers.h"

using namespace numeric;

/**
 * @brief 构建最简单的3x3三对角矩阵用于诊断
 */
CsrMatrix<double> build_tiny_matrix() {
    int n = 3;
    CooMatrix<double> coo(n, n);
    coo.add_value(0, 0, 4.0);
    coo.add_value(0, 1, -1.0);
    coo.add_value(1, 0, -1.0);
    coo.add_value(1, 1, 4.0);
    coo.add_value(1, 2, -1.0);
    coo.add_value(2, 1, -1.0);
    coo.add_value(2, 2, 4.0);
    CsrMatrix<double> csr(n, n);
    csr.build_from_coo(coo);
    return csr;
}

#ifdef HAVE_SUPERLU

TEST(SuperLUDiagnosticTest, StepByStepDecompose) {
    auto A = build_tiny_matrix();
    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);

    solver.set_matrix(A);

    Eigen::VectorXd b(3);
    b << 3.0, 7.0, 3.0;

    auto result = solver.solve(b);
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);

    Eigen::VectorXd expected(3);
    expected << 1.0, 2.0, 1.0;
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(result.x(i), expected(i), 1e-10);
    }
}

#else
TEST(SuperLUDiagnosticTest, Skipped) { GTEST_SKIP() << "SuperLU未启用"; }
#endif
