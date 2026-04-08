/**
 * @file test_superlu_diagnostic.cpp
 * @brief SuperLU_MT 崩溃诊断测试
 * @details 逐步调用SuperLU_MT API，定位SEH异常(0xc0000005)的具体崩溃点
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

#if EM_SOLVER_HAS_SUPERLU
    #include "slu_mt_ddefs.h"
    #include <windows.h>  // 用于 SEH 异常处理
#endif

using namespace numeric;

// 构建最简单的3x3矩阵用于诊断
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

#if EM_SOLVER_HAS_SUPERLU

TEST(SuperLUDiagnosticTest, StepByStepDecompose) {
    auto A = build_tiny_matrix();
    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);

    // 步骤1: set_matrix（内部会调用decompose）
    // 捕获SEH异常以获取更多信息
    __try {
        std::cout << "[DIAG] 调用 set_matrix()..." << std::endl;
        solver.set_matrix(A);
        std::cout << "[DIAG] set_matrix() 成功" << std::endl;
    } __except(EXCEPTION_EXECUTE_HANDLER) {
        std::cout << "[DIAG] SEH异常在set_matrix(): code=0x"
                  << std::hex << GetExceptionCode() << std::endl;
        FAIL() << "set_matrix() 崩溃";
    }

    // 如果到这里，说明分解成功，尝试求解
    Eigen::VectorXd b(3);
    b << 3.0, 7.0, 3.0;

    __try {
        std::cout << "[DIAG] 调用 solve()..." << std::endl;
        auto result = solver.solve(b);
        std::cout << "[DIAG] solve() 成功, status=" << (int)result.status << std::endl;
        EXPECT_EQ(result.status, SolverStatus::SUCCESS);
    } __except(EXCEPTION_EXECUTE_HANDLER) {
        std::cout << "[DIAG] SEH异常在solve(): code=0x"
                  << std::hex << GetExceptionCode() << std::endl;
        FAIL() << "solve() 崩溃";
    }
}

#else
TEST(SuperLUDiagnosticTest, Skipped) { GTEST_SKIP() << "SuperLU未启用"; }
#endif
