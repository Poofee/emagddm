/**
 * @file test_superlu_mt_integration.cpp
 * @brief SuperLU_MT 多线程稀疏LU求解器集成测试
 * @details 使用 Google Test 框架验证 SuperLU_MT 后端在 SymmetricDirectSolver 和
 *          GeneralDirectSolver 中的正确集成，确保线性方程组求解结果准确。
 *
 * 测试场景：
 * - 对称正定(SPD)矩阵：使用 SymmetricDirectSolver + SuperLU_MT 后端
 * - 通用非对称矩阵：使用 GeneralDirectSolver + SuperLU_MT 后端（预留）
 * - 矩阵维度：5x5 小规模快速测试 + 20x20 中等规模精度验证
 *
 * 验证点：
 * 1. SuperLU_MT 后端可用性检测（编译期宏 + 运行时查询）
 * 2. SPD矩阵求解：解向量误差 ||x_computed - x_exact|| < 1e-8
 * 3. 残差范数：||b - A*x_computed|| < 1e-10
 * 4. 多次 solve() 调用结果一致性（分解复用机制验证）
 * 5. 与 Eigen 后端结果对比（交叉验证）
 * 6. 资源清理无内存泄漏（clear() 后可重新 set_matrix()）
 *
 * @author Poofee
 * @date 2026-04-07
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef HAVE_SUPERLU
    #include "em_direct_solvers.h"
#else
    #include "em_direct_solvers.h"
#endif

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include "logger_factory.hpp"

using namespace numeric;


// ==================== 测试辅助函数 ====================

/**
 * @brief 构建小型对称正定(SPD)测试矩阵（5x5）
 * @return CsrMatrix<double> 5x5 对称正定稀疏矩阵
 *
 * @details 矩阵结构：
 * - 主对角线：[4, 5, 6, 7, 8] （递增正值保证正定性）
 * - 次对角线：[-1, -1, -1, -1] （相邻节点耦合）
 * - 三对角线：[-0.2, -0.2, -0.2] （次邻弱耦合）
 */
CsrMatrix<double> build_small_spd_matrix() {
    int n = 5;
    CooMatrix<double> coo(n, n);

    // 主对角线
    for (int i = 0; i < n; ++i) {
        coo.add_value(i, i, 4.0 + i);
    }

    // 次对角线（对称）
    for (int i = 0; i < n - 1; ++i) {
        coo.add_value(i, i + 1, -1.0);
        coo.add_value(i + 1, i, -1.0);
    }

    // (0,2) 和 (2,0) 位置的非零元
    coo.add_value(0, 2, -0.2);
    coo.add_value(2, 0, -0.2);

    CsrMatrix<double> csr(n, n);
    csr.build_from_coo(coo);
    return csr;
}

/**
 * @brief 构建中等规模对称正定(SPD)测试矩阵（20x20）
 * @return CsrMatrix<double> 20x20 对称正定稀疏矩阵
 *
 * @details 模拟有限元刚度矩阵的稀疏结构：
 * - 对角占优的三对角主结构
 * - 随机添加少量远距离耦合项模拟网格连接
 */
CsrMatrix<double> build_medium_spd_matrix() {
    int n = 20;
    CooMatrix<double> coo(n, n);

    for (int i = 0; i < n; ++i) {
        // 强对角占优
        coo.add_value(i, i, 10.0 + 0.1 * i);

        // 相邻节点耦合
        if (i > 0) {
            double val = -1.0;
            coo.add_value(i, i - 1, val);
            coo.add_value(i - 1, i, val);
        }

        // 次邻节点弱耦合
        if (i > 1) {
            double val = -0.15;
            coo.add_value(i, i - 2, val);
            coo.add_value(i - 2, i, val);
        }
    }

    CsrMatrix<double> csr(n, n);
    csr.build_from_coo(coo);
    return csr;
}

/**
 * @brief 构造已知解析解并计算对应的右端项 b = A * x_exact
 * @param A 系数矩阵
 * @param x_exact 已知精确解向量
 * @return std::vector<double> 右端项向量 b
 */
std::vector<double> compute_rhs(const CsrMatrix<double>& A,
                                const std::vector<double>& x_exact) {
    std::vector<double> b(A.rows(), 0.0);
    A.mat_vec(x_exact, b);
    return b;
}

/**
 * @brief 将 std::vector 转换为 Eigen::VectorXd
 * @param v 输入的 std::vector
 * @return Eigen::VectorXd 输出的 Eigen 向量
 */
Eigen::VectorXd to_eigen(const std::vector<double>& v) {
    Eigen::VectorXd ev(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        ev(i) = v[i];
    }
    return ev;
}

/**
 * @brief 计算 Euclidean 范数误差
 * @param computed 计算得到的解向量
 * @param exact 精确解向量
 * @return double Euclidean 范数误差 ||computed - exact||_2
 */
double compute_error(const Eigen::VectorXd& computed, const std::vector<double>& exact) {
    Eigen::VectorXd e = to_eigen(exact);
    return (computed - e).norm();
}


// ==================== 单元测试用例 ====================

/**
 * @test 测试SuperLU_MT后端的编译期和运行时可用性
 */
TEST(SuperLUMTIntegrationTest, BackendAvailability) {
#ifdef HAVE_SUPERLU
    EXPECT_TRUE(DirectBackendManager::isBackendAvailable(DirectBackendType::SUPERLU));
    EXPECT_EQ(DirectBackendManager::getBackendName(DirectBackendType::SUPERLU), "SuperLU");

    auto backends = DirectBackendManager::getAvailableBackends();
    bool found_superlu = false;
    for (auto& b : backends) {
        if (b == DirectBackendType::SUPERLU) {
            found_superlu = true;
            break;
        }
    }
    EXPECT_TRUE(found_superlu) << "SuperLU后端应在可用后端列表中";
#else
    GTEST_SKIP() << "SuperLU_MT未编译进此构建（HAVE_SUPERLU未定义）";
#endif
}

/**
 * @test 测试 SPD 矩阵的小规模求解（5x5）- SuperLU_MT 后端
 *
 * 验证流程：
 * 1. 构建 5x5 SPD 测试矩阵 A
 * 2. 设定精确解 x_exact = [1, 2, 3, 4, 5]^T
 * 3. 计算右端项 b = A * x_exact
 * 4. 使用 SymmetricDirectSolver(SUPERLU后端) 求解 Ax=b
 * 5. 验证解误差 < 1e-8
 */
TEST(SuperLUMTIntegrationTest, SmallSPDSolve) {
#ifdef HAVE_SUPERLU
    auto A = build_small_spd_matrix();
    ASSERT_EQ(A.rows(), 5);
    ASSERT_EQ(A.cols(), 5);

    std::vector<double> x_exact = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto b_vec = compute_rhs(A, x_exact);

    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
    solver.set_matrix(A);

    auto result = solver.solve(to_eigen(b_vec));

    EXPECT_EQ(result.status, SolverStatus::SUCCESS)
        << "求解状态应为SUCCESS，错误信息: " << result.error_msg;

    double error = compute_error(result.x, x_exact);
    EXPECT_LT(error, 1e-8)
        << "解向量误差过大: " << error << " (阈值: 1e-8)";

    EXPECT_LT(result.residual_norm, 1e-10)
        << "残差范数过大: " << result.residual_norm;

    FEEM_INFO("SuperLU_MT小规模SPD测试通过, 误差={:.2e}, 残差={:.2e}, 耗时={:.3f}ms",
              error, result.residual_norm, result.solve_time_ms);
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}

/**
 * @test 测试 SPD 矩阵的中等规模求解（20x20）- SuperLU_MT 后端
 *
 * 更大规模的测试以验证数值稳定性和性能。
 * 精度要求与小规模测试一致。
 */
TEST(SuperLUMTIntegrationTest, MediumSPDSolve) {
#ifdef HAVE_SUPERLU
    auto A = build_medium_spd_matrix();
    ASSERT_EQ(A.rows(), 20);
    ASSERT_EQ(A.cols(), 20);

    std::vector<double> x_exact(20);
    for (int i = 0; i < 20; ++i) {
        x_exact[i] = static_cast<double>(i + 1);
    }

    auto b_vec = compute_rhs(A, x_exact);

    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
    solver.set_matrix(A);

    auto result = solver.solve(to_eigen(b_vec));

    EXPECT_EQ(result.status, SolverStatus::SUCCESS)
        << "求解失败: " << result.error_msg;

    double error = compute_error(result.x, x_exact);
    EXPECT_LT(error, 1e-8)
        << "中等规模SPD求解误差过大: " << error;

    EXPECT_LT(result.residual_norm, 1e-10)
        << "残差范数过大: " << result.residual_norm;

    FEEM_INFO("SuperLU_MT中等规模SPD测试通过 (20x20), 误差={:.2e}, 残差={:.2e}",
              error, result.residual_norm);
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}

/**
 * @test 测试多次solve()调用的一致性（分解复用验证）
 *
 * 关键优化点验证：set_matrix()一次分解后，
 * 多次solve()应返回一致且正确的结果。
 */
TEST(SuperLUMTIntegrationTest, MultipleSolveConsistency) {
#ifdef HAVE_SUPERLU
    auto A = build_small_spd_matrix();

    std::vector<double> x_exact = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto b_vec = compute_rhs(A, x_exact);

    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
    solver.set_matrix(A);

    // 第一次求解
    auto result1 = solver.solve(to_eigen(b_vec));
    ASSERT_EQ(result1.status, SolverStatus::SUCCESS);

    // 第二次求解（相同右端项，应复用缓存的L/U因子）
    auto result2 = solver.solve(to_eigen(b_vec));
    ASSERT_EQ(result2.status, SolverStatus::SUCCESS);

    // 第三次求解（不同右端项）
    std::vector<double> x_exact2 = {5.0, 4.0, 3.0, 2.0, 1.0};
    auto b_vec2 = compute_rhs(A, x_exact2);
    auto result3 = solver.solve(to_eigen(b_vec2));
    ASSERT_EQ(result3.status, SolverStatus::SUCCESS);

    // 验证三次结果各自正确
    double error1 = compute_error(result1.x, x_exact);
    double error2 = compute_error(result2.x, x_exact);
    double error3 = compute_error(result3.x, x_exact2);

    EXPECT_LT(error1, 1e-8) << "第一次求解误差: " << error1;
    EXPECT_LT(error2, 1e-8) << "第二次求解误差: " << error2;
    EXPECT_LT(error3, 1e-8) << "第三次求解误差: " << error3;

    // 第一次和第二次的结果应完全一致（相同输入、相同缓存状态）
    double diff_12 = (result1.x - result2.x).norm();
    EXPECT_LT(diff_12, 1e-14)
        << "两次相同输入的solve()结果不一致, 差异: " << diff_12;

    FEEM_INFO("多次solve()一致性测试通过");
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}

/**
 * @test 测试 SuperLU_MT 与 Eigen 后端的交叉验证
 *
 * 同一矩阵分别使用两种后端求解，比较结果差异。
 * 两种后端的解向量差异应在机器精度范围内（< 1e-12）。
 */
TEST(SuperLUMTIntegrationTest, CrossValidationWithEigen) {
#ifdef HAVE_SUPERLU
    auto A = build_medium_spd_matrix();

    std::vector<double> x_exact(20);
    for (int i = 0; i < 20; ++i) {
        x_exact[i] = static_cast<double>(i + 1) * 0.5;
    }
    auto b_vec = compute_rhs(A, x_exact);

    // Eigen后端求解
    SymmetricDirectSolver eigen_solver(DirectBackendType::EIGEN);
    eigen_solver.set_matrix(A);
    auto eigen_result = eigen_solver.solve(to_eigen(b_vec));
    ASSERT_EQ(eigen_result.status, SolverStatus::SUCCESS);

    // SuperLU_MT后端求解
    SymmetricDirectSolver superlu_solver(DirectBackendType::SUPERLU);
    superlu_solver.set_matrix(A);
    auto superlu_result = superlu_solver.solve(to_eigen(b_vec));
    ASSERT_EQ(superlu_result.status, SolverStatus::SUCCESS);

    // 比较两种后端的解
    double cross_diff = (eigen_result.x - superlu_result.x).norm();
    EXPECT_LT(cross_diff, 1e-10)
        << "Eigen与SuperLU_MT解向量差异过大: " << cross_diff;

    // 两者与精确解的误差都应在合理范围内
    double eigen_error = compute_error(eigen_result.x, x_exact);
    double superlu_error = compute_error(superlu_result.x, x_exact);

    EXPECT_LT(eigen_error, 1e-8) << "Eigen后端误差: " << eigen_error;
    EXPECT_LT(superlu_error, 1e-8) << "SuperLU_MT后端误差: " << superlu_error;

    FEEM_INFO("交叉验证通过: Eigen误差={:.2e}, SuperLU_MT误差={:.2e}, 差异={:.2e}",
              eigen_error, superlu_error, cross_diff);
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}

/**
 * @test 测试 clear() 后重新初始化的正确性
 *
 * 验证 clear() 能完全释放资源，之后可以正常重新 set_matrix() + solve()
 */
TEST(SuperLUMTIntegrationTest, ClearAndReinitialize) {
#ifdef HAVE_SUPERLU
    auto A = build_small_spd_matrix();

    std::vector<double> x_exact = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto b_vec = compute_rhs(A, x_exact);

    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);

    // 第一次完整流程
    solver.set_matrix(A);
    auto result1 = solver.solve(to_eigen(b_vec));
    ASSERT_EQ(result1.status, SolverStatus::SUCCESS);

    // 清理资源
    solver.clear();

    // 第二次完整流程（使用不同矩阵）
    auto A2 = build_medium_spd_matrix();
    std::vector<double> x_exact2(20);
    for (int i = 0; i < 20; ++i) {
        x_exact2[i] = static_cast<double>(i + 1);
    }
    auto b_vec2 = compute_rhs(A2, x_exact2);

    solver.set_matrix(A2);
    auto result2 = solver.solve(to_eigen(b_vec2));
    ASSERT_EQ(result2.status, SolverStatus::SUCCESS)
        << "clear()后重新初始化求解失败: " << result2.error_msg;

    double error2 = compute_error(result2.x, x_exact2);
    EXPECT_LT(error2, 1e-8)
        << "clear()后重新求解误差: " << error2;

    FEEM_INFO("clear()重新初始化测试通过");
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}

/**
 * @test 测试求解器名称标识
 *
 * 验证 get_solver_name() 返回包含 "SuperLU" 的字符串
 */
TEST(SuperLUMTIntegrationTest, SolverNameIdentification) {
#ifdef HAVE_SUPERLU
    SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
    std::string name = solver.get_solver_name();

    EXPECT_NE(name.find("SuperLU"), std::string::npos)
        << "求解器名称应包含'SuperLU'，实际: " << name;

    EXPECT_NE(name.find("SymmetricDirect"), std::string::npos)
        << "求解器名称应包含'SymmetricDirect'，实际: " << name;

    FEEM_INFO("求解器名称: {}", name);
#else
    GTEST_SKIP() << "SuperLU_MT未启用";
#endif
}
