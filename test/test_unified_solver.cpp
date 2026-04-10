/**
 * @file test_unified_solver.cpp
 * @brief 实数/复数统一架构线性求解器完整功能验证
 * @details 使用 Google Test 框架验证统一求解器架构的所有核心功能：
 *          - 实数对称正定矩阵的 UnifiedDirectSolver (Eigen LLT) 求解
 *          - 复数一般矩阵的 UnifiedDirectSolver (Eigen LU) 求解
 *          - 实数矩阵的 CGSolver 和 BiCGSTABSolver 迭代求解
 *          - 复数矩阵的 CGSolver 和 BiCGSTABSolver 迭代求解
 *          - 扩展数据结构接口（merge_duplicates, get_eigen等）
 *
 * @par 测试覆盖范围：
 * 1. 直接求解器正确性（实数SPD + 复数一般矩阵）
 * 2. 迭代求解器收敛性（实数 + 复数）
 * 3. 迭代与直接求解器交叉验证（误差 < 1e-6）
 * 4. 数据结构扩展接口功能验证
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版 - 使用 UnifiedDirectSolver 新架构)
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <memory>

#include "em_direct_solvers.h"
#include "em_iterative_solvers.h"
#include "em_sparse_converter.h"
#include "em_solver_factory.h"
#include "unified_direct_solver.h"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "vector.hpp"
#include "tool/logger_factory.hpp"  // 日志宏（FEEM_INFO/FEEM_ERROR/FEEM_DEBUG）

using namespace numeric;
using namespace emag;


// ==================== 辅助函数定义 ====================

/**
 * @brief 构建5x5对称正定(SPD)测试矩阵
 * @return CsrMatrix<double> 5x5 对称正定稀疏矩阵
 *
 * @details 矩阵构建策略：
 * - 主对角线：渐增正值（4.0 + i），保证正定性
 * - 次对角线：-1.0，模拟相邻节点耦合
 * - 三对角线：-0.2，模拟次邻节点弱耦合
 * - 额外耦合项增强稀疏性和真实性
 */
CsrMatrix<double> build_5x5_spd_matrix() {
    int n = 5;

    // 先构造稠密矩阵再转稀疏（便于控制矩阵结构）
    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    // 主对角线（刚度系数，正值保证正定）
    for (int i = 0; i < n; ++i) {
        K_dense(i, i) = 4.0 + i;  // 渐增的对角线 [4, 5, 6, 7, 8]
    }

    // 次对角线（相邻节点耦合）
    for (int i = 0; i < n - 1; ++i) {
        K_dense(i, i+1) = -1.0;
        K_dense(i+1, i) = -1.0;  // 对称
    }

    // 三对角线（次邻节点弱耦合）
    for (int i = 0; i < n - 2; ++i) {
        K_dense(i, i+2) = -0.2;
        K_dense(i+2, i) = -0.2;  // 对称
    }

    // 额外耦合（增加稀疏结构的复杂性）
    K_dense(0, 3) = -0.1; K_dense(3, 0) = -0.1;
    K_dense(1, 4) = -0.15; K_dense(4, 1) = -0.15;

    // 转换为稀疏矩阵并确保压缩存储
    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    // 转换为项目自定义 CSR 格式
    return SparseConverter::from_eigen(K_sparse);
}


/**
 * @brief 构建4x4复数一般测试矩阵（模拟时谐场刚度矩阵）
 * @return CsrMatrix<std::complex<double>> 4x4 复数稀疏矩阵
 *
 * @details 模拟50Hz时谐场分析中的复数刚度矩阵：
 * - 实部：对称正定结构（类似静磁场刚度矩阵）
 * - 虚部：非对称涡流阻尼项（jωσ 项引入的非对称性）
 */
CsrMatrix<std::complex<double>> build_4x4_complex_matrix() {
    int n = 4;

    // 构造稠密复数矩阵（便于控制结构）
    Eigen::MatrixXcd A_dense = Eigen::MatrixXcd::Zero(n, n);

    // 实部：对称正定三对角 + 额外耦合
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = std::complex<double>(4.0 + i, 0.0);  // 正实部对角线
    }
    for (int i = 0; i < n - 1; ++i) {
        A_dense(i, i+1) = std::complex<double>(-1.0, 0.0);
        A_dense(i+1, i) = std::complex<double>(-1.0, 0.0);   // 对称实部
    }
    A_dense(0, 2) = std::complex<double>(-0.3, 0.0);
    A_dense(2, 0) = std::complex<double>(-0.3, 0.0);

    // 虚部：非对称涡流阻尼项（jωσ 引入的非对称性）
    double omega_sigma = 2.0 * 3.14159 * 50.0 * 1e6;  // ωσ ≈ 314MHz
    A_dense(0, 0) += std::complex<double>(0.0, omega_sigma * 0.1);  // 主对角虚部
    A_dense(1, 1) += std::complex<double>(0.0, omega_sigma * 0.15);
    A_dense(2, 2) += std::complex<double>(0.0, omega_sigma * 0.08);
    A_dense(3, 3) += std::complex<double>(0.0, omega_sigma * 0.12);

    // 非对称虚部耦合（模拟涡流效应的方向性）
    A_dense(0, 1) += std::complex<double>(0.0, -omega_sigma * 0.05);
    A_dense(1, 2) += std::complex<double>(0.0, -omega_sigma * 0.03);
    A_dense(2, 3) += std::complex<double>(0.0, -omega_sigma * 0.04);

    // 转换为稀疏矩阵并压缩存储
    Eigen::SparseMatrix<std::complex<double>> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    // 转换为项目自定义 CSR 格式
    return SparseConverter::from_eigen_complex(A_sparse);
}


/**
 * @brief 构建5x5非对称实数测试矩阵
 * @return CsrMatrix<double> 5x5 非对称稀疏矩阵
 */
CsrMatrix<double> build_5x5_nonsymmetric_matrix() {
    int n = 5;
    Eigen::MatrixXd A_dense = Eigen::MatrixXd::Zero(n, n);

    // 对角线
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = 3.0 + i * 0.5;
    }

    // 下三角元素（不同于上三角 → 非对称）
    A_dense(1, 0) = 2.0;
    A_dense(2, 0) = -1.0;
    A_dense(2, 1) = 3.0;
    A_dense(3, 0) = 0.5;
    A_dense(3, 1) = -0.5;
    A_dense(3, 2) = 2.5;
    A_dense(4, 1) = 1.5;
    A_dense(4, 2) = -1.0;
    A_dense(4, 3) = 3.0;

    // 上三角元素（独立于下三角）
    A_dense(0, 1) = -1.0;
    A_dense(0, 2) = 0.5;
    A_dense(1, 2) = 2.0;
    A_dense(1, 3) = -0.5;
    A_dense(2, 3) = 1.5;
    A_dense(2, 4) = -0.8;
    A_dense(3, 4) = 2.0;

    Eigen::SparseMatrix<double> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen(A_sparse);
}


// ====================================================================
//  测试固件：实数对称正定矩阵直接求解器测试环境
// ====================================================================

/**
 * @class RealSymmetricDirectSolverFixture
 * @brief Google Test 固件类：构建实数SPD矩阵直接求解器测试环境
 */
class RealSymmetricDirectSolverFixture : public ::testing::Test {
protected:
    CsrMatrix<double> A_real_;              ///< 5x5 SPD矩阵
    Eigen::VectorXd x_expected_;            ///< 已知精确解向量
    Eigen::VectorXd b_real_;                ///< 右端项 b = A * x_expected
    bool matrix_valid_;                     ///< 矩阵构建是否成功

    void SetUp() override {
        FEEM_DEBUG("RealSymmetricDirectSolverFixture::SetUp - 开始构建测试环境");

        // 构建5x5 SPD矩阵
        A_real_ = build_5x5_spd_matrix();
        matrix_valid_ = A_real_.is_built();

        if (!matrix_valid_) {
            FEEM_ERROR("SPD矩阵构建失败");
            return;
        }

        // 构造已知精确解 [1, 2, 3, 4, 5]^T
        int n = A_real_.rows();
        x_expected_ = Eigen::VectorXd(n);
        for (int i = 0; i < n; ++i) {
            x_expected_(i) = i + 1.0;
        }

        // 计算右端项 b = A * x_exact
        auto A_eigen = SparseConverter::to_eigen(A_real_);
        b_real_ = A_eigen * x_expected_;

        FEEM_DEBUG("RealSymmetricDirectSolverFixture::SetUp - 测试环境构建完成");
        FEEM_DEBUG("  矩阵尺寸: {}x{}, NNZ: {}", A_real_.rows(), A_real_.cols(), A_real_.nnz());
    }
};


// ====================================================================
//  测试1: 实数对称正定矩阵 UnifiedDirectSolver (LLT) 求解
// ====================================================================

/**
 * @test 测试实数对称正定矩阵的 UnifiedDirectSolver (Eigen LLT) 求解
 * @details 测试场景：手动构建5x5对称正定矩阵，使用新架构的统一求解器
 *
 * 验证要点：
 * - 求解状态为 SUCCESS
 * - 解向量与已知解析解误差 < 1e-6
 * - 残差范数 ||Ax - b|| < 1e-10
 */
TEST_F(RealSymmetricDirectSolverFixture, RealSymmetricDirectSolver_BasicSolve) {
    ASSERT_TRUE(matrix_valid_) << "SPD矩阵未成功构建";

    // 使用工厂创建对称正定直接求解器（Eigen LLT后端）
    auto solver = EMSolverFactory::create_solver(SolverType::EIGEN_SYMMETRIC_DIRECT);
    ASSERT_NE(solver, nullptr) << "工厂创建求解器失败";
    FEEM_INFO("求解器名称: {}", solver->get_solver_name());

    // 设置系数矩阵（触发Cholesky分解）
    ASSERT_NO_THROW(solver->set_matrix(A_real_))
        << "set_matrix() 不应抛出异常";

    // 执行求解
    auto result = solver->solve(b_real_);

    // 验证点1: 求解状态必须为 SUCCESS
    EXPECT_EQ(result.status, SolverStatus::SUCCESS)
        << "求解状态应为SUCCESS，实际: "
        << static_cast<int>(result.status)
        << ", 错误信息: " << result.error_msg;

    // 验证点2: 解向量误差检查
    double solution_error = (result.x - x_expected_).norm();
    EXPECT_LT(solution_error, 1e-6)
        << "解向量误差过大: ||x - x_expected|| = " << solution_error
        << " (应 < 1e-6)";

    // 验证点3: 残差范数检查
    EXPECT_LT(result.residual_norm, 1e-10)
        << "残差范数过大: ||b - Ax|| = " << result.residual_norm
        << " (应 < 1e-10)";

    // 输出详细信息用于调试
    std::cout << "\n========== 实数SPD直接求解结果详情 ==========\n";
    std::printf("  矩阵尺寸:     %dx%d\n", A_real_.rows(), A_real_.cols());
    std::printf("  非零元素数:   %d\n", A_real_.nnz());
    std::printf("  求解状态:     %s\n",
                result.status == SolverStatus::SUCCESS ? "SUCCESS" : "FAILED");
    std::printf("  迭代次数:     %d (直接求解器应为0)\n", result.iterations);
    std::printf("  解向量误差:   %.2e (阈值: 1e-6)\n", solution_error);
    std::printf("  残差范数:     %.2e (阈值: 1e-10)\n", result.residual_norm);
    std::printf("  求解耗时:     %.3f ms\n", result.solve_time_ms);
    std::cout << "==============================================\n";

    FEEM_INFO("实数SPD直接求解完成: error={:.2e}, residual={:.2e}, time={:.3f}ms",
              solution_error, result.residual_norm, result.solve_time_ms);
}


// ====================================================================
//  测试2: 复数一般矩阵 UnifiedDirectSolver (LU) 求解
// ====================================================================

/**
 * @test 测试复数一般矩阵的 UnifiedDirectSolver (Eigen LU) 求解
 * @details 测试场景：模拟时谐场50Hz的复数刚度矩阵，使用LU分解求解
 *
 * 验证要点：
 * - 求解状态为 SUCCESS
 * - 解向量正确（通过残差 ||Ax - b|| 验证）
 * - 残差范数 < 1e-8
 */
TEST(ComplexGeneralDirectSolver, ComplexMatrix_LUSolve) {
    // 构建4x4复数矩阵
    CsrMatrix<std::complex<double>> A_complex = build_4x4_complex_matrix();
    ASSERT_TRUE(A_complex.is_built()) << "复数矩阵构建失败";

    int n = A_complex.rows();

    // 构造已知精确解 [1+j, 2+2j, 3-j, 4]^T
    Eigen::VectorXcd x_expected(n);
    x_expected(0) = std::complex<double>(1.0, 1.0);
    x_expected(1) = std::complex<double>(2.0, 2.0);
    x_expected(2) = std::complex<double>(3.0, -1.0);
    x_expected(3) = std::complex<double>(4.0, 0.0);

    // 计算右端项 b = A * x_expected
    auto A_eigen = SparseConverter::to_eigen_complex(A_complex);
    Eigen::VectorXcd b_complex = A_eigen * x_expected;

    // 使用工厂创建通用直接求解器（LU分解）
    auto solver = EMSolverFactory::create_solver(SolverType::EIGEN_GENERAL_DIRECT);
    ASSERT_NE(solver, nullptr) << "工厂创建求解器失败";

    FEEM_INFO("复数求解器名称: {}", solver->get_solver_name());

    // 设置复数系数矩阵
    ASSERT_NO_THROW(solver->set_matrix(A_complex))
        << "set_matrix() 不应抛出异常";

    // 执行复数求解
    auto result = solver->solve(b_complex);

    // 验证点1: 求解状态必须为 SUCCESS
    EXPECT_EQ(result.status, SolverStatus::SUCCESS)
        << "复数求解状态应为SUCCESS，实际: "
        << static_cast<int>(result.status)
        << ", 错误信息: " << result.error_msg;

    // 验证点2: 通过残差验证解的正确性
    // 注意：复数版本的解存储在 result.x_complex 中
    Eigen::VectorXcd residual = A_eigen * result.x_complex - b_complex;
    double residual_norm = residual.norm();

    EXPECT_LT(residual_norm, 1e-8)
        << "复数残差范数过大: ||Ax - b|| = " << residual_norm
        << " (应 < 1e-8)";

    // 验证点3: 解向量误差检查（与期望解对比）
    double solution_error = (result.x_complex - x_expected).norm();
    EXPECT_LT(solution_error, 1e-6)
        << "复数解向量误差过大: ||x - x_expected|| = " << solution_error
        << " (应 < 1e-6)";

    // 输出详细信息
    std::cout << "\n========== 复数一般矩阵LU求解结果详情 ==========\n";
    std::printf("  矩阵尺寸:     %dx%d\n", A_complex.rows(), A_complex.cols());
    std::printf("  非零元素数:   %d\n", A_complex.nnz());
    std::printf("  求解状态:     %s\n",
                result.status == SolverStatus::SUCCESS ? "SUCCESS" : "FAILED");
    std::printf("  解向量误差:   %.2e (阈值: 1e-6)\n", solution_error);
    std::printf("  残差范数:     %.2e (阈值: 1e-8)\n", residual_norm);
    std::printf("  求解耗时:     %.3f ms\n", result.solve_time_ms);

    // 输出解向量的实部和虚部
    std::cout << "\n  解向量 x:\n";
    for (int i = 0; i < n; ++i) {
        std::printf("    x[%d] = %.6f + j*%.6f\n",
                    i, result.x_complex(i).real(), result.x_complex(i).imag());
    }
    std::cout << "==================================================\n";

    FEEM_INFO("复数LU直接求解完成: error={:.2e}, residual={:.2e}, time={:.3f}ms",
              solution_error, residual_norm, result.solve_time_ms);
}


// ====================================================================
//  测试3: 实数迭代求解器（CGSolver 和 BiCGSTABSolver）
// ====================================================================

/**
 * @test 测试实数矩阵的 CGSolver 和 BiCGSTABSolver 求解
 * @details 测试场景：
 * - 对称正定矩阵使用 CG 求解，并与 UnifiedDirectSolver (LLT) 对比
 * - 非对称矩阵使用 BiCGSTAB 求解，并与 UnifiedDirectSolver (LU) 对比
 *
 * 验证要点：
 * - 两种迭代求解器都能收敛
 * - 结果与直接求解器一致（误差 < 1e-6）
 * - 迭代次数和残差记录正确
 */
TEST(RealIterativeSolvers, CGAndBiCGSTAB_ConvergenceAndAccuracy) {

    // ========== 子测试A: CGSolver 对称正定矩阵 ==========
    {
        std::cout << "\n--- 子测试A: CGSolver 对称正定矩阵 ---\n" << std::endl;

        // 构建SPD矩阵
        CsrMatrix<double> A_spd = build_5x5_spd_matrix();
        ASSERT_TRUE(A_spd.is_built()) << "SPD矩阵构建失败";

        int n = A_spd.rows();

        // 已知解和右端项
        Eigen::VectorXd x_expected(n);
        for (int i = 0; i < n; ++i) x_expected(i) = i + 1.0;

        auto A_eigen = SparseConverter::to_eigen(A_spd);
        Eigen::VectorXd b = A_eigen * x_expected;

        // 直接求解器参考解（使用新架构）
        auto direct_solver = EMSolverFactory::create_solver(SolverType::EIGEN_SYMMETRIC_DIRECT);
        ASSERT_NE(direct_solver, nullptr) << "直接求解器创建失败";
        direct_solver->set_matrix(A_spd);
        auto direct_result = direct_solver->solve(b);
        ASSERT_EQ(direct_result.status, SolverStatus::SUCCESS) << "直接求解器失败";

        // CG迭代求解器配置
        CGConfig cg_config;
        cg_config.tolerance = 1e-10;
        cg_config.max_iterations = 1000;
        cg_config.preconditioner = CGConfig::PreconditionerType::NONE;

        CGSolver cg_solver(cg_config);
        cg_solver.set_matrix(A_spd);
        auto cg_result = cg_solver.solve(b);

        // 验证CG收敛
        EXPECT_EQ(cg_result.status, SolverStatus::SUCCESS)
            << "CG求解器未收敛，状态: " << static_cast<int>(cg_result.status);

        // 与直接求解器对比
        double cg_error = (cg_result.x - direct_result.x).norm();
        EXPECT_LT(cg_error, 1e-6)
            << "CG解与直接解误差过大: " << cg_error;

        // 输出CG结果
        std::cout << "  CG求解器结果:\n";
        std::printf("    状态:       %s\n",
                    cg_result.status == SolverStatus::SUCCESS ? "收敛成功" : "未收敛");
        std::printf("    迭代次数:   %d\n", cg_result.iterations);
        std::printf("    残差范数:   %.2e\n", cg_result.residual_norm);
        std::printf("    与直接解误差: %.2e (阈值: 1e-6)\n", cg_error);
        std::printf("    求解耗时:   %.3f ms\n", cg_result.solve_time_ms);

        // 验证收敛历史
        ConvergenceHistory history = cg_solver.get_convergence_history();
        EXPECT_GT(history.total_iterations(), 0) << "收敛历史应为空";
        std::cout << "    收敛历史记录: " << history.total_iterations() << " 次迭代\n";

        FEEM_INFO("CG求解器测试通过: iterations={}, error={:.2e}",
                  cg_result.iterations, cg_error);
    }

    // ========== 子测试B: BiCGSTABSolver 非对称矩阵 ==========
    {
        std::cout << "\n--- 子测试B: BiCGSTABSolver 非对称矩阵 ---\n" << std::endl;

        // 构建非对称矩阵
        CsrMatrix<double> A_nonsym = build_5x5_nonsymmetric_matrix();
        ASSERT_TRUE(A_nonsym.is_built()) << "非对称矩阵构建失败";

        int n = A_nonsym.rows();

        // 已知解和右端项
        Eigen::VectorXd x_expected(n);
        for (int i = 0; i < n; ++i) x_expected(i) = (i + 1) * 1.5;

        auto A_eigen = SparseConverter::to_eigen(A_nonsym);
        Eigen::VectorXd b = A_eigen * x_expected;

        // 直接求解器参考解（使用新架构）
        auto direct_solver = EMSolverFactory::create_solver(SolverType::EIGEN_GENERAL_DIRECT);
        ASSERT_NE(direct_solver, nullptr) << "直接求解器创建失败";
        direct_solver->set_matrix(A_nonsym);
        auto direct_result = direct_solver->solve(b);
        ASSERT_EQ(direct_result.status, SolverStatus::SUCCESS) << "直接求解器失败";

        // BiCGSTAB迭代求解器配置
        BiCGSTABConfig bicgstab_config;
        bicgstab_config.tolerance = 1e-10;
        bicgstab_config.max_iterations = 2000;
        bicgstab_config.preconditioner = BiCGSTABConfig::PreconditionerType::NONE;

        BiCGSTABSolver bicgstab_solver(bicgstab_config);
        bicgstab_solver.set_matrix(A_nonsym);
        auto bicgstab_result = bicgstab_solver.solve(b);

        // 验证BiCGSTAB收敛
        EXPECT_EQ(bicgstab_result.status, SolverStatus::SUCCESS)
            << "BiCGSTAB求解器未收敛，状态: " << static_cast<int>(bicgstab_result.status);

        // 与直接求解器对比
        double bicgstab_error = (bicgstab_result.x - direct_result.x).norm();
        EXPECT_LT(bicgstab_error, 1e-6)
            << "BiCGSTAB解与直接解误差过大: " << bicgstab_error;

        // 输出BiCGSTAB结果
        std::cout << "  BiCGSTAB求解器结果:\n";
        std::printf("    状态:       %s\n",
                    bicgstab_result.status == SolverStatus::SUCCESS ? "收敛成功" : "未收敛");
        std::printf("    迭代次数:   %d\n", bicgstab_result.iterations);
        std::printf("    残差范数:   %.2e\n", bicgstab_result.residual_norm);
        std::printf("    与直接解误差: %.2e (阈值: 1e-6)\n", bicgstab_error);
        std::printf("    求解耗时:   %.3f ms\n", bicgstab_result.solve_time_ms);

        // 验证收敛历史
        ConvergenceHistory history = bicgstab_solver.get_convergence_history();
        EXPECT_GT(history.total_iterations(), 0) << "收敛历史不应为空";
        std::cout << "    收敛历史记录: " << history.total_iterations() << " 次迭代\n";

        FEEM_INFO("BiCGSTAB求解器测试通过: iterations={}, error={:.2e}",
                  bicgstab_result.iterations, bicgstab_error);
    }
}


// ====================================================================
//  测试4: 复数迭代求解器（CGSolver 和 BiCGSTABSolver）
// ====================================================================

/**
 * @test 测试复数矩阵的 CGSolver 和 BiCGSTABSolver 求解
 * @details 测试场景：复数矩阵使用迭代求解器求解，并与直接求解器对比
 *
 * 验证要点：
 * - 两种迭代求解器都能收敛
 * - 结果与直接求解器一致（误差 < 1e-6）
 */
TEST(ComplexIterativeSolvers, Complex_CGAndBiCGSTAB_Accuracy) {

    // 构建复数矩阵
    CsrMatrix<std::complex<double>> A_complex = build_4x4_complex_matrix();
    ASSERT_TRUE(A_complex.is_built()) << "复数矩阵构建失败";

    int n = A_complex.rows();

    // 已知解和右端项
    Eigen::VectorXcd x_expected(n);
    x_expected(0) = std::complex<double>(1.0, 1.0);
    x_expected(1) = std::complex<double>(2.0, -0.5);
    x_expected(2) = std::complex<double>(3.0, 2.0);
    x_expected(3) = std::complex<double>(4.0, -1.0);

    auto A_eigen = SparseConverter::to_eigen_complex(A_complex);
    Eigen::VectorXcd b = A_eigen * x_expected;

    // 直接求解器参考解（使用新架构）
    auto direct_solver = EMSolverFactory::create_solver(SolverType::EIGEN_GENERAL_DIRECT);
    ASSERT_NE(direct_solver, nullptr) << "复数直接求解器创建失败";
    direct_solver->set_matrix(A_complex);
    auto direct_result = direct_solver->solve(b);
    ASSERT_EQ(direct_result.status, SolverStatus::SUCCESS) << "复数直接求解器失败";

    // ========== 子测试A: 复数CGSolver ==========
    {
        std::cout << "\n--- 子测试A: 复数CGSolver ---\n" << std::endl;

        CGConfig cg_config;
        cg_config.tolerance = 1e-10;
        cg_config.max_iterations = 1000;

        CGSolver cg_solver(cg_config);
        cg_solver.set_matrix(A_complex);
        auto cg_result = cg_solver.solve(b);

        // 验证收敛
        EXPECT_EQ(cg_result.status, SolverStatus::SUCCESS)
            << "复数CG求解器未收敛，状态: " << static_cast<int>(cg_result.status);

        // 与直接求解器对比
        double cg_error = (cg_result.x_complex - direct_result.x_complex).norm();
        EXPECT_LT(cg_error, 1e-6)
            << "复数CG解与直接解误差过大: " << cg_error;

        std::cout << "  复数CG求解器结果:\n";
        std::printf("    状态:       %s\n",
                    cg_result.status == SolverStatus::SUCCESS ? "收敛成功" : "未收敛");
        std::printf("    迭代次数:   %d\n", cg_result.iterations);
        std::printf("    残差范数:   %.2e\n", cg_result.residual_norm);
        std::printf("    与直接解误差: %.2e (阈值: 1e-6)\n", cg_error);

        FEEM_INFO("复数CG求解器测试通过: iterations={}, error={:.2e}",
                  cg_result.iterations, cg_error);
    }

    // ========== 子测试B: 复数BiCGSTABSolver ==========
    {
        std::cout << "\n--- 子测试B: 复数BiCGSTABSolver ---\n" << std::endl;

        BiCGSTABConfig bicgstab_config;
        bicgstab_config.tolerance = 1e-10;
        bicgstab_config.max_iterations = 2000;

        BiCGSTABSolver bicgstab_solver(bicgstab_config);
        bicgstab_solver.set_matrix(A_complex);
        auto bicgstab_result = bicgstab_solver.solve(b);

        // 验证收敛
        EXPECT_EQ(bicgstab_result.status, SolverStatus::SUCCESS)
            << "复数BiCGSTAB求解器未收敛，状态: " << static_cast<int>(bicgstab_result.status);

        // 与直接求解器对比
        double bicgstab_error = (bicgstab_result.x_complex - direct_result.x_complex).norm();
        EXPECT_LT(bicgstab_error, 1e-6)
            << "复数BiCGSTAB解与直接解误差过大: " << bicgstab_error;

        std::cout << "  复数BiCGSTAB求解器结果:\n";
        std::printf("    状态:       %s\n",
                    bicgstab_result.status == SolverStatus::SUCCESS ? "收敛成功" : "未收敛");
        std::printf("    迭代次数:   %d\n", bicgstab_result.iterations);
        std::printf("    残差范数:   %.2e\n", bicgstab_result.residual_norm);
        std::printf("    与直接解误差: %.2e (阈值: 1e-6)\n", bicgstab_error);

        FEEM_INFO("复数BiCGSTAB求解器测试通过: iterations={}, error={:.2e}",
                  bicgstab_result.iterations, bicgstab_error);
    }
}


// ====================================================================
//  测试5: 数据结构扩展接口测试
// ====================================================================

/**
 * @test 测试扩展的数据结构接口（merge_duplicates, get_eigen 等）
 * @details 测试内容：
 * - CooMatrix 的 merge_duplicates() 功能
 * - CooMatrix/CsrMatrix 的 get_eigen_real()/get_eigen_complex() 功能
 * - Vector 的 get_eigen_real()/get_eigen_complex() 功能
 */
TEST(ExtendedDataStructures, MergeDuplicatesAndGetEigenInterfaces) {

    // ========== 子测试A: CooMatrix merge_duplicates() ==========
    {
        std::cout << "\n--- 子测试A: CooMatrix merge_duplicates() ---\n" << std::endl;

        // 创建含重复元素的COO矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 0, 2.0);      // 重复位置 (0,0)，值应为 3.0
        coo.add_value(1, 1, 4.0);
        coo.add_value(1, 1, 5.0);      // 重复位置 (1,1)，值应为 9.0
        coo.add_value(2, 2, 7.0);
        coo.add_value(0, 1, 3.0);
        coo.add_value(1, 0, 3.0);      // 对称位置

        EXPECT_EQ(coo.nnz(), 7) << "合并前应有7个非零元";

        // 执行合并
        coo.merge_duplicates();

        EXPECT_EQ(coo.nnz(), 5) << "合并后应有5个非零元（两个重复被合并）";

        // 验证合并后的值
        const auto& eigen_mat = coo.get_eigen_real();
        EXPECT_NEAR(std::abs(eigen_mat.coeff(0, 0) - 3.0), 0.0, 1e-12)
            << "(0,0)位置合并后值应为3.0";
        EXPECT_NEAR(std::abs(eigen_mat.coeff(1, 1) - 9.0), 0.0, 1e-12)
            << "(1,1)位置合并后值应为9.0";

        std::cout << "  合并前非零元数: 7\n";
        std::cout << "  合并后非零元数: " << coo.nnz() << "\n";
        std::cout << "  (0,0) 值: " << eigen_mat.coeff(0, 0) << " (预期: 3.0)\n";
        std::cout << "  (1,1) 值: " << eigen_mat.coeff(1, 1) << " (预期: 9.0)\n";
        std::cout << "  ✓ merge_duplicates() 测试通过\n";

        FEEM_INFO("CooMatrix merge_duplicates() 测试通过");
    }

    // ========== 子测试B: CooMatrix/CsrMatrix get_eigen 接口 ==========
    {
        std::cout << "\n--- 子测试B: CooMatrix/CsrMatrix get_eigen 接口 ---\n" << std::endl;

        // 实数 COO -> Eigen
        {
            CooMatrixReal coo(2, 2);
            coo.add_value(0, 0, 10.0);
            coo.add_value(0, 1, 20.0);
            coo.add_value(1, 0, 30.0);
            coo.add_value(1, 1, 40.0);

            const auto& eigen_coo = coo.get_eigen_real();
            EXPECT_EQ(eigen_coo.rows(), 2);
            EXPECT_EQ(eigen_coo.cols(), 2);
            EXPECT_EQ(eigen_coo.nonZeros(), 4);
            EXPECT_NEAR(std::abs(eigen_coo.coeff(0, 0) - 10.0), 0.0, 1e-12);
            EXPECT_NEAR(std::abs(eigen_coo.coeff(1, 1) - 40.0), 0.0, 1e-12);
            std::cout << "  ✓ 实数 CooMatrix::get_eigen_real() 测试通过\n";
        }

        // 复数 COO -> Eigen
        {
            CooMatrixComplex coo(2, 2);
            coo.add_value(0, 0, std::complex<double>(1.0, 1.0));
            coo.add_value(1, 1, std::complex<double>(2.0, -2.0));

            const auto& eigen_coo = coo.get_eigen_complex();
            EXPECT_EQ(eigen_coo.rows(), 2);
            EXPECT_EQ(eigen_coo.cols(), 2);
            EXPECT_NEAR(std::abs(eigen_coo.coeff(0, 0) - std::complex<double>(1.0, 1.0)), 0.0, 1e-12);
            std::cout << "  ✓ 复数 CooMatrix::get_eigen_complex() 测试通过\n";
        }

        // 实数 CSR -> Eigen
        {
            CooMatrixReal coo(3, 3);
            coo.add_value(0, 0, 5.0);
            coo.add_value(1, 1, 6.0);
            coo.add_value(2, 2, 7.0);

            CsrMatrixReal csr(3, 3);
            csr.build_from_coo(coo);

            const auto& eigen_csr = csr.get_eigen_real();
            EXPECT_EQ(eigen_csr.rows(), 3);
            EXPECT_EQ(eigen_csr.cols(), 3);
            EXPECT_EQ(eigen_csr.nonZeros(), 3);
            EXPECT_NEAR(std::abs(eigen_csr.coeff(1, 1) - 6.0), 0.0, 1e-12);
            std::cout << "  ✓ 实数 CsrMatrix::get_eigen_real() 测试通过\n";
        }

        // 复数 CSR -> Eigen
        {
            CooMatrixComplex coo(2, 2);
            coo.add_value(0, 0, std::complex<double>(3.0, 4.0));
            coo.add_value(1, 1, std::complex<double>(5.0, -6.0));

            CsrMatrixComplex csr(2, 2);
            csr.build_from_coo(coo);

            const auto& eigen_csr = csr.get_eigen_complex();
            EXPECT_EQ(eigen_csr.rows(), 2);
            EXPECT_EQ(eigen_csr.cols(), 2);
            EXPECT_NEAR(std::abs(eigen_csr.coeff(0, 0) - std::complex<double>(3.0, 4.0)), 0.0, 1e-12);
            std::cout << "  ✓ 复数 CsrMatrix::get_eigen_complex() 测试通过\n";
        }

        FEEM_INFO("CooMatrix/CsrMatrix get_eigen 接口测试通过");
    }

    // ========== 子测试C: Vector get_eigen 接口 ==========
    {
        std::cout << "\n--- 子测试C: Vector get_eigen 接口 ---\n" << std::endl;

        // 实数 Vector -> Eigen
        {
            VectorReal vec({1.0, 2.0, 3.0, 4.0, 5.0});

            Eigen::VectorXd eigen_vec = vec.get_eigen_real();
            EXPECT_EQ(eigen_vec.size(), 5);
            EXPECT_NEAR(std::abs(eigen_vec[0] - 1.0), 0.0, 1e-12);
            EXPECT_NEAR(std::abs(eigen_vec[4] - 5.0), 0.0, 1e-12);
            std::cout << "  ✓ 实数 Vector::get_eigen_real() 测试通过\n";
        }

        // 复数 Vector -> Eigen
        {
            VectorComplex vec({std::complex<double>(1.0, 1.0),
                               std::complex<double>(2.0, -2.0),
                               std::complex<double>(3.0, 3.0)});

            Eigen::VectorXcd eigen_vec = vec.get_eigen_complex();
            EXPECT_EQ(eigen_vec.size(), 3);
            EXPECT_NEAR(std::abs(eigen_vec[0] - std::complex<double>(1.0, 1.0)), 0.0, 1e-12);
            EXPECT_NEAR(std::abs(eigen_vec[1] - std::complex<double>(2.0, -2.0)), 0.0, 1e-12);
            std::cout << "  ✓ 复数 Vector::get_eigen_complex() 测试通过\n";
        }

        FEEM_INFO("Vector get_eigen 接口测试通过");
    }
}


// ==================== 主函数 ====================

int main(int argc, char** argv) {
    // 初始化日志系统
    FEEM_LOG_INIT("test_unified_solver.log", true);
    FEEM_SET_LEVEL(tool::LogLevel::DEBUG);

    ::testing::InitGoogleTest(&argc, argv);

    std::cout << "================================================================" << std::endl;
    std::cout << "  实数/复数统一架构线性求解器完整功能验证测试套件" << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << "测试内容:" << std::endl;
    std::cout << "  1. 实数对称正定矩阵 UnifiedDirectSolver (Eigen LLT)" << std::endl;
    std::cout << "  2. 复数一般矩阵 UnifiedDirectSolver (Eigen LU)" << std::endl;
    std::cout << "  3. 实数迭代求解器 CGSolver + BiCGSTABSolver" << std::endl;
    std::cout << "  4. 复数迭代求解器 CGSolver + BiCGSTABSolver" << std::endl;
    std::cout << "  5. 数据结构扩展接口 (merge_duplicates, get_eigen)" << std::endl;
    std::cout << "================================================================\n" << std::endl;

    FEEM_INFO("================================================================");
    FEEM_INFO("  实数/复数统一架构线性求解器完整功能验证");
    FEEM_INFO("  测试数量: 5个核心测试组");
    FEEM_INFO("================================================================");

    int result = RUN_ALL_TESTS();

    if (result == 0) {
        std::cout << "\n================================================================" << std::endl;
        std::cout << "  测试结果: 全部通过 ✓" << std::endl;
        std::cout << "================================================================" << std::endl;

        FEEM_INFO("================================================================");
        FEEM_INFO("  所有测试用例通过! 统一架构验证完成");
        FEEM_INFO("================================================================");
    } else {
        std::cout << "\n================================================================" << std::endl;
        std::cout << "  测试结果: 存在失败用例 ✗" << std::endl;
        std::cout << "================================================================" << std::endl;

        FEEM_ERROR("部分测试用例失败，返回码: {}", result);
    }

    return result;
}
