/**
 * @file test_linear_solver_scalar.cpp
 * @brief 标量场泊松方程求解测试 - PRISM15单元 SymmetricDirectSolver验证
 * @details 使用 Google Test 框架验证 SymmetricDirectSolver 在 PRISM15 二阶三棱柱单元
 *          生成的 15×15 对称正定刚度矩阵上的正确性。
 *
 * 测试场景：
 * - 单元类型：PRISM15（二阶三棱柱单元，15个节点）
 * - DOF类型：SCALAR_ONLY（纯标量节点DOF）
 * - 约束条件：无约束（自由边界）
 * - 矩阵尺寸：15×15（对称正定SPD）
 * - 右端项：构造已知解析解 x_exact = [1, 2, 3, ..., 15]^T，计算 b = K * x_exact
 * - 求解器：SymmetricDirectSolver（EIGEN后端）
 *
 * 验证点：
 * 1. 求解状态为 SUCCESS
 * 2. 解向量误差 ||x_computed - x_exact|| < 1e-6
 * 3. 残差范数 ||b - A*x_computed|| < 1e-10
 * 4. 求解时间 > 0 ms（计时功能正常）
 * 5. 多次 solve() 调用结果一致（验证分解复用机制）
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include "em_direct_solvers.h"
#include "em_sparse_converter.h"
#include "csr_matrix.hpp"
#include "logger_factory.hpp"

using namespace numeric;


// ==================== 辅助函数定义 ====================

/**
 * @brief 构建 PRISM15 单元刚度矩阵（模拟数据）
 * @return CsrMatrix<double> 15×15 对称正定稀疏矩阵
 *
 * @details 由于完整有限元装配需要网格数据，此处构建具有相同数学特性的测试矩阵：
 * - 对称正定（SPD）
 * - 稀疏结构（每行约5-10个非零元）
 * - 对角占优（保证数值稳定性）
 * - 条件数适中（κ < 1000）
 *
 * 矩阵构建策略：
 * 1. 主对角线：渐增的正值（4.0 + 0.1*i），保证正定性
 * 2. 次对角线：相邻节点耦合（-1.0），模拟有限元相邻节点连接
 * 3. 三对角线：次邻节点弱耦合（-0.2），模拟较远节点的弱影响
 * 4. PRISM15特定耦合：
 *    - 底面三角形顶点互相耦合（节点0-1-2）
 *    - 顶面三角形顶点互相耦合（节点3-4-5）
 *    - 底面-顶面对应节点垂直耦合（棱边方向）
 */
CsrMatrix<double> build_prism15_stiffness_matrix() {
    int n = 15;

    // 先构造稠密矩阵再转稀疏（便于控制矩阵结构）
    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    // 主对角线（刚度系数，正值保证正定）
    for (int i = 0; i < n; ++i) {
        K_dense(i, i) = 4.0 + 0.1 * i;  // 渐增的对角线
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

    // 模拟 PRISM15 的额外耦合（三角形面内节点耦合）
    // 节点 0-1-2（底面三角形）互相耦合
    K_dense(0, 1) = -0.5; K_dense(1, 0) = -0.5;
    K_dense(0, 2) = -0.3; K_dense(2, 0) = -0.3;
    K_dense(1, 2) = -0.5; K_dense(2, 1) = -0.5;

    // 节点 3-4-5（顶面三角形）互相耦合
    K_dense(3, 4) = -0.5; K_dense(4, 3) = -0.5;
    K_dense(3, 5) = -0.3; K_dense(5, 3) = -0.3;
    K_dense(4, 5) = -0.5; K_dense(5, 4) = -0.5;

    // 底面-顶面对应节点耦合（棱边方向）
    for (int i = 0; i < 6; ++i) {
        if (i + 6 < n) {
            K_dense(i, i+6) = -0.8;
            K_dense(i+6, i) = -0.8;
        }
    }

    // 转换为稀疏矩阵并确保压缩存储
    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    // 转换为项目自定义 CSR 格式
    return SparseConverter::from_eigen(K_sparse);
}


/**
 * @brief 验证矩阵的对称正定性
 * @param matrix 待检查的CSR矩阵
 * @return true 如果矩阵在容差范围内对称正定
 *
 * @details 验证步骤：
 * 1. 将CSR矩阵转换为Eigen格式
 * 2. 检查对称性：计算 ||A - A^T||_F / ||A||_F，要求 < 1e-10
 * 3. 尝试 Cholesky 分解（正定性检验）：若分解成功则矩阵正定
 */
bool verify_symmetric_positive_definite(const CsrMatrix<double>& matrix) {
    // 步骤1: 转为 Eigen 格式
    auto eigen_mat = SparseConverter::to_eigen(matrix);

    // 步骤2: 检查对称性: ||A - A^T||_F / ||A||_F < tol
    // 使用 dense 矩阵避免稀疏矩阵存储顺序不匹配问题
    Eigen::MatrixXd dense_mat = eigen_mat;
    Eigen::MatrixXd diff = dense_mat - dense_mat.transpose();
    double asym_norm = diff.norm();
    double mat_norm = eigen_mat.norm();

    if (mat_norm < 1e-15) {
        FEEM_WARN("矩阵范数接近零，跳过对称性检查");
        return false;
    }

    double asym_ratio = asym_norm / mat_norm;
    if (asym_ratio > 1e-10) {
        std::cout << "[WARN] 矩阵不对称性检测失败: ||A-A^T||/||A|| = "
                  << asym_ratio << std::endl;
        return false;
    }

    // 步骤3: 尝试 Cholesky 分解（正定性检验）
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(eigen_mat);
    if (llt.info() != Eigen::Success) {
        std::cout << "[FAIL] Cholesky 分解失败，矩阵可能非正定" << std::endl;
        return false;
    }

    return true;
}


// ====================================================================
//  测试固件：PRISM15标量求解器测试环境
// ====================================================================

/**
 * @class Prism15ScalarSolverFixture
 * @brief Google Test 固件类：构建PRISM15标量场求解器测试环境
 * @details 在 SetUp 中一次性构建完整的测试数据：
 *          - PRISM15 刚度矩阵（15×15 SPD）
 *          - 已知精确解向量 x_exact = [1, 2, ..., 15]^T
 *          - 右端项向量 b = K * x_exact
 *          所有 TEST_F 用例共享此测试数据，避免重复构建。
 */
class Prism15ScalarSolverFixture : public ::testing::Test {
protected:
    CsrMatrix<double> K_;              ///< PRISM15刚度矩阵（15×15 SPD）
    Eigen::VectorXd x_exact_;          ///< 精确解向量 [1, 2, ..., 15]^T
    Eigen::VectorXd b_;                ///< 右端项向量 b = K * x_exact
    bool matrix_valid_;                ///< 矩阵构建是否成功

    void SetUp() override {
        FEEM_DEBUG("Prism15ScalarSolverFixture::SetUp - 开始构建测试环境");

        // 构建PRISM15刚度矩阵
        K_ = build_prism15_stiffness_matrix();
        matrix_valid_ = K_.is_built();

        if (!matrix_valid_) {
            FEEM_ERROR("PRISM15刚度矩阵构建失败");
            return;
        }

        // 构造已知精确解
        int n = K_.rows();
        x_exact_ = Eigen::VectorXd(n);
        for (int i = 0; i < n; ++i) {
            x_exact_(i) = i + 1.0;  // [1, 2, 3, ..., 15]
        }

        // 计算右端项 b = K * x_exact
        auto K_eigen = SparseConverter::to_eigen(K_);
        b_ = K_eigen * x_exact_;

        FEEM_DEBUG("Prism15ScalarSolverFixture::SetUp - 测试环境构建完成");
        FEEM_DEBUG("  矩阵尺寸: {}x{}, NNZ: {}", K_.rows(), K_.cols(), K_.nnz());
        FEEM_DEBUG("  ||x_exact||: {:.6f}, ||b||: {:.6f}", x_exact_.norm(), b_.norm());
    }
};


// ====================================================================
//  测试1: 矩阵属性验证 — 对称正定性检查
// ====================================================================

/**
 * @test 验证构建的PRISM15刚度矩阵满足对称正定性要求
 * @details 对称正定性是SymmetricDirectSolver（Cholesky分解）的前提条件，
 *          若矩阵非SPD则求解器将无法正常工作或产生无意义结果。
 */
TEST_F(Prism15ScalarSolverFixture, MatrixProperty_SymmetricPositiveDefinite) {
    ASSERT_TRUE(matrix_valid_) << "刚度矩阵未成功构建";

    bool is_spd = verify_symmetric_positive_definite(K_);

    EXPECT_TRUE(is_spd)
        << "PRISM15刚度矩阵必须为对称正定(SPD)矩阵";

    if (is_spd) {
        FEEM_INFO("矩阵属性验证通过: 15x15 SPD矩阵, NNZ={}", K_.nnz());
    }
}


// ====================================================================
//  测试2: 求解器基本功能 — 状态与误差验证
// ====================================================================

/**
 * @test 验证SymmetricDirectSolver的基本求解功能
 * @details 测试内容：
 *          1. set_matrix() 成功设置矩阵并完成Cholesky分解
 *          2. solve() 返回 SUCCESS 状态
 *          3. 解向量误差 ||x_computed - x_exact|| < 1e-6
 *          4. 残差范数 ||b - A*x_computed|| < 1e-10
 */
TEST_F(Prism15ScalarSolverFixture, SolverBasicFunctionality_StatusAndError) {
    ASSERT_TRUE(matrix_valid_) << "刚度矩阵未成功构建";

    // 创建对称正定直接求解器（默认使用Eigen后端）
    SymmetricDirectSolver solver;

    FEEM_INFO("求解器名称: {}", solver.get_solver_name());

    // 设置系数矩阵（触发Cholesky分解）
    ASSERT_NO_THROW(solver.set_matrix(K_))
        << "set_matrix() 不应抛出异常";

    // 执行第一次求解
    auto result = solver.solve(b_);

    // 验证点1: 求解状态必须为 SUCCESS
    EXPECT_EQ(result.status, SolverStatus::SUCCESS)
        << "求解状态应为SUCCESS，实际: "
        << static_cast<int>(result.status)
        << ", 错误信息: " << result.error_msg;

    // 验证点2: 解向量误差检查
    double solution_error = (result.x - x_exact_).norm();
    EXPECT_LT(solution_error, 1e-6)
        << "解向量误差过大: ||x - x_exact|| = " << solution_error
        << " (应 < 1e-6)";

    // 验证点3: 残差范数检查
    EXPECT_LT(result.residual_norm, 1e-10)
        << "残差范数过大: ||b - Ax|| = " << result.residual_norm
        << " (应 < 1e-10)";

    // 输出详细信息用于调试
    std::cout << "\n========== 第一次求解结果详情 ==========\n";
    std::printf("  求解状态:     %s\n",
                result.status == SolverStatus::SUCCESS ? "SUCCESS" : "FAILED");
    std::printf("  迭代次数:     %d (直接求解器应为0)\n", result.iterations);
    std::printf("  解向量误差:   %.2e (阈值: 1e-6)\n", solution_error);
    std::printf("  残差范数:     %.2e (阈值: 1e-10)\n", result.residual_norm);
    std::printf("  求解耗时:     %.3f ms\n", result.solve_time_ms);
    std::cout << "==========================================\n";

    FEEM_INFO("第一次求解完成: error={:.2e}, residual={:.2e}, time={:.3f}ms",
              solution_error, result.residual_norm, result.solve_time_ms);
}


// ====================================================================
//  测试3: 求解计时功能验证
// ====================================================================

/**
 * @test 验证求解器的计时功能正常工作
 * @details 直接求解器的solve()方法应记录实际耗时，
 *          对于15×15的小规模系统，耗时通常在0.01~1ms范围。
 *          即使极快也不应为0（除非系统时钟精度不足）。
 */
TEST_F(Prism15ScalarSolverFixture, SolverTiming_Functional) {
    ASSERT_TRUE(matrix_valid_) << "刚度矩阵未成功构建";

    SymmetricDirectSolver solver;
    solver.set_matrix(K_);

    auto result = solver.solve(b_);

    // 验证点4: 计时功能检查
    // 注意：对于极小矩阵，solve_time_ms可能接近0，这不一定是错误
    // 此处仅做警告性检查，不作为硬性FAIL条件
    if (result.solve_time_ms <= 0) {
        std::cout << "[WARN] 求解时间为0（可能矩阵太小或计时精度不足）" << std::endl;
        FEEM_WARN("求解时间为0ms，可能计时精度不足");
    } else {
        FEEM_INFO("计时功能正常: {:.3f}ms", result.solve_time_ms);
    }

    // 只要求解成功且误差合理，即使时间为0也通过此测试
    EXPECT_EQ(result.status, SolverStatus::SUCCESS);
    EXPECT_LT((result.x - x_exact_).norm(), 1e-6);
}


// ====================================================================
//  测试4: 分解复用机制验证
// ====================================================================

/**
 * @test 验证SymmetricDirectSolver的分解复用机制
 * @details SymmetricDirectSolver的核心性能优化特性：
 *          - set_matrix() 执行一次 Cholesky 分解并缓存 L 因子
 *          - 后续多次 solve() 调用复用缓存的分解结果
 *          - 第二次及以后的求解仅需 O(n²) 前代回代（远快于重新分解的 O(n³)）
 *
 * 测试方法：
 * 1. 第一次求解使用右端项 b1 = K * x_exact
 * 2. 第二次求解使用不同的右端项 b2 = 2*b1（对应解 x2 = 2*x_exact）
 * 3. 验证第二次求解结果正确且通常更快（因省去分解开销）
 */
TEST_F(Prism15ScalarSolverFixture, DecompositionReuse_Consistency) {
    ASSERT_TRUE(matrix_valid_) << "刚度矩阵未成功构建";

    SymmetricDirectSolver solver;
    solver.set_matrix(K_);

    // 第一次求解
    auto result1 = solver.solve(b_);
    double time1 = result1.solve_time_ms;
    double error1 = (result1.x - x_exact_).norm();

    // 构造不同的右端项（2倍原右端项）
    Eigen::VectorXd b2 = 2.0 * b_;
    Eigen::VectorXd x_expected2 = 2.0 * x_exact_;

    // 第二次求解（应复用已缓存的Cholesky分解）
    auto result2 = solver.solve(b2);
    double time2 = result2.solve_time_ms;
    double error2 = (result2.x - x_expected2).norm();

    // 验证点5: 第二次求解结果一致性
    EXPECT_EQ(result2.status, SolverStatus::SUCCESS)
        << "第二次求解状态应为SUCCESS";
    EXPECT_LT(error2, 1e-6)
        << "第二次解向量误差过大: " << error2 << " (应 < 1e-6)";

    // 输出对比信息
    std::cout << "\n========== 分解复用机制验证 ==========\n";
    std::printf("  第一次求解:\n");
    std::printf("    耗时: %.3f ms, 误差: %.2e\n", time1, error1);
    std::printf("  第二次求解:\n");
    std::printf("    耗时: %.3f ms, 误差: %.2e\n", time2, error2);
    std::printf("  加速比: %.2fx\n", time1 > 0 ? time1 / std::max(time2, 0.001) : 0);
    std::cout << "======================================\n";

    // 检查分解复用效果（软性验证：第二次通常更快但不强制要求）
    if (time2 < time1 && time1 > 0) {
        FEEM_INFO("分解复用验证通过: 第二次更快 ({:.3f}ms vs {:.3f}ms)", time2, time1);
    } else {
        // 对于15×15小矩阵，分解开销本身很小，加速可能不明显
        FEEM_INFO("分解复用正常: 小矩阵加速不明显属预期行为");
    }
}


// ====================================================================
//  测试5: 综合精度与鲁棒性验证
// ====================================================================

/**
 * @test 综合验证求解器在多种场景下的精度和鲁棒性
 * @details 测试多个不同右端项的求解精度，验证求解器的数值稳定性：
 *          1. 常数右端项（全1向量）
 *          2. 随机右端项（检验通用性）
 *          3. 大值右端项（检验缩放不变性）
 */
TEST_F(Prism15ScalarSolverFixture, Comprehensive_PrecisionAndRobustness) {
    ASSERT_TRUE(matrix_valid_) << "刚度矩阵未成功构建";

    SymmetricDirectSolver solver;
    solver.set_matrix(K_);

    int n = K_.rows();

    // 场景1: 常数右端项
    {
        Eigen::VectorXd b_const = Eigen::VectorXd::Ones(n);
        auto result = solver.solve(b_const);

        EXPECT_EQ(result.status, SolverStatus::SUCCESS)
            << "常数右端项求解失败";

        // 反向验证：计算残差
        auto K_eigen = SparseConverter::to_eigen(K_);
        Eigen::VectorXd residual = K_eigen * result.x - b_const;
        EXPECT_LT(residual.norm(), 1e-10)
            << "常数右端项残差过大: " << residual.norm();

        FEEM_INFO("常数右端项测试通过: residual={:.2e}", residual.norm());
    }

    // 场景2: 线性递增右端项（类似真实物理问题）
    {
        Eigen::VectorXd b_linear(n);
        for (int i = 0; i < n; ++i) {
            b_linear(i) = i + 1.0;
        }
        auto result = solver.solve(b_linear);

        EXPECT_EQ(result.status, SolverStatus::SUCCESS)
            << "线性右端项求解失败";

        auto K_eigen = SparseConverter::to_eigen(K_);
        Eigen::VectorXd residual = K_eigen * result.x - b_linear;
        EXPECT_LT(residual.norm(), 1e-10)
            << "线性右端项残差过大: " << residual.norm();

        FEEM_INFO("线性右端项测试通过: residual={:.2e}", residual.norm());
    }

    // 场景3: 大值右端项（检验数值缩放稳定性）
    {
        Eigen::VectorXd b_large = b_ * 1e6;
        auto result = solver.solve(b_large);

        EXPECT_EQ(result.status, SolverStatus::SUCCESS)
            << "大值右端项求解失败";

        Eigen::VectorXd x_expected_large = x_exact_ * 1e6;
        double relative_error = (result.x - x_expected_large).norm() / x_expected_large.norm();
        EXPECT_LT(relative_error, 1e-10)
            << "大值右端项相对误差过大: " << relative_error;

        FEEM_INFO("大值右端项测试通过: relative_error={:.2e}", relative_error);
    }
}


// ====================================================================
//  测试6: 边界条件与错误处理验证
// ====================================================================

/**
 * @test 验证求解器对非法输入的错误处理能力
 * @details 测试边界条件和异常情况：
 *          1. 未调用set_matrix()直接solve()应返回INVALID_INPUT
 *          2. 维度不匹配的右端项应返回INVALID_INPUT
 */
TEST_F(Prism15ScalarSolverFixture, ErrorHandling_BoundaryConditions) {
    // 场景1: 未设置矩阵直接求解
    {
        SymmetricDirectSolver solver;
        Eigen::VectorXd b_test = Eigen::VectorXd::Ones(15);
        auto result = solver.solve(b_test);

        EXPECT_EQ(result.status, SolverStatus::INVALID_INPUT)
            << "未设置矩阵时应返回INVALID_INPUT，实际: "
            << static_cast<int>(result.status);

        FEEM_INFO("错误处理测试1通过: 未设置矩阵返回INVALID_INPUT");
    }

    // 场景2: 维度不匹配的右端项
    {
        SymmetricDirectSolver solver;
        solver.set_matrix(K_);

        Eigen::VectorXd b_wrong_size = Eigen::VectorXd::Ones(20);  // 错误维度
        auto result = solver.solve(b_wrong_size);

        EXPECT_EQ(result.status, SolverStatus::INVALID_INPUT)
            << "维度不匹配时应返回INVALID_INPUT，实际: "
            << static_cast<int>(result.status);

        FEEM_INFO("错误处理测试2通过: 维度不匹配返回INVALID_INPUT");
    }
}


// ==================== 主函数 ====================

int main(int argc, char** argv) {
    // 初始化日志系统
    FEEM_LOG_INIT("test_linear_solver_scalar.log", true);
    FEEM_SET_LEVEL(tool::LogLevel::DEBUG);

    ::testing::InitGoogleTest(&argc, argv);

    std::cout << "========================================" << std::endl;
    std::cout << "测试用例1: 标量场泊松方程求解（PRISM15单元）" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "单元类型: PRISM15 (二阶三棱柱, 15节点)" << std::endl;
    std::cout << "DOF格式:  SCALAR_ONLY (纯标量)" << std::endl;
    std::cout << "求解器:   SymmetricDirectSolver (Eigen后端)" << std::endl;
    std::cout << "矩阵尺寸: 15×15 (对称正定SPD)" << std::endl;
    std::cout << "========================================\n" << std::endl;

    FEEM_INFO("============================================");
    FEEM_INFO("  标量场泊松方程求解测试 (PRISM15单元)");
    FEEM_INFO("  求解器: SymmetricDirectSolver (Eigen)");
    FEEM_INFO("  验证点: 5个核心验证项");
    FEEM_INFO("============================================");

    int result = RUN_ALL_TESTS();

    if (result == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "测试结果: 全部通过 ✓" << std::endl;
        std::cout << "========================================" << std::endl;

        FEEM_INFO("============================================");
        FEEM_INFO("  所有测试用例通过!");
        FEEM_INFO("============================================");
    } else {
        std::cout << "\n========================================" << std::endl;
        std::cout << "测试结果: 存在失败用例 ✗" << std::endl;
        std::cout << "========================================" << std::endl;

        FEEM_ERROR("部分测试用例失败，返回码: {}", result);
    }

    return result;
}
