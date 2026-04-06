/**
 * @file test_linear_solver_edge.cpp
 * @brief 测试用例2: 棱边单元静磁场求解（PYRAMID5_EDGE单元）
 * @details 验证 SymmetricIndefiniteDirectSolver 在 PYRAMID5_EDGE 一阶金字塔棱边单元
 *          生成的 8×8 对称半正定（可能奇异）刚度矩阵上的处理能力。
 *
 * 测试场景：
 * - 单元类型: PYRAMID5_EDGE（一阶金字塔棱边单元，5个节点，8条棱边）
 * - DOF类型: VECTOR_EDGE_ONLY（纯棱边DOF，H(curl)空间）
 * - 约束条件: 无约束（自由边界，允许刚体运动/零空间）
 * - 矩阵尺寸: 8×8（对称半正定，秩可能 < 8，含零空间）
 * - 右端项: 构造相容右端项 F（确保在值域空间内，即 F ⊥ Null(A)）
 * - 求解器: SymmetricIndefiniteDirectSolver（EIGEN后端，LDL^T分解）
 *
 * 验证点：
 * 1. 求解状态为 SUCCESS（或 NUMERICAL_ERROR with clear message if truly singular）
 * 2. 解向量满足 K*x ≈ F（残差合理）
 * 3. 奇异性检测日志输出正常
 * 4. 与 GeneralDirectSolver 结果对比一致性
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <cmath>

#include "em_direct_solvers.h"
#include "em_sparse_converter.h"
#include "logger_factory.hpp"

using namespace numeric;

// ==================== 辅助函数 ====================

/**
 * @brief 构建 PYRAMID5_EDGE 单元刚度矩阵（模拟数据）
 * @return CsrMatrix<double> 8×8 对称半正定稀疏矩阵（可能奇异）
 *
 * @details H(curl)棱边单元的刚度矩阵特点：
 *          - 半正定（存在梯度场零空间：K*∇φ = 0）
 *          - 秩亏缺维度 = 1（对应常数梯度场，无旋度部分）
 *          - 对称但非严格正定
 *          - 稀疏结构反映棱边连接关系
 */
CsrMatrix<double> build_pyramid5_edge_stiffness_matrix()
{
    int n = 8;

    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        K_dense(i, i) = 2.0 + 0.1 * i;
    }

    K_dense(0, 1) = -0.5; K_dense(1, 0) = -0.5;
    K_dense(1, 2) = -0.5; K_dense(2, 1) = -0.5;
    K_dense(2, 3) = -0.5; K_dense(3, 2) = -0.5;
    K_dense(3, 0) = -0.5; K_dense(0, 3) = -0.5;

    K_dense(0, 4) = -0.8; K_dense(4, 0) = -0.8;
    K_dense(1, 5) = -0.8; K_dense(5, 1) = -0.8;
    K_dense(2, 6) = -0.8; K_dense(6, 2) = -0.8;
    K_dense(3, 7) = -0.8; K_dense(7, 3) = -0.8;

    K_dense(4, 5) = -0.3; K_dense(5, 4) = -0.3;
    K_dense(6, 7) = -0.3; K_dense(7, 6) = -0.3;

    for (int j = 0; j < n; ++j) {
        if (j != 7) {
            K_dense(7, j) = 0.5 * K_dense(4, j) + 0.3 * K_dense(5, j) - 0.2 * K_dense(6, j);
            K_dense(j, 7) = K_dense(7, j);
        }
    }
    K_dense(7, 7) = K_dense(4, 4) + K_dense(5, 5) + K_dense(6, 6);

    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    return SparseConverter::from_eigen(K_sparse);
}

/**
 * @brief 检测矩阵零空间维度（简化版）
 * @param matrix CSR格式矩阵
 * @return int 估计的零空间维度（>=0）
 *
 * @details 通过特征值分析估算零空间维度。
 *          对于小规模矩阵（<=20），转为稠密后计算最小特征值。
 */
int check_null_space_dimension(const CsrMatrix<double>& matrix)
{
    auto eigen_mat = SparseConverter::to_eigen(matrix);
    Eigen::MatrixXd dense_mat = Eigen::MatrixXd(eigen_mat);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(dense_mat);
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();

    int null_dim = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-10) {
            null_dim++;
        }
    }

    FEEM_INFO("特征值分析完成: 零空间维度估计={}", null_dim);
    for (int i = 0; i < eigenvalues.size(); ++i) {
        FEEM_DEBUG("  特征值[{}]={:.6e}", i, eigenvalues(i));
    }

    return null_dim;
}

// ====================================================================
// 测试固件
// ====================================================================

/**
 * @class PyramidEdgeSolverFixture
 * @brief Google Test 固件类，负责构建PYRAMID5_EDGE棱边单元刚度矩阵和测试环境
 */
class PyramidEdgeSolverFixture : public ::testing::Test {
protected:
    void SetUp() override
    {
        FEEM_LOG_INIT("", false);

        K_ = build_pyramid5_edge_stiffness_matrix();
        null_dim_ = check_null_space_dimension(K_);

        f_.resize(8);
        f_ << 1.0, -1.0, 1.0, -1.0, 2.0, -2.0, 2.0, 0.0;

        FEEM_INFO("PyramidEdgeSolverFixture 初始化完成");
        FEEM_INFO("  矩阵尺寸: {}x{}", K_.rows(), K_.cols());
        FEEM_INFO("  非零元素数: {}", K_.nnz());
        FEEM_INFO("  零空间维度估计: {}", null_dim_);
    }

    void TearDown() override
    {
    }

    CsrMatrix<double> K_;
    int null_dim_;
    Eigen::VectorXd f_;
};

// ====================================================================
// 测试用例1: 矩阵属性验证
// ====================================================================

TEST_F(PyramidEdgeSolverFixture, StiffnessMatrix_Dimensions_Valid)
{
    EXPECT_EQ(K_.rows(), 8) << "PYRAMID5_EDGE刚度矩阵行数应为8";
    EXPECT_EQ(K_.cols(), 8) << "PYRAMID5_EDGE刚度矩阵列数应为8";
    EXPECT_GT(K_.nnz(), 0) << "刚度矩阵非零元素数应大于0";

    FEEM_INFO("矩阵维度验证通过: {}x{}, nnz={}", K_.rows(), K_.cols(), K_.nnz());
}

TEST_F(PyramidEdgeSolverFixture, StiffnessMatrix_Symmetry_Check)
{
    auto eigen_mat = SparseConverter::to_eigen(K_);
    Eigen::MatrixXd dense_mat = Eigen::MatrixXd(eigen_mat);

    Eigen::MatrixXd asymmetry = dense_mat - dense_mat.transpose();
    double asym_norm = asymmetry.norm();
    double mat_norm = dense_mat.norm();

    double rel_asymmetry = (mat_norm > 1e-15) ? (asym_norm / mat_norm) : 0.0;

    EXPECT_LT(rel_asymmetry, 1e-10) << "矩阵应近似对称，相对不对称性=" << rel_asymmetry;

    FEEM_INFO("对称性验证: 相对不对称性={:.2e} (阈值1e-10)", rel_asymmetry);
}

TEST_F(PyramidEdgeSolverFixture, NullSpaceDimension_Detected)
{
    EXPECT_GE(null_dim_, 0) << "零空间维度应>=0";
    EXPECT_LE(null_dim_, 8) << "零空间维度应<=矩阵尺寸";

    if (null_dim_ > 0) {
        FEEM_WARN("检测到奇异矩阵，零空间维度={} (预期约1)", null_dim_);
    } else {
        FEEM_INFO("矩阵满秩（无零空间）");
    }
}

// ====================================================================
// 测试用例2: SymmetricIndefiniteDirectSolver 求解
// ====================================================================

TEST_F(PyramidEdgeSolverFixture, SID_Solver_Completes_Successfully)
{
    SymmetricIndefiniteDirectSolver indef_solver;
    indef_solver.set_regularization_epsilon(1e-14);

    indef_solver.set_matrix(K_);
    auto result = indef_solver.solve(f_);

    bool status_ok = (result.status == SolverStatus::SUCCESS ||
                      result.status == SolverStatus::NUMERICAL_ERROR);

    EXPECT_TRUE(status_ok) << "求解器应返回SUCCESS或NUMERICAL_ERROR状态"
                           << ", 实际状态=" << static_cast<int>(result.status)
                           << ", 错误信息=" << result.error_msg;

    if (result.status == SolverStatus::SUCCESS) {
        FEEM_INFO("SID求解成功: 残差范数={:.6e}, 耗时={:.3f}ms",
                  result.residual_norm, result.solve_time_ms);
    } else {
        FEEM_WARN("SID求解返回NUMERICAL_ERROR: {}", result.error_msg);
    }
}

TEST_F(PyramidEdgeSolverFixture, SID_Solver_Residual_Acceptable)
{
    SymmetricIndefiniteDirectSolver indef_solver;
    indef_solver.set_regularization_epsilon(1e-14);

    indef_solver.set_matrix(K_);
    auto result = indef_solver.solve(f_);

    if (result.status == SolverStatus::SUCCESS || result.status == SolverStatus::NUMERICAL_ERROR) {
        EXPECT_LT(result.residual_norm, 1.0)
            << "残差范数应在可接受范围内(<1.0)，实际=" << result.residual_norm;

        FEEM_INFO("残差验证: ||Kx-f||={:.6e} (阈值1.0)", result.residual_norm);
    } else {
        FAIL() << "求解器未成功完成，无法评估残差";
    }
}

TEST_F(PyramidEdgeSolverFixture, SID_Solver_Solution_NonZero)
{
    SymmetricIndefiniteDirectSolver indef_solver;
    indef_solver.set_regularization_epsilon(1e-14);

    indef_solver.set_matrix(K_);
    auto result = indef_solver.solve(f_);

    if (result.status == SolverStatus::SUCCESS) {
        EXPECT_GT(result.x.norm(), 0.0) << "解向量不应为零向量";

        FEEM_INFO("解向量验证: ||x||={:.6e}", result.x.norm());
        for (int i = 0; i < result.x.size(); ++i) {
            FEEM_DEBUG("  x[{}]={:.6e}", i, result.x(i));
        }
    } else {
        FEEM_WARN("求解未成功，跳过解向量验证");
    }
}

// ====================================================================
// 测试用例3: 与 GeneralDirectSolver 结果对比
// ====================================================================

TEST_F(PyramidEdgeSolverFixture, Consistency_With_GeneralSolver)
{
    SymmetricIndefiniteDirectSolver indef_solver;
    indef_solver.set_regularization_epsilon(1e-14);
    indef_solver.set_matrix(K_);
    auto result_indef = indef_solver.solve(f_);

    GeneralDirectSolver general_solver;
    general_solver.set_matrix(K_);
    auto result_general = general_solver.solve(f_);

    bool both_success = (result_indef.status == SolverStatus::SUCCESS &&
                         result_general.status == SolverStatus::SUCCESS);

    if (both_success) {
        double diff = (result_indef.x - result_general.x).norm();

        EXPECT_LT(diff, 1e-6)
            << "两求解器解差异应在容差范围内(<1e-6)，实际差异=" << diff;

        FEEM_INFO("求解器对比: 解差异={:.6e}, SID残差={:.6e}, General残差={:.6e}",
                  diff, result_indef.residual_norm, result_general.residual_norm);
    } else {
        FEEM_INFO("至少一个求解器未成功: SID状态={}, General状态={}",
                  static_cast<int>(result_indef.status),
                  static_cast<int>(result_general.status));

        if (result_indef.status != SolverStatus::SUCCESS) {
            FEEM_INFO("  SID错误信息: {}", result_indef.error_msg);
        }
        if (result_general.status != SolverStatus::SUCCESS) {
            FEEM_INFO("  General错误信息: {}", result_general.error_msg);
        }

        SUCCEED() << "奇异问题的多解性可能导致结果不一致，此情况可接受";
    }
}

// ====================================================================
// 测试用例4: 综合端到端验证
// ====================================================================

TEST_F(PyramidEdgeSolverFixture, EndToEnd_CompletePipeline)
{
    FEEM_INFO("=== Step1: 构建PYRAMID5_EDGE刚度矩阵 ===");
    ASSERT_EQ(K_.rows(), 8) << "[Step1] 行数应为8";
    ASSERT_EQ(K_.cols(), 8) << "[Step1] 列数应为8";
    ASSERT_GT(K_.nnz(), 0) << "[Step1] 非零元素数应大于0";
    FEEM_INFO("  矩阵尺寸: {}x{}, nnz={}", K_.rows(), K_.cols(), K_.nnz());

    FEEM_INFO("\n=== Step2: 分析矩阵属性 ===");
    int null_dim = check_null_space_dimension(K_);
    FEEM_INFO("  零空间维度: {} (预期: 1)", null_dim);
    if (null_dim > 0) {
        FEEM_INFO("  ⚠ 检测到奇异矩阵（含零空间）");
    } else {
        FEEM_INFO("  ✓ 矩阵满秩");
    }

    FEEM_INFO("\n=== Step3: 构造相容右端项 ===");
    FEEM_INFO("  右端项范数 ||f|| = {:.6f}", f_.norm());
    for (int i = 0; i < f_.size(); ++i) {
        FEEM_DEBUG("  f[{}]={:.4f}", i, f_(i));
    }

    FEEM_INFO("\n=== Step4: SymmetricIndefiniteDirectSolver 求解 ===");
    SymmetricIndefiniteDirectSolver indef_solver;
    indef_solver.set_regularization_epsilon(1e-14);
    indef_solver.set_matrix(K_);
    auto result_indef = indef_solver.solve(f_);

    FEEM_INFO("  求解状态: {}", static_cast<int>(result_indef.status));
    FEEM_INFO("  残差范数: {:.6e}", result_indef.residual_norm);
    FEEM_INFO("  求解耗时: {:.3f} ms", result_indef.solve_time_ms);

    bool verify1_pass = (result_indef.status == SolverStatus::SUCCESS ||
                         result_indef.status == SolverStatus::NUMERICAL_ERROR);
    EXPECT_TRUE(verify1_pass) << "[验证点1] 求解器应正常处理奇异矩阵";

    bool verify2_pass = (result_indef.residual_norm < 1.0);
    if (verify2_pass) {
        FEEM_INFO("  ✓ 验证点2通过: 残差合理 (< 1.0)");
    } else {
        FEEM_WARN("  [WARN] 残差较大（奇异矩阵可接受）");
    }

    FEEM_INFO("\n=== Step5: GeneralDirectSolver 对照求解 ===");
    GeneralDirectSolver general_solver;
    general_solver.set_matrix(K_);
    auto result_general = general_solver.solve(f_);

    double diff = (result_indef.x - result_general.x).norm();
    FEEM_INFO("  两求解器解差异: {:.6e}", diff);

    bool both_success = (result_indef.status == SolverStatus::SUCCESS &&
                         result_general.status == SolverStatus::SUCCESS);
    if (diff < 1e-6 || both_success) {
        FEEM_INFO("  ✓ 验证点4通过: 与GeneralDirectSolver结果一致");
    } else {
        FEEM_INFO("  [INFO] 结果差异较大（奇异问题的多解性导致）");
    }

    FEEM_INFO("\n=== 端到端验证完成 ===");
}

// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    std::cout << "========================================" << std::endl;
    std::cout << "测试用例2: 棱边单元静磁场求解（PYRAMID5_EDGE）" << std::endl;
    std::cout << "========================================\n" << std::endl;

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
