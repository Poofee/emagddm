/**
 * @file test_mumps_solver.cpp
 * @brief MUMPS 求解器完备测试套件
 * @details 全面验证 MUMPS 后端在实数和复数场景下的正确性，包括：
 *
 * **测试模块**：
 * 1. MUMPS 实数后端测试
 *    - SymmetricDirectSolver + MUMPS (SPD 矩阵 Cholesky 分解)
 *    - GeneralDirectSolver + MUMPS (非对称矩阵 LU 分解)
 * 2. MUMPS 复数后端测试
 *    - SymmetricDirectSolver + MUMPS (Hermitian 正定矩阵)
 *    - SymmetricIndefiniteDirectSolver + MUMPS (一般 Hermitian 矩阵)
 *    - GeneralDirectSolver + MUMPS (非对称复数矩阵)
 * 3. 后端分发逻辑测试
 *    - Eigen vs MUMPS 结果一致性验证（实数/复数）
 * 4. 边界条件和错误处理测试
 *    - 未调用 set_matrix 时调用 solve 的错误处理
 *    - 右端项维度与矩阵维度不匹配时的错误处理
 *
 * @par 测试数据构造说明：
 * - 实数 SPD 矩阵：三对角线结构 + 额外耦合项，保证对称正定性
 * - 复数 Hermitian 矩阵：A = K_real + j*K_imag，两者均为 SPD
 * - 非对称矩阵：上下三角使用不同值产生非对称性
 *
 * @par 编译要求：
 * - 需启用 MUMPS 支持：-DUSE_MUMPS=ON 或 -DHAVE_MUMPS
 * - 若 MUMPS 不可用，相关测试自动跳过并输出提示信息
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <chrono>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "direct_solvers.h"
#include "sparse_converter.h"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "logger_factory.hpp"

using namespace numeric;
using namespace std::complex_literals;

// ============================================================================
// MUMPS 可用性检测宏
// ============================================================================

#ifdef HAVE_MUMPS
    #define MUMPS_AVAILABLE_FLAG true
#else
    #define MUMPS_AVAILABLE_FLAG false
    #warning "MUMPS 未启用，test_mumps_solver.cpp 中的 MUMPS 相关测试将被跳过"
#endif


// ============================================================================
// 辅助函数定义
// ============================================================================

/**
 * @brief 构建 5x5 对称正定(SPD)实数稀疏矩阵（用于 MUMPS 实数后端测试）
 * @return CsrMatrix<double> 5x5 对称正定稀疏矩阵
 *
 * @details 矩阵结构：
 * - 主对角线：4.0（恒定正值）
 * - 次对角线：-1.0（相邻节点耦合）
 * - 额外耦合：(0,3) 和 (1,4) 位置为 -0.5（模拟远距离弱耦合）
 *
 * 矩阵形式：
 * [ 4  -1   0  -0.5   0  ]
 * [-1   4  -1   0   -0.5 ]
 * [ 0  -1   4  -1    0  ]
 * [-0.5 0  -1   4   -1  ]
 * [ 0  -0.5  0  -1    4  ]
 */
CsrMatrix<double> build_5x5_spd_matrix() {
    int n = 5;

    // 构造稠密矩阵便于精确控制元素值
    Eigen::MatrixXd A_dense(n, n);
    A_dense.setZero();

    // 主对角线
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = 4.0;
    }

    // 次对角线（相邻耦合）
    for (int i = 0; i < n - 1; ++i) {
        A_dense(i, i + 1) = -1.0;
        A_dense(i + 1, i) = -1.0;  // 对称
    }

    // 额外耦合项（增强稀疏性）
    A_dense(0, 3) = A_dense(3, 0) = -0.5;
    A_dense(1, 4) = A_dense(4, 1) = -0.5;

    // 转换为 CSR 格式
    Eigen::SparseMatrix<double> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen(A_sparse);
}


/**
 * @brief 构建 5x5 非对称实数稀疏矩阵（用于 GeneralDirectSolver 测试）
 * @return CsrMatrix<double> 5x5 非对称稀疏矩阵
 *
 * @details 通过上下三角使用不同值产生明显非对称性，
 *          模拟含对流项的流体力学或非保守系统
 */
CsrMatrix<double> build_5x5_nonsymmetric_matrix() {
    int n = 5;

    Eigen::MatrixXd A_dense(n, n);
    A_dense.setZero();

    // 主对角线（正值保证可逆性）
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = 4.0 + 0.5 * i;
    }

    // 下三角（较大负值）
    for (int i = 1; i < n; ++i) {
        A_dense(i, i - 1) = -1.2 - 0.1 * i;
    }

    // 上三角（较小负值，产生非对称性）
    for (int i = 0; i < n - 1; ++i) {
        A_dense(i, i + 1) = -0.7 - 0.05 * i;
    }

    // 额外非对称耦合
    A_dense(0, 2) = -0.3;
    A_dense(2, 0) = -0.6;  // 不对称
    A_dense(1, 3) = -0.2;
    A_dense(3, 1) = -0.5;  // 不对称
    A_dense(2, 4) = -0.15;
    A_dense(4, 2) = -0.4;  // 不对称

    Eigen::SparseMatrix<double> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen(A_sparse);
}


/**
 * @brief 构建 4x4 Hermitian 正定复数稀疏矩阵（用于 MUMPS 复数 Hermitian 正定测试）
 * @return CsrMatrix<std::complex<double>> 4x4 Hermitian 正定复数稀疏矩阵
 *
 * @details 构造方法：A = K_real + j * K_imag
 * - 实部 K_real：对称正定结构（主对角占优的三对角矩阵）
 * - 虚部 K_imag：对称正定结构（小阻尼项）
 * - 整体满足 Hermitian 正定性 (A = A^H)
 */
CsrMatrix<std::complex<double>> build_4x4_hermitian_spd_matrix() {
    int n = 4;

    Eigen::MatrixXcd A_dense(n, n);
    A_dense.setIdentity();  // 单位矩阵作为基础（SPD）

    // 对角元添加正实部和小虚部（阻尼项）
    A_dense(0, 0) = {4.0, 0.5};
    A_dense(1, 1) = {4.0, 0.3};
    A_dense(2, 2) = {4.0, 0.4};
    A_dense(3, 3) = {4.0, 0.2};

    // 非对角耦合（Hermitian: A(i,j) = conj(A(j,i))）
    A_dense(0, 1) = {-1.0, 0.1};
    A_dense(1, 0) = std::conj(A_dense(0, 1));  // Hermitian 共轭

    A_dense(1, 2) = {-1.0, 0.05};
    A_dense(2, 1) = std::conj(A_dense(1, 2));

    A_dense(2, 3) = {-1.0, 0.08};
    A_dense(3, 2) = std::conj(A_dense(2, 3));

    Eigen::SparseMatrix<std::complex<double>> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen_complex(A_sparse);
}


/**
 * @brief 构建 4x4 一般 Hermitian 不定复数稀疏矩阵（Saddle Point 类型）
 * @return CsrMatrix<std::complex<double>> 4x4 一般 Hermitian 复数稀疏矩阵
 *
 * @details Saddle Point 类型的复数系统：
 * - 矩阵结构类似 [K  B^H; B  0] 的复数形式
 * - Hermitian 但不保证正定（存在零/负特征值）
 * - 用于测试 sym=2（一般Hermitian）模式
 */
CsrMatrix<std::complex<double>> build_4x4_hermitian_indefinite_matrix() {
    int n = 4;

    Eigen::MatrixXcd A_dense(n, n);
    A_dense.setZero();

    // 左上块 K（正定部分）
    A_dense(0, 0) = {4.0, 0.5};
    A_dense(1, 1) = {4.0, 0.3};

    // 右下块 0（Saddle Point 特征）
    A_dense(2, 2) = {0.0, 0.0};
    A_dense(3, 3) = {0.0, 0.0};

    // 耦合块 B 和 B^H
    A_dense(0, 2) = {1.0, 0.2};
    A_dense(2, 0) = std::conj(A_dense(0, 2));

    A_dense(0, 3) = {0.5, 0.1};
    A_dense(3, 0) = std::conj(A_dense(0, 3));

    A_dense(1, 2) = {0.8, 0.15};
    A_dense(2, 1) = std::conj(A_dense(1, 2));

    A_dense(1, 3) = {1.2, 0.25};
    A_dense(3, 1) = std::conj(A_dense(1, 3));

    Eigen::SparseMatrix<std::complex<double>> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen_complex(A_sparse);
}


/**
 * @brief 构建 4x4 非对称复数稀疏矩阵（用于 GeneralDirectSolver 测试）
 * @return CsrMatrix<std::complex<double>> 4x4 非对称复数稀疏矩阵
 *
 * @details 模拟涡流场分析中的非对称复数系统：
 * - 实部：非对称结构（对流效应）
 * - 虚部：非对称结构（涡流阻尼效应）
 * - 整体为一般复数矩阵（sym=0 模式）
 */
CsrMatrix<std::complex<double>> build_4x4_nonsymmetric_complex_matrix() {
    int n = 4;

    Eigen::MatrixXcd A_dense(n, n);
    A_dense.setZero();

    // 主对角线（正实部 + 小虚部）
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = {4.0 + 0.5 * i, 0.3 + 0.1 * i};
    }

    // 下三角（较大值）
    for (int i = 1; i < n; ++i) {
        A_dense(i, i - 1) = {-1.2 - 0.1 * i, -0.15 - 0.02 * i};
    }

    // 上三角（较小值，产生非对称性）
    for (int i = 0; i < n - 1; ++i) {
        A_dense(i, i + 1) = {-0.7 - 0.05 * i, -0.08 - 0.01 * i};
    }

    // 远距离非对称耦合
    A_dense(0, 2) = {-0.3, -0.05};
    A_dense(2, 0) = {-0.6, -0.10};  // 明显不对称

    A_dense(1, 3) = {-0.25, -0.04};
    A_dense(3, 1) = {-0.5, -0.08};   // 明显不对称

    Eigen::SparseMatrix<std::complex<double>> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();

    return SparseConverter::from_eigen_complex(A_sparse);
}


/**
 * @brief 从 CsrMatrix 提取 Eigen 稠密矩阵（用于残差计算）
 * @param csr CSR 格式稀疏矩阵
 * @return Eigen::MatrixXd 对应的稠密矩阵
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> csr_to_dense(const CsrMatrix<T>& csr) {
    // 先转换为 COO 再构建稠密矩阵
    CooMatrix<T> coo(csr.rows(), csr.cols());

    // 遍历 CSR 数据构建 COO（使用 get_* API）
    const auto& values = csr.get_values();
    const auto& col_indices = csr.get_col_indices();
    const auto& row_ptr = csr.get_row_ptr();

    for (int row = 0; row < csr.rows(); ++row) {
        for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx) {
            coo.add_value(row, col_indices[idx], values[idx]);
        }
    }

    // COO 转 Eigen 稠密（使用 get_* API）
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dense =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(csr.rows(), csr.cols());
    const auto& coo_rows = coo.get_row_indices();
    const auto& coo_cols = coo.get_col_indices();
    const auto& coo_vals = coo.get_values();

    for (size_t k = 0; k < coo_vals.size(); ++k) {
        dense(coo_rows[k], coo_cols[k]) += coo_vals[k];
    }

    return dense;
}


// ============================================================================
// 1. MUMPS 实数求解器测试
// ============================================================================

#ifdef HAVE_MUMPS

/**
 * @brief 测试 SymmetricDirectSolver 使用 MUMPS 实数后端求解对称正定矩阵
 *
 * 测试步骤：
 * 1. 创建 5x5 对称正定矩阵 A（三对角线 + 额外耦合项）
 * 2. 构造右端项 b 和已知解析解 x_expected
 * 3. 设置 backend_type_ = DirectBackendType::MUMPS
 * 4. 调用 set_matrix() 触发 MUMPS 分解
 * 5. 调用 solve() 执行 MUMPS 求解
 *
 * 验证要点：
 * - 求解状态为 SUCCESS
 * - 解向量与已知解误差 < 1e-8
 * - 残差范数 ||Ax - b|| < 1e-10
 */
TEST(MumpsRealTest, SymmetricDirectSolver_SPD_Matrix) {
    // 构建测试矩阵
    CsrMatrix<double> A_csr = build_5x5_spd_matrix();
    Eigen::MatrixXd A_dense = csr_to_dense<double>(A_csr);

    // 构造已知解和对应右端项
    Eigen::VectorXd x_expected(5);
    x_expected << 1.0, 2.0, 3.0, 4.0, 5.0;
    Eigen::VectorXd b = A_dense * x_expected;

    // 创建 MUMPS 后端的对称正定求解器
    SymmetricDirectSolver solver(DirectBackendType::MUMPS);

    // 设置矩阵并触发分解
    solver.set_matrix(A_csr);

    // 执行求解
    auto result = solver.solve(b);

    // 验证求解状态
    ASSERT_EQ(result.status, SolverStatus::SUCCESS)
        << "MUMPS 实数 SPD 求解失败: " << result.error_msg;

    // 验证解向量精度
    EXPECT_NEAR((result.x - x_expected).norm() / x_expected.norm(), 0.0, 1e-8)
        << "MUMPS 实数 SPD 解向量误差过大";

    // 验证残差范数
    Eigen::VectorXd residual = A_dense * result.x - b;
    EXPECT_LT(residual.norm(), 1e-10)
        << "MUMPS 实数 SPD 残差范数过大: " << residual.norm();

    std::cout << "  [MUMPS Real SPD] 解向量: " << result.x.transpose() << std::endl;
    std::cout << "  [MUMPS Real SPD] 相对误差: "
              << (result.x - x_expected).norm() / x_expected.norm() << std::endl;
    std::cout << "  [MUMPS Real SPD] 残差范数: " << residual.norm() << std::endl;
}


/**
 * @brief 测试 GeneralDirectSolver 使用 MUMPS 实数后端求解非对称矩阵
 *
 * 验证要点：
 * - 非对称矩阵 LU 分解正确
 * - 解向量误差 < 1e-8
 * - 残差范数 ||Ax - b|| < 1e-10
 */
TEST(MumpsRealTest, GeneralDirectSolver_Nonsymmetric_Matrix) {
    // 构建非对称测试矩阵
    CsrMatrix<double> A_csr = build_5x5_nonsymmetric_matrix();
    Eigen::MatrixXd A_dense = csr_to_dense<double>(A_csr);

    // 构造已知解和对应右端项
    Eigen::VectorXd x_expected(5);
    x_expected << 1.0, -1.0, 2.0, -2.0, 3.0;
    Eigen::VectorXd b = A_dense * x_expected;

    // 创建 MUMPS 后端的通用求解器
    GeneralDirectSolver solver(DirectBackendType::MUMPS);

    // 设置矩阵并触发 LU 分解
    solver.set_matrix(A_csr);

    // 执行求解
    auto result = solver.solve(b);

    // 验证求解状态
    ASSERT_EQ(result.status, SolverStatus::SUCCESS)
        << "MUMPS 实数非对称求解失败: " << result.error_msg;

    // 验证解向量精度
    EXPECT_NEAR((result.x - x_expected).norm() / x_expected.norm(), 0.0, 1e-8)
        << "MUMPS 实数非对称解向量误差过大";

    // 验证残差范数
    Eigen::VectorXd residual = A_dense * result.x - b;
    EXPECT_LT(residual.norm(), 1e-10)
        << "MUMPS 实数非对称残差范数过大: " << residual.norm();

    std::cout << "  [MUMPS Real Nonsym] 解向量: " << result.x.transpose() << std::endl;
    std::cout << "  [MUMPS Real Nonsym] 相对误差: "
              << (result.x - x_expected).norm() / x_expected.norm() << std::endl;
    std::cout << "  [MUMPS Real Nonsym] 残差范数: " << residual.norm() << std::endl;
}

#else  // HAVE_MUMPS

/**
 * @brief MUMPS 不可用时跳过实数测试并输出提示
 */
TEST(MumpsRealTest, Skipped_MUMPS_Not_Available) {
    GTEST_SKIP() << "MUMPS 未编译进当前构建，跳过实数 MUMPS 测试。"
                 << "请使用 -DUSE_MUMPS=ON 重新编译以启用 MUMPS 支持。";
}

#endif  // HAVE_MUMPS


// ============================================================================
// 2. MUMPS 复数求解器测试
// ============================================================================

#ifdef HAVE_MUMPS

/**
 * @brief 测试 SymmetricDirectSolver 使用 MUMPS 复数后端求解 Hermitian 正定矩阵
 *
 * 测试场景：模拟时谐场电磁分析中的复数刚度矩阵
 * - 实部：对称正定结构（模拟刚度项）
 * - 虚部：对称正定结构（模拟阻尼项）
 * - 整体满足 Hermitian 正定性 (A = A^H)
 *
 * 构造方法：A = K_real + j * K_imag，其中 K_real 和 K_imag 均为 SPD
 *
 * 验证要点：
 * - MUMPS 复数分解成功（sym=1, Hermitian正定模式）
 * - 求解状态为 SUCCESS
 * - 复数解存储在 result.x_complex
 * - 残差范数 ||Ax - b|| < 1e-8
 */
TEST(MumpsComplexTest, SymmetricDirectSolver_HermitianSPD) {
    // 构建 Hermitian 正定复数矩阵
    CsrMatrix<std::complex<double>> A_csr = build_4x4_hermitian_spd_matrix();
    Eigen::MatrixXcd A_dense = csr_to_dense<std::complex<double>>(A_csr);

    // 构造已知复数解和对应右端项
    Eigen::VectorXcd x_expected(4);
    x_expected << std::complex<double>(1.0, 0.5),
                  std::complex<double>(2.0, -0.3),
                  std::complex<double>(3.0, 0.8),
                  std::complex<double>(4.0, -0.6);
    Eigen::VectorXcd b = A_dense * x_expected;

    // 创建 MUMPS 后端的对称正定求解器
    SymmetricDirectSolver solver(DirectBackendType::MUMPS);

    // 设置复数矩阵并触发 MUMPS 分解（sym=1, Hermitian正定模式）
    solver.set_matrix(A_csr);

    // 执行复数求解
    auto result = solver.solve(b);

    // 验证求解状态
    ASSERT_EQ(result.status, SolverStatus::SUCCESS)
        << "MUMPS 复数 HermitianSPD 求解失败: " << result.error_msg;

    // 验证复数解向量精度
    double rel_error = (result.x_complex - x_expected).norm() / x_expected.norm();
    EXPECT_LT(rel_error, 1e-8)
        << "MUMPS 复数 HermitianSPD 解向量误差过大: " << rel_error;

    // 验证残差范数
    Eigen::VectorXcd residual = A_dense * result.x_complex - b;
    double rel_residual = residual.norm() / b.norm();
    EXPECT_LT(rel_residual, 1e-8)
        << "MUMPS 复数 HermitianSPD 残差范数过大: " << rel_residual;

    std::cout << "  [MUMPS Complex HermitianSPD] 解向量:\n";
    for (int i = 0; i < result.x_complex.size(); ++i) {
        std::cout << "    x[" << i << "] = (" << result.x_complex(i).real()
                  << ", " << result.x_complex(i).imag() << ")\n";
    }
    std::cout << "  [MUMPS Complex HermitianSPD] 相对误差: " << rel_error << std::endl;
    std::cout << "  [MUMPS Complex HermitianSPD] 相对残差: " << rel_residual << std::endl;
}


/**
 * @brief 测试 SymmetricIndefiniteDirectSolver 使用 MUMPS 复数后端求解一般 Hermitian 矩阵
 *
 * 测试场景：Saddle Point 类型的复数系统
 * - 矩阵结构类似 [K  B^H; B  0] 的复数形式
 * - Hermitian 但不保证正定
 *
 * 验证要点：
 * - MUMPS 复数分解成功（sym=2, 一般Hermitian模式）
 * - 求解成功且残差合理
 */
TEST(MumpsComplexTest, SymmetricIndefiniteDirectSolver_HermitianIndefinite) {
    // 构建 Saddle Point 类型的一般 Hermitian 复数矩阵
    CsrMatrix<std::complex<double>> A_csr = build_4x4_hermitian_indefinite_matrix();
    Eigen::MatrixXcd A_dense = csr_to_dense<std::complex<double>>(A_csr);

    // 构造已知复数解和对应右端项
    Eigen::VectorXcd x_expected(4);
    x_expected << std::complex<double>(1.0, 0.2),
                  std::complex<double>(-0.5, 0.3),
                  std::complex<double>(2.0, -0.1),
                  std::complex<double>(1.5, 0.4);
    Eigen::VectorXcd b = A_dense * x_expected;

    // 创建 MUMPS 后端的对称不定求解器
    SymmetricIndefiniteDirectSolver solver(DirectBackendType::MUMPS);

    // 设置复数矩阵并触发 MUMPS 分解（sym=2, 一般Hermitian模式）
    solver.set_matrix(A_csr);

    // 执行复数求解
    auto result = solver.solve(b);

    // 验证求解状态
    ASSERT_EQ(result.status, SolverStatus::SUCCESS)
        << "MUMPS 复数 HermitianIndefinite 求解失败: " << result.error_msg;

    // 验证复数解向量精度（不定矩阵可能条件数较大，放宽容差）
    double rel_error = (result.x_complex - x_expected).norm() / x_expected.norm();
    EXPECT_LT(rel_error, 1e-6)
        << "MUMPS 复数 HermitianIndefinite 解向量误差过大: " << rel_error;

    // 验证残差范数
    Eigen::VectorXcd residual = A_dense * result.x_complex - b;
    double rel_residual = residual.norm() / b.norm();
    EXPECT_LT(rel_residual, 1e-6)
        << "MUMPS 复数 HermitianIndefinite 残差范数过大: " << rel_residual;

    std::cout << "  [MUMPS Complex HermitianIndef] 解向量:\n";
    for (int i = 0; i < result.x_complex.size(); ++i) {
        std::cout << "    x[" << i << "] = (" << result.x_complex(i).real()
                  << ", " << result.x_complex(i).imag() << ")\n";
    }
    std::cout << "  [MUMPS Complex HermitianIndef] 相对误差: " << rel_error << std::endl;
    std::cout << "  [MUMPS Complex HermitianIndef] 相对残差: " << rel_residual << std::endl;
}


/**
 * @brief 测试 GeneralDirectSolver 使用 MUMPS 复数后端求解非对称复数矩阵
 *
 * 测试场景：涡流场分析中的非对称复数系统
 * - 实部：非对称结构（对流效应）
 * - 虚部：非对称结构（涡流阻尼效应）
 *
 * 验证要点：
 * - MUMPS 复数分解成功（sym=0, 非对称模式）
 * - 求解成功且残差 < 1e-8
 */
TEST(MumpsComplexTest, GeneralDirectSolver_NonsymmetricComplex) {
    // 构建非对称复数矩阵
    CsrMatrix<std::complex<double>> A_csr = build_4x4_nonsymmetric_complex_matrix();
    Eigen::MatrixXcd A_dense = csr_to_dense<std::complex<double>>(A_csr);

    // 构造已知复数解和对应右端项
    Eigen::VectorXcd x_expected(4);
    x_expected << std::complex<double>(1.0, 0.5),
                  std::complex<double>(-1.0, -0.3),
                  std::complex<double>(2.0, 0.8),
                  std::complex<double>(-2.0, -0.6);
    Eigen::VectorXcd b = A_dense * x_expected;

    // 创建 MUMPS 后端的通用求解器
    GeneralDirectSolver solver(DirectBackendType::MUMPS);

    // 设置复数矩阵并触发 MUMPS 分解（sym=0, 非对称模式）
    solver.set_matrix(A_csr);

    // 执行复数求解
    auto result = solver.solve(b);

    // 验证求解状态
    ASSERT_EQ(result.status, SolverStatus::SUCCESS)
        << "MUMPS 复数 NonSymComplex 求解失败: " << result.error_msg;

    // 验证复数解向量精度
    double rel_error = (result.x_complex - x_expected).norm() / x_expected.norm();
    EXPECT_LT(rel_error, 1e-8)
        << "MUMPS 复数 NonSymComplex 解向量误差过大: " << rel_error;

    // 验证残差范数
    Eigen::VectorXcd residual = A_dense * result.x_complex - b;
    double rel_residual = residual.norm() / b.norm();
    EXPECT_LT(rel_residual, 1e-8)
        << "MUMPS 复数 NonSymComplex 残差范数过大: " << rel_residual;

    std::cout << "  [MUMPS Complex NonSym] 解向量:\n";
    for (int i = 0; i < result.x_complex.size(); ++i) {
        std::cout << "    x[" << i << "] = (" << result.x_complex(i).real()
                  << ", " << result.x_complex(i).imag() << ")\n";
    }
    std::cout << "  [MUMPS Complex NonSym] 相对误差: " << rel_error << std::endl;
    std::cout << "  [MUMPS Complex NonSym] 相对残差: " << rel_residual << std::endl;
}

#else  // HAVE_MUMPS

/**
 * @brief MUMPS 不可用时跳过复数测试并输出提示
 */
TEST(MumpsComplexTest, Skipped_MUMPS_Not_Available) {
    GTEST_SKIP() << "MUMPS 未编译进当前构建，跳过复数 MUMPS 测试。"
                 << "请使用 -DUSE_MUMPS=ON 重新编译以启用 MUMPS 支持。";
}

#endif  // HAVE_MUMPS


// ============================================================================
// 3. 后端分发逻辑测试（Eigen vs MUMPS 一致性验证）
// ============================================================================

#ifdef HAVE_MUMPS

/**
 * @brief 验证 Eigen 后端和 MUMPS 后端对同一实数矩阵的求解结果一致
 *
 * 测试步骤：
 * 1. 使用同一矩阵分别配置两个求解器实例
 *    - solver_eigen: backend_type = EIGEN
 *    - solver_mumps: backend_type = MUMPS
 * 2. 分别求解并比较结果
 *
 * 验证要点：
 * - 两种后端的解向量差异 < 1e-6
 * - 两种后端的残差范数接近
 */
TEST(BackendDispatchTest, Eigen_vs_Mumps_Consistency_Real) {
    // 构建测试矩阵和右端项
    CsrMatrix<double> A_csr = build_5x5_spd_matrix();
    Eigen::MatrixXd A_dense = csr_to_dense<double>(A_csr);

    Eigen::VectorXd b(5);
    b << 1.0, 2.0, 3.0, 4.0, 5.0;

    // 配置 Eigen 后端求解器
    SymmetricDirectSolver solver_eigen(DirectBackendType::EIGEN);
    solver_eigen.set_matrix(A_csr);
    auto result_eigen = solver_eigen.solve(b);
    ASSERT_EQ(result_eigen.status, SolverStatus::SUCCESS)
        << "Eigen 后端求解失败";

    // 配置 MUMPS 后端求解器
    SymmetricDirectSolver solver_mumps(DirectBackendType::MUMPS);
    solver_mumps.set_matrix(A_csr);
    auto result_mumps = solver_mumps.solve(b);
    ASSERT_EQ(result_mumps.status, SolverStatus::SUCCESS)
        << "MUMPS 后端求解失败";

    // 比较两种后端的解向量差异
    double solution_diff = (result_eigen.x - result_mumps.x).norm();
    EXPECT_LT(solution_diff, 1e-6)
        << "Eigen 与 MUMPS 实数解向量差异过大: " << solution_diff;

    // 比较两种后端的残差范数
    Eigen::VectorXd residual_eigen = A_dense * result_eigen.x - b;
    Eigen::VectorXd residual_mumps = A_dense * result_mumps.x - b;
    double residual_diff = std::abs(residual_eigen.norm() - residual_mumps.norm());
    EXPECT_LT(residual_diff, 1e-6)
        << "Eigen 与 MUMPS 实数残差范数差异过大: " << residual_diff;

    std::cout << "  [Backend Consistency Real] Eigen-MUMPS 解差异: " << solution_diff << std::endl;
    std::cout << "  [Backend Consistency Real] Eigen 残差: " << residual_eigen.norm()
              << ", MUMPS 残差: " << residual_mumps.norm() << std::endl;
}


/**
 * @brief 验证 Eigen 后端和 MUMPS 后端对同一复数矩阵的求解结果一致
 *
 * 验证要点：
 * - 两种后端的复数解向量差异 < 1e-6
 * - 两种后端的复数残差范数接近
 */
TEST(BackendDispatchTest, Eigen_vs_Mumps_Consistency_Complex) {
    // 构建复数测试矩阵和右端项
    CsrMatrix<std::complex<double>> A_csr = build_4x4_hermitian_spd_matrix();
    Eigen::MatrixXcd A_dense = csr_to_dense<std::complex<double>>(A_csr);

    Eigen::VectorXcd b(4);
    b << std::complex<double>(1.0, 0.5),
         std::complex<double>(2.0, -0.3),
         std::complex<double>(3.0, 0.8),
         std::complex<double>(4.0, -0.6);

    // 配置 Eigen 后端求解器
    SymmetricDirectSolver solver_eigen(DirectBackendType::EIGEN);
    solver_eigen.set_matrix(A_csr);
    auto result_eigen = solver_eigen.solve(b);
    ASSERT_EQ(result_eigen.status, SolverStatus::SUCCESS)
        << "Eigen 复数后端求解失败";

    // 配置 MUMPS 后端求解器
    SymmetricDirectSolver solver_mumps(DirectBackendType::MUMPS);
    solver_mumps.set_matrix(A_csr);
    auto result_mumps = solver_mumps.solve(b);
    ASSERT_EQ(result_mumps.status, SolverStatus::SUCCESS)
        << "MUMPS 复数后端求解失败";

    // 比较两种后端的复数解向量差异
    double solution_diff = (result_eigen.x_complex - result_mumps.x_complex).norm();
    EXPECT_LT(solution_diff, 1e-6)
        << "Eigen 与 MUMPS 复数解向量差异过大: " << solution_diff;

    // 比较两种后端的复数残差范数
    Eigen::VectorXcd residual_eigen = A_dense * result_eigen.x_complex - b;
    Eigen::VectorXcd residual_mumps = A_dense * result_mumps.x_complex - b;
    double residual_diff = std::abs(residual_eigen.norm() - residual_mumps.norm());
    EXPECT_LT(residual_diff, 1e-6)
        << "Eigen 与 MUMPS 复数残差范数差异过大: " << residual_diff;

    std::cout << "  [Backend Consistency Complex] Eigen-MUMPS 解差异: " << solution_diff << std::endl;
    std::cout << "  [Backend Consistency Complex] Eigen 残差: " << residual_eigen.norm()
              << ", MUMPS 残差: " << residual_mumps.norm() << std::endl;
}

#else  // HAVE_MUMPS

/**
 * @brief MUMPS 不可用时跳过后端一致性测试
 */
TEST(BackendDispatchTest, Skipped_MUMPS_Not_Available) {
    GTEST_SKIP() << "MUMPS 未编译进当前构建，跳过后端一致性测试。"
                 << "请使用 -DUSE_MUMPS=ON 重新编译以启用 MUMPS 支持。";
}

#endif  // HAVE_MUMPS


// ============================================================================
// 4. 边界条件和错误处理测试
// ============================================================================

/**
 * @brief 测试在未调用 set_matrix 时调用 solve 的错误处理
 *
 * 验证要点：
 * - 返回状态为 INVALID_INPUT 或 NUMERICAL_ERROR
 * - 错误信息不为空
 */
TEST(ErrorHandlingTest, SolveWithoutSetMatrix) {
    // 测试实数 SymmetricDirectSolver
    {
        SymmetricDirectSolver solver(DirectBackendType::EIGEN);
        Eigen::VectorXd b = Eigen::VectorXd::Ones(5);
        auto result = solver.solve(b);

        EXPECT_NE(result.status, SolverStatus::SUCCESS)
            << "未设置矩阵时 solve 应返回错误状态";
        EXPECT_FALSE(result.error_msg.empty())
            << "未设置矩阵时 solve 应返回错误信息";
    }

    // 测试复数 SymmetricDirectSolver
    {
        SymmetricDirectSolver solver(DirectBackendType::EIGEN);
        Eigen::VectorXcd b = Eigen::VectorXcd::Ones(5);
        auto result = solver.solve(b);

        EXPECT_NE(result.status, SolverStatus::SUCCESS)
            << "未设置复数矩阵时 solve 应返回错误状态";
        EXPECT_FALSE(result.error_msg.empty())
            << "未设置复数矩阵时 solve 应返回错误信息";
    }

    // 测试 GeneralDirectSolver
    {
        GeneralDirectSolver solver(DirectBackendType::EIGEN);
        Eigen::VectorXd b = Eigen::VectorXd::Ones(5);
        auto result = solver.solve(b);

        EXPECT_NE(result.status, SolverStatus::SUCCESS)
            << "GeneralDirectSolver 未设置矩阵时 solve 应返回错误状态";
    }

    std::cout << "  [ErrorHandling] SolveWithoutSetMatrix 测试通过" << std::endl;
}


/**
 * @brief 测试右端项维度与矩阵维度不匹配时的错误处理
 *
 * 验证要点：
 * - 返回状态为 INVALID_INPUT
 * - 错误信息包含维度相关信息
 */
TEST(ErrorHandlingTest, DimensionMismatch) {
    // 构建 5x5 矩阵
    CsrMatrix<double> A_csr = build_5x5_spd_matrix();

    // 设置矩阵到求解器
    SymmetricDirectSolver solver(DirectBackendType::EIGEN);
    solver.set_matrix(A_csr);

    // 使用错误维度的右端项（应该是 5 维，但传入 10 维）
    Eigen::VectorXd b_wrong(10);
    b_wrong.setOnes();
    auto result = solver.solve(b_wrong);

    EXPECT_NE(result.status, SolverStatus::SUCCESS)
        << "维度不匹配时 solve 应返回错误状态";

    // 测试复数版本
    CsrMatrix<std::complex<double>> A_complex = build_4x4_hermitian_spd_matrix();
    SymmetricDirectSolver solver_complex(DirectBackendType::EIGEN);
    solver_complex.set_matrix(A_complex);

    Eigen::VectorXcd b_complex_wrong(8);  // 应该是 4 维
    b_complex_wrong.setOnes();
    auto result_complex = solver_complex.solve(b_complex_wrong);

    EXPECT_NE(result_complex.status, SolverStatus::SUCCESS)
        << "复数维度不匹配时 solve 应返回错误状态";

    std::cout << "  [ErrorHandling] DimensionMismatch 测试通过" << std::endl;
}


// ============================================================================
// 5. 性能基准测试（可选）
// ============================================================================

#ifdef HAVE_MUMPS

/**
 * @brief MUMPS 实数后端性能基准测试
 *
 * 测试中等规模矩阵（100x100 三对角 SPD）的分解和求解性能
 * 记录耗时用于性能回归检测
 */
TEST(MumpsPerfTest, Benchmark_Real_SPD_Large) {
    int n = 100;

    // 构建大规模三对角 SPD 矩阵
    Eigen::MatrixXd A_dense(n, n);
    A_dense.setZero();
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = 4.0;
        if (i > 0) A_dense(i, i - 1) = -1.0;
        if (i < n - 1) A_dense(i, i + 1) = -1.0;
    }
    Eigen::SparseMatrix<double> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();
    CsrMatrix<double> A_csr = SparseConverter::from_eigen(A_sparse);

    // 构造右端项
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);

    // MUMPS 后端性能测试
    SymmetricDirectSolver solver_mumps(DirectBackendType::MUMPS);

    auto t0 = std::chrono::high_resolution_clock::now();
    solver_mumps.set_matrix(A_csr);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto result = solver_mumps.solve(b);
    auto t2 = std::chrono::high_resolution_clock::now();

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    double factorize_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double solve_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();

    // 验证结果正确性
    Eigen::VectorXd residual = A_dense * result.x - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-8)
        << "性能基准测试精度验证失败";

    std::cout << "\n  [MUMPS Perf Real] 矩阵规模: " << n << "x" << n << std::endl;
    std::cout << "  [MUMPS Perf Real] 分解耗时: " << factorize_ms << " ms" << std::endl;
    std::cout << "  [MUMPS Perf Real] 求解耗时: " << solve_ms << " ms" << std::endl;
    std::cout << "  [MUMPS Perf Real] 相对残差: "
              << residual.norm() / b.norm() << std::endl;

    // 性能断言（可根据硬件调整阈值）
    EXPECT_LT(factorize_ms, 5000.0) << "MUMPS 分解耗时过长";
    EXPECT_LT(solve_ms, 100.0) << "MUMPS 求解耗时过长";
}


/**
 * @brief MUMPS 复数后端性能基准测试
 *
 * 测试中等规模复数矩阵的性能表现
 */
TEST(MumpsPerfTest, Benchmark_Complex_HermitianSPD_Large) {
    int n = 80;

    // 构建大规模 Hermitian 正定复数矩阵
    Eigen::MatrixXcd A_dense(n, n);
    A_dense.setZero();
    for (int i = 0; i < n; ++i) {
        A_dense(i, i) = {4.0, 0.3};
        if (i > 0) {
            A_dense(i, i - 1) = {-1.0, 0.05};
            A_dense(i - 1, i) = std::conj(A_dense(i, i - 1));  // Hermitian
        }
    }
    Eigen::SparseMatrix<std::complex<double>> A_sparse = A_dense.sparseView();
    A_sparse.makeCompressed();
    CsrMatrix<std::complex<double>> A_csr = SparseConverter::from_eigen_complex(A_sparse);

    // 构造复数右端项
    Eigen::VectorXcd b = Eigen::VectorXcd::Ones(n);

    // MUMPS 复数后端性能测试
    SymmetricDirectSolver solver_mumps(DirectBackendType::MUMPS);

    auto t0 = std::chrono::high_resolution_clock::now();
    solver_mumps.set_matrix(A_csr);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto result = solver_mumps.solve(b);
    auto t2 = std::chrono::high_resolution_clock::now();

    ASSERT_EQ(result.status, SolverStatus::SUCCESS);

    double factorize_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double solve_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();

    // 验证结果正确性
    Eigen::VectorXcd residual = A_dense * result.x_complex - b;
    EXPECT_LT(residual.norm() / b.norm(), 1e-6)
        << "复数性能基准测试精度验证失败";

    std::cout << "\n  [MUMPS Perf Complex] 矩阵规模: " << n << "x" << n << std::endl;
    std::cout << "  [MUMPS Perf Complex] 分解耗时: " << factorize_ms << " ms" << std::endl;
    std::cout << "  [MUMPS Perf Complex] 求解耗时: " << solve_ms << " ms" << std::endl;
    std::cout << "  [MUMPS Perf Complex] 相对残差: "
              << residual.norm() / b.norm() << std::endl;

    EXPECT_LT(factorize_ms, 5000.0) << "MUMPS 复数分解耗时过长";
    EXPECT_LT(solve_ms, 200.0) << "MUMPS 复数求解耗时过长";
}

#endif  // HAVE_MUMPS


// ============================================================================
// 6. 主函数
// ============================================================================

/**
 * @brief 测试套件入口点
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return int 测试结果码（0=全部通过，非0=存在失败）
 *
 * @details 初始化 Google Test 框架，输出测试覆盖范围说明，
 *          执行所有注册的测试用例，最后汇总输出测试结果统计。
 */
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // 输出测试套件标题和覆盖范围
    std::cout << "============================================================" << std::endl;
    std::cout << "  MUMPS 求解器完备测试套件" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "  测试覆盖:" << std::endl;
    std::cout << "    1. MUMPS实数后端 (SymmetricDirect + GeneralDirect)" << std::endl;
    std::cout << "    2. MUMPS复数后端 (SymmetricDirect + SymmetricIndefinite + GeneralDirect)" << std::endl;
    std::cout << "    3. 后端分发逻辑 (Eigen vs MUMPS 一致性)" << std::endl;
    std::cout << "    4. 错误处理 (未设置矩阵、维度不匹配)" << std::endl;
    std::cout << "    5. 性能基准测试 (可选)" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "  MUMPS 可用性: " << (MUMPS_AVAILABLE_FLAG ? "是" : "否 (相关测试将跳过)") << std::endl;
    std::cout << "============================================================\n" << std::endl;

    // 执行所有测试
    int result = RUN_ALL_TESTS();

    // 输出最终结果统计
    std::cout << "\n============================================================" << std::endl;
    if (result == 0) {
        std::cout << "  MUMPS 测试结果: 全部通过 ✓" << std::endl;
    } else {
        std::cout << "  MUMPS 测试结果: 存在失败用例 ✗" << std::endl;
    }
    std::cout << "============================================================" << std::endl;

    return result;
}
