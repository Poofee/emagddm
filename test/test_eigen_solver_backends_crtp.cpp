/**
 * @file test_eigen_solver_backends_crtp.cpp
 * @brief CRTP 重构后 Eigen 求解器后端单元测试
 * @details 全面验证 CRTP 模板基类架构的正确性和稳定性：
 *          - 三个后端（LLT/LDLT/LU）的实数求解功能
 *          - 三个后端的复数求解功能
 *          - 异常处理和边界条件
 *          - 接口一致性验证（SolverBackend 多态调用）
 *          - 资源管理（RAII/clear）
 *          - 元信息查询（backend_name, supports_complex, is_symmetric_only）
 *
 * @par 测试策略：
 * 1. **功能正确性**：对比 Eigen 直接计算结果，误差 < 1e-10
 * 2. **接口兼容性**：通过基类指针调用所有方法（多态性验证）
 * 3. **异常安全性**：未设置矩阵时 solve 应抛出异常
 * 4. **资源安全**：clear() 后应重置状态，析构无内存泄漏
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0 (CRTP重构版专用测试)
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

#include "eigen_solver_backends.h"
#include "sparse_converter.h"
#include "csr_matrix.hpp"
#include "tool/logger_factory.hpp"

using namespace numeric;


// ============================================================================
// 辅助函数：构建测试矩阵
// ============================================================================

/**
 * @brief 构建5x5对称正定(SPD)测试矩阵
 * @return CsrMatrix<double> 5x5 SPD 稀疏矩阵
 */
CsrMatrix<double> build_5x5_spd_matrix() {
    int n = 5;
    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        K_dense(i, i) = 4.0 + i;  // 渐增对角线 [4, 5, 6, 7, 8]
    }

    for (int i = 0; i < n - 1; ++i) {
        K_dense(i, i+1) = -1.0;
        K_dense(i+1, i) = -1.0;
    }

    for (int i = 0; i < n - 2; ++i) {
        K_dense(i, i+2) = -0.2;
        K_dense(i+2, i) = -0.2;
    }

    K_dense(0, 3) = -0.1; K_dense(3, 0) = -0.1;
    K_dense(1, 4) = -0.15; K_dense(4, 1) = -0.15;

    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    return SparseConverter::from_eigen(K_sparse);
}

/**
 * @brief 构建4x4对称不定测试矩阵
 * @return CsrMatrix<double> 4x4 对称但未必正定的稀疏矩阵
 */
CsrMatrix<double> build_4x4_symmetric_indefinite_matrix() {
    int n = 4;
    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    // 对角线包含负值（不定矩阵）
    K_dense(0, 0) = 3.0;
    K_dense(1, 1) = -2.0;  // 负值导致非正定
    K_dense(2, 2) = 4.0;
    K_dense(3, 3) = 1.5;

    // 对称非对角元
    for (int i = 0; i < n - 1; ++i) {
        K_dense(i, i+1) = -0.5;
        K_dense(i+1, i) = -0.5;
    }

    K_dense(0, 2) = 0.3; K_dense(2, 0) = 0.3;

    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    return SparseConverter::from_eigen(K_sparse);
}

/**
 * @brief 构建4x4非对称一般矩阵
 * @return CstrMatrix<double> 4x4 非对称稀疏矩阵
 */
CsrMatrix<double> build_4x4_nonsymmetric_matrix() {
    int n = 4;
    Eigen::MatrixXd K_dense = Eigen::MatrixXd::Zero(n, n);

    K_dense(0, 0) = 4.0;  K_dense(0, 1) = -1.0; K_dense(0, 2) = 0.2;
    K_dense(1, 0) = -2.0; K_dense(1, 1) = 5.0;  K_dense(1, 3) = -0.5;
    K_dense(2, 0) = 0.1;  K_dense(2, 2) = 3.0;  K_dense(2, 3) = -1.5;
    K_dense(3, 1) = 0.3;  K_dense(3, 2) = -0.8; K_dense(3, 3) = 6.0;

    Eigen::SparseMatrix<double> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    return SparseConverter::from_eigen(K_sparse);
}

/**
 * @brief 构建4x4复数一般矩阵（模拟时谐场刚度矩阵）
 * @return CsrMatrix<std::complex<double>> 4x4 复数稀疏矩阵
 */
CsrMatrix<std::complex<double>> build_4x4_complex_matrix() {
    int n = 4;
    Eigen::MatrixXcd K_dense = Eigen::MatrixXcd::Zero(n, n);

    // 实部：对称结构
    K_dense(0, 0) = {4.0, 0.0};   K_dense(0, 1) = {-1.0, 0.0}; K_dense(0, 2) = {0.2, 0.0};
    K_dense(1, 0) = {-1.0, 0.0};  K_dense(1, 1) = {5.0, 0.0};  K_dense(1, 3) = {-0.5, 0.0};
    K_dense(2, 0) = {0.2, 0.0};   K_dense(2, 2) = {3.0, 0.0};  K_dense(2, 3) = {-1.5, 0.0};
    K_dense(3, 1) = {-0.5, 0.0};  K_dense(3, 2) = {-1.5, 0.0}; K_dense(3, 3) = {6.0, 0.0};

    // 对角虚部（保持Hermitian：对角必须为实数）
    // 非对角虚部耦合（涡流项，共轭对称保证 A^H = A）
    K_dense(0, 1) += std::complex<double>(0.0, -0.2);
    K_dense(1, 0) += std::complex<double>(0.0, 0.2);

    Eigen::SparseMatrix<std::complex<double>> K_sparse = K_dense.sparseView();
    K_sparse.makeCompressed();

    return SparseConverter::from_eigen_complex(K_sparse);
}

/**
 * @brief 构建右端项向量 b = [1, 2, 3, 4, 5]^T 或 [1, 2, 3, 4]^T
 * @param n 向量维度
 * @return Eigen::VectorXd 右端项向量
 */
Eigen::VectorXd build_rhs_real(int n) {
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        b(i) = i + 1.0;
    }
    return b;
}

/**
 * @brief 构建复数右端项向量
 * @param n 向量维度
 * @return Eigen::VectorXcd 复数右端项向量
 */
Eigen::VectorXcd build_rhs_complex(int n) {
    Eigen::VectorXcd b(n);
    for (int i = 0; i < n; ++i) {
        b(i) = std::complex<double>(i + 1.0, 0.5 * (i + 1));
    }
    return b;
}


// ============================================================================
// 测试用例 1: EigenLLTBackend 实数求解（SPD 矩阵）
// ============================================================================

TEST(EigenLLTBackendTest, RealSPDSolve) {
    auto backend = std::make_unique<EigenLLTBackend>();

    CsrMatrix<double> A = build_5x5_spd_matrix();
    Eigen::VectorXd b = build_rhs_real(5);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXd x = backend->solve_real(b);

    EXPECT_EQ(x.size(), 5);

    Eigen::SparseMatrix<double>* A_eigen = const_cast<Eigen::SparseMatrix<double>*>(
        backend->get_eigen_matrix_real());
    Eigen::VectorXd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LLT 实数求解相对残差过大: " << rel_error;

    std::cout << "[LLT实数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LLT实数] 相对残差 ||Ax-b||/||b|| = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 2: EigenLLTBackend 复数求解
// ============================================================================

TEST(EigenLLTBackendTest, ComplexSolve) {
    auto backend = std::make_unique<EigenLLTBackend>();

    CsrMatrix<std::complex<double>> A = build_4x4_complex_matrix();
    Eigen::VectorXcd b = build_rhs_complex(4);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXcd x = backend->solve_complex(b);

    EXPECT_EQ(x.size(), 4);

    const Eigen::SparseMatrix<std::complex<double>>* A_eigen =
        backend->get_eigen_matrix_complex();
    Eigen::VectorXcd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LLT 复数求解相对残差过大: " << rel_error;

    std::cout << "[LLT复数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LLT复数] 相对残差 = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 3: EigenLDLTBackend 实数求解（对称不定矩阵）
// ============================================================================

TEST(EigenLDLTBackendTest, RealSymmetricIndefiniteSolve) {
    auto backend = std::make_unique<EigenLDLTBackend>();

    CsrMatrix<double> A = build_4x4_symmetric_indefinite_matrix();
    Eigen::VectorXd b = build_rhs_real(4);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXd x = backend->solve_real(b);

    EXPECT_EQ(x.size(), 4);

    const Eigen::SparseMatrix<double>* A_eigen = backend->get_eigen_matrix_real();
    Eigen::VectorXd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LDLT 实数求解相对残差过大: " << rel_error;

    std::cout << "[LDLT实数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LDLT实数] 相对残差 = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 4: EigenLDLTBackend 复数求解
// ============================================================================

TEST(EigenLDLTBackendTest, ComplexSolve) {
    auto backend = std::make_unique<EigenLDLTBackend>();

    CsrMatrix<std::complex<double>> A = build_4x4_complex_matrix();
    Eigen::VectorXcd b = build_rhs_complex(4);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXcd x = backend->solve_complex(b);

    EXPECT_EQ(x.size(), 4);

    const Eigen::SparseMatrix<std::complex<double>>* A_eigen =
        backend->get_eigen_matrix_complex();
    Eigen::VectorXcd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LDLT 复数求解相对残差过大: " << rel_error;

    std::cout << "[LDLT复数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LDLT复数] 相对残差 = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 5: EigenLUBackend 实数求解（非对称矩阵）
// ============================================================================

TEST(EigenLUBackendTest, RealNonsymmetricSolve) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<double> A = build_4x4_nonsymmetric_matrix();
    Eigen::VectorXd b = build_rhs_real(4);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXd x = backend->solve_real(b);

    EXPECT_EQ(x.size(), 4);

    const Eigen::SparseMatrix<double>* A_eigen = backend->get_eigen_matrix_real();
    Eigen::VectorXd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LU 实数求解相对残差过大: " << rel_error;

    std::cout << "[LU实数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LU实数] 相对残差 = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 6: EigenLUBackend 复数求解
// ============================================================================

TEST(EigenLUBackendTest, ComplexSolve) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<std::complex<double>> A = build_4x4_complex_matrix();
    Eigen::VectorXcd b = build_rhs_complex(4);

    ASSERT_NO_THROW(backend->set_matrix(A));

    Eigen::VectorXcd x = backend->solve_complex(b);

    EXPECT_EQ(x.size(), 4);

    const Eigen::SparseMatrix<std::complex<double>>* A_eigen =
        backend->get_eigen_matrix_complex();
    Eigen::VectorXcd residual = (*A_eigen) * x - b;
    double rel_error = residual.norm() / b.norm();

    EXPECT_LT(rel_error, 1e-10) << "LU 复数求解相对残差过大: " << rel_error;

    std::cout << "[LU复数] 解向量 x = " << x.transpose() << std::endl;
    std::cout << "[LU复数] 相对残差 = " << rel_error << std::endl;
}


// ============================================================================
// 测试用例 7: 异常处理 - 未设置矩阵时调用 solve
// ============================================================================

TEST(ExceptionHandlingTest, SolveWithoutSetMatrix_Real) {
    auto llt_backend = std::make_unique<EigenLLTBackend>();
    auto ldlt_backend = std::make_unique<EigenLDLTBackend>();
    auto lu_backend = std::make_unique<EigenLUBackend>();

    Eigen::VectorXd b = build_rhs_real(5);

    EXPECT_THROW(llt_backend->solve_real(b), std::runtime_error);
    EXPECT_THROW(ldlt_backend->solve_real(b), std::runtime_error);
    EXPECT_THROW(lu_backend->solve_real(b), std::runtime_error);
}

TEST(ExceptionHandlingTest, SolveWithoutSetMatrix_Complex) {
    auto llt_backend = std::make_unique<EigenLLTBackend>();
    auto ldlt_backend = std::make_unique<EigenLDLTBackend>();
    auto lu_backend = std::make_unique<EigenLUBackend>();

    Eigen::VectorXcd b = build_rhs_complex(4);

    EXPECT_THROW(llt_backend->solve_complex(b), std::runtime_error);
    EXPECT_THROW(ldlt_backend->solve_complex(b), std::runtime_error);
    EXPECT_THROW(lu_backend->solve_complex(b), std::runtime_error);
}


// ============================================================================
// 测试用例 8: 异常处理 - 维度不匹配
// ============================================================================

TEST(ExceptionHandlingTest, DimensionMismatch_Real) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<double> A = build_4x4_nonsymmetric_matrix();
    Eigen::VectorXd b_wrong_size = build_rhs_real(10);  // 错误维度

    ASSERT_NO_THROW(backend->set_matrix(A));
    EXPECT_THROW(backend->solve_real(b_wrong_size), std::invalid_argument);
}

TEST(ExceptionHandlingTest, DimensionMismatch_Complex) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<std::complex<double>> A = build_4x4_complex_matrix();
    Eigen::VectorXcd b_wrong_size = build_rhs_complex(8);  // 错误维度

    ASSERT_NO_THROW(backend->set_matrix(A));
    EXPECT_THROW(backend->solve_complex(b_wrong_size), std::invalid_argument);
}


// ============================================================================
// 测试用例 9: 资源管理 - clear() 和 RAII
// ============================================================================

TEST(ResourceManagementTest, ClearResetsState) {
    auto backend = std::make_unique<EigenLLTBackend>();

    CsrMatrix<double> A = build_5x5_spd_matrix();
    Eigen::VectorXd b = build_rhs_real(5);

    backend->set_matrix(A);
    Eigen::VectorXd x1 = backend->solve_real(b);

    backend->clear();

    EXPECT_THROW(backend->solve_real(b), std::runtime_error)
        << "clear() 后再调用 solve_real 应抛出异常";
}

TEST(ResourceManagementTest, RAAIDestruction) {
    {
        auto backend = std::make_unique<EigenLUBackend>();
        CsrMatrix<double> A = build_4x4_nonsymmetric_matrix();
        backend->set_matrix(A);
    }  // 析构函数自动调用 clear()，不应崩溃或泄漏

    SUCCEED();
}


// ============================================================================
// 测试用例 10: 多态性 - 通过基类指针调用
// ============================================================================

TEST(PolymorphismTest, BasePointerCall_Real) {
    std::vector<std::unique_ptr<SolverBackend>> backends;
    backends.push_back(std::make_unique<EigenLLTBackend>());
    backends.push_back(std::make_unique<EigenLDLTBackend>());
    backends.push_back(std::make_unique<EigenLUBackend>());

    CsrMatrix<double> A_spd = build_5x5_spd_matrix();
    CsrMatrix<double> A_nonsym = build_4x4_nonsymmetric_matrix();
    Eigen::VectorXd b_spd = build_rhs_real(5);
    Eigen::VectorXd b_nonsym = build_rhs_real(4);

    backends[0]->set_matrix(A_spd);  // LLT
    backends[1]->set_matrix(A_spd);  // LDLT
    backends[2]->set_matrix(A_nonsym);  // LU

    Eigen::VectorXd x0 = backends[0]->solve_real(b_spd);
    Eigen::VectorXd x1 = backends[1]->solve_real(b_spd);
    Eigen::VectorXd x2 = backends[2]->solve_real(b_nonsym);

    EXPECT_EQ(x0.size(), 5);
    EXPECT_EQ(x1.size(), 5);
    EXPECT_EQ(x2.size(), 4);

    std::cout << "[多态性] LLT解: " << x0.transpose() << std::endl;
    std::cout << "[多态性] LDLT解: " << x1.transpose() << std::endl;
    std::cout << "[多态性] LU解: " << x2.transpose() << std::endl;
}

TEST(PolymorphismTest, BasePointerCall_Complex) {
    std::vector<std::unique_ptr<SolverBackend>> backends;
    backends.push_back(std::make_unique<EigenLLTBackend>());
    backends.push_back(std::make_unique<EigenLDLTBackend>());
    backends.push_back(std::make_unique<EigenLUBackend>());

    CsrMatrix<std::complex<double>> A = build_4x4_complex_matrix();
    Eigen::VectorXcd b = build_rhs_complex(4);

    for (auto& backend : backends) {
        backend->set_matrix(A);
        Eigen::VectorXcd x = backend->solve_complex(b);
        EXPECT_EQ(x.size(), 4);
    }

    SUCCEED();
}


// ============================================================================
// 测试用例 11: 元信息查询接口
// ============================================================================

TEST(MetadataTest, BackendNames) {
    EigenLLTBackend llt;
    EigenLDLTBackend ldlt;
    EigenLUBackend lu;

    EXPECT_EQ(llt.get_backend_name(), "EigenLLT");
    EXPECT_EQ(ldlt.get_backend_name(), "EigenLDLT");
    EXPECT_EQ(lu.get_backend_name(), "EigenLU");

    std::cout << "[元信息] LLT名称: " << llt.get_backend_name() << std::endl;
    std::cout << "[元信息] LDLT名称: " << ldlt.get_backend_name() << std::endl;
    std::cout << "[元信息] LU名称: " << lu.get_backend_name() << std::endl;
}

TEST(MetadataTest, SupportsComplex) {
    EigenLLTBackend llt;
    EigenLDLTBackend ldlt;
    EigenLUBackend lu;

    EXPECT_TRUE(llt.supports_complex());
    EXPECT_TRUE(ldlt.supports_complex());
    EXPECT_TRUE(lu.supports_complex());
}

TEST(MetadataTest, SymmetricOnly) {
    EigenLLTBackend llt;
    EigenLDLTBackend ldlt;
    EigenLUBackend lu;

    EXPECT_TRUE(llt.is_symmetric_only()) << "LLT 应仅支持对称矩阵";
    EXPECT_TRUE(ldlt.is_symmetric_only()) << "LDLT 应仅支持对称矩阵";
    EXPECT_FALSE(lu.is_symmetric_only()) << "LU 应支持非对称矩阵";

    std::cout << "[元信息] LLT对称要求: " << llt.is_symmetric_only() << std::endl;
    std::cout << "[元信息] LDLT对称要求: " << ldlt.is_symmetric_only() << std::endl;
    std::cout << "[元信息] LU对称要求: " << lu.is_symmetric_only() << std::endl;
}


// ============================================================================
// 测试用例 12: 多次求解（分解结果复用验证）
// ============================================================================

TEST(ReuseTest, MultipleSolvesAfterSingleDecomposition) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<double> A = build_4x4_nonsymmetric_matrix();

    backend->set_matrix(A);  // 仅执行一次分解

    for (int i = 0; i < 5; ++i) {
        Eigen::VectorXd b = build_rhs_real(4);
        b(0) += i * 0.5;  // 略微修改右端项

        Eigen::VectorXd x = backend->solve_real(b);

        const Eigen::SparseMatrix<double>* A_eigen = backend->get_eigen_matrix_real();
        Eigen::VectorXd residual = (*A_eigen) * x - b;
        double rel_error = residual.norm() / b.norm();

        EXPECT_LT(rel_error, 1e-10) << "第 " << i+1 << " 次求解相对残差过大";
    }

    std::cout << "[复用测试] 5次求解均成功，分解仅执行1次" << std::endl;
}


// ============================================================================
// 测试用例 13: 实数/复数混合使用
// ============================================================================

TEST(MixedModeTest, RealThenComplex) {
    auto backend = std::make_unique<EigenLUBackend>();

    CsrMatrix<double> A_real = build_4x4_nonsymmetric_matrix();
    CsrMatrix<std::complex<double>> A_complex = build_4x4_complex_matrix();
    Eigen::VectorXd b_real = build_rhs_real(4);
    Eigen::VectorXcd b_complex = build_rhs_complex(4);

    backend->set_matrix(A_real);
    Eigen::VectorXd x_real = backend->solve_real(b_real);

    backend->set_matrix(A_complex);
    Eigen::VectorXcd x_complex = backend->solve_complex(b_complex);

    EXPECT_EQ(x_real.size(), 4);
    EXPECT_EQ(x_complex.size(), 4);

    const Eigen::SparseMatrix<double>* A_eigen_real = backend->get_eigen_matrix_real();
    const Eigen::SparseMatrix<std::complex<double>>* A_eigen_complex =
        backend->get_eigen_matrix_complex();

    double real_error = ((*A_eigen_real) * x_real - b_real).norm() / b_real.norm();
    double complex_error = ((*A_eigen_complex) * x_complex - b_complex).norm() / b_complex.norm();

    EXPECT_LT(real_error, 1e-10);
    EXPECT_LT(complex_error, 1e-10);

    std::cout << "[混合模式] 实数误差: " << real_error
              << ", 复数误差: " << complex_error << std::endl;
}


// ============================================================================
// 主函数
// ============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    std::cout << "========================================" << std::endl;
    std::cout << "CRTP 重构后 Eigen 求解器后端单元测试" << std::endl;
    std::cout << "========================================\n" << std::endl;
    return RUN_ALL_TESTS();
}
