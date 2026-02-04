/**
 * @file test_sparse_matrix_gtest_stage2.cpp
 * @brief 稀疏矩阵模块阶段2功能测试（Google Test版本）
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <vector>
#include <complex>
#include <cmath>
#include <memory>

#include "numeric/vector.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/coo_matrix.hpp"
#include "numeric/preconditioner.hpp"
#include "numeric/preconditioner_impl.hpp"
#include "numeric/sym_csr_matrix.hpp"
#include "numeric/sym_csr_matrix_impl.hpp"

using namespace emag;
using namespace numeric;

/**
 * @brief 测试向量类基本功能
 */
class VectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试向量
        v1_real_ = VectorReal(5);
        v1_real_[0] = 1.0;
        v1_real_[1] = 2.0;
        v1_real_[2] = 3.0;
        v1_real_[3] = 4.0;
        v1_real_[4] = 5.0;
        
        v2_real_ = VectorReal({2.0, 3.0, 4.0, 5.0, 6.0});
        
        v1_complex_ = VectorComplex(3);
        v1_complex_[0] = std::complex<double>(1.0, 2.0);
        v1_complex_[1] = std::complex<double>(3.0, 4.0);
        v1_complex_[2] = std::complex<double>(5.0, 6.0);
    }
    
    VectorReal v1_real_;
    VectorReal v2_real_;
    VectorComplex v1_complex_;
};

TEST_F(VectorTest, BasicProperties) {
    EXPECT_EQ(v1_real_.size(), 5);
    EXPECT_EQ(v1_real_.get_data_type(), VectorDataType::REAL);
    
    EXPECT_EQ(v1_complex_.size(), 3);
    EXPECT_EQ(v1_complex_.get_data_type(), VectorDataType::COMPLEX);
}

TEST_F(VectorTest, ElementAccess) {
    EXPECT_EQ(v1_real_[0], 1.0);
    EXPECT_EQ(v1_real_[4], 5.0);
    
    EXPECT_EQ(v1_complex_[0], std::complex<double>(1.0, 2.0));
    EXPECT_EQ(v1_complex_[2], std::complex<double>(5.0, 6.0));
}

TEST_F(VectorTest, VectorOperations) {
    VectorReal v3 = v1_real_ + v2_real_;
    VectorReal v4 = v1_real_ - v2_real_;
    VectorReal v5 = v1_real_ * 2.0;
    
    EXPECT_EQ(v3[0], 3.0);
    EXPECT_EQ(v4[0], -1.0);
    EXPECT_EQ(v5[0], 2.0);
}

TEST_F(VectorTest, DotProductAndNorm) {
    double dot_result = v1_real_.dot(v2_real_);
    double norm_result = v1_real_.norm();
    
    EXPECT_NEAR(dot_result, 70.0, 1e-10);
    EXPECT_NEAR(norm_result, std::sqrt(55.0), 1e-10);
}

/**
 * @brief 测试矩阵向量乘法
 */
class MatrixVectorMultiplicationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试矩阵
        coo_real_ = CooMatrixReal(3, 3);
        coo_real_.add_value(0, 0, 1.0);
        coo_real_.add_value(0, 1, 2.0);
        coo_real_.add_value(1, 1, 3.0);
        coo_real_.add_value(2, 2, 4.0);
        
        csr_real_ = CsrMatrixReal(coo_real_.rows(), coo_real_.cols());
        csr_real_.build_from_coo(coo_real_);
        
        // 创建测试向量
        x_real_ = {1.0, 2.0, 3.0};
    }
    
    CooMatrixReal coo_real_;
    CsrMatrixReal csr_real_;
    std::vector<double> x_real_;
};

TEST_F(MatrixVectorMultiplicationTest, CsrMatrixVectorMultiplication) {
    std::vector<double> y;
    csr_real_.mat_vec(x_real_, y);
    
    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 5.0, 1e-10);  // 1*1 + 2*2 = 5
    EXPECT_NEAR(y[1], 6.0, 1e-10);  // 3*2 = 6
    EXPECT_NEAR(y[2], 12.0, 1e-10); // 4*3 = 12
}

/**
 * @brief 测试矩阵操作
 */
class MatrixOperationsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试矩阵
        coo_real_ = CooMatrixReal(2, 3);
        coo_real_.add_value(0, 0, 1.0);
        coo_real_.add_value(0, 1, 2.0);
        coo_real_.add_value(1, 1, 3.0);
        coo_real_.add_value(1, 2, 4.0);
        
        csr_real_ = CsrMatrixReal(coo_real_.rows(), coo_real_.cols());
        csr_real_.build_from_coo(coo_real_);
    }
    
    CooMatrixReal coo_real_;
    CsrMatrixReal csr_real_;
};

TEST_F(MatrixOperationsTest, MatrixScaling) {
    csr_real_.scale(2.0);
    
    std::vector<double> x = {1.0, 1.0, 1.0};
    std::vector<double> y;
    csr_real_.mat_vec(x, y);
    
    EXPECT_EQ(y.size(), 2);
    EXPECT_NEAR(y[0], 6.0, 1e-10);  // 2*(1*1 + 2*1) = 6
    EXPECT_NEAR(y[1], 14.0, 1e-10); // 2*(3*1 + 4*1) = 14
}

TEST_F(MatrixOperationsTest, MatrixTranspose) {
    auto csr_trans = csr_real_.transpose();
    
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> y;
    csr_trans.mat_vec(x, y);
    
    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 1.0, 1e-10);  // 1*1 = 1
    EXPECT_NEAR(y[1], 8.0, 1e-10);  // 2*1 + 3*2 = 8
    EXPECT_NEAR(y[2], 8.0, 1e-10);  // 4*2 = 8
}

TEST_F(MatrixOperationsTest, DiagonalOperations) {
    std::vector<double> diag;
    csr_real_.get_diag(diag);
    
    EXPECT_EQ(diag.size(), 2);
    EXPECT_NEAR(diag[0], 1.0, 1e-10);
    EXPECT_NEAR(diag[1], 3.0, 1e-10);
}

/**
 * @brief 测试预条件子
 */
class PreconditionerTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试矩阵
        coo_real_ = CooMatrixReal(3, 3);
        coo_real_.add_value(0, 0, 4.0);
        coo_real_.add_value(0, 1, 1.0);
        coo_real_.add_value(1, 0, 1.0);
        coo_real_.add_value(1, 1, 3.0);
        coo_real_.add_value(1, 2, 1.0);
        coo_real_.add_value(2, 1, 1.0);
        coo_real_.add_value(2, 2, 2.0);
        
        csr_real_ = CsrMatrixReal(coo_real_.rows(), coo_real_.cols());
        csr_real_.build_from_coo(coo_real_);
        
        // 创建测试向量
        x_ = {1.0, 2.0, 3.0};
    }
    
    CooMatrixReal coo_real_;
    CsrMatrixReal csr_real_;
    std::vector<double> x_;
};

TEST_F(PreconditionerTest, JacobiPreconditioner) {
    // 使用正确的构造函数
    JacobiPreconditioner jacobi(csr_real_, 1e-10);
    
    // 创建emag::Vector类型的输入输出向量
    emag::VectorReal x_vec(x_.size());
    for (size_t i = 0; i < x_.size(); ++i) {
        x_vec[i] = x_[i];
    }
    
    emag::VectorReal y_vec;
    jacobi.apply(x_vec, y_vec);
    
    EXPECT_EQ(y_vec.size(), 3);
    EXPECT_NEAR(y_vec[0], 0.25, 1e-10);  // 1.0 / 4.0
    EXPECT_NEAR(y_vec[1], 0.6666666666666666, 1e-10);  // 2.0 / 3.0
    EXPECT_NEAR(y_vec[2], 1.5, 1e-10);  // 3.0 / 2.0
}

TEST_F(PreconditionerTest, ILU0Preconditioner) {
    // 使用正确的构造函数
    ILU0Preconditioner ilu0(csr_real_);
    
    // 创建emag::Vector类型的输入输出向量
    emag::VectorReal x_vec(x_.size());
    for (size_t i = 0; i < x_.size(); ++i) {
        x_vec[i] = x_[i];
    }
    
    emag::VectorReal y_vec;
    ilu0.apply(x_vec, y_vec);
    
    // ILU(0)应用结果应该与原始向量不同
    EXPECT_EQ(y_vec.size(), 3);
    EXPECT_NE(y_vec[0], x_vec[0]);
    EXPECT_NE(y_vec[1], x_vec[1]);
    EXPECT_NE(y_vec[2], x_vec[2]);
}

/**
 * @brief 测试对称CSR矩阵
 */
class SymmetricCsrMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建对称测试矩阵（只包含下三角元素）
        coo_sym_ = CooMatrixReal(3, 3);
        coo_sym_.add_value(0, 0, 4.0);
        coo_sym_.add_value(1, 0, 1.0);  // 对称位置 (0,1) 的值
        coo_sym_.add_value(1, 1, 3.0);
        coo_sym_.add_value(2, 0, 0.0);  // 对称位置 (0,2) 的值（为0）
        coo_sym_.add_value(2, 1, 2.0);  // 对称位置 (1,2) 的值
        coo_sym_.add_value(2, 2, 5.0);
        
        sym_csr_ = SymCsrMatrixReal(3);
        sym_csr_.build_from_coo(coo_sym_);
        
        // 创建测试向量
        x_ = {1.0, 2.0, 3.0};
    }
    
    CooMatrixReal coo_sym_;
    SymCsrMatrixReal sym_csr_;
    std::vector<double> x_;
};

TEST_F(SymmetricCsrMatrixTest, BasicProperties) {
    EXPECT_EQ(sym_csr_.rows(), 3);
    EXPECT_EQ(sym_csr_.cols(), 3);
    EXPECT_EQ(sym_csr_.nnz(), 6);  // 对称矩阵存储下三角部分，包含对角线和所有下三角元素
    EXPECT_EQ(sym_csr_.get_data_type(), MatrixDataType::REAL);
}

TEST_F(SymmetricCsrMatrixTest, MatrixVectorMultiplication) {
    std::vector<double> y;
    sym_csr_.mat_vec(x_, y);
    
    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 6.0, 1e-10);  // 4*1 + 1*2 + 0*3 = 6
    EXPECT_NEAR(y[1], 13.0, 1e-10); // 1*1 + 3*2 + 2*3 = 1 + 6 + 6 = 13
    EXPECT_NEAR(y[2], 19.0, 1e-10); // 0*1 + 2*2 + 5*3 = 0 + 4 + 15 = 19
}

TEST_F(SymmetricCsrMatrixTest, MemoryEfficiency) {
    // 对称矩阵应该比普通CSR矩阵占用更少内存
    CsrMatrixReal csr_regular(coo_sym_.rows(), coo_sym_.cols());
    csr_regular.build_from_coo(coo_sym_);
    
    // 对称矩阵的非零元数应该少于或等于普通矩阵
    EXPECT_LE(sym_csr_.nnz(), csr_regular.nnz());
}

/**
 * @brief 主函数
 */
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}