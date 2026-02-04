/**
 * @file test_sparse_matrix_gtest_stage1.cpp
 * @brief 稀疏矩阵模块阶段1功能测试（Google Test版本）
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <vector>
#include <complex>
#include <cmath>

#include "numeric/sparse_base.hpp"
#include "numeric/coo_matrix.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/matrix_market_io.hpp"

using namespace numeric;

/**
 * @brief 测试COO矩阵基本功能
 */
class CooMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试矩阵
        coo_real_ = CooMatrixReal(3, 3);
        coo_real_.add_value(0, 0, 1.0);
        coo_real_.add_value(0, 1, 2.0);
        coo_real_.add_value(1, 1, 3.0);
        coo_real_.add_value(2, 2, 4.0);
        
        coo_complex_ = CooMatrixComplex(2, 2);
        coo_complex_.add_value(0, 0, std::complex<double>(1.0, 2.0));
        coo_complex_.add_value(1, 0, std::complex<double>(3.0, 4.0));
        coo_complex_.add_value(1, 1, std::complex<double>(5.0, 6.0));
    }
    
    CooMatrixReal coo_real_;
    CooMatrixComplex coo_complex_;
};

TEST_F(CooMatrixTest, BasicProperties) {
    EXPECT_EQ(coo_real_.rows(), 3);
    EXPECT_EQ(coo_real_.cols(), 3);
    EXPECT_EQ(coo_real_.nnz(), 4);
    EXPECT_EQ(coo_real_.get_data_type(), MatrixDataType::REAL);
    
    EXPECT_EQ(coo_complex_.rows(), 2);
    EXPECT_EQ(coo_complex_.cols(), 2);
    EXPECT_EQ(coo_complex_.nnz(), 3);
    EXPECT_EQ(coo_complex_.get_data_type(), MatrixDataType::COMPLEX);
}

TEST_F(CooMatrixTest, ElementAccess) {
    const auto& rows = coo_real_.get_row_indices();
    const auto& cols = coo_real_.get_col_indices();
    const auto& values = coo_real_.get_values();
    
    EXPECT_EQ(rows.size(), 4);
    EXPECT_EQ(cols.size(), 4);
    EXPECT_EQ(values.size(), 4);
    
    // 验证数据正确性（不依赖顺序）
    std::map<std::pair<int, int>, double> expected_data = {
        {{0, 0}, 1.0},
        {{0, 1}, 2.0},
        {{1, 1}, 3.0},
        {{2, 2}, 4.0}
    };
    
    for (int i = 0; i < coo_real_.nnz(); ++i) {
        auto key = std::make_pair(rows[i], cols[i]);
        auto it = expected_data.find(key);
        EXPECT_NE(it, expected_data.end()) << "位置(" << rows[i] << "," << cols[i] << ")的数据不存在";
        if (it != expected_data.end()) {
            EXPECT_NEAR(values[i], it->second, 1e-10) << "位置(" << rows[i] << "," << cols[i] << ")的值不正确";
        }
    }
}

TEST_F(CooMatrixTest, ClearFunction) {
    coo_real_.clear();
    EXPECT_EQ(coo_real_.nnz(), 0);
    // clear()方法只清空非零元素，不改变矩阵尺寸
    EXPECT_EQ(coo_real_.rows(), 3);
    EXPECT_EQ(coo_real_.cols(), 3);
}

/**
 * @brief 测试CSR矩阵基本功能
 */
class CsrMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建COO矩阵并转换为CSR
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 1, 2.0);
        coo.add_value(1, 1, 3.0);
        coo.add_value(2, 2, 4.0);
        
        // 创建正确尺寸的CSR矩阵
        csr_real_ = CsrMatrixReal(coo.rows(), coo.cols());
        csr_real_.build_from_coo(coo);
        
        // 复数矩阵
        CooMatrixComplex coo_complex(2, 2);
        coo_complex.add_value(0, 0, std::complex<double>(1.0, 2.0));
        coo_complex.add_value(1, 0, std::complex<double>(3.0, 4.0));
        coo_complex.add_value(1, 1, std::complex<double>(5.0, 6.0));
        
        csr_complex_ = CsrMatrixComplex(coo_complex.rows(), coo_complex.cols());
        csr_complex_.build_from_coo(coo_complex);
    }
    
    CsrMatrixReal csr_real_;
    CsrMatrixComplex csr_complex_;
};

TEST_F(CsrMatrixTest, BasicProperties) {
    EXPECT_TRUE(csr_real_.is_built());
    EXPECT_EQ(csr_real_.rows(), 3);
    EXPECT_EQ(csr_real_.cols(), 3);
    EXPECT_EQ(csr_real_.nnz(), 4);
    EXPECT_EQ(csr_real_.get_data_type(), MatrixDataType::REAL);
    
    EXPECT_TRUE(csr_complex_.is_built());
    EXPECT_EQ(csr_complex_.rows(), 2);
    EXPECT_EQ(csr_complex_.cols(), 2);
    EXPECT_EQ(csr_complex_.nnz(), 3);
    EXPECT_EQ(csr_complex_.get_data_type(), MatrixDataType::COMPLEX);
}

TEST_F(CsrMatrixTest, MatrixVectorMultiplication) {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y;
    
    csr_real_.mat_vec(x, y);
    
    EXPECT_EQ(y.size(), 3);
    EXPECT_NEAR(y[0], 5.0, 1e-10);  // 1*1 + 2*2 = 5
    EXPECT_NEAR(y[1], 6.0, 1e-10);  // 3*2 = 6
    EXPECT_NEAR(y[2], 12.0, 1e-10); // 4*3 = 12
}

TEST_F(CsrMatrixTest, CSRStructure) {
    const auto& row_ptr = csr_real_.get_row_ptr();
    const auto& col_indices = csr_real_.get_col_indices();
    const auto& values = csr_real_.get_values();
    
    EXPECT_EQ(row_ptr.size(), 4);  // 3行 + 1个结束标记
    EXPECT_EQ(col_indices.size(), 4);
    EXPECT_EQ(values.size(), 4);
    
    // 验证行偏移正确
    EXPECT_EQ(row_ptr[0], 0);
    EXPECT_EQ(row_ptr[1], 2);  // 第0行有2个元素
    EXPECT_EQ(row_ptr[2], 3);  // 第1行有1个元素
    EXPECT_EQ(row_ptr[3], 4);  // 第2行有1个元素
}

/**
 * @brief 测试MatrixMarket I/O功能
 */
class MatrixMarketIOTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试矩阵
        coo_write_real_ = CooMatrixReal(3, 3);
        coo_write_real_.add_value(0, 0, 1.0);
        coo_write_real_.add_value(1, 1, 2.0);
        coo_write_real_.add_value(2, 2, 3.0);
        coo_write_real_.add_value(0, 1, 0.5);
        
        filename_ = "test_matrix_real_gtest.mtx";
    }
    
    void TearDown() override {
        // 清理测试文件
        std::remove(filename_.c_str());
    }
    
    CooMatrixReal coo_write_real_;
    CooMatrixReal coo_read_real_;
    CooMatrixComplex coo_read_complex_;
    std::string filename_;
};

TEST_F(MatrixMarketIOTest, WriteReadRealMatrix) {
    // 写入文件
    MatrixMarketIO::write_coo(filename_, coo_write_real_);
    
    // 读取文件
    MatrixDataType data_type = MatrixMarketIO::read_coo(filename_, coo_read_real_, coo_read_complex_);
    
    EXPECT_EQ(data_type, MatrixDataType::REAL);
    EXPECT_EQ(coo_read_real_.rows(), 3);
    EXPECT_EQ(coo_read_real_.cols(), 3);
    EXPECT_EQ(coo_read_real_.nnz(), 4);
    
    // 验证数据一致性
    const auto& rows_read = coo_read_real_.get_row_indices();
    const auto& cols_read = coo_read_real_.get_col_indices();
    const auto& values_read = coo_read_real_.get_values();
    
    std::map<std::pair<int, int>, double> expected_data = {
        {{0, 0}, 1.0},
        {{0, 1}, 0.5},
        {{1, 1}, 2.0},
        {{2, 2}, 3.0}
    };
    
    for (int i = 0; i < coo_read_real_.nnz(); ++i) {
        auto key = std::make_pair(rows_read[i], cols_read[i]);
        auto it = expected_data.find(key);
        EXPECT_NE(it, expected_data.end()) << "位置(" << rows_read[i] << "," << cols_read[i] << ")的数据不存在";
        if (it != expected_data.end()) {
            EXPECT_NEAR(values_read[i], it->second, 1e-10) << "位置(" << rows_read[i] << "," << cols_read[i] << ")的值不正确";
        }
    }
}

/**
 * @brief 测试稀疏矩阵基类功能
 */
TEST(SparseBaseTest, PolymorphicBehavior) {
    CooMatrixReal coo(2, 2);
    coo.add_value(0, 0, 1.0);
    coo.add_value(1, 1, 2.0);
    
    CsrMatrixReal csr(2, 2);
    csr.build_from_coo(coo);
    
    // 测试多态行为
    SparseMatrixBase* matrix1 = &coo;
    SparseMatrixBase* matrix2 = &csr;
    
    EXPECT_EQ(matrix1->rows(), 2);
    EXPECT_EQ(matrix1->cols(), 2);
    EXPECT_EQ(matrix1->nnz(), 2);
    EXPECT_EQ(matrix1->get_data_type(), MatrixDataType::REAL);
    
    EXPECT_EQ(matrix2->rows(), 2);
    EXPECT_EQ(matrix2->cols(), 2);
    EXPECT_EQ(matrix2->nnz(), 2);
    EXPECT_EQ(matrix2->get_data_type(), MatrixDataType::REAL);
}

/**
 * @brief 主函数
 */
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}