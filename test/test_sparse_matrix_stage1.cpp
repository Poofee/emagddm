/**
 * @file test_sparse_matrix_stage1.cpp
 * @brief 稀疏矩阵模块阶段1单元测试
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#include "numeric/sparse_base.hpp"
#include "numeric/coo_matrix.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/matrix_market_io.hpp"
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>

using namespace numeric;

/**
 * @brief 断言宏，用于简单的测试断言
 */
#define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "测试失败: " << message << " (文件: " << __FILE__ << ", 行: " << __LINE__ << ")" << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 测试COO矩阵基础功能
 */
bool testCooMatrixBasic() {
    std::cout << "=== 测试COO矩阵基础功能 ===" << std::endl;
    
    CooMatrixReal coo(3, 3);
    ASSERT(coo.rows() == 3, "行数不正确");
    ASSERT(coo.cols() == 3, "列数不正确");
    ASSERT(coo.nnz() == 0, "非零元素数量不正确");
    ASSERT(coo.get_data_type() == MatrixDataType::REAL, "数据类型不正确");
    
    std::cout << "✓ COO矩阵基础构造测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试COO矩阵添加元素功能
 */
bool testCooMatrixAddValues() {
    std::cout << "=== 测试COO矩阵添加元素功能 ===" << std::endl;
    
    CooMatrixReal coo(3, 3);
    coo.add_value(0, 0, 1.0);
    coo.add_value(1, 1, 2.0);
    coo.add_value(2, 2, 3.0);
    
    ASSERT(coo.nnz() == 3, "非零元素数量不正确");
    
    const auto& rows = coo.get_row_indices();
    const auto& cols = coo.get_col_indices();
    const auto& values = coo.get_values();
    
    ASSERT(rows[0] == 0, "行索引0不正确");
    ASSERT(cols[0] == 0, "列索引0不正确");
    ASSERT(std::abs(values[0] - 1.0) < 1e-10, "值0不正确");
    
    ASSERT(rows[1] == 1, "行索引1不正确");
    ASSERT(cols[1] == 1, "列索引1不正确");
    ASSERT(std::abs(values[1] - 2.0) < 1e-10, "值1不正确");
    
    ASSERT(rows[2] == 2, "行索引2不正确");
    ASSERT(cols[2] == 2, "列索引2不正确");
    ASSERT(std::abs(values[2] - 3.0) < 1e-10, "值2不正确");
    
    std::cout << "✓ COO矩阵添加元素测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试复数COO矩阵
 */
bool testCooMatrixComplex() {
    std::cout << "=== 测试复数COO矩阵 ===" << std::endl;
    
    CooMatrixComplex coo(2, 2);
    coo.add_value(0, 0, std::complex<double>(1.0, 2.0));
    coo.add_value(1, 1, std::complex<double>(3.0, 4.0));
    
    ASSERT(coo.nnz() == 2, "非零元素数量不正确");
    ASSERT(coo.get_data_type() == MatrixDataType::COMPLEX, "数据类型不正确");
    
    const auto& values = coo.get_values();
    ASSERT(std::abs(values[0].real() - 1.0) < 1e-10, "实部0不正确");
    ASSERT(std::abs(values[0].imag() - 2.0) < 1e-10, "虚部0不正确");
    ASSERT(std::abs(values[1].real() - 3.0) < 1e-10, "实部1不正确");
    ASSERT(std::abs(values[1].imag() - 4.0) < 1e-10, "虚部1不正确");
    
    std::cout << "✓ 复数COO矩阵测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试CSR矩阵构建和矩阵向量乘法
 */
bool testCsrMatrixBuildAndMatVec() {
    std::cout << "=== 测试CSR矩阵构建和矩阵向量乘法 ===" << std::endl;
    
    // 创建COO矩阵
    CooMatrixReal coo(3, 3);
    coo.add_value(0, 0, 1.0);
    coo.add_value(1, 1, 2.0);
    coo.add_value(2, 2, 3.0);
    coo.add_value(0, 1, 0.5);
    coo.add_value(1, 0, 0.5);
    
    // 构建CSR矩阵
    CsrMatrixReal csr(3, 3);
    csr.build_from_coo(coo);
    
    ASSERT(csr.is_built(), "CSR矩阵未构建");
    ASSERT(csr.nnz() == 5, "非零元素数量不正确");
    
    // 验证矩阵向量乘法
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y;
    csr.mat_vec(x, y);
    
    ASSERT(std::abs(y[0] - (1.0 * 1.0 + 0.5 * 2.0)) < 1e-10, "y[0]计算结果不正确");
    ASSERT(std::abs(y[1] - (0.5 * 1.0 + 2.0 * 2.0)) < 1e-10, "y[1]计算结果不正确");
    ASSERT(std::abs(y[2] - (3.0 * 3.0)) < 1e-10, "y[2]计算结果不正确");
    
    std::cout << "✓ CSR矩阵构建和矩阵向量乘法测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试MatrixMarket I/O功能
 */
bool testMatrixMarketIO() {
    std::cout << "=== 测试MatrixMarket I/O功能 ===" << std::endl;
    
    try {
        // 创建测试矩阵
        CooMatrixReal coo_write(3, 3);
        coo_write.add_value(0, 0, 1.0);
        coo_write.add_value(1, 1, 2.0);
        coo_write.add_value(2, 2, 3.0);
        coo_write.add_value(0, 1, 0.5);
        
        // 写入文件
        const std::string filename = "test_matrix_real.mtx";
        std::cout << "写入MatrixMarket文件: " << filename << std::endl;
        MatrixMarketIO::write_coo(filename, coo_write);
        
        // 检查文件是否存在并打印内容
        std::ifstream check_file(filename);
        if (!check_file.is_open()) {
            std::cerr << "错误：文件写入失败，文件不存在: " << filename << std::endl;
            return false;
        }
        
        std::cout << "写入的文件内容:" << std::endl;
        std::string line;
        while (std::getline(check_file, line)) {
            std::cout << "  " << line << std::endl;
        }
        check_file.close();
        
        // 读取文件
        CooMatrixReal coo_read_real;
        CooMatrixComplex coo_read_complex;
        
        std::cout << "读取MatrixMarket文件..." << std::endl;
        MatrixDataType data_type = MatrixMarketIO::read_coo(filename, coo_read_real, coo_read_complex);
        
        ASSERT(data_type == MatrixDataType::REAL, "数据类型不正确");
        ASSERT(coo_read_real.rows() == 3, "行数不正确");
        ASSERT(coo_read_real.cols() == 3, "列数不正确");
        ASSERT(coo_read_real.nnz() == 4, "非零元素数量不正确");
        
        // 验证数据一致性（不依赖顺序）
        const auto& rows_read = coo_read_real.get_row_indices();
        const auto& cols_read = coo_read_real.get_col_indices();
        const auto& values_read = coo_read_real.get_values();
        
        // 创建映射来验证数据
        std::map<std::pair<int, int>, double> expected_data = {
            {{0, 0}, 1.0},
            {{0, 1}, 0.5},
            {{1, 1}, 2.0},
            {{2, 2}, 3.0}
        };
        
        for (int i = 0; i < coo_read_real.nnz(); ++i) {
            auto key = std::make_pair(rows_read[i], cols_read[i]);
            auto it = expected_data.find(key);
            ASSERT(it != expected_data.end(), "位置(" + std::to_string(rows_read[i]) + "," + std::to_string(cols_read[i]) + ")的数据不存在");
            ASSERT(std::abs(values_read[i] - it->second) < 1e-10, "位置(" + std::to_string(rows_read[i]) + "," + std::to_string(cols_read[i]) + ")的值不正确");
        }
        
        // 清理测试文件
        std::remove(filename.c_str());
        
        std::cout << "✓ MatrixMarket I/O功能测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "MatrixMarket I/O测试异常: " << e.what() << std::endl;
        return false;
    } catch (...) {
        std::cerr << "MatrixMarket I/O测试未知异常" << std::endl;
        return false;
    }
}

/**
 * @brief 测试边界条件
 */
bool testEdgeCases() {
    std::cout << "=== 测试边界条件 ===" << std::endl;
    
    // 空矩阵测试
    CooMatrixReal coo_empty(0, 0);
    ASSERT(coo_empty.rows() == 0, "空矩阵行数不正确");
    ASSERT(coo_empty.cols() == 0, "空矩阵列数不正确");
    ASSERT(coo_empty.nnz() == 0, "空矩阵非零元素数量不正确");
    
    CsrMatrixReal csr_empty(0, 0);
    csr_empty.build_from_coo(coo_empty);
    ASSERT(csr_empty.is_built(), "空CSR矩阵未构建");
    ASSERT(csr_empty.nnz() == 0, "空CSR矩阵非零元素数量不正确");
    
    // 单元素矩阵测试
    CooMatrixReal coo_single(1, 1);
    coo_single.add_value(0, 0, 42.0);
    
    CsrMatrixReal csr_single(1, 1);
    csr_single.build_from_coo(coo_single);
    
    std::vector<double> x = {2.0};
    std::vector<double> y;
    csr_single.mat_vec(x, y);
    
    ASSERT(std::abs(y[0] - 84.0) < 1e-10, "单元素矩阵乘法结果不正确");
    
    std::cout << "✓ 边界条件测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试异常处理
 */
bool testExceptionHandling() {
    std::cout << "=== 测试异常处理 ===" << std::endl;
    
    CooMatrixReal coo(2, 2);
    
    // 测试无效索引
    bool caught = false;
    try {
        coo.add_value(-1, 0, 1.0);
    } catch (const std::out_of_range&) {
        caught = true;
    }
    ASSERT(caught, "负行索引异常未正确抛出");
    
    caught = false;
    try {
        coo.add_value(0, -1, 1.0);
    } catch (const std::out_of_range&) {
        caught = true;
    }
    ASSERT(caught, "负列索引异常未正确抛出");
    
    caught = false;
    try {
        coo.add_value(2, 0, 1.0);
    } catch (const std::out_of_range&) {
        caught = true;
    }
    ASSERT(caught, "超出行范围索引异常未正确抛出");
    
    caught = false;
    try {
        coo.add_value(0, 2, 1.0);
    } catch (const std::out_of_range&) {
        caught = true;
    }
    ASSERT(caught, "超出列范围索引异常未正确抛出");
    
    std::cout << "✓ 异常处理测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始稀疏矩阵模块阶段1测试..." << std::endl;
    
    bool all_passed = true;
    
    all_passed = testCooMatrixBasic() && all_passed;
    all_passed = testCooMatrixAddValues() && all_passed;
    all_passed = testCooMatrixComplex() && all_passed;
    all_passed = testCsrMatrixBuildAndMatVec() && all_passed;
    all_passed = testMatrixMarketIO() && all_passed;
    all_passed = testEdgeCases() && all_passed;
    all_passed = testExceptionHandling() && all_passed;
    
    if (all_passed) {
        std::cout << "\n🎉 所有测试通过！稀疏矩阵模块阶段1功能正常" << std::endl;
        return 0;
    } else {
        std::cout << "\n❌ 部分测试失败，请检查实现" << std::endl;
        return 1;
    }
}