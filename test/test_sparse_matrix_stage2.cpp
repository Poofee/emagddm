/**
 * @file test_sparse_matrix_stage2.cpp
 * @brief 稀疏矩阵模块阶段2功能测试
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#include <iostream>
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

// 测试辅助宏
#define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "测试失败: " << message << std::endl; \
            return false; \
        } \
    } while(0)

#define ASSERT_NEAR(a, b, tolerance, message) \
    do { \
        double diff = std::abs((a) - (b)); \
        if (diff > (tolerance)) { \
            std::cerr << "测试失败: " << message << " (差值: " << diff << ")" << std::endl; \
            return false; \
        } \
    } while(0)

/**
 * @brief 测试向量类基本功能
 */
bool testVectorBasic() {
    std::cout << "=== 测试向量类基本功能 ===" << std::endl;
    
    try {
        // 测试实数向量
        VectorReal v1(5);
        ASSERT(v1.size() == 5, "向量大小不正确");
        ASSERT(v1.get_data_type() == VectorDataType::REAL, "向量数据类型不正确");
        
        // 测试向量赋值和访问
        v1[0] = 1.0;
        v1[1] = 2.0;
        v1[2] = 3.0;
        v1[3] = 4.0;
        v1[4] = 5.0;
        
        ASSERT(v1[0] == 1.0, "向量元素访问失败");
        ASSERT(v1[4] == 5.0, "向量元素访问失败");
        
        // 测试向量运算
        VectorReal v2 = {2.0, 3.0, 4.0, 5.0, 6.0};
        VectorReal v3 = v1 + v2;
        VectorReal v4 = v1 - v2;
        VectorReal v5 = v1 * 2.0;
        
        ASSERT(v3[0] == 3.0, "向量加法失败");
        ASSERT(v4[0] == -1.0, "向量减法失败");
        ASSERT(v5[0] == 2.0, "向量数乘失败");
        
        // 测试点积和范数
        double dot_result = v1.dot(v2);
        double norm_result = v1.norm();
        
        ASSERT_NEAR(dot_result, 70.0, 1e-10, "向量点积计算错误");
        ASSERT_NEAR(norm_result, std::sqrt(55.0), 1e-10, "向量范数计算错误");
        
        // 测试复数向量
        VectorComplex vc1(3);
        vc1[0] = std::complex<double>(1.0, 2.0);
        vc1[1] = std::complex<double>(3.0, 4.0);
        vc1[2] = std::complex<double>(5.0, 6.0);
        
        ASSERT(vc1.get_data_type() == VectorDataType::COMPLEX, "复数向量数据类型不正确");
        
        std::cout << "✓ 向量类基本功能测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "向量类测试异常: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 测试矩阵向量乘法
 */
bool testMatrixVectorMultiplication() {
    std::cout << "=== 测试矩阵向量乘法 ===" << std::endl;
    
    try {
        // 创建测试矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 1, 2.0);
        coo.add_value(0, 2, 3.0);
        coo.add_value(1, 0, 4.0);
        coo.add_value(1, 1, 5.0);
        coo.add_value(1, 2, 6.0);
        coo.add_value(2, 0, 7.0);
        coo.add_value(2, 1, 8.0);
        coo.add_value(2, 2, 9.0);
        
        CsrMatrixReal csr(3, 3);
        csr.build_from_coo(coo);
        
        // 测试std::vector版本
        std::vector<double> x_std = {1.0, 2.0, 3.0};
        std::vector<double> y_std;
        csr.mat_vec(x_std, y_std);
        
        ASSERT(y_std.size() == 3, "输出向量尺寸不正确");
        ASSERT_NEAR(y_std[0], 14.0, 1e-10, "矩阵向量乘法结果错误");
        ASSERT_NEAR(y_std[1], 32.0, 1e-10, "矩阵向量乘法结果错误");
        ASSERT_NEAR(y_std[2], 50.0, 1e-10, "矩阵向量乘法结果错误");
        
        // 测试Vector类版本
        VectorReal x_vec = {1.0, 2.0, 3.0};
        VectorReal y_vec;
        csr.mat_vec(x_vec, y_vec);
        
        ASSERT(y_vec.size() == 3, "输出向量尺寸不正确");
        ASSERT_NEAR(y_vec[0], 14.0, 1e-10, "矩阵向量乘法结果错误");
        ASSERT_NEAR(y_vec[1], 32.0, 1e-10, "矩阵向量乘法结果错误");
        ASSERT_NEAR(y_vec[2], 50.0, 1e-10, "矩阵向量乘法结果错误");
        
        std::cout << "✓ 矩阵向量乘法测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "矩阵向量乘法测试异常: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 测试矩阵基础操作
 */
bool testMatrixOperations() {
    std::cout << "=== 测试矩阵基础操作 ===" << std::endl;
    
    try {
        // 创建测试矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 1, 2.0);
        coo.add_value(1, 0, 3.0);
        coo.add_value(1, 1, 4.0);
        coo.add_value(2, 2, 5.0);
        
        CsrMatrixReal csr(3, 3);
        csr.build_from_coo(coo);
        
        // 测试矩阵数乘
        csr.scale(2.0);
        
        std::vector<double> x = {1.0, 1.0, 1.0};
        std::vector<double> y;
        csr.mat_vec(x, y);
        
        ASSERT_NEAR(y[0], 6.0, 1e-10, "矩阵数乘结果错误");
        ASSERT_NEAR(y[1], 14.0, 1e-10, "矩阵数乘结果错误");
        ASSERT_NEAR(y[2], 10.0, 1e-10, "矩阵数乘结果错误");
        
        // 测试矩阵转置
        CsrMatrixReal csr_trans = csr.transpose();
        
        std::vector<double> y_trans;
        csr_trans.mat_vec(x, y_trans);
        
        ASSERT_NEAR(y_trans[0], 8.0, 1e-10, "矩阵转置结果错误");
        ASSERT_NEAR(y_trans[1], 12.0, 1e-10, "矩阵转置结果错误");
        ASSERT_NEAR(y_trans[2], 10.0, 1e-10, "矩阵转置结果错误");
        
        // 测试对角线操作
        std::vector<double> diag;
        csr.get_diag(diag);
        
        ASSERT(diag.size() == 3, "对角线向量尺寸不正确");
        ASSERT_NEAR(diag[0], 2.0, 1e-10, "对角线提取错误");
        ASSERT_NEAR(diag[1], 8.0, 1e-10, "对角线提取错误");
        ASSERT_NEAR(diag[2], 10.0, 1e-10, "对角线提取错误");
        
        std::cout << "✓ 矩阵基础操作测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "矩阵基础操作测试异常: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 测试Jacobi预处理
 */
bool testJacobiPreconditioner() {
    std::cout << "=== 测试Jacobi预处理 ===" << std::endl;
    
    try {
        // 创建对角占优测试矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 4.0);
        coo.add_value(0, 1, 1.0);
        coo.add_value(1, 0, 1.0);
        coo.add_value(1, 1, 5.0);
        coo.add_value(1, 2, 2.0);
        coo.add_value(2, 1, 2.0);
        coo.add_value(2, 2, 6.0);
        
        CsrMatrixReal csr(3, 3);
        csr.build_from_coo(coo);
        
        // 创建Jacobi预处理器
        JacobiPreconditioner jacobi(csr);
        
        // 测试预处理器应用
        VectorReal r = {1.0, 2.0, 3.0};
        VectorReal z;
        jacobi.apply(r, z);
        
        // 验证结果：z ≈ diag(A)^(-1) * r
        ASSERT(z.size() == 3, "预处理输出向量尺寸不正确");
        ASSERT_NEAR(z[0], 1.0/4.0, 1e-10, "Jacobi预处理结果错误");
        ASSERT_NEAR(z[1], 2.0/5.0, 1e-10, "Jacobi预处理结果错误");
        ASSERT_NEAR(z[2], 3.0/6.0, 1e-10, "Jacobi预处理结果错误");
        
        std::cout << "✓ Jacobi预处理测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Jacobi预处理测试异常: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 测试对称矩阵优化
 */
bool testSymmetricMatrix() {
    std::cout << "=== 测试对称矩阵优化 ===" << std::endl;
    
    try {
        // 创建对称测试矩阵（仅下三角部分）
        CooMatrixReal coo_sym(3, 3);
        coo_sym.add_value(0, 0, 1.0);
        coo_sym.add_value(1, 0, 2.0);
        coo_sym.add_value(1, 1, 3.0);
        coo_sym.add_value(2, 0, 4.0);
        coo_sym.add_value(2, 1, 5.0);
        coo_sym.add_value(2, 2, 6.0);
        
        // 创建完整矩阵用于对比
        CooMatrixReal coo_full(3, 3);
        coo_full.add_value(0, 0, 1.0);
        coo_full.add_value(0, 1, 2.0);
        coo_full.add_value(0, 2, 4.0);
        coo_full.add_value(1, 0, 2.0);
        coo_full.add_value(1, 1, 3.0);
        coo_full.add_value(1, 2, 5.0);
        coo_full.add_value(2, 0, 4.0);
        coo_full.add_value(2, 1, 5.0);
        coo_full.add_value(2, 2, 6.0);
        
        // 构建对称矩阵和完整矩阵
        SymCsrMatrixReal sym_csr(3);
        sym_csr.build_from_coo(coo_sym);
        
        CsrMatrixReal full_csr(3, 3);
        full_csr.build_from_coo(coo_full);
        
        // 测试矩阵向量乘法结果一致性
        VectorReal x = {1.0, 2.0, 3.0};
        VectorReal y_sym, y_full;
        
        sym_csr.mat_vec(x, y_sym);
        full_csr.mat_vec(x, y_full);
        
        ASSERT(y_sym.size() == 3, "对称矩阵输出向量尺寸不正确");
        ASSERT(y_full.size() == 3, "完整矩阵输出向量尺寸不正确");
        
        for (int i = 0; i < 3; ++i) {
            ASSERT_NEAR(y_sym[i], y_full[i], 1e-10, 
                        "对称矩阵与完整矩阵乘法结果不一致");
        }
        
        // 测试转置操作（对称矩阵转置等于自身）
        const SymCsrMatrixReal& transposed = sym_csr.transpose();
        VectorReal y_trans;
        transposed.mat_vec(x, y_trans);
        
        for (int i = 0; i < 3; ++i) {
            ASSERT_NEAR(y_sym[i], y_trans[i], 1e-10, 
                        "对称矩阵转置结果不正确");
        }
        
        std::cout << "✓ 对称矩阵优化测试通过" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "对称矩阵测试异常: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始稀疏矩阵模块阶段2功能测试..." << std::endl;
    std::cout << "================================" << std::endl;
    
    bool all_passed = true;
    
    // 运行所有测试
    all_passed &= testVectorBasic();
    all_passed &= testMatrixVectorMultiplication();
    all_passed &= testMatrixOperations();
    all_passed &= testJacobiPreconditioner();
    all_passed &= testSymmetricMatrix();
    
    std::cout << "================================" << std::endl;
    if (all_passed) {
        std::cout << "✓ 所有阶段2功能测试通过！" << std::endl;
        return 0;
    } else {
        std::cout << "✗ 部分测试失败，请检查实现" << std::endl;
        return 1;
    }
}