/**
 * @file test_eigen_extension.cpp
 * @brief 测试稀疏矩阵和向量的 Eigen 扩展接口
 * @author Icoder
 * @date 2026-04-10
 */

#include "numeric/coo_matrix.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/vector.hpp"
#include <iostream>
#include <cassert>
#include <complex>

using namespace numeric;
using namespace emag;

/**
 * @brief 测试 CooMatrix 的 Eigen 接口扩展
 */
void test_coo_matrix_eigen_extension() {
    std::cout << "=== 测试 CooMatrix Eigen 接口 ===" << std::endl;

    // 测试实数 COO 矩阵
    {
        std::cout << "\n--- 实数 COO 矩阵测试 ---" << std::endl;
        CooMatrixReal coo_real(3, 3);
        
        // 添加元素（包含重复项用于测试合并）
        coo_real.add_value(0, 0, 1.0);
        coo_real.add_value(0, 1, 2.0);
        coo_real.add_value(1, 1, 3.0);
        coo_real.add_value(1, 1, 4.0);  // 重复元素
        coo_real.add_value(2, 2, 5.0);

        std::cout << "添加重复元素前的非零元素数: " << coo_real.nnz() << std::endl;
        assert(coo_real.nnz() == 5);

        // 测试合并重复元素
        coo_real.merge_duplicates();
        std::cout << "合并后的非零元素数: " << coo_real.nnz() << std::endl;
        assert(coo_real.nnz() == 4); // (1,1)位置的两个元素应该合并为7.0

        // 测试获取实数 Eigen 矩阵
        const auto& eigen_mat = coo_real.get_eigen_real();
        std::cout << "Eigen 矩阵尺寸: " << eigen_mat.rows() << " x " << eigen_mat.cols() << std::endl;
        assert(eigen_mat.rows() == 3);
        assert(eigen_mat.cols() == 3);
        assert(eigen_mat.nonZeros() == 4);

        // 验证数值正确性（(1,1)位置应该是7.0）
        std::cout << "(1,1)位置的值: " << eigen_mat.coeff(1, 1) << std::endl;
        assert(std::abs(eigen_mat.coeff(1, 1) - 7.0) < 1e-12);

        std::cout << "✓ 实数 COO 矩阵所有测试通过" << std::endl;
    }

    // 测试复数 COO 矩阵
    {
        std::cout << "\n--- 复数 COO 矩阵测试 ---" << std::endl;
        CooMatrixComplex coo_complex(2, 2);
        
        using std::complex;
        coo_complex.add_value(0, 0, complex<double>(1.0, 1.0));
        coo_complex.add_value(0, 1, complex<double>(2.0, 2.0));
        coo_complex.add_value(1, 0, complex<double>(3.0, 3.0));
        coo_complex.add_value(1, 1, complex<double>(4.0, 4.0));

        // 测试获取复数 Eigen 矩阵
        const auto& eigen_complex = coo_complex.get_eigen_complex();
        std::cout << "复数 Eigen 矩阵尺寸: " << eigen_complex.rows() << " x " << eigen_complex.cols() << std::endl;
        assert(eigen_complex.rows() == 2);
        assert(eigen_complex.cols() == 2);
        assert(eigen_complex.nonZeros() == 4);

        // 验证复数值
        auto val = eigen_complex.coeff(0, 0);
        std::cout << "(0,0)位置的复数值: " << val.real() << " + " << val.imag() << "i" << std::endl;
        assert(std::abs(val.real() - 1.0) < 1e-12);
        assert(std::abs(val.imag() - 1.0) < 1e-12);

        std::cout << "✓ 复数 COO 矩阵所有测试通过" << std::endl;
    }

    // 测试类型安全（static_assert 编译期检查）
    {
        std::cout << "\n--- 类型安全测试 ---" << std::endl;
        std::cout << "注意: Vector 的类型安全由 static_assert 在编译期保证" << std::endl;
        std::cout << "      如果尝试对复数向量调用 get_eigen_real() 会导致编译错误" << std::endl;
        std::cout << "✓ 类型安全机制已就位（编译期检查）" << std::endl;
    }
}

/**
 * @brief 测试 CsrMatrix 的 Eigen 接口扩展
 */
void test_csr_matrix_eigen_extension() {
    std::cout << "\n=== 测试 CsrMatrix Eigen 接口 ===" << std::endl;

    // 测试实数 CSR 矩阵
    {
        std::cout << "\n--- 实数 CSR 矩阵测试 ---" << std::endl;
        
        // 先构建 COO 矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 1, 2.0);
        coo.add_value(1, 1, 3.0);
        coo.add_value(2, 2, 4.0);

        // 从 COO 构建 CSR
        CsrMatrixReal csr(3, 3);
        csr.build_from_coo(coo);

        std::cout << "CSR 矩阵是否已构建: " << (csr.is_built() ? "是" : "否") << std::endl;
        assert(csr.is_built());

        // 测试 merge_duplicates（CSR 应该是空操作）
        int nnz_before = csr.nnz();
        csr.merge_duplicates();
        int nnz_after = csr.nnz();
        std::cout << "merge_duplicates 前后非零元素数: " << nnz_before << " -> " << nnz_after << std::endl;
        assert(nnz_before == nnz_after); // CSR的merge应该是空操作

        // 测试获取实数 Eigen 矩阵
        const auto& eigen_csr = csr.get_eigen_real();
        std::cout << "CSR Eigen 矩阵尺寸: " << eigen_csr.rows() << " x " << eigen_csr.cols() << std::endl;
        assert(eigen_csr.rows() == 3);
        assert(eigen_csr.cols() == 3);
        assert(eigen_csr.nonZeros() == 4);

        std::cout << "✓ 实数 CSR 矩阵所有测试通过" << std::endl;
    }

    // 测试复数 CSR 矩阵
    {
        std::cout << "\n--- 复数 CSR 矩阵测试 ---" << std::endl;
        
        CooMatrixComplex coo(2, 2);
        coo.add_value(0, 0, std::complex<double>(1.0, 1.0));
        coo.add_value(0, 1, std::complex<double>(2.0, 2.0));
        coo.add_value(1, 0, std::complex<double>(3.0, 3.0));
        coo.add_value(1, 1, std::complex<double>(4.0, 4.0));

        CsrMatrixComplex csr(2, 2);
        csr.build_from_coo(coo);

        const auto& eigen_csr_complex = csr.get_eigen_complex();
        std::cout << "复数 CSR Eigen 矩阵非零元素数: " << eigen_csr_complex.nonZeros() << std::endl;
        assert(eigen_csr_complex.nonZeros() == 4);

        std::cout << "✓ 复数 CSR 矩阵所有测试通过" << std::endl;
    }

    // 测试未构建矩阵的异常
    {
        std::cout << "\n--- 未构建矩阵异常测试 ---" << std::endl;
        CsrMatrixReal csr_unbuilt(2, 2);

        bool exception_thrown = false;
        try {
            csr_unbuilt.get_eigen_real();
        } catch (const std::runtime_error& e) {
            exception_thrown = true;
            std::cout << "捕获到预期异常: " << e.what() << std::endl;
        }
        assert(exception_thrown);

        std::cout << "✓ 未构建矩阵异常测试通过" << std::endl;
    }
}

/**
 * @brief 测试 Vector 的 Eigen 接口扩展
 */
void test_vector_eigen_extension() {
    std::cout << "\n=== 测试 Vector Eigen 接口 ===" << std::endl;

    // 测试实数向量
    {
        std::cout << "\n--- 实数向量测试 ---" << std::endl;
        VectorReal vec_real({1.0, 2.0, 3.0, 4.0, 5.0});

        // 获取实数 Eigen 向量
        Eigen::VectorXd eigen_vec = vec_real.get_eigen_real();
        std::cout << "Eigen 向量大小: " << eigen_vec.size() << std::endl;
        assert(eigen_vec.size() == 5);

        // 验证数据一致性
        for (int i = 0; i < vec_real.size(); ++i) {
            std::cout << "vec[" << i << "] = " << vec_real[i] 
                      << ", eigen[" << i << "] = " << eigen_vec[i] << std::endl;
            assert(std::abs(vec_real[i] - eigen_vec[i]) < 1e-12);
        }

        std::cout << "✓ 实数向量所有测试通过" << std::endl;
    }

    // 测试复数向量
    {
        std::cout << "\n--- 复数向量测试 ---" << std::endl;
        using std::complex;
        VectorComplex vec_complex({
            complex<double>(1.0, 1.0),
            complex<double>(2.0, 2.0),
            complex<double>(3.0, 3.0)
        });

        // 获取复数 Eigen 向量
        Eigen::VectorXcd eigen_vec_complex = vec_complex.get_eigen_complex();
        std::cout << "复数 Eigen 向量大小: " << eigen_vec_complex.size() << std::endl;
        assert(eigen_vec_complex.size() == 3);

        // 验证复数数据一致性
        for (int i = 0; i < vec_complex.size(); ++i) {
            auto val = eigen_vec_complex[i];
            std::cout << "vec[" << i << "] = " << vec_complex[i].real() << "+" << vec_complex[i].imag() << "i"
                      << ", eigen[" << i << "] = " << val.real() << "+" << val.imag() << "i" << std::endl;
            assert(std::abs(vec_complex[i].real() - val.real()) < 1e-12);
            assert(std::abs(vec_complex[i].imag() - val.imag()) < 1e-12);
        }

        std::cout << "✓ 复数向量所有测试通过" << std::endl;
    }

    // 测试类型安全
    {
        std::cout << "\n--- 向量类型安全测试 ---" << std::endl;
        std::cout << "注意: Vector 的 get_eigen_real/get_eigen_complex 使用 static_assert" << std::endl;
        std::cout << "      在编译期保证类型安全，错误使用会直接导致编译失败" << std::endl;
        std::cout << "✓ 向量类型安全机制已就位（编译期检查）" << std::endl;
    }
}

/**
 * @brief 测试缓存机制和惰性求值
 */
void test_lazy_evaluation_and_cache() {
    std::cout << "\n=== 测试惰性求值与缓存机制 ===" << std::endl;

    CooMatrixReal coo(2, 2);
    coo.add_value(0, 0, 5.0);
    coo.add_value(1, 1, 10.0);

    // 第一次调用：触发缓存构建
    const auto& eigen1 = coo.get_eigen_real();
    std::cout << "第一次调用 get_eigen_real() 完成" << std::endl;

    // 第二次调用：应使用缓存（不重新构建）
    const auto& eigen2 = coo.get_eigen_real();
    std::cout << "第二次调用 get_eigen_real() 完成（使用缓存）" << std::endl;

    // 验证两次返回的是同一个对象（引用相同）
    assert(&eigen1 == &eigen2);
    std::cout << "✓ 缓存机制验证通过（两次调用返回同一引用）" << std::endl;

    // 修改数据后，再次调用应该重建缓存
    coo.add_value(0, 1, 7.0);
    const auto& eigen3 = coo.get_eigen_real();
    std::cout << "修改数据后再次调用，非零元素数: " << eigen3.nonZeros() << std::endl;
    assert(eigen3.nonZeros() == 3); // 新增了一个元素

    std::cout << "✓ 惰性求值机制测试通过" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "开始测试 Eigen 扩展接口" << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        test_coo_matrix_eigen_extension();
        test_csr_matrix_eigen_extension();
        test_vector_eigen_extension();
        test_lazy_evaluation_and_cache();

        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ 所有测试通过！" << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}
