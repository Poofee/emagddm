/**
 * @file test_eigen_extension_simple.cpp
 * @brief 简化版 Eigen 扩展接口测试
 */

#include "numeric/coo_matrix.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/vector.hpp"
#include <iostream>
#include <cassert>
#include <complex>

using namespace numeric;
using namespace emag;

int main() {
    std::cout << "=== 简化版 Eigen 接口测试 ===" << std::endl;

    // 测试1: 实数 COO 矩阵
    {
        std::cout << "\n1. 测试实数 COO 矩阵..." << std::endl;
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(1, 1, 2.0);
        coo.add_value(2, 2, 3.0);

        const auto& eigen_mat = coo.get_eigen_real();
        assert(eigen_mat.rows() == 3);
        assert(eigen_mat.cols() == 3);
        std::cout << "   ✓ 实数 COO 矩阵测试通过" << std::endl;
    }

    // 测试2: 复数 COO 矩阵
    {
        std::cout << "\n2. 测试复数 COO 矩阵..." << std::endl;
        CooMatrixComplex coo(2, 2);
        coo.add_value(0, 0, std::complex<double>(1.0, 1.0));
        coo.add_value(1, 1, std::complex<double>(2.0, 2.0));

        const auto& eigen_mat = coo.get_eigen_complex();
        assert(eigen_mat.rows() == 2);
        assert(eigen_mat.cols() == 2);
        std::cout << "   ✓ 复数 COO 矩阵测试通过" << std::endl;
    }

    // 测试3: 合并重复元素
    {
        std::cout << "\n3. 测试合并重复元素..." << std::endl;
        CooMatrixReal coo(2, 2);
        coo.add_value(0, 0, 1.0);
        coo.add_value(0, 0, 2.0);  // 重复
        coo.add_value(1, 1, 3.0);

        assert(coo.nnz() == 3);
        coo.merge_duplicates();
        assert(coo.nnz() == 2);  // (0,0) 应该合并为 3.0

        const auto& eigen_mat = coo.get_eigen_real();
        assert(std::abs(eigen_mat.coeff(0, 0) - 3.0) < 1e-12);
        std::cout << "   ✓ 合并重复元素测试通过" << std::endl;
    }

    // 测试4: CSR 矩阵
    {
        std::cout << "\n4. 测试实数 CSR 矩阵..." << std::endl;
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 1.0);
        coo.add_value(1, 1, 2.0);
        coo.add_value(2, 2, 3.0);

        CsrMatrixReal csr(3, 3);
        csr.build_from_coo(coo);

        const auto& eigen_csr = csr.get_eigen_real();
        assert(eigen_csr.rows() == 3);
        assert(eigen_csr.cols() == 3);
        assert(eigen_csr.nonZeros() == 3);
        std::cout << "   ✓ 实数 CSR 矩阵测试通过" << std::endl;
    }

    // 测试5: 实数向量
    {
        std::cout << "\n5. 测试实数向量..." << std::endl;
        VectorReal vec({1.0, 2.0, 3.0, 4.0, 5.0});

        Eigen::VectorXd eigen_vec = vec.get_eigen_real();
        assert(eigen_vec.size() == 5);
        assert(std::abs(eigen_vec[0] - 1.0) < 1e-12);
        assert(std::abs(eigen_vec[4] - 5.0) < 1e-12);
        std::cout << "   ✓ 实数向量测试通过" << std::endl;
    }

    // 测试6: 复数向量
    {
        std::cout << "\n6. 测试复数向量..." << std::endl;
        VectorComplex vec({std::complex<double>(1.0, 1.0),
                           std::complex<double>(2.0, 2.0),
                           std::complex<double>(3.0, 3.0)});

        Eigen::VectorXcd eigen_vec = vec.get_eigen_complex();
        assert(eigen_vec.size() == 3);
        assert(std::abs(eigen_vec[0].real() - 1.0) < 1e-12);
        assert(std::abs(eigen_vec[0].imag() - 1.0) < 1e-12);
        std::cout << "   ✓ 复数向量测试通过" << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "✓ 所有测试通过！" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
