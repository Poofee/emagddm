/**
 * @file test_ilu0_debug.cpp
 * @brief ILU0预条件子调试测试
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#include <iostream>
#include <vector>
#include <complex>

#include "numeric/coo_matrix.hpp"
#include "numeric/csr_matrix.hpp"
#include "numeric/preconditioner.hpp"

using namespace emag;
using namespace numeric;

int main() {
    try {
        // 创建测试矩阵
        CooMatrixReal coo(3, 3);
        coo.add_value(0, 0, 4.0);
        coo.add_value(0, 1, 1.0);
        coo.add_value(1, 0, 1.0);
        coo.add_value(1, 1, 3.0);
        coo.add_value(1, 2, 1.0);
        coo.add_value(2, 1, 1.0);
        coo.add_value(2, 2, 2.0);
        
        std::cout << "COO矩阵构建成功，尺寸: " << coo.rows() << "x" << coo.cols() << std::endl;
        
        // 创建CSR矩阵
        CsrMatrixReal csr(coo.rows(), coo.cols());
        csr.build_from_coo(coo);
        
        std::cout << "CSR矩阵构建成功，尺寸: " << csr.rows() << "x" << csr.cols() << std::endl;
        std::cout << "CSR矩阵已构建: " << csr.is_built() << std::endl;
        
        // 测试Jacobi预条件子
        JacobiPreconditioner jacobi(csr, 1e-10);
        std::cout << "Jacobi预条件子构建成功" << std::endl;
        
        // 测试ILU0预条件子
        ILU0Preconditioner ilu0(csr);
        std::cout << "ILU0预条件子构建成功" << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
}