/**
 * @file complex_matrix_ops.hpp
 * @brief 复数矩阵操作工具类
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 * 
 * @details 提供复数矩阵特有的操作功能，如共轭转置、复数矩阵向量乘法等
 *          适配涡流场等复数矩阵场景
 */

#pragma once

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include "sym_csr_matrix.hpp"
#include "matrix_attribute.hpp"
#include <complex>
#include <vector>
#include <functional>

namespace numeric {

/**
 * @class ComplexMatrixOps
 * @brief 复数矩阵操作工具类
 * @details 提供复数矩阵特有的操作功能
 */
class ComplexMatrixOps {
public:
    /**
     * @brief 计算复数矩阵的共轭转置
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @return 共轭转置矩阵
     */
    template<typename T>
    static CsrMatrix<T> conjugate_transpose(const CsrMatrix<T>& matrix) {
        // 对于实数矩阵，转置就是共轭转置
        if constexpr (std::is_same_v<T, double>) {
            return transpose(matrix);
        } else {
            // 复数矩阵的共轭转置
            CsrMatrix<T> result(matrix.cols(), matrix.rows());
            
            // 构建转置矩阵的COO格式
            CooMatrix<T> coo_trans(matrix.cols(), matrix.rows());
            
            const auto& row_ptr = matrix.get_row_ptr();
            const auto& col_indices = matrix.get_col_indices();
            const auto& values = matrix.get_values();
            
            for (int i = 0; i < matrix.rows(); ++i) {
                int start = row_ptr[i];
                int end = row_ptr[i + 1];
                
                for (int j = start; j < end; ++j) {
                    int col = col_indices[j];
                    T value = std::conj(values[j]); // 取共轭
                    coo_trans.add_value(col, i, value); // 行列互换
                }
            }
            
            result.build_from_coo(coo_trans);
            return result;
        }
    }
    
    /**
     * @brief 计算矩阵的转置（实数/复数通用）
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @return 转置矩阵
     */
    template<typename T>
    static CsrMatrix<T> transpose(const CsrMatrix<T>& matrix) {
        CsrMatrix<T> result(matrix.cols(), matrix.rows());
        
        // 构建转置矩阵的COO格式
        CooMatrix<T> coo_trans(matrix.cols(), matrix.rows());
        
        const auto& row_ptr = matrix.get_row_ptr();
        const auto& col_indices = matrix.get_col_indices();
        const auto& values = matrix.get_values();
        
        for (int i = 0; i < matrix.rows(); ++i) {
            int start = row_ptr[i];
            int end = row_ptr[i + 1];
            
            for (int j = start; j < end; ++j) {
                int col = col_indices[j];
                coo_trans.add_value(col, i, values[j]); // 行列互换
            }
        }
        
        result.build_from_coo(coo_trans);
        return result;
    }
    
    /**
     * @brief 复数矩阵向量乘法（支持共轭）
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param x 输入向量
     * @param y 输出向量
     * @param conjugate 是否使用共轭乘法
     */
    template<typename T>
    static void complex_mat_vec(const CsrMatrix<T>& matrix, 
                               const std::vector<T>& x, 
                               std::vector<T>& y, 
                               bool conjugate = false) {
        if (x.size() != static_cast<size_t>(matrix.cols())) {
            throw std::invalid_argument("输入向量尺寸与矩阵列数不匹配");
        }
        
        y.resize(matrix.rows(), T(0));
        
        const auto& row_ptr = matrix.get_row_ptr();
        const auto& col_indices = matrix.get_col_indices();
        const auto& values = matrix.get_values();
        
        for (int i = 0; i < matrix.rows(); ++i) {
            int start = row_ptr[i];
            int end = row_ptr[i + 1];
            
            T sum = T(0);
            for (int j = start; j < end; ++j) {
                int col = col_indices[j];
                T value = values[j];
                
                if constexpr (std::is_same_v<T, std::complex<double>>) {
                    if (conjugate) {
                        sum += std::conj(value) * x[col]; // 共轭乘法
                    } else {
                        sum += value * x[col];
                    }
                } else {
                    sum += value * x[col]; // 实数矩阵
                }
            }
            y[i] = sum;
        }
    }
    
    /**
     * @brief 检查矩阵是否为埃尔米特矩阵
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param tolerance 容差
     * @return 是否为埃尔米特矩阵
     */
    template<typename T>
    static bool is_hermitian(const CsrMatrix<T>& matrix, double tolerance = 1e-10) {
        if (matrix.rows() != matrix.cols()) {
            return false; // 非方阵不可能是埃尔米特矩阵
        }
        
        // 对于实数矩阵，埃尔米特矩阵就是对称矩阵
        if constexpr (std::is_same_v<T, double>) {
            return is_symmetric(matrix, tolerance);
        } else {
            // 复数矩阵：检查 A = A^H
            auto A_H = conjugate_transpose(matrix);
            return is_equal(matrix, A_H, tolerance);
        }
    }
    
    /**
     * @brief 检查矩阵是否为对称矩阵
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param tolerance 容差
     * @return 是否为对称矩阵
     */
    template<typename T>
    static bool is_symmetric(const CsrMatrix<T>& matrix, double tolerance = 1e-10) {
        if (matrix.rows() != matrix.cols()) {
            return false; // 非方阵不可能是对称矩阵
        }
        
        auto A_T = transpose(matrix);
        return is_equal(matrix, A_T, tolerance);
    }
    
    /**
     * @brief 检查两个矩阵是否相等（在容差范围内）
     * @tparam T 矩阵元素类型
     * @param A 第一个矩阵
     * @param B 第二个矩阵
     * @param tolerance 容差
     * @return 是否相等
     */
    template<typename T>
    static bool is_equal(const CsrMatrix<T>& A, const CsrMatrix<T>& B, double tolerance = 1e-10) {
        if (A.rows() != B.rows() || A.cols() != B.cols()) {
            return false;
        }
        
        // 转换为COO格式进行比较
        CooMatrix<T> coo_A(A.rows(), A.cols());
        CooMatrix<T> coo_B(B.rows(), B.cols());
        
        // 将A转换为COO
        const auto& row_ptr_A = A.get_row_ptr();
        const auto& col_indices_A = A.get_col_indices();
        const auto& values_A = A.get_values();
        
        for (int i = 0; i < A.rows(); ++i) {
            int start = row_ptr_A[i];
            int end = row_ptr_A[i + 1];
            
            for (int j = start; j < end; ++j) {
                coo_A.add_value(i, col_indices_A[j], values_A[j]);
            }
        }
        
        // 将B转换为COO
        const auto& row_ptr_B = B.get_row_ptr();
        const auto& col_indices_B = B.get_col_indices();
        const auto& values_B = B.get_values();
        
        for (int i = 0; i < B.rows(); ++i) {
            int start = row_ptr_B[i];
            int end = row_ptr_B[i + 1];
            
            for (int j = start; j < end; ++j) {
                coo_B.add_value(i, col_indices_B[j], values_B[j]);
            }
        }
        
        // 比较两个COO矩阵
        return coo_A.is_equal(coo_B, tolerance);
    }
    
    /**
     * @brief 根据矩阵属性自动选择预处理策略
     * @param attr 矩阵属性
     * @return 推荐的预处理类型字符串
     */
    static std::string recommend_preconditioner(const MatrixAttribute& attr) {
        if (attr.suitable_for_block_preconditioner()) {
            return "块ILU";
        } else if (attr.suitable_for_ilu()) {
            return "ILU(0)";
        } else if (attr.suitable_for_jacobi()) {
            return "Jacobi";
        } else {
            return "无预处理";
        }
    }
    
    /**
     * @brief 根据矩阵属性自动选择求解器
     * @param attr 矩阵属性
     * @return 推荐的求解器类型字符串
     */
    static std::string recommend_solver(const MatrixAttribute& attr) {
        if (attr.suitable_for_cg()) {
            return "CG";
        } else if (attr.data_type == MatrixDataType::COMPLEX) {
            return "GMRES";
        } else if (attr.is_singular) {
            return "MINRES";
        } else {
            return "BiCGSTAB";
        }
    }
    
    /**
     * @brief 创建复数对角矩阵
     * @tparam T 矩阵元素类型
     * @param n 矩阵尺寸
     * @param diag 对角线元素
     * @return 复数对角矩阵
     */
    template<typename T>
    static CsrMatrix<T> create_complex_diagonal(int n, const std::vector<T>& diag) {
        if (diag.size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("对角线元素数量与矩阵尺寸不匹配");
        }
        
        CooMatrix<T> coo(n, n);
        for (int i = 0; i < n; ++i) {
            coo.add_value(i, i, diag[i]);
        }
        
        CsrMatrix<T> result(n, n);
        result.build_from_coo(coo);
        return result;
    }
    
    /**
     * @brief 创建涡流场系统矩阵（示例）
     * @param n 矩阵尺寸
     * @param omega 角频率
     * @param sigma 电导率
     * @param mu 磁导率
     * @return 涡流场系统矩阵
     */
    static CsrMatrix<std::complex<double>> create_eddy_current_matrix(int n, 
                                                                     double omega, 
                                                                     double sigma, 
                                                                     double mu) {
        // 创建示例涡流场矩阵：A = K + iωσM
        // 其中K是刚度矩阵，M是质量矩阵
        
        CooMatrix<std::complex<double>> coo(n, n);
        
        // 创建简单的三对角矩阵作为示例
        std::complex<double> j(0, 1); // 虚数单位
        std::complex<double> k_val(2.0, 0.0); // 刚度项
        std::complex<double> m_val(1.0, 0.0); // 质量项
        
        for (int i = 0; i < n; ++i) {
            // 对角元素
            coo.add_value(i, i, k_val + j * omega * sigma * m_val);
            
            // 非对角元素
            if (i > 0) {
                coo.add_value(i, i-1, -1.0);
            }
            if (i < n-1) {
                coo.add_value(i, i+1, -1.0);
            }
        }
        
        CsrMatrix<std::complex<double>> result(n, n);
        result.build_from_coo(coo);
        return result;
    }
};

} // namespace numeric