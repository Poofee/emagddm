/**
 * @file em_sparse_converter.h
 * @brief 核心数值层 - CSR矩阵与Eigen稀疏矩阵转换工具类
 * @details 提供项目自定义CsrMatrix<double>与Eigen::SparseMatrix<double>之间的双向无缝转换功能，
 *          支持性能统计、维度校验、大型矩阵优化等特性
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

#include "csr_matrix.hpp"
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <chrono>
#include <stdexcept>

namespace numeric {

/**
 * @class SparseConverter
 * @brief CSR矩阵与Eigen稀疏矩阵的双向转换工具类
 * @details 实现numeric::CsrMatrix<double>与Eigen::SparseMatrix<double>之间的高效双向转换，
 *          提供维度校验、性能统计、内存优化等功能。使用静态方法设计，无需实例化即可调用。
 *
 * @note 性能特征：
 *       - 大型矩阵(nnz>10000)自动启用预分配优化
 *       - 使用std::chrono高精度计时
 *       - 支持往返转换精度验证(||A_original - A_roundtrip||_F < 1e-15)
 *
 * @code
 * // CsrMatrix转Eigen示例
 * numeric::CsrMatrix<double> csr_matrix = ...;
 * Eigen::SparseMatrix<double> eigen_mat = numeric::SparseConverter::to_eigen(csr_matrix);
 *
 * // Eigen转CsrMatrix示例
 * Eigen::SparseMatrix<double> eigen_mat = ...;
 * numeric::CsrMatrix<double> csr_matrix = numeric::SparseConverter::from_eigen(eigen_mat);
 * @endcode
 */
class SparseConverter {
public:
    /**
     * @brief 转换结果统计信息结构体
     */
    struct ConversionStats {
        double conversion_time_ms;   ///< 转换耗时（毫秒）
        int source_nnz;              ///< 源矩阵非零元素数
        int target_nnz;              ///< 目标矩阵非零元素数
        bool dimensions_match;       ///< 维度是否匹配
    };

    /**
     * @brief 将CsrMatrix转换为Eigen稀疏矩阵
     * @param csr 输入的CSR格式矩阵（const引用，避免拷贝）
     * @return Eigen::SparseMatrix<double> 转换后的Eigen稀疏矩阵（按行主序存储）
     * @throw std::invalid_argument 矩阵未构建或维度不合法时抛出
     * @note 算法步骤：
     *       1. 调用validate_matrix_dimensions()校验输入合法性
     *       2. 提取CSR内部数据(row_ptr, col_indices, values)
     *       3. 使用Eigen::Triplets构建稀疏矩阵
     *       4. 统计转换耗时并记录到ConversionStats
     * @note 性能优化：
     *       - 对于nnz>10000的大型矩阵，使用reserve()预分配空间
     *       - 使用移动语义减少拷贝开销
     */
    static Eigen::SparseMatrix<double> to_eigen(const CsrMatrix<double>& csr);

    /**
     * @brief 将Eigen稀疏矩阵转换为CsrMatrix
     * @param eigen_mat 输入的Eigen稀疏矩阵（const引用）
     * @return CsrMatrix<double> 转换后的CSR格式矩阵
     * @throw std::invalid_argument Eigen矩阵为空或维度不合法时抛出
     * @note 算法步骤：
     *       1. 校验Eigen矩阵是否已压缩且维度合法
     *       2. 提取Eigen内部CSC格式数据(outer_ptr, inner_idx, values)
     *       3. 将CSC格式转换为CSR格式并构造CsrMatrix对象
     *       4. 保证非零元素数量和维度一致性
     *       5. 统计转换耗时
     * @note Eigen默认使用CSC(压缩列存储)格式，本方法会正确转换为CSR格式
     */
    static CsrMatrix<double> from_eigen(const Eigen::SparseMatrix<double>& eigen_mat);

    /**
     * @brief 校验矩阵维度合法性
     * @param csr 待校验的CSR矩阵
     * @param expected_rows 期望行数（-1表示不检查）
     * @param expected_cols 期望列数（-1表示不检查）
     * @throw std::invalid_argument 以下情况抛出异常：
     *        - 矩阵未构建(is_built()==false)
     *        - 行数或列数<=0
     *        - 提供了expected_rows/expected_cols但维度不匹配
     * @note 建议在to_eigen()/from_eigen()前手动调用以提前发现错误
     */
    static void validate_matrix_dimensions(const CsrMatrix<double>& csr,
                                          int expected_rows = -1,
                                          int expected_cols = -1);

    /**
     * @brief 获取上次转换操作的统计信息
     * @return ConversionStats 上次转换的统计数据（耗时、非零元数、维度匹配状态）
     * @note 在调用to_eigen()或from_eigen()后使用，用于性能监控和调试
     */
    static ConversionStats get_last_conversion_stats();

private:
    static ConversionStats last_stats_;  ///< 内部状态：上次转换统计信息

    /**
     * @brief 将Eigen的CSC格式转换为CSR格式（私有辅助方法）
     * @param eigen_mat 输入的Eigen稀疏矩阵（需已压缩）
     * @param row_ptr 输出行偏移数组
     * @param col_indices 输出列索引数组
     * @param values 输出数值数组
     * @note Eigen默认使用CSC存储，此方法实现CSC到CSR的高效转换
     */
    static void convert_csc_to_csr(const Eigen::SparseMatrix<double>& eigen_mat,
                                   std::vector<int>& row_ptr,
                                   std::vector<int>& col_indices,
                                   std::vector<double>& values);
};

} // namespace numeric
