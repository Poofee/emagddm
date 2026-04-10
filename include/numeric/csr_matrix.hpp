/**
 * @file csr_matrix.hpp
 * @brief CSR格式稀疏矩阵类定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "sparse_base.hpp"
#include "coo_matrix.hpp"
#include <vector>
#include <complex>
#include <stdexcept>

// 前向声明
namespace emag {
    template<typename T>
    class Vector;
}

namespace numeric {

/**
 * @class CsrMatrix
 * @brief CSR格式稀疏矩阵类
 * @details 使用压缩行存储格式，适合矩阵求解阶段
 */
template<typename T>
class CsrMatrix : public SparseMatrixBase {
public:
    /**
     * @brief 默认构造函数
     */
    CsrMatrix();

    /**
     * @brief 构造函数，指定矩阵尺寸
     * @param rows 矩阵行数
     * @param cols 矩阵列数
     */
    CsrMatrix(int rows, int cols);

    /**
     * @brief 析构函数
     */
    ~CsrMatrix() override = default;

    /**
     * @brief 获取矩阵行数
     * @return 矩阵行数
     */
    int rows() const override;

    /**
     * @brief 获取矩阵列数
     * @return 矩阵列数
     */
    int cols() const override;

    /**
     * @brief 获取非零元素数量
     * @return 非零元素数量
     */
    int nnz() const override;

    /**
     * @brief 清空矩阵
     */
    void clear() override;

    /**
     * @brief 打印矩阵信息
     */
    void print_info() const override;

    /**
     * @brief 获取矩阵数据类型
     * @return 矩阵数据类型
     */
    MatrixDataType get_data_type() const override;

    /**
     * @brief 从COO格式构建CSR矩阵
     * @param coo COO格式矩阵
     */
    void build_from_coo(const CooMatrix<T>& coo);

    /**
     * @brief 矩阵向量乘法（使用std::vector）
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const std::vector<T>& x, std::vector<T>& y) const;

    /**
     * @brief 矩阵向量乘法（使用Vector类）
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const emag::Vector<T>& x, emag::Vector<T>& y) const;

    /**
     * @brief 获取行偏移数组
     * @return 行偏移数组引用
     */
    const std::vector<int>& get_row_ptr() const;

    /**
     * @brief 获取列索引数组
     * @return 列索引数组引用
     */
    const std::vector<int>& get_col_indices() const;

    /**
     * @brief 获取数值数组
     * @return 数值数组引用
     */
    const std::vector<T>& get_values() const;

    /**
     * @brief 矩阵数乘
     * @param alpha 缩放因子
     */
    void scale(T alpha);

    /**
     * @brief 获取矩阵转置
     * @return 转置后的CSR矩阵
     */
    CsrMatrix<T> transpose() const;

    /**
     * @brief 提取对角线元素
     * @param diag 对角线元素向量
     */
    void get_diag(std::vector<T>& diag) const;

    /**
     * @brief 设置对角线元素
     * @param diag 对角线元素向量
     */
    void set_diag(const std::vector<T>& diag);

    /**
     * @brief 检查矩阵是否已构建
     * @return 是否已构建
     */
    bool is_built() const;

    /**
     * @brief 合并重复元素
     * @details CSR 格式本身不存储重复元素（构建时已合并），此方法为接口一致性提供空实现
     */
    void merge_duplicates() override;

    /**
     * @brief 获取实数 Eigen 稀疏矩阵
     * @return 实数 Eigen 稀疏矩阵的常量引用
     * @note 如果矩阵为复数类型，将抛出异常
     */
    const Eigen::SparseMatrix<double>& get_eigen_real() const override;

    /**
     * @brief 获取复数 Eigen 稀疏矩阵
     * @return 复数 Eigen 稀疏矩阵的常量引用
     * @note 如果矩阵为实数类型，将抛出异常
     */
    const Eigen::SparseMatrix<std::complex<double>>& get_eigen_complex() const override;

private:
    int rows_;                    ///< 矩阵行数
    int cols_;                    ///< 矩阵列数
    std::vector<int> row_ptr_;    ///< 行偏移数组
    std::vector<int> col_indices_; ///< 列索引数组
    std::vector<T> values_;        ///< 数值数组
    bool built_;                   ///< 矩阵是否已构建
    
    // Eigen 矩阵缓存（mutable 允许在 const 方法中更新）
    mutable Eigen::SparseMatrix<double> eigen_real_cache_;      ///< 实数 Eigen 矩阵缓存
    mutable Eigen::SparseMatrix<std::complex<double>> eigen_complex_cache_; ///< 复数 Eigen 矩阵缓存
    mutable bool eigen_real_dirty_ = true;      ///< 实数缓存是否需要更新
    mutable bool eigen_complex_dirty_ = true;   ///< 复数缓存是否需要更新

    /**
     * @brief 对COO数据进行排序（按行主序，行内按列排序）
     * @param rows 行索引数组
     * @param cols 列索引数组
     * @param values 数值数组
     * @param nnz 非零元素数量
     */
    void sort_coo_data(std::vector<int>& rows, std::vector<int>& cols,
                      std::vector<T>& values, int nnz) const;

    /**
     * @brief 更新实数 Eigen 缓存
     * @details 将当前 CSR 格式数据转换为 Eigen::SparseMatrix<double> 格式。
     */
    void update_eigen_real_cache() const {
        if constexpr (std::is_same_v<T, double>) {
            if (!built_) {
                throw std::runtime_error("CSR 矩阵未构建，无法生成 Eigen 矩阵");
            }

            // 使用三元组格式
            using Triplet = Eigen::Triplet<double>;
            std::vector<Triplet> triplets;
            triplets.reserve(nnz());

            for (int i = 0; i < rows_; ++i) {
                for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
                    triplets.emplace_back(i, col_indices_[j], values_[j]);
                }
            }

            // 构建 Eigen 稀疏矩阵
            eigen_real_cache_.resize(rows_, cols_);
            eigen_real_cache_.setFromTriplets(triplets.begin(), triplets.end());
            eigen_real_cache_.makeCompressed();

            eigen_real_dirty_ = false;
        } else {
            throw std::runtime_error("update_eigen_real_cache 仅支持实数类型 (double)");
        }
    }

    /**
     * @brief 更新复数 Eigen 缓存
     * @details 将当前 CSR 格式数据转换为 Eigen::SparseMatrix<std::complex<double>> 格式。
     */
    void update_eigen_complex_cache() const {
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            if (!built_) {
                throw std::runtime_error("CSR 矩阵未构建，无法生成 Eigen 矩阵");
            }

            // 转换为三元组格式
            using Triplet = Eigen::Triplet<std::complex<double>>;
            std::vector<Triplet> triplets;
            triplets.reserve(nnz());

            for (int i = 0; i < rows_; ++i) {
                for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
                    triplets.emplace_back(i, col_indices_[j], values_[j]);
                }
            }

            // 构建 Eigen 复数稀疏矩阵
            eigen_complex_cache_.resize(rows_, cols_);
            eigen_complex_cache_.setFromTriplets(triplets.begin(), triplets.end());
            eigen_complex_cache_.makeCompressed();

            eigen_complex_dirty_ = false;
        } else {
            throw std::runtime_error("update_eigen_complex_cache 仅支持复数类型 (std::complex<double>)");
        }
    }
};

// 常用类型别名
using CsrMatrixReal = CsrMatrix<double>;
using CsrMatrixComplex = CsrMatrix<std::complex<double>>;

} // namespace numeric

// 包含模板实现
#include "csr_matrix_impl.hpp"