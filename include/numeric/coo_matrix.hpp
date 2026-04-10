/**
 * @file coo_matrix.hpp
 * @brief COO格式稀疏矩阵类定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "sparse_base.hpp"
#include <vector>
#include <complex>
#include <stdexcept>

namespace numeric {

/**
 * @class CooMatrix
 * @brief COO格式稀疏矩阵类
 * @details 使用坐标格式存储稀疏矩阵，适合矩阵组装阶段
 */
template<typename T>
class CooMatrix : public SparseMatrixBase {
public:
    /**
     * @brief 默认构造函数
     */
    CooMatrix();

    /**
     * @brief 构造函数，指定矩阵尺寸
     * @param rows 矩阵行数
     * @param cols 矩阵列数
     */
    CooMatrix(int rows, int cols);

    /**
     * @brief 构造函数，指定矩阵尺寸和非零元素容量
     * @param rows 矩阵行数
     * @param cols 矩阵列数
     * @param capacity 预分配的非零元素容量
     */
    CooMatrix(int rows, int cols, int capacity);

    /**
     * @brief 析构函数
     */
    ~CooMatrix() override = default;

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
     * @brief 添加非零元素
     * @param row 行索引
     * @param col 列索引
     * @param value 元素值
     */
    void add_value(int row, int col, T value);

    /**
     * @brief 批量添加非零元素
     * @param rows 行索引数组
     * @param cols 列索引数组
     * @param values 元素值数组
     * @param count 元素数量
     */
    void add_values(const std::vector<int>& rows, const std::vector<int>& cols, 
                   const std::vector<T>& values, int count);

    /**
     * @brief 设置矩阵尺寸
     * @param rows 行数
     * @param cols 列数
     */
    void set_size(int rows, int cols);

    /**
     * @brief 预留非零元素容量
     * @param capacity 容量大小
     */
    void reserve(int capacity);

    /**
     * @brief 获取行索引数组
     * @return 行索引数组引用
     */
    const std::vector<int>& get_row_indices() const;

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
     * @brief 检查索引是否有效
     * @param row 行索引
     * @param col 列索引
     * @return 是否有效
     */
    bool is_valid_index(int row, int col) const;

    /**
     * @brief 合并重复元素
     * @details 合并相同 (row, col) 位置的元素，值累加
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
    std::vector<int> row_indices_; ///< 行索引数组
    std::vector<int> col_indices_; ///< 列索引数组
    std::vector<T> values_;        ///< 数值数组
    
    // Eigen 矩阵缓存（mutable 允许在 const 方法中更新）
    mutable Eigen::SparseMatrix<double> eigen_real_cache_;      ///< 实数 Eigen 矩阵缓存
    mutable Eigen::SparseMatrix<std::complex<double>> eigen_complex_cache_; ///< 复数 Eigen 矩阵缓存
    mutable bool eigen_real_dirty_ = true;      ///< 实数缓存是否需要更新
    mutable bool eigen_complex_dirty_ = true;   ///< 复数缓存是否需要更新

    /**
     * @brief 更新实数 Eigen 缓存
     * @details 将当前 COO 格式数据转换为 Eigen::SparseMatrix<double> 格式。
     *          使用三元组格式 (Triplets) 构建稀疏矩阵。
     */
    void update_eigen_real_cache() const {
        if constexpr (std::is_same_v<T, double>) {
            // 创建三元组列表用于构建 Eigen 稀疏矩阵
            using Triplet = Eigen::Triplet<double>;
            std::vector<Triplet> triplets;
            triplets.reserve(nnz());

            // 填充三元组：(行, 列, 值)
            for (int i = 0; i < nnz(); ++i) {
                triplets.emplace_back(row_indices_[i], col_indices_[i], values_[i]);
            }

            // 构建 Eigen 稀疏矩阵
            eigen_real_cache_.resize(rows_, cols_);
            eigen_real_cache_.setFromTriplets(triplets.begin(), triplets.end());
            eigen_real_cache_.makeCompressed();

            // 标记缓存已更新
            eigen_real_dirty_ = false;
        } else {
            throw std::runtime_error("update_eigen_real_cache 仅支持实数类型 (double)");
        }
    }

    /**
     * @brief 更新复数 Eigen 缓存
     * @details 将当前 COO 格式数据转换为 Eigen::SparseMatrix<std::complex<double>> 格式。
     */
    void update_eigen_complex_cache() const {
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            // 创建三元组列表用于构建 Eigen 复数稀疏矩阵
            using Triplet = Eigen::Triplet<std::complex<double>>;
            std::vector<Triplet> triplets;
            triplets.reserve(nnz());

            // 填充三元组：(行, 列, 复数值)
            for (int i = 0; i < nnz(); ++i) {
                triplets.emplace_back(row_indices_[i], col_indices_[i], values_[i]);
            }

            // 构建 Eigen 复数稀疏矩阵
            eigen_complex_cache_.resize(rows_, cols_);
            eigen_complex_cache_.setFromTriplets(triplets.begin(), triplets.end());
            eigen_complex_cache_.makeCompressed();

            // 标记缓存已更新
            eigen_complex_dirty_ = false;
        } else {
            throw std::runtime_error("update_eigen_complex_cache 仅支持复数类型 (std::complex<double>)");
        }
    }
};

// 常用类型别名
using CooMatrixReal = CooMatrix<double>;
using CooMatrixComplex = CooMatrix<std::complex<double>>;

} // namespace numeric

// 包含模板实现
#include "coo_matrix_impl.hpp"