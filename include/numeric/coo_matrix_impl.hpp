/**
 * @file coo_matrix_impl.hpp
 * @brief COO格式稀疏矩阵模板实现
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "coo_matrix.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <map>
#include <Eigen/Sparse>

namespace numeric {

template<typename T>
CooMatrix<T>::CooMatrix() : rows_(0), cols_(0),
    eigen_real_cache_(0, 0),
    eigen_complex_cache_(0, 0) {
}

template<typename T>
CooMatrix<T>::CooMatrix(int rows, int cols) : rows_(rows), cols_(cols),
    eigen_real_cache_(rows, cols),
    eigen_complex_cache_(rows, cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
}

template<typename T>
CooMatrix<T>::CooMatrix(int rows, int cols, int capacity) : rows_(rows), cols_(cols),
    eigen_real_cache_(rows, cols),
    eigen_complex_cache_(rows, cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
    if (capacity > 0) {
        reserve(capacity);
    }
}

template<typename T>
int CooMatrix<T>::rows() const {
    return rows_;
}

template<typename T>
int CooMatrix<T>::cols() const {
    return cols_;
}

template<typename T>
int CooMatrix<T>::nnz() const {
    return static_cast<int>(values_.size());
}

template<typename T>
void CooMatrix<T>::clear() {
    row_indices_.clear();
    col_indices_.clear();
    values_.clear();
}

template<typename T>
void CooMatrix<T>::print_info() const {
    std::cout << "COO矩阵信息:" << std::endl;
    std::cout << "  尺寸: " << rows_ << " x " << cols_ << std::endl;
    std::cout << "  非零元素数量: " << nnz() << std::endl;
    std::cout << "  数据类型: " << (get_data_type() == MatrixDataType::REAL ? "实数" : "复数") << std::endl;
    std::cout << "  存储容量: " << values_.capacity() << std::endl;
}

template<typename T>
MatrixDataType CooMatrix<T>::get_data_type() const {
    if constexpr (std::is_same_v<T, double>) {
        return MatrixDataType::REAL;
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return MatrixDataType::COMPLEX;
    } else {
        static_assert(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>, 
                      "仅支持double和std::complex<double>类型");
        return MatrixDataType::REAL; // 默认返回实数，但不会执行到这里
    }
}

template<typename T>
void CooMatrix<T>::add_value(int row, int col, T value) {
    if (!is_valid_index(row, col)) {
        throw std::out_of_range("矩阵索引超出范围");
    }
    
    row_indices_.push_back(row);
    col_indices_.push_back(col);
    values_.push_back(value);
}

template<typename T>
void CooMatrix<T>::add_values(const std::vector<int>& rows, const std::vector<int>& cols, 
                             const std::vector<T>& values, int count) {
    if (rows.size() < count || cols.size() < count || values.size() < count) {
        throw std::invalid_argument("输入数组大小不足");
    }
    
    for (int i = 0; i < count; ++i) {
        add_value(rows[i], cols[i], values[i]);
    }
}

template<typename T>
void CooMatrix<T>::set_size(int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
    rows_ = rows;
    cols_ = cols;
    clear(); // 改变尺寸时清空数据
}

template<typename T>
void CooMatrix<T>::reserve(int capacity) {
    if (capacity <= 0) {
        return;
    }
    row_indices_.reserve(capacity);
    col_indices_.reserve(capacity);
    values_.reserve(capacity);
}

template<typename T>
const std::vector<int>& CooMatrix<T>::get_row_indices() const {
    return row_indices_;
}

template<typename T>
const std::vector<int>& CooMatrix<T>::get_col_indices() const {
    return col_indices_;
}

template<typename T>
const std::vector<T>& CooMatrix<T>::get_values() const {
    return values_;
}

template<typename T>
bool CooMatrix<T>::is_valid_index(int row, int col) const {
    return row >= 0 && row < rows_ && col >= 0 && col < cols_;
}

/**
 * @brief 合并 COO 矩阵中的重复元素
 * @details 使用 std::map 按 (row, col) 键值对进行排序和合并，
 *          相同位置的元素值会累加。合并后更新内部存储并标记 Eigen 缓存为脏。
 * @tparam T 矩阵元素数据类型（double 或 std::complex<double>）
 */
template<typename T>
void CooMatrix<T>::merge_duplicates() {
    if (nnz() == 0) {
        return; // 空矩阵无需处理
    }

    // 使用 map 自动按 (row, col) 排序并合并重复项
    std::map<std::pair<int, int>, T> merged_map;
    
    for (size_t i = 0; i < row_indices_.size(); ++i) {
        auto key = std::make_pair(row_indices_[i], col_indices_[i]);
        merged_map[key] += values_[i];
    }

    // 清空原有数据并用合并后的数据填充
    row_indices_.clear();
    col_indices_.clear();
    values_.clear();

    // 预分配容量以提高性能
    size_t new_nnz = merged_map.size();
    row_indices_.reserve(new_nnz);
    col_indices_.reserve(new_nnz);
    values_.reserve(new_nnz);

    // 将 map 中的数据按顺序填充回向量（map 已自动按键排序）
    for (const auto& [key, value] : merged_map) {
        row_indices_.push_back(key.first);
        col_indices_.push_back(key.second);
        values_.push_back(value);
    }

    // 标记 Eigen 缓存需要重新构建
    eigen_real_dirty_ = true;
    eigen_complex_dirty_ = true;
}

/**
 * @brief 获取实数 Eigen 稀疏矩阵的常量引用
 * @details 返回当前 COO 矩阵对应的 Eigen 实数稀疏矩阵。采用惰性求值策略，
 *          仅在缓存失效时才重新构建 Eigen 矩阵，否则直接返回缓存的引用。
 *          这避免了频繁转换带来的性能开销。
 * @return const Eigen::SparseMatrix<double>& 实数 Eigen 稀疏矩阵的常量引用
 * @throws std::runtime_error 如果调用复数类型的矩阵实例
 * @tparam T 矩阵元素数据类型
 */
template<typename T>
const Eigen::SparseMatrix<double>& CooMatrix<T>::get_eigen_real() const {
    // 类型安全检查：确保只有实数矩阵才能调用此方法
    if constexpr (!std::is_same_v<T, double>) {
        throw std::runtime_error("无法从复数矩阵获取实数 Eigen 矩阵");
    }

    // 惰性求值：仅在缓存脏时才更新
    if (eigen_real_dirty_) {
        update_eigen_real_cache();
    }

    return eigen_real_cache_;
}

/**
 * @brief 获取复数 Eigen 稀疏矩阵的常量引用
 * @details 返回当前 COO 矩阵对应的 Eigen 复数稀疏矩阵。采用惰性求值策略，
 *          仅在缓存失效时才重新构建 Eigen 矩阵，否则直接返回缓存的引用。
 *          这避免了频繁转换带来的性能开销。
 * @return const Eigen::SparseMatrix<std::complex<double>>& 复数 Eigen 稀疏矩阵的常量引用
 * @throws std::runtime_error 如果调用实数类型的矩阵实例
 * @tparam T 矩阵元素数据类型
 */
template<typename T>
const Eigen::SparseMatrix<std::complex<double>>& CooMatrix<T>::get_eigen_complex() const {
    // 类型安全检查：确保只有复数矩阵才能调用此方法
    if constexpr (!std::is_same_v<T, std::complex<double>>) {
        throw std::runtime_error("无法从实数矩阵获取复数 Eigen 矩阵");
    }

    // 惰性求值：仅在缓存脏时才更新
    if (eigen_complex_dirty_) {
        update_eigen_complex_cache();
    }

    return eigen_complex_cache_;
}

} // namespace numeric