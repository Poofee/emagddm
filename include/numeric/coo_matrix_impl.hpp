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

namespace numeric {

template<typename T>
CooMatrix<T>::CooMatrix() : rows_(0), cols_(0) {
}

template<typename T>
CooMatrix<T>::CooMatrix(int rows, int cols) : rows_(rows), cols_(cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
}

template<typename T>
CooMatrix<T>::CooMatrix(int rows, int cols, int capacity) : rows_(rows), cols_(cols) {
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

} // namespace numeric