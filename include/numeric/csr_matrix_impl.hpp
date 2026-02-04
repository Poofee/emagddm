/**
 * @file csr_matrix_impl.hpp
 * @brief CSR格式稀疏矩阵模板实现
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <functional>

namespace numeric {

template<typename T>
CsrMatrix<T>::CsrMatrix() : rows_(0), cols_(0), built_(false) {
    row_ptr_.push_back(0); // 初始化行偏移
}

template<typename T>
CsrMatrix<T>::CsrMatrix(int rows, int cols) : rows_(rows), cols_(cols), built_(false) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
    row_ptr_.resize(rows + 1, 0); // 初始化行偏移数组
}

template<typename T>
int CsrMatrix<T>::rows() const {
    return rows_;
}

template<typename T>
int CsrMatrix<T>::cols() const {
    return cols_;
}

template<typename T>
int CsrMatrix<T>::nnz() const {
    return static_cast<int>(values_.size());
}

template<typename T>
void CsrMatrix<T>::clear() {
    row_ptr_.clear();
    col_indices_.clear();
    values_.clear();
    built_ = false;
    row_ptr_.push_back(0); // 重新初始化
}

template<typename T>
void CsrMatrix<T>::print_info() const {
    std::cout << "CSR矩阵信息:" << std::endl;
    std::cout << "  尺寸: " << rows_ << " x " << cols_ << std::endl;
    std::cout << "  非零元素数量: " << nnz() << std::endl;
    std::cout << "  数据类型: " << (get_data_type() == MatrixDataType::REAL ? "实数" : "复数") << std::endl;
    std::cout << "  构建状态: " << (built_ ? "已构建" : "未构建") << std::endl;
}

template<typename T>
MatrixDataType CsrMatrix<T>::get_data_type() const {
    if constexpr (std::is_same_v<T, double>) {
        return MatrixDataType::REAL;
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return MatrixDataType::COMPLEX;
    } else {
        static_assert(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>, 
                      "仅支持double和std::complex<double>类型");
        return MatrixDataType::REAL;
    }
}

template<typename T>
void CsrMatrix<T>::build_from_coo(const CooMatrix<T>& coo) {
    if (coo.rows() != rows_ || coo.cols() != cols_) {
        throw std::invalid_argument("COO矩阵尺寸与CSR矩阵尺寸不匹配");
    }

    clear();
    
    int nnz = coo.nnz();
    if (nnz == 0) {
        built_ = true;
        return;
    }

    // 获取COO数据
    std::vector<int> rows = coo.get_row_indices();
    std::vector<int> cols = coo.get_col_indices();
    std::vector<T> values = coo.get_values();

    // 对COO数据进行排序
    sort_coo_data(rows, cols, values, nnz);

    // 构建CSR格式
    row_ptr_.resize(rows_ + 1, 0);
    col_indices_.resize(nnz);
    values_.resize(nnz);

    // 统计每行的非零元素数量
    for (int i = 0; i < nnz; ++i) {
        row_ptr_[rows[i] + 1]++;
    }

    // 计算行偏移
    for (int i = 0; i < rows_; ++i) {
        row_ptr_[i + 1] += row_ptr_[i];
    }

    // 填充列索引和数值
    std::vector<int> row_count(rows_, 0); // 每行当前填充位置
    for (int i = 0; i < nnz; ++i) {
        int row = rows[i];
        int pos = row_ptr_[row] + row_count[row];
        col_indices_[pos] = cols[i];
        values_[pos] = values[i];
        row_count[row]++;
    }

    built_ = true;
}

template<typename T>
void CsrMatrix<T>::mat_vec(const std::vector<T>& x, std::vector<T>& y) const {
    if (!built_) {
        throw std::runtime_error("CSR矩阵未构建，无法进行矩阵向量乘法");
    }
    if (x.size() != static_cast<size_t>(cols_)) {
        throw std::invalid_argument("输入向量尺寸与矩阵列数不匹配");
    }

    y.resize(rows_);
    std::fill(y.begin(), y.end(), T(0));

    for (int i = 0; i < rows_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            y[i] += values_[j] * x[col_indices_[j]];
        }
    }
}

template<typename T>
const std::vector<int>& CsrMatrix<T>::get_row_ptr() const {
    return row_ptr_;
}

template<typename T>
const std::vector<int>& CsrMatrix<T>::get_col_indices() const {
    return col_indices_;
}

template<typename T>
const std::vector<T>& CsrMatrix<T>::get_values() const {
    return values_;
}

template<typename T>
void CsrMatrix<T>::scale(T alpha) {
    if (!built_) {
        throw std::runtime_error("CSR矩阵未构建，无法进行数乘");
    }

    for (auto& value : values_) {
        value *= alpha;
    }
}

template<typename T>
CsrMatrix<T> CsrMatrix<T>::transpose() const {
    if (!built_) {
        throw std::runtime_error("CSR矩阵未构建，无法转置");
    }

    CsrMatrix<T> transposed(cols_, rows_);
    
    // 使用COO格式进行转置
    CooMatrix<T> coo(cols_, rows_);
    coo.reserve(nnz());
    
    for (int i = 0; i < rows_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            coo.add_value(col_indices_[j], i, values_[j]);
        }
    }
    
    transposed.build_from_coo(coo);
    return transposed;
}

template<typename T>
void CsrMatrix<T>::get_diag(std::vector<T>& diag) const {
    if (!built_) {
        throw std::runtime_error("CSR矩阵未构建，无法提取对角线");
    }

    diag.resize(std::min(rows_, cols_));
    std::fill(diag.begin(), diag.end(), T(0));

    for (int i = 0; i < rows_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            if (col_indices_[j] == i) {
                diag[i] = values_[j];
                break;
            }
        }
    }
}

template<typename T>
void CsrMatrix<T>::set_diag(const std::vector<T>& diag) {
    if (!built_) {
        throw std::runtime_error("CSR矩阵未构建，无法设置对角线");
    }
    if (diag.size() != static_cast<size_t>(std::min(rows_, cols_))) {
        throw std::invalid_argument("对角线向量尺寸不匹配");
    }

    for (int i = 0; i < rows_; ++i) {
        bool found = false;
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            if (col_indices_[j] == i) {
                values_[j] = diag[i];
                found = true;
                break;
            }
        }
        if (!found) {
            // 如果对角线元素不存在，需要插入（这里简化处理，实际应该重建矩阵）
            throw std::runtime_error("对角线元素不存在，需要重建矩阵");
        }
    }
}

template<typename T>
bool CsrMatrix<T>::is_built() const {
    return built_;
}

template<typename T>
void CsrMatrix<T>::sort_coo_data(std::vector<int>& rows, std::vector<int>& cols, 
                                 std::vector<T>& values, int nnz) const {
    // 创建索引数组
    std::vector<int> indices(nnz);
    for (int i = 0; i < nnz; ++i) {
        indices[i] = i;
    }

    // 按行主序，行内按列排序
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        if (rows[a] != rows[b]) {
            return rows[a] < rows[b];
        }
        return cols[a] < cols[b];
    });

    // 重新排列数据
    std::vector<int> sorted_rows(nnz);
    std::vector<int> sorted_cols(nnz);
    std::vector<T> sorted_values(nnz);

    for (int i = 0; i < nnz; ++i) {
        sorted_rows[i] = rows[indices[i]];
        sorted_cols[i] = cols[indices[i]];
        sorted_values[i] = values[indices[i]];
    }

    rows = std::move(sorted_rows);
    cols = std::move(sorted_cols);
    values = std::move(sorted_values);
}

} // namespace numeric