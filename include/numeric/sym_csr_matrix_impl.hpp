#pragma once

#include "sym_csr_matrix.hpp"
#include <algorithm>
#include <functional>
#include <iostream>

namespace numeric {

template<typename T>
SymCsrMatrix<T>::SymCsrMatrix() : size_(0), built_(false) {
    row_ptr_.push_back(0); // 初始化行偏移
}

template<typename T>
SymCsrMatrix<T>::SymCsrMatrix(int size) : size_(size), built_(false) {
    if (size <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
    row_ptr_.resize(size + 1, 0); // 初始化行偏移数组
}

template<typename T>
int SymCsrMatrix<T>::rows() const {
    return size_;
}

template<typename T>
int SymCsrMatrix<T>::cols() const {
    return size_;
}

template<typename T>
int SymCsrMatrix<T>::nnz() const {
    return static_cast<int>(values_.size());
}

template<typename T>
void SymCsrMatrix<T>::clear() {
    row_ptr_.clear();
    col_indices_.clear();
    values_.clear();
    built_ = false;
    row_ptr_.push_back(0); // 重新初始化
}

template<typename T>
void SymCsrMatrix<T>::print_info() const {
    std::cout << "对称CSR矩阵信息:" << std::endl;
    std::cout << "  尺寸: " << size_ << " x " << size_ << std::endl;
    std::cout << "  非零元素数量: " << nnz() << std::endl;
    std::cout << "  数据类型: " << (get_data_type() == MatrixDataType::REAL ? "实数" : "复数") << std::endl;
    std::cout << "  构建状态: " << (built_ ? "已构建" : "未构建") << std::endl;
}

template<typename T>
MatrixDataType SymCsrMatrix<T>::get_data_type() const {
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
void SymCsrMatrix<T>::build_from_coo(const CooMatrix<T>& coo) {
    if (coo.rows() != size_ || coo.cols() != size_) {
        throw std::invalid_argument("COO矩阵尺寸与对称CSR矩阵尺寸不匹配");
    }
    
    // 验证矩阵对称性
    validate_symmetry(coo);
    
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
    
    // 构建对称CSR格式（仅下三角）
    row_ptr_.resize(size_ + 1, 0);
    col_indices_.resize(nnz);
    values_.resize(nnz);
    
    // 统计每行的非零元素数量
    for (int i = 0; i < nnz; ++i) {
        row_ptr_[rows[i] + 1]++;
    }
    
    // 计算行偏移
    for (int i = 0; i < size_; ++i) {
        row_ptr_[i + 1] += row_ptr_[i];
    }
    
    // 填充列索引和数值
    std::vector<int> row_count(size_, 0); // 每行当前填充位置
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
void SymCsrMatrix<T>::build_from_full_csr(const CsrMatrix<T>& csr) {
    if (!csr.is_built()) {
        throw std::runtime_error("输入CSR矩阵未构建");
    }
    if (csr.rows() != size_ || csr.cols() != size_) {
        throw std::invalid_argument("CSR矩阵尺寸与对称CSR矩阵尺寸不匹配");
    }
    
    clear();
    
    // 提取下三角部分
    CooMatrix<T> coo(size_, size_);
    const auto& row_ptr = csr.get_row_ptr();
    const auto& col_indices = csr.get_col_indices();
    const auto& csr_values = csr.get_values();
    
    for (int i = 0; i < size_; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = col_indices[j];
            if (col <= i) { // 仅取下三角部分（包括对角线）
                coo.add_value(i, col, csr_values[j]);
            }
        }
    }
    
    build_from_coo(coo);
}

template<typename T>
void SymCsrMatrix<T>::mat_vec(const emag::Vector<T>& x, emag::Vector<T>& y) const {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法进行矩阵向量乘法");
    }
    if (x.size() != size_) {
        throw std::invalid_argument("输入向量尺寸与矩阵尺寸不匹配");
    }
    
    y.resize(size_);
    y.set_zero();
    
    // 对称矩阵向量乘法：y = A * x
    // 由于只存储下三角，需要利用对称性
    for (int i = 0; i < size_; ++i) {
        // 处理下三角部分（包括对角线）
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            int col = col_indices_[j];
            T value = values_[j];
            
            if (col == i) {
                // 对角线元素
                y[i] += value * x[col];
            } else {
                // 下三角非对角线元素
                y[i] += value * x[col];
                // 利用对称性：A[i][col] = A[col][i]
                y[col] += value * x[i];
            }
        }
    }
}

template<typename T>
void SymCsrMatrix<T>::mat_vec(const std::vector<T>& x, std::vector<T>& y) const {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法进行矩阵向量乘法");
    }
    if (x.size() != static_cast<size_t>(size_)) {
        throw std::invalid_argument("输入向量尺寸与矩阵尺寸不匹配");
    }
    
    y.resize(size_);
    std::fill(y.begin(), y.end(), T(0));
    
    // 对称矩阵向量乘法：y = A * x
    for (int i = 0; i < size_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            int col = col_indices_[j];
            T value = values_[j];
            
            if (col == i) {
                // 对角线元素
                y[i] += value * x[col];
            } else {
                // 下三角非对角线元素
                y[i] += value * x[col];
                y[col] += value * x[i];
            }
        }
    }
}

template<typename T>
void SymCsrMatrix<T>::scale(T alpha) {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法进行数乘");
    }
    
    for (auto& value : values_) {
        value *= alpha;
    }
}

template<typename T>
const SymCsrMatrix<T>& SymCsrMatrix<T>::transpose() const {
    // 对称矩阵的转置等于自身
    return *this;
}

template<typename T>
void SymCsrMatrix<T>::get_diag(emag::Vector<T>& diag) const {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法提取对角线");
    }
    
    diag.resize(size_);
    diag.set_zero();
    
    for (int i = 0; i < size_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            if (col_indices_[j] == i) {
                diag[i] = values_[j];
                break;
            }
        }
    }
}

template<typename T>
void SymCsrMatrix<T>::set_diag(const emag::Vector<T>& diag) {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法设置对角线");
    }
    if (diag.size() != size_) {
        throw std::invalid_argument("对角线向量尺寸不匹配");
    }
    
    for (int i = 0; i < size_; ++i) {
        bool found = false;
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            if (col_indices_[j] == i) {
                values_[j] = diag[i];
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("对角线元素不存在，需要重建矩阵");
        }
    }
}

template<typename T>
bool SymCsrMatrix<T>::is_built() const {
    return built_;
}

template<typename T>
const std::vector<int>& SymCsrMatrix<T>::get_row_ptr() const {
    return row_ptr_;
}

template<typename T>
const std::vector<int>& SymCsrMatrix<T>::get_col_indices() const {
    return col_indices_;
}

template<typename T>
const std::vector<T>& SymCsrMatrix<T>::get_values() const {
    return values_;
}

template<typename T>
CsrMatrix<T> SymCsrMatrix<T>::to_full_csr() const {
    if (!built_) {
        throw std::runtime_error("对称CSR矩阵未构建，无法转换为完整CSR");
    }
    
    // 创建完整COO矩阵
    CooMatrix<T> coo_full(size_, size_);
    
    for (int i = 0; i < size_; ++i) {
        for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
            int col = col_indices_[j];
            T value = values_[j];
            
            // 添加下三角元素
            coo_full.add_value(i, col, value);
            
            // 利用对称性添加上三角元素（对角线只加一次）
            if (col != i) {
                coo_full.add_value(col, i, value);
            }
        }
    }
    
    // 构建完整CSR矩阵
    CsrMatrix<T> full_csr(size_, size_);
    full_csr.build_from_coo(coo_full);
    
    return full_csr;
}

template<typename T>
void SymCsrMatrix<T>::sort_coo_data(std::vector<int>& rows, std::vector<int>& cols, 
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

template<typename T>
void SymCsrMatrix<T>::validate_symmetry(const CooMatrix<T>& coo) const {
    // 简化实现：检查所有元素是否都在下三角（包括对角线）
    const auto& rows = coo.get_row_indices();
    const auto& cols = coo.get_col_indices();
    
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i] < cols[i]) {
            // 发现上三角元素，抛出异常
            throw std::invalid_argument("输入矩阵包含上三角元素，不是对称矩阵的下三角表示");
        }
    }
}

// 显式实例化
template class SymCsrMatrix<double>;
template class SymCsrMatrix<std::complex<double>>;

} // namespace numeric