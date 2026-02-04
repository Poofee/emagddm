/**
 * @file block_csr_matrix_impl.hpp
 * @brief 块CSR格式稀疏矩阵模板实现
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "block_csr_matrix.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace numeric {

template<typename T>
BlockCsrMatrix<T>::BlockCsrMatrix() 
    : rows_(0), cols_(0), block_size_(BlockSize::BLOCK_1x1), block_dim_(1), built_(false) {
    block_row_ptr_.push_back(0); // 初始化块行偏移
}

template<typename T>
BlockCsrMatrix<T>::BlockCsrMatrix(int rows, int cols, BlockSize block_size) 
    : rows_(rows), cols_(cols), block_size_(block_size), built_(false) {
    
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("矩阵尺寸必须为正数");
    }
    
    // 计算块尺寸
    block_dim_ = static_cast<int>(block_size);
    
    // 初始化块行偏移数组
    block_row_ptr_.resize(rows + 1, 0);
}

template<typename T>
int BlockCsrMatrix<T>::rows() const {
    return rows_;
}

template<typename T>
int BlockCsrMatrix<T>::cols() const {
    return cols_;
}

template<typename T>
int BlockCsrMatrix<T>::nnz() const {
    return static_cast<int>(block_values_.size());
}

template<typename T>
int BlockCsrMatrix<T>::block_nnz() const {
    return static_cast<int>(block_col_indices_.size());
}

template<typename T>
void BlockCsrMatrix<T>::clear() {
    block_row_ptr_.clear();
    block_col_indices_.clear();
    block_values_.clear();
    
    // 重新初始化块行偏移
    block_row_ptr_.push_back(0);
    built_ = false;
}

template<typename T>
void BlockCsrMatrix<T>::print_info() const {
    std::cout << "块CSR矩阵信息:" << std::endl;
    std::cout << "  尺寸: " << rows_ << " × " << cols_ << " (块)" << std::endl;
    std::cout << "  块大小: " << block_dim_ << " × " << block_dim_ << std::endl;
    std::cout << "  块非零元数: " << block_nnz() << std::endl;
    std::cout << "  总非零元数: " << nnz() << std::endl;
    std::cout << "  数据类型: " << (get_data_type() == MatrixDataType::REAL ? "实数" : "复数") << std::endl;
    std::cout << "  已构建: " << (built_ ? "是" : "否") << std::endl;
}

template<typename T>
MatrixDataType BlockCsrMatrix<T>::get_data_type() const {
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
BlockSize BlockCsrMatrix<T>::get_block_size() const {
    return block_size_;
}

template<typename T>
int BlockCsrMatrix<T>::get_block_dim() const {
    return block_dim_;
}

template<typename T>
void BlockCsrMatrix<T>::build_from_coo(const CooMatrix<T>& coo) {
    if (coo.rows() != rows_ * block_dim_ || coo.cols() != cols_ * block_dim_) {
        throw std::invalid_argument("COO矩阵尺寸与块CSR矩阵尺寸不匹配");
    }
    
    clear();
    
    int total_nnz = coo.nnz();
    if (total_nnz == 0) {
        built_ = true;
        return;
    }
    
    // 获取COO矩阵数据
    const auto& coo_rows = coo.get_row_indices();
    const auto& coo_cols = coo.get_col_indices();
    const auto& coo_values = coo.get_values();
    
    // 将元素按块行主序排序
    std::vector<int> indices(total_nnz);
    for (int i = 0; i < total_nnz; ++i) {
        indices[i] = i;
    }
    
    // 按块行、块列排序
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        int block_row_a = coo_rows[a] / block_dim_;
        int block_row_b = coo_rows[b] / block_dim_;
        if (block_row_a != block_row_b) {
            return block_row_a < block_row_b;
        }
        int block_col_a = coo_cols[a] / block_dim_;
        int block_col_b = coo_cols[b] / block_dim_;
        return block_col_a < block_col_b;
    });
    
    // 构建块CSR结构
    block_row_ptr_.resize(rows_ + 1, 0);
    int current_block_row = -1;
    
    for (int i = 0; i < total_nnz; ++i) {
        int idx = indices[i];
        int block_row = coo_rows[idx] / block_dim_;
        int block_col = coo_cols[idx] / block_dim_;
        
        if (block_row != current_block_row) {
            for (int r = current_block_row + 1; r <= block_row; ++r) {
                block_row_ptr_[r] = block_col_indices_.size();
            }
            current_block_row = block_row;
        }
        
        // 检查是否为新块
        if (block_col_indices_.empty() || block_col_indices_.back() != block_col) {
            block_col_indices_.push_back(block_col);
            // 为新区块预留空间
            block_values_.resize(block_values_.size() + block_dim_ * block_dim_, T(0));
        }
        
        // 设置块内元素值
        int block_index = block_col_indices_.size() - 1;
        int elem_row = coo_rows[idx] % block_dim_;
        int elem_col = coo_cols[idx] % block_dim_;
        int elem_index = compute_element_index(block_index, elem_row, elem_col);
        block_values_[elem_index] = coo_values[idx];
    }
    
    // 设置剩余行的偏移
    for (int r = current_block_row + 1; r <= rows_; ++r) {
        block_row_ptr_[r] = block_col_indices_.size();
    }
    
    built_ = true;
}

template<typename T>
void BlockCsrMatrix<T>::build_from_csr(const CsrMatrix<T>& csr, BlockSize block_size) {
    if (csr.rows() % static_cast<int>(block_size) != 0 || 
        csr.cols() % static_cast<int>(block_size) != 0) {
        throw std::invalid_argument("CSR矩阵尺寸与块大小不兼容");
    }
    
    rows_ = csr.rows() / static_cast<int>(block_size);
    cols_ = csr.cols() / static_cast<int>(block_size);
    block_size_ = block_size;
    block_dim_ = static_cast<int>(block_size);
    
    // 先转换为COO，再构建块CSR
    CooMatrix<T> coo(csr.rows(), csr.cols());
    
    const auto& row_ptr = csr.get_row_ptr();
    const auto& col_indices = csr.get_col_indices();
    const auto& values = csr.get_values();
    
    for (int i = 0; i < csr.rows(); ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            coo.add_value(i, col_indices[j], values[j]);
        }
    }
    
    build_from_coo(coo);
}

template<typename T>
void BlockCsrMatrix<T>::mat_vec(const emag::Vector<T>& x, emag::Vector<T>& y) const {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法执行矩阵向量乘法");
    }
    
    if (x.size() != cols_ * block_dim_) {
        throw std::invalid_argument("输入向量尺寸与矩阵列数不匹配");
    }
    
    y.resize(rows_ * block_dim_);
    
    // 使用std::vector进行块级计算
    std::vector<T> x_vec(x.size());
    std::vector<T> y_vec(y.size(), T(0));
    
    for (int i = 0; i < x.size(); ++i) {
        x_vec[i] = x[i];
    }
    
    block_mat_vec(x_vec, y_vec);
    
    for (int i = 0; i < y.size(); ++i) {
        y[i] = y_vec[i];
    }
}

template<typename T>
void BlockCsrMatrix<T>::mat_vec(const std::vector<T>& x, std::vector<T>& y) const {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法执行矩阵向量乘法");
    }
    
    if (x.size() != cols_ * block_dim_) {
        throw std::invalid_argument("输入向量尺寸与矩阵列数不匹配");
    }
    
    y.resize(rows_ * block_dim_, T(0));
    block_mat_vec(x, y);
}

template<typename T>
void BlockCsrMatrix<T>::block_mat_vec(const std::vector<T>& x_block, std::vector<T>& y_block) const {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法执行块级矩阵向量乘法");
    }
    
    // 块级矩阵向量乘法
    for (int block_row = 0; block_row < rows_; ++block_row) {
        int block_start = block_row_ptr_[block_row];
        int block_end = block_row_ptr_[block_row + 1];
        
        for (int block_idx = block_start; block_idx < block_end; ++block_idx) {
            int block_col = block_col_indices_[block_idx];
            
            // 块乘法：y_block += A_block * x_block
            for (int i = 0; i < block_dim_; ++i) {
                int y_index = block_row * block_dim_ + i;
                
                for (int j = 0; j < block_dim_; ++j) {
                    int x_index = block_col * block_dim_ + j;
                    int elem_index = compute_element_index(block_idx, i, j);
                    
                    y_block[y_index] += block_values_[elem_index] * x_block[x_index];
                }
            }
        }
    }
}

template<typename T>
void BlockCsrMatrix<T>::get_block_diag(std::vector<T>& diag_blocks) const {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法获取块对角元素");
    }
    
    diag_blocks.clear();
    diag_blocks.reserve(rows_ * block_dim_ * block_dim_);
    
    for (int block_row = 0; block_row < rows_; ++block_row) {
        int block_start = block_row_ptr_[block_row];
        int block_end = block_row_ptr_[block_row + 1];
        
        bool found_diag = false;
        for (int block_idx = block_start; block_idx < block_end; ++block_idx) {
            if (block_col_indices_[block_idx] == block_row) {
                // 找到对角块
                for (int i = 0; i < block_dim_; ++i) {
                    for (int j = 0; j < block_dim_; ++j) {
                        int elem_index = compute_element_index(block_idx, i, j);
                        diag_blocks.push_back(block_values_[elem_index]);
                    }
                }
                found_diag = true;
                break;
            }
        }
        
        if (!found_diag) {
            // 对角块为零，添加零块
            for (int i = 0; i < block_dim_ * block_dim_; ++i) {
                diag_blocks.push_back(T(0));
            }
        }
    }
}

template<typename T>
void BlockCsrMatrix<T>::set_block_diag(const std::vector<T>& diag_blocks) {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法设置块对角元素");
    }
    
    if (diag_blocks.size() != static_cast<size_t>(rows_ * block_dim_ * block_dim_)) {
        throw std::invalid_argument("对角块尺寸与矩阵尺寸不匹配");
    }
    
    for (int block_row = 0; block_row < rows_; ++block_row) {
        int block_start = block_row_ptr_[block_row];
        int block_end = block_row_ptr_[block_row + 1];
        
        for (int block_idx = block_start; block_idx < block_end; ++block_idx) {
            if (block_col_indices_[block_idx] == block_row) {
                // 设置对角块
                for (int i = 0; i < block_dim_; ++i) {
                    for (int j = 0; j < block_dim_; ++j) {
                        int elem_index = compute_element_index(block_idx, i, j);
                        int diag_index = block_row * block_dim_ * block_dim_ + i * block_dim_ + j;
                        block_values_[elem_index] = diag_blocks[diag_index];
                    }
                }
                break;
            }
        }
    }
}

template<typename T>
void BlockCsrMatrix<T>::scale(T alpha) {
    if (!built_) {
        throw std::runtime_error("矩阵未构建，无法执行缩放");
    }
    
    for (auto& value : block_values_) {
        value *= alpha;
    }
}

template<typename T>
const std::vector<int>& BlockCsrMatrix<T>::get_block_row_ptr() const {
    return block_row_ptr_;
}

template<typename T>
const std::vector<int>& BlockCsrMatrix<T>::get_block_col_indices() const {
    return block_col_indices_;
}

template<typename T>
const std::vector<T>& BlockCsrMatrix<T>::get_block_values() const {
    return block_values_;
}

template<typename T>
bool BlockCsrMatrix<T>::is_built() const {
    return built_;
}

template<typename T>
void BlockCsrMatrix<T>::validate_block_indices(int block_row, int block_col) const {
    if (block_row < 0 || block_row >= rows_) {
        throw std::out_of_range("块行索引越界");
    }
    if (block_col < 0 || block_col >= cols_) {
        throw std::out_of_range("块列索引越界");
    }
}

template<typename T>
int BlockCsrMatrix<T>::compute_element_index(int block_index, int i, int j) const {
    if (i < 0 || i >= block_dim_ || j < 0 || j >= block_dim_) {
        throw std::out_of_range("块内元素索引越界");
    }
    
    return block_index * block_dim_ * block_dim_ + i * block_dim_ + j;
}

} // namespace numeric