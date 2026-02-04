#pragma once

#include "preconditioner.hpp"
#include "coo_matrix.hpp"
#include "vector.hpp"
#include <memory>
#include <iostream>
#include <cmath>

namespace numeric {

// JacobiPreconditioner实现

JacobiPreconditioner::JacobiPreconditioner(const CsrMatrixReal& matrix, double epsilon) 
    : is_complex_(false) {
    build_jacobi(matrix, epsilon);
}

JacobiPreconditioner::JacobiPreconditioner(const CsrMatrixComplex& matrix, double epsilon) 
    : is_complex_(true) {
    build_jacobi(matrix, epsilon);
}

template<typename T>
void JacobiPreconditioner::build_jacobi(const CsrMatrix<T>& matrix, double epsilon) {
    if (!matrix.is_built()) {
        throw std::runtime_error("输入矩阵未构建，无法构建Jacobi预处理器");
    }
    
    int size = matrix.rows();
    if (size != matrix.cols()) {
        throw std::invalid_argument("Jacobi预处理器仅适用于方阵");
    }
    
    // 提取对角线元素
    std::vector<T> diag;
    matrix.get_diag(diag);
    
    // 构建对角逆矩阵
    if constexpr (std::is_same_v<T, double>) {
        diag_inv_real_.resize(size);
        for (int i = 0; i < size; ++i) {
            double diag_val = diag[i];
            if (std::abs(diag_val) < epsilon) {
                // 避免除零，使用微小值
                diag_val = (diag_val >= 0) ? epsilon : -epsilon;
            }
            diag_inv_real_[i] = 1.0 / diag_val;
        }
    } else {
        diag_inv_complex_.resize(size);
        for (int i = 0; i < size; ++i) {
            std::complex<double> diag_val = diag[i];
            double magnitude = std::abs(diag_val);
            if (magnitude < epsilon) {
                // 避免除零，使用微小值
                diag_val = std::complex<double>(epsilon, 0.0);
            }
            diag_inv_complex_[i] = 1.0 / diag_val;
        }
    }
}

void JacobiPreconditioner::apply(const emag::VectorReal& r, emag::VectorReal& z) const {
    if (is_complex_) {
        throw std::runtime_error("Jacobi预处理器为复数版本，无法应用于实数向量");
    }
    
    if (r.size() != diag_inv_real_.size()) {
        throw std::invalid_argument("输入向量尺寸与预处理器尺寸不匹配");
    }
    
    z.resize(r.size());
    for (int i = 0; i < r.size(); ++i) {
        z[i] = diag_inv_real_[i] * r[i];
    }
}

void JacobiPreconditioner::apply(const emag::VectorComplex& r, emag::VectorComplex& z) const {
    if (!is_complex_) {
        throw std::runtime_error("Jacobi预处理器为实数版本，无法应用于复数向量");
    }
    
    if (r.size() != diag_inv_complex_.size()) {
        throw std::invalid_argument("输入向量尺寸与预处理器尺寸不匹配");
    }
    
    z.resize(r.size());
    for (int i = 0; i < r.size(); ++i) {
        z[i] = diag_inv_complex_[i] * r[i];
    }
}

void JacobiPreconditioner::print_info() const {
    std::cout << "Jacobi预处理器信息:" << std::endl;
    std::cout << "  类型: " << (is_complex_ ? "复数" : "实数") << std::endl;
    std::cout << "  尺寸: " << (is_complex_ ? diag_inv_complex_.size() : diag_inv_real_.size()) << std::endl;
}

// ILU0Preconditioner实现

ILU0Preconditioner::ILU0Preconditioner(const CsrMatrixReal& matrix) 
    : is_complex_(false) {
    build_ilu0(matrix);
}

ILU0Preconditioner::ILU0Preconditioner(const CsrMatrixComplex& matrix) 
    : is_complex_(true) {
    build_ilu0(matrix);
}

template<typename T>
void ILU0Preconditioner::build_ilu0(const CsrMatrix<T>& matrix) {
    if (!matrix.is_built()) {
        throw std::runtime_error("输入矩阵未构建，无法构建ILU(0)预处理器");
    }
    
    int size = matrix.rows();
    if (size != matrix.cols()) {
        throw std::invalid_argument("ILU(0)预处理器仅适用于方阵");
    }
    
    // 获取矩阵数据
    const auto& row_ptr = matrix.get_row_ptr();
    const auto& col_indices = matrix.get_col_indices();
    const auto& values = matrix.get_values();
    
    // 创建L和U矩阵的COO格式
    CooMatrix<T> L_coo(size, size);
    CooMatrix<T> U_coo(size, size);
    
    // 为L和U矩阵预留空间
    L_coo.reserve(matrix.nnz());
    U_coo.reserve(matrix.nnz());
    
    // 临时存储用于ILU分解
    std::vector<T> diag(size, T(0));
    
    // ILU(0)分解算法
    for (int i = 0; i < size; ++i) {
        // 处理第i行
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = col_indices[j];
            T value = values[j];
            
            if (col < i) {
                // 下三角部分，属于L矩阵
                L_coo.add_value(i, col, value);
            } else if (col == i) {
                // 对角线元素
                diag[i] = value;
                // 设置L和U的对角线为1
                L_coo.add_value(i, i, T(1));
                U_coo.add_value(i, i, value);
            } else {
                // 上三角部分，属于U矩阵
                U_coo.add_value(i, col, value);
            }
        }
        
        // 简单的ILU(0)处理：保持原矩阵的非零模式
        // 实际ILU分解需要更复杂的算法，这里简化实现
    }
    
    // 构建L和U矩阵
    if constexpr (std::is_same_v<T, double>) {
        L_real_.build_from_coo(L_coo);
        U_real_.build_from_coo(U_coo);
    } else {
        L_complex_.build_from_coo(L_coo);
        U_complex_.build_from_coo(U_coo);
    }
}

template<typename T>
void ILU0Preconditioner::forward_substitution(const CsrMatrix<T>& L, const emag::Vector<T>& r, emag::Vector<T>& z) const {
    int size = r.size();
    z.resize(size);
    
    const auto& L_row_ptr = L.get_row_ptr();
    const auto& L_col_indices = L.get_col_indices();
    const auto& L_values = L.get_values();
    
    for (int i = 0; i < size; ++i) {
        T sum = r[i];
        for (int j = L_row_ptr[i]; j < L_row_ptr[i + 1]; ++j) {
            int col = L_col_indices[j];
            if (col < i) {
                sum -= L_values[j] * z[col];
            }
        }
        // L的对角线为1
        z[i] = sum;
    }
}

template<typename T>
void ILU0Preconditioner::backward_substitution(const CsrMatrix<T>& U, const emag::Vector<T>& z, emag::Vector<T>& x) const {
    int size = z.size();
    x.resize(size);
    
    const auto& U_row_ptr = U.get_row_ptr();
    const auto& U_col_indices = U.get_col_indices();
    const auto& U_values = U.get_values();
    
    for (int i = size - 1; i >= 0; --i) {
        T sum = z[i];
        for (int j = U_row_ptr[i]; j < U_row_ptr[i + 1]; ++j) {
            int col = U_col_indices[j];
            if (col > i) {
                sum -= U_values[j] * x[col];
            }
        }
        // 除以对角线元素
        for (int j = U_row_ptr[i]; j < U_row_ptr[i + 1]; ++j) {
            if (U_col_indices[j] == i) {
                x[i] = sum / U_values[j];
                break;
            }
        }
    }
}

void ILU0Preconditioner::apply(const emag::VectorReal& r, emag::VectorReal& z) const {
    if (is_complex_) {
        throw std::runtime_error("ILU(0)预处理器为复数版本，无法应用于实数向量");
    }
    
    emag::VectorReal y;
    forward_substitution(L_real_, r, y);
    backward_substitution(U_real_, y, z);
}

void ILU0Preconditioner::apply(const emag::VectorComplex& r, emag::VectorComplex& z) const {
    if (!is_complex_) {
        throw std::runtime_error("ILU(0)预处理器为实数版本，无法应用于复数向量");
    }
    
    emag::VectorComplex y;
    forward_substitution(L_complex_, r, y);
    backward_substitution(U_complex_, y, z);
}

void ILU0Preconditioner::print_info() const {
    std::cout << "ILU(0)预处理器信息:" << std::endl;
    std::cout << "  类型: " << (is_complex_ ? "复数" : "实数") << std::endl;
    if (is_complex_) {
        std::cout << "  L矩阵非零元: " << L_complex_.nnz() << std::endl;
        std::cout << "  U矩阵非零元: " << U_complex_.nnz() << std::endl;
    } else {
        std::cout << "  L矩阵非零元: " << L_real_.nnz() << std::endl;
        std::cout << "  U矩阵非零元: " << U_real_.nnz() << std::endl;
    }
}

// 预处理器工厂函数

std::unique_ptr<Preconditioner> create_preconditioner(PreconditionerType type, const CsrMatrixReal& matrix) {
    switch (type) {
        case PreconditionerType::JACOBI:
            return std::make_unique<JacobiPreconditioner>(matrix);
        case PreconditionerType::ILU0:
            return std::make_unique<ILU0Preconditioner>(matrix);
        default:
            throw std::invalid_argument("不支持的预处理器类型");
    }
}

std::unique_ptr<Preconditioner> create_preconditioner(PreconditionerType type, const CsrMatrixComplex& matrix) {
    switch (type) {
        case PreconditionerType::JACOBI:
            return std::make_unique<JacobiPreconditioner>(matrix);
        case PreconditionerType::ILU0:
            return std::make_unique<ILU0Preconditioner>(matrix);
        default:
            throw std::invalid_argument("不支持的预处理器类型");
    }
}

// 显式实例化模板函数
template void JacobiPreconditioner::build_jacobi<double>(const CsrMatrix<double>&, double);
template void JacobiPreconditioner::build_jacobi<std::complex<double>>(const CsrMatrix<std::complex<double>>&, double);

template void ILU0Preconditioner::build_ilu0<double>(const CsrMatrix<double>&);
template void ILU0Preconditioner::build_ilu0<std::complex<double>>(const CsrMatrix<std::complex<double>>&);

template void ILU0Preconditioner::forward_substitution<double>(const CsrMatrix<double>&, const emag::Vector<double>&, emag::Vector<double>&) const;
template void ILU0Preconditioner::forward_substitution<std::complex<double>>(const CsrMatrix<std::complex<double>>&, const emag::Vector<std::complex<double>>&, emag::Vector<std::complex<double>>&) const;

template void ILU0Preconditioner::backward_substitution<double>(const CsrMatrix<double>&, const emag::Vector<double>&, emag::Vector<double>&) const;
template void ILU0Preconditioner::backward_substitution<std::complex<double>>(const CsrMatrix<std::complex<double>>&, const emag::Vector<std::complex<double>>&, emag::Vector<std::complex<double>>&) const;

} // namespace numeric