/**
 * @file preconditioner_impl.hpp
 * @brief 预处理器模板实现文件
 * @author Poofee
 * @date 2026-02-06
 * @version 2.0
 */

#pragma once

#include "preconditioner.hpp"
#include "coo_matrix.hpp"
#include <iostream>
#include <cmath>

namespace numeric {

template<typename T>
JacobiPreconditioner<T>::JacobiPreconditioner(const CsrMatrix<T>& matrix, double epsilon)
    : epsilon_(epsilon) {
    std::vector<T> diag;
    matrix.get_diag(diag);

    int n = static_cast<int>(diag.size());
    diag_inv_.resize(n);
    for (int i = 0; i < n; ++i) {
        T abs_val = std::abs(diag[i]) < epsilon ? T(epsilon) : diag[i];
        diag_inv_[i] = T(1.0) / abs_val;
    }
}

template<typename T>
void JacobiPreconditioner<T>::apply(const vector_type& r, vector_type& z) const {
    if (r.size() != z.size()) {
        z.resize(r.size());
    }

    int n = static_cast<int>(r.size());
    for (int i = 0; i < n; ++i) {
        z[i] = r[i] * diag_inv_[i];
    }
}

template<typename T>
void JacobiPreconditioner<T>::print_info() const {
    std::cout << "Jacobi Preconditioner" << std::endl;
    std::cout << "  Diagonal inverse size: " << diag_inv_.size() << std::endl;
}

template<typename T>
ILU0Preconditioner<T>::ILU0Preconditioner(const CsrMatrix<T>& matrix) {
    build_ilu0(matrix);
}

template<typename T>
void ILU0Preconditioner<T>::build_ilu0(const CsrMatrix<T>& matrix) {
    int n = matrix.rows();

    CooMatrix<T> L_coo(n, n);
    CooMatrix<T> U_coo(n, n);
    L_coo.reserve(n * n);
    U_coo.reserve(n * n);

    const auto& row_ptr = matrix.get_row_ptr();
    const auto& col_idx = matrix.get_col_indices();
    const auto& values = matrix.get_values();

    std::vector<std::vector<int>> U_cols(n);
    std::vector<std::vector<T>> U_vals(n);

    for (int i = 0; i < n; ++i) {
        std::vector<T> row_val(n, T(0));
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            row_val[col_idx[j]] = values[j];
        }

        for (int k = 0; k < i; ++k) {
            double mag = std::abs(row_val[k]);
            if (mag < 1e-15) continue;
            T factor = row_val[k] / U_vals[k][0];
            row_val[k] = factor;
            for (size_t p = 1; p < U_vals[k].size(); ++p) {
                int col = U_cols[k][p];
                if (col > k) row_val[col] -= factor * U_vals[k][p];
            }
        }

        L_coo.add_value(i, i, T(1));
        for (int k = 0; k < i; ++k) {
            if (std::abs(row_val[k]) > 1e-15) {
                L_coo.add_value(i, k, row_val[k]);
            }
        }

        U_cols[i].clear();
        U_vals[i].clear();
        for (int k = i; k < n; ++k) {
            if (std::abs(row_val[k]) > 1e-15) {
                U_cols[i].push_back(k);
                U_vals[i].push_back(row_val[k]);
                U_coo.add_value(i, k, row_val[k]);
            }
        }
    }

    L_ = CsrMatrix<T>(n, n);
    U_ = CsrMatrix<T>(n, n);
    L_.build_from_coo(L_coo);
    U_.build_from_coo(U_coo);
}

template<typename T>
void ILU0Preconditioner<T>::forward_substitution(const CsrMatrix<T>& L,
                                                   const vector_type& r,
                                                   vector_type& z) const {
    int n = static_cast<int>(r.size());
    z.resize(n);

    const auto& row_ptr = L.get_row_ptr();
    const auto& col_idx = L.get_col_indices();
    const auto& values = L.get_values();

    for (int i = 0; i < n; ++i) {
        T sum = r[i];
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = col_idx[j];
            if (col < i) {
                sum -= values[j] * z[col];
            }
        }
        z[i] = sum;
    }
}

template<typename T>
void ILU0Preconditioner<T>::backward_substitution(const CsrMatrix<T>& U,
                                                    const vector_type& z,
                                                    vector_type& x) const {
    int n = static_cast<int>(z.size());
    x.resize(n);

    const auto& row_ptr = U.get_row_ptr();
    const auto& col_idx = U.get_col_indices();
    const auto& values = U.get_values();

    for (int i = n - 1; i >= 0; --i) {
        T sum = z[i];
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = col_idx[j];
            if (col > i) {
                sum -= values[j] * x[col];
            }
        }
        x[i] = sum / values[row_ptr[i]];
    }
}

template<typename T>
void ILU0Preconditioner<T>::apply(const vector_type& r, vector_type& z) const {
    vector_type temp;
    forward_substitution(L_, r, temp);
    backward_substitution(U_, temp, z);
}

template<typename T>
void ILU0Preconditioner<T>::print_info() const {
    std::cout << "ILU(0) Preconditioner" << std::endl;
    std::cout << "  L matrix: " << L_.rows() << "x" << L_.cols() << std::endl;
    std::cout << "  U matrix: " << U_.rows() << "x" << U_.cols() << std::endl;
}

template<typename T>
std::unique_ptr<Preconditioner<T>> create_preconditioner(PreconditionerType type, 
                                                            const CsrMatrix<T>& matrix) {
    switch (type) {
        case PreconditionerType::JACOBI:
            return std::make_unique<JacobiPreconditioner<T>>(matrix);
        case PreconditionerType::ILU0:
            return std::make_unique<ILU0Preconditioner<T>>(matrix);
        default:
            return std::make_unique<JacobiPreconditioner<T>>(matrix);
    }
}

template class JacobiPreconditioner<double>;
template class JacobiPreconditioner<std::complex<double>>;
template class ILU0Preconditioner<double>;
template class ILU0Preconditioner<std::complex<double>>;

} // namespace numeric
