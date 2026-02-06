/**
 * @file preconditioner_impl.hpp
 * @brief 预处理器模板实现文件
 * @author Poofee
 * @date 2026-02-06
 * @version 2.0
 */

#pragma once

#include "preconditioner.hpp"
#include <iostream>

namespace numeric {

template<typename T>
JacobiPreconditioner<T>::JacobiPreconditioner(const CsrMatrix<T>& matrix, double epsilon)
    : epsilon_(epsilon) {
}

template<typename T>
void JacobiPreconditioner<T>::apply(const vector_type& r, vector_type& z) const {
    if ((int)r.size() != (int)z.size()) {
        z.resize(r.size());
    }
    
    int n = (int)r.size();
    for (int i = 0; i < n; ++i) {
        z[i] = r[i];
    }
}

template<typename T>
void JacobiPreconditioner<T>::print_info() const {
    std::cout << "Jacobi Preconditioner" << std::endl;
    std::cout << "  Diagonal inverse size: " << diag_inv_.size() << std::endl;
}

template<typename T>
ILU0Preconditioner<T>::ILU0Preconditioner(const CsrMatrix<T>& matrix) {
}

template<typename T>
void ILU0Preconditioner<T>::build_ilu0(const CsrMatrix<T>& matrix) {
    int n = matrix.rows();
    L_ = CsrMatrix<T>(n, n);
    U_ = CsrMatrix<T>(n, n);
}

template<typename T>
void ILU0Preconditioner<T>::forward_substitution(const CsrMatrix<T>& L, 
                                                   const vector_type& r, 
                                                   vector_type& z) const {
    int n = (int)r.size();
    z.resize(n);
    
    for (int i = 0; i < n; ++i) {
        z[i] = r[i];
    }
}

template<typename T>
void ILU0Preconditioner<T>::backward_substitution(const CsrMatrix<T>& U, 
                                                    const vector_type& z, 
                                                    vector_type& x) const {
    int n = (int)z.size();
    x.resize(n);
    
    for (int i = 0; i < n; ++i) {
        x[i] = z[i];
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
