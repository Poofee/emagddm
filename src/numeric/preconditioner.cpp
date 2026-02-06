/**
 * @file preconditioner.cpp
 * @brief 预处理器源文件
 * @author Poofee
 * @date 2026-02-06
 * @version 2.0
 */

#include "numeric/preconditioner.hpp"

namespace numeric {

// 模板类的显式实例化
template class JacobiPreconditioner<double>;
template class JacobiPreconditioner<std::complex<double>>;
template class ILU0Preconditioner<double>;
template class ILU0Preconditioner<std::complex<double>>;

} // namespace numeric
