/**
 * @file preconditioner.hpp
 * @brief 预处理器模块头文件
 * @author Poofee
 * @date 2026-02-06
 * @version 2.0
 */

#pragma once

#include "csr_matrix.hpp"
#include "vector.hpp"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <memory>

namespace numeric {

/**
 * @brief 预处理器类型枚举
 */
enum class PreconditionerType {
    JACOBI,     ///< Jacobi预处理（对角预处理）
    ILU0        ///< ILU(0)预处理（无填充不完全LU分解）
};

/**
 * @class Preconditioner
 * @brief 预处理器抽象基类（模板版本）
 * @details 定义预处理器的通用接口，支持实数和复数类型
 */
template<typename T>
class Preconditioner {
public:
    using value_type = T;
    using vector_type = emag::Vector<T>;
    
    /**
     * @brief 虚析构函数
     */
    virtual ~Preconditioner() = default;
    
    /**
     * @brief 应用预处理器（求解 Mz = r）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    virtual void apply(const vector_type& r, vector_type& z) const = 0;
    
    /**
     * @brief 获取预处理器类型
     * @return 预处理器类型
     */
    virtual PreconditionerType get_type() const = 0;
    
    /**
     * @brief 打印预处理器信息
     */
    virtual void print_info() const = 0;
};

/**
 * @brief Jacobi预处理器别名（实数）
 */
using PreconditionerReal = Preconditioner<double>;

/**
 * @brief Jacobi预处理器别名（复数）
 */
using PreconditionerComplex = Preconditioner<std::complex<double>>;

/**
 * @class JacobiPreconditioner
 * @brief Jacobi预处理器（对角预处理，模板版本）
 * @details 使用矩阵对角元的倒数作为预处理矩阵
 */
template<typename T>
class JacobiPreconditioner : public Preconditioner<T> {
public:
    using value_type = T;
    using vector_type = emag::Vector<T>;
    
    /**
     * @brief 构造函数
     * @param matrix 输入矩阵
     * @param epsilon 对角元最小值（避免除零）
     */
    JacobiPreconditioner(const CsrMatrix<T>& matrix, double epsilon = 1e-12);
    
    /**
     * @brief 应用预处理器
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const vector_type& r, vector_type& z) const override;
    
    /**
     * @brief 获取预处理器类型
     * @return 预处理器类型
     */
    PreconditionerType get_type() const override { return PreconditionerType::JACOBI; }
    
    /**
     * @brief 打印预处理器信息
     */
    void print_info() const override;
    
    /**
     * @brief 获取对角逆矩阵
     * @return 对角逆矩阵向量
     */
    const vector_type& get_diag_inv() const { return diag_inv_; }

private:
    vector_type diag_inv_;  ///< 对角逆矩阵
    double epsilon_;        ///< 避免除零的小量
};

/**
 * @class ILU0Preconditioner
 * @brief ILU(0)预处理器（无填充不完全LU分解，模板版本）
 * @details 对矩阵进行不完全LU分解，不引入填充元
 */
template<typename T>
class ILU0Preconditioner : public Preconditioner<T> {
public:
    using value_type = T;
    using vector_type = emag::Vector<T>;
    
    /**
     * @brief 构造函数
     * @param matrix 输入矩阵
     */
    ILU0Preconditioner(const CsrMatrix<T>& matrix);
    
    /**
     * @brief 析构函数
     */
    ~ILU0Preconditioner() override = default;
    
    /**
     * @brief 应用预处理器（前向+后向替换）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const vector_type& r, vector_type& z) const override;
    
    /**
     * @brief 获取预处理器类型
     * @return 预处理器类型
     */
    PreconditionerType get_type() const override { return PreconditionerType::ILU0; }
    
    /**
     * @brief 打印预处理器信息
     */
    void print_info() const override;
    
    /**
     * @brief 获取L矩阵
     * @return L矩阵（下三角）
     */
    const CsrMatrix<T>& get_L() const { return L_; }
    
    /**
     * @brief 获取U矩阵
     * @return U矩阵（上三角）
     */
    const CsrMatrix<T>& get_U() const { return U_; }

private:
    CsrMatrix<T> L_;  ///< L矩阵（下三角）
    CsrMatrix<T> U_;  ///< U矩阵（上三角）
    
    /**
     * @brief 构建ILU(0)分解
     * @param matrix 输入矩阵
     */
    void build_ilu0(const CsrMatrix<T>& matrix);
    
    /**
     * @brief 前向替换（求解 Lz = r）
     * @param L L矩阵
     * @param r 输入向量
     * @param z 输出向量
     */
    void forward_substitution(const CsrMatrix<T>& L, const vector_type& r, vector_type& z) const;
    
    /**
     * @brief 后向替换（求解 Ux = z）
     * @param U U矩阵
     * @param z 输入向量
     * @param x 输出向量
     */
    void backward_substitution(const CsrMatrix<T>& U, const vector_type& z, vector_type& x) const;
};

// 类型别名
using JacobiPreconditionerReal = JacobiPreconditioner<double>;
using JacobiPreconditionerComplex = JacobiPreconditioner<std::complex<double>>;
using ILU0PreconditionerReal = ILU0Preconditioner<double>;
using ILU0PreconditionerComplex = ILU0Preconditioner<std::complex<double>>;

/**
 * @brief 创建预处理器
 * @param type 预处理器类型
 * @param matrix 输入矩阵
 * @return 预处理器指针
 */
template<typename T>
std::unique_ptr<Preconditioner<T>> create_preconditioner(PreconditionerType type, const CsrMatrix<T>& matrix);

// 显式实例化声明
extern template class JacobiPreconditioner<double>;
extern template class JacobiPreconditioner<std::complex<double>>;
extern template class ILU0Preconditioner<double>;
extern template class ILU0Preconditioner<std::complex<double>>;

} // namespace numeric

// 模板实现
#include "preconditioner_impl.hpp"
