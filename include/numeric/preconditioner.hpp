#pragma once

#include "csr_matrix.hpp"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <memory>

// 前向声明
namespace emag {
    template<typename T>
    class Vector;
}

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
 * @brief 预处理器抽象基类
 * @details 定义预处理器的通用接口
 */
class Preconditioner {
public:
    /**
     * @brief 虚析构函数
     */
    virtual ~Preconditioner() = default;
    
    /**
     * @brief 应用预处理器（求解 Mz = r）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    virtual void apply(const emag::VectorReal& r, emag::VectorReal& z) const = 0;
    
    /**
     * @brief 应用预处理器（复数版本）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    virtual void apply(const emag::VectorComplex& r, emag::VectorComplex& z) const = 0;
    
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
 * @class JacobiPreconditioner
 * @brief Jacobi预处理器（对角预处理）
 * @details 使用矩阵对角元的倒数作为预处理矩阵
 */
class JacobiPreconditioner : public Preconditioner {
public:
    /**
     * @brief 构造函数
     * @param matrix 输入矩阵
     * @param epsilon 对角元最小值（避免除零）
     */
    JacobiPreconditioner(const CsrMatrixReal& matrix, double epsilon = 1e-12);
    
    /**
     * @brief 构造函数（复数版本）
     * @param matrix 输入矩阵
     * @param epsilon 对角元最小值（避免除零）
     */
    JacobiPreconditioner(const CsrMatrixComplex& matrix, double epsilon = 1e-12);
    
    /**
     * @brief 应用预处理器
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const emag::VectorReal& r, emag::VectorReal& z) const override;
    
    /**
     * @brief 应用预处理器（复数版本）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const emag::VectorComplex& r, emag::VectorComplex& z) const override;
    
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
    const emag::VectorReal& get_diag_inv_real() const { return diag_inv_real_; }
    
    /**
     * @brief 获取对角逆矩阵（复数版本）
     * @return 对角逆矩阵向量
     */
    const emag::VectorComplex& get_diag_inv_complex() const { return diag_inv_complex_; }

private:
    emag::VectorReal diag_inv_real_;      ///< 实数对角逆矩阵
    emag::VectorComplex diag_inv_complex_; ///< 复数对角逆矩阵
    bool is_complex_;                     ///< 是否为复数矩阵
    
    /**
     * @brief 构建Jacobi预处理器
     * @param matrix 输入矩阵
     * @param epsilon 对角元最小值
     */
    template<typename T>
    void build_jacobi(const CsrMatrix<T>& matrix, double epsilon);
};

/**
 * @class ILU0Preconditioner
 * @brief ILU(0)预处理器（无填充不完全LU分解）
 * @details 对矩阵进行不完全LU分解，不引入填充元
 */
class ILU0Preconditioner : public Preconditioner {
public:
    /**
     * @brief 构造函数
     * @param matrix 输入矩阵
     */
    ILU0Preconditioner(const CsrMatrixReal& matrix);
    
    /**
     * @brief 构造函数（复数版本）
     * @param matrix 输入矩阵
     */
    ILU0Preconditioner(const CsrMatrixComplex& matrix);
    
    /**
     * @brief 析构函数
     */
    ~ILU0Preconditioner() override = default;
    
    /**
     * @brief 应用预处理器（前向+后向替换）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const emag::VectorReal& r, emag::VectorReal& z) const override;
    
    /**
     * @brief 应用预处理器（复数版本）
     * @param r 输入向量（残差）
     * @param z 输出向量
     */
    void apply(const emag::VectorComplex& r, emag::VectorComplex& z) const override;
    
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
    const CsrMatrixReal& get_L_real() const { return L_real_; }
    
    /**
     * @brief 获取U矩阵
     * @return U矩阵（上三角）
     */
    const CsrMatrixReal& get_U_real() const { return U_real_; }
    
    /**
     * @brief 获取L矩阵（复数版本）
     * @return L矩阵（下三角）
     */
    const CsrMatrixComplex& get_L_complex() const { return L_complex_; }
    
    /**
     * @brief 获取U矩阵（复数版本）
     * @return U矩阵（上三角）
     */
    const CsrMatrixComplex& get_U_complex() const { return U_complex_; }

private:
    CsrMatrixReal L_real_;        ///< 实数L矩阵（下三角）
    CsrMatrixReal U_real_;        ///< 实数U矩阵（上三角）
    CsrMatrixComplex L_complex_;  ///< 复数L矩阵（下三角）
    CsrMatrixComplex U_complex_;  ///< 复数U矩阵（上三角）
    bool is_complex_;             ///< 是否为复数矩阵
    
    /**
     * @brief 构建ILU(0)分解
     * @param matrix 输入矩阵
     */
    template<typename T>
    void build_ilu0(const CsrMatrix<T>& matrix);
    
    /**
     * @brief 前向替换（求解 Lz = r）
     * @param L L矩阵
     * @param r 输入向量
     * @param z 输出向量
     */
    template<typename T>
    void forward_substitution(const CsrMatrix<T>& L, const emag::Vector<T>& r, emag::Vector<T>& z) const;
    
    /**
     * @brief 后向替换（求解 Ux = z）
     * @param U U矩阵
     * @param z 输入向量
     * @param x 输出向量
     */
    template<typename T>
    void backward_substitution(const CsrMatrix<T>& U, const emag::Vector<T>& z, emag::Vector<T>& x) const;
};

/**
 * @brief 创建预处理器
 * @param type 预处理器类型
 * @param matrix 输入矩阵
 * @return 预处理器指针
 */
std::unique_ptr<Preconditioner> create_preconditioner(PreconditionerType type, const CsrMatrixReal& matrix);

/**
 * @brief 创建预处理器（复数版本）
 * @param type 预处理器类型
 * @param matrix 输入矩阵
 * @return 预处理器指针
 */
std::unique_ptr<Preconditioner> create_preconditioner(PreconditionerType type, const CsrMatrixComplex& matrix);

} // namespace numeric