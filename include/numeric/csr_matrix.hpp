/**
 * @file csr_matrix.hpp
 * @brief CSR格式稀疏矩阵类定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "sparse_base.hpp"
#include <vector>
#include <complex>

namespace numeric {

/**
 * @class CsrMatrix
 * @brief CSR格式稀疏矩阵类
 * @details 使用压缩行存储格式，适合矩阵求解阶段
 */
template<typename T>
class CsrMatrix : public SparseMatrixBase {
public:
    /**
     * @brief 默认构造函数
     */
    CsrMatrix();

    /**
     * @brief 构造函数，指定矩阵尺寸
     * @param rows 矩阵行数
     * @param cols 矩阵列数
     */
    CsrMatrix(int rows, int cols);

    /**
     * @brief 析构函数
     */
    ~CsrMatrix() override = default;

    /**
     * @brief 获取矩阵行数
     * @return 矩阵行数
     */
    int rows() const override;

    /**
     * @brief 获取矩阵列数
     * @return 矩阵列数
     */
    int cols() const override;

    /**
     * @brief 获取非零元素数量
     * @return 非零元素数量
     */
    int nnz() const override;

    /**
     * @brief 清空矩阵
     */
    void clear() override;

    /**
     * @brief 打印矩阵信息
     */
    void print_info() const override;

    /**
     * @brief 获取矩阵数据类型
     * @return 矩阵数据类型
     */
    MatrixDataType get_data_type() const override;

    /**
     * @brief 从COO格式构建CSR矩阵
     * @param coo COO格式矩阵
     */
    void build_from_coo(const CooMatrix<T>& coo);

    /**
     * @brief 矩阵向量乘法
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const std::vector<T>& x, std::vector<T>& y) const;

    /**
     * @brief 获取行偏移数组
     * @return 行偏移数组引用
     */
    const std::vector<int>& get_row_ptr() const;

    /**
     * @brief 获取列索引数组
     * @return 列索引数组引用
     */
    const std::vector<int>& get_col_indices() const;

    /**
     * @brief 获取数值数组
     * @return 数值数组引用
     */
    const std::vector<T>& get_values() const;

    /**
     * @brief 矩阵数乘
     * @param alpha 缩放因子
     */
    void scale(T alpha);

    /**
     * @brief 获取矩阵转置
     * @return 转置后的CSR矩阵
     */
    CsrMatrix<T> transpose() const;

    /**
     * @brief 提取对角线元素
     * @param diag 对角线元素向量
     */
    void get_diag(std::vector<T>& diag) const;

    /**
     * @brief 设置对角线元素
     * @param diag 对角线元素向量
     */
    void set_diag(const std::vector<T>& diag);

    /**
     * @brief 检查矩阵是否已构建
     * @return 是否已构建
     */
    bool is_built() const;

private:
    int rows_;                    ///< 矩阵行数
    int cols_;                    ///< 矩阵列数
    std::vector<int> row_ptr_;    ///< 行偏移数组
    std::vector<int> col_indices_; ///< 列索引数组
    std::vector<T> values_;        ///< 数值数组
    bool built_;                   ///< 矩阵是否已构建

    /**
     * @brief 对COO数据进行排序（按行主序，行内按列排序）
     * @param rows 行索引数组
     * @param cols 列索引数组
     * @param values 数值数组
     * @param nnz 非零元素数量
     */
    void sort_coo_data(std::vector<int>& rows, std::vector<int>& cols, 
                      std::vector<T>& values, int nnz) const;
};

// 常用类型别名
using CsrMatrixReal = CsrMatrix<double>;
using CsrMatrixComplex = CsrMatrix<std::complex<double>>;

} // namespace numeric

// 包含模板实现
#include "csr_matrix_impl.hpp"