#pragma once

#include "sparse_base.hpp"
#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include <vector>
#include <complex>
#include <stdexcept>

// 前向声明
namespace emag {
    template<typename T>
    class Vector;
}

namespace numeric {

/**
 * @class SymCsrMatrix
 * @brief 对称CSR格式稀疏矩阵类
 * @details 仅存储下三角部分，内存占用减半，适合对称矩阵
 */
template<typename T>
class SymCsrMatrix : public SparseMatrixBase {
public:
    /**
     * @brief 默认构造函数
     */
    SymCsrMatrix();
    
    /**
     * @brief 构造函数，指定矩阵尺寸
     * @param size 矩阵尺寸（方阵）
     */
    explicit SymCsrMatrix(int size);
    
    /**
     * @brief 析构函数
     */
    ~SymCsrMatrix() override = default;
    
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
     * @brief 从COO格式构建对称CSR矩阵
     * @param coo COO格式矩阵（仅需下三角部分）
     */
    void build_from_coo(const CooMatrix<T>& coo);
    
    /**
     * @brief 从完整CSR矩阵构建对称CSR矩阵
     * @param csr 完整CSR矩阵
     */
    void build_from_full_csr(const CsrMatrix<T>& csr);
    
    /**
     * @brief 矩阵向量乘法
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const emag::Vector<T>& x, emag::Vector<T>& y) const;
    
    /**
     * @brief 矩阵向量乘法（使用std::vector）
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const std::vector<T>& x, std::vector<T>& y) const;
    
    /**
     * @brief 矩阵数乘
     * @param alpha 缩放因子
     */
    void scale(T alpha);
    
    /**
     * @brief 获取矩阵转置（对称矩阵转置等于自身）
     * @return 转置后的对称CSR矩阵（返回自身引用）
     */
    const SymCsrMatrix<T>& transpose() const;
    
    /**
     * @brief 提取对角线元素
     * @param diag 对角线元素向量
     */
    void get_diag(emag::Vector<T>& diag) const;
    
    /**
     * @brief 设置对角线元素
     * @param diag 对角线元素向量
     */
    void set_diag(const emag::Vector<T>& diag);
    
    /**
     * @brief 检查矩阵是否已构建
     * @return 是否已构建
     */
    bool is_built() const;
    
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
     * @brief 转换为完整CSR矩阵
     * @return 完整CSR矩阵
     */
    CsrMatrix<T> to_full_csr() const;

private:
    int size_;                    ///< 矩阵尺寸（方阵）
    std::vector<int> row_ptr_;    ///< 行偏移数组
    std::vector<int> col_indices_; ///< 列索引数组（仅下三角）
    std::vector<T> values_;        ///< 数值数组（仅下三角）
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
    
    /**
     * @brief 验证矩阵对称性
     * @param coo COO格式矩阵
     */
    void validate_symmetry(const CooMatrix<T>& coo) const;
};

// 常用类型别名
using SymCsrMatrixReal = SymCsrMatrix<double>;
using SymCsrMatrixComplex = SymCsrMatrix<std::complex<double>>;

} // namespace numeric