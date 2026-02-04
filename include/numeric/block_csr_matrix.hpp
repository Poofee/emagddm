/**
 * @file block_csr_matrix.hpp
 * @brief 块CSR格式稀疏矩阵类定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 * 
 * @details 适配低频电磁矢量元（Nedelec边元）的块稀疏矩阵，支持2×2/3×3块结构
 *          提升矢量元计算效率，避免单元素循环
 */

#pragma once

#include "sparse_base.hpp"
#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include <vector>
#include <complex>
#include <stdexcept>
#include <memory>

// 前向声明
namespace emag {
    template<typename T>
    class Vector;
}

namespace numeric {

/**
 * @brief 块大小枚举
 */
enum class BlockSize {
    BLOCK_1x1 = 1,  ///< 标量元（1×1块）
    BLOCK_2x2 = 2,  ///< 二维矢量元（2×2块）
    BLOCK_3x3 = 3   ///< 三维矢量元（3×3块）
};

/**
 * @class BlockCsrMatrix
 * @brief 块CSR格式稀疏矩阵类
 * @details 适配矢量元的块稀疏结构，提升计算效率
 */
template<typename T>
class BlockCsrMatrix : public SparseMatrixBase {
public:
    /**
     * @brief 默认构造函数
     */
    BlockCsrMatrix();
    
    /**
     * @brief 构造函数，指定矩阵尺寸和块大小
     * @param rows 矩阵行数（块行数）
     * @param cols 矩阵列数（块列数）
     * @param block_size 块大小
     */
    BlockCsrMatrix(int rows, int cols, BlockSize block_size);
    
    /**
     * @brief 析构函数
     */
    ~BlockCsrMatrix() override = default;
    
    /**
     * @brief 获取矩阵行数（块行数）
     * @return 矩阵行数
     */
    int rows() const override;
    
    /**
     * @brief 获取矩阵列数（块列数）
     * @return 矩阵列数
     */
    int cols() const override;
    
    /**
     * @brief 获取非零元素数量（块非零元数）
     * @return 非零元素数量
     */
    int nnz() const override;
    
    /**
     * @brief 获取块非零元数
     * @return 块非零元数
     */
    int block_nnz() const;
    
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
     * @brief 获取块大小
     * @return 块大小
     */
    BlockSize get_block_size() const;
    
    /**
     * @brief 获取块尺寸
     * @return 块尺寸（块大小×块大小）
     */
    int get_block_dim() const;
    
    /**
     * @brief 从COO格式构建块CSR矩阵
     * @param coo COO格式矩阵
     */
    void build_from_coo(const CooMatrix<T>& coo);
    
    /**
     * @brief 从普通CSR矩阵构建块CSR矩阵
     * @param csr 普通CSR矩阵
     * @param block_size 块大小
     */
    void build_from_csr(const CsrMatrix<T>& csr, BlockSize block_size);
    
    /**
     * @brief 矩阵向量乘法（块级）
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const emag::Vector<T>& x, emag::Vector<T>& y) const;
    
    /**
     * @brief 矩阵向量乘法（块级，使用std::vector）
     * @param x 输入向量
     * @param y 输出向量
     */
    void mat_vec(const std::vector<T>& x, std::vector<T>& y) const;
    
    /**
     * @brief 块级矩阵向量乘法
     * @param x_block 输入块向量
     * @param y_block 输出块向量
     */
    void block_mat_vec(const std::vector<T>& x_block, std::vector<T>& y_block) const;
    
    /**
     * @brief 获取块对角元素
     * @param diag_blocks 块对角元素向量
     */
    void get_block_diag(std::vector<T>& diag_blocks) const;
    
    /**
     * @brief 设置块对角元素
     * @param diag_blocks 块对角元素向量
     */
    void set_block_diag(const std::vector<T>& diag_blocks);
    
    /**
     * @brief 矩阵数乘（块级）
     * @param alpha 缩放因子
     */
    void scale(T alpha);
    
    /**
     * @brief 获取块行偏移数组
     * @return 块行偏移数组引用
     */
    const std::vector<int>& get_block_row_ptr() const;
    
    /**
     * @brief 获取块列索引数组
     * @return 块列索引数组引用
     */
    const std::vector<int>& get_block_col_indices() const;
    
    /**
     * @brief 获取块数值数组
     * @return 块数值数组引用
     */
    const std::vector<T>& get_block_values() const;
    
    /**
     * @brief 检查矩阵是否已构建
     * @return 是否已构建
     */
    bool is_built() const;

private:
    int rows_;                          ///< 矩阵行数（块行数）
    int cols_;                          ///< 矩阵列数（块列数）
    BlockSize block_size_;              ///< 块大小
    int block_dim_;                     ///< 块尺寸（块大小×块大小）
    
    std::vector<int> block_row_ptr_;    ///< 块行偏移数组
    std::vector<int> block_col_indices_; ///< 块列索引数组
    std::vector<T> block_values_;       ///< 块数值数组
    
    bool built_;                        ///< 矩阵是否已构建
    
    /**
     * @brief 验证块索引范围
     * @param block_row 块行索引
     * @param block_col 块列索引
     */
    void validate_block_indices(int block_row, int block_col) const;
    
    /**
     * @brief 计算块内元素索引
     * @param block_index 块索引
     * @param i 块内行索引
     * @param j 块内列索引
     * @return 元素在块数值数组中的索引
     */
    int compute_element_index(int block_index, int i, int j) const;
};

// 常用类型别名
using BlockCsrMatrixReal = BlockCsrMatrix<double>;
using BlockCsrMatrixComplex = BlockCsrMatrix<std::complex<double>>;

} // namespace numeric