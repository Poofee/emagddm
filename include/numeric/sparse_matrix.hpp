/**
 * @file sparse_matrix.hpp
 * @brief 核心数值层 - 稀疏矩阵管理模块头文件
 * @details 提供稀疏矩阵的基本操作接口，基于新的稀疏矩阵架构
 * @author Poofee
 * @date 2026-02-04
 * @version 2.0
 * 
 * @note 此文件为兼容性接口，推荐使用新的稀疏矩阵类：
 * - CooMatrix: COO格式稀疏矩阵，适合矩阵组装
 * - CsrMatrix: CSR格式稀疏矩阵，适合矩阵求解
 * - MatrixMarketIO: MatrixMarket格式I/O功能
 */

#pragma once

#include "sparse_base.hpp"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "matrix_market_io.hpp"

namespace numeric {

/**
 * @class SparseMatrix
 * @brief 稀疏矩阵管理类（兼容性接口）
 * @deprecated 推荐使用新的CooMatrix或CsrMatrix类
 */
class SparseMatrix {
public:
    /**
     * @brief 构造函数
     */
    SparseMatrix();

    /**
     * @brief 析构函数
     */
    ~SparseMatrix();

    /**
     * @brief 设置矩阵尺寸
     * @param rows 行数
     * @param cols 列数
     */
    void setSize(int rows, int cols);

    /**
     * @brief 插入矩阵元素
     * @param row 行索引
     * @param col 列索引
     * @param value 元素值
     */
    void insert(int row, int col, double value);

    /**
     * @brief 获取矩阵行数
     * @return int 行数
     */
    int rows() const;

    /**
     * @brief 获取矩阵列数
     * @return int 列数
     */
    int cols() const;

private:
    int rows_;
    int cols_;
    std::vector<double> values_;
};

} // namespace numeric