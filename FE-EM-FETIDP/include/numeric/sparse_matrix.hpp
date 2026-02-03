/**
 * @file sparse_matrix.hpp
 * @brief 核心数值层 - 稀疏矩阵管理模块头文件
 * @details 提供稀疏矩阵的基本操作接口，封装Eigen库的稀疏矩阵功能
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include <vector>

namespace numeric {

/**
 * @class SparseMatrix
 * @brief 稀疏矩阵管理类
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