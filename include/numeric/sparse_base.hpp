/**
 * @file sparse_base.hpp
 * @brief 稀疏矩阵抽象基类定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include <vector>
#include <complex>

namespace numeric {

/**
 * @brief 矩阵数据类型枚举
 */
enum class MatrixDataType {
    REAL,    ///< 实数矩阵
    COMPLEX  ///< 复数矩阵
};

/**
 * @class SparseMatrixBase
 * @brief 稀疏矩阵抽象基类
 * @details 定义稀疏矩阵的统一接口，支持实数和复数矩阵
 */
class SparseMatrixBase {
public:
    /**
     * @brief 虚析构函数
     */
    virtual ~SparseMatrixBase() = default;

    /**
     * @brief 获取矩阵行数
     * @return 矩阵行数
     */
    virtual int rows() const = 0;

    /**
     * @brief 获取矩阵列数
     * @return 矩阵列数
     */
    virtual int cols() const = 0;

    /**
     * @brief 获取非零元素数量
     * @return 非零元素数量
     */
    virtual int nnz() const = 0;

    /**
     * @brief 清空矩阵
     */
    virtual void clear() = 0;

    /**
     * @brief 打印矩阵信息
     */
    virtual void print_info() const = 0;

    /**
     * @brief 获取矩阵数据类型
     * @return 矩阵数据类型
     */
    virtual MatrixDataType get_data_type() const = 0;
};

} // namespace numeric