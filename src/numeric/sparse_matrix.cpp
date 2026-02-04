/**
 * @file sparse_matrix.cpp
 * @brief 核心数值层 - 稀疏矩阵管理模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "sparse_matrix.hpp"

namespace numeric {

SparseMatrix::SparseMatrix() : rows_(0), cols_(0) {
}

SparseMatrix::~SparseMatrix() {
}

void SparseMatrix::setSize(int rows, int cols) {
    rows_ = rows;
    cols_ = cols;
    // TODO: 实现稀疏矩阵存储结构初始化
}

void SparseMatrix::insert(int row, int col, double value) {
    // TODO: 实现稀疏矩阵元素插入逻辑
}

int SparseMatrix::rows() const {
    return rows_;
}

int SparseMatrix::cols() const {
    return cols_;
}

} // namespace numeric