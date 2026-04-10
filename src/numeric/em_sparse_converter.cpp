/**
 * @file em_sparse_converter.cpp
 * @brief 核心数值层 - CSR矩阵与Eigen稀疏矩阵转换工具类实现
 * @details 实现SparseConverter类的所有方法，包括双向转换、维度校验、性能统计等功能
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "em_sparse_converter.h"
#include "coo_matrix.hpp"
#include "logger_factory.hpp"

namespace numeric {

// 静态成员变量初始化
SparseConverter::ConversionStats SparseConverter::last_stats_ = {0.0, 0, 0, false};

Eigen::SparseMatrix<double> SparseConverter::to_eigen(const CsrMatrix<double>& csr) {
    // 记录转换开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1: 校验输入矩阵合法性
    validate_matrix_dimensions(csr);

    // 步骤2: 提取CSR内部数据
    const std::vector<int>& row_ptr = csr.get_row_ptr();
    const std::vector<int>& col_indices = csr.get_col_indices();
    const std::vector<double>& values = csr.get_values();

    int num_rows = csr.rows();
    int num_cols = csr.cols();
    int nnz = csr.nnz();

    // 记录源矩阵非零元素数
    last_stats_.source_nnz = nnz;

    // 大型矩阵转换日志提示
    if (nnz > 10000) {
        FEEM_DEBUG("开始大型CSR->Eigen转换，矩阵尺寸: {}x{}, 非零元数: {}", 
                   num_rows, num_cols, nnz);
    }

    // 步骤3: 使用Eigen::Triplets构建稀疏矩阵
    // 预分配三元组空间以避免多次内存重分配（性能优化）
    std::vector<Eigen::Triplet<double>> triplets;
    if (nnz > 10000) {
        triplets.reserve(nnz);
    }

    // 按行主序遍历CSR数据并构建三元组
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            triplets.emplace_back(i, col_indices[j], values[j]);
        }
    }

    // 构建Eigen稀疏矩阵并设置数据
    Eigen::SparseMatrix<double> eigen_mat(num_rows, num_cols);
    eigen_mat.setFromTriplets(triplets.begin(), triplets.end());

    // 确保矩阵压缩存储以优化后续操作
    eigen_mat.makeCompressed();

    // 步骤4: 统计转换耗时和结果信息
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    last_stats_.conversion_time_ms = duration.count() / 1000.0;
    last_stats_.target_nnz = eigen_mat.nonZeros();
    last_stats_.dimensions_match = true;

    FEEM_INFO("CSR->Eigen转换完成，耗时: {:.3f}ms, 非零元数: {}->{}", 
              last_stats_.conversion_time_ms, nnz, eigen_mat.nonZeros());

    // 使用移动语义返回，避免拷贝开销
    return eigen_mat;
}

CsrMatrix<double> SparseConverter::from_eigen(const Eigen::SparseMatrix<double>& eigen_mat) {
    // 记录转换开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1: 校验Eigen矩阵合法性
    if (eigen_mat.size() == 0) {
        throw std::invalid_argument("SparseConverter::from_eigen: Eigen矩阵为空");
    }

    int num_rows = eigen_mat.rows();
    int num_cols = eigen_mat.cols();

    if (num_rows <= 0 || num_cols <= 0) {
        throw std::invalid_argument("SparseConverter::from_eigen: 矩阵维度不合法，rows=" + 
                                    std::to_string(num_rows) + ", cols=" + std::to_string(num_cols));
    }

    // 记录源矩阵非零元素数
    last_stats_.source_nnz = eigen_mat.nonZeros();

    // 大型矩阵转换日志提示
    if (eigen_mat.nonZeros() > 10000) {
        FEEM_DEBUG("开始大型Eigen->CSR转换，矩阵尺寸: {}x{}, 非零元数: {}", 
                   num_rows, num_cols, eigen_mat.nonZeros());
    }

    // 步骤2: 确保Eigen矩阵为压缩存储格式
    // 注意：const引用无法直接调用makeCompressed()，需通过副本操作
    // 但由于输入是const的，我们假设调用者已确保矩阵压缩或使用内部指针安全访问

    // 步骤3: 将Eigen的CSC格式转换为CSR格式
    std::vector<int> row_ptr;
    std::vector<int> col_indices;
    std::vector<double> values;

    convert_csc_to_csr(eigen_mat, row_ptr, col_indices, values);

    // 步骤4: 通过COO中间格式构建CsrMatrix
    // 创建COO矩阵并预分配容量（性能优化）
    CooMatrix<double> coo(num_rows, num_cols, static_cast<int>(values.size()));

    // 将CSR格式数据展开为COO格式（按行主序遍历）
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            coo.add_value(i, col_indices[j], values[j]);
        }
    }

    // 构造CsrMatrix并通过build_from_coo填充数据
    CsrMatrix<double> csr(num_rows, num_cols);
    csr.build_from_coo(coo);

    // 步骤5: 统计转换耗时和结果信息
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    last_stats_.conversion_time_ms = duration.count() / 1000.0;
    last_stats_.target_nnz = csr.nnz();
    last_stats_.dimensions_match = true;

    FEEM_INFO("Eigen->CSR转换完成，耗时: {:.3f}ms, 非零元数: {}->{}", 
              last_stats_.conversion_time_ms, eigen_mat.nonZeros(), csr.nnz());

    // 使用移动语义返回，避免拷贝开销
    return csr;
}

void SparseConverter::validate_matrix_dimensions(const CsrMatrix<double>& csr,
                                                 int expected_rows,
                                                 int expected_cols) {
    // 检查矩阵是否已构建
    if (!csr.is_built()) {
        throw std::invalid_argument("SparseConverter: 输入矩阵未构建(is_built()==false)");
    }

    // 检查行数合法性
    int actual_rows = csr.rows();
    if (actual_rows <= 0) {
        throw std::invalid_argument("SparseConverter: 矩阵行数不合法，rows=" + 
                                    std::to_string(actual_rows));
    }

    // 检查列数合法性
    int actual_cols = csr.cols();
    if (actual_cols <= 0) {
        throw std::invalid_argument("SparseConverter: 矩阵列数不合法，cols=" + 
                                    std::to_string(actual_cols));
    }

    // 检查期望维度匹配（如果提供了期望值）
    if (expected_rows >= 0 && actual_rows != expected_rows) {
        throw std::invalid_argument("SparseConverter: 行数不匹配，期望=" + 
                                    std::to_string(expected_rows) + 
                                    ", 实际=" + std::to_string(actual_rows));
    }

    if (expected_cols >= 0 && actual_cols != expected_cols) {
        throw std::invalid_argument("SparseConverter: 列数不匹配，期望=" +
                                    std::to_string(expected_cols) +
                                    ", 实际=" + std::to_string(actual_cols));
    }
}

void SparseConverter::validate_matrix_dimensions(const CsrMatrix<std::complex<double>>& csr,
                                                 int expected_rows,
                                                 int expected_cols) {
    // 检查矩阵是否已构建
    if (!csr.is_built()) {
        throw std::invalid_argument("SparseConverter: 输入复数矩阵未构建(is_built()==false)");
    }

    // 检查行数合法性
    int actual_rows = csr.rows();
    if (actual_rows <= 0) {
        throw std::invalid_argument("SparseConverter: 复数矩阵行数不合法，rows=" +
                                    std::to_string(actual_rows));
    }

    // 检查列数合法性
    int actual_cols = csr.cols();
    if (actual_cols <= 0) {
        throw std::invalid_argument("SparseConverter: 复数矩阵列数不合法，cols=" +
                                    std::to_string(actual_cols));
    }

    // 检查期望维度匹配（如果提供了期望值）
    if (expected_rows >= 0 && actual_rows != expected_rows) {
        throw std::invalid_argument("SparseConverter: 复数矩阵行数不匹配，期望=" +
                                    std::to_string(expected_rows) +
                                    ", 实际=" + std::to_string(actual_rows));
    }

    if (expected_cols >= 0 && actual_cols != expected_cols) {
        throw std::invalid_argument("SparseConverter: 复数矩阵列数不匹配，期望=" +
                                    std::to_string(expected_cols) +
                                    ", 实际=" + std::to_string(actual_cols));
    }
}

SparseConverter::ConversionStats SparseConverter::get_last_conversion_stats() {
    return last_stats_;
}

void SparseConverter::convert_csc_to_csr(const Eigen::SparseMatrix<double>& eigen_mat,
                                         std::vector<int>& row_ptr,
                                         std::vector<int>& col_indices,
                                         std::vector<double>& values) {
    // 获取Eigen矩阵基本信息
    int num_rows = eigen_mat.rows();
    int num_cols = eigen_mat.cols();
    int nnz = eigen_mat.nonZeros();

    // 预分配输出数组空间（性能优化）
    row_ptr.resize(num_rows + 1, 0);
    col_indices.resize(nnz);
    values.resize(nnz);

    // 初始化行计数器数组
    std::vector<int> row_counts(num_rows, 0);

    // 第一遍扫描：统计每行的非零元素数量
    // Eigen的innerIndexPtr()返回行索引（对于CSC格式）
    const int* inner_idx = eigen_mat.innerIndexPtr();
    for (int i = 0; i < nnz; ++i) {
        row_counts[inner_idx[i]]++;
    }

    // 构建行偏移数组（前缀和）
    row_ptr[0] = 0;
    for (int i = 0; i < num_rows; ++i) {
        row_ptr[i + 1] = row_ptr[i] + row_counts[i];
    }

    // 重置行计数器用于第二遍扫描
    std::fill(row_counts.begin(), row_counts.end(), 0);

    // 第二遍扫描：按行主序填充列索引和数值
    const int* outer_ptr = eigen_mat.outerIndexPtr();  // CSC的列指针
    const double* eigen_values = eigen_mat.valuePtr();

    // 遍历每一列（CSC的外层循环）
    for (int col = 0; col < num_cols; ++col) {
        // 遍历当前列的所有非零元素
        for (int j = outer_ptr[col]; j < outer_ptr[col + 1]; ++j) {
            int row = inner_idx[j];

            // 计算该元素在CSR格式中的位置
            int pos = row_ptr[row] + row_counts[row];

            // 填充列索引和数值
            col_indices[pos] = col;
            values[pos] = eigen_values[j];

            // 更新行计数器
            row_counts[row]++;
        }
    }
}

Eigen::SparseMatrix<std::complex<double>> SparseConverter::to_eigen_complex(const CsrMatrix<std::complex<double>>& csr) {
    // 记录转换开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1: 校验输入矩阵合法性
    validate_matrix_dimensions(csr);

    // 步骤2: 提取CSR内部数据
    const std::vector<int>& row_ptr = csr.get_row_ptr();
    const std::vector<int>& col_indices = csr.get_col_indices();
    const std::vector<std::complex<double>>& values = csr.get_values();

    int num_rows = csr.rows();
    int num_cols = csr.cols();
    int nnz = csr.nnz();

    // 记录源矩阵非零元素数
    last_stats_.source_nnz = nnz;

    // 大型矩阵转换日志提示
    if (nnz > 10000) {
        FEEM_DEBUG("开始大型CSR->Eigen转换（复数），矩阵尺寸: {}x{}, 非零元数: {}", 
                   num_rows, num_cols, nnz);
    }

    // 步骤3: 使用Eigen::Triplets构建稀疏矩阵
    // 预分配三元组空间以避免多次内存重分配（性能优化）
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    if (nnz > 10000) {
        triplets.reserve(nnz);
    }

    // 按行主序遍历CSR数据并构建三元组
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            triplets.emplace_back(i, col_indices[j], values[j]);
        }
    }

    // 构建Eigen稀疏矩阵并设置数据
    Eigen::SparseMatrix<std::complex<double>> eigen_mat(num_rows, num_cols);
    eigen_mat.setFromTriplets(triplets.begin(), triplets.end());

    // 确保矩阵压缩存储以优化后续操作
    eigen_mat.makeCompressed();

    // 步骤4: 统计转换耗时和结果信息
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    last_stats_.conversion_time_ms = duration.count() / 1000.0;
    last_stats_.target_nnz = eigen_mat.nonZeros();
    last_stats_.dimensions_match = true;

    FEEM_INFO("CSR->Eigen转换完成（复数），耗时: {:.3f}ms, 非零元数: {}->{}", 
              last_stats_.conversion_time_ms, nnz, eigen_mat.nonZeros());

    // 使用移动语义返回，避免拷贝开销
    return eigen_mat;
}

CsrMatrix<std::complex<double>> SparseConverter::from_eigen_complex(const Eigen::SparseMatrix<std::complex<double>>& eigen_mat) {
    // 记录转换开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 步骤1: 校验Eigen矩阵合法性
    if (eigen_mat.size() == 0) {
        throw std::invalid_argument("SparseConverter::from_eigen_complex: Eigen矩阵为空");
    }

    int num_rows = eigen_mat.rows();
    int num_cols = eigen_mat.cols();

    if (num_rows <= 0 || num_cols <= 0) {
        throw std::invalid_argument("SparseConverter::from_eigen_complex: 矩阵维度不合法，rows=" + 
                                    std::to_string(num_rows) + ", cols=" + std::to_string(num_cols));
    }

    // 记录源矩阵非零元素数
    last_stats_.source_nnz = eigen_mat.nonZeros();

    // 大型矩阵转换日志提示
    if (eigen_mat.nonZeros() > 10000) {
        FEEM_DEBUG("开始大型Eigen->CSR转换（复数），矩阵尺寸: {}x{}, 非零元数: {}", 
                   num_rows, num_cols, eigen_mat.nonZeros());
    }

    // 步骤2: 确保Eigen矩阵为压缩存储格式
    // 注意：const引用无法直接调用makeCompressed()，需通过副本操作
    // 但由于输入是const的，我们假设调用者已确保矩阵压缩或使用内部指针安全访问

    // 步骤3: 将Eigen的CSC格式转换为CSR格式
    std::vector<int> row_ptr;
    std::vector<int> col_indices;
    std::vector<std::complex<double>> values;

    convert_csc_to_csr_complex(eigen_mat, row_ptr, col_indices, values);

    // 步骤4: 通过COO中间格式构建CsrMatrix
    // 创建COO矩阵并预分配容量（性能优化）
    CooMatrix<std::complex<double>> coo(num_rows, num_cols, static_cast<int>(values.size()));

    // 将CSR格式数据展开为COO格式（按行主序遍历）
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            coo.add_value(i, col_indices[j], values[j]);
        }
    }

    // 构造CsrMatrix并通过build_from_coo填充数据
    CsrMatrix<std::complex<double>> csr(num_rows, num_cols);
    csr.build_from_coo(coo);

    // 步骤5: 统计转换耗时和结果信息
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    last_stats_.conversion_time_ms = duration.count() / 1000.0;
    last_stats_.target_nnz = csr.nnz();
    last_stats_.dimensions_match = true;

    FEEM_INFO("Eigen->CSR转换完成（复数），耗时: {:.3f}ms, 非零元数: {}->{}", 
              last_stats_.conversion_time_ms, eigen_mat.nonZeros(), csr.nnz());

    // 使用移动语义返回，避免拷贝开销
    return csr;
}

void SparseConverter::convert_csc_to_csr_complex(const Eigen::SparseMatrix<std::complex<double>>& eigen_mat,
                                                 std::vector<int>& row_ptr,
                                                 std::vector<int>& col_indices,
                                                 std::vector<std::complex<double>>& values) {
    // 获取Eigen矩阵基本信息
    int num_rows = eigen_mat.rows();
    int num_cols = eigen_mat.cols();
    int nnz = eigen_mat.nonZeros();

    // 预分配输出数组空间（性能优化）
    row_ptr.resize(num_rows + 1, 0);
    col_indices.resize(nnz);
    values.resize(nnz);

    // 初始化行计数器数组
    std::vector<int> row_counts(num_rows, 0);

    // 第一遍扫描：统计每行的非零元素数量
    // Eigen的innerIndexPtr()返回行索引（对于CSC格式）
    const int* inner_idx = eigen_mat.innerIndexPtr();
    for (int i = 0; i < nnz; ++i) {
        row_counts[inner_idx[i]]++;
    }

    // 构建行偏移数组（前缀和）
    row_ptr[0] = 0;
    for (int i = 0; i < num_rows; ++i) {
        row_ptr[i + 1] = row_ptr[i] + row_counts[i];
    }

    // 重置行计数器用于第二遍扫描
    std::fill(row_counts.begin(), row_counts.end(), 0);

    // 第二遍扫描：按行主序填充列索引和数值
    const int* outer_ptr = eigen_mat.outerIndexPtr();  // CSC的列指针
    const std::complex<double>* eigen_values = eigen_mat.valuePtr();

    // 遍历每一列（CSC的外层循环）
    for (int col = 0; col < num_cols; ++col) {
        // 遍历当前列的所有非零元素
        for (int j = outer_ptr[col]; j < outer_ptr[col + 1]; ++j) {
            int row = inner_idx[j];

            // 计算该元素在CSR格式中的位置
            int pos = row_ptr[row] + row_counts[row];

            // 填充列索引和数值
            col_indices[pos] = col;
            values[pos] = eigen_values[j];

            // 更新行计数器
            row_counts[row]++;
        }
    }
}

} // namespace numeric
