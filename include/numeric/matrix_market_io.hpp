/**
 * @file matrix_market_io.hpp
 * @brief MatrixMarket格式稀疏矩阵I/O功能
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "sparse_base.hpp"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

namespace numeric {

/**
 * @brief MatrixMarket格式支持的数据类型
 */
enum class MatrixMarketDataType {
    REAL,           ///< 实数
    COMPLEX,        ///< 复数
    INTEGER,        ///< 整数（不支持）
    PATTERN         ///< 模式（不支持）
};

/**
 * @brief MatrixMarket格式支持的矩阵结构
 */
enum class MatrixMarketStructure {
    COORDINATE,     ///< 坐标格式（稀疏）
    ARRAY           ///< 数组格式（稠密，不支持）
};

/**
 * @brief MatrixMarket格式支持的对称性
 */
enum class MatrixMarketSymmetry {
    GENERAL,        ///< 一般矩阵
    SYMMETRIC,      ///< 对称矩阵
    SKEW_SYMMETRIC, ///< 反对称矩阵
    HERMITIAN       ///< 埃尔米特矩阵
};

/**
 * @class MatrixMarketIO
 * @brief MatrixMarket格式稀疏矩阵读写类
 */
class MatrixMarketIO {
public:
    /**
     * @brief 从MatrixMarket文件读取COO矩阵
     * @param filename 文件名
     * @param coo_real 实数COO矩阵（输出）
     * @param coo_complex 复数COO矩阵（输出）
     * @return 读取的矩阵数据类型
     */
    static MatrixDataType read_coo(const std::string& filename, 
                                  CooMatrixReal& coo_real, 
                                  CooMatrixComplex& coo_complex);

    /**
     * @brief 从MatrixMarket文件读取CSR矩阵
     * @param filename 文件名
     * @param csr_real 实数CSR矩阵（输出）
     * @param csr_complex 复数CSR矩阵（输出）
     * @return 读取的矩阵数据类型
     */
    static MatrixDataType read_csr(const std::string& filename, 
                                  CsrMatrixReal& csr_real, 
                                  CsrMatrixComplex& csr_complex);

    /**
     * @brief 将COO矩阵写入MatrixMarket文件
     * @param filename 文件名
     * @param coo COO矩阵
     */
    static void write_coo(const std::string& filename, const CooMatrixReal& coo);

    /**
     * @brief 将COO矩阵写入MatrixMarket文件
     * @param filename 文件名
     * @param coo 复数COO矩阵
     */
    static void write_coo(const std::string& filename, const CooMatrixComplex& coo);

    /**
     * @brief 将CSR矩阵写入MatrixMarket文件
     * @param filename 文件名
     * @param csr CSR矩阵
     */
    static void write_csr(const std::string& filename, const CsrMatrixReal& csr);

    /**
     * @brief 将CSR矩阵写入MatrixMarket文件
     * @param filename 文件名
     * @param csr 复数CSR矩阵
     */
    static void write_csr(const std::string& filename, const CsrMatrixComplex& csr);

private:
    /**
     * @brief 解析MatrixMarket文件头
     * @param line 文件头行
     * @param data_type 数据类型（输出）
     * @param structure 矩阵结构（输出）
     * @param symmetry 对称性（输出）
     */
    static void parse_header(const std::string& line, 
                            MatrixMarketDataType& data_type,
                            MatrixMarketStructure& structure,
                            MatrixMarketSymmetry& symmetry);

    /**
     * @brief 跳过注释行
     * @param file 文件流
     */
    static void skip_comments(std::ifstream& file);

    /**
     * @brief 读取尺寸信息行
     * @param file 文件流
     * @param rows 行数（输出）
     * @param cols 列数（输出）
     * @param nnz 非零元素数量（输出）
     */
    static void read_size_line(std::ifstream& file, int& rows, int& cols, int& nnz);

    /**
     * @brief 读取数据行（实数）
     * @param file 文件流
     * @param row 行索引（输出）
     * @param col 列索引（输出）
     * @param value 数值（输出）
     */
    static void read_data_line(std::ifstream& file, int& row, int& col, double& value);

    /**
     * @brief 读取数据行（复数）
     * @param file 文件流
     * @param row 行索引（输出）
     * @param col 列索引（输出）
     * @param real 实部（输出）
     * @param imag 虚部（输出）
     */
    static void read_data_line(std::ifstream& file, int& row, int& col, 
                              double& real, double& imag);

    /**
     * @brief 写入文件头
     * @param file 文件流
     * @param rows 行数
     * @param cols 列数
     * @param nnz 非零元素数量
     * @param is_complex 是否为复数
     */
    static void write_header(std::ofstream& file, int rows, int cols, int nnz, bool is_complex);

    /**
     * @brief 写入数据行（实数）
     * @param file 文件流
     * @param row 行索引
     * @param col 列索引
     * @param value 数值
     */
    static void write_data_line(std::ofstream& file, int row, int col, double value);

    /**
     * @brief 写入数据行（复数）
     * @param file 文件流
     * @param row 行索引
     * @param col 列索引
     * @param real 实部
     * @param imag 虚部
     */
    static void write_data_line(std::ofstream& file, int row, int col, double real, double imag);
};

} // namespace numeric

// 包含实现
#include "matrix_market_io_impl.hpp"