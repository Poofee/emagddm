/**
 * @file matrix_market_io_impl.hpp
 * @brief MatrixMarket格式稀疏矩阵I/O实现
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 */

#pragma once

#include "matrix_market_io.hpp"
#include <algorithm>
#include <cctype>

namespace numeric {

MatrixDataType MatrixMarketIO::read_coo(const std::string& filename, 
                                       CooMatrixReal& coo_real, 
                                       CooMatrixComplex& coo_complex) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件: " + filename);
    }

    std::string header_line;
    std::getline(file, header_line);

    MatrixMarketDataType data_type;
    MatrixMarketStructure structure;
    MatrixMarketSymmetry symmetry;
    parse_header(header_line, data_type, structure, symmetry);

    if (structure != MatrixMarketStructure::COORDINATE) {
        throw std::runtime_error("仅支持坐标格式稀疏矩阵");
    }

    skip_comments(file);

    int rows, cols, nnz;
    read_size_line(file, rows, cols, nnz);

    if (data_type == MatrixMarketDataType::REAL) {
        coo_real.set_size(rows, cols);
        coo_real.reserve(nnz);

        for (int i = 0; i < nnz; ++i) {
            int row, col;
            double value;
            read_data_line(file, row, col, value);
            coo_real.add_value(row - 1, col - 1, value); // MatrixMarket使用1-based索引
        }

        return MatrixDataType::REAL;
    } else if (data_type == MatrixMarketDataType::COMPLEX) {
        coo_complex.set_size(rows, cols);
        coo_complex.reserve(nnz);

        for (int i = 0; i < nnz; ++i) {
            int row, col;
            double real, imag;
            read_data_line(file, row, col, real, imag);
            coo_complex.add_value(row - 1, col - 1, std::complex<double>(real, imag));
        }

        return MatrixDataType::COMPLEX;
    } else {
        throw std::runtime_error("不支持的数据类型");
    }
}

MatrixDataType MatrixMarketIO::read_csr(const std::string& filename, 
                                       CsrMatrixReal& csr_real, 
                                       CsrMatrixComplex& csr_complex) {
    CooMatrixReal coo_real;
    CooMatrixComplex coo_complex;
    
    MatrixDataType data_type = read_coo(filename, coo_real, coo_complex);
    
    if (data_type == MatrixDataType::REAL) {
        csr_real.build_from_coo(coo_real);
    } else {
        csr_complex.build_from_coo(coo_complex);
    }
    
    return data_type;
}

void MatrixMarketIO::write_coo(const std::string& filename, const CooMatrixReal& coo) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法创建文件: " + filename);
    }

    int rows = coo.rows();
    int cols = coo.cols();
    int nnz = coo.nnz();

    write_header(file, rows, cols, nnz, false);

    const auto& row_indices = coo.get_row_indices();
    const auto& col_indices = coo.get_col_indices();
    const auto& values = coo.get_values();

    for (int i = 0; i < nnz; ++i) {
        write_data_line(file, row_indices[i] + 1, col_indices[i] + 1, values[i]);
    }
}

void MatrixMarketIO::write_coo(const std::string& filename, const CooMatrixComplex& coo) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法创建文件: " + filename);
    }

    int rows = coo.rows();
    int cols = coo.cols();
    int nnz = coo.nnz();

    write_header(file, rows, cols, nnz, true);

    const auto& row_indices = coo.get_row_indices();
    const auto& col_indices = coo.get_col_indices();
    const auto& values = coo.get_values();

    for (int i = 0; i < nnz; ++i) {
        write_data_line(file, row_indices[i] + 1, col_indices[i] + 1, 
                       values[i].real(), values[i].imag());
    }
}

void MatrixMarketIO::write_csr(const std::string& filename, const CsrMatrixReal& csr) {
    if (!csr.is_built()) {
        throw std::runtime_error("CSR矩阵未构建，无法写入");
    }

    // 转换为COO格式写入
    CooMatrixReal coo(csr.rows(), csr.cols());
    const auto& row_ptr = csr.get_row_ptr();
    const auto& col_indices = csr.get_col_indices();
    const auto& values = csr.get_values();

    for (int i = 0; i < csr.rows(); ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            coo.add_value(i, col_indices[j], values[j]);
        }
    }

    write_coo(filename, coo);
}

void MatrixMarketIO::write_csr(const std::string& filename, const CsrMatrixComplex& csr) {
    if (!csr.is_built()) {
        throw std::runtime_error("CSR矩阵未构建，无法写入");
    }

    // 转换为COO格式写入
    CooMatrixComplex coo(csr.rows(), csr.cols());
    const auto& row_ptr = csr.get_row_ptr();
    const auto& col_indices = csr.get_col_indices();
    const auto& values = csr.get_values();

    for (int i = 0; i < csr.rows(); ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            coo.add_value(i, col_indices[j], values[j]);
        }
    }

    write_coo(filename, coo);
}

void MatrixMarketIO::parse_header(const std::string& line, 
                                MatrixMarketDataType& data_type,
                                MatrixMarketStructure& structure,
                                MatrixMarketSymmetry& symmetry) {
    std::istringstream iss(line);
    std::string token;
    
    // 读取标记
    iss >> token;
    if (token != "%%MatrixMarket") {
        throw std::runtime_error("无效的MatrixMarket文件头");
    }

    // 读取矩阵类型
    iss >> token;
    if (token == "matrix") {
        // 继续读取
    } else {
        throw std::runtime_error("仅支持矩阵类型");
    }

    // 读取结构
    iss >> token;
    if (token == "coordinate") {
        structure = MatrixMarketStructure::COORDINATE;
    } else if (token == "array") {
        throw std::runtime_error("不支持数组格式");
    } else {
        throw std::runtime_error("无效的矩阵结构");
    }

    // 读取数据类型
    iss >> token;
    if (token == "real") {
        data_type = MatrixMarketDataType::REAL;
    } else if (token == "complex") {
        data_type = MatrixMarketDataType::COMPLEX;
    } else if (token == "integer") {
        throw std::runtime_error("不支持整数类型");
    } else if (token == "pattern") {
        throw std::runtime_error("不支持模式类型");
    } else {
        throw std::runtime_error("无效的数据类型");
    }

    // 读取对称性
    iss >> token;
    if (token == "general") {
        symmetry = MatrixMarketSymmetry::GENERAL;
    } else if (token == "symmetric") {
        symmetry = MatrixMarketSymmetry::SYMMETRIC;
    } else if (token == "skew-symmetric") {
        symmetry = MatrixMarketSymmetry::SKEW_SYMMETRIC;
    } else if (token == "hermitian") {
        symmetry = MatrixMarketSymmetry::HERMITIAN;
    } else {
        throw std::runtime_error("无效的对称性");
    }
}

void MatrixMarketIO::skip_comments(std::ifstream& file) {
    std::string line;
    
    // 获取当前位置
    std::streampos current_pos = file.tellg();
    
    while (std::getline(file, line)) {
        // 去除行首空格
        std::string trimmed_line = line;
        trimmed_line.erase(0, trimmed_line.find_first_not_of(" \t"));
        
        // 如果是注释行，继续读取下一行
        if (!trimmed_line.empty() && trimmed_line[0] == '%') {
            current_pos = file.tellg(); // 更新位置
            continue;
        }
        
        // 非注释行，回退到这一行的开始位置
        file.seekg(current_pos);
        break;
    }
}

void MatrixMarketIO::read_size_line(std::ifstream& file, int& rows, int& cols, int& nnz) {
    std::string line;
    
    // 读取尺寸信息行（跳过注释后应该就是这一行）
    if (!std::getline(file, line)) {
        throw std::runtime_error("文件结束，未找到尺寸信息行");
    }
    
    // 去除Windows换行符\r
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    
    std::istringstream iss(line);
    if (!(iss >> rows >> cols >> nnz)) {
        throw std::runtime_error("无效的尺寸信息行: " + line);
    }
    
    if (rows <= 0 || cols <= 0 || nnz < 0) {
        throw std::runtime_error("无效的矩阵尺寸");
    }
}

void MatrixMarketIO::read_data_line(std::ifstream& file, int& row, int& col, double& value) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    if (!(iss >> row >> col >> value)) {
        throw std::runtime_error("无效的数据行");
    }
}

void MatrixMarketIO::read_data_line(std::ifstream& file, int& row, int& col, 
                                   double& real, double& imag) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    if (!(iss >> row >> col >> real >> imag)) {
        throw std::runtime_error("无效的复数数据行");
    }
}

void MatrixMarketIO::write_header(std::ofstream& file, int rows, int cols, int nnz, bool is_complex) {
    file << "%%MatrixMarket matrix coordinate " << (is_complex ? "complex" : "real") << " general\n";
    file << "% Generated by Elmer electromagnetic FEM solver\n";
    file << rows << " " << cols << " " << nnz << "\n";
    file.flush(); // 确保数据写入文件
}

void MatrixMarketIO::write_data_line(std::ofstream& file, int row, int col, double value) {
    file << row << " " << col << " " << value << "\n";
}

void MatrixMarketIO::write_data_line(std::ofstream& file, int row, int col, double real, double imag) {
    file << row << " " << col << " " << real << " " << imag << "\n";
}

} // namespace numeric