/**
 * @file test_image_extraction.cpp
 * @brief 测试从AEDT文件中提取Base64图像
 * @author AI Developer
 * @date 2026-04-03
 */

#include "tool/maxwell_parser.hpp"
#include "tool/base64_utils.hpp"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace fe_em::tool::maxwell_parser;
using namespace fe_em::tool;

/**
 * @brief 从AEDT文件中提取Base64图像数据
 * @param file_path AEDT文件路径
 * @return Base64编码的图像数据
 */
std::string extractImage64FromFile(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << file_path << std::endl;
        return "";
    }
    
    std::string line;
    std::string base64_data;
    bool in_image_section = false;
    
    while (std::getline(file, line)) {
        if (line.find("Image64='") != std::string::npos) {
            in_image_section = true;
            size_t start_pos = line.find("Image64='") + 9;
            base64_data = line.substr(start_pos);
            if (base64_data.back() == '\'') {
                base64_data.pop_back();
                break;
            }
            continue;
        }
        
        if (in_image_section) {
            for (char c : line) {
                if (c == '\'') {
                    in_image_section = false;
                    break;
                }
                if (c != '\\' && c != '\n' && c != '\r') {
                    base64_data += c;
                }
            }
            if (!in_image_section) {
                break;
            }
        }
    }
    
    file.close();
    return base64_data;
}

/**
 * @brief 保存二进制数据到文件
 * @param data 二进制数据
 * @param file_path 输出文件路径
 * @return 是否成功
 */
bool saveBinaryFile(const std::vector<uint8_t>& data, const std::string& file_path) {
    std::ofstream file(file_path, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法创建文件: " << file_path << std::endl;
        return false;
    }
    
    file.write(reinterpret_cast<const char*>(data.data()), data.size());
    file.close();
    
    return true;
}

/**
 * @brief 测试Base64解码功能
 */
void test_base64_decode() {
    std::cout << "=== 测试Base64解码功能 ===" << std::endl;
    
    std::string test_encoded = "SGVsbG8gV29ybGQh";
    std::vector<uint8_t> decoded = Base64Utils::decode(test_encoded);
    
    std::string result(decoded.begin(), decoded.end());
    std::cout << "解码结果: " << result << std::endl;
    
    assert(result == "Hello World!" && "Base64解码测试失败");
    
    std::cout << "Base64解码测试通过" << std::endl;
}

/**
 * @brief 测试图像类型检测
 */
void test_image_type_detection() {
    std::cout << "\n=== 测试图像类型检测 ===" << std::endl;
    
    std::vector<uint8_t> jpeg_header = {0xFF, 0xD8, 0xFF, 0xE0};
    std::string jpeg_type = Base64Utils::detectImageType(jpeg_header);
    std::cout << "JPEG类型检测: " << jpeg_type << std::endl;
    assert(jpeg_type == ".jpg" && "JPEG类型检测失败");
    
    std::vector<uint8_t> png_header = {0x89, 0x50, 0x4E, 0x47};
    std::string png_type = Base64Utils::detectImageType(png_header);
    std::cout << "PNG类型检测: " << png_type << std::endl;
    assert(png_type == ".png" && "PNG类型检测失败");
    
    std::cout << "图像类型检测测试通过" << std::endl;
}

/**
 * @brief 测试从Project49.aedt提取图像
 */
void test_extract_image_from_project49() {
    std::cout << "\n=== 测试从Project49.aedt提取图像 ===" << std::endl;
    
    std::string aedt_file = "docs/project/Project49.aedt";
    
    std::cout << "正在读取文件: " << aedt_file << std::endl;
    std::string base64_data = extractImage64FromFile(aedt_file);
    
    if (base64_data.empty()) {
        std::cerr << "未找到Image64数据" << std::endl;
        return;
    }
    
    std::cout << "Base64数据长度: " << base64_data.length() << " 字符" << std::endl;
    std::cout << "Base64数据前50字符: " << base64_data.substr(0, 50) << "..." << std::endl;
    
    std::cout << "\n正在解码Base64数据..." << std::endl;
    std::vector<uint8_t> image_data = Base64Utils::decode(base64_data);
    
    std::cout << "解码后数据大小: " << image_data.size() << " 字节" << std::endl;
    
    std::string image_type = Base64Utils::detectImageType(image_data);
    std::cout << "检测到图像类型: " << image_type << std::endl;
    
    std::string output_file = "build/output_image" + image_type;
    std::cout << "\n正在保存图像到: " << output_file << std::endl;
    
    if (saveBinaryFile(image_data, output_file)) {
        std::cout << "✅ 图像保存成功!" << std::endl;
        std::cout << "图像文件路径: " << output_file << std::endl;
    } else {
        std::cerr << "❌ 图像保存失败" << std::endl;
    }
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始图像提取测试..." << std::endl;
    
    try {
        test_base64_decode();
        test_image_type_detection();
        test_extract_image_from_project49();
        
        std::cout << "\n✅ 所有测试通过!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}
