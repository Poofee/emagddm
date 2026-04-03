/**
 * @file base64_utils.hpp
 * @brief Base64编码/解码工具
 * @author AI Developer
 * @date 2026-04-03
 * @version 1.0
 */

#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace fe_em {
namespace tool {

/**
 * @brief Base64编解码工具类
 */
class Base64Utils {
public:
    /**
     * @brief Base64解码
     * @param encoded Base64编码的字符串
     * @return 解码后的二进制数据
     */
    static std::vector<uint8_t> decode(const std::string& encoded) {
        static const std::string base64_chars = 
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        
        std::vector<uint8_t> decoded;
        
        if (encoded.empty()) {
            return decoded;
        }
        
        std::string clean_encoded;
        for (char c : encoded) {
            if (c != '\n' && c != '\r' && c != ' ' && c != '\\') {
                clean_encoded += c;
            }
        }
        
        size_t padding = 0;
        if (clean_encoded.length() >= 1 && clean_encoded[clean_encoded.length() - 1] == '=') padding++;
        if (clean_encoded.length() >= 2 && clean_encoded[clean_encoded.length() - 2] == '=') padding++;
        
        if (clean_encoded.length() % 4 != 0) {
            throw std::runtime_error("Invalid Base64 string length");
        }
        
        decoded.reserve((clean_encoded.length() / 4) * 3);
        
        uint32_t val = 0;
        int valb = -8;
        
        for (size_t i = 0; i < clean_encoded.length(); i++) {
            if (clean_encoded[i] == '=') {
                valb -= 6;
                continue;
            }
            
            size_t pos = base64_chars.find(clean_encoded[i]);
            if (pos == std::string::npos) {
                continue;
            }
            
            val = (val << 6) | static_cast<uint32_t>(pos);
            valb += 6;
            
            if (valb >= 0) {
                decoded.push_back(static_cast<uint8_t>((val >> valb) & 0xFF));
                valb -= 8;
            }
        }
        
        if (padding > 0 && decoded.size() >= padding) {
            decoded.resize(decoded.size() - padding);
        }
        
        return decoded;
    }
    
    /**
     * @brief Base64编码
     * @param data 二进制数据
     * @return Base64编码的字符串
     */
    static std::string encode(const std::vector<uint8_t>& data) {
        static const std::string base64_chars = 
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        
        std::string encoded;
        encoded.reserve(((data.size() + 2) / 3) * 4);
        
        uint32_t val = 0;
        int valb = -6;
        
        for (uint8_t c : data) {
            val = (val << 8) | c;
            valb += 8;
            while (valb >= 0) {
                encoded.push_back(base64_chars[(val >> valb) & 0x3F]);
                valb -= 6;
            }
        }
        
        if (valb > -6) {
            encoded.push_back(base64_chars[((val << 8) >> (valb + 8)) & 0x3F]);
        }
        
        while (encoded.size() % 4) {
            encoded.push_back('=');
        }
        
        return encoded;
    }
    
    /**
     * @brief 检测图像类型
     * @param data 图像二进制数据
     * @return 图像文件扩展名（如".jpg", ".png"等）
     */
    static std::string detectImageType(const std::vector<uint8_t>& data) {
        if (data.size() < 4) {
            return ".bin";
        }
        
        if (data[0] == 0xFF && data[1] == 0xD8 && data[2] == 0xFF) {
            return ".jpg";
        }
        
        if (data[0] == 0x89 && data[1] == 0x50 && data[2] == 0x4E && data[3] == 0x47) {
            return ".png";
        }
        
        if (data[0] == 0x47 && data[1] == 0x49 && data[2] == 0x46) {
            return ".gif";
        }
        
        if (data[0] == 0x42 && data[1] == 0x4D) {
            return ".bmp";
        }
        
        return ".bin";
    }
};

} // namespace tool
} // namespace fe_em
