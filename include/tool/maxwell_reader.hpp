/**
 * @file maxwell_reader.hpp
 * @brief Maxwell数据读取模块 - 核心类型定义
 * @details 定义Maxwell文件解析所需的错误码、异常类、文件信息结构和解析器接口
 * @author AI Developer
 * @date 2026-02-06
 * @version 2.0
 */

#pragma once

#include "project_data.hpp"
#include "em_enums.hpp"
#include "file_utils.hpp"
#include "../../lib/nlohmann/json.hpp"
#include <memory>
#include <vector>
#include <string>
#include <exception>

namespace tool {

/**
 * @brief Maxwell错误代码枚举
 */
enum class MaxwellErrorCode {
    FILE_NOT_FOUND,           // 文件未找到
    INVALID_FORMAT,           // 无效的文件格式
    VERSION_NOT_SUPPORTED,    // 版本不支持
    DATA_CORRUPTED,           // 数据损坏
    CONVERSION_FAILED,        // 数据转换失败
    VALIDATION_FAILED         // 数据验证失败
};

/**
 * @brief Maxwell异常类
 */
class MaxwellException : public std::exception {
private:
    MaxwellErrorCode error_code_;
    std::string message_;

public:
    /**
     * @brief 构造函数
     * @param code 错误代码
     * @param message 错误消息
     */
    MaxwellException(MaxwellErrorCode code, const std::string& message)
        : error_code_(code), message_(message) {}

    /**
     * @brief 获取错误代码
     * @return 错误代码
     */
    MaxwellErrorCode getErrorCode() const { return error_code_; }

    /**
     * @brief 获取错误消息
     * @return 错误消息
     */
    const char* what() const noexcept override { return message_.c_str(); }
};

/**
 * @brief Maxwell文件信息结构
 */
struct MaxwellFileInfo {
    std::string file_path;                    // 文件路径
    std::string file_format;                  // 文件格式
    std::string maxwell_version;              // Maxwell版本
    std::string simulation_type;              // 仿真类型
    std::string dimension;                    // 维度
    std::string project_name;                 // 项目名称
    std::string created_date;                 // 创建日期
    std::string modified_date;                // 修改日期
    uint64_t file_size;                       // 文件大小
    bool is_compressed;                       // 是否压缩

    /**
     * @brief 转换为JSON格式
     * @return JSON对象
     */
    nlohmann::json toJson() const {
        return {
            {"file_path", file_path},
            {"file_format", file_format},
            {"maxwell_version", maxwell_version},
            {"simulation_type", simulation_type},
            {"dimension", dimension},
            {"project_name", project_name},
            {"created_date", created_date},
            {"modified_date", modified_date},
            {"file_size", file_size},
            {"is_compressed", is_compressed}
        };
    }

    /**
     * @brief 从JSON格式解析
     * @param json JSON对象
     * @return 是否成功
     */
    bool fromJson(const nlohmann::json& json) {
        try {
            file_path = json.value("file_path", "");
            file_format = json.value("file_format", "");
            maxwell_version = json.value("maxwell_version", "");
            simulation_type = json.value("simulation_type", "");
            dimension = json.value("dimension", "");
            project_name = json.value("project_name", "");
            created_date = json.value("created_date", "");
            modified_date = json.value("modified_date", "");
            file_size = json.value("file_size", 0ULL);
            is_compressed = json.value("is_compressed", false);
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }
};

/**
 * @brief Maxwell文件解析器接口
 * @details 定义AEDT文件解析的统一接口，由具体实现类（如MaxwellParserImpl）继承
 */
class IMaxwellParser {
public:
    virtual ~IMaxwellParser() = default;

    /**
     * @brief 检查是否可以解析指定文件
     * @param file_path 文件路径
     * @return 是否可以解析
     */
    virtual bool canParse(const std::string& file_path) = 0;

    /**
     * @brief 解析文件信息
     * @return 文件信息
     */
    virtual MaxwellFileInfo parseFileInfo() = 0;

    /**
     * @brief 解析材料数据（返回JSON中间格式）
     * @return 材料数据列表
     */
    virtual std::vector<nlohmann::json> parseMaterials() = 0;

    /**
     * @brief 解析边界条件数据（返回JSON中间格式）
     * @return 边界条件数据列表
     */
    virtual std::vector<nlohmann::json> parseBoundaries() = 0;

    /**
     * @brief 解析激励源数据（返回JSON中间格式）
     * @return 激励源数据列表
     */
    virtual std::vector<nlohmann::json> parseExcitations() = 0;

    /**
     * @brief 解析求解设置数据（返回JSON中间格式）
     * @return 求解设置数据
     */
    virtual nlohmann::json parseSolutionSetup() = 0;

    /**
     * @brief 解析几何数据（返回JSON中间格式）
     * @return 几何数据
     */
    virtual nlohmann::json parseGeometry() = 0;

    /**
     * @brief 获取所有解析的数据
     * @return 完整项目数据
     */
    virtual nlohmann::json parseAllData() = 0;
};

} // namespace tool
