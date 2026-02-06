/**
 * @file maxwell_reader.hpp
 * @brief Maxwell数据读取模块 - 核心接口定义
 * @details 定义Maxwell文件解析、数据转换和验证的核心接口
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "project_data.hpp"
#include "em_enums.hpp"
#include "file_utils.hpp"
#include "../../lib/nlohmann/json.hpp"
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
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
 * @brief 验证结果结构
 */
struct ValidationResult {
    bool is_valid;                            // 是否有效
    std::vector<std::string> warnings;        // 警告信息
    std::vector<std::string> errors;          // 错误信息
    std::string summary;                      // 摘要信息

    /**
     * @brief 构造函数
     */
    ValidationResult() : is_valid(true) {}

    /**
     * @brief 添加警告
     * @param warning 警告信息
     */
    void addWarning(const std::string& warning) {
        warnings.push_back(warning);
    }

    /**
     * @brief 添加错误
     * @param error 错误信息
     */
    void addError(const std::string& error) {
        errors.push_back(error);
        is_valid = false;
    }

    /**
     * @brief 生成摘要
     */
    void generateSummary() {
        summary = "验证结果: " + std::string(is_valid ? "通过" : "失败") + 
                  "，警告: " + std::to_string(warnings.size()) + 
                  "，错误: " + std::to_string(errors.size());
    }

    /**
     * @brief 转换为JSON格式
     * @return JSON对象
     */
    nlohmann::json toJson() const {
        return {
            {"is_valid", is_valid},
            {"warnings", warnings},
            {"errors", errors},
            {"summary", summary}
        };
    }
};

/**
 * @brief Maxwell文件解析器接口
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
     * @brief 解析材料数据
     * @return 材料数据列表
     */
    virtual std::vector<nlohmann::json> parseMaterials() = 0;

    /**
     * @brief 解析边界条件数据
     * @return 边界条件数据列表
     */
    virtual std::vector<nlohmann::json> parseBoundaries() = 0;

    /**
     * @brief 解析激励源数据
     * @return 激励源数据列表
     */
    virtual std::vector<nlohmann::json> parseExcitations() = 0;

    /**
     * @brief 解析求解设置数据
     * @return 求解设置数据
     */
    virtual nlohmann::json parseSolutionSetup() = 0;

    /**
     * @brief 解析几何数据
     * @return 几何数据
     */
    virtual nlohmann::json parseGeometry() = 0;

    /**
     * @brief 获取所有解析的数据
     * @return 完整项目数据
     */
    virtual nlohmann::json parseAllData() = 0;
};

/**
 * @brief Maxwell数据转换器接口
 */
class IMaxwellConverter {
public:
    virtual ~IMaxwellConverter() = default;

    /**
     * @brief 转换材料数据
     * @param material_data 原始材料数据
     * @return 转换后的材料对象
     */
    virtual std::shared_ptr<Material> convertMaterial(const nlohmann::json& material_data) = 0;

    /**
     * @brief 转换边界条件数据
     * @param boundary_data 原始边界条件数据
     * @return 转换后的边界条件对象
     */
    virtual std::shared_ptr<Boundary> convertBoundary(const nlohmann::json& boundary_data) = 0;

    /**
     * @brief 转换激励源数据
     * @param excitation_data 原始激励源数据
     * @return 转换后的激励源对象
     */
    virtual std::shared_ptr<Excitation> convertExcitation(const nlohmann::json& excitation_data) = 0;

    /**
     * @brief 转换求解设置数据
     * @param setup_data 原始求解设置数据
     * @return 转换后的求解设置对象
     */
    virtual std::shared_ptr<SolutionSetup> convertSolutionSetup(const nlohmann::json& setup_data) = 0;

    /**
     * @brief 转换几何数据
     * @param geometry_data 原始几何数据
     * @return 转换后的几何对象（预留接口）
     */
    virtual std::shared_ptr<void> convertGeometry(const nlohmann::json& geometry_data) = 0;

    /**
     * @brief 批量转换材料数据
     * @param materials_data 原始材料数据列表
     * @return 转换后的材料对象列表
     */
    virtual std::vector<std::shared_ptr<Material>> convertMaterials(const std::vector<nlohmann::json>& materials_data) = 0;

    /**
     * @brief 批量转换边界条件数据
     * @param boundaries_data 原始边界条件数据列表
     * @return 转换后的边界条件对象列表
     */
    virtual std::vector<std::shared_ptr<Boundary>> convertBoundaries(const std::vector<nlohmann::json>& boundaries_data) = 0;

    /**
     * @brief 批量转换激励源数据
     * @param excitations_data 原始激励源数据列表
     * @return 转换后的激励源对象列表
     */
    virtual std::vector<std::shared_ptr<Excitation>> convertExcitations(const std::vector<nlohmann::json>& excitations_data) = 0;
};

/**
 * @brief Maxwell数据验证器接口
 */
class IMaxwellValidator {
public:
    virtual ~IMaxwellValidator() = default;

    /**
     * @brief 验证材料数据
     * @param materials 材料对象列表
     * @return 验证结果
     */
    virtual ValidationResult validateMaterials(const std::vector<std::shared_ptr<Material>>& materials) = 0;

    /**
     * @brief 验证边界条件数据
     * @param boundaries 边界条件对象列表
     * @return 验证结果
     */
    virtual ValidationResult validateBoundaries(const std::vector<std::shared_ptr<Boundary>>& boundaries) = 0;

    /**
     * @brief 验证激励源数据
     * @param excitations 激励源对象列表
     * @return 验证结果
     */
    virtual ValidationResult validateExcitations(const std::vector<std::shared_ptr<Excitation>>& excitations) = 0;

    /**
     * @brief 验证求解设置数据
     * @param setup 求解设置对象
     * @return 验证结果
     */
    virtual ValidationResult validateSolutionSetup(const std::shared_ptr<SolutionSetup>& setup) = 0;

    /**
     * @brief 验证几何数据
     * @param geometry 几何对象
     * @return 验证结果（预留接口）
     */
    virtual ValidationResult validateGeometry(const std::shared_ptr<void>& geometry) = 0;

    /**
     * @brief 验证完整项目数据
     * @param materials 材料对象列表
     * @param boundaries 边界条件对象列表
     * @param excitations 激励源对象列表
     * @param setup 求解设置对象
     * @return 验证结果
     */
    virtual ValidationResult validateProjectData(
        const std::vector<std::shared_ptr<Material>>& materials,
        const std::vector<std::shared_ptr<Boundary>>& boundaries,
        const std::vector<std::shared_ptr<Excitation>>& excitations,
        const std::shared_ptr<SolutionSetup>& setup) = 0;
};

/**
 * @brief Maxwell数据读取器主类
 */
class MaxwellReader {
private:
    std::unique_ptr<IMaxwellParser> parser_;
    std::unique_ptr<IMaxwellConverter> converter_;
    std::unique_ptr<IMaxwellValidator> validator_;
    std::string file_path_;

public:
    /**
     * @brief 构造函数
     * @param parser 文件解析器
     * @param converter 数据转换器
     * @param validator 数据验证器
     */
    MaxwellReader(std::unique_ptr<IMaxwellParser> parser,
                  std::unique_ptr<IMaxwellConverter> converter,
                  std::unique_ptr<IMaxwellValidator> validator);

    /**
     * @brief 析构函数
     */
    ~MaxwellReader() = default;

    /**
     * @brief 设置要读取的文件路径
     * @param file_path 文件路径
     */
    void setFilePath(const std::string& file_path);

    /**
     * @brief 读取并转换整个项目数据
     * @return 是否成功
     */
    bool readProject();

    /**
     * @brief 获取材料数据
     * @return 材料对象列表
     */
    std::vector<std::shared_ptr<Material>> getMaterials() const;

    /**
     * @brief 获取边界条件数据
     * @return 边界条件对象列表
     */
    std::vector<std::shared_ptr<Boundary>> getBoundaries() const;

    /**
     * @brief 获取激励源数据
     * @return 激励源对象列表
     */
    std::vector<std::shared_ptr<Excitation>> getExcitations() const;

    /**
     * @brief 获取求解设置数据
     * @return 求解设置对象
     */
    std::shared_ptr<SolutionSetup> getSolutionSetup() const;

    /**
     * @brief 获取验证结果
     * @return 验证结果
     */
    ValidationResult getValidationResult() const;

    /**
     * @brief 获取文件信息
     * @return 文件信息
     */
    MaxwellFileInfo getFileInfo() const;

    /**
     * @brief 导出为JSON格式
     * @return JSON对象
     */
    nlohmann::json exportToJson() const;

    /**
     * @brief 保存为项目文件
     * @param output_path 输出路径
     * @return 是否成功
     */
    bool saveToProject(const std::string& output_path) const;

private:
    std::vector<std::shared_ptr<Material>> materials_;
    std::vector<std::shared_ptr<Boundary>> boundaries_;
    std::vector<std::shared_ptr<Excitation>> excitations_;
    std::shared_ptr<SolutionSetup> solution_setup_;
    MaxwellFileInfo file_info_;
    ValidationResult validation_result_;
};

} // namespace tool