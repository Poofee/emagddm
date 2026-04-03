/**
 * @file aedt_project_loader.hpp
 * @brief 基础工具层 - AEDT格式项目加载器
 * @details 实现IProjectLoader接口，负责将ANSYS Maxwell的AEDT文件解析为统一的项目数据结构
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include "project_loader.hpp"
#include "maxwell_parser_impl.hpp"
#include "base64_utils.hpp"

namespace tool {

/**
 * @class AedtProjectLoader
 * @brief AEDT格式项目加载器
 * 
 * 职责：读取 .aedt 文件 → 使用 MaxwellParserImpl 解析 → 转换为 ProjectManager 统一数据结构
 */
class AedtProjectLoader : public IProjectLoader {
public:
    AedtProjectLoader();
    ~AedtProjectLoader() override = default;

    void setFilePath(const std::string& file_path) override { file_path_ = file_path; }
    bool canLoad(const std::string& file_path) const override;
    bool load(ProjectManager& manager) override;
    InputFormat getSupportedFormat() const override { return InputFormat::AEDT; }
    std::string getFormatDescription() const override { return "ANSYS Maxwell AEDT"; }
    std::string getLastError() const override { return last_error_; }

private:
    /**
     * @brief 将JSON材料数据转换为Material对象并添加到项目管理器
     */
    void convertAndAddMaterials(ProjectManager& manager, 
                                 const std::vector<nlohmann::json>& materials_json);

    /**
     * @brief 将JSON边界条件数据转换为Boundary对象并添加到项目管理器
     */
    void convertAndAddBoundaries(ProjectManager& manager,
                                   const std::vector<nlohmann::json>& boundaries_json);

    /**
     * @brief 将JSON激励源数据转换为Excitation对象并添加到项目管理器
     */
    void convertAndAddExcitations(ProjectManager& manager,
                                   const std::vector<nlohmann::json>& excitations_json);

    /**
     * @brief 将JSON求解设置数据转换为SolutionSetup对象并添加到项目管理器
     */
    void convertAndAddSolutionSetup(ProjectManager& manager,
                                    const nlohmann::json& setup_json);

    /**
     * @brief 提取并保存预览图像（可选操作）
     */
    void extractPreviewImage(const std::string& file_path);

    std::string file_path_;
    std::string last_error_;
    MaxwellParserImpl parser_;
};

} // namespace tool
