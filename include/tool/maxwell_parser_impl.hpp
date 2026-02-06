#pragma once

#include "tool/maxwell_reader.hpp"
#include "tool/maxwell_parser.hpp"
#include <memory>
#include <string>

namespace tool {

/**
 * @brief Maxwell解析器实现类
 * @details 基于maxwell_parser::MaxwellParser实现IMaxwellParser接口
 */
class MaxwellParserImpl : public IMaxwellParser {
private:
    fe_em::tool::maxwell_parser::MaxwellParser parser_;
    std::string file_path_;
    
    /**
     * @brief 从解析树中提取文件信息
     */
    MaxwellFileInfo extractFileInfo() const;
    
    /**
     * @brief 从解析树中提取材料数据
     */
    std::vector<std::shared_ptr<Material>> extractMaterials() const;
    
    /**
     * @brief 从解析树中提取边界条件数据
     */
    std::vector<std::shared_ptr<Boundary>> extractBoundaries() const;
    
    /**
     * @brief 从解析树中提取激励源数据
     */
    std::vector<std::shared_ptr<Excitation>> extractExcitations() const;
    
    /**
     * @brief 从解析树中提取求解设置数据
     */
    std::shared_ptr<SolutionSetup> extractSolutionSetup() const;
    
    /**
     * @brief 从解析树中提取几何数据
     */
    nlohmann::json extractGeometry() const;
    
    /**
     * @brief 检查解析树是否有效
     */
    bool isParseTreeValid() const;
    
    /**
     * @brief 查找指定名称的块
     */
    std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode> findBlock(
        const std::string& block_name,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& start_node) const;
    
    /**
     * @brief 转换属性值为JSON
     */
    nlohmann::json propertyToJson(const fe_em::tool::maxwell_parser::Property& prop) const;
    
    /**
     * @brief 转换块为JSON
     */
    nlohmann::json blockToJson(const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block) const;

public:
    /**
     * @brief 构造函数
     */
    MaxwellParserImpl() = default;
    
    /**
     * @brief 析构函数
     */
    virtual ~MaxwellParserImpl() = default;

    /**
     * @brief 检查是否可以解析指定文件
     * @param file_path 文件路径
     * @return 是否可以解析
     */
    bool canParse(const std::string& file_path) override;

    /**
     * @brief 解析文件信息
     * @return 文件信息
     */
    MaxwellFileInfo parseFileInfo() override;

    /**
     * @brief 解析材料数据
     * @return 材料数据列表
     */
    std::vector<nlohmann::json> parseMaterials() override;

    /**
     * @brief 解析边界条件数据
     * @return 边界条件数据列表
     */
    std::vector<nlohmann::json> parseBoundaries() override;

    /**
     * @brief 解析激励源数据
     * @return 激励源数据列表
     */
    std::vector<nlohmann::json> parseExcitations() override;

    /**
     * @brief 解析求解设置数据
     * @return 求解设置数据
     */
    nlohmann::json parseSolutionSetup() override;

    /**
     * @brief 解析几何数据
     * @return 几何数据
     */
    nlohmann::json parseGeometry() override;

    /**
     * @brief 获取所有解析的数据
     * @return 完整项目数据
     */
    nlohmann::json parseAllData() override;
    
    /**
     * @brief 获取解析器实例（用于调试）
     */
    fe_em::tool::maxwell_parser::MaxwellParser& getParser() { return parser_; }
};

} // namespace tool