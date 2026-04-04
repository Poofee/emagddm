#pragma once

#include "tool/maxwell_reader.hpp"
#include "tool/maxwell_parser.hpp"
#include <memory>
#include <string>

namespace tool {

/**
 * @brief Maxwell解析器实现类
 * @details 基于maxwell_parser::MaxwellParser实现IMaxwellParser接口，
 *          负责AEDT文件的全部解析工作，包括项目数据、预览图像等
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
    std::vector<nlohmann::json> extractMaterials() const;
    
    /**
     * @brief 从单个材料块提取完整材料数据
     * @param material_block 材料块节点
     * @return 材料JSON数据
     */
    nlohmann::json extractSingleMaterial(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) const;
    
    /**
     * @brief 从材料块提取磁导率数据（线性或非线性）
     * @param material_block 材料块节点
     * @param material 输出JSON对象
     */
    void extractPermeability(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;
    
    /**
     * @brief 从材料块提取矫顽力数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象
     */
    void extractMagneticCoercivity(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;
    
    /**
     * @brief 从属性中提取双精度数值
     * @param prop 属性节点
     * @return 双精度数值
     */
    double propertyValueToDouble(const fe_em::tool::maxwell_parser::Property& prop) const;
    
    /**
     * @brief 从属性中提取字符串值
     * @param prop 属性节点
     * @return 字符串值
     */
    std::string propertyValueToString(const fe_em::tool::maxwell_parser::Property& prop) const;
    
    /**
     * @brief 从Points属性提取BH曲线数据
     * @param prop Points属性节点
     * @return BH曲线数据 [[H1,B1], [H2,B2], ...]
     */
    std::vector<std::vector<double>> extractBHCurve(
        const fe_em::tool::maxwell_parser::Property& prop) const;
    
    /**
     * @brief 从解析树中提取边界条件数据
     */
    std::vector<nlohmann::json> extractBoundaries() const;
    
    /**
     * @brief 从解析树中提取激励源数据
     */
    std::vector<nlohmann::json> extractExcitations() const;
    
    /**
     * @brief 从解析树中提取求解设置数据
     */
    nlohmann::json extractSolutionSetup() const;
    
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

    /**
     * @brief 从材料块提取温度相关B-H曲线数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象（将写入 temperature_bh_curves 字段）
     */
    void extractTemperatureBHCurves(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;

    /**
     * @brief 从材料块提取CoreLoss多频率曲线数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象（将写入 core_loss_multicurve 字段）
     */
    void extractCoreLossMulticurve(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;

    /**
     * @brief 从材料块提取铁损系数设置数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象（将写入 coefficient_setup 字段）
     */
    void extractCoeffSetup(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;

    /**
     * @brief 从材料块提取热修正器数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象（将写入 thermal_modifiers 字段）
     */
    void extractThermalModifiers(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;

    /**
     * @brief 从材料块提取外观数据
     * @param material_block 材料块节点
     * @param material 输出JSON对象（将写入 appearance 字段）
     */
    void extractMatAppearance(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
        nlohmann::json& material) const;

    /**
     * @brief 从属性安全提取双精度数值（带默认值）
     * @param block 块节点
     * @param key 属性名
     * @param default_val 默认值
     * @return 提取到的双精度值或默认值
     */
    double safeGetDouble(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& key,
        double default_val = 0.0) const;

    /**
     * @brief 从属性安全提取整型数值（带默认值）
     * @param block 块节点
     * @param key 属性名
     * @param default_val 默认值
     * @return 提取到的整型值或默认值
     */
    int safeGetInt(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& key,
        int default_val = 0) const;

    /**
     * @brief 从属性安全提取字符串值（带默认值）
     * @param block 块节点
     * @param key 属性名
     * @param default_val 默认值
     * @return 提取到的字符串值或默认值
     */
    std::string safeGetString(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& key,
        const std::string& default_val = "") const;

    /**
     * @brief 从属性安全提取布尔值（带默认值）
     * @param block 块节点
     * @param key 属性名
     * @param default_val 默认值
     * @return 提取到的布尔值或默认值
     */
    bool safeGetBool(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& key,
        bool default_val = false) const;

    /**
     * @brief 从属性安全提取int64_t数值（带默认值）
     * @param block 块节点
     * @param key 属性名
     * @param default_val 默认值
     * @return 提取到的int64_t值或默认值
     */
    int64_t safeGetInt64(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& key,
        int64_t default_val = 0) const;

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
     * @brief 提取AEDT文件中的预览图像（Base64编码）
     *
     * AEDT文件中的图像数据存储在 Image64='...' 属性中，可能跨越多行。
     * 本函数负责定位并提取完整的Base64字符串。
     *
     * @return Base64编码的图像数据，失败返回空字符串
     */
    std::string extractPreviewImage() const;

    /**
     * @brief 提取绕组组边界数据（BoundType='Winding Group'）
     * @return 绕组组JSON列表
     */
    std::vector<nlohmann::json> extractWindings() const;

    /**
     * @brief 提取运动设置数据（MotionSetupList区域）
     * @return 运动设置JSON列表
     */
    std::vector<nlohmann::json> extractMotionSetups() const;

    /**
     * @brief 提取全局网格设置（MeshSettings子块）
     * @return 全局网格设置JSON对象
     */
    nlohmann::json extractGlobalMeshSettings() const;

    /**
     * @brief 提取网格操作列表（MeshOperations子块）
     * @return 网格操作JSON列表
     */
    std::vector<nlohmann::json> extractMeshOperations() const;

    /**
     * @brief 提取设计变量（ModelSetup > Properties > VariableProp）
     * @return 设计变量JSON列表
     */
    std::vector<nlohmann::json> extractDesignVariables() const;

    /**
     * @brief 提取非索引变量（NonIndexedVariables区域）
     * @return 非索引变量JSON列表
     */
    std::vector<nlohmann::json> extractNonIndexedVariables() const;

    /**
     * @brief 提取输出变量（OutputVariable > OutputVariables区域）
     * @return 输出变量JSON列表
     */
    std::vector<nlohmann::json> extractOutputVariables() const;

    /**
     * @brief 提取温度设置（TemperatureSettings区域）
     * @return 温度设置JSON对象
     */
    nlohmann::json extractTemperatureSettings() const;

    /**
     * @brief 提取全局边界数据（BoundarySetup顶层属性）
     * @return 全局边界数据JSON对象
     */
    nlohmann::json extractGlobalBoundData() const;
    
    /**
     * @brief 获取解析器实例（用于调试）
     */
    fe_em::tool::maxwell_parser::MaxwellParser& getParser() { return parser_; }
};

} // namespace tool
