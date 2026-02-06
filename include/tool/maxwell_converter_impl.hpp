/**
 * @file maxwell_converter_impl.hpp
 * @brief Maxwell数据转换器实现
 * @details 提供从Maxwell解析数据到内部数据模型的直接转换功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include "tool/maxwell_parser.hpp"
#include "tool/maxwell_parser_impl.hpp"
#include <memory>
#include <string>
#include <unordered_map>

namespace tool {

/**
 * @brief Maxwell数据转换器实现类
 * @details 提供从Maxwell解析数据到内部数据模型的直接转换功能
 */
class MaxwellConverterImpl {
public:
    MaxwellConverterImpl() = default;
    ~MaxwellConverterImpl() = default;

    /**
     * @brief 将Maxwell材料数据转换为内部Material对象
     * @param material_block Maxwell材料块节点
     * @return 转换后的Material对象指针
     */
    std::shared_ptr<Material> convertMaterialDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 将Maxwell边界条件数据转换为内部Boundary对象
     * @param boundary_block Maxwell边界条件块节点
     * @return 转换后的Boundary对象指针
     */
    std::shared_ptr<Boundary> convertBoundaryDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block);

    /**
     * @brief 将Maxwell激励源数据转换为内部Excitation对象
     * @param excitation_block Maxwell激励源块节点
     * @return 转换后的Excitation对象指针
     */
    std::shared_ptr<Excitation> convertExcitationDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 将Maxwell几何数据转换为内部Geometry对象
     * @param geometry_block Maxwell几何块节点
     * @return 转换后的Geometry对象指针
     */
    std::shared_ptr<Geometry> convertGeometryDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& geometry_block);

    /**
     * @brief 将Maxwell求解设置数据转换为内部SolutionSetup对象
     * @param solution_block Maxwell求解设置块节点
     * @return 转换后的SolutionSetup对象指针
     */
    std::shared_ptr<SolutionSetup> convertSolutionSetupDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& solution_block);

    /**
     * @brief 将Maxwell项目文件信息转换为内部ProjectData对象
     * @param project_block Maxwell项目块节点
     * @return 转换后的ProjectData对象指针
     */
    std::shared_ptr<ProjectData> convertProjectDataDirect(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& project_block);



    /**
     * @brief 转换坐标系类型枚举
     * @param maxwell_coord_str Maxwell坐标系类型字符串
     * @return 对应的CoordinateSystemType枚举
     */
    CoordinateSystemType convertCoordinateSystemType(const std::string& maxwell_coord_str);

public:
    /**
     * @brief 转换材料基本属性
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialBasicProperties(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换材料电磁属性
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialElectromagneticProperties(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换材料B-H曲线数据
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialBHCurve(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换材料磁芯损耗参数
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialCoreLossParameters(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换材料各向异性参数
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialAnisotropicProperties(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换材料温度相关参数
     * @param material 目标Material对象
     * @param material_block Maxwell材料块节点
     */
    void convertMaterialTemperatureProperties(
        Material& material,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block);

    /**
     * @brief 转换边界条件基本属性
     * @param boundary 目标Boundary对象
     * @param boundary_block Maxwell边界条件块节点
     */
    void convertBoundaryBasicProperties(
        Boundary& boundary,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block);

    /**
     * @brief 转换边界条件参数
     * @param boundary 目标Boundary对象
     * @param boundary_block Maxwell边界条件块节点
     */
    void convertBoundaryParameters(
        Boundary& boundary,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block);

    /**
     * @brief 转换边界条件关联几何对象
     * @param boundary 目标Boundary对象
     * @param boundary_block Maxwell边界条件块节点
     */
    void convertBoundaryGeometryLinks(
        Boundary& boundary,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block);

    /**
     * @brief 转换Maxwell专属边界数据
     * @param boundary 目标Boundary对象
     * @param boundary_block Maxwell边界条件块节点
     */
    void convertBoundaryMaxwellSpecificProperties(
        Boundary& boundary,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block);

    /**
     * @brief 转换激励源基本属性
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationBasicProperties(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 转换激励源参数
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationParameters(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 转换激励源线圈参数
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationCoilParameters(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 转换激励源波形参数
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationWaveformParameters(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 转换激励源运动参数
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationMotionParameters(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 转换Maxwell专属激励数据
     * @param excitation 目标Excitation对象
     * @param excitation_block Maxwell激励源块节点
     */
    void convertExcitationMaxwellSpecificProperties(
        Excitation& excitation,
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block);

    /**
     * @brief 从属性值字符串转换为数值
     * @param value_str 属性值字符串
     * @return 转换后的数值
     */
    double parseNumericValue(const std::string& value_str);

    /**
     * @brief 从属性值字符串转换为布尔值
     * @param value_str 属性值字符串
     * @return 转换后的布尔值
     */
    bool parseBooleanValue(const std::string& value_str);

    /**
     * @brief 从属性值字符串转换为字符串数组
     * @param value_str 属性值字符串
     * @return 转换后的字符串数组
     */
    std::vector<std::string> parseStringArray(const std::string& value_str);

    /**
     * @brief 从属性值字符串转换为数值数组
     * @param value_str 属性值字符串
     * @return 转换后的数值数组
     */
    std::vector<double> parseNumericArray(const std::string& value_str);

    /**
     * @brief 查找属性值
     * @param block 块节点
     * @param property_name 属性名称
     * @return 属性值字符串，如果不存在则返回空字符串
     */
    std::string findPropertyValue(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& property_name);

    /**
     * @brief 判断属性是否存在
     * @param block 块节点
     * @param property_name 属性名称
     * @return 是否存在
     */
    bool hasProperty(
        const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
        const std::string& property_name);

    /**
     * @brief 转换材料类型枚举
     * @param maxwell_type_str Maxwell材料类型字符串
     * @return 对应的MatType枚举
     */
    MatType convertMaterialType(const std::string& maxwell_type_str);

    /**
     * @brief 转换B-H曲线类型枚举
     * @param maxwell_curve_type_str Maxwell曲线类型字符串
     * @return 对应的BHCurveType枚举
     */
    BHCurveType convertBHCurveType(const std::string& maxwell_curve_type_str);

    /**
     * @brief 转换磁芯损耗模型类型枚举
     * @param maxwell_model_str Maxwell模型类型字符串
     * @return 对应的CoreLossModelType枚举
     */
    CoreLossModelType convertCoreLossModelType(const std::string& maxwell_model_str);

    /**
     * @brief 转换边界条件类型枚举
     * @param maxwell_type_str Maxwell边界条件类型字符串
     * @return 对应的BndType枚举
     */
    BndType convertBoundaryType(const std::string& maxwell_type_str);

    /**
     * @brief 转换边界条件子类型枚举
     * @param maxwell_subtype_str Maxwell边界条件子类型字符串
     * @return 对应的BoundarySubType枚举
     */
    BoundarySubType convertBoundarySubType(const std::string& maxwell_subtype_str);

    /**
     * @brief 转换周期性映射类型枚举
     * @param maxwell_mapping_str Maxwell映射类型字符串
     * @return 对应的PeriodicMappingType枚举
     */
    PeriodicMappingType convertPeriodicMappingType(const std::string& maxwell_mapping_str);

    /**
     * @brief 转换激励源类型枚举
     * @param maxwell_type_str Maxwell激励源类型字符串
     * @return 对应的ExcitationType枚举
     */
    ExcitationType convertExcitationType(const std::string& maxwell_type_str);

    /**
     * @brief 转换激励源波形类型枚举
     * @param maxwell_waveform_str Maxwell波形类型字符串
     * @return 对应的ExcitationWaveformType枚举
     */
    ExcitationWaveformType convertExcitationWaveformType(const std::string& maxwell_waveform_str);

    /**
     * @brief 转换线圈连接类型枚举
     * @param maxwell_conn_str Maxwell连接类型字符串
     * @return 对应的CoilConnectionType枚举
     */
    CoilConnectionType convertCoilConnectionType(const std::string& maxwell_conn_str);

    /**
     * @brief 转换绕组类型枚举
     * @param maxwell_winding_str Maxwell绕组类型字符串
     * @return 对应的WindingType枚举
     */
    WindingType convertWindingType(const std::string& maxwell_winding_str);

    /**
     * @brief 转换运动类型枚举
     * @param maxwell_motion_str Maxwell运动类型字符串
     * @return 对应的MotionType枚举
     */
    MotionType convertMotionType(const std::string& maxwell_motion_str);
};

} // namespace tool