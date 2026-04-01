#pragma once

#include "maxwell_data_structures.hpp"
#include "maxwell_parser.hpp"
#include "logger_factory.hpp"
#include <memory>
#include <unordered_map>
#include <functional>

namespace fe_em {
namespace tool {

/**
 * @brief Maxwell数据转换器命名空间
 * 负责将解析器解析的原始数据转换为结构化数据模型
 */
namespace maxwell_converter {

// 日志宏
#define LOG_CONVERTER(level, ...) \\
    LOGGER_FACTORY_CREATE("MaxwellConverter")->level(__VA_ARGS__)

/**
 * @brief 单位字符串到UnitType的转换函数
 */
UnitType string_to_unit_type(const std::string& unit_str);

/**
 * @brief 物理量字符串解析
 * 解析如"10mA"、"1.5V"等带单位的物理量字符串
 */
PhysicalQuantity parse_physical_quantity(const std::string& quantity_str);

/**
 * @brief 坐标系类型转换
 */
CoordinateSystemType string_to_coord_system_type(const std::string& coord_str);

/**
 * @brief 物理类型转换
 */
PhysicsType string_to_physics_type(const std::string& physics_str);

/**
 * @brief 边界条件类型转换
 */
BoundaryType string_to_boundary_type(const std::string& boundary_str);

/**
 * @brief 求解器类型转换
 */
SolverType string_to_solver_type(const std::string& solver_str);

/**
 * @brief 几何操作类型转换
 */
GeometryOperationType string_to_geometry_operation_type(const std::string& op_str);

/**
 * @brief 材料属性转换器
 */
class MaterialConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> material_node_;
    
public:
    MaterialConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : material_node_(node) {}
    
    /**
     * @brief 转换材料数据
     */
    maxwell_data::Material convert();
    
private:
    /**
     * @brief 解析材料属性
     */
    void parse_material_properties(maxwell_data::Material& material);
    
    /**
     * @brief 解析外观数据
     */
    void parse_appearance_data(maxwell_data::Material& material);
    
    /**
     * @brief 解析物理类型
     */
    void parse_physics_types(maxwell_data::Material& material);
};

/**
 * @brief 几何操作转换器
 */
class GeometryOperationConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> operation_node_;
    
public:
    GeometryOperationConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : operation_node_(node) {}
    
    /**
     * @brief 转换几何操作数据
     */
    maxwell_data::GeometryOperation convert();
    
private:
    /**
     * @brief 解析长方体参数
     */
    maxwell_data::BoxParameters parse_box_parameters();
    
    /**
     * @brief 解析操作身份信息
     */
    maxwell_data::OperationIdentity parse_operation_identity();
    
    /**
     * @brief 解析拓扑信息
     */
    maxwell_data::TopologyInfo parse_topology_info();
};

/**
 * @brief 几何部件转换器
 */
class GeometryPartConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> part_node_;
    
public:
    GeometryPartConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : part_node_(node) {}
    
    /**
     * @brief 转换几何部件数据
     */
    maxwell_data::GeometryPart convert();
    
private:
    /**
     * @brief 解析几何属性
     */
    maxwell_data::GeometryAttributes parse_attributes();
    
    /**
     * @brief 解析几何操作
     */
    std::vector<maxwell_data::GeometryOperation> parse_operations();
};

/**
 * @brief 边界条件转换器
 */
class BoundaryConditionConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> boundary_node_;
    
public:
    BoundaryConditionConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : boundary_node_(node) {}
    
    /**
     * @brief 转换边界条件数据
     */
    maxwell_data::BoundaryCondition convert();
    
private:
    /**
     * @brief 解析几何位置
     */
    maxwell_data::GeometryPosition parse_geometry_position(
        std::shared_ptr<maxwell_parser::BlockNode> pos_node);
    
    /**
     * @brief 解析坐标系向量
     */
    maxwell_data::CoordinateSystemVector parse_coord_system_vector();
    
    /**
     * @brief 解析激励源参数
     */
    void parse_excitation_parameters(maxwell_data::BoundaryCondition& boundary);
};

/**
 * @brief 求解器设置转换器
 */
class SolverSetupConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> setup_node_;
    
public:
    SolverSetupConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : setup_node_(node) {}
    
    /**
     * @brief 转换求解器设置数据
     */
    maxwell_data::SolverSetup convert();
};

/**
 * @brief 网格设置转换器
 */
class MeshSettingsConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> mesh_node_;
    
public:
    MeshSettingsConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : mesh_node_(node) {}
    
    /**
     * @brief 转换网格设置数据
     */
    maxwell_data::MeshSettings convert();
};

/**
 * @brief 设计模型转换器
 */
class DesignModelConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> model_node_;
    
public:
    DesignModelConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : model_node_(node) {}
    
    /**
     * @brief 转换设计模型数据
     */
    maxwell_data::DesignModel convert();
    
private:
    /**
     * @brief 解析温度设置
     */
    maxwell_data::TemperatureSettings parse_temperature_settings();
    
    /**
     * @brief 解析几何部件
     */
    std::vector<maxwell_data::GeometryPart> parse_geometry_parts();
    
    /**
     * @brief 解析边界条件
     */
    std::vector<maxwell_data::BoundaryCondition> parse_boundary_conditions();
    
    /**
     * @brief 解析求解结果
     */
    std::vector<maxwell_data::Solution> parse_solutions();
};

/**
 * @brief 项目元数据转换器
 */
class ProjectMetadataConverter {
private:
    std::shared_ptr<maxwell_parser::BlockNode> project_node_;
    
public:
    ProjectMetadataConverter(std::shared_ptr<maxwell_parser::BlockNode> node) 
        : project_node_(node) {}
    
    /**
     * @brief 转换项目元数据
     */
    maxwell_data::ProjectMetadata convert();
};

/**
 * @brief Maxwell数据转换器主类
 * 负责协调各个子转换器，完成完整的数据转换
 */
class MaxwellDataConverter {
private:
    std::shared_ptr<maxwell_parser::MaxwellParser> parser_;
    
public:
    MaxwellDataConverter(std::shared_ptr<maxwell_parser::MaxwellParser> parser) 
        : parser_(parser) {}
    
    /**
     * @brief 转换完整的Maxwell项目数据
     */
    maxwell_data::MaxwellProjectData convert();
    
    /**
     * @brief 从文件路径直接转换
     */
    static std::optional<maxwell_data::MaxwellProjectData> 
    convert_from_file(const std::string& file_path);
    
private:
    /**
     * @brief 转换项目元数据
     */
    maxwell_data::ProjectMetadata convert_project_metadata(
        std::shared_ptr<maxwell_parser::BlockNode> root_node);
    
    /**
     * @brief 转换材料库
     */
    std::unordered_map<std::string, maxwell_data::Material> convert_materials(
        std::shared_ptr<maxwell_parser::BlockNode> definitions_node);
    
    /**
     * @brief 转换设计模型
     */
    std::vector<maxwell_data::DesignModel> convert_design_models(
        std::shared_ptr<maxwell_parser::BlockNode> root_node);
    
    /**
     * @brief 查找指定名称的块节点
     */
    std::shared_ptr<maxwell_parser::BlockNode> find_block_node(
        std::shared_ptr<maxwell_parser::BlockNode> parent, 
        const std::string& block_name);
    
    /**
     * @brief 查找所有指定名称的子块节点
     */
    std::vector<std::shared_ptr<maxwell_parser::BlockNode>> find_all_block_nodes(
        std::shared_ptr<maxwell_parser::BlockNode> parent, 
        const std::string& block_name);
};

} // namespace maxwell_converter
} // namespace tool
} // namespace fe_em