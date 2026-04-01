#pragma once

#include "maxwell_data_converter.hpp"
#include <regex>
#include <sstream>
#include <algorithm>

namespace fe_em {
namespace tool {
namespace maxwell_converter {

// 单位字符串映射表
static const std::unordered_map<std::string, UnitType> UNIT_MAP = {
    {"m", UnitType::METER},
    {"meter", UnitType::METER},
    {"mm", UnitType::MILLIMETER},
    {"millimeter", UnitType::MILLIMETER},
    {"cm", UnitType::CENTIMETER},
    {"centimeter", UnitType::CENTIMETER},
    {"A", UnitType::AMPERE},
    {"Ampere", UnitType::AMPERE},
    {"mA", UnitType::MILLIAMPER},
    {"milliampere", UnitType::MILLIAMPER},
    {"V", UnitType::VOLT},
    {"Volt", UnitType::VOLT},
    {"mV", UnitType::MILLIVOLT},
    {"millivolt", UnitType::MILLIVOLT},
    {"T", UnitType::TESLA},
    {"Tesla", UnitType::TESLA},
    {"A_per_meter", UnitType::AMPERE_PER_METER},
    {"S_per_meter", UnitType::SIEMENS_PER_METER},
    {"H_per_meter", UnitType::HENRY_PER_METER},
    {"F_per_meter", UnitType::FARAD_PER_METER},
    {"deg", UnitType::DEGREE},
    {"degree", UnitType::DEGREE},
    {"rad", UnitType::RADIAN},
    {"radian", UnitType::RADIAN},
    {"cel", UnitType::CELSIUS},
    {"celsius", UnitType::CELSIUS},
    {"K", UnitType::KELVIN},
    {"kelvin", UnitType::KELVIN},
    {"Pa", UnitType::PASCAL},
    {"pascal", UnitType::PASCAL},
    {"N", UnitType::NEWTON},
    {"newton", UnitType::NEWTON},
    {"W", UnitType::WATT},
    {"watt", UnitType::WATT},
    {"Hz", UnitType::HERTZ},
    {"hertz", UnitType::HERTZ}
};

// 坐标系类型映射表
static const std::unordered_map<std::string, CoordinateSystemType> COORD_SYSTEM_MAP = {
    {"Cartesian", CoordinateSystemType::CARTESIAN},
    {"Cylindrical", CoordinateSystemType::CYLINDRICAL},
    {"Spherical", CoordinateSystemType::SPHERICAL}
};

// 物理类型映射表
static const std::unordered_map<std::string, PhysicsType> PHYSICS_TYPE_MAP = {
    {"Electromagnetic", PhysicsType::ELECTROMAGNETIC},
    {"Thermal", PhysicsType::THERMAL},
    {"Structural", PhysicsType::STRUCTURAL}
};

// 边界条件类型映射表
static const std::unordered_map<std::string, BoundaryType> BOUNDARY_TYPE_MAP = {
    {"Insulating", BoundaryType::INSULATING},
    {"Zero Tangential H Field", BoundaryType::ZERO_TANGENTIAL_H_FIELD},
    {"Zero Integrated Tangential H Field", BoundaryType::INTEGRATED_ZERO_TANGENTIAL_H_FIELD},
    {"Symmetry", BoundaryType::SYMMETRY},
    {"Independent", BoundaryType::INDEPENDENT},
    {"Dependent", BoundaryType::DEPENDENT},
    {"Tangential H Field", BoundaryType::TANGENTIAL_H_FIELD},
    {"Current", BoundaryType::CURRENT},
    {"Current Density", BoundaryType::CURRENT_DENSITY},
    {"Current Density Terminal", BoundaryType::CURRENT_DENSITY_TERMINAL},
    {"Voltage", BoundaryType::VOLTAGE},
    {"Voltage Drop", BoundaryType::VOLTAGE_DROP}
};

// 求解器类型映射表
static const std::unordered_map<std::string, SolverType> SOLVER_TYPE_MAP = {
    {"Magnetostatic", SolverType::MAGNETOSTATIC},
    {"EddyCurrent", SolverType::EDDY_CURRENT},
    {"Transient", SolverType::TRANSIENT}
};

// 几何操作类型映射表
static const std::unordered_map<std::string, GeometryOperationType> GEOMETRY_OPERATION_MAP = {
    {"Box", GeometryOperationType::BOX},
    {"Cylinder", GeometryOperationType::CYLINDER},
    {"Sphere", GeometryOperationType::SPHERE},
    {"Extrude", GeometryOperationType::EXTRUDE},
    {"Revolve", GeometryOperationType::REVOLVE},
    {"Sweep", GeometryOperationType::SWEEP}
};

UnitType string_to_unit_type(const std::string& unit_str) {
    auto it = UNIT_MAP.find(unit_str);
    if (it != UNIT_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知单位类型: {}", unit_str);
    return UnitType::UNKNOWN;
}

PhysicalQuantity parse_physical_quantity(const std::string& quantity_str) {
    PhysicalQuantity quantity;
    quantity.raw_string = quantity_str;
    
    // 正则表达式匹配数值和单位
    std::regex pattern(R"(([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*([a-zA-Z_]*))");
    std::smatch matches;
    
    if (std::regex_search(quantity_str, matches, pattern)) {
        if (matches.size() >= 4) {
            // 解析数值
            try {
                quantity.value = std::stod(matches[1].str());
            } catch (const std::exception& e) {
                LOG_CONVERTER(error, "解析物理量数值失败: {}, 原始字符串: {}", e.what(), quantity_str);
                return quantity;
            }
            
            // 解析单位
            std::string unit_str = matches[3].str();
            if (!unit_str.empty()) {
                quantity.unit = string_to_unit_type(unit_str);
            }
        }
    }
    
    return quantity;
}

CoordinateSystemType string_to_coord_system_type(const std::string& coord_str) {
    auto it = COORD_SYSTEM_MAP.find(coord_str);
    if (it != COORD_SYSTEM_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知坐标系类型: {}", coord_str);
    return CoordinateSystemType::UNKNOWN;
}

PhysicsType string_to_physics_type(const std::string& physics_str) {
    auto it = PHYSICS_TYPE_MAP.find(physics_str);
    if (it != PHYSICS_TYPE_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知物理类型: {}", physics_str);
    return PhysicsType::UNKNOWN;
}

BoundaryType string_to_boundary_type(const std::string& boundary_str) {
    auto it = BOUNDARY_TYPE_MAP.find(boundary_str);
    if (it != BOUNDARY_TYPE_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知边界条件类型: {}", boundary_str);
    return BoundaryType::UNKNOWN;
}

SolverType string_to_solver_type(const std::string& solver_str) {
    auto it = SOLVER_TYPE_MAP.find(solver_str);
    if (it != SOLVER_TYPE_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知求解器类型: {}", solver_str);
    return SolverType::UNKNOWN;
}

GeometryOperationType string_to_geometry_operation_type(const std::string& op_str) {
    auto it = GEOMETRY_OPERATION_MAP.find(op_str);
    if (it != GEOMETRY_OPERATION_MAP.end()) {
        return it->second;
    }
    LOG_CONVERTER(warn, "未知几何操作类型: {}", op_str);
    return GeometryOperationType::UNKNOWN;
}

// MaterialConverter 实现
maxwell_data::Material MaterialConverter::convert() {
    maxwell_data::Material material;
    
    // 解析材料名称
    material.name = material_node_->name;
    
    // 解析坐标系类型
    auto coord_prop = material_node_->find_property("CoordinateSystemType");
    if (coord_prop) {
        if (std::holds_alternative<std::string>(coord_prop->value)) {
            material.coord_system = string_to_coord_system_type(
                std::get<std::string>(coord_prop->value));
        }
    }
    
    // 解析体/面材料类型
    auto bulk_prop = material_node_->find_property("BulkOrSurfaceType");
    if (bulk_prop) {
        if (std::holds_alternative<double>(bulk_prop->value)) {
            material.is_bulk_material = (std::get<double>(bulk_prop->value) == 1);
        }
    }
    
    // 解析物理类型
    parse_physics_types(material);
    
    // 解析材料属性
    parse_material_properties(material);
    
    // 解析外观数据
    parse_appearance_data(material);
    
    // 解析库信息
    auto lib_prop = material_node_->find_property("Library");
    if (lib_prop && std::holds_alternative<std::string>(lib_prop->value)) {
        material.library = std::get<std::string>(lib_prop->value);
    }
    
    auto loc_prop = material_node_->find_property("LibLocation");
    if (loc_prop && std::holds_alternative<std::string>(loc_prop->value)) {
        material.lib_location = std::get<std::string>(loc_prop->value);
    }
    
    auto mod_prop = material_node_->find_property("ModSinceLib");
    if (mod_prop) {
        if (std::holds_alternative<bool>(mod_prop->value)) {
            material.modified_since_lib = std::get<bool>(mod_prop->value);
        } else if (std::holds_alternative<std::string>(mod_prop->value)) {
            material.modified_since_lib = (std::get<std::string>(mod_prop->value) == "true");
        }
    }
    
    return material;
}

void MaterialConverter::parse_physics_types(maxwell_data::Material& material) {
    auto physics_node = material_node_->find_children("PhysicsTypes");
    if (!physics_node.empty()) {
        auto set_node = physics_node[0]->find_children("set");
        if (!set_node.empty()) {
            for (const auto& prop : set_node[0]->properties) {
                if (std::holds_alternative<std::string>(prop.value)) {
                    auto physics_type = string_to_physics_type(
                        std::get<std::string>(prop.value));
                    if (physics_type != PhysicsType::UNKNOWN) {
                        material.physics_types.push_back(physics_type);
                    }
                }
            }
        }
    }
}

void MaterialConverter::parse_material_properties(maxwell_data::Material& material) {
    // 解析电磁属性
    auto permittivity_prop = material_node_->find_property("permittivity");
    if (permittivity_prop && std::holds_alternative<std::string>(permittivity_prop->value)) {
        material.properties.permittivity = parse_physical_quantity(
            std::get<std::string>(permittivity_prop->value));
    }
    
    auto permeability_prop = material_node_->find_property("permeability");
    if (permeability_prop && std::holds_alternative<std::string>(permeability_prop->value)) {
        material.properties.permeability = parse_physical_quantity(
            std::get<std::string>(permeability_prop->value));
    }
    
    auto conductivity_prop = material_node_->find_property("conductivity");
    if (conductivity_prop && std::holds_alternative<std::string>(conductivity_prop->value)) {
        material.properties.conductivity = parse_physical_quantity(
            std::get<std::string>(conductivity_prop->value));
    }
    
    // 解析矫顽力（向量属性）
    auto coercivity_node = material_node_->find_children("magnetic_coercivity");
    if (!coercivity_node.empty()) {
        maxwell_data::Vector3D coercivity = {0, 0, 0};
        
        auto mag_prop = coercivity_node[0]->find_property("Magnitude");
        if (mag_prop && std::holds_alternative<std::string>(mag_prop->value)) {
            auto magnitude = parse_physical_quantity(std::get<std::string>(mag_prop->value));
            
            auto dir1_prop = coercivity_node[0]->find_property("DirComp1");
            auto dir2_prop = coercivity_node[0]->find_property("DirComp2");
            auto dir3_prop = coercivity_node[0]->find_property("DirComp3");
            
            double dir1 = 0, dir2 = 0, dir3 = 0;
            if (dir1_prop && std::holds_alternative<double>(dir1_prop->value)) {
                dir1 = std::get<double>(dir1_prop->value);
            }
            if (dir2_prop && std::holds_alternative<double>(dir2_prop->value)) {
                dir2 = std::get<double>(dir2_prop->value);
            }
            if (dir3_prop && std::holds_alternative<double>(dir3_prop->value)) {
                dir3 = std::get<double>(dir3_prop->value);
            }
            
            coercivity[0] = magnitude.value * dir1;
            coercivity[1] = magnitude.value * dir2;
            coercivity[2] = magnitude.value * dir3;
            
            material.properties.magnetic_coercivity = coercivity;
        }
    }
    
    // 解析热学属性（类似方式）
    // ... 其他属性的解析实现
}

void MaterialConverter::parse_appearance_data(maxwell_data::Material& material) {
    auto appearance_node = material_node_->find_children("MatAppearanceData");
    if (!appearance_node.empty()) {
        maxwell_data::MaterialAppearance appearance;
        
        auto red_prop = appearance_node[0]->find_property("Red");
        auto green_prop = appearance_node[0]->find_property("Green");
        auto blue_prop = appearance_node[0]->find_property("Blue");
        auto trans_prop = appearance_node[0]->find_property("Transparency");
        
        if (red_prop && std::holds_alternative<double>(red_prop->value)) {
            appearance.color[0] = static_cast<int>(std::get<double>(red_prop->value));
        }
        if (green_prop && std::holds_alternative<double>(green_prop->value)) {
            appearance.color[1] = static_cast<int>(std::get<double>(green_prop->value));
        }
        if (blue_prop && std::holds_alternative<double>(blue_prop->value)) {
            appearance.color[2] = static_cast<int>(std::get<double>(blue_prop->value));
        }
        if (trans_prop && std::holds_alternative<double>(trans_prop->value)) {
            appearance.transparency = std::get<double>(trans_prop->value);
        }
        
        material.properties.appearance = appearance;
    }
}

// GeometryOperationConverter 实现
maxwell_data::GeometryOperation GeometryOperationConverter::convert() {
    maxwell_data::GeometryOperation operation;
    
    // 解析操作类型
    auto type_prop = operation_node_->find_property("OperationType");
    if (type_prop && std::holds_alternative<std::string>(type_prop->value)) {
        operation.type = string_to_geometry_operation_type(
            std::get<std::string>(type_prop->value));
    }
    
    // 解析操作ID
    auto id_prop = operation_node_->find_property("ID");
    if (id_prop && std::holds_alternative<double>(id_prop->value)) {
        operation.id = static_cast<int>(std::get<double>(id_prop->value));
    }
    
    // 解析参考坐标系ID
    auto ref_prop = operation_node_->find_property("ReferenceCoordSystemID");
    if (ref_prop && std::holds_alternative<double>(ref_prop->value)) {
        operation.reference_coord_system_id = static_cast<int>(std::get<double>(ref_prop->value));
    }
    
    // 解析是否抑制
    auto suppress_prop = operation_node_->find_property("IsSuppressed");
    if (suppress_prop) {
        if (std::holds_alternative<bool>(suppress_prop->value)) {
            operation.is_suppressed = std::get<bool>(suppress_prop->value);
        } else if (std::holds_alternative<std::string>(suppress_prop->value)) {
            operation.is_suppressed = (std::get<std::string>(suppress_prop->value) == "true");
        }
    }
    
    // 根据操作类型解析参数
    switch (operation.type) {
        case GeometryOperationType::BOX:
            operation.parameters = parse_box_parameters();
            break;
        default:
            LOG_CONVERTER(warn, "暂不支持的操作类型: {}", static_cast<int>(operation.type));
            break;
    }
    
    // 解析操作身份信息
    auto identity_node = operation_node_->find_children("OperationIdentity");
    if (!identity_node.empty()) {
        operation.identity = parse_operation_identity();
    }
    
    return operation;
}

maxwell_data::BoxParameters GeometryOperationConverter::parse_box_parameters() {
    maxwell_data::BoxParameters params;
    
    auto box_node = operation_node_->find_children("BoxParameters");
    if (!box_node.empty()) {
        auto kernel_prop = box_node[0]->find_property("KernelVersion");
        if (kernel_prop && std::holds_alternative<double>(kernel_prop->value)) {
            params.kernel_version = static_cast<int>(std::get<double>(kernel_prop->value));
        }
        
        // 解析位置和尺寸参数
        // ... 实现具体的参数解析
    }
    
    return params;
}

maxwell_data::OperationIdentity GeometryOperationConverter::parse_operation_identity() {
    maxwell_data::OperationIdentity identity;
    
    auto identity_node = operation_node_->find_children("OperationIdentity");
    if (!identity_node.empty()) {
        // 解析拓扑信息
        auto topology_node = identity_node[0]->find_children("Topology");
        if (!topology_node.empty()) {
            identity.topology = parse_topology_info();
        }
        
        // 解析其他身份信息
        // ... 实现具体的身份信息解析
    }
    
    return identity;
}

maxwell_data::TopologyInfo GeometryOperationConverter::parse_topology_info() {
    maxwell_data::TopologyInfo topology;
    
    auto topology_node = operation_node_->find_children("Topology");
    if (!topology_node.empty()) {
        // 解析拓扑参数
        // ... 实现具体的拓扑信息解析
    }
    
    return topology;
}

// 其他转换器的实现类似，这里省略详细实现
// 实际项目中需要完整实现所有转换器

// MaxwellDataConverter 主类实现
maxwell_data::MaxwellProjectData MaxwellDataConverter::convert() {
    maxwell_data::MaxwellProjectData project_data;
    
    auto root_node = parser_->get_root_node();
    if (!root_node) {
        LOG_CONVERTER(error, "解析器未成功解析文件");
        return project_data;
    }
    
    // 转换项目元数据
    project_data.metadata = convert_project_metadata(root_node);
    
    // 转换材料库
    auto definitions_node = find_block_node(root_node, "Definitions");
    if (definitions_node) {
        project_data.materials = convert_materials(definitions_node);
    }
    
    // 转换设计模型
    project_data.design_models = convert_design_models(root_node);
    
    LOG_CONVERTER(info, "成功转换Maxwell项目数据，包含{}个材料，{}个设计模型", 
                 project_data.materials.size(), project_data.design_models.size());
    
    return project_data;
}

std::optional<maxwell_data::MaxwellProjectData> 
MaxwellDataConverter::convert_from_file(const std::string& file_path) {
    try {
        auto parser = std::make_shared<maxwell_parser::MaxwellParser>();
        parser->parse_file(file_path);
        
        MaxwellDataConverter converter(parser);
        return converter.convert();
    } catch (const std::exception& e) {
        LOG_CONVERTER(error, "转换文件失败: {}, 文件路径: {}", e.what(), file_path);
        return std::nullopt;
    }
}

maxwell_data::ProjectMetadata MaxwellDataConverter::convert_project_metadata(
    std::shared_ptr<maxwell_parser::BlockNode> root_node) {
    maxwell_data::ProjectMetadata metadata;
    
    // 解析项目元数据属性
    auto created_prop = root_node->find_property("Created");
    if (created_prop && std::holds_alternative<std::string>(created_prop->value)) {
        metadata.created_time = std::get<std::string>(created_prop->value);
    }
    
    auto product_prop = root_node->find_property("Product");
    if (product_prop && std::holds_alternative<std::string>(product_prop->value)) {
        metadata.product = std::get<std::string>(product_prop->value);
    }
    
    // 解析其他元数据属性
    // ... 实现具体的元数据解析
    
    return metadata;
}

std::unordered_map<std::string, maxwell_data::Material> 
MaxwellDataConverter::convert_materials(
    std::shared_ptr<maxwell_parser::BlockNode> definitions_node) {
    std::unordered_map<std::string, maxwell_data::Material> materials;
    
    auto materials_node = find_block_node(definitions_node, "Materials");
    if (materials_node) {
        for (const auto& child : materials_node->children) {
            MaterialConverter converter(child);
            auto material = converter.convert();
            materials[material.name] = material;
        }
    }
    
    return materials;
}

std::vector<maxwell_data::DesignModel> MaxwellDataConverter::convert_design_models(
    std::shared_ptr<maxwell_parser::BlockNode> root_node) {
    std::vector<maxwell_data::DesignModel> design_models;
    
    // 查找所有设计模型节点
    auto maxwell3d_nodes = find_all_block_nodes(root_node, "Maxwell3DModel");
    auto maxwell2d_nodes = find_all_block_nodes(root_node, "Maxwell2DModel");
    
    // 转换3D设计模型
    for (const auto& node : maxwell3d_nodes) {
        DesignModelConverter converter(node);
        design_models.push_back(converter.convert());
    }
    
    // 转换2D设计模型
    for (const auto& node : maxwell2d_nodes) {
        DesignModelConverter converter(node);
        design_models.push_back(converter.convert());
    }
    
    return design_models;
}

std::shared_ptr<maxwell_parser::BlockNode> MaxwellDataConverter::find_block_node(
    std::shared_ptr<maxwell_parser::BlockNode> parent, 
    const std::string& block_name) {
    auto children = parent->find_children(block_name);
    if (!children.empty()) {
        return children[0];
    }
    return nullptr;
}

std::vector<std::shared_ptr<maxwell_parser::BlockNode>> 
MaxwellDataConverter::find_all_block_nodes(
    std::shared_ptr<maxwell_parser::BlockNode> parent, 
    const std::string& block_name) {
    return parent->find_children(block_name);
}

} // namespace maxwell_converter
} // namespace tool
} // namespace fe_em