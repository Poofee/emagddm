/**
 * @file project_data.cpp
 * @brief 项目数据结构实现文件
 * @details 实现所有数据类的序列化与反序列化功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/project_data.hpp"
#include "tool/em_enums.hpp"
#include "tool/id_generator.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace tool {

// ============================================================================
// Material类实现
// ============================================================================

Material::Material(const std::string& name) : name_(name) {
    id_ = IDGenerator::getInstance().generateID(IDCategory::MATERIAL);
}

nlohmann::json Material::toJson() const {
    nlohmann::json json;
    
    // 基础信息
    json["name"] = name_;
    json["id"] = id_;
    json["type"] = static_cast<int>(type_);
    
    // 基本属性
    json["relative_permeability"] = relative_permeability_;
    json["conductivity"] = conductivity_;
    json["mass_density"] = mass_density_;
    
    // 磁芯损耗设置
    json["core_loss_enabled"] = core_loss_enabled_;
    json["core_loss_model"] = static_cast<int>(core_loss_model_);
    json["core_loss_ks"] = core_loss_ks_;
    json["core_loss_alpha"] = core_loss_alpha_;
    json["core_loss_beta"] = core_loss_beta_;
    json["core_loss_kn"] = core_loss_kn_;
    
    // B-H曲线
    json["bh_curve_type"] = static_cast<int>(bh_curve_type_);
    nlohmann::json bh_json;
    for (const auto& point : bh_curve_) {
        nlohmann::json point_json;
        point_json["h"] = point.h;
        point_json["b"] = point.b;
        bh_json.push_back(point_json);
    }
    json["bh_curve"] = bh_json;
    
    // 自定义属性
    nlohmann::json properties_json;
    for (const auto& property : properties_) {
        nlohmann::json prop_json;
        prop_json["name"] = property.name;
        prop_json["type"] = property.type;
        prop_json["unit"] = property.unit;
        prop_json["is_temperature_dependent"] = property.is_temperature_dependent;
        
        // 处理variant类型值
        if (std::holds_alternative<double>(property.value)) {
            prop_json["value"] = std::get<double>(property.value);
        } else if (std::holds_alternative<std::vector<double>>(property.value)) {
            prop_json["value"] = std::get<std::vector<double>>(property.value);
        } else if (std::holds_alternative<std::string>(property.value)) {
            prop_json["value"] = std::get<std::string>(property.value);
        }
        
        // 温度相关数据
        nlohmann::json temp_data_json;
        for (const auto& temp_data : property.temp_dependent_data) {
            nlohmann::json temp_json;
            temp_json["temperature"] = temp_data.first;
            temp_json["values"] = temp_data.second;
            temp_data_json.push_back(temp_json);
        }
        prop_json["temp_dependent_data"] = temp_data_json;
        
        properties_json.push_back(prop_json);
    }
    json["properties"] = properties_json;
    
    // Maxwell专属数据
    json["maxwell_material_id"] = maxwell_material_id_;
    json["temperature_coefficient"] = temperature_coefficient_;
    json["bh_custom_curve_file"] = bh_custom_curve_file_;
    json["coreloss_user_data_file"] = coreloss_user_data_file_;
    json["anisotropic_permeability"] = anisotropic_permeability_;
    json["anisotropic_conductivity"] = anisotropic_conductivity_;
    json["maxwell_specific_params"] = maxwell_specific_params_;
    
    return json;
}

bool Material::fromJson(const nlohmann::json& json) {
    try {
        // 基础信息
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        type_ = static_cast<MatType>(json.value("type", 0));
        
        // 基本属性
        relative_permeability_ = json.value("relative_permeability", 1.0);
        conductivity_ = json.value("conductivity", 0.0);
        mass_density_ = json.value("mass_density", 0.0);
        
        // 磁芯损耗设置
        core_loss_enabled_ = json.value("core_loss_enabled", false);
        core_loss_model_ = static_cast<CoreLossModelType>(json.value("core_loss_model", 0));
        core_loss_ks_ = json.value("core_loss_ks", 0.0);
        core_loss_alpha_ = json.value("core_loss_alpha", 0.0);
        core_loss_beta_ = json.value("core_loss_beta", 0.0);
        core_loss_kn_ = json.value("core_loss_kn", 0.0);
        
        // B-H曲线
        bh_curve_type_ = static_cast<BHCurveType>(json.value("bh_curve_type", 0));
        bh_curve_.clear();
        if (json.contains("bh_curve")) {
            for (const auto& point_json : json["bh_curve"]) {
                BHDataPoint point;
                point.h = point_json.value("h", 0.0);
                point.b = point_json.value("b", 0.0);
                bh_curve_.push_back(point);
            }
        }
        
        // 自定义属性
        properties_.clear();
        if (json.contains("properties")) {
            for (const auto& prop_json : json["properties"]) {
                MaterialProperty property;
                property.name = prop_json.value("name", "");
                property.type = prop_json.value("type", "");
                property.unit = prop_json.value("unit", "");
                property.is_temperature_dependent = prop_json.value("is_temperature_dependent", false);
                
                // 处理variant类型值
                if (prop_json.contains("value")) {
                    if (prop_json["value"].is_number()) {
                        property.value = prop_json["value"].get<double>();
                    } else if (prop_json["value"].is_array()) {
                        property.value = prop_json["value"].get<std::vector<double>>();
                    } else if (prop_json["value"].is_string()) {
                        property.value = prop_json["value"].get<std::string>();
                    }
                }
                
                // 温度相关数据
                if (prop_json.contains("temp_dependent_data")) {
                    for (const auto& temp_json : prop_json["temp_dependent_data"]) {
                        double temperature = temp_json.value("temperature", 0.0);
                        std::vector<double> values = temp_json.value("values", std::vector<double>());
                        property.temp_dependent_data.emplace_back(temperature, values);
                    }
                }
                
                properties_.push_back(property);
            }
        }
        
        // Maxwell专属数据
        maxwell_material_id_ = json.value("maxwell_material_id", "");
        temperature_coefficient_ = json.value("temperature_coefficient", 0.0);
        bh_custom_curve_file_ = json.value("bh_custom_curve_file", "");
        coreloss_user_data_file_ = json.value("coreloss_user_data_file", "");
        anisotropic_permeability_ = json.value("anisotropic_permeability", std::vector<double>());
        anisotropic_conductivity_ = json.value("anisotropic_conductivity", std::vector<double>());
        maxwell_specific_params_ = json.value("maxwell_specific_params", std::unordered_map<std::string, std::string>());
        
        return validate();
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Material::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 序列化为JSON，然后转换为二进制
        nlohmann::json json = toJson();
        std::string json_str = json.dump();
        
        // 添加版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(version));
        
        // 添加数据长度
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size), 
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(data_size));
        
        // 添加JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Material::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取数据长度
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        
        // 解析JSON
        nlohmann::json json = nlohmann::json::parse(json_str);
        return fromJson(json);
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Material::validate() const {
    // 验证材料名称
    if (name_.empty()) return false;
    
    // 验证材料类型
    if (type_ < MatType::LINEAR_ISOTROPIC || type_ > MatType::SUPERCONDUCTOR) return false;
    
    // 验证基本属性范围
    if (relative_permeability_ < 0.0) return false;
    if (conductivity_ < 0.0) return false;
    if (mass_density_ < 0.0) return false;
    
    // 验证磁芯损耗系数
    if (core_loss_enabled_) {
        if (core_loss_ks_ < 0.0 || core_loss_alpha_ < 0.0 || 
            core_loss_beta_ < 0.0 || core_loss_kn_ < 0.0) return false;
    }
    
    // 验证B-H曲线数据
    for (const auto& point : bh_curve_) {
        if (point.h < 0.0 || point.b < 0.0) return false;
    }
    
    return true;
}

// ============================================================================
// Boundary类实现
// ============================================================================

Boundary::Boundary(const std::string& name) : name_(name) {
    id_ = IDGenerator::getInstance().generateID(IDCategory::BOUNDARY);
}

nlohmann::json Boundary::toJson() const {
    nlohmann::json json;
    
    // 基础信息
    json["name"] = name_;
    json["id"] = id_;
    json["type"] = static_cast<int>(type_);
    
    // 几何关联
    json["faces"] = faces_;
    json["edges"] = edges_;
    json["objects"] = objects_;
    
    // 边界参数
    json["impedance_value"] = impedance_value_;
    json["vector_potential"] = vector_potential_;
    json["voltage"] = voltage_;
    json["current"] = current_;
    
    // 主从边界
    json["master_name"] = master_name_;
    json["slave_name"] = slave_name_;
    
    // Maxwell专属边界数据
    json["boundary_sub_type"] = static_cast<int>(boundary_sub_type_);
    json["periodic_mapping_type"] = static_cast<int>(periodic_mapping_type_);
    json["radiation_distance"] = radiation_distance_;
    json["perfect_e_symmetry"] = perfect_e_symmetry_;
    json["perfect_h_symmetry"] = perfect_h_symmetry_;
    json["infinite_sphere_radius"] = infinite_sphere_radius_;
    json["maxwell_boundary_id"] = maxwell_boundary_id_;
    json["boundary_subdivision_params"] = boundary_subdivision_params_;
    json["maxwell_specific_params"] = maxwell_specific_params_;
    
    return json;
}

bool Boundary::fromJson(const nlohmann::json& json) {
    try {
        // 基础信息
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        type_ = static_cast<BndType>(json.value("type", 0));
        
        // 几何关联
        faces_ = json.value("faces", std::vector<std::string>());
        edges_ = json.value("edges", std::vector<std::string>());
        objects_ = json.value("objects", std::vector<std::string>());
        
        // 边界参数
        impedance_value_ = json.value("impedance_value", 0.0);
        vector_potential_ = json.value("vector_potential", 0.0);
        voltage_ = json.value("voltage", 0.0);
        current_ = json.value("current", 0.0);
        
        // 主从边界
        master_name_ = json.value("master_name", "");
        slave_name_ = json.value("slave_name", "");
        
        // Maxwell专属边界数据
        boundary_sub_type_ = static_cast<BoundarySubType>(json.value("boundary_sub_type", 0));
        periodic_mapping_type_ = static_cast<PeriodicMappingType>(json.value("periodic_mapping_type", 0));
        radiation_distance_ = json.value("radiation_distance", 0.0);
        perfect_e_symmetry_ = json.value("perfect_e_symmetry", false);
        perfect_h_symmetry_ = json.value("perfect_h_symmetry", false);
        infinite_sphere_radius_ = json.value("infinite_sphere_radius", 0.0);
        maxwell_boundary_id_ = json.value("maxwell_boundary_id", "");
        boundary_subdivision_params_ = json.value("boundary_subdivision_params", std::vector<double>());
        maxwell_specific_params_ = json.value("maxwell_specific_params", std::unordered_map<std::string, std::string>());
        
        return validate();
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Boundary::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 将对象转换为JSON
        nlohmann::json j = toJson();
        
        // 序列化JSON为字符串
        std::string json_str = j.dump();
        
        // 写入版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(uint32_t));
        
        // 写入数据长度
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size), 
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(uint32_t));
        
        // 写入JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Boundary::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取数据长度
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        
        // 解析JSON
        nlohmann::json j = nlohmann::json::parse(json_str);
        return fromJson(j);
    } catch (const std::exception& e) {
        return false;
    }
}

bool Boundary::validate() const {
    // 验证边界名称
    if (name_.empty()) return false;
    
    // 验证边界类型
    if (type_ < BndType::DIRICHLET || type_ > BndType::CONTACTS) return false;
    
    // 验证参数范围
    if (impedance_value_ < 0.0) return false;
    if (voltage_ < 0.0) return false;
    if (current_ < 0.0) return false;
    
    // 验证主从边界一致性
    if (!master_name_.empty() && slave_name_.empty()) return false;
    if (master_name_.empty() && !slave_name_.empty()) return false;
    
    return true;
}

// ============================================================================
// Excitation类实现
// ============================================================================

Excitation::Excitation(const std::string& name) : name_(name) {
    id_ = IDGenerator::getInstance().generateID(IDCategory::EXCITATION);
}

nlohmann::json Excitation::toJson() const {
    nlohmann::json json;
    
    // 基础信息
    json["name"] = name_;
    json["id"] = id_;
    json["type"] = static_cast<int>(type_);
    
    // 激励参数
    json["value"] = value_;
    json["phase"] = phase_;
    json["is_solid"] = is_solid_;
    
    // 线圈设置
    json["coil_group"] = coil_group_;
    json["connection_type"] = static_cast<int>(connection_type_);
    json["number_of_turns"] = number_of_turns_;
    
    // 多边形点
    nlohmann::json points_json;
    for (const auto& point : polygon_points_) {
        nlohmann::json point_json;
        point_json["x"] = point.first;
        point_json["y"] = point.second;
        points_json.push_back(point_json);
    }
    json["polygon_points"] = points_json;
    
    json["direction"] = direction_;
    
    // Maxwell专属激励数据
    json["waveform_type"] = static_cast<int>(waveform_type_);
    json["frequency"] = frequency_;
    json["duty_cycle"] = duty_cycle_;
    json["winding_type"] = static_cast<int>(winding_type_);
    json["motion_type"] = static_cast<int>(motion_type_);
    json["rotation_speed"] = rotation_speed_;
    json["translation_speed"] = translation_speed_;
    json["external_circuit_file"] = external_circuit_file_;
    json["custom_waveform_file"] = custom_waveform_file_;
    json["maxwell_excitation_id"] = maxwell_excitation_id_;
    json["waveform_params"] = waveform_params_;
    json["maxwell_specific_params"] = maxwell_specific_params_;
    
    return json;
}

bool Excitation::fromJson(const nlohmann::json& json) {
    try {
        // 基础信息
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        type_ = static_cast<ExcitationType>(json.value("type", 0));
        
        // 激励参数
        value_ = json.value("value", 0.0);
        phase_ = json.value("phase", 0.0);
        is_solid_ = json.value("is_solid", false);
        
        // 线圈设置
        coil_group_ = json.value("coil_group", "");
        connection_type_ = static_cast<CoilConnectionType>(json.value("connection_type", 0));
        number_of_turns_ = json.value("number_of_turns", 1);
        
        // 多边形点
        polygon_points_.clear();
        if (json.contains("polygon_points")) {
            for (const auto& point_json : json["polygon_points"]) {
                double x = point_json.value("x", 0.0);
                double y = point_json.value("y", 0.0);
                polygon_points_.emplace_back(x, y);
            }
        }
        
        direction_ = json.value("direction", 1);
        
        // Maxwell专属激励数据
        waveform_type_ = static_cast<ExcitationWaveformType>(json.value("waveform_type", 0));
        frequency_ = json.value("frequency", 0.0);
        duty_cycle_ = json.value("duty_cycle", 0.5);
        winding_type_ = static_cast<WindingType>(json.value("winding_type", 0));
        motion_type_ = static_cast<MotionType>(json.value("motion_type", 0));
        rotation_speed_ = json.value("rotation_speed", 0.0);
        translation_speed_ = json.value("translation_speed", 0.0);
        external_circuit_file_ = json.value("external_circuit_file", "");
        custom_waveform_file_ = json.value("custom_waveform_file", "");
        maxwell_excitation_id_ = json.value("maxwell_excitation_id", "");
        waveform_params_ = json.value("waveform_params", std::vector<double>());
        maxwell_specific_params_ = json.value("maxwell_specific_params", std::unordered_map<std::string, std::string>());
        
        return validate();
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Excitation::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 将对象转换为JSON
        nlohmann::json j = toJson();
        
        // 序列化JSON为字符串
        std::string json_str = j.dump();
        
        // 写入版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(uint32_t));
        
        // 写入数据长度
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size), 
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(uint32_t));
        
        // 写入JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Excitation::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取数据长度
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        
        // 解析JSON
        nlohmann::json j = nlohmann::json::parse(json_str);
        return fromJson(j);
    } catch (const std::exception& e) {
        return false;
    }
}

bool Excitation::validate() const {
    // 验证激励名称
    if (name_.empty()) return false;
    
    // 验证激励类型
    if (type_ < ExcitationType::CURRENT_DENSITY || type_ > ExcitationType::DISPLACEMENT_CURRENT) return false;
    
    // 验证参数范围
    if (value_ < 0.0) return false;
    if (phase_ < -360.0 || phase_ > 360.0) return false;
    if (number_of_turns_ < 1) return false;
    if (direction_ != -1 && direction_ != 1) return false;
    
    // 验证波形参数
    if (frequency_ < 0.0) return false;
    if (duty_cycle_ < 0.0 || duty_cycle_ > 1.0) return false;
    if (rotation_speed_ < 0.0) return false;
    if (translation_speed_ < 0.0) return false;
    
    return true;
}

// ============================================================================
// SolutionSetup类实现
// ============================================================================

SolutionSetup::SolutionSetup(const std::string& name) : name_(name) {
    id_ = IDGenerator::getInstance().generateID(IDCategory::PROJECT);
}

nlohmann::json SolutionSetup::toJson() const {
    nlohmann::json json;
    
    // 基础信息
    json["name"] = name_;
    json["id"] = id_;
    json["solution_type"] = static_cast<int>(solution_type_);
    json["solver_type"] = static_cast<int>(solver_type_);
    json["convergence_type"] = static_cast<int>(convergence_type_);
    
    // 求解参数
    json["convergence_value"] = convergence_value_;
    json["maximum_iterations"] = maximum_iterations_;
    json["frequency"] = frequency_;
    json["stator_frequency"] = stator_frequency_;
    
    // 网格设置
    json["mesh_refinement"] = mesh_refinement_;
    json["mesh_refinement_percent"] = mesh_refinement_percent_;
    json["percent_error"] = percent_error_;
    
    // Maxwell专属求解设置
    json["hpc_parallel_mode"] = static_cast<int>(hpc_parallel_mode_);
    json["hpc_solver_mode"] = static_cast<int>(hpc_solver_mode_);
    json["num_cores"] = num_cores_;
    json["domain_decomposition_type"] = static_cast<int>(domain_decomposition_type_);
    json["adaptive_mesh_refinement"] = adaptive_mesh_refinement_;
    json["adaptive_depth"] = adaptive_depth_;
    json["skin_depth_refinement"] = skin_depth_refinement_;
    json["coreloss_refinement"] = coreloss_refinement_;
    json["maxwell_solver_id"] = maxwell_solver_id_;
    json["hpc_params"] = hpc_params_;
    json["maxwell_specific_params"] = maxwell_specific_params_;
    
    return json;
}

bool SolutionSetup::fromJson(const nlohmann::json& json) {
    try {
        // 基础信息
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        solution_type_ = static_cast<SimulationType>(json.value("solution_type", 0));
        solver_type_ = static_cast<SolverType>(json.value("solver_type", 0));
        convergence_type_ = static_cast<ConvergenceType>(json.value("convergence_type", 0));
        
        // 求解参数
        convergence_value_ = json.value("convergence_value", 0.001);
        maximum_iterations_ = json.value("maximum_iterations", 100);
        frequency_ = json.value("frequency", 0.0);
        stator_frequency_ = json.value("stator_frequency", 0.0);
        
        // 网格设置
        mesh_refinement_ = json.value("mesh_refinement", false);
        mesh_refinement_percent_ = json.value("mesh_refinement_percent", 0.0);
        percent_error_ = json.value("percent_error", 1.0);
        
        // Maxwell专属求解设置
        hpc_parallel_mode_ = static_cast<HPCParallelMode>(json.value("hpc_parallel_mode", 0));
        hpc_solver_mode_ = static_cast<HPCSolverMode>(json.value("hpc_solver_mode", 0));
        num_cores_ = json.value("num_cores", 1);
        domain_decomposition_type_ = static_cast<DomainDecompositionType>(json.value("domain_decomposition_type", 0));
        adaptive_mesh_refinement_ = json.value("adaptive_mesh_refinement", false);
        adaptive_depth_ = json.value("adaptive_depth", 0);
        skin_depth_refinement_ = json.value("skin_depth_refinement", false);
        coreloss_refinement_ = json.value("coreloss_refinement", false);
        maxwell_solver_id_ = json.value("maxwell_solver_id", "");
        hpc_params_ = json.value("hpc_params", std::unordered_map<std::string, std::string>());
        maxwell_specific_params_ = json.value("maxwell_specific_params", std::unordered_map<std::string, std::string>());
        
        return validate();
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool SolutionSetup::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 将对象转换为JSON
        nlohmann::json j = toJson();
        
        // 序列化JSON为字符串
        std::string json_str = j.dump();
        
        // 写入版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(uint32_t));
        
        // 写入数据长度
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size), 
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(uint32_t));
        
        // 写入JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool SolutionSetup::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取数据长度
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        
        // 解析JSON
        nlohmann::json j = nlohmann::json::parse(json_str);
        return fromJson(j);
    } catch (const std::exception& e) {
        return false;
    }
}

bool SolutionSetup::validate() const {
    // 验证求解设置名称
    if (name_.empty()) return false;
    
    // 验证求解类型
    if (solution_type_ < SimulationType::ELECTROSTATIC || solution_type_ > SimulationType::HARMONIC) return false;
    
    // 验证求解器类型
    if (solver_type_ < SolverType::AUTO || solver_type_ > SolverType::MOONEY_RIVLIN) return false;
    
    // 验证收敛参数
    if (convergence_value_ <= 0.0) return false;
    if (maximum_iterations_ < 1) return false;
    if (frequency_ < 0.0) return false;
    if (stator_frequency_ < 0.0) return false;
    
    // 验证网格设置
    if (mesh_refinement_percent_ < 0.0 || mesh_refinement_percent_ > 100.0) return false;
    if (percent_error_ < 0.0) return false;
    
    // 验证HPC设置
    if (num_cores_ < 1) return false;
    if (adaptive_depth_ < 0) return false;
    
    return true;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

void Material::setAnisotropicPermeability(const std::vector<double>& permeability) {
    anisotropic_permeability_ = permeability;
}

void Material::setAnisotropicConductivity(const std::vector<double>& conductivity) {
    anisotropic_conductivity_ = conductivity;
}

void Material::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

void Boundary::setBoundarySubdivisionParameters(const std::vector<double>& params) {
    boundary_subdivision_params_ = params;
}

void Boundary::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

void Excitation::setWaveformParameters(const std::vector<double>& params) {
    waveform_params_ = params;
}

void Excitation::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

void SolutionSetup::setHPCParameters(const std::unordered_map<std::string, std::string>& params) {
    hpc_params_ = params;
}

void SolutionSetup::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

} // namespace tool