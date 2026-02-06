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





// Excitation类缺失的方法实现
void Excitation::setWaveformParameters(const std::vector<double>& params) {
    waveform_params_ = params;
}

void Excitation::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

// SolutionSetup类缺失的方法实现
void SolutionSetup::setHPCParameters(const std::unordered_map<std::string, std::string>& params) {
    hpc_params_ = params;
}

void SolutionSetup::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

} // namespace tool