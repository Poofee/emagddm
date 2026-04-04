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
        (void)*reinterpret_cast<const uint32_t*>(&data[offset]);
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

    // 瞬态求解扩展字段
    json["enabled"] = enabled_;
    json["mesh_link_import"] = mesh_link_import_;
    json["time_integration_method"] = time_integration_method_;
    json["smooth_bh_curve"] = smooth_bh_curve_;
    json["output_error"] = output_error_;
    json["output_per_object_core_loss"] = output_per_object_core_loss_;
    json["output_per_object_solid_loss"] = output_per_object_solid_loss_;
    json["use_control_program"] = use_control_program_;
    json["control_program_name"] = control_program_name_;
    json["control_program_arg"] = control_program_arg_;
    json["fast_reach_steady_state"] = fast_reach_steady_state_;
    json["auto_detect_steady_state"] = auto_detect_steady_state_;
    json["frequency_of_added_voltage_source"] = frequency_of_added_voltage_source_;
    json["stop_criterion"] = stop_criterion_;
    json["is_general_transient"] = is_general_transient_;
    json["is_half_periodic_transient"] = is_half_periodic_transient_;
    json["save_fields_type"] = save_fields_type_;
    json["n_steps"] = n_steps_;
    json["steps_from"] = steps_from_;
    json["steps_to"] = steps_to_;
    json["use_nonlinear_iter_num"] = use_nonlinear_iter_num_;
    json["cache_save_kind"] = cache_save_kind_;
    json["number_solve_steps"] = number_solve_steps_;
    json["range_start"] = range_start_;
    json["range_end"] = range_end_;
    json["use_adaptive_time_step"] = use_adaptive_time_step_;
    json["initial_time_step"] = initial_time_step_;
    json["min_time_step"] = min_time_step_;
    json["max_time_step"] = max_time_step_;
    json["time_step_err_tolerance"] = time_step_err_tolerance_;
    
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

        // 瞬态求解扩展字段
        enabled_ = json.value("enabled", true);
        mesh_link_import_ = json.value("mesh_link_import", false);
        time_integration_method_ = json.value("time_integration_method", 0);
        smooth_bh_curve_ = json.value("smooth_bh_curve", false);
        output_error_ = json.value("output_error", false);
        output_per_object_core_loss_ = json.value("output_per_object_core_loss", false);
        output_per_object_solid_loss_ = json.value("output_per_object_solid_loss", false);
        use_control_program_ = json.value("use_control_program", false);
        control_program_name_ = json.value("control_program_name", "");
        control_program_arg_ = json.value("control_program_arg", "");
        fast_reach_steady_state_ = json.value("fast_reach_steady_state", true);
        auto_detect_steady_state_ = json.value("auto_detect_steady_state", true);
        frequency_of_added_voltage_source_ = json.value("frequency_of_added_voltage_source", "");
        stop_criterion_ = json.value("stop_criterion", 0.005);
        is_general_transient_ = json.value("is_general_transient", true);
        is_half_periodic_transient_ = json.value("is_half_periodic_transient", false);
        save_fields_type_ = json.value("save_fields_type", "");
        n_steps_ = json.value("n_steps", "");
        steps_from_ = json.value("steps_from", "");
        steps_to_ = json.value("steps_to", "");
        use_nonlinear_iter_num_ = json.value("use_nonlinear_iter_num", false);
        cache_save_kind_ = json.value("cache_save_kind", "");
        number_solve_steps_ = json.value("number_solve_steps", 0);
        range_start_ = json.value("range_start", "");
        range_end_ = json.value("range_end", "");
        use_adaptive_time_step_ = json.value("use_adaptive_time_step", false);
        initial_time_step_ = json.value("initial_time_step", "");
        min_time_step_ = json.value("min_time_step", "");
        max_time_step_ = json.value("max_time_step", "");
        time_step_err_tolerance_ = json.value("time_step_err_tolerance", 0.0001);
        
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
        (void)*reinterpret_cast<const uint32_t*>(&data[offset]);
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

    // 验证瞬态求解扩展字段
    if (stop_criterion_ <= 0.0 || stop_criterion_ >= 1.0) return false;
    if (time_step_err_tolerance_ <= 0.0) return false;
    if (number_solve_steps_ < 0) return false;

    // 引用所有新字段以消除编译器警告
    (void)enabled_;
    (void)mesh_link_import_;
    (void)time_integration_method_;
    (void)smooth_bh_curve_;
    (void)output_error_;
    (void)output_per_object_core_loss_;
    (void)output_per_object_solid_loss_;
    (void)use_control_program_;
    (void)control_program_name_;
    (void)control_program_arg_;
    (void)fast_reach_steady_state_;
    (void)auto_detect_steady_state_;
    (void)frequency_of_added_voltage_source_;
    (void)is_general_transient_;
    (void)is_half_periodic_transient_;
    (void)save_fields_type_;
    (void)n_steps_;
    (void)steps_from_;
    (void)steps_to_;
    (void)use_nonlinear_iter_num_;
    (void)cache_save_kind_;
    (void)range_start_;
    (void)range_end_;
    (void)use_adaptive_time_step_;
    (void)initial_time_step_;
    (void)min_time_step_;
    (void)max_time_step_;

    return true;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

// ============================================================================
// Material类新增方法实现
// ============================================================================

void Material::setLibraryInfo(const std::string& lib, const std::string& location, bool mod_since) {
    library_ = lib;
    lib_location_ = location;
    mod_since_lib_ = mod_since;
}

void Material::setAppearance(int r, int g, int b, double transparency) {
    appearance_data_ = MaterialAppearance{r, g, b, transparency};
}

void Material::addThermalModifier(const ThermalModifier& mod) {
    thermal_modifiers_.push_back(mod);
}

void Material::setTemperatureBHCurves(const std::map<double, std::vector<BHDataPoint>>& curves) {
    temperature_bh_curves_ = curves;
}

void Material::setCoreLossFreqCurves(const std::map<double, std::vector<std::pair<double, double>>>& curves) {
    core_loss_freq_curves_ = curves;
}

// ============================================================================
// Winding类实现
// ============================================================================

Winding::Winding(const std::string& name) : name_(name) {}

nlohmann::json Winding::toJson() const {
    nlohmann::json json;
    json["name"] = name_;
    json["id"] = id_;
    json["excitation_type"] = static_cast<int>(excitation_type_);
    json["is_solid"] = is_solid_;
    json["current_expression"] = current_expression_;
    json["voltage_expression"] = voltage_expression_;
    json["resistance"] = resistance_;
    json["inductance"] = inductance_;
    json["parallel_branches_num"] = parallel_branches_num_;
    json["coils"] = coils_;
    return json;
}

bool Winding::fromJson(const nlohmann::json& json) {
    try {
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        excitation_type_ = static_cast<WindingExcitationType>(json.value("excitation_type", 0));
        is_solid_ = json.value("is_solid", false);
        current_expression_ = json.value("current_expression", "");
        voltage_expression_ = json.value("voltage_expression", "");
        resistance_ = json.value("resistance", "");
        inductance_ = json.value("inductance", "");
        parallel_branches_num_ = json.value("parallel_branches_num", "");
        coils_ = json.value("coils", std::vector<std::string>());
        return validate();
    } catch (...) { return false; }
}

bool Winding::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string json_str = j.dump();
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version),
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(uint32_t));
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size),
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(uint32_t));
        data.insert(data.end(), json_str.begin(), json_str.end());
        return true;
    } catch (...) { return false; }
}

bool Winding::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        return fromJson(nlohmann::json::parse(json_str));
    } catch (...) { return false; }
}

bool Winding::validate() const { return !name_.empty(); }

void Winding::addCoil(const std::string& coil_name) { coils_.push_back(coil_name); }

// ============================================================================
// MotionSetup类实现
// ============================================================================

MotionSetup::MotionSetup(const std::string& name) : name_(name) {}

nlohmann::json MotionSetup::toJson() const {
    nlohmann::json json;
    json["name"] = name_;
    json["id"] = id_;
    json["motion_type"] = static_cast<int>(motion_type_);
    json["move_type"] = static_cast<int>(move_type_);
    json["coordinate_system"] = coordinate_system_;
    json["axis"] = static_cast<int>(axis_);
    json["is_positive"] = is_positive_;
    json["initial_position"] = initial_position_;
    json["has_rotate_limit"] = has_rotate_limit_;
    json["non_cylindrical"] = non_cylindrical_;
    json["consider_mechanical_transient"] = consider_mechanical_transient_;
    json["angular_velocity"] = angular_velocity_;
    json["linear_velocity"] = linear_velocity_;
    json["band_name_ref"] = band_name_ref_;
    json["objects"] = objects_;
    json["moving_objects"] = moving_objects_;
    return json;
}

bool MotionSetup::fromJson(const nlohmann::json& json) {
    try {
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        motion_type_ = static_cast<MotionSetupType>(json.value("motion_type", 0));
        move_type_ = static_cast<MoveType>(json.value("move_type", 0));
        coordinate_system_ = json.value("coordinate_system", 1);
        axis_ = static_cast<MotionAxis>(json.value("axis", 2));
        is_positive_ = json.value("is_positive", true);
        initial_position_ = json.value("initial_position", "");
        angular_velocity_ = json.value("angular_velocity", "");
        band_name_ref_ = json.value("band_name_ref", -1);
        objects_ = json.value("objects", std::vector<int>());
        moving_objects_ = json.value("moving_objects", std::vector<int>());
        return validate();
    } catch (...) { return false; }
}

bool MotionSetup::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string json_str = j.dump();
        uint32_t v = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&v), reinterpret_cast<const uint8_t*>(&v) + sizeof(uint32_t));
        uint32_t sz = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&sz), reinterpret_cast<const uint8_t*>(&sz) + sizeof(uint32_t));
        data.insert(data.end(), json_str.begin(), json_str.end());
        return true;
    } catch (...) { return false; }
}

bool MotionSetup::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t sz = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + sz > data.size()) return false;
        std::string s(data.begin() + offset, data.begin() + offset + sz);
        offset += sz;
        return fromJson(nlohmann::json::parse(s));
    } catch (...) { return false; }
}

bool MotionSetup::validate() const { return !name_.empty(); }

void MotionSetup::addObject(int obj_id) { objects_.push_back(obj_id); }
void MotionSetup::addMovingObject(int obj_id) { moving_objects_.push_back(obj_id); }

// ============================================================================
// MeshOperation类实现
// ============================================================================

MeshOperation::MeshOperation(const std::string& name) : name_(name) {}

nlohmann::json MeshOperation::toJson() const {
    nlohmann::json json;
    json["name"] = name_;
    json["id"] = id_;
    json["type"] = static_cast<int>(type_);
    json["enabled"] = enabled_;
    json["is_component"] = is_component_;
    json["is_global"] = is_global_;
    json["objects"] = objects_;
    json["refine_inside"] = refine_inside_;
    json["restrict_elem"] = restrict_elem_;
    json["num_max_elem"] = num_max_elem_;
    json["restrict_length"] = restrict_length_;
    json["max_length"] = max_length_;
    json["apply_to_initial_mesh"] = apply_to_initial_mesh_;
    // SurfApproxBased 特有字段
    json["surf_approx_mode"] = surf_approx_mode_;
    json["surf_dev"] = surf_dev_;
    json["normal_dev"] = normal_dev_;
    // CylindricalGap 特有字段
    json["use_band_mapping_angle"] = use_band_mapping_angle_;
    json["band_mapping_angle"] = band_mapping_angle_;
    return json;
}

bool MeshOperation::fromJson(const nlohmann::json& json) {
    try {
        name_ = json.value("name", "");
        id_ = json.value("id", 0);
        type_ = static_cast<MeshOperationType>(json.value("type", 0));
        enabled_ = json.value("enabled", true);
        is_component_ = json.value("is_component", false);
        is_global_ = json.value("is_global", false);
        objects_ = json.value("objects", std::vector<int>());
        refine_inside_ = json.value("refine_inside", true);
        restrict_elem_ = json.value("restrict_elem", false);
        num_max_elem_ = json.value("num_max_elem", "");
        restrict_length_ = json.value("restrict_length", true);
        max_length_ = json.value("max_length", "");
        apply_to_initial_mesh_ = json.value("apply_to_initial_mesh", false);
        surf_approx_mode_ = json.value("surf_approx_mode", "ManualSettings");
        surf_dev_ = json.value("surf_dev", "");
        normal_dev_ = json.value("normal_dev", "");
        use_band_mapping_angle_ = json.value("use_band_mapping_angle", false);
        band_mapping_angle_ = json.value("band_mapping_angle", "");
        return validate();
    } catch (...) { return false; }
}

bool MeshOperation::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string s = j.dump();
        uint32_t v = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&v), reinterpret_cast<const uint8_t*>(&v) + sizeof(uint32_t));
        uint32_t sz = static_cast<uint32_t>(s.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&sz), reinterpret_cast<const uint8_t*>(&sz) + sizeof(uint32_t));
        data.insert(data.end(), s.begin(), s.end());
        return true;
    } catch (...) { return false; }
}

bool MeshOperation::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t sz = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + sz > data.size()) return false;
        std::string str(data.begin() + offset, data.begin() + offset + sz);
        offset += sz;
        return fromJson(nlohmann::json::parse(str));
    } catch (...) { return false; }
}

bool MeshOperation::validate() const {
    if (name_.empty()) return false;
    // 引用所有字段以消除编译器警告
    (void)surf_dev_choice_;
    (void)normal_dev_choice_;
    (void)aspect_ratio_choice_;
    return true;
}

void MeshOperation::addObject(int obj_id) { objects_.push_back(obj_id); }

// ============================================================================
// DesignVariable类实现
// ============================================================================

DesignVariable::DesignVariable(const std::string& name, const std::string& value,
                                 const std::string& unit)
    : name_(name), value_(value), unit_(unit) {}

nlohmann::json DesignVariable::toJson() const {
    nlohmann::json json;
    json["name"] = name_;
    json["value"] = value_;
    json["unit"] = unit_;
    json["is_indexed"] = is_indexed_;
    json["expression"] = expression_;
    return json;
}

bool DesignVariable::fromJson(const nlohmann::json& json) {
    try {
        name_ = json.value("name", "");
        value_ = json.value("value", "");
        unit_ = json.value("unit", "");
        is_indexed_ = json.value("is_indexed", false);
        expression_ = json.value("expression", "");
        return validate();
    } catch (...) { return false; }
}

bool DesignVariable::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string s = j.dump();
        uint32_t v = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&v), reinterpret_cast<const uint8_t*>(&v) + sizeof(uint32_t));
        uint32_t sz = static_cast<uint32_t>(s.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&sz), reinterpret_cast<const uint8_t*>(&sz) + sizeof(uint32_t));
        data.insert(data.end(), s.begin(), s.end());
        return true;
    } catch (...) { return false; }
}

bool DesignVariable::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t sz = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + sz > data.size()) return false;
        std::string str(data.begin() + offset, data.begin() + offset + sz);
        offset += sz;
        return fromJson(nlohmann::json::parse(str));
    } catch (...) { return false; }
}

bool DesignVariable::validate() const { return !name_.empty(); }

// ============================================================================
// OutputVariable类实现
// ============================================================================

OutputVariable::OutputVariable(const std::string& name, int id, const std::string& expression,
                                 const std::string& result_unit, const std::string& display_unit)
    : name_(name), id_(id), expression_(expression), result_unit_(result_unit), display_unit_(display_unit) {}

nlohmann::json OutputVariable::toJson() const {
    nlohmann::json json;
    json["name"] = name_;
    json["id"] = id_;
    json["expression"] = expression_;
    json["result_unit"] = result_unit_;
    json["display_unit"] = display_unit_;
    return json;
}

bool OutputVariable::fromJson(const nlohmann::json& json) {
    try {
        name_ = json.value("name", "");
        id_ = json.value("id", 0);
        expression_ = json.value("expression", "");
        result_unit_ = json.value("result_unit", "");
        display_unit_ = json.value("display_unit", "");
        return validate();
    } catch (...) { return false; }
}

bool OutputVariable::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string s = j.dump();
        uint32_t v = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&v), reinterpret_cast<const uint8_t*>(&v) + sizeof(uint32_t));
        uint32_t sz = static_cast<uint32_t>(s.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&sz), reinterpret_cast<const uint8_t*>(&sz) + sizeof(uint32_t));
        data.insert(data.end(), s.begin(), s.end());
        return true;
    } catch (...) { return false; }
}

bool OutputVariable::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t sz = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + sz > data.size()) return false;
        std::string str(data.begin() + offset, data.begin() + offset + sz);
        offset += sz;
        return fromJson(nlohmann::json::parse(str));
    } catch (...) { return false; }
}

bool OutputVariable::validate() const { return !name_.empty(); }

// ============================================================================
// TemperatureSettings类实现
// ============================================================================

TemperatureSettings::TemperatureSettings() = default;

nlohmann::json TemperatureSettings::toJson() const {
    nlohmann::json json;
    json["include_temperature_dependent"] = include_temperature_dependent_;
    json["enable_feedback"] = enable_feedback_;
    nlohmann::json temp_map;
    for (const auto& [obj_id, temp_ref] : object_temperature_map_) {
        temp_map[std::to_string(obj_id)] = temp_ref;
    }
    json["object_temperatures"] = temp_map;
    return json;
}

bool TemperatureSettings::fromJson(const nlohmann::json& json) {
    try {
        include_temperature_dependent_ = json.value("include_temperature_dependent", false);
        enable_feedback_ = json.value("enable_feedback", false);
        object_temperature_map_.clear();
        if (json.contains("object_temperatures") && json["object_temperatures"].is_object()) {
            for (const auto& [key, val] : json["object_temperatures"].items()) {
                try { object_temperature_map_[std::stoi(key)] = val.get<std::string>(); } catch (...) {}
            }
        }
        return validate();
    } catch (...) { return false; }
}

bool TemperatureSettings::toBinary(std::vector<uint8_t>& data) const {
    try {
        auto j = toJson();
        std::string s = j.dump();
        uint32_t v = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&v), reinterpret_cast<const uint8_t*>(&v) + sizeof(uint32_t));
        uint32_t sz = static_cast<uint32_t>(s.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&sz), reinterpret_cast<const uint8_t*>(&sz) + sizeof(uint32_t));
        data.insert(data.end(), s.begin(), s.end());
        return true;
    } catch (...) { return false; }
}

bool TemperatureSettings::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        if (offset + sizeof(uint32_t) > data.size()) return false;
        offset += sizeof(uint32_t);
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t sz = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        if (offset + sz > data.size()) return false;
        std::string str(data.begin() + offset, data.begin() + offset + sz);
        offset += sz;
        return fromJson(nlohmann::json::parse(str));
    } catch (...) { return false; }
}

bool TemperatureSettings::validate() const { return true; }

void TemperatureSettings::setObjectTemperature(int obj_id, const std::string& temp_ref) {
    object_temperature_map_[obj_id] = temp_ref;
}





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