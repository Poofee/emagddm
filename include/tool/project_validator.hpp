/**
 * @file project_validator.hpp
 * @brief 基础工具层 - 项目数据校验器头文件
 * @details 实现数据合法性、业务规则、关联关系的全面校验
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include "tool/project_manager.hpp"
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <memory>
#include <chrono>

namespace tool {

struct ValidationError;

using ValidationCallback = std::function<void(const ValidationError&)>;

struct ValidationError {
    std::string error_code;
    std::string error_message;
    std::string data_type;
    std::string entity_id;
    std::string field_name;
    
    enum class SeverityLevel {
        INFO,
        WARNING,
        ERROR,
        FATAL
    };
    
    SeverityLevel severity = SeverityLevel::ERROR;
    
    std::string to_string() const;
};

class ProjectValidator {
public:
    ProjectValidator();
    ~ProjectValidator() = default;
    
    ValidationResult validate_project(const ProjectManager& project);
    
    ValidationResult validate_material(const Material& material);
    ValidationResult validate_geometry(const Geometry& geometry);
    ValidationResult validate_boundary(const Boundary& boundary);
    ValidationResult validate_excitation(const Excitation& excitation);
    ValidationResult validate_mesh(const Mesh& mesh);
    ValidationResult validate_solution_setup(const SolutionSetup& setup);
    
    void set_strict_mode(bool strict);
    void set_warning_as_error(bool treat_warnings_as_errors);
    void set_validation_callback(ValidationCallback callback);
    
    static std::string format_validation_error(const ValidationError& error);
    
private:
    bool strict_mode_ = false;
    bool warnings_as_errors_ = false;
    ValidationCallback callback_;
    
    ValidationResult validate_material_properties(const Material& material);
    ValidationResult validate_bh_curve(const Material& material);
    ValidationResult validate_core_loss_parameters(const Material& material);
    ValidationResult validate_boundary_association(const Boundary& boundary, 
                                                   const std::unordered_map<std::string, GeometryPtr>& geometries);
    ValidationResult validate_excitation_association(const Excitation& excitation,
                                                     const std::unordered_map<std::string, GeometryPtr>& geometries);
    ValidationResult validate_geometry_topology(const Geometry& geometry);
    void collect_errors(ValidationResult& result, const std::vector<ValidationError>& errors);
};

class ValidationResult {
public:
    ValidationResult();
    ~ValidationResult() = default;
    
    bool is_valid() const { return errors_.empty(); }
    
    void add_error(const ValidationError& error);
    void add_warning(const ValidationError& warning);
    void add_info(const ValidationError& info);
    
    const std::vector<ValidationError>& get_errors() const { return errors_; }
    const std::vector<ValidationError>& get_warnings() const { return warnings_; }
    const std::vector<ValidationError>& get_info() const { return infos_; }
    
    size_t get_error_count() const { return errors_.size(); }
    size_t get_warning_count() const { return warnings_.size(); }
    size_t get_info_count() const { return infos_.size(); }
    
    std::string get_summary() const;
    
    void merge(const ValidationResult& other);
    
private:
    std::vector<ValidationError> errors_;
    std::vector<ValidationError> warnings_;
    std::vector<ValidationError> infos_;
};

class ValidationReport {
public:
    ValidationReport() = default;
    ~ValidationReport() = default;
    
    void add_result(const std::string& data_type, const ValidationResult& result);
    
    void set_validation_time(const std::chrono::system_clock::time_point& time);
    std::chrono::system_clock::time_point get_validation_time() const;
    
    std::string to_string() const;
    bool save_to_file(const std::string& file_path) const;
    
private:
    std::unordered_map<std::string, ValidationResult> results_;
    std::chrono::system_clock::time_point validation_time_;
};

} // namespace tool
