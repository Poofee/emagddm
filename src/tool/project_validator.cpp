/**
 * @file project_validator.cpp
 * @brief 基础工具层 - 项目数据校验器源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/project_validator.hpp"
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace tool {

ValidationError::SeverityLevel string_to_severity(const std::string& str) {
    if (str == "INFO") return ValidationError::SeverityLevel::INFO;
    if (str == "WARNING") return ValidationError::SeverityLevel::WARNING;
    if (str == "ERROR") return ValidationError::SeverityLevel::ERROR;
    if (str == "FATAL") return ValidationError::SeverityLevel::FATAL;
    return ValidationError::SeverityLevel::ERROR;
}

std::string ValidationError::to_string() const {
    std::ostringstream oss;
    oss << "[" << error_code << "] ";
    switch (severity) {
        case SeverityLevel::INFO: oss << "INFO: "; break;
        case SeverityLevel::WARNING: oss << "WARNING: "; break;
        case SeverityLevel::ERROR: oss << "ERROR: "; break;
        case SeverityLevel::FATAL: oss << "FATAL: "; break;
    }
    oss << error_message;
    if (!data_type.empty()) oss << " [Type: " << data_type << "]";
    if (!entity_id.empty()) oss << " [ID: " << entity_id << "]";
    if (!field_name.empty()) oss << " [Field: " << field_name << "]";
    return oss.str();
}

std::string ProjectValidator::format_validation_error(const ValidationError& error) {
    return error.to_string();
}

ProjectValidator::ProjectValidator() = default;

ValidationResult ProjectValidator::validate_project(const ProjectData& project) {
    ValidationResult result;
    
    auto& materials = const_cast<ProjectData&>(project).getAllMaterials();
    auto& geometries = const_cast<ProjectData&>(project).getAllGeometries();
    auto& boundaries = const_cast<ProjectData&>(project).getAllBoundaries();
    auto& excitations = const_cast<ProjectData&>(project).getAllExcitations();
    
    for (const auto& material : materials) {
        ValidationResult material_result = validate_material(*material.second);
        result.merge(material_result);
    }
    
    for (const auto& geometry : geometries) {
        ValidationResult geometry_result = validate_geometry(*geometry.second);
        result.merge(geometry_result);
        
        ValidationResult topology_result = validate_geometry_topology(*geometry.second);
        result.merge(topology_result);
    }
    
    for (const auto& boundary : boundaries) {
        ValidationResult boundary_result = validate_boundary(*boundary.second);
        result.merge(boundary_result);
        
        ValidationResult association_result = validate_boundary_association(*boundary.second, geometries);
        result.merge(association_result);
    }
    
    for (const auto& excitation : excitations) {
        ValidationResult excitation_result = validate_excitation(*excitation.second);
        result.merge(excitation_result);
        
        ValidationResult association_result = validate_excitation_association(*excitation.second, geometries);
        result.merge(association_result);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_material(const Material& material) {
    ValidationResult result;
    
    if (material.getName().empty()) {
        ValidationError error;
        error.error_code = "MAT_001";
        error.error_message = "材料名称不能为空";
        error.data_type = "Material";
        error.entity_id = std::to_string(material.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (material.getID() == 0) {
        ValidationError error;
        error.error_code = "MAT_002";
        error.error_message = "材料ID无效（必须大于0）";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "id";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    ValidationResult properties_result = validate_material_properties(material);
    result.merge(properties_result);
    
    if (material.isCoreLossEnabled()) {
        ValidationResult core_loss_result = validate_core_loss_parameters(material);
        result.merge(core_loss_result);
    }
    
    auto bh_curve = material.getBHCurve();
    if (!bh_curve.empty()) {
        ValidationResult bh_result = validate_bh_curve(material);
        result.merge(bh_result);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_material_properties(const Material& material) {
    ValidationResult result;
    
    double mu_r = material.getRelativePermeability();
    if (mu_r <= 0) {
        ValidationError error;
        error.error_code = "MAT_P_001";
        error.error_message = "相对磁导率必须大于0";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "relative_permeability";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    } else if (mu_r < 0.001 || mu_r > 1000000) {
        ValidationError warning;
        warning.error_code = "MAT_P_W001";
        warning.error_message = "相对磁导率超出常规范围（0.001~1,000,000），请确认输入正确";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "relative_permeability";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    double sigma = material.getConductivity();
    if (sigma < 0) {
        ValidationError error;
        error.error_code = "MAT_P_002";
        error.error_message = "电导率不能为负数";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "conductivity";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    } else if (sigma > 1e10) {
        ValidationError warning;
        warning.error_code = "MAT_P_W002";
        warning.error_message = "电导率超出常规范围（>1e10 S/m），请确认输入正确";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "conductivity";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    double rho = material.getMassDensity();
    if (rho < 0) {
        ValidationError error;
        error.error_code = "MAT_P_003";
        error.error_message = "质量密度不能为负数";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "mass_density";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    } else if (rho > 50000) {
        ValidationError warning;
        warning.error_code = "MAT_P_W003";
        warning.error_message = "质量密度超出常规范围（>50,000 kg/m³），请确认输入正确";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "mass_density";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_bh_curve(const Material& material) {
    ValidationResult result;
    
    auto bh_curve = material.getBHCurve();
    
    if (bh_curve.empty()) {
        return result;
    }
    
    if (bh_curve.size() < 2) {
        ValidationError error;
        error.error_code = "MAT_BH_001";
        error.error_message = "BH曲线至少需要2个数据点";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "bh_curve";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
        return result;
    }
    
    double prev_h = bh_curve[0].h;
    double prev_b = bh_curve[0].b;
    
    for (size_t i = 1; i < bh_curve.size(); ++i) {
        if (bh_curve[i].h < prev_h) {
            ValidationError error;
            error.error_code = "MAT_BH_002";
            error.error_message = "BH曲线的H值必须单调递增（第" + std::to_string(i) + "点不满足）";
            error.data_type = "Material";
            error.entity_id = material.getName();
            error.field_name = "bh_curve";
            error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(error);
        }
        
        if (bh_curve[i].b < prev_b && material.getBHCurveType() == BHCurveType::LINEAR) {
            ValidationError warning;
            warning.error_code = "MAT_BH_W001";
            warning.error_message = "BH曲线的B值出现下降，可能导致求解不收敛";
            warning.data_type = "Material";
            warning.entity_id = material.getName();
            warning.field_name = "bh_curve";
            warning.severity = ValidationError::SeverityLevel::WARNING;
            result.add_warning(warning);
        }
        
        if (bh_curve[i].b <= 0) {
            ValidationError warning;
            warning.error_code = "MAT_BH_W002";
            warning.error_message = "BH曲线包含非正值B值";
            warning.data_type = "Material";
            warning.entity_id = material.getName();
            warning.field_name = "bh_curve";
            warning.severity = ValidationError::SeverityLevel::WARNING;
            result.add_warning(warning);
        }
        
        prev_h = bh_curve[i].h;
        prev_b = bh_curve[i].b;
    }
    
    if (bh_curve.front().h < 0) {
        ValidationError warning;
        warning.error_code = "MAT_BH_W003";
        warning.error_message = "BH曲线首点H值为负，可能导致初始磁导率计算异常";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "bh_curve";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_core_loss_parameters(const Material& material) {
    ValidationResult result;
    
    double ks, alpha, beta, kn;
    material.getCoreLossCoefficients(ks, alpha, beta, kn);
    
    if (ks < 0) {
        ValidationError error;
        error.error_code = "MAT_CL_001";
        error.error_message = "铁损系数Ks不能为负";
        error.data_type = "Material";
        error.entity_id = material.getName();
        error.field_name = "core_loss_coefficients";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (alpha < 0 || alpha > 3) {
        ValidationError warning;
        warning.error_code = "MAT_CL_W001";
        warning.error_message = "铁损系数Alpha超出常规范围（0~3），请确认输入正确";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "core_loss_coefficients";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    if (beta < 0 || beta > 4) {
        ValidationError warning;
        warning.error_code = "MAT_CL_W002";
        warning.error_message = "铁损系数Beta超出常规范围（0~4），请确认输入正确";
        warning.data_type = "Material";
        warning.entity_id = material.getName();
        warning.field_name = "core_loss_coefficients";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_geometry(const Geometry& geometry) {
    ValidationResult result;
    
    if (geometry.getName().empty()) {
        ValidationError error;
        error.error_code = "GEO_001";
        error.error_message = "几何名称不能为空";
        error.data_type = "Geometry";
        error.entity_id = std::to_string(geometry.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (geometry.getID() == 0) {
        ValidationError error;
        error.error_code = "GEO_002";
        error.error_message = "几何ID无效（必须大于0）";
        error.data_type = "Geometry";
        error.entity_id = geometry.getName();
        error.field_name = "id";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    auto& file_path = geometry.getFilePath();
    if (!file_path.empty()) {
        if (file_path.length() > 260) {
            ValidationError warning;
            warning.error_code = "GEO_W001";
            warning.error_message = "文件路径过长，可能在某些系统上无法访问";
            warning.data_type = "Geometry";
            warning.entity_id = geometry.getName();
            warning.field_name = "file_path";
            warning.severity = ValidationError::SeverityLevel::WARNING;
            result.add_warning(warning);
        }
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_geometry_topology(const Geometry& geometry) {
    ValidationResult result;
    
    auto& object_material_map = geometry.getObjectMaterialMap();
    std::unordered_set<std::string> defined_materials;
    
    for (const auto& obj : object_material_map) {
        if (obj.first.empty()) {
            ValidationError error;
            error.error_code = "GEO_TOPO_001";
            error.error_message = "发现空名称的几何对象";
            error.data_type = "Geometry";
            error.entity_id = geometry.getName();
            error.field_name = "object_material_map";
            error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(error);
        }
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_boundary(const Boundary& boundary) {
    ValidationResult result;
    
    if (boundary.getName().empty()) {
        ValidationError error;
        error.error_code = "BND_001";
        error.error_message = "边界名称不能为空";
        error.data_type = "Boundary";
        error.entity_id = std::to_string(boundary.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    auto type = boundary.getType();
    switch (type) {
        case BndType::IMPEDANCE:
            if (boundary.getImpedanceValue() <= 0) {
                ValidationError error;
                error.error_code = "BND_002";
                error.error_message = "阻抗边界阻抗值必须大于0";
                error.data_type = "Boundary";
                error.entity_id = boundary.getName();
                error.field_name = "impedance_value";
                error.severity = ValidationError::SeverityLevel::ERROR;
                result.add_error(error);
            }
            break;
        case BndType::VECTOR_POTENTIAL:
            if (std::abs(boundary.getVectorPotential()) > 1e6) {
                ValidationError warning;
                warning.error_code = "BND_W001";
                warning.error_message = "矢量位值超出常规范围（|Az| > 1e6），请确认输入正确";
                warning.data_type = "Boundary";
                warning.entity_id = boundary.getName();
                warning.field_name = "vector_potential";
                warning.severity = ValidationError::SeverityLevel::WARNING;
                result.add_warning(warning);
            }
            break;
        default:
            break;
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_boundary_association(
    const Boundary& boundary,
    const std::unordered_map<std::string, GeometryPtr>& geometries) {
    
    ValidationResult result;
    
    for (const auto& obj : boundary.getObjects()) {
        if (geometries.find(obj) == geometries.end()) {
            ValidationError error;
            error.error_code = "BND_ASSOC_001";
            error.error_message = "边界绑定的几何对象不存在: " + obj;
            error.data_type = "Boundary";
            error.entity_id = boundary.getName();
            error.field_name = "objects";
            error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(error);
        }
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_excitation(const Excitation& excitation) {
    ValidationResult result;
    
    if (excitation.getName().empty()) {
        ValidationError error;
        error.error_code = "EXT_001";
        error.error_message = "激励名称不能为空";
        error.data_type = "Excitation";
        error.entity_id = std::to_string(excitation.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (excitation.getID() == 0) {
        ValidationError error;
        error.error_code = "EXT_002";
        error.error_message = "激励ID无效（必须大于0）";
        error.data_type = "Excitation";
        error.entity_id = excitation.getName();
        error.field_name = "id";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (excitation.getValue() > 1e12) {
        ValidationError warning;
        warning.error_code = "EXT_W001";
        warning.error_message = "激励值超出常规范围（>1e12），请确认输入正确";
        warning.data_type = "Excitation";
        warning.entity_id = excitation.getName();
        warning.field_name = "value";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    if (excitation.getPhase() < -180 || excitation.getPhase() > 180) {
        ValidationError error;
        error.error_code = "EXT_003";
        error.error_message = "激励相位超出范围（-180° ~ 180°）";
        error.data_type = "Excitation";
        error.entity_id = excitation.getName();
        error.field_name = "phase";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (excitation.getNumberOfTurns() <= 0) {
        ValidationError error;
        error.error_code = "EXT_004";
        error.error_message = "线圈匝数必须大于0";
        error.data_type = "Excitation";
        error.entity_id = excitation.getName();
        error.field_name = "number_of_turns";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_excitation_association(
    const Excitation& excitation,
    const std::unordered_map<std::string, GeometryPtr>& geometries) {
    
    ValidationResult result;
    
    const auto& coil_group = excitation.getCoilGroup();
    if (!coil_group.empty()) {
        bool found = false;
        for (const auto& geo : geometries) {
            auto& objects = geo.second->getObjectMaterialMap();
            if (objects.find(coil_group) != objects.end()) {
                found = true;
                break;
            }
        }
        
        if (!found) {
            ValidationError warning;
            warning.error_code = "EXT_ASSOC_W001";
            warning.error_message = "激励绑定的线圈组在几何中未找到: " + coil_group;
            warning.data_type = "Excitation";
            warning.entity_id = excitation.getName();
            warning.field_name = "coil_group";
            warning.severity = ValidationError::SeverityLevel::WARNING;
            result.add_warning(warning);
        }
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_mesh(const Mesh& mesh) {
    ValidationResult result;
    
    if (mesh.getName().empty()) {
        ValidationError error;
        error.error_code = "MSH_001";
        error.error_message = "网格配置名称不能为空";
        error.data_type = "Mesh";
        error.entity_id = std::to_string(mesh.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (mesh.getMaxElementSize() <= 0) {
        ValidationError error;
        error.error_code = "MSH_002";
        error.error_message = "最大单元尺寸必须大于0";
        error.data_type = "Mesh";
        error.entity_id = mesh.getName();
        error.field_name = "max_element_size";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (mesh.getMinElementSize() >= mesh.getMaxElementSize() && mesh.getMinElementSize() > 0) {
        ValidationError warning;
        warning.error_code = "MSH_W001";
        warning.error_message = "最小单元尺寸大于等于最大单元尺寸，可能导致网格生成问题";
        warning.data_type = "Mesh";
        warning.entity_id = mesh.getName();
        warning.field_name = "min_element_size";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    if (mesh.isBoundaryLayerEnabled()) {
        if (mesh.getBoundaryLayerNumberOfLayers() <= 0) {
            ValidationError error;
            error.error_code = "MSH_003";
            error.error_message = "边界层数量必须大于0";
            error.data_type = "Mesh";
            error.entity_id = mesh.getName();
            error.field_name = "boundary_layer_number_of_layers";
            error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(error);
        }
        
        if (mesh.getBoundaryLayerThickness() <= 0) {
            ValidationError error;
            error.error_code = "MSH_004";
            error.error_message = "边界层厚度必须大于0";
            error.data_type = "Mesh";
            error.entity_id = mesh.getName();
            error.field_name = "boundary_layer_thickness";
            error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(error);
        }
    }
    
    return result;
}

ValidationResult ProjectValidator::validate_solution_setup(const SolutionSetup& setup) {
    ValidationResult result;
    
    if (setup.getName().empty()) {
        ValidationError error;
        error.error_code = "SOL_001";
        error.error_message = "求解设置名称不能为空";
        error.data_type = "SolutionSetup";
        error.entity_id = std::to_string(setup.getID());
        error.field_name = "name";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (setup.getConvergenceValue() <= 0) {
        ValidationError error;
        error.error_code = "SOL_002";
        error.error_message = "收敛阈值必须大于0";
        error.data_type = "SolutionSetup";
        error.entity_id = setup.getName();
        error.field_name = "convergence_value";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (setup.getMaximumIterations() <= 0) {
        ValidationError error;
        error.error_code = "SOL_003";
        error.error_message = "最大迭代次数必须大于0";
        error.data_type = "SolutionSetup";
        error.entity_id = setup.getName();
        error.field_name = "maximum_iterations";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (setup.getFrequency() < 0) {
        ValidationError error;
        error.error_code = "SOL_004";
        error.error_message = "频率不能为负";
        error.data_type = "SolutionSetup";
        error.entity_id = setup.getName();
        error.field_name = "frequency";
        error.severity = ValidationError::SeverityLevel::ERROR;
        result.add_error(error);
    }
    
    if (setup.getFrequency() > 1e12) {
        ValidationError warning;
        warning.error_code = "SOL_W001";
        warning.error_message = "频率超出常规范围（>1e12 Hz），请确认输入正确";
        warning.data_type = "SolutionSetup";
        warning.entity_id = setup.getName();
        warning.field_name = "frequency";
        warning.severity = ValidationError::SeverityLevel::WARNING;
        result.add_warning(warning);
    }
    
    return result;
}

void ProjectValidator::set_strict_mode(bool strict) {
    strict_mode_ = strict;
}

void ProjectValidator::set_warning_as_error(bool treat_warnings_as_errors) {
    warnings_as_errors_ = treat_warnings_as_errors;
}

void ProjectValidator::set_validation_callback(ValidationCallback callback) {
    callback_ = callback;
}

void ProjectValidator::collect_errors(ValidationResult& result, const std::vector<ValidationError>& errors) {
    for (const auto& error : errors) {
        if (warnings_as_errors_ && error.severity == ValidationError::SeverityLevel::WARNING) {
            ValidationError modified_error = error;
            modified_error.severity = ValidationError::SeverityLevel::ERROR;
            result.add_error(modified_error);
        } else {
            result.add_error(error);
        }
    }
}

ValidationResult::ValidationResult() = default;

void ValidationResult::add_error(const ValidationError& error) {
    errors_.push_back(error);
}

void ValidationResult::add_warning(const ValidationError& warning) {
    warnings_.push_back(warning);
}

void ValidationResult::add_info(const ValidationError& info) {
    infos_.push_back(info);
}

std::string ValidationResult::get_summary() const {
    std::ostringstream oss;
    oss << "验证结果: ";
    oss << errors_.size() << " 个错误";
    oss << ", " << warnings_.size() << " 个警告";
    oss << ", " << infos_.size() << " 个提示";
    
    if (errors_.empty()) {
        oss << " [通过]";
    } else {
        oss << " [失败]";
    }
    
    return oss.str();
}

void ValidationResult::merge(const ValidationResult& other) {
    errors_.insert(errors_.end(), other.errors_.begin(), other.errors_.end());
    warnings_.insert(warnings_.end(), other.warnings_.begin(), other.warnings_.end());
    infos_.insert(infos_.end(), other.infos_.begin(), other.infos_.end());
}

ValidationReport::ValidationReport() {
    validation_time_ = std::chrono::system_clock::now();
}

void ValidationReport::add_result(const std::string& data_type, const ValidationResult& result) {
    results_[data_type] = result;
}

void ValidationReport::set_validation_time(const std::chrono::system_clock::time_point& time) {
    validation_time_ = time;
}

std::chrono::system_clock::time_point ValidationReport::get_validation_time() const {
    return validation_time_;
}

std::string ValidationReport::to_string() const {
    std::ostringstream oss;
    
    auto time_t_now = std::chrono::system_clock::to_time_t(validation_time_);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        validation_time_.time_since_epoch()) % 1000;
    
    oss << "=== 验证报告 ===\n";
    oss << "验证时间: " << std::put_time(std::localtime(&time_t_now), "%Y-%m-%d %H:%M:%S");
    oss << "." << std::setfill('0') << std::setw(3) << ms.count() << "\n\n";
    
    size_t total_errors = 0, total_warnings = 0;
    
    for (const auto& [data_type, result] : results_) {
        total_errors += result.get_error_count();
        total_warnings += result.get_warning_count();
        
        oss << "[" << data_type << "]\n";
        oss << "  状态: " << (result.is_valid() ? "通过" : "失败") << "\n";
        oss << "  错误: " << result.get_error_count() << "\n";
        oss << "  警告: " << result.get_warning_count() << "\n";
        oss << "  提示: " << result.get_info_count() << "\n\n";
    }
    
    oss << "=== 汇总 ===\n";
    oss << "总错误: " << total_errors << "\n";
    oss << "总警告: " << total_warnings << "\n";
    
    return oss.str();
}

bool ValidationReport::save_to_file(const std::string& file_path) const {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return false;
    }
    
    file << to_string();
    file.close();
    
    return true;
}

} // namespace tool
