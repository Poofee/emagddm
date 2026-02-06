/**
 * @file material.cpp
 * @brief Material类实现文件
 * @details 实现Material类的序列化、反序列化及属性管理功能
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
                        double temp = temp_json.value("temperature", 0.0);
                        std::vector<double> values = temp_json.value("values", std::vector<double>());
                        property.temp_dependent_data.push_back({temp, values});
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
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Material::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 序列化为JSON，然后转换为二进制
        nlohmann::json json = toJson();
        std::string json_str = json.dump();
        
        // 写入版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(version));
        
        // 写入JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Material::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        std::string json_str(data.begin() + offset, data.end());
        nlohmann::json json = nlohmann::json::parse(json_str);
        
        bool result = fromJson(json);
        if (result) {
            offset += json_str.length();
        }
        
        return result;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Material::validate() const {
    // 验证材料名称
    if (name_.empty()) return false;
    
    // 验证材料类型
    if (type_ < MatType::LINEAR_ISOTROPIC || type_ > MatType::SUPERCONDUCTOR) return false;
    
    // 验证基本属性
    if (relative_permeability_ <= 0) return false;
    if (conductivity_ < 0) return false;
    if (mass_density_ < 0) return false;
    
    // 验证磁芯损耗系数
    if (core_loss_enabled_) {
        if (core_loss_ks_ < 0 || core_loss_alpha_ < 0 || 
            core_loss_beta_ < 0 || core_loss_kn_ < 0) {
            return false;
        }
    }
    
    // 验证B-H曲线
    if (bh_curve_type_ != BHCurveType::NONE && bh_curve_.empty()) {
        return false;
    }
    
    return true;
}

void Material::setRelativePermeability(double mu_r) {
    relative_permeability_ = mu_r;
}

void Material::setConductivity(double sigma) {
    conductivity_ = sigma;
}

void Material::setMassDensity(double rho) {
    mass_density_ = rho;
}

void Material::setBHCurve(const std::vector<BHDataPoint>& curve) {
    bh_curve_ = curve;
}

void Material::setCoreLossCoefficients(double ks, double alpha, double beta, double kn) {
    core_loss_ks_ = ks;
    core_loss_alpha_ = alpha;
    core_loss_beta_ = beta;
    core_loss_kn_ = kn;
}

void Material::addProperty(const MaterialProperty& property) {
    properties_.push_back(property);
}

std::optional<MaterialProperty> Material::getProperty(const std::string& name) const {
    for (const auto& prop : properties_) {
        if (prop.name == name) {
            return prop;
        }
    }
    return std::nullopt;
}

void Material::setAnisotropicPermeability(const std::vector<double>& permeability) {
    anisotropic_permeability_ = permeability;
}

void Material::setAnisotropicConductivity(const std::vector<double>& conductivity) {
    anisotropic_conductivity_ = conductivity;
}

void Material::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

} // namespace tool