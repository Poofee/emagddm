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

    // 材料元信息
    json["library"] = library_;
    json["lib_location"] = lib_location_;
    json["mod_since_lib"] = mod_since_lib_;
    json["mod_time"] = mod_time_;
    json["physics_types"] = physics_types_;
    json["coordinate_system_type"] = coordinate_system_type_;
    json["bulk_or_surface_type"] = bulk_or_surface_type_;

    // 扩展物理属性
    json["permittivity"] = permittivity_;
    json["youngs_modulus"] = youngs_modulus_;
    json["poisons_ratio"] = poisons_ratio_;
    json["thermal_expansion_coefficient"] = thermal_expansion_coefficient_;

    // B-H曲线增强信息
    json["bh_curve_data_type"] = bh_curve_data_type_;
    json["h_unit"] = h_unit_;
    json["b_unit"] = b_unit_;
    json["is_temperature_dependent"] = is_temperature_dependent_;

    // 温度相关B-H曲线
    nlohmann::json temp_bh_json;
    for (const auto& [temp, points] : temperature_bh_curves_) {
        std::string key = std::to_string(temp);
        nlohmann::json points_arr;
        for (const auto& pt : points) {
            nlohmann::json pt_json;
            pt_json[0] = pt.h;
            pt_json[1] = pt.b;
            points_arr.push_back(pt_json);
        }
        temp_bh_json[key] = points_arr;
    }
    json["temperature_bh_curves"] = temp_bh_json;

    // 铁损增强
    json["core_loss_type_choice"] = core_loss_type_choice_;
    json["stacking_type_str"] = stacking_type_str_;

    // 多频率铁损曲线
    nlohmann::json freq_curves_json;
    for (const auto& [freq, bp_list] : core_loss_freq_curves_) {
        std::string key = std::to_string(freq);
        nlohmann::json bp_arr;
        for (const auto& [b_val, p_val] : bp_list) {
            nlohmann::json bp_json;
            bp_json[0] = b_val;
            bp_json[1] = p_val;
            bp_arr.push_back(bp_json);
        }
        freq_curves_json[key] = bp_arr;
    }
    json["core_loss_freq_curves"] = freq_curves_json;
    json["core_loss_unit"] = core_loss_unit_;

    // 铁损系数设置数据
    if (core_loss_coeff_setup_.has_value()) {
        nlohmann::json coeff_json;
        coeff_json["mode"] = core_loss_coeff_setup_->mode;
        coeff_json["frequency"] = core_loss_coeff_setup_->frequency;
        coeff_json["thickness"] = core_loss_coeff_setup_->thickness;
        coeff_json["conductivity"] = core_loss_coeff_setup_->conductivity;
        nlohmann::json bp_arr;
        for (const auto& [b_val, p_val] : core_loss_coeff_setup_->bp_curve) {
            nlohmann::json bp_json;
            bp_json[0] = b_val;
            bp_json[1] = p_val;
            bp_arr.push_back(bp_json);
        }
        coeff_json["bp_curve"] = bp_arr;
        json["core_loss_coeff_setup"] = coeff_json;
    }

    // 热修正器
    nlohmann::json tm_json;
    for (const auto& mod : thermal_modifiers_) {
        nlohmann::json m;
        m["property_name"] = mod.property_name;
        m["index"] = mod.index;
        m["formula_string"] = mod.formula_string;
        tm_json.push_back(m);
    }
    json["thermal_modifiers"] = tm_json;

    // 外观数据
    if (appearance_data_.has_value()) {
        nlohmann::json app_json;
        app_json["red"] = appearance_data_->red;
        app_json["green"] = appearance_data_->green;
        app_json["blue"] = appearance_data_->blue;
        app_json["transparency"] = appearance_data_->transparency;
        json["appearance_data"] = app_json;
    }

    // 备注和关键词
    json["notes"] = notes_;
    json["keywords"] = keywords_;
    
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

        // 材料元信息
        library_ = json.value("library", "");
        lib_location_ = json.value("lib_location", "");
        mod_since_lib_ = json.value("mod_since_lib", false);
        mod_time_ = json.value("mod_time", 0LL);
        physics_types_ = json.value("physics_types", std::vector<std::string>());
        coordinate_system_type_ = json.value("coordinate_system_type", "Cartesian");
        bulk_or_surface_type_ = json.value("bulk_or_surface_type", 1);

        // 扩展物理属性
        permittivity_ = json.value("permittivity", 1.0);
        youngs_modulus_ = json.value("youngs_modulus", 0.0);
        poisons_ratio_ = json.value("poisons_ratio", 0.0);
        thermal_expansion_coefficient_ = json.value("thermal_expansion_coefficient", 0.0);

        // B-H曲线增强信息
        bh_curve_data_type_ = json.value("bh_curve_data_type", "normal");
        h_unit_ = json.value("h_unit", "A_per_meter");
        b_unit_ = json.value("b_unit", "tesla");
        is_temperature_dependent_ = json.value("is_temperature_dependent", false);

        // 温度相关B-H曲线
        temperature_bh_curves_.clear();
        if (json.contains("temperature_bh_curves") && json["temperature_bh_curves"].is_object()) {
            for (const auto& [key, points_arr] : json["temperature_bh_curves"].items()) {
                double temp = 0;
                try { temp = std::stod(key); } catch (...) { continue; }
                std::vector<BHDataPoint> bh_points;
                for (const auto& pt : points_arr) {
                    BHDataPoint dp;
                    dp.h = pt[0].get<double>();
                    dp.b = pt[1].get<double>();
                    bh_points.push_back(dp);
                }
                temperature_bh_curves_[temp] = bh_points;
            }
        }

        // 铁损增强
        core_loss_type_choice_ = json.value("core_loss_type_choice", "");
        stacking_type_str_ = json.value("stacking_type_str", "Solid");

        // 多频率铁损曲线
        core_loss_freq_curves_.clear();
        if (json.contains("core_loss_freq_curves") && json["core_loss_freq_curves"].is_object()) {
            for (const auto& [key, bp_arr] : json["core_loss_freq_curves"].items()) {
                double freq = 0;
                try { freq = std::stod(key); } catch (...) { continue; }
                std::vector<std::pair<double, double>> bp_list;
                for (const auto& pt : bp_arr) {
                    bp_list.push_back({pt[0].get<double>(), pt[1].get<double>()});
                }
                core_loss_freq_curves_[freq] = bp_list;
            }
        }
        core_loss_unit_ = json.value("core_loss_unit", "w_per_kg");

        // 铁损系数设置数据
        core_loss_coeff_setup_.reset();
        if (json.contains("core_loss_coeff_setup") && json["core_loss_coeff_setup"].is_object()) {
            const auto& coeff_json = json["core_loss_coeff_setup"];
            CoreLossCoefficientSetup coeff;
            coeff.mode = coeff_json.value("mode", "");
            coeff.frequency = coeff_json.value("frequency", 0.0);
            coeff.thickness = coeff_json.value("thickness", "");
            coeff.conductivity = coeff_json.value("conductivity", "");
            if (coeff_json.contains("bp_curve") && coeff_json["bp_curve"].is_array()) {
                for (const auto& pt : coeff_json["bp_curve"]) {
                    coeff.bp_curve.push_back({pt[0].get<double>(), pt[1].get<double>()});
                }
            }
            core_loss_coeff_setup_ = coeff;
        }

        // 热修正器
        thermal_modifiers_.clear();
        if (json.contains("thermal_modifiers") && json["thermal_modifiers"].is_array()) {
            for (const auto& tm : json["thermal_modifiers"]) {
                ThermalModifier mod;
                mod.property_name = tm.value("property_name", "");
                mod.index = tm.value("index", 0);
                mod.formula_string = tm.value("formula_string", "");
                thermal_modifiers_.push_back(mod);
            }
        }

        // 外观数据
        appearance_data_.reset();
        if (json.contains("appearance_data") && json["appearance_data"].is_object()) {
            const auto& app_json = json["appearance_data"];
            appearance_data_ = MaterialAppearance{
                app_json.value("red", 200),
                app_json.value("green", 200),
                app_json.value("blue", 200),
                app_json.value("transparency", 0.0)
            };
        }

        // 备注和关键词
        notes_ = json.value("notes", "");
        keywords_ = json.value("keywords", "");
        
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
        (void)*reinterpret_cast<const uint32_t*>(&data[offset]);
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