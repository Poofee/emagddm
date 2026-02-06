/**
 * @file maxwell_converter_impl.cpp
 * @brief Maxwell数据转换器实现
 * @details 实现从Maxwell解析数据到内部数据模型的直接转换功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */
#include "logger_factory.hpp"
#include "maxwell_converter_impl.hpp"
#include <sstream>
#include <regex>
#include <algorithm>
#include <cctype>

namespace tool {

// ============================================================================
// 公共转换方法实现
// ============================================================================

std::shared_ptr<Material> MaxwellConverterImpl::convertMaterialDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    if (!material_block) {
        FEEM_ERROR("MaxwellConverterImpl: 材料块节点为空");
        return nullptr;
    }

    FEEM_DEBUG("MaxwellConverterImpl: 开始转换材料数据");
    
    // 创建Material对象
    auto material = std::make_shared<Material>();
    
    // 转换材料基本属性
    convertMaterialBasicProperties(*material, material_block);
    
    // 转换电磁属性
    convertMaterialElectromagneticProperties(*material, material_block);
    
    // 转换B-H曲线
    convertMaterialBHCurve(*material, material_block);
    
    // 转换磁芯损耗参数
    convertMaterialCoreLossParameters(*material, material_block);
    
    // 转换各向异性参数
    convertMaterialAnisotropicProperties(*material, material_block);
    
    // 转换温度相关参数
    convertMaterialTemperatureProperties(*material, material_block);
    
    FEEM_DEBUG("MaxwellConverterImpl: 材料数据转换完成");
    return material;
}

std::shared_ptr<Boundary> MaxwellConverterImpl::convertBoundaryDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block) {
    
    if (!boundary_block) {
        FEEM_ERROR("MaxwellConverterImpl: 边界条件块节点为空");
        return nullptr;
    }

    FEEM_DEBUG("MaxwellConverterImpl: 开始转换边界条件数据");
    
    // 创建Boundary对象
    auto boundary = std::make_shared<Boundary>();
    
    // 转换边界条件基本属性
    convertBoundaryBasicProperties(*boundary, boundary_block);
    
    // 转换边界条件参数
    convertBoundaryParameters(*boundary, boundary_block);
    
    // 转换关联几何对象
    convertBoundaryGeometryLinks(*boundary, boundary_block);
    
    // 转换Maxwell专属边界数据
    convertBoundaryMaxwellSpecificProperties(*boundary, boundary_block);
    
    FEEM_DEBUG("MaxwellConverterImpl: 边界条件数据转换完成");
    return boundary;
}

std::shared_ptr<Excitation> MaxwellConverterImpl::convertExcitationDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    if (!excitation_block) {
        FEEM_ERROR("MaxwellConverterImpl: 激励源块节点为空");
        return nullptr;
    }

    FEEM_DEBUG("MaxwellConverterImpl: 开始转换激励源数据");
    
    // 创建Excitation对象
    auto excitation = std::make_shared<Excitation>();
    
    // 转换激励源基本属性
    convertExcitationBasicProperties(*excitation, excitation_block);
    
    // 转换激励源参数
    convertExcitationParameters(*excitation, excitation_block);
    
    // 转换线圈参数
    convertExcitationCoilParameters(*excitation, excitation_block);
    
    // 转换波形参数
    convertExcitationWaveformParameters(*excitation, excitation_block);
    
    // 转换运动参数
    convertExcitationMotionParameters(*excitation, excitation_block);
    
    // 转换Maxwell专属激励数据
    convertExcitationMaxwellSpecificProperties(*excitation, excitation_block);
    
    FEEM_DEBUG("MaxwellConverterImpl: 激励源数据转换完成");
    return excitation;
}

std::shared_ptr<Geometry> MaxwellConverterImpl::convertGeometryDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& geometry_block) {
    
    FEEM_WARN("MaxwellConverterImpl: 几何数据转换功能待实现");
    return nullptr;
}

std::shared_ptr<SolutionSetup> MaxwellConverterImpl::convertSolutionSetupDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& solution_block) {
    
    FEEM_WARN("MaxwellConverterImpl: 求解设置转换功能待实现");
    return nullptr;
}

std::shared_ptr<ProjectData> MaxwellConverterImpl::convertProjectDataDirect(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& project_block) {
    
    FEEM_WARN("MaxwellConverterImpl: 项目数据转换功能待实现");
    return nullptr;
}

// ============================================================================
// 私有转换方法实现
// ============================================================================

void MaxwellConverterImpl::convertMaterialBasicProperties(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 获取材料名称
    std::string name = findPropertyValue(material_block, "Name");
    if (!name.empty()) {
        // 移除可能的引号
        if (name.front() == '\'' && name.back() == '\'') {
            name = name.substr(1, name.length() - 2);
        }
        material.setName(name);
    }
    
    // 获取材料类型
    std::string type_str = findPropertyValue(material_block, "Type");
    if (!type_str.empty()) {
        MatType mat_type = convertMaterialType(type_str);
        material.setType(mat_type);
    }
    
    // 设置Maxwell材料ID
    material.setMaxwellMaterialID(material_block->name);
}

void MaxwellConverterImpl::convertMaterialElectromagneticProperties(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 转换相对磁导率
    std::string mu_r_str = findPropertyValue(material_block, "RelativePermeability");
    if (!mu_r_str.empty()) {
        double mu_r = parseNumericValue(mu_r_str);
        material.setRelativePermeability(mu_r);
    }
    
    // 转换电导率
    std::string sigma_str = findPropertyValue(material_block, "Conductivity");
    if (!sigma_str.empty()) {
        double sigma = parseNumericValue(sigma_str);
        material.setConductivity(sigma);
    }
    
    // 转换质量密度
    std::string density_str = findPropertyValue(material_block, "MassDensity");
    if (!density_str.empty()) {
        double density = parseNumericValue(density_str);
        material.setMassDensity(density);
    }
}

void MaxwellConverterImpl::convertMaterialBHCurve(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 检查是否有B-H曲线数据
    if (hasProperty(material_block, "BHCurve")) {
        std::string bh_data_str = findPropertyValue(material_block, "BHCurve");
        if (!bh_data_str.empty()) {
            // 解析B-H曲线数据
            std::vector<double> bh_data = parseNumericArray(bh_data_str);
            if (bh_data.size() % 2 == 0) {
                std::vector<BHDataPoint> bh_curve;
                for (size_t i = 0; i < bh_data.size(); i += 2) {
                    BHDataPoint point;
                    point.h = bh_data[i];     // H值
                    point.b = bh_data[i + 1]; // B值
                    bh_curve.push_back(point);
                }
                material.setBHCurve(bh_curve);
                
                // 设置B-H曲线类型
                material.setBHCurveType(BHCurveType::SINGLE_CURVE);
            }
        }
    }
    
    // 检查是否有自定义B-H曲线文件
    std::string bh_file_str = findPropertyValue(material_block, "BHCustomCurveFile");
    if (!bh_file_str.empty()) {
        material.setBHCustomCurveFile(bh_file_str);
    }
}

void MaxwellConverterImpl::convertMaterialCoreLossParameters(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 检查是否启用磁芯损耗
    std::string core_loss_enabled_str = findPropertyValue(material_block, "CoreLossEnabled");
    if (!core_loss_enabled_str.empty()) {
        bool enabled = parseBooleanValue(core_loss_enabled_str);
        material.setCoreLossEnabled(enabled);
    }
    
    // 转换磁芯损耗模型类型
    std::string model_str = findPropertyValue(material_block, "CoreLossModel");
    if (!model_str.empty()) {
        CoreLossModelType model = convertCoreLossModelType(model_str);
        material.setCoreLossModel(model);
    }
    
    // 转换磁芯损耗系数
    if (hasProperty(material_block, "CoreLossCoefficients")) {
        std::string coeffs_str = findPropertyValue(material_block, "CoreLossCoefficients");
        if (!coeffs_str.empty()) {
            std::vector<double> coeffs = parseNumericArray(coeffs_str);
            if (coeffs.size() >= 4) {
                material.setCoreLossCoefficients(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
            }
        }
    }
    
    // 检查是否有自定义损耗数据文件
    std::string loss_file_str = findPropertyValue(material_block, "CoreLossUserDataFile");
    if (!loss_file_str.empty()) {
        material.setCoreLossUserDataFile(loss_file_str);
    }
}

void MaxwellConverterImpl::convertMaterialAnisotropicProperties(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 转换各向异性磁导率
    if (hasProperty(material_block, "AnisotropicPermeability")) {
        std::string perm_str = findPropertyValue(material_block, "AnisotropicPermeability");
        if (!perm_str.empty()) {
            std::vector<double> perm_data = parseNumericArray(perm_str);
            material.setAnisotropicPermeability(perm_data);
        }
    }
    
    // 转换各向异性电导率
    if (hasProperty(material_block, "AnisotropicConductivity")) {
        std::string cond_str = findPropertyValue(material_block, "AnisotropicConductivity");
        if (!cond_str.empty()) {
            std::vector<double> cond_data = parseNumericArray(cond_str);
            material.setAnisotropicConductivity(cond_data);
        }
    }
}

void MaxwellConverterImpl::convertMaterialTemperatureProperties(
    Material& material,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) {
    
    // 转换温度系数
    std::string temp_coeff_str = findPropertyValue(material_block, "TemperatureCoefficient");
    if (!temp_coeff_str.empty()) {
        double temp_coeff = parseNumericValue(temp_coeff_str);
        material.setTemperatureCoefficient(temp_coeff);
    }
}

// ============================================================================
// 辅助方法实现
// ============================================================================

double MaxwellConverterImpl::parseNumericValue(const std::string& value_str) {
    try {
        return std::stod(value_str);
    } catch (const std::exception& e) {
        FEEM_ERROR("MaxwellConverterImpl: 数值转换失败: {}", value_str);
        return 0.0;
    }
}

bool MaxwellConverterImpl::parseBooleanValue(const std::string& value_str) {
    std::string lower_str = value_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str == "true" || lower_str == "1" || lower_str == "yes") {
        return true;
    } else if (lower_str == "false" || lower_str == "0" || lower_str == "no") {
        return false;
    }
    
    FEEM_WARN("MaxwellConverterImpl: 无法识别的布尔值: {}", value_str);
    return false;
}

std::vector<std::string> MaxwellConverterImpl::parseStringArray(const std::string& value_str) {
    std::vector<std::string> result;
    std::regex pattern("\\s*(?:'([^']*)'|(\\S+))\\s*");
    std::sregex_iterator it(value_str.begin(), value_str.end(), pattern);
    std::sregex_iterator end;
    
    while (it != end) {
        std::smatch match = *it;
        if (match[1].matched) {
            result.push_back(match[1].str());
        } else if (match[2].matched) {
            result.push_back(match[2].str());
        }
        ++it;
    }
    
    return result;
}

std::vector<double> MaxwellConverterImpl::parseNumericArray(const std::string& value_str) {
    std::vector<double> result;
    std::regex pattern("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?");
    std::sregex_iterator it(value_str.begin(), value_str.end(), pattern);
    std::sregex_iterator end;
    
    while (it != end) {
        std::smatch match = *it;
        try {
            result.push_back(std::stod(match.str()));
        } catch (const std::exception& e) {
            FEEM_ERROR("MaxwellConverterImpl: 数组元素转换失败: {}", match.str());
        }
        ++it;
    }
    
    return result;
}

std::string MaxwellConverterImpl::findPropertyValue(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& property_name) {
    
    auto prop = block->find_property(property_name);
    if (prop.has_value()) {
        // 将Property的value转换为字符串表示
        auto& value = prop->value;
        
        // 使用std::visit来处理variant类型
        return std::visit([](auto&& arg) -> std::string {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, std::string>) {
                return arg;
            } else if constexpr (std::is_same_v<T, double>) {
                return std::to_string(arg);
            } else if constexpr (std::is_same_v<T, bool>) {
                return arg ? "true" : "false";
            } else if constexpr (std::is_same_v<T, std::vector<fe_em::tool::maxwell_parser::Value>>) {
                // 数组类型需要特殊处理
                std::string result;
                for (const auto& elem : arg) {
                    result += std::visit([](auto&& inner_arg) -> std::string {
                        using InnerT = std::decay_t<decltype(inner_arg)>;
                        if constexpr (std::is_same_v<InnerT, std::string>) {
                            return inner_arg + " ";
                        } else if constexpr (std::is_same_v<InnerT, double>) {
                            return std::to_string(inner_arg) + " ";
                        } else if constexpr (std::is_same_v<InnerT, bool>) {
                            return (inner_arg ? "true" : "false") + " ";
                        } else {
                            return "";
                        }
                    }, elem);
                }
                return result;
            } else {
                return "";
            }
        }, value);
    }
    
    return "";
}

bool MaxwellConverterImpl::hasProperty(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& property_name) {
    
    return block->find_property(property_name).has_value();
}

MatType MaxwellConverterImpl::convertMaterialType(const std::string& maxwell_type_str) {
    std::string lower_str = maxwell_type_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("linear") != std::string::npos) {
        if (lower_str.find("anisotropic") != std::string::npos) {
            return MatType::LINEAR_ANISOTROPIC;
        } else {
            return MatType::LINEAR_ISOTROPIC;
        }
    } else if (lower_str.find("nonlinear") != std::string::npos) {
        if (lower_str.find("anisotropic") != std::string::npos) {
            return MatType::NONLINEAR_ANISOTROPIC;
        } else {
            return MatType::NONLINEAR_ISOTROPIC;
        }
    }
    
    // 默认返回线性各向同性
    return MatType::LINEAR_ISOTROPIC;
}

BHCurveType MaxwellConverterImpl::convertBHCurveType(const std::string& maxwell_curve_type_str) {
    std::string lower_str = maxwell_curve_type_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("linear") != std::string::npos || lower_str.find("single") != std::string::npos) {
        return BHCurveType::SINGLE_CURVE;
    } else if (lower_str.find("spline") != std::string::npos || lower_str.find("temp") != std::string::npos) {
        return BHCurveType::TEMP_DEPENDENT;
    } else if (lower_str.find("freq") != std::string::npos) {
        return BHCurveType::FREQ_DEPENDENT;
    } else if (lower_str.find("custom") != std::string::npos) {
        return BHCurveType::CUSTOM_CURVE;
    }
    
    return BHCurveType::NONE;
}

CoreLossModelType MaxwellConverterImpl::convertCoreLossModelType(const std::string& maxwell_model_str) {
    std::string lower_str = maxwell_model_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("steinmetz") != std::string::npos) {
        return CoreLossModelType::STEINMETZ;
    } else if (lower_str.find("bertotti") != std::string::npos) {
        return CoreLossModelType::Bertotti;
    } else if (lower_str.find("custom") != std::string::npos) {
        return CoreLossModelType::CUSTOM;
    }
    
    return CoreLossModelType::NONE;
}

// ============================================================================
// 边界条件转换方法实现
// ============================================================================

void MaxwellConverterImpl::convertBoundaryBasicProperties(
    Boundary& boundary,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block) {
    
    // 获取边界条件名称
    std::string name = findPropertyValue(boundary_block, "Name");
    if (!name.empty()) {
        // 移除可能的引号
        if (name.front() == '\'' && name.back() == '\'') {
            name = name.substr(1, name.length() - 2);
        }
        boundary.setName(name);
    }
    
    // 获取边界条件类型
    std::string type_str = findPropertyValue(boundary_block, "Type");
    if (!type_str.empty()) {
        BndType bnd_type = convertBoundaryType(type_str);
        boundary.setType(bnd_type);
    }
    
    // 设置Maxwell边界条件ID
    boundary.setMaxwellBoundaryID(boundary_block->name);
}

void MaxwellConverterImpl::convertBoundaryParameters(
    Boundary& boundary,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block) {
    
    // 转换阻抗值
    std::string impedance_str = findPropertyValue(boundary_block, "Impedance");
    if (!impedance_str.empty()) {
        double impedance = parseNumericValue(impedance_str);
        boundary.setImpedanceValue(impedance);
    }
    
    // 转换矢量势值
    std::string vector_potential_str = findPropertyValue(boundary_block, "VectorPotential");
    if (!vector_potential_str.empty()) {
        double vector_potential = parseNumericValue(vector_potential_str);
        boundary.setVectorPotential(vector_potential);
    }
    
    // 转换电压值
    std::string voltage_str = findPropertyValue(boundary_block, "Voltage");
    if (!voltage_str.empty()) {
        double voltage = parseNumericValue(voltage_str);
        boundary.setVoltage(voltage);
    }
    
    // 转换电流值
    std::string current_str = findPropertyValue(boundary_block, "Current");
    if (!current_str.empty()) {
        double current = parseNumericValue(current_str);
        boundary.setCurrent(current);
    }
}

void MaxwellConverterImpl::convertBoundaryGeometryLinks(
    Boundary& boundary,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block) {
    
    // 转换关联面
    std::string faces_str = findPropertyValue(boundary_block, "Faces");
    if (!faces_str.empty()) {
        std::vector<std::string> faces = parseStringArray(faces_str);
        for (const auto& face : faces) {
            boundary.addFace(face);
        }
    }
    
    // 转换关联边
    std::string edges_str = findPropertyValue(boundary_block, "Edges");
    if (!edges_str.empty()) {
        std::vector<std::string> edges = parseStringArray(edges_str);
        for (const auto& edge : edges) {
            boundary.addEdge(edge);
        }
    }
    
    // 转换关联对象
    std::string objects_str = findPropertyValue(boundary_block, "Objects");
    if (!objects_str.empty()) {
        std::vector<std::string> objects = parseStringArray(objects_str);
        for (const auto& object : objects) {
            boundary.addObject(object);
        }
    }
    
    // 转换主从边界名称
    std::string master_name = findPropertyValue(boundary_block, "MasterName");
    if (!master_name.empty()) {
        boundary.setMasterName(master_name);
    }
    
    std::string slave_name = findPropertyValue(boundary_block, "SlaveName");
    if (!slave_name.empty()) {
        boundary.setSlaveName(slave_name);
    }
}

void MaxwellConverterImpl::convertBoundaryMaxwellSpecificProperties(
    Boundary& boundary,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& boundary_block) {
    
    // 转换边界条件子类型
    std::string subtype_str = findPropertyValue(boundary_block, "SubType");
    if (!subtype_str.empty()) {
        BoundarySubType subtype = convertBoundarySubType(subtype_str);
        boundary.setBoundarySubType(subtype);
    }
    
    // 转换周期性映射类型
    std::string mapping_str = findPropertyValue(boundary_block, "PeriodicMappingType");
    if (!mapping_str.empty()) {
        PeriodicMappingType mapping_type = convertPeriodicMappingType(mapping_str);
        boundary.setPeriodicMappingType(mapping_type);
    }
    
    // 转换辐射距离
    std::string radiation_str = findPropertyValue(boundary_block, "RadiationDistance");
    if (!radiation_str.empty()) {
        double radiation_distance = parseNumericValue(radiation_str);
        boundary.setRadiationDistance(radiation_distance);
    }
    
    // 转换完美E对称性
    std::string e_symmetry_str = findPropertyValue(boundary_block, "PerfectESymmetry");
    if (!e_symmetry_str.empty()) {
        bool e_symmetry = parseBooleanValue(e_symmetry_str);
        boundary.setPerfectESymmetry(e_symmetry);
    }
    
    // 转换完美H对称性
    std::string h_symmetry_str = findPropertyValue(boundary_block, "PerfectHSymmetry");
    if (!h_symmetry_str.empty()) {
        bool h_symmetry = parseBooleanValue(h_symmetry_str);
        boundary.setPerfectHSymmetry(h_symmetry);
    }
    
    // 转换无限球半径
    std::string sphere_radius_str = findPropertyValue(boundary_block, "InfiniteSphereRadius");
    if (!sphere_radius_str.empty()) {
        double sphere_radius = parseNumericValue(sphere_radius_str);
        boundary.setInfiniteSphereRadius(sphere_radius);
    }
    
    // 转换边界细分参数
    std::string subdivision_str = findPropertyValue(boundary_block, "BoundarySubdivisionParameters");
    if (!subdivision_str.empty()) {
        std::vector<double> subdivision_params = parseNumericArray(subdivision_str);
        boundary.setBoundarySubdivisionParameters(subdivision_params);
    }
}

// ============================================================================
// 激励源转换方法实现
// ============================================================================

void MaxwellConverterImpl::convertExcitationBasicProperties(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 获取激励源名称
    std::string name = findPropertyValue(excitation_block, "Name");
    if (!name.empty()) {
        // 移除可能的引号
        if (name.front() == '\'' && name.back() == '\'') {
            name = name.substr(1, name.length() - 2);
        }
        excitation.setName(name);
    }
    
    // 获取激励源类型
    std::string type_str = findPropertyValue(excitation_block, "Type");
    if (!type_str.empty()) {
        ExcitationType excitation_type = convertExcitationType(type_str);
        excitation.setType(excitation_type);
    }
    
    // 设置Maxwell激励源ID
    excitation.setMaxwellExcitationID(excitation_block->name);
}

void MaxwellConverterImpl::convertExcitationParameters(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 转换激励值
    std::string value_str = findPropertyValue(excitation_block, "Value");
    if (!value_str.empty()) {
        double value = parseNumericValue(value_str);
        excitation.setValue(value);
    }
    
    // 转换相位
    std::string phase_str = findPropertyValue(excitation_block, "Phase");
    if (!phase_str.empty()) {
        double phase = parseNumericValue(phase_str);
        excitation.setPhase(phase);
    }
    
    // 转换是否为实心导体
    std::string is_solid_str = findPropertyValue(excitation_block, "IsSolid");
    if (!is_solid_str.empty()) {
        bool is_solid = parseBooleanValue(is_solid_str);
        excitation.setIsSolid(is_solid);
    }
    
    // 转换方向
    std::string direction_str = findPropertyValue(excitation_block, "Direction");
    if (!direction_str.empty()) {
        int direction = static_cast<int>(parseNumericValue(direction_str));
        excitation.setDirection(direction);
    }
}

void MaxwellConverterImpl::convertExcitationCoilParameters(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 转换线圈组
    std::string coil_group = findPropertyValue(excitation_block, "CoilGroup");
    if (!coil_group.empty()) {
        excitation.setCoilGroup(coil_group);
    }
    
    // 转换连接类型
    std::string conn_type_str = findPropertyValue(excitation_block, "ConnectionType");
    if (!conn_type_str.empty()) {
        CoilConnectionType conn_type = convertCoilConnectionType(conn_type_str);
        excitation.setConnectionType(conn_type);
    }
    
    // 转换匝数
    std::string turns_str = findPropertyValue(excitation_block, "NumberOfTurns");
    if (!turns_str.empty()) {
        int turns = static_cast<int>(parseNumericValue(turns_str));
        excitation.setNumberOfTurns(turns);
    }
    
    // 转换多边形点
    std::string polygon_str = findPropertyValue(excitation_block, "PolygonPoints");
    if (!polygon_str.empty()) {
        std::vector<double> points_data = parseNumericArray(polygon_str);
        if (points_data.size() % 2 == 0) {
            std::vector<std::pair<double, double>> polygon_points;
            for (size_t i = 0; i < points_data.size(); i += 2) {
                polygon_points.emplace_back(points_data[i], points_data[i + 1]);
            }
            excitation.setPolygonPoints(polygon_points);
        }
    }
}

void MaxwellConverterImpl::convertExcitationWaveformParameters(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 转换波形类型
    std::string waveform_str = findPropertyValue(excitation_block, "WaveformType");
    if (!waveform_str.empty()) {
        ExcitationWaveformType waveform = convertExcitationWaveformType(waveform_str);
        excitation.setWaveformType(waveform);
    }
    
    // 转换频率
    std::string frequency_str = findPropertyValue(excitation_block, "Frequency");
    if (!frequency_str.empty()) {
        double frequency = parseNumericValue(frequency_str);
        excitation.setFrequency(frequency);
    }
    
    // 转换占空比
    std::string duty_str = findPropertyValue(excitation_block, "DutyCycle");
    if (!duty_str.empty()) {
        double duty_cycle = parseNumericValue(duty_str);
        excitation.setDutyCycle(duty_cycle);
    }
    
    // 转换波形参数
    std::string waveform_params_str = findPropertyValue(excitation_block, "WaveformParameters");
    if (!waveform_params_str.empty()) {
        std::vector<double> waveform_params = parseNumericArray(waveform_params_str);
        excitation.setWaveformParameters(waveform_params);
    }
    
    // 转换外部电路文件
    std::string circuit_file = findPropertyValue(excitation_block, "ExternalCircuitFile");
    if (!circuit_file.empty()) {
        excitation.setExternalCircuitFile(circuit_file);
    }
    
    // 转换自定义波形文件
    std::string custom_waveform_file = findPropertyValue(excitation_block, "CustomWaveformFile");
    if (!custom_waveform_file.empty()) {
        excitation.setCustomWaveformFile(custom_waveform_file);
    }
}

void MaxwellConverterImpl::convertExcitationMotionParameters(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 转换绕组类型
    std::string winding_str = findPropertyValue(excitation_block, "WindingType");
    if (!winding_str.empty()) {
        WindingType winding = convertWindingType(winding_str);
        excitation.setWindingType(winding);
    }
    
    // 转换运动类型
    std::string motion_str = findPropertyValue(excitation_block, "MotionType");
    if (!motion_str.empty()) {
        MotionType motion = convertMotionType(motion_str);
        excitation.setMotionType(motion);
    }
    
    // 转换旋转速度
    std::string rotation_str = findPropertyValue(excitation_block, "RotationSpeed");
    if (!rotation_str.empty()) {
        double rotation_speed = parseNumericValue(rotation_str);
        excitation.setRotationSpeed(rotation_speed);
    }
    
    // 转换平移速度
    std::string translation_str = findPropertyValue(excitation_block, "TranslationSpeed");
    if (!translation_str.empty()) {
        double translation_speed = parseNumericValue(translation_str);
        excitation.setTranslationSpeed(translation_speed);
    }
}

void MaxwellConverterImpl::convertExcitationMaxwellSpecificProperties(
    Excitation& excitation,
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& excitation_block) {
    
    // 这里可以添加Maxwell专属激励参数的转换逻辑
    // 例如特定的参数映射或特殊处理
}

// ============================================================================
// 枚举类型转换方法实现
// ============================================================================

BndType MaxwellConverterImpl::convertBoundaryType(const std::string& maxwell_type_str) {
    std::string lower_str = maxwell_type_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("dirichlet") != std::string::npos) {
        return BndType::DIRICHLET;
    } else if (lower_str.find("neumann") != std::string::npos) {
        return BndType::NEUMANN;
    } else if (lower_str.find("periodic") != std::string::npos) {
        return BndType::PERIODIC;
    } else if (lower_str.find("even") != std::string::npos && lower_str.find("symmetry") != std::string::npos) {
        return BndType::EVEN_SYMMETRY;
    } else if (lower_str.find("odd") != std::string::npos && lower_str.find("symmetry") != std::string::npos) {
        return BndType::ODD_SYMMETRY;
    } else if (lower_str.find("radiation") != std::string::npos) {
        return BndType::RADIATION;
    } else if (lower_str.find("impedance") != std::string::npos) {
        return BndType::IMPEDANCE;
    }
    
    // 默认返回Dirichlet边界
    return BndType::DIRICHLET;
}

BoundarySubType MaxwellConverterImpl::convertBoundarySubType(const std::string& maxwell_subtype_str) {
    std::string lower_str = maxwell_subtype_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("skin") != std::string::npos) {
        return BoundarySubType::SKIN_DEPTH;
    } else if (lower_str.find("eddy") != std::string::npos) {
        return BoundarySubType::EDDY_CURRENT;
    } else if (lower_str.find("proximity") != std::string::npos) {
        return BoundarySubType::PROXIMITY_EFFECT;
    } else if (lower_str.find("edge") != std::string::npos) {
        return BoundarySubType::EDGE_BASED;
    } else if (lower_str.find("face") != std::string::npos) {
        return BoundarySubType::FACE_BASED;
    }
    
    return BoundarySubType::NONE;
}

PeriodicMappingType MaxwellConverterImpl::convertPeriodicMappingType(const std::string& maxwell_mapping_str) {
    std::string lower_str = maxwell_mapping_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("translational") != std::string::npos) {
        return PeriodicMappingType::TRANSLATIONAL;
    } else if (lower_str.find("rotational") != std::string::npos) {
        return PeriodicMappingType::ROTATIONAL;
    }
    
    return PeriodicMappingType::NONE;
}

ExcitationType MaxwellConverterImpl::convertExcitationType(const std::string& maxwell_type_str) {
    std::string lower_str = maxwell_type_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("current") != std::string::npos) {
        return ExcitationType::CURRENT_DENSITY;
    } else if (lower_str.find("voltage") != std::string::npos) {
        return ExcitationType::VOLTAGE_SOURCE;
    } else if (lower_str.find("coil") != std::string::npos) {
        return ExcitationType::COIL;
    } else if (lower_str.find("winding") != std::string::npos) {
        return ExcitationType::WINDING;
    }
    
    // 默认返回电流密度激励
    return ExcitationType::CURRENT_DENSITY;
}

ExcitationWaveformType MaxwellConverterImpl::convertExcitationWaveformType(const std::string& maxwell_waveform_str) {
    std::string lower_str = maxwell_waveform_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("dc") != std::string::npos) {
        return ExcitationWaveformType::DC;
    } else if (lower_str.find("sinusoidal") != std::string::npos) {
        return ExcitationWaveformType::SINUSOIDAL;
    } else if (lower_str.find("cosine") != std::string::npos) {
        return ExcitationWaveformType::COSINE;
    } else if (lower_str.find("square") != std::string::npos) {
        return ExcitationWaveformType::SQUARE;
    } else if (lower_str.find("pulse") != std::string::npos) {
        return ExcitationWaveformType::PULSE;
    }
    
    return ExcitationWaveformType::DC;
}

CoilConnectionType MaxwellConverterImpl::convertCoilConnectionType(const std::string& maxwell_conn_str) {
    std::string lower_str = maxwell_conn_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("series") != std::string::npos) {
        return CoilConnectionType::SERIES;
    } else if (lower_str.find("parallel") != std::string::npos) {
        return CoilConnectionType::PARALLEL;
    }
    
    return CoilConnectionType::SERIES;
}

WindingType MaxwellConverterImpl::convertWindingType(const std::string& maxwell_winding_str) {
    std::string lower_str = maxwell_winding_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("solid") != std::string::npos) {
        return WindingType::SOLID;
    } else if (lower_str.find("stranded") != std::string::npos) {
        return WindingType::STRANDED;
    }
    
    return WindingType::SOLID;
}

MotionType MaxwellConverterImpl::convertMotionType(const std::string& maxwell_motion_str) {
    std::string lower_str = maxwell_motion_str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str.find("rotation") != std::string::npos) {
        return MotionType::ROTATION;
    } else if (lower_str.find("translation") != std::string::npos) {
        return MotionType::TRANSLATION;
    }
    
    return MotionType::NONE;
}

} // namespace tool