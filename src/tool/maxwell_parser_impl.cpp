// TODO: 材料数据提取功能已完成
// 实现了从Maxwell .aedt文件中提取完整材料数据的功能
// 包括：
// 1. 线性材料（如copper）：电导率、磁导率、质量密度、热导率、比热容
// 2. 非线性材料（如B27AV1400）：BH曲线、磁芯损耗系数、矫顽力
// 3. 自动识别线性/非线性磁导率
// 4. 完整的测试用例验证

#include "tool/maxwell_parser_impl.hpp"
#include "tool/logger_factory.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>

namespace tool {

bool MaxwellParserImpl::canParse(const std::string& file_path) {
    // 检查文件扩展名
    std::filesystem::path path(file_path);
    std::string extension = path.extension().string();
    
    // 支持的文件扩展名
    std::vector<std::string> supported_extensions = {
        ".aedt", ".aedtz", ".amat", ".xml"
    };
    
    bool supported = std::find(supported_extensions.begin(), 
                              supported_extensions.end(), 
                              extension) != supported_extensions.end();
    
    if (supported) {
        FEEM_DEBUG("文件格式支持: " + file_path);
        // 设置文件路径，供后续解析使用
        file_path_ = file_path;
    } else {
        FEEM_DEBUG("文件格式不支持: " + file_path);
    }
    
    return supported;
}

MaxwellFileInfo MaxwellParserImpl::parseFileInfo() {
    if (!parser_.parse_file(file_path_)) {
        FEEM_ERROR("文件解析失败: " + file_path_);
        throw MaxwellException(MaxwellErrorCode::INVALID_FORMAT, 
                              "文件解析失败: " + file_path_);
    }
    
    return extractFileInfo();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseMaterials() {
    // 如果解析树无效且有文件路径，则先解析文件
    if (!isParseTreeValid()) {
        if (!file_path_.empty()) {
            if (!parser_.parse_file(file_path_)) {
                FEEM_ERROR("文件解析失败: " + file_path_);
                throw MaxwellException(MaxwellErrorCode::INVALID_FORMAT, 
                                      "文件解析失败: " + file_path_);
            }
        } else {
            FEEM_ERROR("解析树无效，请先调用parseFileInfo");
            throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                                  "解析树无效，请先调用parseFileInfo");
        }
    }
    
    FEEM_DEBUG("开始解析材料数据");
    return extractMaterials();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseBoundaries() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析边界条件数据");
    return extractBoundaries();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseExcitations() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析激励源数据");
    return extractExcitations();
}

nlohmann::json MaxwellParserImpl::parseSolutionSetup() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析求解设置数据");
    return extractSolutionSetup();
}

nlohmann::json MaxwellParserImpl::parseGeometry() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析几何数据");
    return extractGeometry();
}

nlohmann::json MaxwellParserImpl::parseAllData() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析所有数据");
    
    nlohmann::json result;
    result["file_info"] = extractFileInfo().toJson();
    result["materials"] = extractMaterials();
    result["boundaries"] = extractBoundaries();
    result["excitations"] = extractExcitations();
    result["solution_setup"] = extractSolutionSetup();
    result["geometry"] = extractGeometry();
    
    return result;
}

// 私有方法实现

MaxwellFileInfo MaxwellParserImpl::extractFileInfo() const {
    MaxwellFileInfo info;
    
    auto root = parser_.get_root();
    if (!root) {
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, "根节点为空");
    }
    
    info.file_path = file_path_;
    info.file_format = "Maxwell .aedt";
    
    // 从根块提取基本信息
    if (auto created_prop = root->find_property("Created")) {
        info.created_date = std::get<std::string>(created_prop->value);
    }
    
    if (auto product_prop = root->find_property("Product")) {
        info.maxwell_version = std::get<std::string>(product_prop->value);
    }
    
    // 从Desktop块提取版本信息
    auto desktop_block = findBlock("Desktop", root);
    if (desktop_block) {
        // 可以进一步提取版本信息
    }
    
    // 设置文件大小
    std::filesystem::path path(file_path_);
    if (std::filesystem::exists(path)) {
        info.file_size = std::filesystem::file_size(path);
    }
    
    return info;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractMaterials() const {
    std::vector<nlohmann::json> materials;
    
    auto root = parser_.get_root();
    if (!root) return materials;
    
    // 查找Definitions块
    auto definitions_block = findBlock("Definitions", root);
    if (!definitions_block) {
        FEEM_DEBUG("未找到Definitions块");
        return materials;
    }
    
    // 从Definitions中查找Materials块
    auto materials_block = findBlock("Materials", definitions_block);
    if (!materials_block) {
        FEEM_DEBUG("未找到Materials块");
        return materials;
    }
    
    FEEM_DEBUG("找到材料块，包含 " + std::to_string(materials_block->children.size()) + " 个材料");
    
    // 提取所有材料子块
    for (const auto& material_block : materials_block->children) {
        nlohmann::json material = extractSingleMaterial(material_block);
        if (!material.empty()) {
            materials.push_back(material);
        }
    }
    
    FEEM_DEBUG("成功提取 " + std::to_string(materials.size()) + " 个材料");
    return materials;
}

/**
 * @brief 从单个材料块提取完整材料数据
 * @param material_block 材料块节点
 * @return 材料JSON数据
 */
nlohmann::json MaxwellParserImpl::extractSingleMaterial(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block) const {
    
    nlohmann::json material;
    
    if (!material_block) return material;
    
    // 材料名称
    std::string material_name = material_block->name;
    material["name"] = material_name;
    
    // 提取电导率 (conductivity)
    if (auto prop = material_block->find_property("conductivity")) {
        material["conductivity"] = propertyValueToDouble(*prop);
    }
    
    // 提取质量密度 (mass_density)
    if (auto prop = material_block->find_property("mass_density")) {
        material["mass_density"] = propertyValueToDouble(*prop);
    }
    
    // 提取热导率 (thermal_conductivity)
    if (auto prop = material_block->find_property("thermal_conductivity")) {
        material["thermal_conductivity"] = propertyValueToDouble(*prop);
    }
    
    // 提取比热容 (specific_heat)
    if (auto prop = material_block->find_property("specific_heat")) {
        material["specific_heat"] = propertyValueToDouble(*prop);
    }
    
    // 提取磁芯损耗系数
    if (auto prop = material_block->find_property("core_loss_kh")) {
        material["core_loss_kh"] = propertyValueToDouble(*prop);
    }
    if (auto prop = material_block->find_property("core_loss_kc")) {
        material["core_loss_kc"] = propertyValueToDouble(*prop);
    }
    if (auto prop = material_block->find_property("core_loss_ke")) {
        material["core_loss_ke"] = propertyValueToDouble(*prop);
    }
    
    // 提取磁导率 - 需要检查是否为非线性
    extractPermeability(material_block, material);
    
    // 提取矫顽力 (magnetic_coercivity)
    extractMagneticCoercivity(material_block, material);
    
    return material;
}

/**
 * @brief 从材料块提取磁导率数据（线性或非线性）
 * @param material_block 材料块节点
 * @param material 输出JSON对象
 */
void MaxwellParserImpl::extractPermeability(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {
    
    // 首先检查是否为简单的permeability属性（线性材料）
    if (auto prop = material_block->find_property("permeability")) {
        // 检查属性值是否为数字（线性磁导率）
        const auto& value = prop->value;
        if (std::holds_alternative<double>(value)) {
            material["permeability"] = std::get<double>(value);
            material["permeability_type"] = "linear";
            return;
        } else if (std::holds_alternative<std::string>(value)) {
            // 尝试转换为double
            std::string str_val = std::get<std::string>(value);
            try {
                double val = std::stod(str_val);
                material["permeability"] = val;
                material["permeability_type"] = "linear";
                return;
            } catch (...) {
                // 不是数字，可能是特殊标记
            }
        }
    }
    
    // 检查是否有permeability子块（非线性材料）
    auto perm_block = findBlock("permeability", material_block);
    if (perm_block) {
        material["permeability_type"] = "nonlinear";
        
        // 提取BH曲线数据
        auto bh_block = findBlock("BHCoordinates", perm_block);
        if (bh_block) {
            // 从Points属性提取BH曲线
            if (auto prop = bh_block->find_property("Points")) {
                auto bh_curve = extractBHCurve(*prop);
                if (!bh_curve.empty()) {
                    material["bh_curve"] = bh_curve;
                }
            }
        }
    }
}

/**
 * @brief 从材料块提取矫顽力数据
 * @param material_block 材料块节点
 * @param material 输出JSON对象
 */
void MaxwellParserImpl::extractMagneticCoercivity(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {
    
    // 查找magnetic_coercivity子块
    auto coercivity_block = findBlock("magnetic_coercivity", material_block);
    if (coercivity_block) {
        nlohmann::json coercivity;
        
        // 提取Magnitude
        if (auto prop = coercivity_block->find_property("Magnitude")) {
            coercivity["magnitude"] = propertyValueToString(*prop);
        }
        
        // 提取方向分量
        std::vector<double> direction;
        if (auto prop = coercivity_block->find_property("DirComp1")) {
            direction.push_back(propertyValueToDouble(*prop));
        }
        if (auto prop = coercivity_block->find_property("DirComp2")) {
            direction.push_back(propertyValueToDouble(*prop));
        }
        if (auto prop = coercivity_block->find_property("DirComp3")) {
            direction.push_back(propertyValueToDouble(*prop));
        }
        
        if (!direction.empty()) {
            coercivity["direction"] = direction;
        }
        
        if (!coercivity.empty()) {
            material["coercivity"] = coercivity;
        }
    }
}

/**
 * @brief 从属性中提取双精度数值
 * @param prop 属性节点
 * @return 双精度数值
 */
double MaxwellParserImpl::propertyValueToDouble(
    const fe_em::tool::maxwell_parser::Property& prop) const {
    
    const auto& value = prop.value;
    
    if (std::holds_alternative<double>(value)) {
        return std::get<double>(value);
    } else if (std::holds_alternative<std::string>(value)) {
        std::string str_val = std::get<std::string>(value);
        try {
            return std::stod(str_val);
        } catch (...) {
            return 0.0;
        }
    }
    
    return 0.0;
}

/**
 * @brief 从属性中提取字符串值
 * @param prop 属性节点
 * @return 字符串值
 */
std::string MaxwellParserImpl::propertyValueToString(
    const fe_em::tool::maxwell_parser::Property& prop) const {
    
    const auto& value = prop.value;
    
    if (std::holds_alternative<std::string>(value)) {
        return std::get<std::string>(value);
    }
    
    return "";
}

/**
 * @brief 从Points属性提取BH曲线数据
 * @param prop Points属性节点
 * @return BH曲线数据 [[H1,B1], [H2,B2], ...]
 */
std::vector<std::vector<double>> MaxwellParserImpl::extractBHCurve(
    const fe_em::tool::maxwell_parser::Property& prop) const {
    
    std::vector<std::vector<double>> bh_curve;
    
    const auto& value = prop.value;
    
    // 处理数组类型的值
    if (std::holds_alternative<std::vector<double>>(value)) {
        const auto& points = std::get<std::vector<double>>(value);
        
        // BH曲线是交替的 H, B, H, B, ... 序列
        for (size_t i = 0; i + 1 < points.size(); i += 2) {
            std::vector<double> point = {points[i], points[i + 1]};
            bh_curve.push_back(point);
        }
    }
    // 注意：Value类型不包含std::vector<Value>，只有std::vector<double>和std::vector<std::string>
    
    return bh_curve;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractBoundaries() const {
    std::vector<nlohmann::json> boundaries;
    
    auto root = parser_.get_root();
    if (!root) return boundaries;
    
    // 查找 Maxwell2DModel 块
    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) {
        FEEM_WARN("未找到 Maxwell2DModel 块");
        return boundaries;
    }
    
    // 查找 BoundarySetup 块
    auto boundary_setup = findBlock("BoundarySetup", maxwell_model);
    if (!boundary_setup) {
        FEEM_WARN("未找到 BoundarySetup 块");
        return boundaries;
    }
    
    // 查找 Boundaries 块
    auto boundaries_block = findBlock("Boundaries", boundary_setup);
    if (!boundaries_block) {
        FEEM_WARN("未找到 Boundaries 块");
        return boundaries;
    }
    
    // 遍历所有边界条件子块
    for (const auto& boundary_block : boundaries_block->children) {
        nlohmann::json boundary_json;
        
        // 边界条件名称
        boundary_json["name"] = boundary_block->name;
        
        // 提取 ID
        if (auto id_prop = boundary_block->find_property("ID")) {
            if (std::holds_alternative<double>(id_prop->value)) {
                boundary_json["id"] = static_cast<int>(std::get<double>(id_prop->value));
            }
        }
        
        // 提取 BoundType
        if (auto bound_type_prop = boundary_block->find_property("BoundType")) {
            if (std::holds_alternative<std::string>(bound_type_prop->value)) {
                boundary_json["bound_type"] = std::get<std::string>(bound_type_prop->value);
            }
        }
        
        // 提取 Edges（函数调用形式）
        if (auto edges_prop = boundary_block->find_property("Edges")) {
            if (std::holds_alternative<std::vector<double>>(edges_prop->value)) {
                auto edges = std::get<std::vector<double>>(edges_prop->value);
                std::vector<int> edges_int;
                for (double edge : edges) {
                    edges_int.push_back(static_cast<int>(edge));
                }
                boundary_json["edges"] = edges_int;
            } else if (std::holds_alternative<double>(edges_prop->value)) {
                boundary_json["edges"] = std::vector<int>{static_cast<int>(std::get<double>(edges_prop->value))};
            }
        }
        
        // 提取 Objects（函数调用形式）
        if (auto objects_prop = boundary_block->find_property("Objects")) {
            if (std::holds_alternative<std::vector<double>>(objects_prop->value)) {
                auto objects = std::get<std::vector<double>>(objects_prop->value);
                std::vector<int> objects_int;
                for (double obj : objects) {
                    objects_int.push_back(static_cast<int>(obj));
                }
                boundary_json["objects"] = objects_int;
            } else if (std::holds_alternative<double>(objects_prop->value)) {
                boundary_json["objects"] = std::vector<int>{static_cast<int>(std::get<double>(objects_prop->value))};
            }
        }
        
        // 提取 Value
        if (auto value_prop = boundary_block->find_property("Value")) {
            if (std::holds_alternative<std::string>(value_prop->value)) {
                boundary_json["value"] = std::get<std::string>(value_prop->value);
            }
        }
        
        // 提取 Winding
        if (auto winding_prop = boundary_block->find_property("Winding")) {
            if (std::holds_alternative<double>(winding_prop->value)) {
                boundary_json["winding"] = static_cast<int>(std::get<double>(winding_prop->value));
            }
        }
        
        // 提取 PolarityType
        if (auto polarity_prop = boundary_block->find_property("PolarityType")) {
            if (std::holds_alternative<std::string>(polarity_prop->value)) {
                boundary_json["polarity_type"] = std::get<std::string>(polarity_prop->value);
            }
        }
        
        // 提取 Conductor number
        if (auto conductor_prop = boundary_block->find_property("Conductor number")) {
            if (std::holds_alternative<std::string>(conductor_prop->value)) {
                boundary_json["conductor_number"] = std::get<std::string>(conductor_prop->value);
            }
        }
        
        boundaries.push_back(boundary_json);
    }
    
    FEEM_INFO("成功提取 " + std::to_string(boundaries.size()) + " 个边界条件");
    return boundaries;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractExcitations() const {
    std::vector<nlohmann::json> excitations;
    
    // 从边界条件中识别激励源（BoundType='Coil' 即为激励源）
    auto boundaries = extractBoundaries();
    
    for (const auto& boundary : boundaries) {
        // 检查是否为 Coil 类型
        std::string bound_type = "";
        if (boundary.contains("bound_type") && boundary["bound_type"].is_string()) {
            bound_type = boundary["bound_type"].get<std::string>();
        }
        if (bound_type != "Coil") {
            continue;
        }
        
        nlohmann::json excitation;
        
        // 激励源名称
        excitation["name"] = boundary.value("name", "");
        
        // 激励源类型
        excitation["type"] = "Coil";
        
        // 提取 Winding
        if (boundary.contains("winding")) {
            excitation["winding"] = boundary["winding"];
        }
        
        // 提取极性
        if (boundary.contains("polarity_type")) {
            excitation["polarity"] = boundary["polarity_type"];
        }
        
        // 提取 Objects
        if (boundary.contains("objects")) {
            excitation["objects"] = boundary["objects"];
        }
        
        // 提取 Conductor number
        if (boundary.contains("conductor_number")) {
            excitation["conductor_number"] = boundary["conductor_number"];
        }
        
        excitations.push_back(excitation);
    }
    
    FEEM_INFO("成功提取 " + std::to_string(excitations.size()) + " 个激励源");
    return excitations;
}

nlohmann::json MaxwellParserImpl::extractSolutionSetup() const {
    nlohmann::json setup;
    
    auto root = parser_.get_root();
    if (!root) return setup;
    
    // 查找 Maxwell2DModel 块
    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) {
        FEEM_WARN("未找到 Maxwell2DModel 块");
        return setup;
    }
    
    // 查找 AnalysisSetup 块
    auto analysis_setup = findBlock("AnalysisSetup", maxwell_model);
    if (!analysis_setup) {
        FEEM_WARN("未找到 AnalysisSetup 块");
        return setup;
    }
    
    // 查找 SolveSetups 块
    auto solve_setups = findBlock("SolveSetups", analysis_setup);
    if (!solve_setups) {
        FEEM_WARN("未找到 SolveSetups 块");
        return setup;
    }
    
    // 获取第一个求解设置（通常是 Setup1）
    if (solve_setups->children.empty()) {
        FEEM_WARN("SolveSetups 块中没有求解设置");
        return setup;
    }
    
    auto setup_block = solve_setups->children[0];
    
    // 求解设置名称
    setup["name"] = setup_block->name;
    
    // 提取 SetupType
    if (auto setup_type_prop = setup_block->find_property("SetupType")) {
        if (std::holds_alternative<std::string>(setup_type_prop->value)) {
            setup["setup_type"] = std::get<std::string>(setup_type_prop->value);
        }
    }
    
    // 提取 StopTime
    if (auto stop_time_prop = setup_block->find_property("StopTime")) {
        if (std::holds_alternative<std::string>(stop_time_prop->value)) {
            setup["stop_time"] = std::get<std::string>(stop_time_prop->value);
        }
    }
    
    // 提取 TimeStep
    if (auto time_step_prop = setup_block->find_property("TimeStep")) {
        if (std::holds_alternative<std::string>(time_step_prop->value)) {
            setup["time_step"] = std::get<std::string>(time_step_prop->value);
        }
    }
    
    // 提取 NonlinearSolverResidual
    if (auto residual_prop = setup_block->find_property("NonlinearSolverResidual")) {
        if (std::holds_alternative<std::string>(residual_prop->value)) {
            setup["nonlinear_residual"] = std::get<std::string>(residual_prop->value);
        }
    }
    
    // 提取 FrequencyOfAddedVoltageSource
    if (auto frequency_prop = setup_block->find_property("FrequencyOfAddedVoltageSource")) {
        if (std::holds_alternative<std::string>(frequency_prop->value)) {
            setup["frequency"] = std::get<std::string>(frequency_prop->value);
        }
    }
    
    // 提取 UseAdaptiveTimeStep
    if (auto adaptive_prop = setup_block->find_property("UseAdaptiveTimeStep")) {
        if (std::holds_alternative<bool>(adaptive_prop->value)) {
            setup["adaptive_time_step"] = std::get<bool>(adaptive_prop->value);
        }
    }
    
    // 提取 InitialTimeStep
    if (auto initial_step_prop = setup_block->find_property("InitialTimeStep")) {
        if (std::holds_alternative<std::string>(initial_step_prop->value)) {
            setup["initial_time_step"] = std::get<std::string>(initial_step_prop->value);
        }
    }
    
    // 提取 MinTimeStep
    if (auto min_step_prop = setup_block->find_property("MinTimeStep")) {
        if (std::holds_alternative<std::string>(min_step_prop->value)) {
            setup["min_time_step"] = std::get<std::string>(min_step_prop->value);
        }
    }
    
    // 提取 MaxTimeStep
    if (auto max_step_prop = setup_block->find_property("MaxTimeStep")) {
        if (std::holds_alternative<std::string>(max_step_prop->value)) {
            setup["max_time_step"] = std::get<std::string>(max_step_prop->value);
        }
    }
    
    // 提取 StopCriterion
    if (auto criterion_prop = setup_block->find_property("StopCriterion")) {
        if (std::holds_alternative<double>(criterion_prop->value)) {
            setup["stop_criterion"] = std::get<double>(criterion_prop->value);
        }
    }
    
    FEEM_INFO("成功提取求解设置: {}", setup["name"].get<std::string>());
    return setup;
}

nlohmann::json MaxwellParserImpl::extractGeometry() const {
    nlohmann::json geometry;
    
    auto root = parser_.get_root();
    if (!root) return geometry;
    
    // 查找 Maxwell2DModel 块
    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) {
        FEEM_WARN("未找到 Maxwell2DModel 块");
        return geometry;
    }
    
    // 提取 Name（设计名称）
    if (auto name_prop = maxwell_model->find_property("Name")) {
        if (std::holds_alternative<std::string>(name_prop->value)) {
            geometry["design_name"] = std::get<std::string>(name_prop->value);
        }
    }
    
    // 提取 SolutionType
    if (auto solution_type_prop = maxwell_model->find_property("SolutionType")) {
        if (std::holds_alternative<std::string>(solution_type_prop->value)) {
            geometry["solution_type"] = std::get<std::string>(solution_type_prop->value);
        }
    }
    
    // 提取 GeometryMode
    if (auto geometry_mode_prop = maxwell_model->find_property("GeometryMode")) {
        if (std::holds_alternative<std::string>(geometry_mode_prop->value)) {
            geometry["geometry_mode"] = std::get<std::string>(geometry_mode_prop->value);
        }
    }
    
    // 提取 ModelDepth
    if (auto model_depth_prop = maxwell_model->find_property("ModelDepth")) {
        if (std::holds_alternative<std::string>(model_depth_prop->value)) {
            geometry["model_depth"] = std::get<std::string>(model_depth_prop->value);
        }
    }
    
    // 提取 BackgroundMaterialName
    if (auto background_material_prop = maxwell_model->find_property("BackgroundMaterialName")) {
        if (std::holds_alternative<std::string>(background_material_prop->value)) {
            geometry["background_material"] = std::get<std::string>(background_material_prop->value);
        }
    }
    
    // 查找 GeometryCore 块以提取 Units
    auto geometry_core = findBlock("GeometryCore", maxwell_model);
    if (geometry_core) {
        if (auto units_prop = geometry_core->find_property("Units")) {
            if (std::holds_alternative<std::string>(units_prop->value)) {
                geometry["units"] = std::get<std::string>(units_prop->value);
            }
        }
    }
    
    FEEM_INFO("成功提取几何信息");
    return geometry;
}

bool MaxwellParserImpl::isParseTreeValid() const {
    return parser_.get_root() != nullptr && parser_.validate();
}

std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode> 
MaxwellParserImpl::findBlock(const std::string& block_name,
                            const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& start_node) const {
    if (!start_node) return nullptr;
    
    // 广度优先搜索
    if (start_node->name == block_name) {
        return start_node;
    }
    
    for (const auto& child : start_node->children) {
        if (child->name == block_name) {
            return child;
        }
    }
    
    // 递归搜索子块
    for (const auto& child : start_node->children) {
        auto result = findBlock(block_name, child);
        if (result) {
            return result;
        }
    }
    
    return nullptr;
}

nlohmann::json MaxwellParserImpl::propertyToJson(const fe_em::tool::maxwell_parser::Property& prop) const {
    nlohmann::json result;
    
    result["name"] = prop.name;
    result["type"] = static_cast<int>(prop.type);
    result["line_number"] = prop.line_number;
    
    // 根据类型转换值
    std::visit([&](const auto& value) {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, std::string>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, double>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, bool>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, std::vector<fe_em::tool::maxwell_parser::Value>>) {
            nlohmann::json array;
            for (const auto& item : value) {
                // 递归处理数组元素
                if (std::holds_alternative<double>(item)) {
                    array.push_back(std::get<double>(item));
                } else if (std::holds_alternative<std::string>(item)) {
                    array.push_back(std::get<std::string>(item));
                } else if (std::holds_alternative<bool>(item)) {
                    array.push_back(std::get<bool>(item));
                }
            }
            result["value"] = array;
        } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
            result["value"] = value;
        }
    }, prop.value);
    
    return result;
}

nlohmann::json MaxwellParserImpl::blockToJson(const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block) const {
    nlohmann::json result;
    
    if (!block) return result;
    
    result["name"] = block->name;
    result["start_line"] = block->start_line;
    result["end_line"] = block->end_line;
    
    // 转换属性
    nlohmann::json properties;
    for (const auto& prop : block->properties) {
        properties[prop.name] = propertyToJson(prop);
    }
    result["properties"] = properties;
    
    // 递归转换子块
    nlohmann::json children;
    for (const auto& child : block->children) {
        children.push_back(blockToJson(child));
    }
    result["children"] = children;
    
    return result;
}

std::string MaxwellParserImpl::extractPreviewImage() const {
    std::ifstream file(file_path_);
    if (!file.is_open()) {
        FEEM_WARN("无法打开AEDT文件: {}", file_path_);
        return "";
    }
    
    std::string line;
    std::string base64_data;
    bool in_image_section = false;  // 标记是否正在读取图像数据区域
    
    while (std::getline(file, line)) {
        // 检测Image64字段的起始位置
        if (line.find("Image64='") != std::string::npos) {
            in_image_section = true;
            // 提取起始行中的数据部分（跳过"Image64='"前缀）
            size_t start_pos = line.find("Image64='") + 9;
            base64_data = line.substr(start_pos);
            
            // 如果起始行以单引号结尾，说明数据在同一行内完成
            if (!base64_data.empty() && base64_data.back() == '\'') {
                base64_data.pop_back();
                break;
            }
            continue;
        }
        
        // 在图像数据区域内继续收集后续行的数据
        if (in_image_section) {
            for (char c : line) {
                // 单引号表示数据结束
                if (c == '\'') {
                    in_image_section = false;
                    break;
                }
                // 过滤掉转义字符和换行符
                if (c != '\\' && c != '\n' && c != '\r') {
                    base64_data += c;
                }
            }
            // 数据已完整提取
            if (!in_image_section) {
                break;
            }
        }
    }
    
    file.close();
    
    if (base64_data.empty()) {
        FEEM_DEBUG("AEDT文件中未找到预览图像数据");
    } else {
        FEEM_DEBUG("成功提取预览图像Base64数据，长度: {} 字符", base64_data.size());
    }
    
    return base64_data;
}

} // namespace tool