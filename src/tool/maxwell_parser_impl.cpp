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

    // 提取材料元信息
    material["library"] = safeGetString(material_block, "Library", "");
    material["lib_location"] = safeGetString(material_block, "LibLocation", "");
    material["mod_since_lib"] = safeGetBool(material_block, "ModSinceLib", false);
    material["mod_time"] = safeGetInt64(material_block, "ModTime", 0);

    // 提取坐标系统类型
    material["coordinate_system_type"] = safeGetString(material_block, "CoordinateSystemType", "");
    material["bulk_or_surface_type"] = safeGetInt(material_block, "BulkOrSurfaceType", 0);

    // 提取物理类型列表（PhysicsTypes函数调用格式）
    if (auto physics_prop = material_block->find_property("PhysicsTypes")) {
        std::vector<std::string> physics_types;
        if (std::holds_alternative<std::vector<std::string>>(physics_prop->value)) {
            physics_types = std::get<std::vector<std::string>>(physics_prop->value);
        }
        material["physics_types"] = physics_types;
    }

    // 提取B-H曲线元信息
    material["bh_curve_data_type"] = safeGetString(material_block, "permeability.BTypeForSingleCurve", "normal");
    material["h_unit"] = safeGetString(material_block, "permeability.HUnit", "A_per_meter");
    material["b_unit"] = safeGetString(material_block, "permeability.BUnit", "tesla");
    material["is_temperature_dependent"] = safeGetBool(material_block, "permeability.IsTemperatureDependent", false);

    // 提取铁损类型和堆叠类型
    material["core_loss_type_choice"] = safeGetString(material_block, "core_loss_type.Choice", "");
    material["stacking_type_str"] = safeGetString(material_block, "stacking_type.Choice", "");

    // 提取附加物理属性
    material["permittivity"] = safeGetDouble(material_block, "permittivity", 1.0);
    material["youngs_modulus"] = safeGetDouble(material_block, "youngs_modulus", 0.0);
    material["poisons_ratio"] = safeGetDouble(material_block, "poisons_ratio", 0.0);
    material["thermal_expansion_coefficient"] = safeGetDouble(material_block, "thermal_expansion_coefficient", 0.0);

    // 提取温度相关B-H曲线
    extractTemperatureBHCurves(material_block, material);

    // 提取CoreLoss多频率曲线数据
    extractCoreLossMulticurve(material_block, material);

    // 提取铁损系数设置数据
    extractCoeffSetup(material_block, material);

    // 提取热修正器数据
    extractThermalModifiers(material_block, material);

    // 提取外观数据
    extractMatAppearance(material_block, material);

    // 提取备注和关键词（在AttachedData子块内）
    auto attached_data = findBlock("AttachedData", material_block);
    if (attached_data) {
        auto notes_block = findBlock("MatNotesData", attached_data);
        if (notes_block) {
            material["notes"] = safeGetString(notes_block, "Notes", "");
        }
        auto keywords_block = findBlock("MatKeywordsData", attached_data);
        if (keywords_block) {
            material["keywords"] = safeGetString(keywords_block, "Keywords", "");
        }
    }

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

double MaxwellParserImpl::safeGetDouble(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& key,
    double default_val) const {

    if (!block) return default_val;
    auto prop = block->find_property(key);
    if (!prop) return default_val;
    return propertyValueToDouble(*prop);
}

int MaxwellParserImpl::safeGetInt(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& key,
    int default_val) const {

    if (!block) return default_val;
    auto prop = block->find_property(key);
    if (!prop) return default_val;
    const auto& value = prop->value;
    if (std::holds_alternative<double>(value)) {
        return static_cast<int>(std::get<double>(value));
    }
    return default_val;
}

std::string MaxwellParserImpl::safeGetString(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& key,
    const std::string& default_val) const {

    if (!block) return default_val;
    auto prop = block->find_property(key);
    if (!prop) return default_val;
    return propertyValueToString(*prop);
}

bool MaxwellParserImpl::safeGetBool(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& key,
    bool default_val) const {

    if (!block) return default_val;
    auto prop = block->find_property(key);
    if (!prop) return default_val;
    const auto& value = prop->value;
    if (std::holds_alternative<bool>(value)) {
        return std::get<bool>(value);
    }
    if (std::holds_alternative<double>(value)) {
        return std::get<double>(value) != 0.0;
    }
    if (std::holds_alternative<std::string>(value)) {
        std::string str = std::get<std::string>(value);
        if (str == "true") return true;
        if (str == "false") return false;
    }
    return default_val;
}

int64_t MaxwellParserImpl::safeGetInt64(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block,
    const std::string& key,
    int64_t default_val) const {

    if (!block) return default_val;
    auto prop = block->find_property(key);
    if (!prop) return default_val;
    const auto& value = prop->value;
    if (std::holds_alternative<double>(value)) {
        return static_cast<int64_t>(std::get<double>(value));
    }
    return default_val;
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

void MaxwellParserImpl::extractTemperatureBHCurves(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {

    auto temps_block = findBlock("Temperatures", material_block);
    if (!temps_block) return;

    nlohmann::json temp_curves = nlohmann::json::array();
    for (const auto& temp_block : temps_block->children) {
        nlohmann::json temp_entry;
        temp_entry["temp_ref"] = safeGetInt(temp_block, "TempRef", 20);

        // 提取该温度下的BHCoordinates曲线
        auto bh_coords = findBlock("BHCoordinates", temp_block);
        if (bh_coords) {
            if (auto points_prop = bh_coords->find_property("Points")) {
                auto curve = extractBHCurve(*points_prop);
                if (!curve.empty()) {
                    temp_entry["bh_coordinates"] = curve;
                }
            }
        }

        if (!temp_entry.empty()) {
            temp_curves.push_back(temp_entry);
        }
    }

    if (!temp_curves.empty()) {
        material["temperature_bh_curves"] = temp_curves;
    }
}

void MaxwellParserImpl::extractCoreLossMulticurve(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {

    auto attached_data = findBlock("AttachedData", material_block);
    if (!attached_data) return;

    auto multicurve_block = findBlock("CoreLossMultiCurveData", attached_data);
    if (!multicurve_block) return;

    nlohmann::json multicurve;
    multicurve["coreloss_unit"] = safeGetString(multicurve_block, "coreloss_unit", "w_per_kg");

    auto all_curves = findBlock("AllCurves", multicurve_block);
    if (all_curves) {
        nlohmann::json curves = nlohmann::json::array();
        for (const auto& one_curve : all_curves->children) {
            nlohmann::json curve_entry;
            curve_entry["frequency"] = safeGetInt(one_curve, "Frequency", 0);

            auto coords_block = findBlock("Coordinates", one_curve);
            if (coords_block) {
                if (auto points_prop = coords_block->find_property("Points")) {
                    auto bp_curve = extractBHCurve(*points_prop);
                    if (!bp_curve.empty()) {
                        curve_entry["coordinates"] = bp_curve;
                    }
                }
            }

            curves.push_back(curve_entry);
        }
        multicurve["curves"] = curves;
    }

    material["core_loss_multicurve"] = multicurve;
}

void MaxwellParserImpl::extractCoeffSetup(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {

    auto attached_data = findBlock("AttachedData", material_block);
    if (!attached_data) return;

    auto coeff_block = findBlock("CoefficientSetupData", attached_data);
    if (!coeff_block) return;

    nlohmann::json coeff_setup;
    coeff_setup["coefficient_setup"] = safeGetString(coeff_block, "coefficient_setup", "");
    coeff_setup["frequency"] = safeGetInt(coeff_block, "Frequency", 0);
    coeff_setup["thickness"] = safeGetString(coeff_block, "Thickness", "");
    coeff_setup["conductivity"] = safeGetString(coeff_block, "Conductivity", "");

    auto coords_block = findBlock("Coordinates", coeff_block);
    if (coords_block) {
        if (auto points_prop = coords_block->find_property("Points")) {
            auto bp_data = extractBHCurve(*points_prop);
            if (!bp_data.empty()) {
                coeff_setup["coordinates"] = bp_data;
            }
        }
    }

    material["coefficient_setup"] = coeff_setup;
}

void MaxwellParserImpl::extractThermalModifiers(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {

    auto modifier_data = findBlock("ModifierData", material_block);
    if (!modifier_data) return;

    nlohmann::json modifiers = nlohmann::json::array();
    for (const auto& mod_block : modifier_data->children) {
        nlohmann::json modifier;
        modifier["property"] = safeGetString(mod_block, "Property", "");
        modifier["index"] = safeGetInt(mod_block, "Index", 0);
        modifier["free_form_value"] = safeGetString(mod_block, "free_form_value", "");

        if (!modifier["property"].get<std::string>().empty()) {
            modifiers.push_back(modifier);
        }
    }

    if (!modifiers.empty()) {
        material["thermal_modifiers"] = modifiers;
    }
}

void MaxwellParserImpl::extractMatAppearance(
    const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& material_block,
    nlohmann::json& material) const {

    auto attached_data = findBlock("AttachedData", material_block);
    if (!attached_data) return;

    auto appearance_block = findBlock("MatAppearanceData", attached_data);
    if (!appearance_block) return;

    nlohmann::json appearance;
    appearance["red"] = safeGetInt(appearance_block, "Red", 200);
    appearance["green"] = safeGetInt(appearance_block, "Green", 200);
    appearance["blue"] = safeGetInt(appearance_block, "Blue", 200);
    appearance["transparency"] = safeGetDouble(appearance_block, "Transparency", 0.0);

    material["appearance"] = appearance;
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

        // 提取 ParentBndID（Coil边界指向Winding Group的ID）
        if (auto parent_prop = boundary_block->find_property("ParentBndID")) {
            if (std::holds_alternative<double>(parent_prop->value)) {
                boundary_json["parent_bnd_id"] = static_cast<int>(std::get<double>(parent_prop->value));
            }
        }

        // 提取 IsComponent
        if (auto comp_prop = boundary_block->find_property("IsComponent")) {
            if (std::holds_alternative<bool>(comp_prop->value)) {
                boundary_json["is_component"] = std::get<bool>(comp_prop->value);
            } else if (std::holds_alternative<double>(comp_prop->value)) {
                boundary_json["is_component"] = std::get<double>(comp_prop->value) != 0.0;
            }
        }

        // 提取 CoordinateSystem
        if (auto coord_prop = boundary_block->find_property("CoordinateSystem")) {
            if (std::holds_alternative<double>(coord_prop->value)) {
                boundary_json["coordinate_system"] = static_cast<int>(std::get<double>(coord_prop->value));
            }
        }

        // 提取 Type（绕组组特有：Current/Voltage）
        if (auto type_prop = boundary_block->find_property("Type")) {
            if (std::holds_alternative<std::string>(type_prop->value)) {
                boundary_json["type"] = std::get<std::string>(type_prop->value);
            }
        }

        // 提取 IsSolid（绕组组特有）
        if (auto solid_prop = boundary_block->find_property("IsSolid")) {
            if (std::holds_alternative<bool>(solid_prop->value)) {
                boundary_json["is_solid"] = std::get<bool>(solid_prop->value);
            } else if (std::holds_alternative<double>(solid_prop->value)) {
                boundary_json["is_solid"] = std::get<double>(solid_prop->value) != 0.0;
            }
        }

        // 提取 Current（绕组组电流表达式）
        if (auto current_prop = boundary_block->find_property("Current")) {
            if (std::holds_alternative<std::string>(current_prop->value)) {
                boundary_json["current"] = std::get<std::string>(current_prop->value);
            }
        }

        // 提取 Voltage（绕组组电压表达式）
        if (auto voltage_prop = boundary_block->find_property("Voltage")) {
            if (std::holds_alternative<std::string>(voltage_prop->value)) {
                boundary_json["voltage"] = std::get<std::string>(voltage_prop->value);
            }
        }

        // 提取 Resistance（绕组组电阻）
        if (auto resist_prop = boundary_block->find_property("Resistance")) {
            if (std::holds_alternative<std::string>(resist_prop->value)) {
                boundary_json["resistance"] = std::get<std::string>(resist_prop->value);
            }
        }

        // 提取 Inductance（绕组组电感）
        if (auto induct_prop = boundary_block->find_property("Inductance")) {
            if (std::holds_alternative<std::string>(induct_prop->value)) {
                boundary_json["inductance"] = std::get<std::string>(induct_prop->value);
            }
        }

        // 提取 ParallelBranchesNum（并联支路数）
        if (auto branch_prop = boundary_block->find_property("ParallelBranchesNum")) {
            if (std::holds_alternative<std::string>(branch_prop->value)) {
                boundary_json["parallel_branches_num"] = std::get<std::string>(branch_prop->value);
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

    // 补充瞬态求解器字段
    if (auto enabled_prop = setup_block->find_property("Enabled")) {
        if (std::holds_alternative<bool>(enabled_prop->value)) {
            setup["enabled"] = std::get<bool>(enabled_prop->value);
        }
    }

    // 时间积分方法
    if (auto time_int_prop = setup_block->find_property("TimeIntegrationMethod")) {
        if (std::holds_alternative<std::string>(time_int_prop->value)) {
            setup["time_integration_method"] = std::get<std::string>(time_int_prop->value);
        }
    }

    // 平滑BH曲线
    if (auto smooth_bh_prop = setup_block->find_property("SmoothBHCurve")) {
        if (std::holds_alternative<bool>(smooth_bh_prop->value)) {
            setup["smooth_bh_curve"] = std::get<bool>(smooth_bh_prop->value);
        } else if (std::holds_alternative<double>(smooth_bh_prop->value)) {
            setup["smooth_bh_curve"] = std::get<double>(smooth_bh_prop->value) != 0.0;
        }
    }

    // 保存场类型
    if (auto save_fields_prop = setup_block->find_property("SaveFieldsType")) {
        if (std::holds_alternative<std::string>(save_fields_prop->value)) {
            setup["save_fields_type"] = std::get<std::string>(save_fields_prop->value);
        }
    }

    // 自适应时间步相关字段
    if (auto max_iter_prop = setup_block->find_property("MaxAdaptiveStepIterations")) {
        if (std::holds_alternative<double>(max_iter_prop->value)) {
            setup["max_adaptive_step_iterations"] = static_cast<int>(std::get<double>(max_iter_prop->value));
        }
    }

    if (auto error_tol_prop = setup_block->find_property("AdaptiveStepErrorTolerance")) {
        if (std::holds_alternative<std::string>(error_tol_prop->value)) {
            setup["adaptive_step_error_tolerance"] = std::get<std::string>(error_tol_prop->value);
        }
    }

    if (auto min_num_steps_prop = setup_block->find_property("MinNumberOfStepsPerCycle")) {
        if (std::holds_alternative<double>(min_num_steps_prop->value)) {
            setup["min_steps_per_cycle"] = static_cast<int>(std::get<double>(min_num_steps_prop->value));
        }
    }

    if (auto max_num_steps_prop = setup_block->find_property("MaxNumberOfStepsPerCycle")) {
        if (std::holds_alternative<double>(max_num_steps_prop->value)) {
            setup["max_steps_per_cycle"] = static_cast<int>(std::get<double>(max_num_steps_prop->value));
        }
    }

    // 输出设置
    if (auto save_each_prop = setup_block->find_property("SaveEachNthStep")) {
        if (std::holds_alternative<double>(save_each_prop->value)) {
            setup["save_each_nth_step"] = static_cast<int>(std::get<double>(save_each_prop->value));
        }
    }

    if (auto save_at_time_prop = setup_block->find_property("SaveSpecificTimePoints")) {
        if (std::holds_alternative<std::string>(save_at_time_prop->value)) {
            setup["save_specific_time_points"] = std::get<std::string>(save_at_time_prop->value);
        }
    }

    // 非线性求解器设置
    if (auto max_nonlinear_prop = setup_block->find_property("MaxNonlinearIterations")) {
        if (std::holds_alternative<double>(max_nonlinear_prop->value)) {
            setup["max_nonlinear_iterations"] = static_cast<int>(std::get<double>(max_nonlinear_prop->value));
        }
    }

    if (auto derating_prop = setup_block->find_property("DerivativeMaxNormForConvergence")) {
        if (std::holds_alternative<std::string>(derating_prop->value)) {
            setup["derivative_max_norm"] = std::get<std::string>(derating_prop->value);
        }
    }

    if (auto conv_order_prop = setup_block->find_property("ConvergenceOrder")) {
        if (std::holds_alternative<double>(conv_order_prop->value)) {
            setup["convergence_order"] = std::get<double>(conv_order_prop->value);
        }
    }

    // 运动耦合设置
    if (auto motion_coupling_prop = setup_block->find_property("EnableMotionCoupling")) {
        if (std::holds_alternative<bool>(motion_coupling_prop->value)) {
            setup["enable_motion_coupling"] = std::get<bool>(motion_coupling_prop->value);
        } else if (std::holds_alternative<double>(motion_coupling_prop->value)) {
            setup["enable_motion_coupling"] = std::get<double>(motion_coupling_prop->value) != 0.0;
        }
    }

    // 温度耦合设置
    if (auto temp_coupling_prop = setup_block->find_property("TemperatureCoupling")) {
        if (std::holds_alternative<std::string>(temp_coupling_prop->value)) {
            setup["temperature_coupling"] = std::get<std::string>(temp_coupling_prop->value);
        }
    }

    // 求解顺序
    if (auto solve_order_prop = setup_block->find_property("SolveOrder")) {
        if (std::holds_alternative<std::string>(solve_order_prop->value)) {
            setup["solve_order"] = std::get<std::string>(solve_order_prop->value);
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

std::vector<nlohmann::json> MaxwellParserImpl::extractWindings() const {
    std::vector<nlohmann::json> windings;

    auto root = parser_.get_root();
    if (!root) return windings;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return windings;

    auto boundary_setup = findBlock("BoundarySetup", maxwell_model);
    if (!boundary_setup) return windings;

    auto boundaries_block = findBlock("Boundaries", boundary_setup);
    if (!boundaries_block) return windings;

    for (const auto& bnd_block : boundaries_block->children) {
        // 筛选 BoundType == "Winding Group" 的条目
        std::string bound_type = safeGetString(bnd_block, "BoundType", "");
        if (bound_type != "Winding Group") continue;

        nlohmann::json winding;
        winding["name"] = bnd_block->name;
        winding["id"] = safeGetInt(bnd_block, "ID", 0);
        winding["bound_type"] = bound_type;
        winding["type"] = safeGetString(bnd_block, "Type", "");
        winding["is_solid"] = safeGetBool(bnd_block, "IsSolid", false);
        winding["current"] = safeGetString(bnd_block, "Current", "");
        winding["voltage"] = safeGetString(bnd_block, "Voltage", "");
        winding["resistance"] = safeGetString(bnd_block, "Resistance", "");
        winding["inductance"] = safeGetString(bnd_block, "Inductance", "");
        winding["parallel_branches_num"] = safeGetString(bnd_block, "ParallelBranchesNum", "");

        windings.push_back(winding);
    }

    FEEM_INFO("成功提取 {} 个绕组组", windings.size());
    return windings;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractMotionSetups() const {
    std::vector<nlohmann::json> motion_setups;

    auto root = parser_.get_root();
    if (!root) return motion_setups;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return motion_setups;

    auto motion_list = findBlock("MotionSetupList", maxwell_model);
    if (!motion_list) return motion_setups;

    for (const auto& child : motion_list->children) {
        const std::string& child_name = child->name;

        // MotionSetup* 条目：运动参数定义（排除Moving*条目）
        if (child_name.find("MotionSetup") != std::string::npos &&
            child_name.find("Moving") == std::string::npos) {

            nlohmann::json motion;
            motion["name"] = child_name;
            motion["motion_type"] = safeGetString(child, "MotionType", "");
            motion["move_type"] = safeGetString(child, "MoveType", "");
            motion["coordinate_system"] = safeGetInt(child, "CoordinateSystem", 0);
            motion["axis"] = safeGetString(child, "Axis", "");
            motion["is_positive"] = safeGetBool(child, "IsPositive", true);
            motion["init_pos"] = safeGetString(child, "InitPos", "");
            motion["has_rotate_limit"] = safeGetBool(child, "HasRotateLimit", false);
            motion["non_cylindrical"] = safeGetBool(child, "NonCylindrical", false);
            motion["consider_mechanical_transient"] = safeGetBool(
                child, "ConsiderMechanicalTransient", false);
            motion["angular_velocity"] = safeGetString(child, "AngularVelocity", "");

            // 提取运动对象列表
            if (auto obj_prop = child->find_property("Objects")) {
                if (std::holds_alternative<std::vector<double>>(obj_prop->value)) {
                    auto objs = std::get<std::vector<double>>(obj_prop->value);
                    std::vector<int> obj_ids;
                    for (double o : objs) obj_ids.push_back(static_cast<int>(o));
                    motion["objects"] = obj_ids;
                }
            }

            motion_setups.push_back(motion);
        }

        // Moving* 条目：运动对象关联信息
        if (child_name.find("Moving") != std::string::npos &&
            child_name.find("MotionSetup") == std::string::npos) {

            nlohmann::json moving;
            moving["name"] = child_name;
            moving["motion_type"] = safeGetString(child, "MotionType", "");
            moving["band_name_ref"] = safeGetInt(child, "BandName", 0);

            if (auto obj_prop = child->find_property("Objects")) {
                if (std::holds_alternative<std::vector<double>>(obj_prop->value)) {
                    auto objs = std::get<std::vector<double>>(obj_prop->value);
                    std::vector<int> obj_ids;
                    for (double o : objs) obj_ids.push_back(static_cast<int>(o));
                    moving["moving_objects"] = obj_ids;
                }
            }

            motion_setups.push_back(moving);
        }
    }

    FEEM_INFO("成功提取 {} 个运动设置项", motion_setups.size());
    return motion_setups;
}

nlohmann::json MaxwellParserImpl::extractGlobalMeshSettings() const {
    nlohmann::json settings;

    auto root = parser_.get_root();
    if (!root) return settings;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return settings;

    auto mesh_setup = findBlock("MeshSetup", maxwell_model);
    if (!mesh_setup) return settings;

    auto mesh_settings = findBlock("MeshSettings", mesh_setup);
    if (!mesh_settings) return settings;

    // 全局曲面近似设置
    auto surf_approx = findBlock("GlobalSurfApproximation", mesh_settings);
    if (surf_approx) {
        settings["surf_approx_choice"] = safeGetString(surf_approx, "SurfApproxChoice", "");
        settings["slider_mesh_settings"] = safeGetInt(surf_approx, "SliderMeshSettings", 5);
    }

    // 全局模型分辨率
    auto model_res = findBlock("GlobalModelRes", mesh_settings);
    if (model_res) {
        settings["use_auto_length"] = safeGetBool(model_res, "UseAutoLength", true);
    }

    FEEM_DEBUG("成功提取全局网格设置");
    return settings;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractMeshOperations() const {
    std::vector<nlohmann::json> operations;

    auto root = parser_.get_root();
    if (!root) return operations;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return operations;

    auto mesh_setup = findBlock("MeshSetup", maxwell_model);
    if (!mesh_setup) return operations;

    auto mesh_ops = findBlock("MeshOperations", mesh_setup);
    if (!mesh_ops) return operations;

    for (const auto& op_block : mesh_ops->children) {
        nlohmann::json op;
        op["name"] = op_block->name;
        op["id"] = safeGetInt(op_block, "ID", 0);
        op["type"] = safeGetString(op_block, "Type", "");
        op["enabled"] = safeGetBool(op_block, "Enabled", true);
        op["is_component"] = safeGetBool(op_block, "IsComponent", false);
        op["is_global"] = safeGetBool(op_block, "IsGlobal", false);

        // 提取对象列表
        if (auto obj_prop = op_block->find_property("Objects")) {
            if (std::holds_alternative<std::vector<double>>(obj_prop->value)) {
                auto objs = std::get<std::vector<double>>(obj_prop->value);
                std::vector<int> obj_ids;
                for (double o : objs) obj_ids.push_back(static_cast<int>(o));
                op["objects"] = obj_ids;
            }
        }

        // 根据类型提取特有字段
        std::string op_type = op["type"].get<std::string>();

        if (op_type == "SurfApproxBased") {
            op["refine_inside"] = safeGetBool(op_block, "RefineInside", true);
            op["surf_approx_mode"] = safeGetString(op_block, "SurfApproxMode", "");
            op["surf_dev_choice"] = safeGetInt(op_block, "SurfDevChoice", 0);
            op["surf_dev"] = safeGetString(op_block, "SurfDev", "");
            op["normal_dev_choice"] = safeGetInt(op_block, "NormalDevChoice", 0);
            op["normal_dev"] = safeGetString(op_block, "NormalDev", "");
            op["aspect_ratio_choice"] = safeGetInt(op_block, "AspectRatioChoice", 0);

        } else if (op_type == "LengthBased") {
            op["restrict_length"] = safeGetBool(op_block, "RestrictLength", true);
            op["max_length"] = safeGetString(op_block, "MaxLength", "");
            op["apply_to_initial_mesh"] = safeGetBool(op_block, "ApplyToInitialMesh", false);

        } else if (op_type == "CylindricalGap") {
            op["use_band_mapping_angle"] = safeGetBool(op_block, "UseBandMappingAngle", true);
            op["band_mapping_angle"] = safeGetString(op_block, "BandMappingAngle", "");
        }

        operations.push_back(op);
    }

    FEEM_INFO("成功提取 {} 个网格操作", operations.size());
    return operations;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractDesignVariables() const {
    std::vector<nlohmann::json> variables;

    auto root = parser_.get_root();
    if (!root) return variables;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return variables;

    auto model_setup = findBlock("ModelSetup", maxwell_model);
    if (!model_setup) return variables;

    auto properties = findBlock("Properties", model_setup);
    if (!properties) return variables;

    // VariableProp 是一个函数调用类型的属性，其值包含所有变量定义
    if (auto varprop = properties->find_property("VariableProp")) {
        const auto& value = varprop->value;
        if (std::holds_alternative<std::vector<std::string>>(value)) {
            auto entries = std::get<std::vector<std::string>>(value);
            for (size_t i = 0; i + 1 < entries.size(); i += 2) {
                nlohmann::json var;
                var["name"] = entries[i];
                var["expression"] = entries[i + 1];
                variables.push_back(var);
            }
        }
    }

    FEEM_INFO("成功提取 {} 个设计变量", variables.size());
    return variables;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractNonIndexedVariables() const {
    std::vector<nlohmann::json> variables;

    auto root = parser_.get_root();
    if (!root) return variables;

    auto nonindexed = findBlock("NonIndexedVariables", root);
    if (!nonindexed) return variables;

    for (const auto& var_block : nonindexed->children) {
        nlohmann::json var;
        var["name"] = var_block->name;
        var["prop_type"] = safeGetString(var_block, "PropType", "");
        var["expression"] = safeGetString(var_block, "Expression", "");

        variables.push_back(var);
    }

    FEEM_INFO("成功提取 {} 个非索引变量", variables.size());
    return variables;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractOutputVariables() const {
    std::vector<nlohmann::json> output_vars;

    auto root = parser_.get_root();
    if (!root) return output_vars;

    auto output_var = findBlock("OutputVariable", root);
    if (!output_var) return output_vars;

    auto output_variables = findBlock("OutputVariables", output_var);
    if (!output_variables) return output_vars;

    for (const auto& var_block : output_variables->children) {
        nlohmann::json var;
        var["name"] = var_block->name;
        var["id"] = safeGetInt(var_block, "ID", 0);
        var["expression"] = safeGetString(var_block, "Expression", "");
        var["result_unit"] = safeGetString(var_block, "ResultUnit", "");

        output_vars.push_back(var);
    }

    FEEM_INFO("成功提取 {} 个输出变量", output_vars.size());
    return output_vars;
}

nlohmann::json MaxwellParserImpl::extractTemperatureSettings() const {
    nlohmann::json temp_settings;

    auto root = parser_.get_root();
    if (!root) return temp_settings;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return temp_settings;

    auto temp_settings_block = findBlock("TemperatureSettings", maxwell_model);
    if (!temp_settings_block) return temp_settings;

    temp_settings["include_temperature_dependence"] = safeGetBool(
        temp_settings_block, "IncludeTemperatureDependence", false);
    temp_settings["enable_feedback"] = safeGetBool(
        temp_settings_block, "EnableFeedback", false);

    // 提取温度映射表 Temperatures(1='22cel', 2='22cel', ...)
    if (auto temps_prop = temp_settings_block->find_property("Temperatures")) {
        const auto& value = temps_prop->value;
        if (std::holds_alternative<std::vector<std::string>>(value)) {
            auto entries = std::get<std::vector<std::string>>(value);
            nlohmann::json mapping;
            for (size_t i = 0; i + 1 < entries.size(); i += 2) {
                mapping[entries[i]] = entries[i + 1];
            }
            temp_settings["temperature_mapping"] = mapping;
        }
    }

    FEEM_DEBUG("成功提取温度设置");
    return temp_settings;
}

nlohmann::json MaxwellParserImpl::extractGlobalBoundData() const {
    nlohmann::json global_data;

    auto root = parser_.get_root();
    if (!root) return global_data;

    auto maxwell_model = findBlock("Maxwell2DModel", root);
    if (!maxwell_model) return global_data;

    auto boundary_setup = findBlock("BoundarySetup", maxwell_model);
    if (!boundary_setup) return global_data;

    // 提取 BoundarySetup 顶层属性
    global_data["name"] = safeGetString(boundary_setup, "Name", "");
    global_data["id"] = safeGetInt(boundary_setup, "ID", 0);

    FEEM_DEBUG("成功提取全局边界数据");
    return global_data;
}

} // namespace tool