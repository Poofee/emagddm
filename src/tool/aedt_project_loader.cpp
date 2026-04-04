/**
 * @file aedt_project_loader.cpp
 * @brief AEDT格式项目加载器实现
 * @details 将AEDT文件解析为统一的项目数据结构（ProjectManager）
 *          转换流程：MaxwellParserImpl解析 → JSON中间数据 → Material/Boundary/Excitation/SolutionSetup对象 → ProjectManager
 * @author Poofee
 * @date 2026-04-04
 */

#include "aedt_project_loader.hpp"
#include "project_manager.hpp"
#include "logger_factory.hpp"

#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace tool {

namespace {

/**
 * @brief 从JSON值安全提取double数值
 * @details 兼容number和string两种类型（AEDT文件中数值常存储为字符串）
 * @param json_val JSON值
 * @param default_val 解析失败时的默认值
 * @return 提取的double值
 */
double safeGetDouble(const nlohmann::json& json_val, double default_val = 0.0) {
    if (json_val.is_number()) {
        return json_val.get<double>();
    } else if (json_val.is_string()) {
        try {
            return std::stod(json_val.get<std::string>());
        } catch (...) {
            return default_val;
        }
    }
    return default_val;
}

} // anonymous namespace

AedtProjectLoader::AedtProjectLoader() {
}

bool AedtProjectLoader::canLoad(const std::string& file_path) const {
    if (file_path.empty()) return false;
    
    size_t pos = file_path.find_last_of('.');
    if (pos == std::string::npos) return false;
    
    std::string ext = file_path.substr(pos);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    return (ext == ".aedt");
}

bool AedtProjectLoader::load(ProjectManager& manager) {
    FEEM_INFO("========================================");
    FEEM_INFO("开始加载AEDT项目文件: {}", file_path_);
    FEEM_INFO("========================================");

    if (file_path_.empty()) {
        last_error_ = "文件路径未设置";
        FEEM_ERROR(last_error_);
        return false;
    }

    if (!fs::exists(file_path_)) {
        last_error_ = "文件不存在: " + file_path_;
        FEEM_ERROR(last_error_);
        return false;
    }

    if (!parser_.canParse(file_path_)) {
        last_error_ = "不支持的AEDT文件格式";
        FEEM_ERROR(last_error_);
        return false;
    }

    FEEM_INFO("");
    FEEM_INFO("[1/6] 解析文件信息...");
    auto file_info = parser_.parseFileInfo();
    manager.setProjectName(file_info.project_name);
    FEEM_INFO("  项目名称: {}", file_info.project_name);
    FEEM_INFO("  Maxwell版本: {}", file_info.maxwell_version);
    FEEM_INFO("  仿真类型: {}", file_info.simulation_type);
    FEEM_INFO("  维度: {}", file_info.dimension);

    FEEM_INFO("");
    FEEM_INFO("[2/6] 解析材料数据...");
    auto materials_json = parser_.parseMaterials();
    convertAndAddMaterials(manager, materials_json);

    FEEM_INFO("");
    FEEM_INFO("[3/6] 解析边界条件...");
    auto boundaries_json = parser_.parseBoundaries();
    convertAndAddBoundaries(manager, boundaries_json);

    FEEM_INFO("");
    FEEM_INFO("[4/6] 解析激励源...");
    auto excitations_json = parser_.parseExcitations();
    convertAndAddExcitations(manager, excitations_json);

    FEEM_INFO("");
    FEEM_INFO("[5/13] 解析求解设置...");
    auto setup_json = parser_.parseSolutionSetup();
    convertAndAddSolutionSetup(manager, setup_json);

    FEEM_INFO("");
    FEEM_INFO("[6/13] 解析绕组组...");
    auto windings_json = parser_.extractWindings();
    convertAndAddWindings(manager, windings_json);

    FEEM_INFO("");
    FEEM_INFO("[7/13] 解析运动设置...");
    auto motions_json = parser_.extractMotionSetups();
    convertAndAddMotionSetups(manager, motions_json);

    FEEM_INFO("");
    FEEM_INFO("[8/13] 解析全局网格设置...");
    auto mesh_settings_json = parser_.extractGlobalMeshSettings();

    FEEM_INFO("");
    FEEM_INFO("[9/13] 解析网格操作...");
    auto mesh_ops_json = parser_.extractMeshOperations();
    convertAndAddMeshOperations(manager, mesh_ops_json);

    FEEM_INFO("");
    FEEM_INFO("[10/13] 解析设计变量和输出变量...");
    auto design_vars_json = parser_.extractDesignVariables();
    auto nonindexed_vars_json = parser_.extractNonIndexedVariables();
    convertAndAddDesignVariables(manager, design_vars_json);
    convertAndAddDesignVariables(manager, nonindexed_vars_json);
    auto output_vars_json = parser_.extractOutputVariables();
    convertAndAddOutputVariables(manager, output_vars_json);

    FEEM_INFO("");
    FEEM_INFO("[11/13] 解析温度设置...");
    auto temp_settings_json = parser_.extractTemperatureSettings();
    convertAndAddTemperatureSettings(manager, temp_settings_json);

    FEEM_INFO("");
    FEEM_INFO("[12/13] 解析全局边界数据...");
    auto global_bound_data = parser_.extractGlobalBoundData();

    FEEM_INFO("");
    FEEM_INFO("[13/13] 提取预览图像...");
    extractPreviewImage(file_path_);

    FEEM_INFO("");
    FEEM_INFO("========================================");
    FEEM_INFO("AEDT项目数据加载完成");
    FEEM_INFO("========================================");

    return true;
}

void AedtProjectLoader::convertAndAddMaterials(ProjectManager& manager,
                                                 const std::vector<nlohmann::json>& materials_json) {
    for (const auto& mat_json : materials_json) {
        std::string name = mat_json.value("name", "Unknown");

        auto material = std::make_shared<Material>(name);

        if (mat_json.contains("permeability") && (mat_json["permeability"].is_number() || mat_json["permeability"].is_string())) {
            double mu_r = safeGetDouble(mat_json["permeability"]);
            material->setRelativePermeability(mu_r);
            material->setType(MatType::LINEAR_ISOTROPIC);
        }

        if (mat_json.contains("permeability_type") && mat_json["permeability_type"] == "nonlinear") {
            material->setType(MatType::NONLINEAR_ISOTROPIC);
            material->setBHCurveType(BHCurveType::SINGLE_CURVE);

            if (mat_json.contains("bh_curve") && mat_json["bh_curve"].is_array()) {
                std::vector<BHDataPoint> bh_points;
                for (const auto& point : mat_json["bh_curve"]) {
                    if (point.is_array() && point.size() >= 2) {
                        BHDataPoint dp;
                        dp.h = point[0].get<double>();
                        dp.b = point[1].get<double>();
                        bh_points.push_back(dp);
                    }
                }
                if (!bh_points.empty()) {
                    material->setBHCurve(bh_points);
                }
            }
        }

        if (mat_json.contains("conductivity") && (mat_json["conductivity"].is_number() || mat_json["conductivity"].is_string())) {
            double sigma = safeGetDouble(mat_json["conductivity"]);
            material->setConductivity(sigma);
        }

        if (mat_json.contains("mass_density") && (mat_json["mass_density"].is_number() || mat_json["mass_density"].is_string())) {
            double rho = safeGetDouble(mat_json["mass_density"]);
            material->setMassDensity(rho);
        }

        if (mat_json.contains("coercivity") && mat_json["coercivity"].is_object()) {
            const auto& coercivity = mat_json["coercivity"];
            if (coercivity.contains("magnitude")) {
                MaterialProperty prop;
                prop.name = "coercivity_magnitude";
                prop.type = "scalar";
                prop.value = safeGetDouble(coercivity["magnitude"]);
                prop.unit = "A/m";
                material->addProperty(prop);
            }
        }

        // 元信息（材料库来源等）
        if (mat_json.contains("library")) {
            material->setLibraryInfo(
                mat_json.value("library", ""),
                mat_json.value("lib_location", ""),
                mat_json.value("mod_since_lib", false));
        }

        // 扩展物理属性
        if (mat_json.contains("permittivity"))
            material->setPermittivity(safeGetDouble(mat_json["permittivity"], 1.0));
        if (mat_json.contains("youngs_modulus"))
            material->setYoungsModulus(safeGetDouble(mat_json["youngs_modulus"], 0.0));
        if (mat_json.contains("poisons_ratio"))
            material->setPoisonsRatio(safeGetDouble(mat_json["poisons_ratio"], 0.0));
        if (mat_json.contains("thermal_expansion_coefficient"))
            material->setThermalExpansionCoefficient(
                safeGetDouble(mat_json["thermal_expansion_coefficient"], 0.0));

        // B-H曲线元信息
        if (mat_json.contains("bh_curve_data_type"))
            material->setBHCurveDataType(mat_json["bh_curve_data_type"].get<std::string>());
        if (mat_json.contains("is_temperature_dependent"))
            material->setIsTemperatureDependent(mat_json["is_temperature_dependent"].get<bool>());

        // 铁损/堆叠类型
        if (mat_json.contains("core_loss_type_choice"))
            material->setCoreLossTypeChoice(mat_json["core_loss_type_choice"].get<std::string>());
        if (mat_json.contains("stacking_type_str"))
            material->setStackingTypeStr(mat_json["stacking_type_str"].get<std::string>());

        // 温度相关B-H曲线
        if (mat_json.contains("temperature_bh_curves") && mat_json["temperature_bh_curves"].is_object()) {
            std::map<double, std::vector<BHDataPoint>> temp_curves;
            for (auto& [temp_key, points_arr] : mat_json["temperature_bh_curves"].items()) {
                double temp = 0;
                try {
                    temp = std::stod(temp_key);
                } catch (const std::exception& e) {
                    FEEM_WARN("Invalid temperature key '{}': {}", temp_key, e.what());
                    continue;
                }
                std::vector<BHDataPoint> bh_points;
                for (const auto& pt : points_arr) {
                    BHDataPoint dp;
                    dp.h = pt[0].get<double>();
                    dp.b = pt[1].get<double>();
                    bh_points.push_back(dp);
                }
                temp_curves[temp] = bh_points;
            }
            material->setTemperatureBHCurves(temp_curves);
        }

        // 多频率铁损曲线
        if (mat_json.contains("core_loss_freq_curves") && mat_json["core_loss_freq_curves"].is_object()) {
            std::map<double, std::vector<std::pair<double, double>>> freq_curves;
            for (auto& [freq_key, bp_arr] : mat_json["core_loss_freq_curves"].items()) {
                double freq = 0;
                try {
                    freq = std::stod(freq_key);
                } catch (const std::exception& e) {
                    FEEM_WARN("Invalid frequency key '{}': {}", freq_key, e.what());
                    continue;
                }
                std::vector<std::pair<double, double>> bp_points;
                for (const auto& pt : bp_arr) {
                    bp_points.push_back({pt[0].get<double>(), pt[1].get<double>()});
                }
                freq_curves[freq] = bp_points;
            }
            material->setCoreLossFreqCurves(freq_curves);
        }

        // 热修正器
        if (mat_json.contains("thermal_modifiers") && mat_json["thermal_modifiers"].is_array()) {
            for (const auto& tm : mat_json["thermal_modifiers"]) {
                ThermalModifier mod;
                mod.property_name = tm.value("property_name", "");
                mod.index = tm.value("index", 0);
                mod.formula_string = tm.value("formula_string", "");
                material->addThermalModifier(mod);
            }
        }

        // 外观颜色数据
        if (mat_json.contains("appearance_data") && mat_json["appearance_data"].is_object()) {
            auto& app = mat_json["appearance_data"];
            material->setAppearance(app.value("red", 200), app.value("green", 200),
                                    app.value("blue", 200), app.value("transparency", 0.0));
        }

        manager.addMaterial(material);
        FEEM_INFO("  已加载材料: {} (μr={})", name, material->getRelativePermeability());
    }

    FEEM_INFO("  共加载 {} 个材料", materials_json.size());
}

void AedtProjectLoader::convertAndAddBoundaries(ProjectManager& manager,
                                                  const std::vector<nlohmann::json>& boundaries_json) {
    for (const auto& bnd_json : boundaries_json) {
        std::string name = bnd_json.value("name", "Unknown");

        auto boundary = std::make_shared<Boundary>(name);

        std::string bound_type_str = bnd_json.value("bound_type", "");
        BndType bnd_type = BndType::DIRICHLET;

        if (bound_type_str == "Balloon" || bound_type_str == "Balloon") {
            bnd_type = BndType::BALLOON;
        } else if (bound_type_str == "Perfect E") {
            bnd_type = BndType::PERFECT_E;
        } else if (bound_type_str == "Perfect H") {
            bnd_type = BndType::PERFECT_H;
        } else if (bound_type_str == "Symmetry") {
            bnd_type = BndType::EVEN_SYMMETRY;
        } else if (bound_type_str == "MasterSlave" || bound_type_str == "Independent") {
            bnd_type = BndType::MASTER_SLAVE;
        } else if (bound_type_str == "Slave" || bound_type_str == "Dependent") {
            bnd_type = BndType::MASTER_SLAVE;
        }

        boundary->setType(bnd_type);

        if (bnd_json.contains("edges") && bnd_json["edges"].is_array()) {
            for (const auto& edge : bnd_json["edges"]) {
                boundary->addEdge(std::to_string(edge.get<int>()));
            }
        }

        if (bnd_json.contains("objects") && bnd_json["objects"].is_array()) {
            for (const auto& obj : bnd_json["objects"]) {
                boundary->addObject(std::to_string(obj.get<int>()));
            }
        }

        if (bnd_json.contains("value") && bnd_json["value"].is_string()) {
            std::string value_str = bnd_json["value"].get<std::string>();
            try {
                boundary->setCurrent(std::stod(value_str));
            } catch (...) {
                boundary->setVoltage(std::stod(value_str));
            }
        }

        // 边界扩展字段
        if (bnd_json.contains("parent_bnd_id"))
            boundary->setParentBndID(bnd_json["parent_bnd_id"].get<int>());
        if (bnd_json.contains("is_component"))
            boundary->setIsComponent(bnd_json["is_component"].get<bool>());
        if (bnd_json.contains("coordinate_system"))
            boundary->setCoordinateSystem(bnd_json["coordinate_system"].get<int>());
        if (bnd_json.contains("conductor_number_str"))
            boundary->setConductorNumberStr(bnd_json["conductor_number_str"].get<std::string>());

        manager.addBoundary(boundary);
        FEEM_INFO("  已加载边界条件: {} ({})", name, bound_type_str);
    }

    FEEM_INFO("  共加载 {} 个边界条件", boundaries_json.size());
}

void AedtProjectLoader::convertAndAddExcitations(ProjectManager& manager,
                                                   const std::vector<nlohmann::json>& excitations_json) {
    for (const auto& exc_json : excitations_json) {
        std::string name = exc_json.value("name", "Unknown");

        auto excitation = std::make_shared<Excitation>(name);

        std::string type_str = exc_json.value("type", "Coil");

        if (type_str == "Current" || type_str == "Coil") {
            excitation->setType(ExcitationType::CURRENT_SOURCE);
        } else if (type_str == "Voltage") {
            excitation->setType(ExcitationType::VOLTAGE_SOURCE);
        } else if (type_str == "Current Density") {
            excitation->setType(ExcitationType::CURRENT_DENSITY);
        } else {
            excitation->setType(ExcitationType::COIL);
        }

        if (exc_json.contains("winding") && exc_json["winding"].is_number()) {
            int winding_num = exc_json["winding"].get<int>();
            excitation->setCoilGroup(std::to_string(winding_num));
        }

        if (exc_json.contains("conductor_number") && exc_json["conductor_number"].is_string()) {
            std::string conductor_str = exc_json["conductor_number"].get<std::string>();
            try {
                double current_val = std::stod(conductor_str);
                excitation->setValue(current_val);
            } catch (...) {
            }
        }

        manager.addExcitation(excitation);
        FEEM_INFO("  已加载激励源: {} ({}, {} A)", name, type_str, excitation->getValue());
    }

    FEEM_INFO("  共加载 {} 个激励源", excitations_json.size());
}

void AedtProjectLoader::convertAndAddSolutionSetup(ProjectManager& manager,
                                                     const nlohmann::json& setup_json) {
    std::string name = setup_json.value("name", "Setup1");

    auto setup = std::make_shared<SolutionSetup>(name);

    std::string setup_type_str = setup_json.value("setup_type", "");

    SimulationType sim_type = SimulationType::MAGNETOSTATIC;
    if (setup_type_str.find("EddyCurrent") != std::string::npos ||
        setup_type_str.find("Eddy Current") != std::string::npos) {
        sim_type = SimulationType::EDDYCURRENT;
    } else if (setup_type_str.find("Transient") != std::string::npos) {
        sim_type = SimulationType::TRANSIENT;
    } else if (setup_type_str.find("Electrostatic") != std::string::npos) {
        sim_type = SimulationType::ELECTROSTATIC;
    }
    setup->setSolutionType(sim_type);

    if (setup_json.contains("frequency") && (setup_json["frequency"].is_number() || setup_json["frequency"].is_string())) {
        double freq = safeGetDouble(setup_json["frequency"]);
        setup->setFrequency(freq);
    }

    if (setup_json.contains("nonlinear_residual") && (setup_json["nonlinear_residual"].is_number() || setup_json["nonlinear_residual"].is_string())) {
        double residual = safeGetDouble(setup_json["nonlinear_residual"]);
        setup->setConvergenceValue(residual);
    }

    if (setup_json.contains("stop_criterion") && (setup_json["stop_criterion"].is_number() || setup_json["stop_criterion"].is_string())) {
        double error_val = safeGetDouble(setup_json["stop_criterion"]);
        setup->setPercentError(error_val);
    }

    // 新增瞬态求解字段
    setup->setEnabled(setup_json.value("enabled", true));
    setup->setTimeIntegrationMethod(setup_json.value("time_integration_method", 0));
    setup->setSmoothBHCurve(setup_json.value("smooth_bh_curve", false));
    setup->setErrorOutput(setup_json.value("output_error", false));
    setup->setFastReachSteadyState(setup_json.value("fast_reach_steady_state", true));
    setup->setAutoDetectSteadyState(setup_json.value("auto_detect_steady_state", true));
    setup->setStopCriterion(setup_json.value("stop_criterion", 0.005));
    setup->setSaveFieldsType(setup_json.value("save_fields_type", ""));
    setup->setNSteps(setup_json.value("n_steps", ""));
    setup->setStepsFrom(setup_json.value("steps_from", ""));
    setup->setStepsTo(setup_json.value("steps_to", ""));
    setup->setUseAdaptiveTimeStep(setup_json.value("use_adaptive_time_step", false));
    setup->setInitialTimeStep(setup_json.value("initial_time_step", ""));
    setup->setMinTimeStep(setup_json.value("min_time_step", ""));
    setup->setMaxTimeStep(setup_json.value("max_time_step", ""));
    setup->setTimeStepErrTolerance(setup_json.value("time_step_err_tolerance", 0.0001));

    manager.addSolutionSetup(setup);

    FEEM_INFO("  求解器名称: {}", name);
    FEEM_INFO("  求解类型: {}", setup_type_str);
}

void AedtProjectLoader::convertAndAddWindings(ProjectManager& manager,
                                               const std::vector<nlohmann::json>& windings_json) {
    for (const auto& w_json : windings_json) {
        std::string name = w_json.value("name", "Unknown");
        auto winding = std::make_shared<Winding>(name);

        // 设置激励类型
        std::string type_str = w_json.value("type", "Current");
        if (type_str == "Voltage") {
            winding->setExcitationType(WindingExcitationType::VOLTAGE);
        } else {
            winding->setExcitationType(WindingExcitationType::CURRENT);
        }

        // 设置其他字段
        winding->setIsSolid(w_json.value("is_solid", false));
        winding->setCurrentExpression(w_json.value("current", ""));
        winding->setVoltageExpression(w_json.value("voltage", ""));
        winding->setResistance(w_json.value("resistance", ""));
        winding->setInductance(w_json.value("inductance", ""));
        winding->setParallelBranchesNum(w_json.value("parallel_branches_num", ""));

        manager.addWinding(winding);
        FEEM_INFO("  已加载绕组: {} ({})", name, type_str);
    }
    FEEM_INFO("  共加载 {} 个绕组组", windings_json.size());
}

void AedtProjectLoader::convertAndAddMotionSetups(ProjectManager& manager,
                                                    const std::vector<nlohmann::json>& motions_json) {
    for (const auto& m_json : motions_json) {
        std::string name = m_json.value("name", "MotionSetup1");
        auto motion = std::make_shared<MotionSetup>(name);

        // 运动类型映射
        std::string motion_type = m_json.value("motion_type", "Band");
        if (motion_type == "Band") motion->setMotionType(MotionSetupType::BAND);
        else if (motion_type == "Rotation") motion->setMotionType(MotionSetupType::ROTATION);
        else if (motion_type == "Translation") motion->setMotionType(MotionSetupType::TRANSLATION);

        std::string move_type = m_json.value("move_type", "Rotate");
        if (move_type == "Linear") motion->setMoveType(MoveType::LINEAR);
        else motion->setMoveType(MoveType::ROTATE);

        // 轴映射
        std::string axis = m_json.value("axis", "Z");
        if (axis == "X") motion->setAxis(MotionAxis::X);
        else if (axis == "Y") motion->setAxis(MotionAxis::Y);
        else motion->setAxis(MotionAxis::Z);

        motion->setInitialPosition(m_json.value("initial_position", ""));
        motion->setAngularVelocity(m_json.value("angular_velocity", ""));
        motion->setBandNameRef(m_json.value("band_name_ref", -1));

        // 关联对象
        if (m_json.contains("objects") && m_json["objects"].is_array()) {
            for (const auto& obj : m_json["objects"]) {
                motion->addObject(obj.get<int>());
            }
        }
        // Moving子项对象
        if (m_json.contains("moving_objects") && m_json["moving_objects"].is_array()) {
            for (const auto& obj : m_json["moving_objects"]) {
                motion->addMovingObject(obj.get<int>());
            }
        }

        manager.addMotionSetup(motion);
        FEEM_INFO("  已加载运动设置: {} ({}, {})", name, motion_type, move_type);
    }
}

void AedtProjectLoader::convertAndAddMeshOperations(ProjectManager& manager,
                                                      const std::vector<nlohmann::json>& ops_json) {
    for (const auto& op_json : ops_json) {
        std::string name = op_json.value("name", "Unknown");
        auto op = std::make_shared<MeshOperation>(name);

        // 类型映射
        std::string type_str = op_json.value("type", "LengthBased");
        if (type_str == "SurfApproxBased" || type_str == "_skin_depth_based") {
            op->setType(MeshOperationType::SURF_APPROX_BASED);
        } else if (type_str == "CylindricalGap") {
            op->setType(MeshOperationType::CYLINDRICAL_GAP);
        } else {
            op->setType(MeshOperationType::LENGTH_BASED);
        }

        op->setEnabled(op_json.value("enabled", true));

        // 对象列表
        if (op_json.contains("objects") && op_json["objects"].is_array()) {
            for (const auto& obj : op_json["objects"]) {
                op->addObject(obj.get<int>());
            }
        }

        // LengthBased特有字段
        if (op_json.contains("max_length")) op->setMaxLength(op_json["max_length"].get<std::string>());
        if (op_json.contains("num_max_elem")) op->setNumMaxElem(op_json["num_max_elem"].get<std::string>());

        // SurfApproxBased特有字段
        if (op_json.contains("surf_dev")) op->setSurfDev(op_json["surf_dev"].get<std::string>());
        if (op_json.contains("normal_dev")) op->setNormalDev(op_json["normal_dev"].get<std::string>());

        // CylindricalGap特有字段
        if (op_json.contains("band_mapping_angle")) op->setBandMappingAngle(op_json["band_mapping_angle"].get<std::string>());

        manager.addMeshOperation(op);
        FEEM_INFO("  已加载网格操作: {} ({})", name, type_str);
    }
    FEEM_INFO("  共加载 {} 个网格操作", ops_json.size());
}

void AedtProjectLoader::convertAndAddDesignVariables(ProjectManager& manager,
                                                       const std::vector<nlohmann::json>& vars_json) {
    for (const auto& v_json : vars_json) {
        std::string name = v_json.value("name", "Unknown");
        std::string value = v_json.value("value", "");
        std::string unit = v_json.value("unit", "");

        auto var = std::make_shared<DesignVariable>(name, value, unit);
        var->setIndexed(v_json.value("is_indexed", false));
        if (v_json.contains("expression")) {
            var->setExpression(v_json["expression"].get<std::string>());
        }

        manager.addDesignVariable(var);
        FEEM_DEBUG("  已加载设计变量: {} = {}", name, value);
    }
    FEEM_INFO("  共加载 {} 个设计变量", vars_json.size());
}

void AedtProjectLoader::convertAndAddOutputVariables(ProjectManager& manager,
                                                       const std::vector<nlohmann::json>& vars_json) {
    for (const auto& v_json : vars_json) {
        std::string name = v_json.value("name", "Unknown");
        int id = v_json.value("id", 0);
        std::string expr = v_json.value("expression", "");
        std::string result_unit = v_json.value("result_unit", "");
        std::string display_unit = v_json.value("display_unit", "");

        auto var = std::make_shared<OutputVariable>(name, id, expr, result_unit, display_unit);
        manager.addOutputVariable(var);
        FEEM_DEBUG("  已加载输出变量: {} [{}]", name, id);
    }
    FEEM_INFO("  共加载 {} 个输出变量", vars_json.size());
}

void AedtProjectLoader::convertAndAddTemperatureSettings(ProjectManager& manager,
                                                           const nlohmann::json& temp_json) {
    if (temp_json.empty()) return;

    auto settings = std::make_shared<TemperatureSettings>();
    settings->setIncludeTemperatureDependence(temp_json.value("include_temperature_dependence", false));
    settings->setEnableFeedback(temp_json.value("enable_feedback", false));

    if (temp_json.contains("object_temperatures") && temp_json["object_temperatures"].is_object()) {
        for (auto& [key, val] : temp_json["object_temperatures"].items()) {
            try {
                int obj_id = std::stoi(key);
                std::string temp_ref = val.get<std::string>();
                settings->setObjectTemperature(obj_id, temp_ref);
            } catch (...) {}
        }
    }

    manager.setTemperatureSettings(settings);
    FEEM_INFO("  已加载温度设置 ({}个物体映射)",
              settings->getObjectTemperatureMap().size());
}

void AedtProjectLoader::extractPreviewImage(const std::string& file_path) {
    std::string base64_image = parser_.extractPreviewImage();

    if (base64_image.empty()) {
        FEEM_WARN("未找到预览图像");
        return;
    }

    std::vector<uint8_t> image_data = fe_em::tool::Base64Utils::decode(base64_image);

    std::string image_type = fe_em::tool::Base64Utils::detectImageType(image_data);

    fs::path output_dir = fs::path(file_path).parent_path() / "output";
    fs::create_directories(output_dir);

    std::string base_name = fs::path(file_path).stem().string();
    std::string image_saved_path = (output_dir / (base_name + "_preview" + image_type)).string();

    std::ofstream img_file(image_saved_path, std::ios::binary);
    if (img_file.is_open()) {
        img_file.write(reinterpret_cast<const char*>(image_data.data()), image_data.size());
        img_file.close();
        FEEM_INFO("  图像已保存: {} ({} 字节)", image_saved_path, image_data.size());
    } else {
        FEEM_WARN("无法创建图像文件: {}", image_saved_path);
    }
}

} // namespace tool
