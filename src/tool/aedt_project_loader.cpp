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
    FEEM_INFO("[5/6] 解析求解设置...");
    auto setup_json = parser_.parseSolutionSetup();
    convertAndAddSolutionSetup(manager, setup_json);

    FEEM_INFO("");
    FEEM_INFO("[6/6] 提取预览图像...");
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

    manager.addSolutionSetup(setup);

    FEEM_INFO("  求解器名称: {}", name);
    FEEM_INFO("  求解类型: {}", setup_type_str);
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
