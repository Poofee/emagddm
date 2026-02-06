/**
 * @file maxwell_reader.cpp
 * @brief Maxwell数据读取模块 - 核心实现
 * @details 实现Maxwell文件解析、数据转换和验证的核心功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/maxwell_reader.hpp"
#include "tool/logger.hpp"
#include <fstream>
#include <sstream>

namespace tool {

// MaxwellReader类实现

MaxwellReader::MaxwellReader(std::unique_ptr<IMaxwellParser> parser,
                             std::unique_ptr<IMaxwellConverter> converter,
                             std::unique_ptr<IMaxwellValidator> validator)
    : parser_(std::move(parser))
    , converter_(std::move(converter))
    , validator_(std::move(validator))
    , solution_setup_(nullptr) {
    
    // 初始化日志
    auto logger = LoggerFactory::getLogger("MaxwellReader");
    logger->info("MaxwellReader initialized");
}

void MaxwellReader::setFilePath(const std::string& file_path) {
    file_path_ = file_path;
    
    auto logger = LoggerFactory::getLogger("MaxwellReader");
    logger->info("Set file path: {}", file_path);
}

bool MaxwellReader::readProject() {
    auto logger = LoggerFactory::getLogger("MaxwellReader");
    
    try {
        // 1. 检查文件路径
        if (file_path_.empty()) {
            logger->error("File path is empty");
            return false;
        }

        // 2. 检查文件是否存在
        if (!file_utils::fileExists(file_path_)) {
            logger->error("File not found: {}", file_path_);
            throw MaxwellException(MaxwellErrorCode::FILE_NOT_FOUND, 
                                  "File not found: " + file_path_);
        }

        // 3. 检查解析器是否支持该文件
        if (!parser_->canParse(file_path_)) {
            logger->error("File format not supported: {}", file_path_);
            throw MaxwellException(MaxwellErrorCode::INVALID_FORMAT, 
                                  "File format not supported: " + file_path_);
        }

        logger->info("Starting to read Maxwell project: {}", file_path_);

        // 4. 解析文件信息
        logger->debug("Parsing file info");
        file_info_ = parser_->parseFileInfo();
        
        // 5. 解析材料数据
        logger->debug("Parsing materials");
        auto materials_data = parser_->parseMaterials();
        materials_ = converter_->convertMaterials(materials_data);
        logger->info("Parsed {} materials", materials_.size());

        // 6. 解析边界条件数据
        logger->debug("Parsing boundaries");
        auto boundaries_data = parser_->parseBoundaries();
        boundaries_ = converter_->convertBoundaries(boundaries_data);
        logger->info("Parsed {} boundaries", boundaries_.size());

        // 7. 解析激励源数据
        logger->debug("Parsing excitations");
        auto excitations_data = parser_->parseExcitations();
        excitations_ = converter_->convertExcitations(excitations_data);
        logger->info("Parsed {} excitations", excitations_.size());

        // 8. 解析求解设置数据
        logger->debug("Parsing solution setup");
        auto setup_data = parser_->parseSolutionSetup();
        solution_setup_ = converter_->convertSolutionSetup(setup_data);
        logger->info("Parsed solution setup");

        // 9. 验证数据
        logger->debug("Validating project data");
        validation_result_ = validator_->validateProjectData(
            materials_, boundaries_, excitations_, solution_setup_);
        validation_result_.generateSummary();
        
        if (!validation_result_.is_valid) {
            logger->warn("Project validation failed: {}", validation_result_.summary);
            for (const auto& error : validation_result_.errors) {
                logger->error("Validation error: {}", error);
            }
        } else {
            logger->info("Project validation passed: {}", validation_result_.summary);
        }

        logger->info("Successfully read Maxwell project: {}", file_path_);
        return true;

    } catch (const MaxwellException& e) {
        logger->error("MaxwellException during reading: {} (code: {})", 
                     e.what(), static_cast<int>(e.getErrorCode()));
        return false;
    } catch (const std::exception& e) {
        logger->error("Exception during reading: {}", e.what());
        return false;
    }
}

std::vector<std::shared_ptr<Material>> MaxwellReader::getMaterials() const {
    return materials_;
}

std::vector<std::shared_ptr<Boundary>> MaxwellReader::getBoundaries() const {
    return boundaries_;
}

std::vector<std::shared_ptr<Excitation>> MaxwellReader::getExcitations() const {
    return excitations_;
}

std::shared_ptr<SolutionSetup> MaxwellReader::getSolutionSetup() const {
    return solution_setup_;
}

ValidationResult MaxwellReader::getValidationResult() const {
    return validation_result_;
}

MaxwellFileInfo MaxwellReader::getFileInfo() const {
    return file_info_;
}

nlohmann::json MaxwellReader::exportToJson() const {
    nlohmann::json result;

    try {
        // 导出文件信息
        result["file_info"] = file_info_.toJson();

        // 导出材料数据
        nlohmann::json materials_json;
        for (const auto& material : materials_) {
            materials_json.push_back(material->toJson());
        }
        result["materials"] = materials_json;

        // 导出边界条件数据
        nlohmann::json boundaries_json;
        for (const auto& boundary : boundaries_) {
            boundaries_json.push_back(boundary->toJson());
        }
        result["boundaries"] = boundaries_json;

        // 导出激励源数据
        nlohmann::json excitations_json;
        for (const auto& excitation : excitations_) {
            excitations_json.push_back(excitation->toJson());
        }
        result["excitations"] = excitations_json;

        // 导出求解设置数据
        if (solution_setup_) {
            result["solution_setup"] = solution_setup_->toJson();
        }

        // 导出验证结果
        result["validation_result"] = validation_result_.toJson();

    } catch (const std::exception& e) {
        auto logger = LoggerFactory::getLogger("MaxwellReader");
        logger->error("Error during JSON export: {}", e.what());
    }

    return result;
}

bool MaxwellReader::saveToProject(const std::string& output_path) const {
    auto logger = LoggerFactory::getLogger("MaxwellReader");
    
    try {
        // 1. 验证数据
        if (!validation_result_.is_valid) {
            logger->warn("Cannot save invalid project data");
            return false;
        }

        // 2. 导出为JSON
        auto project_json = exportToJson();

        // 3. 保存到文件
        std::ofstream file(output_path);
        if (!file.is_open()) {
            logger->error("Failed to open output file: {}", output_path);
            return false;
        }

        file << project_json.dump(4); // 缩进4个空格，便于阅读
        file.close();

        logger->info("Project saved to: {}", output_path);
        return true;

    } catch (const std::exception& e) {
        logger->error("Error during project save: {}", e.what());
        return false;
    }
}

// 基础解析器实现（示例）

class BasicMaxwellParser : public IMaxwellParser {
private:
    std::string file_path_;
    nlohmann::json project_data_;

public:
    BasicMaxwellParser() = default;

    bool canParse(const std::string& file_path) override {
        // 检查文件扩展名
        std::string ext = file_utils::getFileExtension(file_path);
        return ext == ".aedt" || ext == ".aedtz" || ext == ".xml" || ext == ".amat";
    }

    MaxwellFileInfo parseFileInfo() override {
        MaxwellFileInfo info;
        
        // 解析文件基本信息
        info.file_path = file_path_;
        info.file_format = file_utils::getFileExtension(file_path_);
        
        // 这里应该解析Maxwell文件的实际内容
        // 暂时返回示例数据
        info.maxwell_version = "R2024.1";
        info.simulation_type = "Magnetostatic";
        info.dimension = "2D";
        info.project_name = "Example Project";
        info.created_date = "2026-02-06T10:00:00Z";
        info.modified_date = "2026-02-06T10:00:00Z";
        info.file_size = 1024; // 示例大小
        info.is_compressed = false;

        return info;
    }

    std::vector<nlohmann::json> parseMaterials() override {
        std::vector<nlohmann::json> materials;
        
        // 这里应该解析Maxwell文件中的材料数据
        // 暂时返回示例数据
        nlohmann::json material1 = {
            {"name", "Copper"},
            {"type", "Conductor"},
            {"conductivity", 5.8e7},
            {"relative_permeability", 1.0},
            {"maxwell_specific_params", {
                {"material_id", "mat_copper_001"},
                {"library_name", "Built-in"}
            }}
        };

        nlohmann::json material2 = {
            {"name", "Air"},
            {"type", "NonMagnetic"},
            {"conductivity", 0.0},
            {"relative_permeability", 1.0},
            {"maxwell_specific_params", {
                {"material_id", "mat_air_001"},
                {"library_name", "Built-in"}
            }}
        };

        materials.push_back(material1);
        materials.push_back(material2);

        return materials;
    }

    std::vector<nlohmann::json> parseBoundaries() override {
        std::vector<nlohmann::json> boundaries;
        
        // 这里应该解析Maxwell文件中的边界条件数据
        // 暂时返回示例数据
        nlohmann::json boundary1 = {
            {"name", "Outer Boundary"},
            {"type", "Balloon"},
            {"impedance_value", 0.0},
            {"voltage", 0.0},
            {"current", 0.0},
            {"maxwell_specific_params", {
                {"boundary_id", "bnd_outer_001"},
                {"geometry_reference", "outer_edge"}
            }}
        };

        boundaries.push_back(boundary1);

        return boundaries;
    }

    std::vector<nlohmann::json> parseExcitations() override {
        std::vector<nlohmann::json> excitations;
        
        // 这里应该解析Maxwell文件中的激励源数据
        // 暂时返回示例数据
        nlohmann::json excitation1 = {
            {"name", "Coil Excitation"},
            {"type", "Current"},
            {"magnitude", 10.0},
            {"phase", 0.0},
            {"frequency", 50.0},
            {"maxwell_specific_params", {
                {"excitation_id", "exc_coil_001"},
                {"winding_type", "Solid"}
            }}
        };

        excitations.push_back(excitation1);

        return excitations;
    }

    nlohmann::json parseSolutionSetup() override {
        // 这里应该解析Maxwell文件中的求解设置数据
        // 暂时返回示例数据
        nlohmann::json setup = {
            {"name", "Setup1"},
            {"solution_type", "Magnetostatic"},
            {"solver_type", "Iterative"},
            {"convergence_type", "Relative"},
            {"convergence_value", 0.001},
            {"maximum_iterations", 100},
            {"frequency", 0.0},
            {"stator_frequency", 0.0},
            {"maxwell_specific_params", {
                {"setup_id", "setup_001"},
                {"adaptive_meshing", true}
            }}
        };

        return setup;
    }

    nlohmann::json parseGeometry() override {
        // 这里应该解析Maxwell文件中的几何数据
        // 暂时返回示例数据
        nlohmann::json geometry = {
            {"dimension", "2D"},
            {"units", "mm"},
            {"maxwell_specific_params", {
                {"geometry_id", "geom_001"},
                {"coordinate_system", "Cartesian"}
            }}
        };

        return geometry;
    }

    nlohmann::json parseAllData() override {
        // 解析所有数据并返回完整JSON
        nlohmann::json all_data;
        
        all_data["file_info"] = parseFileInfo().toJson();
        
        // 材料数据
        auto materials = parseMaterials();
        all_data["materials"] = materials;
        
        // 边界条件数据
        auto boundaries = parseBoundaries();
        all_data["boundaries"] = boundaries;
        
        // 激励源数据
        auto excitations = parseExcitations();
        all_data["excitations"] = excitations;
        
        // 求解设置数据
        all_data["solution_setup"] = parseSolutionSetup();
        
        // 几何数据
        all_data["geometry"] = parseGeometry();

        return all_data;
    }
};

// 基础转换器实现

class BasicMaxwellConverter : public IMaxwellConverter {
public:
    BasicMaxwellConverter() = default;

    std::shared_ptr<Material> convertMaterial(const nlohmann::json& material_data) override {
        try {
            auto material = std::make_shared<Material>(
                material_data.value("name", "Unknown Material")
            );
            
            // 设置材料类型
            std::string type_str = material_data.value("type", "");
            if (type_str == "Conductor") {
                material->setType(MatType::CONDUCTOR);
            } else if (type_str == "NonMagnetic") {
                material->setType(MatType::NON_MAGNETIC);
            } else {
                material->setType(MatType::LINEAR_MAGNETIC);
            }
            
            // 设置材料参数
            material->setConductivity(material_data.value("conductivity", 0.0));
            material->setRelativePermeability(material_data.value("relative_permeability", 1.0));
            
            // 设置Maxwell特定参数
            auto maxwell_params = material_data.value("maxwell_specific_params", nlohmann::json::object());
            for (auto& [key, value] : maxwell_params.items()) {
                material->setMaxwellSpecificParam(key, value);
            }
            
            return material;
            
        } catch (const std::exception& e) {
            auto logger = LoggerFactory::getLogger("MaxwellConverter");
            logger->error("Error converting material: {}", e.what());
            return nullptr;
        }
    }

    std::shared_ptr<Boundary> convertBoundary(const nlohmann::json& boundary_data) override {
        try {
            auto boundary = std::make_shared<Boundary>(
                boundary_data.value("name", "Unknown Boundary")
            );
            
            // 设置边界类型
            std::string type_str = boundary_data.value("type", "");
            if (type_str == "Balloon") {
                boundary->setType(BndType::BALLOON);
            } else if (type_str == "Symmetry") {
                boundary->setType(BndType::SYMMETRY);
            } else {
                boundary->setType(BndType::DIRICHLET);
            }
            
            // 设置边界参数
            boundary->setImpedanceValue(boundary_data.value("impedance_value", 0.0));
            boundary->setVoltage(boundary_data.value("voltage", 0.0));
            boundary->setCurrent(boundary_data.value("current", 0.0));
            
            // 设置Maxwell特定参数
            auto maxwell_params = boundary_data.value("maxwell_specific_params", nlohmann::json::object());
            for (auto& [key, value] : maxwell_params.items()) {
                boundary->setMaxwellSpecificParam(key, value);
            }
            
            return boundary;
            
        } catch (const std::exception& e) {
            auto logger = LoggerFactory::getLogger("MaxwellConverter");
            logger->error("Error converting boundary: {}", e.what());
            return nullptr;
        }
    }

    std::shared_ptr<Excitation> convertExcitation(const nlohmann::json& excitation_data) override {
        try {
            auto excitation = std::make_shared<Excitation>(
                excitation_data.value("name", "Unknown Excitation")
            );
            
            // 设置激励类型
            std::string type_str = excitation_data.value("type", "");
            if (type_str == "Current") {
                excitation->setType(ExcType::CURRENT);
            } else if (type_str == "Voltage") {
                excitation->setType(ExcType::VOLTAGE);
            } else {
                excitation->setType(ExcType::CURRENT);
            }
            
            // 设置激励参数
            excitation->setMagnitude(excitation_data.value("magnitude", 0.0));
            excitation->setPhase(excitation_data.value("phase", 0.0));
            excitation->setFrequency(excitation_data.value("frequency", 0.0));
            
            // 设置Maxwell特定参数
            auto maxwell_params = excitation_data.value("maxwell_specific_params", nlohmann::json::object());
            for (auto& [key, value] : maxwell_params.items()) {
                excitation->setMaxwellSpecificParam(key, value);
            }
            
            return excitation;
            
        } catch (const std::exception& e) {
            auto logger = LoggerFactory::getLogger("MaxwellConverter");
            logger->error("Error converting excitation: {}", e.what());
            return nullptr;
        }
    }

    std::shared_ptr<SolutionSetup> convertSolutionSetup(const nlohmann::json& setup_data) override {
        try {
            auto setup = std::make_shared<SolutionSetup>(
                setup_data.value("name", "Unknown Setup")
            );
            
            // 设置求解类型
            std::string solution_type_str = setup_data.value("solution_type", "");
            if (solution_type_str == "Magnetostatic") {
                setup->setSolutionType(SimulationType::MAGNETOSTATIC);
            } else {
                setup->setSolutionType(SimulationType::MAGNETOSTATIC);
            }
            
            // 设置求解器类型
            std::string solver_type_str = setup_data.value("solver_type", "");
            if (solver_type_str == "Iterative") {
                setup->setSolverType(SolverType::ITERATIVE);
            } else {
                setup->setSolverType(SolverType::DIRECT);
            }
            
            // 设置收敛类型
            std::string convergence_type_str = setup_data.value("convergence_type", "");
            if (convergence_type_str == "Relative") {
                setup->setConvergenceType(ConvergenceType::RELATIVE);
            } else {
                setup->setConvergenceType(ConvergenceType::ABSOLUTE);
            }
            
            // 设置求解参数
            setup->setConvergenceValue(setup_data.value("convergence_value", 0.001));
            setup->setMaximumIterations(setup_data.value("maximum_iterations", 100));
            setup->setFrequency(setup_data.value("frequency", 0.0));
            setup->setStatorFrequency(setup_data.value("stator_frequency", 0.0));
            
            // 设置Maxwell特定参数
            auto maxwell_params = setup_data.value("maxwell_specific_params", nlohmann::json::object());
            for (auto& [key, value] : maxwell_params.items()) {
                setup->setMaxwellSpecificParam(key, value);
            }
            
            return setup;
            
        } catch (const std::exception& e) {
            auto logger = LoggerFactory::getLogger("MaxwellConverter");
            logger->error("Error converting solution setup: {}", e.what());
            return nullptr;
        }
    }

    std::shared_ptr<void> convertGeometry(const nlohmann::json& geometry_data) override {
        // 几何数据转换（预留接口）
        // 返回空指针，表示功能尚未实现
        return nullptr;
    }

    std::vector<std::shared_ptr<Material>> convertMaterials(const std::vector<nlohmann::json>& materials_data) override {
        std::vector<std::shared_ptr<Material>> materials;
        
        for (const auto& material_data : materials_data) {
            auto material = convertMaterial(material_data);
            if (material) {
                materials.push_back(material);
            }
        }
        
        return materials;
    }

    std::vector<std::shared_ptr<Boundary>> convertBoundaries(const std::vector<nlohmann::json>& boundaries_data) override {
        std::vector<std::shared_ptr<Boundary>> boundaries;
        
        for (const auto& boundary_data : boundaries_data) {
            auto boundary = convertBoundary(boundary_data);
            if (boundary) {
                boundaries.push_back(boundary);
            }
        }
        
        return boundaries;
    }

    std::vector<std::shared_ptr<Excitation>> convertExcitations(const std::vector<nlohmann::json>& excitations_data) override {
        std::vector<std::shared_ptr<Excitation>> excitations;
        
        for (const auto& excitation_data : excitations_data) {
            auto excitation = convertExcitation(excitation_data);
            if (excitation) {
                excitations.push_back(excitation);
            }
        }
        
        return excitations;
    }
};

// 基础验证器实现

class BasicMaxwellValidator : public IMaxwellValidator {
public:
    BasicMaxwellValidator() = default;

    ValidationResult validateMaterials(const std::vector<std::shared_ptr<Material>>& materials) override {
        ValidationResult result;
        
        for (const auto& material : materials) {
            if (!material) {
                result.addError("Null material pointer found");
                continue;
            }
            
            if (!material->validate()) {
                result.addError("Material validation failed: " + material->getName());
            }
            
            // 检查材料参数范围
            if (material->getConductivity() < 0.0) {
                result.addWarning("Material conductivity is negative: " + material->getName());
            }
            
            if (material->getRelativePermeability() <= 0.0) {
                result.addError("Material relative permeability must be positive: " + material->getName());
            }
        }
        
        return result;
    }

    ValidationResult validateBoundaries(const std::vector<std::shared_ptr<Boundary>>& boundaries) override {
        ValidationResult result;
        
        for (const auto& boundary : boundaries) {
            if (!boundary) {
                result.addError("Null boundary pointer found");
                continue;
            }
            
            if (!boundary->validate()) {
                result.addError("Boundary validation failed: " + boundary->getName());
            }
            
            // 检查边界参数范围
            if (boundary->getImpedanceValue() < 0.0) {
                result.addWarning("Boundary impedance is negative: " + boundary->getName());
            }
            
            if (boundary->getVoltage() < 0.0) {
                result.addWarning("Boundary voltage is negative: " + boundary->getName());
            }
        }
        
        return result;
    }

    ValidationResult validateExcitations(const std::vector<std::shared_ptr<Excitation>>& excitations) override {
        ValidationResult result;
        
        for (const auto& excitation : excitations) {
            if (!excitation) {
                result.addError("Null excitation pointer found");
                continue;
            }
            
            if (!excitation->validate()) {
                result.addError("Excitation validation failed: " + excitation->getName());
            }
            
            // 检查激励参数范围
            if (excitation->getMagnitude() < 0.0) {
                result.addWarning("Excitation magnitude is negative: " + excitation->getName());
            }
            
            if (excitation->getFrequency() < 0.0) {
                result.addError("Excitation frequency must be non-negative: " + excitation->getName());
            }
        }
        
        return result;
    }

    ValidationResult validateSolutionSetup(const std::shared_ptr<SolutionSetup>& setup) override {
        ValidationResult result;
        
        if (!setup) {
            result.addError("Null solution setup pointer found");
            return result;
        }
        
        if (!setup->validate()) {
            result.addError("Solution setup validation failed: " + setup->getName());
        }
        
        // 检查求解参数范围
        if (setup->getConvergenceValue() <= 0.0) {
            result.addError("Convergence value must be positive: " + setup->getName());
        }
        
        if (setup->getMaximumIterations() <= 0) {
            result.addError("Maximum iterations must be positive: " + setup->getName());
        }
        
        if (setup->getFrequency() < 0.0) {
            result.addError("Frequency must be non-negative: " + setup->getName());
        }
        
        return result;
    }

    ValidationResult validateGeometry(const std::shared_ptr<void>& geometry) override {
        // 几何数据验证（预留接口）
        ValidationResult result;
        result.addWarning("Geometry validation not implemented yet");
        return result;
    }

    ValidationResult validateProjectData(
        const std::vector<std::shared_ptr<Material>>& materials,
        const std::vector<std::shared_ptr<Boundary>>& boundaries,
        const std::vector<std::shared_ptr<Excitation>>& excitations,
        const std::shared_ptr<SolutionSetup>& setup) override {
        
        ValidationResult result;
        
        // 验证各个组件
        auto materials_result = validateMaterials(materials);
        auto boundaries_result = validateBoundaries(boundaries);
        auto excitations_result = validateExcitations(excitations);
        auto setup_result = validateSolutionSetup(setup);
        
        // 合并结果
        result.is_valid = materials_result.is_valid && boundaries_result.is_valid && 
                         excitations_result.is_valid && setup_result.is_valid;
        
        result.warnings.insert(result.warnings.end(), materials_result.warnings.begin(), materials_result.warnings.end());
        result.warnings.insert(result.warnings.end(), boundaries_result.warnings.begin(), boundaries_result.warnings.end());
        result.warnings.insert(result.warnings.end(), excitations_result.warnings.begin(), excitations_result.warnings.end());
        result.warnings.insert(result.warnings.end(), setup_result.warnings.begin(), setup_result.warnings.end());
        
        result.errors.insert(result.errors.end(), materials_result.errors.begin(), materials_result.errors.end());
        result.errors.insert(result.errors.end(), boundaries_result.errors.begin(), boundaries_result.errors.end());
        result.errors.insert(result.errors.end(), excitations_result.errors.begin(), excitations_result.errors.end());
        result.errors.insert(result.errors.end(), setup_result.errors.begin(), setup_result.errors.end());
        
        // 检查数据一致性
        if (materials.empty()) {
            result.addWarning("No materials defined in project");
        }
        
        if (boundaries.empty()) {
            result.addWarning("No boundaries defined in project");
        }
        
        if (excitations.empty()) {
            result.addWarning("No excitations defined in project");
        }
        
        if (!setup) {
            result.addError("No solution setup defined in project");
        }
        
        return result;
    }
};

} // namespace tool