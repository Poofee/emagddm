/**
 * @file project_manager.cpp
 * @brief 基础工具层 - 项目管理器源文件
 */

#include "tool/project_manager.hpp"
#include "tool/logger_factory.hpp"
#include <fstream>
#include <sstream>

namespace tool {

ProjectManager::ProjectManager()
    : project_name_("Untitled"),
      project_file_path_(""),
      state_(ProjectState::CREATED),
      is_open_(false),
      is_modified_(false),
      design_type_(DimType::D3),
      solution_type_(SimulationType::MAGNETOSTATIC),
      file_version_(MaxwellVersion::UNKNOWN) {
}

ProjectManager::~ProjectManager() {
    if (is_open_) {
        closeProject();
    }
}

ProjectManager& ProjectManager::getInstance() {
    static ProjectManager instance;
    return instance;
}

bool ProjectManager::createProject(const std::string& name, const std::string& file_path) {
    if (is_open_) {
        if (!closeProject()) {
            last_error_ = "Failed to close current project";
            return false;
        }
    }

    project_name_ = name;
    project_file_path_ = file_path;
    state_ = ProjectState::CREATED;
    is_open_ = true;
    is_modified_ = false;

    materials_.clear();
    geometry_.reset();
    boundaries_.clear();
    excitations_.clear();
    mesh_.reset();
    solution_setups_.clear();

    notifyListeners(ProjectState::CLOSING, state_);
    FEEM_INFO("Project '{}' created successfully", name);
    return true;
}

bool ProjectManager::newProject() {
    return createProject("Untitled", "");
}

bool ProjectManager::openProject(const std::string& file_path) {
    if (!file_utils::exists(file_path)) {
        last_error_ = "Project file not found: " + file_path;
        FEEM_ERROR(last_error_);
        return false;
    }

    if (is_open_) {
        if (!closeProject()) {
            last_error_ = "Failed to close current project";
            return false;
        }
    }

    std::string ext = file_utils::getExtension(file_path);
    bool success = false;

    if (ext == ".emat" || ext == ".EMAT") {
        success = loadFromXML(file_path);
    } else if (ext == ".aedt" || ext == ".AEDT") {
        success = importAEDTFile(file_path);
    } else if (ext == ".emf" || ext == ".EMF") {
        success = importEMFFile(file_path);
    } else {
        last_error_ = "Unsupported project file format: " + ext;
        FEEM_ERROR(last_error_);
        return false;
    }

    if (success) {
        is_open_ = true;
        is_modified_ = false;
        state_ = ProjectState::LOADED;
        notifyListeners(ProjectState::CLOSING, state_);
        FEEM_INFO("Project loaded from '{}'", file_path);
    }

    return success;
}

bool ProjectManager::saveProject(const std::string& file_path) {
    if (!is_open_) {
        last_error_ = "No project is currently open";
        FEEM_ERROR(last_error_);
        return false;
    }

    std::string target_path = file_path.empty() ? project_file_path_ : file_path;
    if (target_path.empty()) {
        last_error_ = "No file path specified for saving";
        FEEM_ERROR(last_error_);
        return false;
    }

    if (!validateProject()) {
        FEEM_ERROR("Project validation failed: {}", last_error_);
        return false;
    }

    std::string ext = file_utils::getExtension(target_path);
    bool success = false;

    if (ext == ".emat" || ext == ".EMAT") {
        success = saveToXML(target_path);
    } else {
        last_error_ = "Unsupported file format for saving: " + ext;
        FEEM_ERROR(last_error_);
        return false;
    }

    if (success) {
        project_file_path_ = target_path;
        is_modified_ = false;
        state_ = ProjectState::SAVED;
        notifyListeners(ProjectState::MODIFIED, state_);
        FEEM_INFO("Project saved to '{}'", target_path);
    }

    return success;
}

bool ProjectManager::closeProject() {
    if (!is_open_) {
        return true;
    }

    if (is_modified_) {
        FEEM_WARN("Project '{}' has unsaved changes", project_name_);
    }

    state_ = ProjectState::CLOSING;
    notifyListeners(state_, ProjectState::CLOSING);

    materials_.clear();
    geometry_.reset();
    boundaries_.clear();
    excitations_.clear();
    mesh_.reset();
    solution_setups_.clear();

    project_name_ = "Untitled";
    project_file_path_ = "";
    is_open_ = false;
    is_modified_ = false;
    state_ = ProjectState::CLOSING;

    FEEM_INFO("Project closed");
    return true;
}

bool ProjectManager::loadFromXML(const std::string& file_path) {
    XmlDocument doc;
    if (!doc.loadFromFile(file_path)) {
        last_error_ = "Failed to load XML file: " + file_path;
        return false;
    }

    XmlNode root = doc.getRootNode();
    if (!root.isValid()) {
        last_error_ = "Invalid XML document structure";
        return false;
    }

    std::string root_name = root.getName();
    if (root_name != "MaxwellProject" && root_name != "Project") {
        last_error_ = "Invalid project file format";
        return false;
    }

    auto version_attr = root.getAttribute("Version");
    if (version_attr.has_value()) {
        file_version_ = stringToMaxwellVersion(version_attr.value());
    }

    auto name_node = root.getChild("Name");
    if (name_node.isValid()) {
        project_name_ = name_node.getText();
    }

    auto design_node = root.getChild("Design");
    if (design_node.isValid()) {
        auto dim_node = design_node.getChild("Dimension");
        if (dim_node.isValid()) {
            std::string dim_str = dim_node.getText();
            try {
                design_type_ = stringToDimType(dim_str);
            } catch (...) {
                design_type_ = DimType::D3;
            }
        }
    }

    auto materials_node = root.getChild("Materials");
    if (materials_node.isValid()) {
        auto material_nodes = materials_node.getChildren("Material");
        for (const auto& mat_node : material_nodes) {
            std::string mat_name = mat_node.getName();
            auto name_attr = mat_node.getAttribute("Name");
            if (name_attr.has_value()) {
                mat_name = name_attr.value();
            }

            auto material = std::make_shared<Material>(mat_name);

            auto type_node = mat_node.getChild("Type");
            if (type_node.isValid()) {
                std::string type_str = type_node.getText();
                try {
                    material->setType(stringToMatType(type_str));
                } catch (...) {}
            }

            auto mu_node = mat_node.getChild("RelativePermeability");
            if (mu_node.isValid()) {
                try {
                    double mu_r = std::stod(mu_node.getText());
                    material->setRelativePermeability(mu_r);
                } catch (...) {}
            }

            auto sigma_node = mat_node.getChild("Conductivity");
            if (sigma_node.isValid()) {
                try {
                    double sigma = std::stod(sigma_node.getText());
                    material->setConductivity(sigma);
                } catch (...) {}
            }

            materials_[mat_name] = material;
        }
    }

    return true;
}

bool ProjectManager::saveToXML(const std::string& file_path) {
    XmlDocument doc;
    doc.createNew("MaxwellProject");
    XmlNode root = doc.getRootNode();
    root.setAttribute("Version", getVersion());

    XmlNode name_node = root.appendChild("Name");
    name_node.setText(project_name_);

    XmlNode design_node = root.appendChild("Design");
    XmlNode dim_node = design_node.appendChild("Dimension");
    dim_node.setText(dimTypeToString(design_type_));

    XmlNode materials_node = root.appendChild("Materials");
    for (const auto& pair : materials_) {
        XmlNode mat_node = materials_node.appendChild("Material");
        mat_node.setAttribute("Name", pair.first);
        mat_node.setAttribute("ID", std::to_string(pair.second->getID()));

        XmlNode type_node = mat_node.appendChild("Type");
        type_node.setText(matTypeToString(pair.second->getType()));

        XmlNode mu_node = mat_node.appendChild("RelativePermeability");
        mu_node.setText(std::to_string(pair.second->getRelativePermeability()));

        XmlNode sigma_node = mat_node.appendChild("Conductivity");
        sigma_node.setText(std::to_string(pair.second->getConductivity()));
    }

    if (!doc.saveToFile(file_path, true)) {
        last_error_ = "Failed to save XML file: " + file_path;
        return false;
    }

    return true;
}

bool ProjectManager::validateProject() {
    if (materials_.empty()) {
        last_error_ = "Project has no materials defined";
        return false;
    }

    if (!geometry_.has_value()) {
        last_error_ = "Project has no geometry defined";
        return false;
    }

    return true;
}

void ProjectManager::notifyListeners(ProjectState old_state, ProjectState new_state) {
    for (const auto& listener : listeners_) {
        try {
            listener(old_state, new_state);
        } catch (...) {
            FEEM_ERROR("Exception in project listener");
        }
    }
}

void ProjectManager::addMaterial(MaterialPtr material) {
    materials_[material->getName()] = material;
    is_modified_ = true;
}

bool ProjectManager::removeMaterial(const std::string& name) {
    auto it = materials_.find(name);
    if (it != materials_.end()) {
        materials_.erase(it);
        is_modified_ = true;
        return true;
    }
    return false;
}

std::optional<MaterialPtr> ProjectManager::getMaterial(const std::string& name) const {
    auto it = materials_.find(name);
    if (it != materials_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void ProjectManager::setGeometry(GeometryPtr geometry) {
    geometry_ = geometry;
    is_modified_ = true;
}

void ProjectManager::addBoundary(BoundaryPtr boundary) {
    boundaries_[boundary->getName()] = boundary;
    is_modified_ = true;
}

bool ProjectManager::removeBoundary(const std::string& name) {
    auto it = boundaries_.find(name);
    if (it != boundaries_.end()) {
        boundaries_.erase(it);
        is_modified_ = true;
        return true;
    }
    return false;
}

std::optional<BoundaryPtr> ProjectManager::getBoundary(const std::string& name) const {
    auto it = boundaries_.find(name);
    if (it != boundaries_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void ProjectManager::addExcitation(ExcitationPtr excitation) {
    excitations_[excitation->getName()] = excitation;
    is_modified_ = true;
}

bool ProjectManager::removeExcitation(const std::string& name) {
    auto it = excitations_.find(name);
    if (it != excitations_.end()) {
        excitations_.erase(it);
        is_modified_ = true;
        return true;
    }
    return false;
}

std::optional<ExcitationPtr> ProjectManager::getExcitation(const std::string& name) const {
    auto it = excitations_.find(name);
    if (it != excitations_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void ProjectManager::setMesh(MeshPtr mesh) {
    mesh_ = mesh;
    is_modified_ = true;
}

void ProjectManager::addSolutionSetup(SolutionSetupPtr setup) {
    solution_setups_[setup->getName()] = setup;
    is_modified_ = true;
}

bool ProjectManager::removeSolutionSetup(const std::string& name) {
    auto it = solution_setups_.find(name);
    if (it != solution_setups_.end()) {
        solution_setups_.erase(it);
        is_modified_ = true;
        return true;
    }
    return false;
}

std::optional<SolutionSetupPtr> ProjectManager::getSolutionSetup(const std::string& name) const {
    auto it = solution_setups_.find(name);
    if (it != solution_setups_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void ProjectManager::registerProjectListener(std::function<void(ProjectState, ProjectState)> listener) {
    listeners_.push_back(listener);
}

void ProjectManager::unregisterProjectListener(std::function<void(ProjectState, ProjectState)> listener) {
    listeners_.erase(
        std::remove_if(listeners_.begin(), listeners_.end(),
            [&listener](const auto& l) { return l.target_type() == listener.target_type(); }),
        listeners_.end()
    );
}

bool ProjectManager::importMaxwellFile(const std::string& file_path) {
    std::string ext = file_utils::getExtension(file_path);
    if (ext == ".aedt" || ext == ".AEDT") {
        return importAEDTFile(file_path);
    } else if (ext == ".emf" || ext == ".EMF") {
        return importEMFFile(file_path);
    }
    last_error_ = "Unsupported Maxwell file format: " + ext;
    return false;
}

bool ProjectManager::importAEDTFile(const std::string& file_path) {
    FEEM_INFO("AEDT file import not fully implemented yet");
    last_error_ = "AEDT import requires additional implementation";
    return false;
}

bool ProjectManager::importEMFFile(const std::string& file_path) {
    FEEM_INFO("EMF file import not fully implemented yet");
    last_error_ = "EMF import requires additional implementation";
    return false;
}

std::string ProjectManager::getVersionString() {
    return "EMAGDDM Project Version " + std::string(CURRENT_VERSION);
}

std::string ProjectManager::getVersion() {
    return CURRENT_VERSION;
}

bool createProjectFromTemplate(const std::string& template_path, const std::string& output_path,
                              const std::unordered_map<std::string, std::string>& replacements) {
    if (!file_utils::exists(template_path)) {
        FEEM_ERROR("Template file not found: {}", template_path);
        return false;
    }

    std::string content;
    if (!file_utils::readText(template_path, content)) {
        FEEM_ERROR("Failed to read template file: {}", template_path);
        return false;
    }

    for (const auto& pair : replacements) {
        size_t pos = 0;
        while ((pos = content.find(pair.first, pos)) != std::string::npos) {
            content.replace(pos, pair.first.length(), pair.second);
            pos += pair.second.length();
        }
    }

    if (!file_utils::writeText(output_path, content)) {
        FEEM_ERROR("Failed to write project file: {}", output_path);
        return false;
    }

    FEEM_INFO("Project created from template: {}", output_path);
    return true;
}

} // namespace tool
