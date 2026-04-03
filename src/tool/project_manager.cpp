/**
 * @file project_manager.cpp
 * @brief 基础工具层 - 项目管理器源文件
 */

#include "tool/project_manager.hpp"
#include "tool/project_loader.hpp"
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
    loader_.reset();

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

    loader_ = createProjectLoader(file_path);
    if (!loader_) {
        last_error_ = "Unsupported file format: " + file_utils::getExtension(file_path);
        FEEM_ERROR(last_error_);
        return false;
    }

    loader_->setFilePath(file_path);
    bool success = loader_->load(*this);

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

    if (ext == ".emat" || ext == ".EMAT") {
        last_error_ = "EMAT save not supported yet (requires XmlProjectLoader)";
        FEEM_ERROR(last_error_);
        return false;
    } else {
        last_error_ = "Unsupported file format for saving: " + ext;
        FEEM_ERROR(last_error_);
        return false;
    }
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
