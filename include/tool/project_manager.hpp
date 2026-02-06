/**
 * @file project_manager.hpp
 * @brief 基础工具层 - 项目管理器头文件
 * @details 提供项目管理功能，包括创建、打开、保存、关闭项目等
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include "tool/file_utils.hpp"
#include "tool/xml_interface.hpp"
#include "tool/em_exception.hpp"
#include <string>
#include <memory>
#include <unordered_map>
#include <functional>
#include <vector>

namespace tool {

class ProjectManager {
public:
    static ProjectManager& getInstance();

    bool createProject(const std::string& name, const std::string& file_path = "");
    bool openProject(const std::string& file_path);
    bool saveProject(const std::string& file_path = "");
    bool closeProject();
    bool newProject();

    bool isProjectOpen() const { return is_open_; }
    bool isProjectModified() const { return is_modified_; }

    const std::string& getProjectName() const { return project_name_; }
    const std::string& getProjectFilePath() const { return project_file_path_; }
    ProjectState getProjectState() const { return state_; }

    void setProjectName(const std::string& name) { project_name_ = name; is_modified_ = true; }

    void addMaterial(MaterialPtr material);
    bool removeMaterial(const std::string& name);
    std::optional<MaterialPtr> getMaterial(const std::string& name) const;
    const std::unordered_map<std::string, MaterialPtr>& getAllMaterials() const { return materials_; }

    void setGeometry(GeometryPtr geometry);
    std::optional<GeometryPtr> getGeometry() const { return geometry_; }

    void addBoundary(BoundaryPtr boundary);
    bool removeBoundary(const std::string& name);
    std::optional<BoundaryPtr> getBoundary(const std::string& name) const;
    const std::unordered_map<std::string, BoundaryPtr>& getAllBoundaries() const { return boundaries_; }

    void addExcitation(ExcitationPtr excitation);
    bool removeExcitation(const std::string& name);
    std::optional<ExcitationPtr> getExcitation(const std::string& name) const;
    const std::unordered_map<std::string, ExcitationPtr>& getAllExcitations() const { return excitations_; }

    void setMesh(MeshPtr mesh);
    std::optional<MeshPtr> getMesh() const { return mesh_; }

    void addSolutionSetup(SolutionSetupPtr setup);
    bool removeSolutionSetup(const std::string& name);
    std::optional<SolutionSetupPtr> getSolutionSetup(const std::string& name) const;
    const std::unordered_map<std::string, SolutionSetupPtr>& getAllSolutionSetups() const { return solution_setups_; }

    void setDesignType(DimType type) { design_type_ = type; is_modified_ = true; }
    DimType getDesignType() const { return design_type_; }

    void setSolutionType(SimulationType type) { solution_type_ = type; is_modified_ = true; }
    SimulationType getSolutionType() const { return solution_type_; }

    void registerProjectListener(std::function<void(ProjectState, ProjectState)> listener);
    void unregisterProjectListener(std::function<void(ProjectState, ProjectState)> listener);

    MaxwellVersion getFileVersion() const { return file_version_; }
    void setFileVersion(MaxwellVersion version) { file_version_ = version; }

    bool importMaxwellFile(const std::string& file_path);
    bool importAEDTFile(const std::string& file_path);
    bool importEMFFile(const std::string& file_path);

    std::string getLastError() const { return last_error_; }

    static std::string getVersionString();
    static std::string getVersion();

private:
    ProjectManager();
    ~ProjectManager();
    ProjectManager(const ProjectManager&) = delete;
    ProjectManager& operator=(const ProjectManager&) = delete;

    bool loadFromXML(const std::string& file_path);
    bool saveToXML(const std::string& file_path);

    bool validateProject();
    void notifyListeners(ProjectState old_state, ProjectState new_state);

    std::string project_name_;
    std::string project_file_path_;
    ProjectState state_ = ProjectState::CREATED;
    bool is_open_ = false;
    bool is_modified_ = false;

    DimType design_type_ = DimType::D3;
    SimulationType solution_type_ = SimulationType::MAGNETOSTATIC;
    MaxwellVersion file_version_ = MaxwellVersion::UNKNOWN;

    std::unordered_map<std::string, MaterialPtr> materials_;
    std::optional<GeometryPtr> geometry_;
    std::unordered_map<std::string, BoundaryPtr> boundaries_;
    std::unordered_map<std::string, ExcitationPtr> excitations_;
    std::optional<MeshPtr> mesh_;
    std::unordered_map<std::string, SolutionSetupPtr> solution_setups_;

    std::vector<std::function<void(ProjectState, ProjectState)>> listeners_;
    std::string last_error_;

    static constexpr const char* CURRENT_VERSION = "1.0";
};

bool createProjectFromTemplate(const std::string& template_path, const std::string& output_path,
                               const std::unordered_map<std::string, std::string>& replacements);

} // namespace tool
