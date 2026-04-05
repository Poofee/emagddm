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
#include "tool/project_loader.hpp"
#include "tool/em_exception.hpp"
#include "fe_em/em_mesh_data.hpp"
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

    // ========== 电磁网格拓扑数据管理（方案C架构） ==========
    /**
     * @brief 设置电磁网格拓扑数据
     * @param data 网格拓扑数据（转移所有权）
     * 
     * @note 架构说明：
     *       - mesh_ (tool::Mesh): 存储**网格剖分配置参数**（max_size, boundary_layer等）
     *       - em_mesh_data_ (fe_em::EMMeshData): 存储**网格拓扑数据**（节点坐标、单元连接关系）
     *       - 两者职责分离，统一由ProjectManager管理
     */
    void setEMMeshData(std::unique_ptr<fe_em::EMMeshData> data);
    
    /**
     * @brief 获取电磁网格拓扑数据
     * @return fe_em::EMMeshData* 网格拓扑数据指针（不转移所有权，可为nullptr）
     */
    fe_em::EMMeshData* getEMMeshData() const;
    
    /**
     * @brief 检查是否有网格拓扑数据
     * @return bool 如果存在返回true
     */
    bool hasEMMeshData() const;

    void addSolutionSetup(SolutionSetupPtr setup);
    bool removeSolutionSetup(const std::string& name);
    std::optional<SolutionSetupPtr> getSolutionSetup(const std::string& name) const;
    const std::unordered_map<std::string, SolutionSetupPtr>& getAllSolutionSetups() const { return solution_setups_; }

    void setDesignType(DimType type) { design_type_ = type; is_modified_ = true; }
    DimType getDesignType() const { return design_type_; }

    void setSolutionType(SimulationType type) { solution_type_ = type; is_modified_ = true; }
    SimulationType getSolutionType() const { return solution_type_; }

    // ========== 绕组组管理 ==========
    void addWinding(WindingPtr winding);
    const std::unordered_map<std::string, WindingPtr>& getWindings() const { return windings_; }
    size_t getWindingCount() const { return windings_.size(); }
    bool hasWinding(const std::string& name) const;
    WindingPtr getWinding(const std::string& name) const;

    // ========== 运动设置管理 ==========
    void addMotionSetup(MotionSetupPtr setup);
    const std::vector<MotionSetupPtr>& getMotionSetups() const { return motion_setups_; }
    size_t getMotionSetupCount() const { return motion_setups_.size(); }

    // ========== 网格操作管理 ==========
    void addMeshOperation(MeshOperationPtr operation);
    const std::unordered_map<std::string, MeshOperationPtr>& getMeshOperations() const { return mesh_operations_; }
    size_t getMeshOperationCount() const { return mesh_operations_.size(); }
    bool hasMeshOperation(const std::string& name) const;
    MeshOperationPtr getMeshOperation(const std::string& name) const;

    // ========== 设计变量管理 ==========
    void addDesignVariable(DesignVariablePtr variable);
    const std::unordered_map<std::string, DesignVariablePtr>& getDesignVariables() const { return design_variables_; }
    size_t getDesignVariableCount() const { return design_variables_.size(); }

    // ========== 输出变量管理 ==========
    void addOutputVariable(OutputVariablePtr variable);
    const std::unordered_map<std::string, OutputVariablePtr>& getOutputVariables() const { return output_variables_; }
    size_t getOutputVariableCount() const { return output_variables_.size(); }

    // ========== 温度设置管理 ==========
    void setTemperatureSettings(TemperatureSettingsPtr settings);
    TemperatureSettingsPtr getTemperatureSettings() const;
    bool hasTemperatureSettings() const { return temperature_settings_.has_value(); }

    void registerProjectListener(std::function<void(ProjectState, ProjectState)> listener);
    void unregisterProjectListener(std::function<void(ProjectState, ProjectState)> listener);

    MaxwellVersion getFileVersion() const { return file_version_; }
    void setFileVersion(MaxwellVersion version) { file_version_ = version; }

    std::string getLastError() const { return last_error_; }

    static std::string getVersionString();
    static std::string getVersion();

private:
    ProjectManager();
    ~ProjectManager();
    ProjectManager(const ProjectManager&) = delete;
    ProjectManager& operator=(const ProjectManager&) = delete;

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
    std::optional<MeshPtr> mesh_;                    ///< 网格剖分配置参数（tool::Mesh对象）
    std::unique_ptr<fe_em::EMMeshData> em_mesh_data_; ///< 网格拓扑数据（fe_em::EMMeshData对象，方案C架构）
    std::unordered_map<std::string, SolutionSetupPtr> solution_setups_;

    // 绕组组集合 (Winding Group)
    std::unordered_map<std::string, WindingPtr> windings_;

    // 运动设置列表 (通常只有1个)
    std::vector<MotionSetupPtr> motion_setups_;

    // 网格操作集合
    std::unordered_map<std::string, MeshOperationPtr> mesh_operations_;

    // 设计变量集合
    std::unordered_map<std::string, DesignVariablePtr> design_variables_;

    // 输出变量集合
    std::unordered_map<std::string, OutputVariablePtr> output_variables_;

    // 温度设置 (可选)
    std::optional<TemperatureSettingsPtr> temperature_settings_;

    std::vector<std::function<void(ProjectState, ProjectState)>> listeners_;
    std::string last_error_;
    std::unique_ptr<IProjectLoader> loader_;       ///< 内部加载器（格式相关，openProject时自动创建）

    static constexpr const char* CURRENT_VERSION = "1.0";
};

} // namespace tool
