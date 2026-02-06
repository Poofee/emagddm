/**
 * @file project_data.hpp
 * @brief 基础工具层 - 项目数据结构头文件
 * @details 定义项目管理所需的全部数据结构
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/em_enums.hpp"
#include "tool/id_generator.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <optional>
#include <variant>
#include <functional>

namespace tool {

class ProjectData;

using MaterialPtr = std::shared_ptr<class Material>;
using GeometryPtr = std::shared_ptr<class Geometry>;
using BoundaryPtr = std::shared_ptr<class Boundary>;
using ExcitationPtr = std::shared_ptr<class Excitation>;
using MeshPtr = std::shared_ptr<class Mesh>;
using ResultPtr = std::shared_ptr<class Result>;
using SolutionSetupPtr = std::shared_ptr<class SolutionSetup>;

struct MaterialProperty {
    std::string name;
    std::string type;
    std::variant<double, std::vector<double>, std::string> value;
    std::string unit;
    bool is_temperature_dependent = false;
    std::vector<std::pair<double, std::vector<double>>> temp_dependent_data;
};

struct BHDataPoint {
    double h; 
    double b;
};

class Material {
public:
    Material(const std::string& name);
    ~Material() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }
    MatType getType() const { return type_; }
    void setType(MatType type) { type_ = type; }

    void addProperty(const MaterialProperty& property);
    std::optional<MaterialProperty> getProperty(const std::string& name) const;
    const std::vector<MaterialProperty>& getAllProperties() const { return properties_; }

    void setCoreLossEnabled(bool enabled) { core_loss_enabled_ = enabled; }
    bool isCoreLossEnabled() const { return core_loss_enabled_; }

    void setBHCurve(const std::vector<BHDataPoint>& curve);
    const std::vector<BHDataPoint>& getBHCurve() const { return bh_curve_; }

    void setBHCurveType(BHCurveType type) { bh_curve_type_ = type; }
    BHCurveType getBHCurveType() const { return bh_curve_type_; }

    void setCoreLossModel(CoreLossModelType model) { core_loss_model_ = model; }
    CoreLossModelType getCoreLossModel() const { return core_loss_model_; }

    void setCoreLossCoefficients(double ks, double alpha, double beta, double kn);
    void getCoreLossCoefficients(double& ks, double& alpha, double& beta, double& kn) const;

    void setRelativePermeability(double mu_r);
    double getRelativePermeability() const { return relative_permeability_; }

    void setConductivity(double sigma);
    double getConductivity() const { return conductivity_; }

    void setMassDensity(double rho);
    double getMassDensity() const { return mass_density_; }

private:
    std::string name_;
    uint64_t id_;
    MatType type_ = MatType::LINEAR_ISOTROPIC;
    std::vector<MaterialProperty> properties_;
    std::unordered_map<std::string, int> property_index_;

    bool core_loss_enabled_ = false;
    CoreLossModelType core_loss_model_ = CoreLossModelType::NONE;
    BHCurveType bh_curve_type_ = BHCurveType::NONE;
    std::vector<BHDataPoint> bh_curve_;

    double relative_permeability_ = 1.0;
    double conductivity_ = 0.0;
    double mass_density_ = 0.0;
    double core_loss_ks_ = 0.0;
    double core_loss_alpha_ = 0.0;
    double core_loss_beta_ = 0.0;
    double core_loss_kn_ = 0.0;
};

class Geometry {
public:
    Geometry(const std::string& name);
    ~Geometry() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }
    DimType getDimension() const { return dimension_; }
    void setDimension(DimType dim) { dimension_ = dim; }

    void addSubGeometry(const std::string& name, GeometryPtr sub_geo);
    std::optional<GeometryPtr> getSubGeometry(const std::string& name) const;
    const std::unordered_map<std::string, GeometryPtr>& getAllSubGeometries() const { return sub_geometries_; }

    void addObject(const std::string& name, const std::string& material_name);
    void assignMaterial(const std::string& object_name, const std::string& material_name);

    const std::unordered_map<std::string, std::string>& getObjectMaterialMap() const { return object_material_map_; }

    void setFilePath(const std::string& path) { file_path_ = path; }
    const std::string& getFilePath() const { return file_path_; }

    void setGeometryType(const std::string& type) { geometry_type_ = type; }
    const std::string& getGeometryType() const { return geometry_type_; }

private:
    std::string name_;
    uint64_t id_;
    DimType dimension_ = DimType::D3;
    std::string file_path_;
    std::string geometry_type_;
    std::unordered_map<std::string, GeometryPtr> sub_geometries_;
    std::unordered_map<std::string, std::string> object_material_map_;
};

class Boundary {
public:
    Boundary(const std::string& name);
    ~Boundary() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }
    BndType getType() const { return type_; }
    void setType(BndType type) { type_ = type; }

    void addFace(const std::string& face_id);
    void addEdge(const std::string& edge_id);
    void addObject(const std::string& object_name);

    const std::vector<std::string>& getFaces() const { return faces_; }
    const std::vector<std::string>& getEdges() const { return edges_; }
    const std::vector<std::string>& getObjects() const { return objects_; }

    void setImpedanceValue(double z) { impedance_value_ = z; }
    double getImpedanceValue() const { return impedance_value_; }

    void setVectorPotential(double az) { vector_potential_ = az; }
    double getVectorPotential() const { return vector_potential_; }

    void setVoltage(double v) { voltage_ = v; }
    double getVoltage() const { return voltage_; }

    void setCurrent(double i) { current_ = i; }
    double getCurrent() const { return current_; }

    void setMasterName(const std::string& name) { master_name_ = name; }
    const std::string& getMasterName() const { return master_name_; }

    void setSlaveName(const std::string& name) { slave_name_ = name; }
    const std::string& getSlaveName() const { return slave_name_; }

private:
    std::string name_;
    uint64_t id_;
    BndType type_ = BndType::DIRICHLET;
    std::vector<std::string> faces_;
    std::vector<std::string> edges_;
    std::vector<std::string> objects_;

    double impedance_value_ = 0.0;
    double vector_potential_ = 0.0;
    double voltage_ = 0.0;
    double current_ = 0.0;
    std::string master_name_;
    std::string slave_name_;
};

class Excitation {
public:
    Excitation(const std::string& name);
    ~Excitation() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }
    ExcitationType getType() const { return type_; }
    void setType(ExcitationType type) { type_ = type; }

    void setValue(double value) { value_ = value; }
    double getValue() const { return value_; }

    void setPhase(double phase) { phase_ = phase; }
    double getPhase() const { return phase_; }

    void setIsSolid(bool is_solid) { is_solid_ = is_solid; }
    bool isSolid() const { return is_solid_; }

    void setCoilGroup(const std::string& group) { coil_group_ = group; }
    const std::string& getCoilGroup() const { return coil_group_; }

    void setConnectionType(CoilConnectionType conn_type) { connection_type_ = conn_type; }
    CoilConnectionType getConnectionType() const { return connection_type_; }

    void setNumberOfTurns(int turns) { number_of_turns_ = turns; }
    int getNumberOfTurns() const { return number_of_turns_; }

    void setPolygonPoints(const std::vector<std::pair<double, double>>& points) { polygon_points_ = points; }
    const std::vector<std::pair<double, double>>& getPolygonPoints() const { return polygon_points_; }

    void setDirection(int direction) { direction_ = direction; }
    int getDirection() const { return direction_; }

private:
    std::string name_;
    uint64_t id_;
    ExcitationType type_ = ExcitationType::CURRENT_DENSITY;
    double value_ = 0.0;
    double phase_ = 0.0;
    bool is_solid_ = false;
    std::string coil_group_;
    CoilConnectionType connection_type_ = CoilConnectionType::SERIES;
    int number_of_turns_ = 1;
    std::vector<std::pair<double, double>> polygon_points_;
    int direction_ = 1;
};

class Mesh {
public:
    Mesh(const std::string& name);
    ~Mesh() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }

    void setGenerationType(MeshGenerationType type) { generation_type_ = type; }
    MeshGenerationType getGenerationType() const { return generation_type_; }

    void setMaxElementSize(double size) { max_element_size_ = size; }
    double getMaxElementSize() const { return max_element_size_; }

    void setMinElementSize(double size) { min_element_size_ = size; }
    double getMinElementSize() const { return min_element_size_; }

    void setSurfaceApproximation(double approx) { surface_approximation_ = approx; }
    double getSurfaceApproximation() const { return surface_approximation_; }

    void enableCurvatureRefinement(bool enable) { curvature_refinement_ = enable; }
    bool isCurvatureRefinementEnabled() const { return curvature_refinement_; }

    void enableBoundaryLayer(bool enable) { boundary_layer_ = enable; }
    bool isBoundaryLayerEnabled() const { return boundary_layer_; }

    void setBoundaryLayerNumberOfLayers(int layers) { boundary_layer_num_ = layers; }
    int getBoundaryLayerNumberOfLayers() const { return boundary_layer_num_; }

    void setBoundaryLayerThickness(double thickness) { boundary_layer_thickness_ = thickness; }
    double getBoundaryLayerThickness() const { return boundary_layer_thickness_; }

    void setAdaptiveDepth(int depth) { adaptive_depth_ = depth; }
    int getAdaptiveDepth() const { return adaptive_depth_; }

    void addObjectMeshSettings(const std::string& object_name, double min, double max);
    const std::unordered_map<std::string, std::pair<double, double>>& getObjectMeshSettings() const { return object_mesh_settings_; }

private:
    std::string name_;
    uint64_t id_;
    MeshGenerationType generation_type_ = MeshGenerationType::AUTOMATIC;
    double max_element_size_ = 0.0;
    double min_element_size_ = 0.0;
    double surface_approximation_ = 0.0;
    bool curvature_refinement_ = false;
    bool boundary_layer_ = false;
    int boundary_layer_num_ = 0;
    double boundary_layer_thickness_ = 0.0;
    int adaptive_depth_ = 0;
    std::unordered_map<std::string, std::pair<double, double>> object_mesh_settings_;
};

class SolutionSetup {
public:
    SolutionSetup(const std::string& name);
    ~SolutionSetup() = default;

    const std::string& getName() const { return name_; }
    uint64_t getID() const { return id_; }

    void setSolutionType(SimulationType type) { solution_type_ = type; }
    SimulationType getSolutionType() const { return solution_type_; }

    void setSolverType(SolverType type) { solver_type_ = type; }
    SolverType getSolverType() const { return solver_type_; }

    void setConvergenceType(ConvergenceType type) { convergence_type_ = type; }
    ConvergenceType getConvergenceType() const { return convergence_type_; }

    void setConvergenceValue(double value) { convergence_value_ = value; }
    double getConvergenceValue() const { return convergence_value_; }

    void setMaximumIterations(int max_iter) { maximum_iterations_ = max_iter; }
    int getMaximumIterations() const { return maximum_iterations_; }

    void setFrequency(double freq) { frequency_ = freq; }
    double getFrequency() const { return frequency_; }

    void setStatorFrequency(double freq) { stator_frequency_ = freq; }
    double getStatorFrequency() const { return stator_frequency_; }

    void setMeshRefinementEnabled(bool enable) { mesh_refinement_ = enable; }
    bool isMeshRefinementEnabled() const { return mesh_refinement_; }

    void setMeshRefinementPercent(double percent) { mesh_refinement_percent_ = percent; }
    double getMeshRefinementPercent() const { return mesh_refinement_percent_; }

    void setPercentError(double error) { percent_error_ = error; }
    double getPercentError() const { return percent_error_; }

private:
    std::string name_;
    uint64_t id_;
    SimulationType solution_type_ = SimulationType::MAGNETOSTATIC;
    SolverType solver_type_ = SolverType::AUTO;
    ConvergenceType convergence_type_ = ConvergenceType::RESIDUAL;
    double convergence_value_ = 0.001;
    int maximum_iterations_ = 100;
    double frequency_ = 0.0;
    double stator_frequency_ = 0.0;
    bool mesh_refinement_ = false;
    double mesh_refinement_percent_ = 0.0;
    double percent_error_ = 1.0;
};

} // namespace tool
