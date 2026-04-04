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
#include <map>
#include <memory>
#include <optional>
#include <variant>
#include <functional>
#include "../../lib/nlohmann/json.hpp"

namespace tool {

/**
 * @brief 序列化接口
 * @details 提供JSON和二进制序列化功能
 */
class ISerializable {
public:
    virtual ~ISerializable() = default;

    /**
     * @brief 序列化为JSON对象
     * @return JSON对象
     */
    virtual nlohmann::json toJson() const = 0;

    /**
     * @brief 从JSON对象反序列化
     * @param json JSON对象
     * @return 是否成功
     */
    virtual bool fromJson(const nlohmann::json& json) = 0;

    /**
     * @brief 序列化为二进制数据
     * @param data 二进制数据输出
     * @return 是否成功
     */
    virtual bool toBinary(std::vector<uint8_t>& data) const = 0;

    /**
     * @brief 从二进制数据反序列化
     * @param data 二进制数据
     * @param offset 数据偏移量
     * @return 是否成功
     */
    virtual bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) = 0;

    /**
     * @brief 获取序列化版本
     * @return 版本号
     */
    virtual uint32_t getSerializationVersion() const = 0;

    /**
     * @brief 验证数据完整性
     * @return 是否有效
     */
    virtual bool validate() const = 0;
};

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

class Material : public ISerializable {
public:
    // 嵌套结构体定义（供公共接口使用，必须在方法声明之前）
    struct CoreLossCoefficientSetup {
        std::string mode;
        double frequency = 0;
        std::string thickness;
        std::string conductivity;
        std::vector<std::pair<double, double>> bp_curve;
    };
    struct ThermalModifier {
        std::string property_name;
        int index = 0;
        std::string formula_string;
    };
    struct MaterialAppearance {
        int red = 200, green = 200, blue = 200;
        double transparency = 0.0;
    };

    Material(const std::string& name);
    Material() = default; // 默认构造函数用于反序列化
    ~Material() override = default;

    // ISerializable接口实现
    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }
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

    // Maxwell专属数据接口
    void setMaxwellMaterialID(const std::string& id) { maxwell_material_id_ = id; }
    const std::string& getMaxwellMaterialID() const { return maxwell_material_id_; }

    void setTemperatureCoefficient(double tc) { temperature_coefficient_ = tc; }
    double getTemperatureCoefficient() const { return temperature_coefficient_; }

    void setBHCustomCurveFile(const std::string& file_path) { bh_custom_curve_file_ = file_path; }
    const std::string& getBHCustomCurveFile() const { return bh_custom_curve_file_; }

    void setCoreLossUserDataFile(const std::string& file_path) { coreloss_user_data_file_ = file_path; }
    const std::string& getCoreLossUserDataFile() const { return coreloss_user_data_file_; }

    void setAnisotropicPermeability(const std::vector<double>& permeability);
    const std::vector<double>& getAnisotropicPermeability() const { return anisotropic_permeability_; }

    void setAnisotropicConductivity(const std::vector<double>& conductivity);
    const std::vector<double>& getAnisotropicConductivity() const { return anisotropic_conductivity_; }

    void setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params);
    const std::unordered_map<std::string, std::string>& getMaxwellSpecificParameters() const { return maxwell_specific_params_; }

    // 材料元信息接口
    void setLibraryInfo(const std::string& lib, const std::string& location, bool mod_since);
    const std::string& getLibrary() const { return library_; }
    const std::string& getLibLocation() const { return lib_location_; }
    bool isModSinceLib() const { return mod_since_lib_; }
    int64_t getModTime() const { return mod_time_; }
    void setModTime(int64_t time) { mod_time_ = time; }
    void setPhysicsTypes(const std::vector<std::string>& types) { physics_types_ = types; }
    const std::vector<std::string>& getPhysicsTypes() const { return physics_types_; }
    void setCoordinateSystemType(const std::string& type) { coordinate_system_type_ = type; }
    const std::string& getCoordinateSystemType() const { return coordinate_system_type_; }
    void setBulkOrSurfaceType(int type) { bulk_or_surface_type_ = type; }
    int getBulkOrSurfaceType() const { return bulk_or_surface_type_; }

    // 扩展物理属性接口
    void setPermittivity(double perm) { permittivity_ = perm; }
    double getPermittivity() const { return permittivity_; }
    void setYoungsModulus(double e) { youngs_modulus_ = e; }
    double getYoungsModulus() const { return youngs_modulus_; }
    void setPoisonsRatio(double v) { poisons_ratio_ = v; }
    double getPoisonsRatio() const { return poisons_ratio_; }
    void setThermalExpansionCoefficient(double alpha) { thermal_expansion_coefficient_ = alpha; }
    double getThermalExpansionCoefficient() const { return thermal_expansion_coefficient_; }

    // B-H曲线增强信息接口
    void setBHCurveDataType(const std::string& type) { bh_curve_data_type_ = type; }
    const std::string& getBHCurveDataType() const { return bh_curve_data_type_; }
    void setHUnit(const std::string& unit) { h_unit_ = unit; }
    const std::string& getHUnit() const { return h_unit_; }
    void setBUnit(const std::string& unit) { b_unit_ = unit; }
    const std::string& getBUnit() const { return b_unit_; }
    void setIsTemperatureDependent(bool dep) { is_temperature_dependent_ = dep; }
    bool isTemperatureDependent() const { return is_temperature_dependent_; }
    void setTemperatureBHCurves(const std::map<double, std::vector<BHDataPoint>>& curves);
    const std::map<double, std::vector<BHDataPoint>>& getTemperatureBHCurves() const { return temperature_bh_curves_; }

    // 铁损增强接口
    void setCoreLossTypeChoice(const std::string& choice) { core_loss_type_choice_ = choice; }
    const std::string& getCoreLossTypeChoice() const { return core_loss_type_choice_; }
    void setStackingTypeStr(const std::string& type) { stacking_type_str_ = type; }
    const std::string& getStackingTypeStr() const { return stacking_type_str_; }
    void setCoreLossFreqCurves(const std::map<double, std::vector<std::pair<double, double>>>& curves);
    const std::map<double, std::vector<std::pair<double, double>>>& getCoreLossFreqCurves() const { return core_loss_freq_curves_; }
    void setCoreLossUnit(const std::string& unit) { core_loss_unit_ = unit; }
    const std::string& getCoreLossUnit() const { return core_loss_unit_; }
    void setCoreLossCoeffSetup(const CoreLossCoefficientSetup& setup);
    const std::optional<CoreLossCoefficientSetup>& getCoreLossCoeffSetup() const { return core_loss_coeff_setup_; }

    // 热修正器接口
    void addThermalModifier(const ThermalModifier& modifier);
    const std::vector<ThermalModifier>& getThermalModifiers() const { return thermal_modifiers_; }

    // 外观数据接口
    void setAppearance(int r, int g, int b, double transparency);
    const std::optional<MaterialAppearance>& getAppearance() const { return appearance_data_; }

    // 备注和关键词接口
    void setNotes(const std::string& notes) { notes_ = notes; }
    const std::string& getNotes() const { return notes_; }
    void setKeywords(const std::string& keywords) { keywords_ = keywords; }
    const std::string& getKeywords() const { return keywords_; }

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

    // Maxwell专属数据成员
    std::string maxwell_material_id_;
    double temperature_coefficient_ = 0.0;
    std::string bh_custom_curve_file_;
    std::string coreloss_user_data_file_;
    std::vector<double> anisotropic_permeability_;
    std::vector<double> anisotropic_conductivity_;
    std::unordered_map<std::string, std::string> maxwell_specific_params_;

    // 材料元信息
    std::string library_;
    std::string lib_location_;
    bool mod_since_lib_ = false;
    int64_t mod_time_ = 0;
    std::vector<std::string> physics_types_;
    std::string coordinate_system_type_ = "Cartesian";
    int bulk_or_surface_type_ = 1;

    // 扩展物理属性
    double permittivity_ = 1.0;
    double youngs_modulus_ = 0.0;
    double poisons_ratio_ = 0.0;
    double thermal_expansion_coefficient_ = 0.0;

    // B-H曲线增强信息
    std::string bh_curve_data_type_ = "normal";
    std::string h_unit_ = "A_per_meter";
    std::string b_unit_ = "tesla";
    bool is_temperature_dependent_ = false;
    std::map<double, std::vector<BHDataPoint>> temperature_bh_curves_;

    // 铁损增强
    std::string core_loss_type_choice_;
    std::string stacking_type_str_ = "Solid";
    std::map<double, std::vector<std::pair<double, double>>> core_loss_freq_curves_;
    std::string core_loss_unit_ = "w_per_kg";
    std::optional<CoreLossCoefficientSetup> core_loss_coeff_setup_;

    // 热修正器
    std::vector<ThermalModifier> thermal_modifiers_;

    // 外观数据
    std::optional<MaterialAppearance> appearance_data_;

    // 备注和关键词
    std::string notes_;
    std::string keywords_;
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

class Boundary : public ISerializable {
public:
    Boundary(const std::string& name);
    Boundary() = default; // 默认构造函数用于反序列化
    ~Boundary() override = default;

    // ISerializable接口实现
    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }
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

    // Maxwell专属边界数据接口
    void setBoundarySubType(BoundarySubType sub_type) { boundary_sub_type_ = sub_type; }
    BoundarySubType getBoundarySubType() const { return boundary_sub_type_; }

    void setPeriodicMappingType(PeriodicMappingType mapping_type) { periodic_mapping_type_ = mapping_type; }
    PeriodicMappingType getPeriodicMappingType() const { return periodic_mapping_type_; }

    void setRadiationDistance(double distance) { radiation_distance_ = distance; }
    double getRadiationDistance() const { return radiation_distance_; }

    void setPerfectESymmetry(bool symmetry) { perfect_e_symmetry_ = symmetry; }
    bool getPerfectESymmetry() const { return perfect_e_symmetry_; }

    void setPerfectHSymmetry(bool symmetry) { perfect_h_symmetry_ = symmetry; }
    bool getPerfectHSymmetry() const { return perfect_h_symmetry_; }

    void setInfiniteSphereRadius(double radius) { infinite_sphere_radius_ = radius; }
    double getInfiniteSphereRadius() const { return infinite_sphere_radius_; }

    void setMaxwellBoundaryID(const std::string& id) { maxwell_boundary_id_ = id; }
    const std::string& getMaxwellBoundaryID() const { return maxwell_boundary_id_; }

    void setBoundarySubdivisionParameters(const std::vector<double>& params);
    const std::vector<double>& getBoundarySubdivisionParameters() const { return boundary_subdivision_params_; }

    void setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params);
    const std::unordered_map<std::string, std::string>& getMaxwellSpecificParameters() const { return maxwell_specific_params_; }

    // 边界扩展接口
    void setParentBndID(int id) { parent_bnd_id_ = id; }
    int getParentBndID() const { return parent_bnd_id_; }
    void setIsComponent(bool is_comp) { is_component_ = is_comp; }
    bool isComponent() const { return is_component_; }
    void setCoordinateSystem(int coord_sys) { coordinate_system_ = coord_sys; }
    int getCoordinateSystem() const { return coordinate_system_; }
    void setConductorNumberStr(const std::string& num) { conductor_number_str_ = num; }
    const std::string& getConductorNumberStr() const { return conductor_number_str_; }

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

    // Maxwell专属边界数据成员
    BoundarySubType boundary_sub_type_ = BoundarySubType::NONE;
    PeriodicMappingType periodic_mapping_type_ = PeriodicMappingType::NONE;
    double radiation_distance_ = 0.0;
    bool perfect_e_symmetry_ = false;
    bool perfect_h_symmetry_ = false;
    double infinite_sphere_radius_ = 0.0;
    std::string maxwell_boundary_id_;
    std::vector<double> boundary_subdivision_params_;
    std::unordered_map<std::string, std::string> maxwell_specific_params_;

    // 边界扩展成员
    int parent_bnd_id_ = -1;
    bool is_component_ = false;
    int coordinate_system_ = -1;
    std::string conductor_number_str_;
};

class Excitation : public ISerializable {
public:
    Excitation(const std::string& name);
    Excitation() = default; // 默认构造函数用于反序列化
    ~Excitation() override = default;

    // ISerializable接口实现
    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }
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

    // Maxwell专属激励数据接口
    void setWaveformType(ExcitationWaveformType waveform) { waveform_type_ = waveform; }
    ExcitationWaveformType getWaveformType() const { return waveform_type_; }

    void setFrequency(double freq) { frequency_ = freq; }
    double getFrequency() const { return frequency_; }

    void setDutyCycle(double duty) { duty_cycle_ = duty; }
    double getDutyCycle() const { return duty_cycle_; }

    void setWindingType(WindingType winding) { winding_type_ = winding; }
    WindingType getWindingType() const { return winding_type_; }

    void setMotionType(MotionType motion) { motion_type_ = motion; }
    MotionType getMotionType() const { return motion_type_; }

    void setRotationSpeed(double speed) { rotation_speed_ = speed; }
    double getRotationSpeed() const { return rotation_speed_; }

    void setTranslationSpeed(double speed) { translation_speed_ = speed; }
    double getTranslationSpeed() const { return translation_speed_; }

    void setExternalCircuitFile(const std::string& file_path) { external_circuit_file_ = file_path; }
    const std::string& getExternalCircuitFile() const { return external_circuit_file_; }

    void setCustomWaveformFile(const std::string& file_path) { custom_waveform_file_ = file_path; }
    const std::string& getCustomWaveformFile() const { return custom_waveform_file_; }

    void setMaxwellExcitationID(const std::string& id) { maxwell_excitation_id_ = id; }
    const std::string& getMaxwellExcitationID() const { return maxwell_excitation_id_; }

    void setWaveformParameters(const std::vector<double>& params);
    const std::vector<double>& getWaveformParameters() const { return waveform_params_; }

    void setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params);
    const std::unordered_map<std::string, std::string>& getMaxwellSpecificParameters() const { return maxwell_specific_params_; }

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

    // Maxwell专属激励数据成员
    ExcitationWaveformType waveform_type_ = ExcitationWaveformType::DC;
    double frequency_ = 0.0;
    double duty_cycle_ = 0.5;
    WindingType winding_type_ = WindingType::SOLID;
    MotionType motion_type_ = MotionType::NONE;
    double rotation_speed_ = 0.0;
    double translation_speed_ = 0.0;
    std::string external_circuit_file_;
    std::string custom_waveform_file_;
    std::string maxwell_excitation_id_;
    std::vector<double> waveform_params_;
    std::unordered_map<std::string, std::string> maxwell_specific_params_;
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

    // 网格扩展接口
    void setSurfApproxChoice(const std::string& choice) { surf_approx_choice_ = choice; }
    const std::string& getSurfApproxChoice() const { return surf_approx_choice_; }
    void setSliderMeshSettings(int settings) { slider_mesh_settings_ = settings; }
    int getSliderMeshSettings() const { return slider_mesh_settings_; }
    void setUseAutoLength(bool use) { use_auto_length_ = use; }
    bool isUseAutoLength() const { return use_auto_length_; }

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

    // 网格扩展字段
    std::string surf_approx_choice_ = "UseSlider";
    int slider_mesh_settings_ = 5;
    bool use_auto_length_ = true;
};

class SolutionSetup : public ISerializable {
public:
    SolutionSetup(const std::string& name);
    SolutionSetup() = default; // 默认构造函数用于反序列化
    ~SolutionSetup() override = default;

    // ISerializable接口实现
    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

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

    // Maxwell专属求解设置接口
    void setHPCParallelMode(HPCParallelMode mode) { hpc_parallel_mode_ = mode; }
    HPCParallelMode getHPCParallelMode() const { return hpc_parallel_mode_; }

    void setHPCSolverMode(HPCSolverMode mode) { hpc_solver_mode_ = mode; }
    HPCSolverMode getHPCSolverMode() const { return hpc_solver_mode_; }

    void setNumCores(int cores) { num_cores_ = cores; }
    int getNumCores() const { return num_cores_; }

    void setDomainDecompositionType(DomainDecompositionType type) { domain_decomposition_type_ = type; }
    DomainDecompositionType getDomainDecompositionType() const { return domain_decomposition_type_; }

    void setAdaptiveMeshRefinement(bool enable) { adaptive_mesh_refinement_ = enable; }
    bool isAdaptiveMeshRefinementEnabled() const { return adaptive_mesh_refinement_; }

    void setAdaptiveDepth(int depth) { adaptive_depth_ = depth; }
    int getAdaptiveDepth() const { return adaptive_depth_; }

    void setSkinDepthRefinement(bool enable) { skin_depth_refinement_ = enable; }
    bool isSkinDepthRefinementEnabled() const { return skin_depth_refinement_; }

    void setCoreLossRefinement(bool enable) { coreloss_refinement_ = enable; }
    bool isCoreLossRefinementEnabled() const { return coreloss_refinement_; }

    void setMaxwellSolverID(const std::string& id) { maxwell_solver_id_ = id; }
    const std::string& getMaxwellSolverID() const { return maxwell_solver_id_; }

    void setHPCParameters(const std::unordered_map<std::string, std::string>& params);
    const std::unordered_map<std::string, std::string>& getHPCParameters() const { return hpc_params_; }

    void setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params);
    const std::unordered_map<std::string, std::string>& getMaxwellSpecificParameters() const { return maxwell_specific_params_; }

    // 瞬态求解扩展接口
    void setEnabled(bool val) { enabled_ = val; }
    bool isEnabled() const { return enabled_; }
    void setMeshLinkImport(bool import) { mesh_link_import_ = import; }
    bool isMeshLinkImport() const { return mesh_link_import_; }
    void setTimeIntegrationMethod(int method) { time_integration_method_ = method; }
    int getTimeIntegrationMethod() const { return time_integration_method_; }
    void setSmoothBHCurve(bool smooth) { smooth_bh_curve_ = smooth; }
    bool isSmoothBHCurve() const { return smooth_bh_curve_; }
    void setErrorOutput(bool output) { output_error_ = output; }
    bool isErrorOutput() const { return output_error_; }
    void setOutputPerObjectCoreLoss(bool output) { output_per_object_core_loss_ = output; }
    bool isOutputPerObjectCoreLoss() const { return output_per_object_core_loss_; }
    void setOutputPerObjectSolidLoss(bool output) { output_per_object_solid_loss_ = output; }
    bool isOutputPerObjectSolidLoss() const { return output_per_object_solid_loss_; }
    void setUseControlProgram(bool use) { use_control_program_ = use; }
    bool isUseControlProgram() const { return use_control_program_; }
    void setControlProgramName(const std::string& name) { control_program_name_ = name; }
    const std::string& getControlProgramName() const { return control_program_name_; }
    void setControlProgramArg(const std::string& arg) { control_program_arg_ = arg; }
    const std::string& getControlProgramArg() const { return control_program_arg_; }
    void setFastReachSteadyState(bool fast) { fast_reach_steady_state_ = fast; }
    bool isFastReachSteadyState() const { return fast_reach_steady_state_; }
    void setAutoDetectSteadyState(bool auto_detect) { auto_detect_steady_state_ = auto_detect; }
    bool isAutoDetectSteadyState() const { return auto_detect_steady_state_; }
    void setFrequencyOfAddedVoltageSource(const std::string& freq) { frequency_of_added_voltage_source_ = freq; }
    const std::string& getFrequencyOfAddedVoltageSource() const { return frequency_of_added_voltage_source_; }
    void setStopCriterion(double criterion) { stop_criterion_ = criterion; }
    double getStopCriterion() const { return stop_criterion_; }
    void setIsGeneralTransient(bool is_general) { is_general_transient_ = is_general; }
    bool isGeneralTransient() const { return is_general_transient_; }
    void setIsHalfPeriodicTransient(bool is_half) { is_half_periodic_transient_ = is_half; }
    bool isHalfPeriodicTransient() const { return is_half_periodic_transient_; }
    void setSaveFieldsType(const std::string& type) { save_fields_type_ = type; }
    const std::string& getSaveFieldsType() const { return save_fields_type_; }
    void setNSteps(const std::string& n) { n_steps_ = n; }
    const std::string& getNSteps() const { return n_steps_; }
    void setStepsFrom(const std::string& from) { steps_from_ = from; }
    const std::string& getStepsFrom() const { return steps_from_; }
    void setStepsTo(const std::string& to) { steps_to_ = to; }
    const std::string& getStepsTo() const { return steps_to_; }
    void setUseNonlinearIterNum(bool use) { use_nonlinear_iter_num_ = use; }
    bool isUseNonlinearIterNum() const { return use_nonlinear_iter_num_; }
    void setCacheSaveKind(const std::string& kind) { cache_save_kind_ = kind; }
    const std::string& getCacheSaveKind() const { return cache_save_kind_; }
    void setNumberSolveSteps(int n) { number_solve_steps_ = n; }
    int getNumberSolveSteps() const { return number_solve_steps_; }
    void setRangeStart(const std::string& start) { range_start_ = start; }
    const std::string& getRangeStart() const { return range_start_; }
    void setRangeEnd(const std::string& end) { range_end_ = end; }
    const std::string& getRangeEnd() const { return range_end_; }
    void setUseAdaptiveTimeStep(bool use) { use_adaptive_time_step_ = use; }
    bool isUseAdaptiveTimeStep() const { return use_adaptive_time_step_; }
    void setInitialTimeStep(const std::string& step) { initial_time_step_ = step; }
    const std::string& getInitialTimeStep() const { return initial_time_step_; }
    void setMinTimeStep(const std::string& step) { min_time_step_ = step; }
    const std::string& getMinTimeStep() const { return min_time_step_; }
    void setMaxTimeStep(const std::string& step) { max_time_step_ = step; }
    const std::string& getMaxTimeStep() const { return max_time_step_; }
    void setTimeStepErrTolerance(double tol) { time_step_err_tolerance_ = tol; }
    double getTimeStepErrTolerance() const { return time_step_err_tolerance_; }

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

    // Maxwell专属求解设置成员
    HPCParallelMode hpc_parallel_mode_ = HPCParallelMode::SERIAL;
    HPCSolverMode hpc_solver_mode_ = HPCSolverMode::SHARED_MEMORY;
    int num_cores_ = 1;
    DomainDecompositionType domain_decomposition_type_ = DomainDecompositionType::GEOMETRIC;
    bool adaptive_mesh_refinement_ = false;
    int adaptive_depth_ = 0;
    bool skin_depth_refinement_ = false;
    bool coreloss_refinement_ = false;
    std::string maxwell_solver_id_;
    std::unordered_map<std::string, std::string> hpc_params_;
    std::unordered_map<std::string, std::string> maxwell_specific_params_;

    // 瞬态求解扩展成员
    bool enabled_ = true;
    bool mesh_link_import_ = false;
    int time_integration_method_ = 0;
    bool smooth_bh_curve_ = false;
    bool output_error_ = false;
    bool output_per_object_core_loss_ = false;
    bool output_per_object_solid_loss_ = false;
    bool use_control_program_ = false;
    std::string control_program_name_;
    std::string control_program_arg_;
    bool fast_reach_steady_state_ = true;
    bool auto_detect_steady_state_ = true;
    std::string frequency_of_added_voltage_source_;
    double stop_criterion_ = 0.005;
    bool is_general_transient_ = true;
    bool is_half_periodic_transient_ = false;
    std::string save_fields_type_;
    std::string n_steps_;
    std::string steps_from_;
    std::string steps_to_;
    bool use_nonlinear_iter_num_ = false;
    std::string cache_save_kind_;
    int number_solve_steps_ = 0;
    std::string range_start_;
    std::string range_end_;
    bool use_adaptive_time_step_ = false;
    std::string initial_time_step_;
    std::string min_time_step_;
    std::string max_time_step_;
    double time_step_err_tolerance_ = 0.0;
};

// ==================== 绕组设置类 ====================

class Winding : public ISerializable {
public:
    Winding(const std::string& name);
    Winding() = default;
    ~Winding() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    void setExcitationType(WindingExcitationType type) { excitation_type_ = type; }
    WindingExcitationType getExcitationType() const { return excitation_type_; }

    void setIsSolid(bool is_solid) { is_solid_ = is_solid; }
    bool isSolid() const { return is_solid_; }

    void setCurrentExpression(const std::string& expr) { current_expression_ = expr; }
    const std::string& getCurrentExpression() const { return current_expression_; }
    void setVoltageExpression(const std::string& expr) { voltage_expression_ = expr; }
    const std::string& getVoltageExpression() const { return voltage_expression_; }

    void setResistance(const std::string& r) { resistance_ = r; }
    const std::string& getResistance() const { return resistance_; }
    void setInductance(const std::string& l) { inductance_ = l; }
    const std::string& getInductance() const { return inductance_; }
    void setParallelBranchesNum(const std::string& n) { parallel_branches_num_ = n; }
    const std::string& getParallelBranchesNum() const { return parallel_branches_num_; }

    void addCoil(const std::string& coil_name);
    const std::vector<std::string>& getCoils() const { return coils_; }

private:
    std::string name_;
    uint64_t id_ = 0;
    WindingExcitationType excitation_type_ = WindingExcitationType::CURRENT;
    bool is_solid_ = false;
    std::string current_expression_;
    std::string voltage_expression_;
    std::string resistance_;
    std::string inductance_;
    std::string parallel_branches_num_;
    std::vector<std::string> coils_;
};
using WindingPtr = std::shared_ptr<Winding>;

// ==================== 运动设置类 ====================

class MotionSetup : public ISerializable {
public:
    MotionSetup(const std::string& name);
    MotionSetup() = default;
    ~MotionSetup() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setMotionType(MotionSetupType type) { motion_type_ = type; }
    MotionSetupType getMotionType() const { return motion_type_; }
    void setMoveType(MoveType type) { move_type_ = type; }
    MoveType getMoveType() const { return move_type_; }
    void setAxis(MotionAxis axis) { axis_ = axis; }
    MotionAxis getAxis() const { return axis_; }
    void setInitialPosition(const std::string& pos) { initial_position_ = pos; }
    const std::string& getInitialPosition() const { return initial_position_; }
    void setAngularVelocity(const std::string& vel) { angular_velocity_ = vel; }
    const std::string& getAngularVelocity() const { return angular_velocity_; }
    void setBandNameRef(int ref) { band_name_ref_ = ref; }
    int getBandNameRef() const { return band_name_ref_; }
    void addMovingObject(int obj_id);
    void addObject(int obj_id);
    const std::vector<int>& getObjects() const { return objects_; }
    const std::vector<int>& getMovingObjects() const { return moving_objects_; }

private:
    std::string name_;
    uint64_t id_ = 0;
    MotionSetupType motion_type_ = MotionSetupType::BAND;
    MoveType move_type_ = MoveType::ROTATE;
    int coordinate_system_ = 1;
    MotionAxis axis_ = MotionAxis::Z;
    bool is_positive_ = true;
    std::string initial_position_;
    bool has_rotate_limit_ = false;
    bool non_cylindrical_ = false;
    bool consider_mechanical_transient_ = false;
    std::string angular_velocity_;
    std::string linear_velocity_;
    std::vector<int> objects_;
    int band_name_ref_ = -1;
    std::vector<int> moving_objects_;
};
using MotionSetupPtr = std::shared_ptr<MotionSetup>;

// ==================== 网格操作类 ====================

class MeshOperation : public ISerializable {
public:
    MeshOperation(const std::string& name);
    MeshOperation() = default;
    ~MeshOperation() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    void setType(MeshOperationType type) { type_ = type; }
    MeshOperationType getType() const { return type_; }
    void setEnabled(bool enabled) { enabled_ = enabled; }
    bool isEnabled() const { return enabled_; }
    void addObject(int obj_id);
    const std::vector<int>& getObjects() const { return objects_; }

    // LengthBased特有方法
    void setMaxLength(const std::string& len) { max_length_ = len; }
    void setNumMaxElem(const std::string& num) { num_max_elem_ = num; }

    // SurfApproxBased特有方法
    void setSurfDev(const std::string& dev) { surf_dev_ = dev; }
    void setNormalDev(const std::string& dev) { normal_dev_ = dev; }

    // CylindricalGap特有方法
    void setBandMappingAngle(const std::string& angle) { band_mapping_angle_ = angle; }

private:
    std::string name_;
    int id_ = 0;
    MeshOperationType type_ = MeshOperationType::LENGTH_BASED;
    bool enabled_ = true;
    bool is_component_ = false;
    bool is_global_ = false;
    std::vector<int> objects_;
    bool refine_inside_ = true;
    bool restrict_elem_ = false;
    std::string num_max_elem_;
    bool restrict_length_ = true;
    std::string max_length_;
    bool apply_to_initial_mesh_ = false;

    // SurfApproxBased 特有字段
    std::string surf_approx_mode_ = "ManualSettings";
    int surf_dev_choice_ = 2;
    std::string surf_dev_;
    int normal_dev_choice_ = 2;
    std::string normal_dev_;
    int aspect_ratio_choice_ = 1;

    // CylindricalGap 特有字段
    bool use_band_mapping_angle_ = false;
    std::string band_mapping_angle_;
};
using MeshOperationPtr = std::shared_ptr<MeshOperation>;

// ==================== 设计变量类 ====================

class DesignVariable : public ISerializable {
public:
    DesignVariable(const std::string& name, const std::string& value,
                   const std::string& unit = "");
    DesignVariable() = default;
    ~DesignVariable() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    const std::string& getValue() const { return value_; }
    const std::string& getUnit() const { return unit_; }
    void setIndexed(bool indexed) { is_indexed_ = indexed; }
    bool isIndexed() const { return is_indexed_; }
    void setExpression(const std::string& expr) { expression_ = expr; }
    const std::string& getExpression() const { return expression_; }

private:
    std::string name_;
    std::string value_;
    std::string unit_;
    bool is_indexed_ = false;
    std::string expression_;
};
using DesignVariablePtr = std::shared_ptr<DesignVariable>;

// ==================== 输出变量类 ====================

class OutputVariable : public ISerializable {
public:
    OutputVariable(const std::string& name, int id, const std::string& expression,
                   const std::string& result_unit = "", const std::string& display_unit = "");
    OutputVariable() = default;
    ~OutputVariable() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    const std::string& getName() const { return name_; }
    int getID() const { return id_; }
    const std::string& getExpression() const { return expression_; }
    const std::string& getResultUnit() const { return result_unit_; }
    const std::string& getDisplayUnit() const { return display_unit_; }

private:
    std::string name_;
    int id_ = 0;
    std::string expression_;
    std::string result_unit_;
    std::string display_unit_;
};
using OutputVariablePtr = std::shared_ptr<OutputVariable>;

// ==================== 温度设置类 ====================

class TemperatureSettings : public ISerializable {
public:
    TemperatureSettings();
    ~TemperatureSettings() override = default;

    nlohmann::json toJson() const override;
    bool fromJson(const nlohmann::json& json) override;
    bool toBinary(std::vector<uint8_t>& data) const override;
    bool fromBinary(const std::vector<uint8_t>& data, size_t& offset) override;
    uint32_t getSerializationVersion() const override { return 1; }
    bool validate() const override;

    void setIncludeTemperatureDependence(bool include) { include_temperature_dependent_ = include; }
    bool isIncludeTemperatureDependence() const { return include_temperature_dependent_; }
    void setEnableFeedback(bool enable) { enable_feedback_ = enable; }
    bool isEnableFeedback() const { return enable_feedback_; }
    void setObjectTemperature(int object_id, const std::string& temp_ref);
    const std::unordered_map<int, std::string>& getObjectTemperatureMap() { return object_temperature_map_; }

private:
    bool include_temperature_dependent_ = false;
    bool enable_feedback_ = false;
    std::unordered_map<int, std::string> object_temperature_map_;
};
using TemperatureSettingsPtr = std::shared_ptr<TemperatureSettings>;

} // namespace tool
