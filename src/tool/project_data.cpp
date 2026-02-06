/**
 * @file project_data.cpp
 * @brief 基础工具层 - 项目数据结构源文件
 */

#include "tool/project_data.hpp"
#include "tool/em_exception.hpp"
#include <algorithm>

namespace tool {

Material::Material(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateMaterialID(name)) {
}

void Material::addProperty(const MaterialProperty& property) {
    auto it = property_index_.find(property.name);
    if (it != property_index_.end()) {
        properties_[it->second] = property;
    } else {
        property_index_[property.name] = static_cast<int>(properties_.size());
        properties_.push_back(property);
    }
}

std::optional<MaterialProperty> Material::getProperty(const std::string& name) const {
    auto it = property_index_.find(name);
    if (it != property_index_.end()) {
        return properties_[it->second];
    }
    return std::nullopt;
}

void Material::setBHCurve(const std::vector<BHDataPoint>& curve) {
    bh_curve_ = curve;
    bh_curve_type_ = BHCurveType::SINGLE_CURVE;
}

void Material::setCoreLossCoefficients(double ks, double alpha, double beta, double kn) {
    core_loss_ks_ = ks;
    core_loss_alpha_ = alpha;
    core_loss_beta_ = beta;
    core_loss_kn_ = kn;
    core_loss_enabled_ = true;
}

void Material::getCoreLossCoefficients(double& ks, double& alpha, double& beta, double& kn) const {
    ks = core_loss_ks_;
    alpha = core_loss_alpha_;
    beta = core_loss_beta_;
    kn = core_loss_kn_;
}

void Material::setRelativePermeability(double mu_r) {
    relative_permeability_ = mu_r;
}

void Material::setConductivity(double sigma) {
    conductivity_ = sigma;
}

void Material::setMassDensity(double rho) {
    mass_density_ = rho;
}

Geometry::Geometry(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateEntityID(name)) {
}

void Geometry::addSubGeometry(const std::string& name, GeometryPtr sub_geo) {
    sub_geometries_[name] = sub_geo;
}

std::optional<GeometryPtr> Geometry::getSubGeometry(const std::string& name) const {
    auto it = sub_geometries_.find(name);
    if (it != sub_geometries_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void Geometry::addObject(const std::string& name, const std::string& material_name) {
    object_material_map_[name] = material_name;
}

void Geometry::assignMaterial(const std::string& object_name, const std::string& material_name) {
    object_material_map_[object_name] = material_name;
}

Boundary::Boundary(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateBoundaryID(name)) {
}

void Boundary::addFace(const std::string& face_id) {
    faces_.push_back(face_id);
}

void Boundary::addEdge(const std::string& edge_id) {
    edges_.push_back(edge_id);
}

void Boundary::addObject(const std::string& object_name) {
    objects_.push_back(object_name);
}

Excitation::Excitation(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateExcitationID(name)) {
}

Mesh::Mesh(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateEntityID(name)) {
}

void Mesh::addObjectMeshSettings(const std::string& object_name, double min, double max) {
    object_mesh_settings_[object_name] = std::make_pair(min, max);
}

SolutionSetup::SolutionSetup(const std::string& name)
    : name_(name), id_(EntityIDGenerator::getInstance().generateEntityID(name)) {
}

} // namespace tool
