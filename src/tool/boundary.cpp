/**
 * @file boundary.cpp
 * @brief Boundary类实现文件
 * @details 实现Boundary类的序列化、反序列化及边界条件管理功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/project_data.hpp"
#include "tool/em_enums.hpp"
#include "tool/id_generator.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace tool {

Boundary::Boundary(const std::string& name) : name_(name) {
    id_ = IDGenerator::getInstance().generateID(IDCategory::BOUNDARY);
}

nlohmann::json Boundary::toJson() const {
    nlohmann::json json;
    
    // 基础信息
    json["name"] = name_;
    json["id"] = id_;
    json["type"] = static_cast<int>(type_);
    
    // 几何关联
    json["faces"] = faces_;
    json["edges"] = edges_;
    json["objects"] = objects_;
    
    // 边界参数
    json["impedance_value"] = impedance_value_;
    json["vector_potential"] = vector_potential_;
    json["voltage"] = voltage_;
    json["current"] = current_;
    
    // 主从边界
    json["master_name"] = master_name_;
    json["slave_name"] = slave_name_;
    
    // Maxwell专属边界数据
    json["boundary_sub_type"] = static_cast<int>(boundary_sub_type_);
    json["periodic_mapping_type"] = static_cast<int>(periodic_mapping_type_);
    json["radiation_distance"] = radiation_distance_;
    json["perfect_e_symmetry"] = perfect_e_symmetry_;
    json["perfect_h_symmetry"] = perfect_h_symmetry_;
    json["infinite_sphere_radius"] = infinite_sphere_radius_;
    json["maxwell_boundary_id"] = maxwell_boundary_id_;
    json["boundary_subdivision_params"] = boundary_subdivision_params_;
    json["maxwell_specific_params"] = maxwell_specific_params_;
    
    return json;
}

bool Boundary::fromJson(const nlohmann::json& json) {
    try {
        // 基础信息
        name_ = json.value("name", "");
        id_ = json.value("id", 0ULL);
        type_ = static_cast<BndType>(json.value("type", 0));
        
        // 几何关联
        faces_ = json.value("faces", std::vector<std::string>());
        edges_ = json.value("edges", std::vector<std::string>());
        objects_ = json.value("objects", std::vector<std::string>());
        
        // 边界参数
        impedance_value_ = json.value("impedance_value", 0.0);
        vector_potential_ = json.value("vector_potential", 0.0);
        voltage_ = json.value("voltage", 0.0);
        current_ = json.value("current", 0.0);
        
        // 主从边界
        master_name_ = json.value("master_name", "");
        slave_name_ = json.value("slave_name", "");
        
        // Maxwell专属边界数据
        boundary_sub_type_ = static_cast<BoundarySubType>(json.value("boundary_sub_type", 0));
        periodic_mapping_type_ = static_cast<PeriodicMappingType>(json.value("periodic_mapping_type", 0));
        radiation_distance_ = json.value("radiation_distance", 0.0);
        perfect_e_symmetry_ = json.value("perfect_e_symmetry", false);
        perfect_h_symmetry_ = json.value("perfect_h_symmetry", false);
        infinite_sphere_radius_ = json.value("infinite_sphere_radius", 0.0);
        maxwell_boundary_id_ = json.value("maxwell_boundary_id", "");
        boundary_subdivision_params_ = json.value("boundary_subdivision_params", std::vector<double>());
        maxwell_specific_params_ = json.value("maxwell_specific_params", std::unordered_map<std::string, std::string>());
        
        return validate();
    } catch (const std::exception& e) {
        // 记录错误日志
        return false;
    }
}

bool Boundary::toBinary(std::vector<uint8_t>& data) const {
    try {
        // 将对象转换为JSON
        nlohmann::json j = toJson();
        
        // 序列化JSON为字符串
        std::string json_str = j.dump();
        
        // 写入版本信息
        uint32_t version = getSerializationVersion();
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&version), 
                   reinterpret_cast<const uint8_t*>(&version) + sizeof(uint32_t));
        
        // 写入数据长度
        uint32_t data_size = static_cast<uint32_t>(json_str.size());
        data.insert(data.end(), reinterpret_cast<const uint8_t*>(&data_size), 
                   reinterpret_cast<const uint8_t*>(&data_size) + sizeof(uint32_t));
        
        // 写入JSON数据
        data.insert(data.end(), json_str.begin(), json_str.end());
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

bool Boundary::fromBinary(const std::vector<uint8_t>& data, size_t& offset) {
    try {
        // 读取版本信息
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t version = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取数据长度
        if (offset + sizeof(uint32_t) > data.size()) return false;
        uint32_t data_size = *reinterpret_cast<const uint32_t*>(&data[offset]);
        offset += sizeof(uint32_t);
        
        // 读取JSON数据
        if (offset + data_size > data.size()) return false;
        std::string json_str(data.begin() + offset, data.begin() + offset + data_size);
        offset += data_size;
        
        // 解析JSON
        nlohmann::json j = nlohmann::json::parse(json_str);
        return fromJson(j);
    } catch (const std::exception& e) {
        return false;
    }
}

bool Boundary::validate() const {
    // 验证边界名称
    if (name_.empty()) return false;
    
    // 验证边界类型
    if (type_ < BndType::DIRICHLET || type_ > BndType::CONTACTS) return false;
    
    // 验证参数范围
    if (impedance_value_ < 0.0) return false;
    if (voltage_ < 0.0) return false;
    if (current_ < 0.0) return false;
    
    // 验证Maxwell专属参数
    if (radiation_distance_ < 0.0) return false;
    if (infinite_sphere_radius_ < 0.0) return false;
    
    return true;
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

void Boundary::setBoundarySubdivisionParameters(const std::vector<double>& params) {
    boundary_subdivision_params_ = params;
}

void Boundary::setMaxwellSpecificParameters(const std::unordered_map<std::string, std::string>& params) {
    maxwell_specific_params_ = params;
}

} // namespace tool