/**
 * @file association_map.hpp
 * @brief 基础工具层 - 关联映射表头文件
 * @details 实现几何实体与材料/边界/激励之间的关联关系管理
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <optional>

namespace tool {

class AssociationMap {
public:
    AssociationMap();
    ~AssociationMap() = default;
    
    bool link_entity_material(const std::string& entity_id, const std::string& material_id);
    bool link_entity_boundary(const std::string& entity_id, const std::string& boundary_id);
    bool link_entity_excitation(const std::string& entity_id, const std::string& excitation_id);
    
    bool unlink_entity_material(const std::string& entity_id);
    bool unlink_entity_boundary(const std::string& entity_id, const std::string& boundary_id);
    bool unlink_entity_excitation(const std::string& entity_id);
    
    std::optional<std::string> get_material_id(const std::string& entity_id) const;
    std::vector<std::string> get_boundary_ids(const std::string& entity_id) const;
    std::optional<std::string> get_excitation_id(const std::string& entity_id) const;
    
    std::vector<std::string> get_entities_with_material(const std::string& material_id) const;
    std::vector<std::string> get_entities_with_boundary(const std::string& boundary_id) const;
    std::vector<std::string> get_entities_with_excitation(const std::string& excitation_id) const;
    
    bool check_dangling_references(const ProjectManager& project, 
                                  std::vector<std::string>& dangling_materials,
                                  std::vector<std::string>& dangling_boundaries,
                                  std::vector<std::string>& dangling_excitations) const;
    
    void clear();
    void clear_entity(const std::string& entity_id);
    
    size_t get_material_link_count() const { return entity_to_material_.size(); }
    size_t get_boundary_link_count() const {
        size_t count = 0;
        for (const auto& pair : entity_to_boundaries_) {
            count += pair.second.size();
        }
        return count;
    }
    size_t get_excitation_link_count() const { return entity_to_excitation_.size(); }
    
    bool serialize(const std::string& file_path) const;
    bool deserialize(const std::string& file_path);
    
    std::string to_string() const;
    
private:
    std::unordered_map<std::string, std::string> entity_to_material_;
    std::unordered_map<std::string, std::unordered_set<std::string>> entity_to_boundaries_;
    std::unordered_map<std::string, std::string> entity_to_excitation_;
    
    std::unordered_map<std::string, std::unordered_set<std::string>> material_to_entities_;
    std::unordered_map<std::string, std::unordered_set<std::string>> boundary_to_entities_;
    std::unordered_map<std::string, std::unordered_set<std::string>> excitation_to_entities_;
};

} // namespace tool
