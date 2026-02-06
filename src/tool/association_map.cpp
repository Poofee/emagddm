/**
 * @file association_map.cpp
 * @brief 基础工具层 - 关联映射表源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/association_map.hpp"
#include "tool/project_manager.hpp"
#include "tool/file_utils.hpp"
#include <fstream>
#include <sstream>

namespace tool {

AssociationMap::AssociationMap() = default;

bool AssociationMap::link_entity_material(const std::string& entity_id, const std::string& material_id) {
    if (entity_id.empty() || material_id.empty()) {
        return false;
    }
    
    auto it = entity_to_material_.find(entity_id);
    if (it != entity_to_material_.end()) {
        material_to_entities_[it->second].erase(entity_id);
    }
    
    entity_to_material_[entity_id] = material_id;
    material_to_entities_[material_id].insert(entity_id);
    
    return true;
}

bool AssociationMap::link_entity_boundary(const std::string& entity_id, const std::string& boundary_id) {
    if (entity_id.empty() || boundary_id.empty()) {
        return false;
    }
    
    entity_to_boundaries_[entity_id].insert(boundary_id);
    boundary_to_entities_[boundary_id].insert(entity_id);
    
    return true;
}

bool AssociationMap::link_entity_excitation(const std::string& entity_id, const std::string& excitation_id) {
    if (entity_id.empty() || excitation_id.empty()) {
        return false;
    }
    
    auto it = entity_to_excitation_.find(entity_id);
    if (it != entity_to_excitation_.end()) {
        excitation_to_entities_[it->second].erase(entity_id);
    }
    
    entity_to_excitation_[entity_id] = excitation_id;
    excitation_to_entities_[excitation_id].insert(entity_id);
    
    return true;
}

bool AssociationMap::unlink_entity_material(const std::string& entity_id) {
    auto it = entity_to_material_.find(entity_id);
    if (it == entity_to_material_.end()) {
        return false;
    }
    
    material_to_entities_[it->second].erase(entity_id);
    entity_to_material_.erase(it);
    
    return true;
}

bool AssociationMap::unlink_entity_boundary(const std::string& entity_id, const std::string& boundary_id) {
    auto it = entity_to_boundaries_.find(entity_id);
    if (it == entity_to_boundaries_.end()) {
        return false;
    }
    
    it->second.erase(boundary_id);
    boundary_to_entities_[boundary_id].erase(entity_id);
    
    return true;
}

bool AssociationMap::unlink_entity_excitation(const std::string& entity_id) {
    auto it = entity_to_excitation_.find(entity_id);
    if (it == entity_to_excitation_.end()) {
        return false;
    }
    
    excitation_to_entities_[it->second].erase(entity_id);
    entity_to_excitation_.erase(it);
    
    return true;
}

std::optional<std::string> AssociationMap::get_material_id(const std::string& entity_id) const {
    auto it = entity_to_material_.find(entity_id);
    if (it != entity_to_material_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::vector<std::string> AssociationMap::get_boundary_ids(const std::string& entity_id) const {
    auto it = entity_to_boundaries_.find(entity_id);
    if (it != entity_to_boundaries_.end()) {
        return std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return {};
}

std::optional<std::string> AssociationMap::get_excitation_id(const std::string& entity_id) const {
    auto it = entity_to_excitation_.find(entity_id);
    if (it != entity_to_excitation_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::vector<std::string> AssociationMap::get_entities_with_material(const std::string& material_id) const {
    auto it = material_to_entities_.find(material_id);
    if (it != material_to_entities_.end()) {
        return std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return {};
}

std::vector<std::string> AssociationMap::get_entities_with_boundary(const std::string& boundary_id) const {
    auto it = boundary_to_entities_.find(boundary_id);
    if (it != boundary_to_entities_.end()) {
        return std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return {};
}

std::vector<std::string> AssociationMap::get_entities_with_excitation(const std::string& excitation_id) const {
    auto it = excitation_to_entities_.find(excitation_id);
    if (it != excitation_to_entities_.end()) {
        return std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return {};
}

bool AssociationMap::check_dangling_references(const ProjectManager& project,
                                              std::vector<std::string>& dangling_materials,
                                              std::vector<std::string>& dangling_boundaries,
                                              std::vector<std::string>& dangling_excitations) const {
    dangling_materials.clear();
    dangling_boundaries.clear();
    dangling_excitations.clear();
    
    auto all_materials = project.getAllMaterials();
    for (const auto& mat_pair : all_materials) {
        const std::string& mat_id = mat_pair.first;
        if (material_to_entities_.find(mat_id) == material_to_entities_.end()) {
            dangling_materials.push_back(mat_id);
        }
    }
    
    auto all_boundaries = project.getAllBoundaries();
    for (const auto& bnd_pair : all_boundaries) {
        const std::string& bnd_id = bnd_pair.first;
        if (boundary_to_entities_.find(bnd_id) == boundary_to_entities_.end()) {
            dangling_boundaries.push_back(bnd_id);
        }
    }
    
    auto all_excitations = project.getAllExcitations();
    for (const auto& exc_pair : all_excitations) {
        const std::string& exc_id = exc_pair.first;
        if (excitation_to_entities_.find(exc_id) == excitation_to_entities_.end()) {
            dangling_excitations.push_back(exc_id);
        }
    }
    
    return dangling_materials.empty() && dangling_boundaries.empty() && dangling_excitations.empty();
}

void AssociationMap::clear() {
    entity_to_material_.clear();
    entity_to_boundaries_.clear();
    entity_to_excitation_.clear();
    material_to_entities_.clear();
    boundary_to_entities_.clear();
    excitation_to_entities_.clear();
}

void AssociationMap::clear_entity(const std::string& entity_id) {
    auto mat_it = entity_to_material_.find(entity_id);
    if (mat_it != entity_to_material_.end()) {
        material_to_entities_[mat_it->second].erase(entity_id);
        entity_to_material_.erase(mat_it);
    }
    
    auto bnd_it = entity_to_boundaries_.find(entity_id);
    if (bnd_it != entity_to_boundaries_.end()) {
        for (const auto& bnd_id : bnd_it->second) {
            boundary_to_entities_[bnd_id].erase(entity_id);
        }
        entity_to_boundaries_.erase(bnd_it);
    }
    
    auto exc_it = entity_to_excitation_.find(entity_id);
    if (exc_it != entity_to_excitation_.end()) {
        excitation_to_entities_[exc_it->second].erase(entity_id);
        entity_to_excitation_.erase(exc_it);
    }
}

bool AssociationMap::serialize(const std::string& file_path) const {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return false;
    }
    
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<AssociationMap>\n";
    
    file << "  <EntityMaterialLinks>\n";
    for (const auto& pair : entity_to_material_) {
        file << "    <Link entity=\"" << pair.first << "\" material=\"" << pair.second << "\"/>\n";
    }
    file << "  </EntityMaterialLinks>\n";
    
    file << "  <EntityBoundaryLinks>\n";
    for (const auto& pair : entity_to_boundaries_) {
        for (const auto& bnd_id : pair.second) {
            file << "    <Link entity=\"" << pair.first << "\" boundary=\"" << bnd_id << "\"/>\n";
        }
    }
    file << "  </EntityBoundaryLinks>\n";
    
    file << "  <EntityExcitationLinks>\n";
    for (const auto& pair : entity_to_excitation_) {
        file << "    <Link entity=\"" << pair.first << "\" excitation=\"" << pair.second << "\"/>\n";
    }
    file << "  </EntityExcitationLinks>\n";
    
    file << "</AssociationMap>\n";
    
    file.close();
    return true;
}

bool AssociationMap::deserialize(const std::string& file_path) {
    clear();
    
    std::ifstream file(file_path);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    std::string current_section;
    
    while (std::getline(file, line)) {
        if (line.find("<EntityMaterialLinks>") != std::string::npos) {
            current_section = "material";
            continue;
        } else if (line.find("<EntityBoundaryLinks>") != std::string::npos) {
            current_section = "boundary";
            continue;
        } else if (line.find("<EntityExcitationLinks>") != std::string::npos) {
            current_section = "excitation";
            continue;
        } else if (line.find("</") != std::string::npos) {
            current_section = "";
            continue;
        }
        
        if (line.find("<Link") != std::string::npos) {
            size_t entity_pos = line.find("entity=\"");
            size_t material_pos = line.find("material=\"");
            size_t boundary_pos = line.find("boundary=\"");
            size_t excitation_pos = line.find("excitation=\"");
            
            if (entity_pos == std::string::npos) continue;
            
            size_t entity_end = line.find("\"", entity_pos + 8);
            std::string entity_id = line.substr(entity_pos + 8, entity_end - entity_pos - 8);
            
            if (current_section == "material" && material_pos != std::string::npos) {
                size_t material_end = line.find("\"", material_pos + 9);
                std::string material_id = line.substr(material_pos + 9, material_end - material_pos - 9);
                link_entity_material(entity_id, material_id);
            } else if (current_section == "boundary" && boundary_pos != std::string::npos) {
                size_t boundary_end = line.find("\"", boundary_pos + 10);
                std::string boundary_id = line.substr(boundary_pos + 10, boundary_end - boundary_pos - 10);
                link_entity_boundary(entity_id, boundary_id);
            } else if (current_section == "excitation" && excitation_pos != std::string::npos) {
                size_t excitation_end = line.find("\"", excitation_pos + 12);
                std::string excitation_id = line.substr(excitation_pos + 12, excitation_end - excitation_pos - 12);
                link_entity_excitation(entity_id, excitation_id);
            }
        }
    }
    
    file.close();
    return true;
}

std::string AssociationMap::to_string() const {
    std::ostringstream oss;
    oss << "AssociationMap:\n";
    oss << "  Entity-Material links: " << entity_to_material_.size() << "\n";
    oss << "  Entity-Boundary links: " << get_boundary_link_count() << "\n";
    oss << "  Entity-Excitation links: " << entity_to_excitation_.size() << "\n";
    return oss.str();
}

} // namespace tool
