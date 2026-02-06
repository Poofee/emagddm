/**
 * @file version_manager.cpp
 * @brief 基础工具层 - 版本管理源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/version_manager.hpp"
#include "tool/file_utils.hpp"
#include "tool/xml_interface.hpp"
#include <chrono>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <random>

namespace tool {

std::string ProjectVersion::to_string() const {
    std::ostringstream oss;
    oss << "Version[" << version_id << "] " << version_name;
    oss << " (v" << project_version << ")";
    oss << " - " << version_description;
    return oss.str();
}

void VersionDiff::add_entry(const DiffEntry& entry) {
    entries_.push_back(entry);
}

void VersionDiff::add_new_entity(const std::string& data_type, const std::string& entity_id) {
    DiffEntry entry;
    entry.data_type = data_type;
    entry.entity_id = entity_id;
    entry.is_new = true;
    entries_.push_back(entry);
}

void VersionDiff::add_deleted_entity(const std::string& data_type, const std::string& entity_id) {
    DiffEntry entry;
    entry.data_type = data_type;
    entry.entity_id = entity_id;
    entry.is_deleted = true;
    entries_.push_back(entry);
}

void VersionDiff::add_modified_field(const std::string& data_type, const std::string& entity_id,
                                    const std::string& field_name, const std::string& old_value,
                                    const std::string& new_value) {
    DiffEntry entry;
    entry.data_type = data_type;
    entry.entity_id = entity_id;
    entry.field_name = field_name;
    entry.old_value = old_value;
    entry.new_value = new_value;
    entries_.push_back(entry);
}

size_t VersionDiff::get_new_count() const {
    size_t count = 0;
    for (const auto& entry : entries_) {
        if (entry.is_new) count++;
    }
    return count;
}

size_t VersionDiff::get_deleted_count() const {
    size_t count = 0;
    for (const auto& entry : entries_) {
        if (entry.is_deleted) count++;
    }
    return count;
}

size_t VersionDiff::get_modified_count() const {
    size_t count = 0;
    for (const auto& entry : entries_) {
        if (!entry.is_new && !entry.is_deleted && !entry.field_name.empty()) count++;
    }
    return count;
}

std::string VersionDiff::to_string() const {
    std::ostringstream oss;
    oss << "VersionDiff:\n";
    oss << "  Total changes: " << entries_.size() << "\n";
    oss << "  New entities: " << get_new_count() << "\n";
    oss << "  Deleted entities: " << get_deleted_count() << "\n";
    oss << "  Modified fields: " << get_modified_count() << "\n";
    return oss.str();
}

VersionManager::VersionManager(const std::string& storage_directory)
    : storage_directory_(storage_directory) {}

VersionManager::~VersionManager() {
    shutdown();
}

bool VersionManager::initialize() {
    if (!file_utils::exists(storage_directory_)) {
        if (!file_utils::createDirectories(storage_directory_)) {
            return false;
        }
    }
    
    load_version_metadata();
    
    return true;
}

void VersionManager::shutdown() {
}

uint64_t VersionManager::create_version(ProjectManager& project, const std::string& name,
                                      const std::string& description, VersionType type) {
    current_project_version_++;
    
    ProjectVersion version;
    version.version_id = versions_.size() + 1;
    version.version_name = name;
    version.version_description = description;
    version.type = type;
    version.creation_time = std::chrono::system_clock::now();
    version.modification_time = version.creation_time;
    version.project_version = current_project_version_;
    version.created_by = "";
    
    version.checksum = calculate_project_checksum(project);
    
    std::string version_path = generate_version_path(version.version_id);
    FragmentedWriter writer(version_path);
    if (!writer.create()) {
        return 0;
    }
    
    writer.write_string("version_name", name);
    writer.write_string("version_description", description);
    writer.write_string("type", type == VersionType::MANUAL ? "MANUAL" : 
                       type == VersionType::AUTOMATIC ? "AUTOMATIC" : "SNAPSHOT");
    writer.finalize();
    
    versions_.push_back(version);
    save_version_metadata(version);
    
    if (type == VersionType::AUTOMATIC) {
        size_t auto_count = get_automatic_backup_count();
        if (auto_count > max_automatic_backups_) {
            delete_old_versions(max_automatic_backups_);
        }
    }
    
    return version.version_id;
}

uint64_t VersionManager::create_automatic_backup(ProjectManager& project) {
    std::string name = "AutoBackup_v" + std::to_string(current_project_version_ + 1);
    std::string description = "Automatic backup created at " + 
        std::to_string(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
    return create_version(project, name, description, VersionType::AUTOMATIC);
}

std::optional<ProjectVersion> VersionManager::get_version(uint64_t version_id) const {
    for (const auto& version : versions_) {
        if (version.version_id == version_id) {
            return version;
        }
    }
    return std::nullopt;
}

std::optional<ProjectVersion> VersionManager::get_version_by_name(const std::string& name) const {
    for (const auto& version : versions_) {
        if (version.version_name == name) {
            return version;
        }
    }
    return std::nullopt;
}

std::vector<ProjectVersion> VersionManager::get_all_versions() const {
    return versions_;
}

std::vector<ProjectVersion> VersionManager::get_versions_by_type(VersionType type) const {
    std::vector<ProjectVersion> result;
    for (const auto& version : versions_) {
        if (version.type == type) {
            result.push_back(version);
        }
    }
    return result;
}

std::vector<ProjectVersion> VersionManager::get_recent_versions(size_t count) const {
    std::vector<ProjectVersion> result;
    size_t start = count >= versions_.size() ? 0 : versions_.size() - count;
    for (size_t i = start; i < versions_.size(); ++i) {
        result.push_back(versions_[i]);
    }
    return result;
}

bool VersionManager::restore_version(uint64_t version_id, ProjectManager& project) {
    auto version = get_version(version_id);
    if (!version.has_value()) {
        return false;
    }
    
    std::string version_path = generate_version_path(version_id);
    if (!file_utils::exists(version_path)) {
        return false;
    }
    
    return true;
}

bool VersionManager::rollback_to_version(uint64_t version_id, ProjectManager& project) {
    return restore_version(version_id, project);
}

std::unique_ptr<VersionDiff> VersionManager::compare_versions(uint64_t version1_id, 
                                                           uint64_t version2_id) const {
    auto diff = std::make_unique<VersionDiff>();
    
    auto version1 = get_version(version1_id);
    auto version2 = get_version(version2_id);
    
    if (!version1.has_value() || !version2.has_value()) {
        return diff;
    }
    
    return diff;
}

std::unique_ptr<VersionDiff> VersionManager::compare_with_current(uint64_t version_id,
                                                                 const ProjectManager& project) const {
    return compare_versions(version_id, versions_.back().version_id);
}

bool VersionManager::delete_version(uint64_t version_id) {
    auto it = versions_.begin();
    while (it != versions_.end()) {
        if (it->version_id == version_id) {
            std::string version_path = generate_version_path(version_id);
            if (file_utils::exists(version_path)) {
                file_utils::remove(version_path);
            }
            versions_.erase(it);
            return true;
        }
        ++it;
    }
    return false;
}

bool VersionManager::delete_old_versions(size_t keep_count) {
    if (versions_.size() <= keep_count) {
        return true;
    }
    
    std::vector<ProjectVersion> to_delete;
    for (size_t i = 0; i < versions_.size() - keep_count; ++i) {
        if (versions_[i].type == VersionType::AUTOMATIC) {
            to_delete.push_back(versions_[i]);
        }
    }
    
    for (const auto& version : to_delete) {
        delete_version(version.version_id);
    }
    
    return true;
}

bool VersionManager::delete_versions_older_than(const std::chrono::system_clock::time_point& cutoff_time) {
    std::vector<ProjectVersion> to_delete;
    for (const auto& version : versions_) {
        if (version.creation_time < cutoff_time) {
            to_delete.push_back(version);
        }
    }
    
    for (const auto& version : to_delete) {
        delete_version(version.version_id);
    }
    
    return true;
}

bool VersionManager::export_version(uint64_t version_id, const std::string& file_path) {
    std::string version_path = generate_version_path(version_id);
    return file_utils::copyFile(version_path, file_path);
}

bool VersionManager::import_version(const std::string& file_path, ProjectManager& project) {
    if (!file_utils::exists(file_path)) {
        return false;
    }
    
    std::string file_name = file_utils::getFilename(file_path);
    
    return true;
}

size_t VersionManager::get_automatic_backup_count() const {
    size_t count = 0;
    for (const auto& version : versions_) {
        if (version.type == VersionType::AUTOMATIC) {
            count++;
        }
    }
    return count;
}

std::string VersionManager::generate_version_path(uint64_t version_id) const {
    return file_utils::combinePath(storage_directory_, "v" + std::to_string(version_id) + ".ver");
}

bool VersionManager::save_version_metadata(const ProjectVersion& version) {
    std::string meta_path = file_utils::combinePath(storage_directory_, "versions.xml");
    
    std::ofstream file(meta_path);
    if (!file.is_open()) {
        return false;
    }
    
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<VersionMetadata>\n";
    
    for (const auto& v : versions_) {
        file << "  <Version id=\"" << v.version_id << "\">\n";
        file << "    <Name>" << v.version_name << "</Name>\n";
        file << "    <Description>" << v.version_description << "</Description>\n";
        file << "    <Type>" << (v.type == VersionType::MANUAL ? "MANUAL" : 
                                  v.type == VersionType::AUTOMATIC ? "AUTOMATIC" : "SNAPSHOT") << "</Type>\n";
        file << "    <ProjectVersion>" << v.project_version << "</ProjectVersion>\n";
        file << "  </Version>\n";
    }
    
    file << "</VersionMetadata>\n";
    
    file.close();
    return true;
}

bool VersionManager::load_version_metadata() {
    std::string meta_path = file_utils::combinePath(storage_directory_, "versions.xml");
    
    if (!file_utils::exists(meta_path)) {
        return true;
    }
    
    return true;
}

std::string VersionManager::calculate_project_checksum(const ProjectManager& project) {
    std::ostringstream oss;
    oss << project.getProjectName();
    oss << project.getProjectFilePath();
    oss << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
    return std::to_string(std::hash<std::string>{}(oss.str()));
}

VersionManager& VersionManagerSingleton::getInstance() {
    static VersionManager instance;
    return instance;
}

} // namespace tool
