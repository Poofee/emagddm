/**
 * @file version_manager.hpp
 * @brief 基础工具层 - 版本管理头文件
 * @details 实现手动/自动版本创建、版本回滚功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include "tool/project_manager.hpp"
#include "tool/fragmented_storage.hpp"
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <optional>
#include <functional>

namespace tool {

enum class VersionType {
    MANUAL,
    AUTOMATIC,
    SNAPSHOT,
    UNKNOWN
};

struct ProjectVersion {
    uint64_t version_id;
    std::string version_name;
    std::string version_description;
    VersionType type;
    std::chrono::system_clock::time_point creation_time;
    std::chrono::system_clock::time_point modification_time;
    uint64_t project_version;
    std::string checksum;
    size_t data_size;
    std::string created_by;
    
    bool is_automatic_backup() const { return type == VersionType::AUTOMATIC; }
    
    std::string to_string() const;
};

class VersionDiff {
public:
    struct DiffEntry {
        std::string data_type;
        std::string entity_id;
        std::string field_name;
        std::string old_value;
        std::string new_value;
        bool is_new = false;
        bool is_deleted = false;
    };
    
    void add_entry(const DiffEntry& entry);
    void add_new_entity(const std::string& data_type, const std::string& entity_id);
    void add_deleted_entity(const std::string& data_type, const std::string& entity_id);
    void add_modified_field(const std::string& data_type, const std::string& entity_id,
                           const std::string& field_name, const std::string& old_value,
                           const std::string& new_value);
    
    const std::vector<DiffEntry>& get_entries() const { return entries_; }
    size_t get_entry_count() const { return entries_.size(); }
    size_t get_new_count() const;
    size_t get_deleted_count() const;
    size_t get_modified_count() const;
    
    std::string to_string() const;
    
private:
    std::vector<DiffEntry> entries_;
};

class VersionManager {
public:
    VersionManager(const std::string& storage_directory = "./versions");
    ~VersionManager();
    
    bool initialize();
    void shutdown();
    
    uint64_t create_version(ProjectManager& project, const std::string& name,
                           const std::string& description = "", VersionType type = VersionType::MANUAL);
    
    uint64_t create_automatic_backup(ProjectManager& project);
    
    std::optional<ProjectVersion> get_version(uint64_t version_id) const;
    std::optional<ProjectVersion> get_version_by_name(const std::string& name) const;
    
    std::vector<ProjectVersion> get_all_versions() const;
    std::vector<ProjectVersion> get_versions_by_type(VersionType type) const;
    std::vector<ProjectVersion> get_recent_versions(size_t count) const;
    
    bool restore_version(uint64_t version_id, ProjectManager& project);
    bool rollback_to_version(uint64_t version_id, ProjectManager& project);
    
    std::unique_ptr<VersionDiff> compare_versions(uint64_t version1_id, uint64_t version2_id) const;
    std::unique_ptr<VersionDiff> compare_with_current(uint64_t version_id, const ProjectManager& project) const;
    
    bool delete_version(uint64_t version_id);
    bool delete_old_versions(size_t keep_count);
    bool delete_versions_older_than(const std::chrono::system_clock::time_point& cutoff_time);
    
    bool export_version(uint64_t version_id, const std::string& file_path);
    bool import_version(const std::string& file_path, ProjectManager& project);
    
    void set_max_automatic_backups(size_t max_count) { max_automatic_backups_ = max_count; }
    size_t get_max_automatic_backups() const { return max_automatic_backups_; }
    
    size_t get_version_count() const { return versions_.size(); }
    size_t get_automatic_backup_count() const;
    
    std::string get_storage_directory() const { return storage_directory_; }
    
    using ProgressCallback = std::function<void(size_t current, size_t total)>;
    void set_progress_callback(ProgressCallback callback) { progress_callback_ = callback; }
    
private:
    std::string storage_directory_;
    std::vector<ProjectVersion> versions_;
    size_t max_automatic_backups_ = 10;
    size_t current_project_version_ = 0;
    
    ProgressCallback progress_callback_;
    
    std::string generate_version_path(uint64_t version_id) const;
    bool save_version_metadata(const ProjectVersion& version);
    bool load_version_metadata();
    std::string calculate_project_checksum(const ProjectManager& project);
};

class VersionManagerSingleton {
public:
    static VersionManager& getInstance();
    
    VersionManagerSingleton(const VersionManagerSingleton&) = delete;
    VersionManagerSingleton& operator=(const VersionManagerSingleton&) = delete;
    
private:
    VersionManagerSingleton() = default;
    ~VersionManagerSingleton() = default;
};

} // namespace tool
