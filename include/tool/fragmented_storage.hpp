/**
 * @file fragmented_storage.hpp
 * @brief 基础工具层 - 分片存储头文件
 * @details 实现大体积数据的分片存储与按需加载
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>
#include <cstdint>
#include <optional>

namespace tool {

constexpr size_t DEFAULT_FRAGMENT_SIZE = 64 * 1024;
constexpr size_t MIN_FRAGMENT_SIZE = 1024;
constexpr size_t MAX_FRAGMENT_SIZE = 1024 * 1024;

struct FragmentInfo {
    std::string fragment_id;
    size_t fragment_index;
    size_t offset;
    size_t size;
    uint64_t checksum;
    std::string data_type;
    std::string description;
};

struct StorageHeader {
    uint32_t magic_number;
    uint32_t version;
    uint64_t total_size;
    size_t fragment_count;
    size_t fragment_size;
    uint64_t creation_time;
    uint64_t modification_time;
    std::string data_type;
    std::string description;
};

class FragmentedStorage {
public:
    explicit FragmentedStorage(size_t fragment_size = DEFAULT_FRAGMENT_SIZE);
    ~FragmentedStorage();
    
    bool create(const std::string& file_path);
    bool open(const std::string& file_path, bool read_only = false);
    bool close();
    
    bool write_fragment(const std::string& fragment_id, const void* data, size_t size);
    bool read_fragment(const std::string& fragment_id, void* buffer, size_t buffer_size);
    std::vector<uint8_t> read_fragment(const std::string& fragment_id);
    
    bool append_data(const void* data, size_t size);
    bool append_string(const std::string& str);
    
    std::vector<FragmentInfo> get_all_fragments() const;
    std::optional<FragmentInfo> get_fragment_info(const std::string& fragment_id) const;
    
    bool has_fragment(const std::string& fragment_id) const;
    bool remove_fragment(const std::string& fragment_id);
    bool remove_fragments(const std::vector<std::string>& fragment_ids);
    
    size_t get_fragment_size() const { return fragment_size_; }
    size_t get_fragment_count() const { return fragments_.size(); }
    uint64_t get_total_size() const { return header_.total_size; }
    
    static uint64_t calculate_checksum(const void* data, size_t size);
    
    bool is_open() const { return is_open_; }
    const std::string& get_file_path() const { return file_path_; }
    
    void set_fragment_size(size_t size);
    void set_data_type(const std::string& type) { header_.data_type = type; }
    void set_description(const std::string& desc) { header_.description = desc; }
    
    using ProgressCallback = std::function<void(size_t current, size_t total)>;
    void set_read_progress_callback(ProgressCallback callback) { read_progress_callback_ = callback; }
    void set_write_progress_callback(ProgressCallback callback) { write_progress_callback_ = callback; }
    
private:
    bool write_header();
    bool read_header();
    bool write_fragment_table();
    bool read_fragment_table();
    
    std::string file_path_;
    FILE* file_ = nullptr;
    bool is_open_ = false;
    bool read_only_ = false;
    size_t fragment_size_ = DEFAULT_FRAGMENT_SIZE;
    
    StorageHeader header_;
    std::unordered_map<std::string, FragmentInfo> fragments_;
    std::vector<std::string> fragment_order_;
    
    ProgressCallback read_progress_callback_;
    ProgressCallback write_progress_callback_;
    
    static constexpr uint32_t MAGIC_NUMBER = 0x46475354;
    static constexpr uint32_t VERSION = 1;
};

class FragmentedReader {
public:
    explicit FragmentedReader(const std::string& file_path);
    ~FragmentedReader();
    
    bool open();
    void close();
    
    std::vector<FragmentInfo> list_fragments() const;
    bool read_fragment(const std::string& fragment_id, void* buffer, size_t buffer_size);
    std::vector<uint8_t> read_fragment(const std::string& fragment_id);
    std::string read_string(const std::string& fragment_id);
    
    size_t get_fragment_size() const { return storage_ ? storage_->get_fragment_size() : 0; }
    size_t get_fragment_count() const { return storage_ ? storage_->get_fragment_count() : 0; }
    
    bool is_open() const { return is_open_; }
    
private:
    std::string file_path_;
    std::unique_ptr<FragmentedStorage> storage_;
    bool is_open_ = false;
};

class FragmentedWriter {
public:
    FragmentedWriter(const std::string& file_path, size_t fragment_size = DEFAULT_FRAGMENT_SIZE);
    ~FragmentedWriter();
    
    bool create();
    
    bool write_fragment(const std::string& fragment_id, const void* data, size_t size);
    bool write_string(const std::string& fragment_id, const std::string& str);
    
    void set_data_type(const std::string& type) { if (storage_) storage_->set_data_type(type); }
    void set_description(const std::string& desc) { if (storage_) storage_->set_description(desc); }
    
    bool finalize();
    
    bool is_open() const { return storage_ && storage_->is_open(); }
    
private:
    std::string file_path_;
    size_t fragment_size_;
    std::unique_ptr<FragmentedStorage> storage_;
};

} // namespace tool
