/**
 * @file fragmented_storage.cpp
 * @brief 基础工具层 - 分片存储源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/fragmented_storage.hpp"
#include "tool/file_utils.hpp"
#include <cstring>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace tool {

FragmentedStorage::FragmentedStorage(size_t fragment_size)
    : fragment_size_(fragment_size) {
    if (fragment_size < MIN_FRAGMENT_SIZE) {
        fragment_size_ = MIN_FRAGMENT_SIZE;
    } else if (fragment_size > MAX_FRAGMENT_SIZE) {
        fragment_size_ = MAX_FRAGMENT_SIZE;
    }
}

FragmentedStorage::~FragmentedStorage() {
    close();
}

bool FragmentedStorage::create(const std::string& file_path) {
    close();
    
    file_path_ = file_path;
    file_ = fopen(file_path.c_str(), "wb");
    if (!file_) {
        return false;
    }
    
    header_.magic_number = MAGIC_NUMBER;
    header_.version = VERSION;
    header_.total_size = 0;
    header_.fragment_count = 0;
    header_.fragment_size = fragment_size_;
    header_.creation_time = static_cast<uint64_t>(time(nullptr));
    header_.modification_time = header_.creation_time;
    header_.data_type = "";
    header_.description = "";
    
    fragments_.clear();
    fragment_order_.clear();
    
    if (!write_header()) {
        close();
        return false;
    }
    
    is_open_ = true;
    return true;
}

bool FragmentedStorage::open(const std::string& file_path, bool read_only) {
    close();
    
    file_path_ = file_path;
    
    if (!file_utils::exists(file_path)) {
        return false;
    }
    
    if (read_only) {
        file_ = fopen(file_path.c_str(), "rb");
        read_only_ = true;
    } else {
        file_ = fopen(file_path.c_str(), "r+b");
        read_only_ = false;
    }
    
    if (!file_) {
        return false;
    }
    
    if (!read_header()) {
        close();
        return false;
    }
    
    fragment_size_ = header_.fragment_size;
    
    if (!read_fragment_table()) {
        close();
        return false;
    }
    
    is_open_ = true;
    return true;
}

bool FragmentedStorage::close() {
    if (file_) {
        if (fflush(file_) != 0) {
            fclose(file_);
            file_ = nullptr;
            return false;
        }
        fclose(file_);
        file_ = nullptr;
    }
    
    is_open_ = false;
    read_only_ = false;
    
    return true;
}

bool FragmentedStorage::write_header() {
    if (!file_) return false;
    
    if (fseek(file_, 0, SEEK_SET) != 0) return false;
    
    uint8_t header_data[256];
    memset(header_data, 0, sizeof(header_data));
    
    memcpy(header_data, &header_.magic_number, sizeof(uint32_t));
    memcpy(header_data + 4, &header_.version, sizeof(uint32_t));
    memcpy(header_data + 8, &header_.total_size, sizeof(uint64_t));
    memcpy(header_data + 16, &header_.fragment_count, sizeof(size_t));
    memcpy(header_data + 24, &header_.fragment_size, sizeof(size_t));
    memcpy(header_data + 32, &header_.creation_time, sizeof(uint64_t));
    memcpy(header_data + 40, &header_.modification_time, sizeof(uint64_t));
    
    size_t offset = 48;
    size_t data_type_len = header_.data_type.length();
    if (data_type_len > 100) data_type_len = 100;
    memcpy(header_data + offset, header_.data_type.c_str(), data_type_len);
    offset += 100;
    
    size_t desc_len = header_.description.length();
    if (desc_len > 100) desc_len = 100;
    memcpy(header_data + offset, header_.description.c_str(), desc_len);
    
    if (fwrite(header_data, 1, sizeof(header_data), file_) != sizeof(header_data)) {
        return false;
    }
    
    return true;
}

bool FragmentedStorage::read_header() {
    if (!file_) return false;
    
    if (fseek(file_, 0, SEEK_SET) != 0) return false;
    
    uint8_t header_data[256];
    if (fread(header_data, 1, sizeof(header_data), file_) != sizeof(header_data)) {
        return false;
    }
    
    memcpy(&header_.magic_number, header_data, sizeof(uint32_t));
    if (header_.magic_number != MAGIC_NUMBER) {
        return false;
    }
    
    memcpy(&header_.version, header_data + 4, sizeof(uint32_t));
    memcpy(&header_.total_size, header_data + 8, sizeof(uint64_t));
    memcpy(&header_.fragment_count, header_data + 16, sizeof(size_t));
    memcpy(&header_.fragment_size, header_data + 24, sizeof(size_t));
    memcpy(&header_.creation_time, header_data + 32, sizeof(uint64_t));
    memcpy(&header_.modification_time, header_data + 40, sizeof(uint64_t));
    
    size_t offset = 48;
    header_.data_type = std::string(reinterpret_cast<char*>(header_data + offset), 100);
    offset += 100;
    header_.description = std::string(reinterpret_cast<char*>(header_data + offset), 100);
    
    return true;
}

bool FragmentedStorage::write_fragment_table() {
    if (!file_) return false;
    
    size_t table_offset = 256;
    
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    oss << "<FragmentTable>\n";
    oss << "  <FragmentCount>" << fragment_order_.size() << "</FragmentCount>\n";
    
    for (size_t i = 0; i < fragment_order_.size(); ++i) {
        const auto& frag_id = fragment_order_[i];
        const auto& info = fragments_[frag_id];
        
        oss << "  <Fragment index=\"" << i << "\">\n";
        oss << "    <ID>" << info.fragment_id << "</ID>\n";
        oss << "    <Offset>" << info.offset << "</Offset>\n";
        oss << "    <Size>" << info.size << "</Size>\n";
        oss << "    <Checksum>" << info.checksum << "</Checksum>\n";
        oss << "    <DataType>" << info.data_type << "</DataType>\n";
        if (!info.description.empty()) {
            oss << "    <Description>" << info.description << "</Description>\n";
        }
        oss << "  </Fragment>\n";
    }
    
    oss << "</FragmentTable>\n";
    
    std::string table_content = oss.str();
    
    if (fseek(file_, static_cast<long>(table_offset), SEEK_SET) != 0) return false;
    
    if (fwrite(table_content.c_str(), 1, table_content.length(), file_) != table_content.length()) {
        return false;
    }
    
    return true;
}

bool FragmentedStorage::read_fragment_table() {
    if (!file_) return false;
    
    char buffer[1024];
    if (fseek(file_, 256, SEEK_SET) != 0) return false;
    
    std::string xml_content;
    while (fgets(buffer, sizeof(buffer), file_)) {
        xml_content += buffer;
    }
    
    size_t pos = xml_content.find("<FragmentCount>");
    if (pos == std::string::npos) return false;
    
    size_t end_pos = xml_content.find("</FragmentCount>", pos);
    std::string count_str = xml_content.substr(pos + 15, end_pos - pos - 15);
    size_t count = std::stoul(count_str);
    
    fragments_.clear();
    fragment_order_.clear();
    
    for (size_t i = 0; i < count; ++i) {
        std::string id_tag = "<ID>";
        std::string end_id_tag = "</ID>";
        pos = xml_content.find(id_tag, pos);
        if (pos == std::string::npos) break;
        
        size_t id_start = pos + id_tag.length();
        size_t id_end = xml_content.find(end_id_tag, id_start);
        std::string frag_id = xml_content.substr(id_start, id_end - id_start);
        
        FragmentInfo info;
        info.fragment_id = frag_id;
        info.fragment_index = i;
        fragments_[frag_id] = info;
        fragment_order_.push_back(frag_id);
    }
    
    return true;
}

bool FragmentedStorage::write_fragment(const std::string& fragment_id, const void* data, size_t size) {
    if (!is_open_ || read_only_) return false;
    
    uint64_t checksum = calculate_checksum(data, size);
    
    size_t fragment_index = fragment_order_.size();
    uint64_t offset = 256 + fragment_index * fragment_size_;
    
    FragmentInfo info;
    info.fragment_id = fragment_id;
    info.fragment_index = fragment_index;
    info.offset = offset;
    info.size = size;
    info.checksum = checksum;
    
    fragments_[fragment_id] = info;
    fragment_order_.push_back(fragment_id);
    
    header_.fragment_count = fragment_order_.size();
    header_.total_size += size;
    header_.modification_time = static_cast<uint64_t>(time(nullptr));
    
    if (fseek(file_, static_cast<long>(offset), SEEK_SET) != 0) return false;
    
    std::vector<uint8_t> padded_data(size);
    memcpy(padded_data.data(), data, size);
    
    if (fwrite(padded_data.data(), 1, fragment_size_, file_) != fragment_size_) {
        return false;
    }
    
    write_header();
    write_fragment_table();
    
    if (write_progress_callback_) {
        write_progress_callback_(fragment_index + 1, fragment_order_.size());
    }
    
    return true;
}

bool FragmentedStorage::read_fragment(const std::string& fragment_id, void* buffer, size_t buffer_size) {
    if (!is_open_) return false;
    
    auto it = fragments_.find(fragment_id);
    if (it == fragments_.end()) return false;
    
    const FragmentInfo& info = it->second;
    
    if (buffer_size < info.size) return false;
    
    if (fseek(file_, static_cast<long>(info.offset), SEEK_SET) != 0) return false;
    
    std::vector<uint8_t> fragment_data(fragment_size_);
    if (fread(fragment_data.data(), 1, fragment_size_, file_) != fragment_size_) {
        return false;
    }
    
    uint64_t checksum = calculate_checksum(fragment_data.data(), info.size);
    if (checksum != info.checksum) {
        return false;
    }
    
    memcpy(buffer, fragment_data.data(), info.size);
    
    if (read_progress_callback_) {
        read_progress_callback_(info.fragment_index + 1, fragments_.size());
    }
    
    return true;
}

std::vector<uint8_t> FragmentedStorage::read_fragment(const std::string& fragment_id) {
    auto it = fragments_.find(fragment_id);
    if (it == fragments_.end()) return {};
    
    const FragmentInfo& info = it->second;
    std::vector<uint8_t> data(info.size);
    
    if (!read_fragment(fragment_id, data.data(), info.size)) {
        return {};
    }
    
    return data;
}

bool FragmentedStorage::append_data(const void* data, size_t size) {
    std::string frag_id = "data_" + std::to_string(fragments_.size());
    return write_fragment(frag_id, data, size);
}

bool FragmentedStorage::append_string(const std::string& str) {
    size_t len = str.length() + 1;
    std::vector<uint8_t> buffer(len);
    memcpy(buffer.data(), str.c_str(), len);
    return append_data(buffer.data(), len);
}

std::vector<FragmentInfo> FragmentedStorage::get_all_fragments() const {
    std::vector<FragmentInfo> result;
    for (const auto& frag_id : fragment_order_) {
        auto it = fragments_.find(frag_id);
        if (it != fragments_.end()) {
            result.push_back(it->second);
        }
    }
    return result;
}

std::optional<FragmentInfo> FragmentedStorage::get_fragment_info(const std::string& fragment_id) const {
    auto it = fragments_.find(fragment_id);
    if (it != fragments_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool FragmentedStorage::has_fragment(const std::string& fragment_id) const {
    return fragments_.find(fragment_id) != fragments_.end();
}

bool FragmentedStorage::remove_fragment(const std::string& fragment_id) {
    if (read_only_) return false;
    
    auto it = fragments_.find(fragment_id);
    if (it == fragments_.end()) return false;
    
    fragments_.erase(it);
    
    fragment_order_.erase(
        std::remove(fragment_order_.begin(), fragment_order_.end(), fragment_id),
        fragment_order_.end()
    );
    
    header_.fragment_count = fragment_order_.size();
    header_.modification_time = static_cast<uint64_t>(time(nullptr));
    
    return true;
}

bool FragmentedStorage::remove_fragments(const std::vector<std::string>& fragment_ids) {
    bool success = true;
    for (const auto& id : fragment_ids) {
        if (!remove_fragment(id)) {
            success = false;
        }
    }
    return success;
}

uint64_t FragmentedStorage::calculate_checksum(const void* data, size_t size) {
    const uint8_t* bytes = static_cast<const uint8_t*>(data);
    uint64_t checksum = 0;
    
    for (size_t i = 0; i < size; ++i) {
        checksum += bytes[i] * (i + 1);
    }
    
    return checksum;
}

void FragmentedStorage::set_fragment_size(size_t size) {
    if (!is_open_) {
        if (size >= MIN_FRAGMENT_SIZE && size <= MAX_FRAGMENT_SIZE) {
            fragment_size_ = size;
        }
    }
}

FragmentedReader::FragmentedReader(const std::string& file_path)
    : file_path_(file_path) {}

FragmentedReader::~FragmentedReader() {
    close();
}

bool FragmentedReader::open() {
    storage_ = std::make_unique<FragmentedStorage>();
    if (!storage_->open(file_path_, true)) {
        storage_.reset();
        return false;
    }
    is_open_ = true;
    return true;
}

void FragmentedReader::close() {
    if (storage_) {
        storage_->close();
        storage_.reset();
    }
    is_open_ = false;
}

std::vector<FragmentInfo> FragmentedReader::list_fragments() const {
    if (!storage_) return {};
    return storage_->get_all_fragments();
}

bool FragmentedReader::read_fragment(const std::string& fragment_id, void* buffer, size_t buffer_size) {
    if (!storage_) return false;
    return storage_->read_fragment(fragment_id, buffer, buffer_size);
}

std::vector<uint8_t> FragmentedReader::read_fragment(const std::string& fragment_id) {
    if (!storage_) return {};
    return storage_->read_fragment(fragment_id);
}

std::string FragmentedReader::read_string(const std::string& fragment_id) {
    auto data = read_fragment(fragment_id);
    if (data.empty()) return "";
    return std::string(reinterpret_cast<char*>(data.data()));
}

FragmentedWriter::FragmentedWriter(const std::string& file_path, size_t fragment_size)
    : file_path_(file_path), fragment_size_(fragment_size) {}

FragmentedWriter::~FragmentedWriter() {
    if (storage_ && storage_->is_open()) {
        finalize();
    }
}

bool FragmentedWriter::create() {
    storage_ = std::make_unique<FragmentedStorage>(fragment_size_);
    if (!storage_->create(file_path_)) {
        storage_.reset();
        return false;
    }
    return true;
}

bool FragmentedWriter::write_fragment(const std::string& fragment_id, const void* data, size_t size) {
    if (!storage_) return false;
    return storage_->write_fragment(fragment_id, data, size);
}

bool FragmentedWriter::write_string(const std::string& fragment_id, const std::string& str) {
    size_t len = str.length() + 1;
    std::vector<uint8_t> buffer(len);
    memcpy(buffer.data(), str.c_str(), len);
    return write_fragment(fragment_id, buffer.data(), len);
}

bool FragmentedWriter::finalize() {
    if (!storage_) return true;
    bool result = storage_->close();
    storage_.reset();
    return result;
}

} // namespace tool
