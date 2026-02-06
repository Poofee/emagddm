/**
 * @file format_adapter.cpp
 * @brief 基础工具层 - 格式适配器源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/format_adapter.hpp"
#include "tool/xml_interface.hpp"
#include "tool/file_utils.hpp"
#include <algorithm>
#include <cctype>

namespace tool {

FormatType detect_format(const std::string& file_path) {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    if (ext == ".xml" || ext == ".maxwell") {
        return FormatType::XML;
    } else if (ext == ".json") {
        return FormatType::JSON;
    } else if (ext == ".yaml" || ext == ".yml") {
        return FormatType::YAML;
    } else if (ext == ".emat") {
        return FormatType::EMAT;
    } else if (ext == ".amat") {
        return FormatType::AMAT;
    }
    
    return FormatType::UNKNOWN;
}

std::string format_type_to_string(FormatType type) {
    switch (type) {
        case FormatType::XML: return "XML";
        case FormatType::JSON: return "JSON";
        case FormatType::YAML: return "YAML";
        case FormatType::EMAT: return "EMAT";
        case FormatType::MAXWELL: return "Maxwell";
        case FormatType::AMAT: return "AMAT";
        default: return "Unknown";
    }
}

FormatType string_to_format_type(const std::string& str) {
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str == "xml" || lower_str == "maxwell") {
        return FormatType::XML;
    } else if (lower_str == "json") {
        return FormatType::JSON;
    } else if (lower_str == "yaml" || lower_str == "yml") {
        return FormatType::YAML;
    } else if (lower_str == "emat") {
        return FormatType::EMAT;
    } else if (lower_str == "amat") {
        return FormatType::AMAT;
    }
    
    return FormatType::UNKNOWN;
}

FormatAdapterFactory& FormatAdapterFactory::getInstance() {
    static FormatAdapterFactory instance;
    return instance;
}

FormatAdapterFactory::FormatAdapterFactory() {
    format_infos_[FormatType::XML] = {
        FormatType::XML,
        ".xml",
        "application/xml",
        "XML格式（通用配置格式）",
        false,
        true,
        true
    };
    
    format_infos_[FormatType::JSON] = {
        FormatType::JSON,
        ".json",
        "application/json",
        "JSON格式（明文格式，方便编辑）",
        false,
        true,
        true
    };
    
    format_infos_[FormatType::YAML] = {
        FormatType::YAML,
        ".yaml",
        "application/x-yaml",
        "YAML格式（人类可读的配置格式）",
        false,
        true,
        true
    };
    
    format_infos_[FormatType::EMAT] = {
        FormatType::EMAT,
        ".emat",
        "application/octet-stream",
        "EMAT格式（项目包格式）",
        false,
        true,
        true
    };
    
    format_infos_[FormatType::AMAT] = {
        FormatType::AMAT,
        ".amat",
        "application/octet-stream",
        "AMAT格式（Maxwell材料库格式）",
        false,
        true,
        false
    };
}

void FormatAdapterFactory::register_adapter(FormatType type, 
                                         std::function<std::shared_ptr<FormatAdapter>()> creator) {
    adapters_[type] = creator;
}

void FormatAdapterFactory::unregister_adapter(FormatType type) {
    adapters_.erase(type);
}

std::shared_ptr<FormatAdapter> FormatAdapterFactory::get_adapter(FormatType type) {
    auto it = adapters_.find(type);
    if (it != adapters_.end()) {
        return it->second();
    }
    return nullptr;
}

std::shared_ptr<FormatAdapter> FormatAdapterFactory::get_adapter(const std::string& file_path) {
    FormatType type = detect_format(file_path);
    return get_adapter(type);
}

std::shared_ptr<FormatAdapter> FormatAdapterFactory::get_adapter_by_extension(const std::string& extension) {
    for (const auto& pair : format_infos_) {
        if (pair.second.extension == extension) {
            return get_adapter(pair.first);
        }
    }
    return nullptr;
}

std::vector<FormatType> FormatAdapterFactory::get_supported_read_formats() const {
    std::vector<FormatType> formats;
    for (const auto& pair : adapters_) {
        auto it = format_infos_.find(pair.first);
        if (it != format_infos_.end() && it->second.support_read) {
            formats.push_back(pair.first);
        }
    }
    return formats;
}

std::vector<FormatType> FormatAdapterFactory::get_supported_write_formats() const {
    std::vector<FormatType> formats;
    for (const auto& pair : adapters_) {
        auto it = format_infos_.find(pair.first);
        if (it != format_infos_.end() && it->second.support_write) {
            formats.push_back(pair.first);
        }
    }
    return formats;
}

std::vector<FormatType> FormatAdapterFactory::get_all_formats() const {
    std::vector<FormatType> formats;
    for (const auto& pair : format_infos_) {
        formats.push_back(pair.first);
    }
    return formats;
}

FormatInfo FormatAdapterFactory::get_format_info(FormatType type) const {
    auto it = format_infos_.find(type);
    if (it != format_infos_.end()) {
        return it->second;
    }
    return FormatInfo();
}

XmlFormatAdapter::XmlFormatAdapter() = default;

FormatInfo XmlFormatAdapter::get_format_info() const {
    FormatInfo info;
    info.type = FormatType::XML;
    info.extension = ".xml";
    info.mime_type = "application/xml";
    info.description = "XML格式（通用配置格式）";
    info.is_binary = false;
    info.support_read = true;
    info.support_write = true;
    return info;
}

bool XmlFormatAdapter::can_read(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".xml" || ext == ".maxwell";
}

bool XmlFormatAdapter::can_write(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".xml";
}

std::shared_ptr<ProjectManager> XmlFormatAdapter::read(const std::string& file_path) {
    clear_error();
    
    auto project = std::make_shared<ProjectManager>();
    
    if (!project->openProject(file_path)) {
        last_error_ = project->getLastError();
        return nullptr;
    }
    
    return project;
}

bool XmlFormatAdapter::write(const std::string& file_path, const ProjectManager& project) {
    clear_error();
    
    if (!project.saveProject(file_path)) {
        last_error_ = project.getLastError();
        return false;
    }
    
    return true;
}

JsonFormatAdapter::JsonFormatAdapter() = default;

FormatInfo JsonFormatAdapter::get_format_info() const {
    FormatInfo info;
    info.type = FormatType::JSON;
    info.extension = ".json";
    info.mime_type = "application/json";
    info.description = "JSON格式（明文格式，方便编辑）";
    info.is_binary = false;
    info.support_read = true;
    info.support_write = true;
    return info;
}

bool JsonFormatAdapter::can_read(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".json";
}

bool JsonFormatAdapter::can_write(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".json";
}

std::shared_ptr<ProjectManager> JsonFormatAdapter::read(const std::string& file_path) {
    clear_error();
    last_error_ = "JSON格式读取尚未实现";
    return nullptr;
}

bool JsonFormatAdapter::write(const std::string& file_path, const ProjectManager& project) {
    clear_error();
    last_error_ = "JSON格式写入尚未实现";
    return false;
}

EmatFormatAdapter::EmatFormatAdapter() = default;

FormatInfo EmatFormatAdapter::get_format_info() const {
    FormatInfo info;
    info.type = FormatType::EMAT;
    info.extension = ".emat";
    info.mime_type = "application/octet-stream";
    info.description = "EMAT格式（项目包格式，支持分片存储）";
    info.is_binary = false;
    info.support_read = true;
    info.support_write = true;
    return info;
}

bool EmatFormatAdapter::can_read(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".emat";
}

bool EmatFormatAdapter::can_write(const std::string& file_path) const {
    std::string ext = file_utils::getExtension(file_path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".emat";
}

std::shared_ptr<ProjectManager> EmatFormatAdapter::read(const std::string& file_path) {
    clear_error();
    last_error_ = "EMAT格式读取尚未实现";
    return nullptr;
}

bool EmatFormatAdapter::write(const std::string& file_path, const ProjectManager& project) {
    clear_error();
    last_error_ = "EMAT格式写入尚未实现";
    return false;
}

} // namespace tool
