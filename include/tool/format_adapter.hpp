/**
 * @file format_adapter.hpp
 * @brief 基础工具层 - 格式适配器头文件
 * @details 提供多格式兼容与无缝转换的适配层
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include "tool/project_data.hpp"
#include "tool/project_manager.hpp"
#include <string>
#include <memory>
#include <functional>

namespace tool {

enum class FormatType {
    UNKNOWN,
    XML,
    JSON,
    YAML,
    EMAT,
    MAXWELL,
    AMAT
};

struct FormatInfo {
    FormatType type = FormatType::UNKNOWN;
    std::string extension;
    std::string mime_type;
    std::string description;
    bool is_binary = false;
    bool support_read = false;
    bool support_write = false;
};

class FormatAdapter {
public:
    virtual ~FormatAdapter() = default;
    
    virtual FormatInfo get_format_info() const = 0;
    virtual bool can_read(const std::string& file_path) const = 0;
    virtual bool can_write(const std::string& file_path) const = 0;
    
    virtual std::shared_ptr<ProjectManager> read(const std::string& file_path) = 0;
    virtual bool write(const std::string& file_path, const ProjectManager& project) = 0;
    
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;
};

class FormatAdapterFactory {
public:
    static FormatAdapterFactory& getInstance();
    
    void register_adapter(FormatType type, std::function<std::shared_ptr<FormatAdapter>()> creator);
    void unregister_adapter(FormatType type);
    
    std::shared_ptr<FormatAdapter> get_adapter(FormatType type);
    std::shared_ptr<FormatAdapter> get_adapter(const std::string& file_path);
    std::shared_ptr<FormatAdapter> get_adapter_by_extension(const std::string& extension);
    
    std::vector<FormatType> get_supported_read_formats() const;
    std::vector<FormatType> get_supported_write_formats() const;
    std::vector<FormatType> get_all_formats() const;
    
    FormatInfo get_format_info(FormatType type) const;
    
private:
    FormatAdapterFactory();
    ~FormatAdapterFactory() = default;
    
    std::unordered_map<FormatType, std::function<std::shared_ptr<FormatAdapter>()>> adapters_;
    std::unordered_map<FormatType, FormatInfo> format_infos_;
};

class XmlFormatAdapter : public FormatAdapter {
public:
    XmlFormatAdapter();
    ~XmlFormatAdapter() override = default;
    
    FormatInfo get_format_info() const override;
    bool can_read(const std::string& file_path) const override;
    bool can_write(const std::string& file_path) const override;
    
    std::shared_ptr<ProjectManager> read(const std::string& file_path) override;
    bool write(const std::string& file_path, const ProjectManager& project) override;
    
    std::string get_last_error() const override { return last_error_; }
    void clear_error() override { last_error_.clear(); }
    
private:
    std::string last_error_;
};

class JsonFormatAdapter : public FormatAdapter {
public:
    JsonFormatAdapter();
    ~JsonFormatAdapter() override = default;
    
    FormatInfo get_format_info() const override;
    bool can_read(const std::string& file_path) const override;
    bool can_write(const std::string& file_path) const override;
    
    std::shared_ptr<ProjectManager> read(const std::string& file_path) override;
    bool write(const std::string& file_path, const ProjectManager& project) override;
    
    std::string get_last_error() const override { return last_error_; }
    void clear_error() override { last_error_.clear(); }
    
private:
    std::string last_error_;
};

class EmatFormatAdapter : public FormatAdapter {
public:
    EmatFormatAdapter();
    ~EmatFormatAdapter() override = default;
    
    FormatInfo get_format_info() const override;
    bool can_read(const std::string& file_path) const override;
    bool can_write(const std::string& file_path) const override;
    
    std::shared_ptr<ProjectManager> read(const std::string& file_path) override;
    bool write(const std::string& file_path, const ProjectManager& project) override;
    
    std::string get_last_error() const override { return last_error_; }
    void clear_error() override { last_error_.clear(); }
    
private:
    std::string last_error_;
};

FormatType detect_format(const std::string& file_path);
std::string format_type_to_string(FormatType type);
FormatType string_to_format_type(const std::string& str);

} // namespace tool
