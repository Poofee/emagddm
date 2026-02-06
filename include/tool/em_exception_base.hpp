/**
 * @file em_exception_base.hpp
 * @brief 基础工具层 - 电磁模块异常基类头文件
 */

#pragma once

#include <string>
#include <stdexcept>

namespace tool {

class EmException : public std::runtime_error {
public:
    explicit EmException(const std::string& message)
        : std::runtime_error(message), error_code_(0), module_name_("Unknown") {}

    EmException(const std::string& module_name, const std::string& message, int error_code)
        : std::runtime_error("[" + module_name + "] " + message), 
          error_code_(error_code), 
          module_name_(module_name) {}

    virtual ~EmException() = default;

    virtual int getErrorCode() const { return error_code_; }
    virtual std::string getModuleName() const { return module_name_; }
    virtual std::string getDetail() const { return what(); }
    virtual std::string getSuggestion() const { return ""; }

protected:
    int error_code_;
    std::string module_name_;
};

} // namespace tool

namespace project {

class ProjectException : public tool::EmException {
public:
    explicit ProjectException(const std::string& message)
        : tool::EmException("Project", message, 1000), suggestion_("") {}

    ProjectException(const std::string& message, int error_code)
        : tool::EmException("Project", message, error_code), suggestion_("") {}

    std::string getSuggestion() const override { return suggestion_; }

protected:
    std::string suggestion_;
};

} // namespace project

namespace data {

class DataReadException : public tool::EmException {
public:
    explicit DataReadException(const std::string& message)
        : tool::EmException("DataRead", message, 2000), node_path_("") {}

    std::string getNodePath() const { return node_path_; }

protected:
    std::string node_path_;
};

class DataValidationException : public tool::EmException {
public:
    explicit DataValidationException(const std::string& message)
        : tool::EmException("DataValidation", message, 4000), field_name_("") {}

    std::string getFieldName() const { return field_name_; }

protected:
    std::string field_name_;
};

} // namespace data

namespace material {

class MaterialException : public tool::EmException {
public:
    explicit MaterialException(const std::string& message)
        : tool::EmException("Material", message, 7000), material_name_("") {}

protected:
    std::string material_name_;
};

} // namespace material

namespace boundary {

class BoundaryException : public tool::EmException {
public:
    explicit BoundaryException(const std::string& message)
        : tool::EmException("Boundary", message, 8000), boundary_name_("") {}

protected:
    std::string boundary_name_;
};

} // namespace boundary

namespace excitation {

class ExcitationException : public tool::EmException {
public:
    explicit ExcitationException(const std::string& message)
        : tool::EmException("Excitation", message, 9000), excitation_name_("") {}

protected:
    std::string excitation_name_;
};

} // namespace excitation

namespace hpc {

class HPCConfigException : public tool::EmException {
public:
    explicit HPCConfigException(const std::string& message)
        : tool::EmException("HPC", message, 10000), config_param_("") {}

protected:
    std::string config_param_;
};

} // namespace hpc

namespace resource {

class ResourceException : public tool::EmException {
public:
    explicit ResourceException(const std::string& message)
        : tool::EmException("Resource", message, 11000), resource_type_(""), file_path_("") {}

protected:
    std::string resource_type_;
    std::string file_path_;
};

} // namespace resource

namespace geometry {

class GeometryException : public tool::EmException {
public:
    explicit GeometryException(const std::string& message)
        : tool::EmException("Geometry", message, 12000), entity_name_("") {}

protected:
    std::string entity_name_;
};

} // namespace geometry
