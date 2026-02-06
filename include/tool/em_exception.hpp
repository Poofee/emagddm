/**
 * @file em_exception.hpp
 * @brief 基础工具层 - 电磁模块异常处理工具头文件（简化版）
 */

#pragma once

#include <stdexcept>
#include <string>

namespace tool {

class EmException : public std::runtime_error {
public:
    explicit EmException(const std::string& msg) : std::runtime_error(msg), code_(0), module_("Unknown") {}
    EmException(const std::string& module, const std::string& msg, int code) 
        : std::runtime_error("[" + module + "] " + msg), code_(code), module_(module) {}
    virtual ~EmException() = default;
    virtual int getErrorCode() const { return code_; }
    virtual std::string getModuleName() const { return module_; }
protected:
    int code_;
    std::string module_;
};

namespace project {

class ProjectException : public EmException {
public:
    explicit ProjectException(const std::string& msg) : EmException("Project", msg, 1000) {}
};

class ProjectNotFoundException : public ProjectException {
public:
    explicit ProjectNotFoundException(const std::string& file_path) 
        : ProjectException("Project file not found: " + file_path), file_path_(file_path) {}
    std::string getFilePath() const { return file_path_; }
private:
    std::string file_path_;
};

} // namespace project

namespace data {

class DataReadException : public EmException {
public:
    explicit DataReadException(const std::string& msg) : EmException("DataRead", msg, 2000) {}
};

class DataValidationException : public EmException {
public:
    explicit DataValidationException(const std::string& msg) : EmException("DataValidation", msg, 4000) {}
};

} // namespace data

namespace format {

class XMLParseException : public data::DataReadException {
public:
    XMLParseException(int line, int col, const std::string& msg) 
        : data::DataReadException("XML parse error at line " + std::to_string(line) + 
                                   ", col " + std::to_string(col) + ": " + msg),
          line_(line), col_(col) {}
    int getLineNumber() const { return line_; }
    int getColumnNumber() const { return col_; }
private:
    int line_, col_;
};

} // namespace format

namespace material {

class MaterialException : public EmException {
public:
    explicit MaterialException(const std::string& msg) : EmException("Material", msg, 7000) {}
};

class BHCurveException : public MaterialException {
public:
    BHCurveException(const std::string& mat, const std::string& msg) 
        : MaterialException(mat + ": BH curve error: " + msg) {}
};

class CorelossFileException : public data::DataReadException {
public:
    explicit CorelossFileException(const std::string& file, const std::string& msg) 
        : data::DataReadException("coreloss_user.data parse error in " + file + ": " + msg) {}
};

} // namespace material

namespace boundary {

class BoundaryException : public EmException {
public:
    explicit BoundaryException(const std::string& msg) : EmException("Boundary", msg, 8000) {}
};

class BoundaryConflictException : public BoundaryException {
public:
    BoundaryConflictException(const std::string& entity, const std::string& b1, const std::string& b2)
        : BoundaryException("Boundary conflict on " + entity + ": " + b1 + " vs " + b2) {}
};

} // namespace boundary

namespace excitation {

class ExcitationException : public EmException {
public:
    explicit ExcitationException(const std::string& msg) : EmException("Excitation", msg, 9000) {}
};

class WaveformException : public ExcitationException {
public:
    WaveformException(const std::string& exc, const std::string& msg)
        : ExcitationException(exc + ": Waveform error: " + msg) {}
};

} // namespace excitation

namespace hpc {

class HPCConfigException : public EmException {
public:
    explicit HPCConfigException(const std::string& msg) : EmException("HPC", msg, 10000) {}
};

} // namespace hpc

namespace geometry {

class GeometryException : public EmException {
public:
    explicit GeometryException(const std::string& msg) : EmException("Geometry", msg, 12000) {}
};

} // namespace geometry

} // namespace tool
