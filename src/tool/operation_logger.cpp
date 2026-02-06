/**
 * @file operation_logger.cpp
 * @brief 基础工具层 - 操作日志源文件
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include "tool/operation_logger.hpp"
#include <chrono>
#include <sstream>
#include <iomanip>
#include <fstream>

namespace tool {

static std::string log_type_to_string(OperationType type) {
    switch (type) {
        case OperationType::CREATE: return "CREATE";
        case OperationType::READ: return "READ";
        case OperationType::WRITE: return "WRITE";
        case OperationType::MODIFY: return "MODIFY";
        case OperationType::DELETE: return "DELETE";
        case OperationType::IMPORT: return "IMPORT";
        case OperationType::EXPORT: return "EXPORT";
        case OperationType::VALIDATE: return "VALIDATE";
        case OperationType::FORMAT_CONVERT: return "FORMAT_CONVERT";
        case OperationType::VERSION_CREATE: return "VERSION_CREATE";
        case OperationType::VERSION_ROLLBACK: return "VERSION_ROLLBACK";
        default: return "UNKNOWN";
    }
}

static std::string log_status_to_string(OperationStatus status) {
    switch (status) {
        case OperationStatus::SUCCESS: return "SUCCESS";
        case OperationStatus::FAILED: return "FAILED";
        case OperationStatus::CANCELLED: return "CANCELLED";
        case OperationStatus::IN_PROGRESS: return "IN_PROGRESS";
        default: return "UNKNOWN";
    }
}

std::string OperationLog::to_string() const {
    std::ostringstream oss;
    oss << "[" << log_id << "] ";
    oss << log_type_to_string(type) << " ";
    oss << log_status_to_string(status) << " ";
    
    auto time_t_now = std::chrono::system_clock::to_time_t(timestamp);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        timestamp.time_since_epoch()) % 1000;
    oss << std::put_time(std::localtime(&time_t_now), "%Y-%m-%d %H:%M:%S");
    oss << "." << std::setfill('0') << std::setw(3) << ms.count() << " ";
    
    oss << target_type << ":" << target_id << " ";
    oss << description;
    return oss.str();
}

OperationLogger::OperationLogger() = default;

uint64_t OperationLogger::log_operation(OperationType type, const std::string& description,
                                     const std::string& target_type,
                                     const std::string& target_id) {
    OperationLog log;
    log.log_id = next_log_id_++;
    log.type = type;
    log.status = OperationStatus::SUCCESS;
    log.operator_name = operator_name_;
    log.timestamp = std::chrono::system_clock::now();
    log.target_type = target_type;
    log.target_id = target_id;
    log.description = description;
    log.session_id = session_id_;
    log.project_name = project_name_;
    
    logs_.push_back(log);
    
    return log.log_id;
}

uint64_t OperationLogger::log_create(const std::string& target_type, const std::string& target_id,
                                   const std::string& description) {
    return log_operation(OperationType::CREATE, description, target_type, target_id);
}

uint64_t OperationLogger::log_read(const std::string& target_type, const std::string& target_id,
                                 const std::string& description) {
    return log_operation(OperationType::READ, description, target_type, target_id);
}

uint64_t OperationLogger::log_write(const std::string& target_type, const std::string& target_id,
                                  const std::string& description) {
    return log_operation(OperationType::WRITE, description, target_type, target_id);
}

uint64_t OperationLogger::log_modify(const std::string& target_type, const std::string& target_id,
                                    const std::string& description, const std::string& before_value,
                                    const std::string& after_value) {
    OperationLog log;
    log.log_id = next_log_id_++;
    log.type = OperationType::MODIFY;
    log.status = OperationStatus::SUCCESS;
    log.operator_name = operator_name_;
    log.timestamp = std::chrono::system_clock::now();
    log.target_type = target_type;
    log.target_id = target_id;
    log.description = description;
    log.before_value = before_value;
    log.after_value = after_value;
    log.session_id = session_id_;
    log.project_name = project_name_;
    
    logs_.push_back(log);
    
    return log.log_id;
}

uint64_t OperationLogger::log_delete(const std::string& target_type, const std::string& target_id,
                                   const std::string& description) {
    return log_operation(OperationType::DELETE, description, target_type, target_id);
}

uint64_t OperationLogger::log_import(const std::string& source_path, const std::string& target_type,
                                    const std::string& description) {
    OperationLog log;
    log.log_id = next_log_id_++;
    log.type = OperationType::IMPORT;
    log.status = OperationStatus::SUCCESS;
    log.operator_name = operator_name_;
    log.timestamp = std::chrono::system_clock::now();
    log.target_type = target_type;
    log.target_id = source_path;
    log.description = description;
    log.additional_info = "Source: " + source_path;
    log.session_id = session_id_;
    log.project_name = project_name_;
    
    logs_.push_back(log);
    
    return log.log_id;
}

uint64_t OperationLogger::log_export(const std::string& target_type, const std::string& target_id,
                                   const std::string& destination_path, const std::string& description) {
    OperationLog log;
    log.log_id = next_log_id_++;
    log.type = OperationType::EXPORT;
    log.status = OperationStatus::SUCCESS;
    log.operator_name = operator_name_;
    log.timestamp = std::chrono::system_clock::now();
    log.target_type = target_type;
    log.target_id = target_id;
    log.description = description;
    log.additional_info = "Destination: " + destination_path;
    log.session_id = session_id_;
    log.project_name = project_name_;
    
    logs_.push_back(log);
    
    return log.log_id;
}

uint64_t OperationLogger::log_validate(const std::string& target_type, const std::string& target_id,
                                     const std::string& description) {
    return log_operation(OperationType::VALIDATE, description, target_type, target_id);
}

uint64_t OperationLogger::log_format_convert(const std::string& source_format, const std::string& target_format,
                                           const std::string& description) {
    OperationLog log;
    log.log_id = next_log_id_++;
    log.type = OperationType::FORMAT_CONVERT;
    log.status = OperationStatus::SUCCESS;
    log.operator_name = operator_name_;
    log.timestamp = std::chrono::system_clock::now();
    log.target_type = "Format";
    log.target_id = source_format + " -> " + target_format;
    log.description = description;
    log.additional_info = "Source: " + source_format + ", Target: " + target_format;
    log.session_id = session_id_;
    log.project_name = project_name_;
    
    logs_.push_back(log);
    
    return log.log_id;
}

bool OperationLogger::update_operation_status(uint64_t log_id, OperationStatus status,
                                           const std::string& error_message) {
    for (auto& log : logs_) {
        if (log.log_id == log_id) {
            log.status = status;
            if (!error_message.empty()) {
                log.additional_info = "Error: " + error_message;
            }
            return true;
        }
    }
    return false;
}

std::vector<OperationLog> OperationLogger::get_all_logs() const {
    return logs_;
}

std::vector<OperationLog> OperationLogger::get_logs_by_type(OperationType type) const {
    std::vector<OperationLog> result;
    for (const auto& log : logs_) {
        if (log.type == type) {
            result.push_back(log);
        }
    }
    return result;
}

std::vector<OperationLog> OperationLogger::get_logs_by_time_range(
    const std::chrono::system_clock::time_point& start,
    const std::chrono::system_clock::time_point& end) const {
    
    std::vector<OperationLog> result;
    for (const auto& log : logs_) {
        if (log.timestamp >= start && log.timestamp <= end) {
            result.push_back(log);
        }
    }
    return result;
}

std::vector<OperationLog> OperationLogger::get_logs_by_target(const std::string& target_type,
                                                             const std::string& target_id) const {
    std::vector<OperationLog> result;
    for (const auto& log : logs_) {
        if (log.target_type == target_type && log.target_id == target_id) {
            result.push_back(log);
        }
    }
    return result;
}

std::optional<OperationLog> OperationLogger::get_log(uint64_t log_id) const {
    for (const auto& log : logs_) {
        if (log.log_id == log_id) {
            return log;
        }
    }
    return std::nullopt;
}

bool OperationLogger::export_to_file(const std::string& file_path) const {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return false;
    }
    
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<OperationLogs>\n";
    file << "  <Count>" << logs_.size() << "</Count>\n";
    
    for (const auto& log : logs_) {
        auto time_t_now = std::chrono::system_clock::to_time_t(log.timestamp);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            log.timestamp.time_since_epoch()) % 1000;
        
        file << "  <Log id=\"" << log.log_id << "\">\n";
        file << "    <Type>" << log_type_to_string(log.type) << "</Type>\n";
        file << "    <Status>" << log_status_to_string(log.status) << "</Status>\n";
        file << "    <Timestamp>" << std::put_time(std::localtime(&time_t_now), "%Y-%m-%d %H:%M:%S");
        file << "." << std::setfill('0') << std::setw(3) << ms.count() << "</Timestamp>\n";
        file << "    <Operator>" << log.operator_name << "</Operator>\n";
        file << "    <TargetType>" << log.target_type << "</TargetType>\n";
        file << "    <TargetID>" << log.target_id << "</TargetID>\n";
        file << "    <Description>" << log.description << "</Description>\n";
        if (!log.before_value.empty()) {
            file << "    <BeforeValue>" << log.before_value << "</BeforeValue>\n";
        }
        if (!log.after_value.empty()) {
            file << "    <AfterValue>" << log.after_value << "</AfterValue>\n";
        }
        file << "  </Log>\n";
    }
    
    file << "</OperationLogs>\n";
    
    file.close();
    return true;
}

bool OperationLogger::clear_logs() {
    logs_.clear();
    next_log_id_ = 1;
    return true;
}

std::string OperationLogger::operation_type_to_string(OperationType type) {
    return log_type_to_string(type);
}

std::string OperationLogger::operation_status_to_string(OperationStatus status) {
    return log_status_to_string(status);
}

OperationType OperationLogger::string_to_operation_type(const std::string& str) {
    if (str == "CREATE") return OperationType::CREATE;
    if (str == "READ") return OperationType::READ;
    if (str == "WRITE") return OperationType::WRITE;
    if (str == "MODIFY") return OperationType::MODIFY;
    if (str == "DELETE") return OperationType::DELETE;
    if (str == "IMPORT") return OperationType::IMPORT;
    if (str == "EXPORT") return OperationType::EXPORT;
    if (str == "VALIDATE") return OperationType::VALIDATE;
    if (str == "FORMAT_CONVERT") return OperationType::FORMAT_CONVERT;
    return OperationType::UNKNOWN;
}

OperationLogger& OperationLoggerSingleton::getInstance() {
    static OperationLogger instance;
    return instance;
}

} // namespace tool
