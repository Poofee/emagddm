/**
 * @file operation_logger.hpp
 * @brief 基础工具层 - 操作日志头文件
 * @details 记录所有数据操作的完整日志
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#pragma once

#include <string>
#include <vector>
#include <chrono>
#include <memory>
#include <functional>
#include <optional>

namespace tool {

enum class OperationType {
    CREATE,
    READ,
    WRITE,
    MODIFY,
    DELETE,
    IMPORT,
    EXPORT,
    VALIDATE,
    FORMAT_CONVERT,
    VERSION_CREATE,
    VERSION_ROLLBACK,
    UNKNOWN
};

enum class OperationStatus {
    SUCCESS,
    FAILED,
    CANCELLED,
    IN_PROGRESS,
    UNKNOWN
};

struct OperationLog {
    uint64_t log_id;
    OperationType type;
    OperationStatus status;
    std::string operator_name;
    std::chrono::system_clock::time_point timestamp;
    std::string target_type;
    std::string target_id;
    std::string description;
    std::string before_value;
    std::string after_value;
    std::string ip_address;
    std::string session_id;
    std::string project_name;
    std::string additional_info;
    
    std::string to_string() const;
};

class OperationLogger {
public:
    OperationLogger();
    ~OperationLogger() = default;
    
    uint64_t log_operation(OperationType type, const std::string& description,
                          const std::string& target_type = "", 
                          const std::string& target_id = "");
    
    uint64_t log_create(const std::string& target_type, const std::string& target_id,
                       const std::string& description);
    uint64_t log_read(const std::string& target_type, const std::string& target_id,
                     const std::string& description);
    uint64_t log_write(const std::string& target_type, const std::string& target_id,
                      const std::string& description);
    uint64_t log_modify(const std::string& target_type, const std::string& target_id,
                       const std::string& description, const std::string& before_value = "",
                       const std::string& after_value = "");
    uint64_t log_delete(const std::string& target_type, const std::string& target_id,
                       const std::string& description);
    uint64_t log_import(const std::string& source_path, const std::string& target_type,
                        const std::string& description);
    uint64_t log_export(const std::string& target_type, const std::string& target_id,
                       const std::string& destination_path, const std::string& description);
    uint64_t log_validate(const std::string& target_type, const std::string& target_id,
                         const std::string& description);
    uint64_t log_format_convert(const std::string& source_format, const std::string& target_format,
                                const std::string& description);
    
    bool update_operation_status(uint64_t log_id, OperationStatus status, 
                               const std::string& error_message = "");
    
    std::vector<OperationLog> get_all_logs() const;
    std::vector<OperationLog> get_logs_by_type(OperationType type) const;
    std::vector<OperationLog> get_logs_by_time_range(
        const std::chrono::system_clock::time_point& start,
        const std::chrono::system_clock::time_point& end) const;
    std::vector<OperationLog> get_logs_by_target(const std::string& target_type,
                                                  const std::string& target_id) const;
    
    std::optional<OperationLog> get_log(uint64_t log_id) const;
    
    bool export_to_file(const std::string& file_path) const;
    bool clear_logs();
    
    void set_operator_name(const std::string& name) { operator_name_ = name; }
    void set_session_id(const std::string& session_id) { session_id_ = session_id; }
    void set_project_name(const std::string& project_name) { project_name_ = project_name; }
    
    static std::string operation_type_to_string(OperationType type);
    static std::string operation_status_to_string(OperationStatus status);
    static OperationType string_to_operation_type(const std::string& str);
    
private:
    std::vector<OperationLog> logs_;
    uint64_t next_log_id_ = 1;
    std::string operator_name_ = "Unknown";
    std::string session_id_ = "";
    std::string project_name_ = "";
};

class OperationLoggerSingleton {
public:
    static OperationLogger& getInstance();
    
    OperationLoggerSingleton(const OperationLoggerSingleton&) = delete;
    OperationLoggerSingleton& operator=(const OperationLoggerSingleton&) = delete;
    
private:
    OperationLoggerSingleton() = default;
    ~OperationLoggerSingleton() = default;
};

} // namespace tool
