/**
 * @file em_exception.cpp
 * @brief 基础工具层 - 电磁模块异常处理工具源文件
 */

#include "tool/em_exception.hpp"

namespace exception_utils {

std::string getExceptionChain(const std::exception& e) {
    return std::string(e.what());
}

bool isRecoverable(const std::exception& e) {
    std::string msg = e.what();
    if (msg.find("not found") != std::string::npos || 
        msg.find("corrupt") != std::string::npos) {
        return false;
    }
    return true;
}

std::string getRecoverySuggestion(const std::exception& e) {
    std::string msg = e.what();
    if (msg.find("Project file not found") != std::string::npos) {
        return "Please verify the file path is correct.";
    }
    if (msg.find("XML parsing error") != std::string::npos) {
        return "Please check the XML file structure.";
    }
    return "Please check the log file for more details.";
}

} // namespace exception_utils
