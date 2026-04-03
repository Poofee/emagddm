/**
 * @file project_loader.cpp
 * @brief 项目加载器工厂函数实现
 * @details 根据文件扩展名创建对应格式的加载器实例
 * @author Poofee
 * @date 2026-04-04
 */

#include "project_loader.hpp"
#include "project_manager.hpp"
#include "aedt_project_loader.hpp"
#include <algorithm>

namespace tool {

std::unique_ptr<IProjectLoader> createProjectLoader(const std::string& file_path) {
    if (file_path.empty()) {
        return nullptr;
    }
    
    // 提取文件扩展名
    size_t pos = file_path.find_last_of('.');
    if (pos == std::string::npos) {
        return nullptr;
    }
    
    std::string ext = file_path.substr(pos);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    // 根据扩展名创建对应的加载器
    if (ext == ".aedt") {
        return std::make_unique<AedtProjectLoader>();
    } else if (ext == ".json") {
        // TODO: 实现JsonProjectLoader
        // return std::make_unique<JsonProjectLoader>();
        return nullptr;
    } else if (ext == ".xml") {
        // TODO: 实现XmlProjectLoader
        return nullptr;
    }
    
    return nullptr;
}

std::unique_ptr<IProjectLoader> createProjectLoader(InputFormat format) {
    switch (format) {
        case InputFormat::AEDT:
            return std::make_unique<AedtProjectLoader>();
        case InputFormat::JSON:
            // TODO: 实现JsonProjectLoader
            return nullptr;
        case InputFormat::XML:
            // TODO: 实现XmlProjectLoader
            return nullptr;
        default:
            return nullptr;
    }
}

} // namespace tool
