/**
 * @file id_generator.cpp
 * @brief 基础工具层 - 全局唯一ID生成工具源文件
 */

#include "tool/id_generator.hpp"
#include <stdexcept>

namespace tool {

IDGenerator::IDGenerator() {
    category_counters_[static_cast<int>(IDCategory::MATERIAL)] = 1000;
    category_counters_[static_cast<int>(IDCategory::GEOMETRY)] = 10000;
    category_counters_[static_cast<int>(IDCategory::BOUNDARY)] = 20000;
    category_counters_[static_cast<int>(IDCategory::EXCITATION)] = 30000;
}

IDGenerator& IDGenerator::getInstance() {
    static IDGenerator instance;
    return instance;
}

uint64_t IDGenerator::generateID(IDCategory category) {
    std::lock_guard<std::mutex> lock(mutex_);
    int key = static_cast<int>(category);
    uint64_t id = ++category_counters_[key];
    return id;
}

std::string IDGenerator::generateIDString(IDCategory category, uint64_t id) {
    static const char* prefixes[] = {
        "Mat", "Geom", "Bnd", "Exc", "Mesh", "Unit", "DOF",
        "Proj", "Result", "UDP", "Wire", "Coil", "Winding", "Scale", "Param", "Var"
    };
    int idx = static_cast<int>(category);
    if (idx < 0 || idx >= 16) idx = 0;
    return std::string(prefixes[idx]) + "_" + std::to_string(id);
}

void IDGenerator::resetAll() {
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto& pair : category_counters_) {
        pair.second = 0;
    }
}

EntityIDGenerator::EntityIDGenerator() : id_gen_(IDGenerator::getInstance()) {
}

EntityIDGenerator& EntityIDGenerator::getInstance() {
    static EntityIDGenerator instance;
    return instance;
}

uint64_t EntityIDGenerator::generateMaterialID(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = name_map_.find(name);
    if (it != name_map_.end()) return it->second;
    uint64_t id = id_gen_.generateID(IDCategory::MATERIAL);
    name_map_[name] = id;
    return id;
}

uint64_t EntityIDGenerator::generateBoundaryID(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = name_map_.find(name);
    if (it != name_map_.end()) return it->second;
    uint64_t id = id_gen_.generateID(IDCategory::BOUNDARY);
    name_map_[name] = id;
    return id;
}

uint64_t EntityIDGenerator::generateExcitationID(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = name_map_.find(name);
    if (it != name_map_.end()) return it->second;
    uint64_t id = id_gen_.generateID(IDCategory::EXCITATION);
    name_map_[name] = id;
    return id;
}

uint64_t EntityIDGenerator::generateEntityID(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = name_map_.find(name);
    if (it != name_map_.end()) return it->second;
    uint64_t id = id_gen_.generateID(IDCategory::GEOMETRY);
    name_map_[name] = id;
    return id;
}

void EntityIDGenerator::reset() {
    std::lock_guard<std::mutex> lock(mutex_);
    id_gen_.resetAll();
    name_map_.clear();
}

} // namespace tool
