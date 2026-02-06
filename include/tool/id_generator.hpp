/**
 * @file id_generator.hpp
 * @brief 基础工具层 - 全局唯一ID生成工具头文件（简化版）
 */

#pragma once
#include "em_enums.hpp"
#include <cstdint>
#include <string>
#include <mutex>
#include <atomic>
#include <unordered_map>

namespace tool {

class IDGenerator {
public:
    static IDGenerator& getInstance();
    uint64_t generateID(IDCategory category);
    std::string generateIDString(IDCategory category, uint64_t id);
    void resetAll();

private:
    IDGenerator();
    ~IDGenerator() = default;
    IDGenerator(const IDGenerator&) = delete;
    IDGenerator& operator=(const IDGenerator&) = delete;

    std::mutex mutex_;
    std::unordered_map<int, std::atomic<uint64_t>> category_counters_;
};

class EntityIDGenerator {
public:
    static EntityIDGenerator& getInstance();
    uint64_t generateMaterialID(const std::string& name);
    uint64_t generateBoundaryID(const std::string& name);
    uint64_t generateExcitationID(const std::string& name);
    uint64_t generateEntityID(const std::string& name);
    void reset();

private:
    EntityIDGenerator();
    ~EntityIDGenerator() = default;
    EntityIDGenerator(const EntityIDGenerator&) = delete;
    EntityIDGenerator& operator=(const EntityIDGenerator&) = delete;

    IDGenerator& id_gen_;
    std::mutex mutex_;
    std::unordered_map<std::string, uint64_t> name_map_;
};

} // namespace tool
