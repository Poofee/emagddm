/**
 * @file project_loader.hpp
 * @brief 基础工具层 - 项目加载器抽象接口
 * @details 定义统一的项目数据加载规范，支持多种输入格式（AEDT/JSON/XML）
 *          所有格式的加载器实现此接口，将外部文件数据转换为统一的内部数据结构
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include "project_data.hpp"
#include <string>
#include <memory>

namespace tool {

class ProjectManager;

/**
 * @enum InputFormat
 * @brief 支持的输入文件格式枚举
 */
enum class InputFormat {
    AEDT,    ///< ANSYS Maxwell AEDT项目文件
    JSON,    ///< JSON配置文件
    XML,     ///< XML配置文件（预留）
    UNKNOWN  ///< 未知格式
};

/**
 * @class IProjectLoader
 * @brief 项目加载器抽象接口
 * 
 * 所有格式的文件加载器都必须实现此接口。
 * 核心职责：读取外部文件 → 转换为统一数据结构 → 存入 ProjectManager
 * 
 * 使用方式：
 * @code
 * std::unique_ptr<IProjectLoader> loader = createLoader("project.aedt");
 * if (loader->load(project_manager)) {
 *     // 数据已存入 project_manager，可直接用于求解
 * }
 * @endcode
 */
class IProjectLoader {
public:
    virtual ~IProjectLoader() = default;

    /**
     * @brief 设置要加载的文件路径
     * @param file_path 文件路径
     */
    virtual void setFilePath(const std::string& file_path) = 0;

    /**
     * @brief 检查是否支持指定格式的文件
     * @param file_path 文件路径
     * @return bool 支持返回true
     */
    virtual bool canLoad(const std::string& file_path) const = 0;

    /**
     * @brief 从文件加载数据到项目管理器
     * 
     * 这是核心方法，负责：
     * 1. 解析外部文件（AEDT/JSON/XML等）
     * 2. 将解析的数据转换为统一的内部结构（Material/Boundary/Excitation等）
     * 3. 将转换后的数据写入 ProjectManager
     * 
     * @param manager 项目管理器指针，用于存储加载的数据
     * @return bool 加载成功返回true
     */
    virtual bool load(ProjectManager& manager) = 0;

    /**
     * @brief 获取支持的文件格式
     * @return InputFormat 文件格式枚举值
     */
    virtual InputFormat getSupportedFormat() const = 0;

    /**
     * @brief 获取格式描述字符串
     * @return 格式描述（如 "ANSYS Maxwell AEDT"）
     */
    virtual std::string getFormatDescription() const = 0;

    /**
     * @brief 获取最后一次错误信息
     * @return 错误信息字符串
     */
    virtual std::string getLastError() const = 0;
};

/**
 * @brief 工厂函数：根据文件扩展名创建对应格式的加载器
 * 
 * 自动检测文件扩展名并返回对应的加载器实例：
 * - .aedt → AedtProjectLoader
 * - .json → JsonProjectLoader
 * - .xml  → XmlProjectLoader（预留）
 * 
 * @param file_path 输入文件路径
 * @return std::unique_ptr<IProjectLoader> 加载器实例，不支持则返回nullptr
 */
std::unique_ptr<IProjectLoader> createProjectLoader(const std::string& file_path);

/**
 * @brief 工厂函数：根据格式枚举创建对应格式的加载器
 * @param format 目标格式
 * @return std::unique_ptr<IProjectLoader> 加载器实例
 */
std::unique_ptr<IProjectLoader> createProjectLoader(InputFormat format);

} // namespace tool
