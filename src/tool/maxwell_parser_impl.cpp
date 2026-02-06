#include "tool/maxwell_parser_impl.hpp"
#include "tool/logger_factory.hpp"
#include <algorithm>
#include <filesystem>

namespace tool {

bool MaxwellParserImpl::canParse(const std::string& file_path) {
    // 检查文件扩展名
    std::filesystem::path path(file_path);
    std::string extension = path.extension().string();
    
    // 支持的文件扩展名
    std::vector<std::string> supported_extensions = {
        ".aedt", ".aedtz", ".amat", ".xml"
    };
    
    bool supported = std::find(supported_extensions.begin(), 
                              supported_extensions.end(), 
                              extension) != supported_extensions.end();
    
    if (supported) {
        FEEM_DEBUG("文件格式支持: {}", file_path);
    } else {
        FEEM_DEBUG("文件格式不支持: {}", file_path);
    }
    
    return supported;
}

MaxwellFileInfo MaxwellParserImpl::parseFileInfo() {
    if (!parser_.parse_file(file_path_)) {
        FEEM_ERROR("文件解析失败: {}", file_path_);
        throw MaxwellException(MaxwellErrorCode::INVALID_FORMAT, 
                              "文件解析失败: " + file_path_);
    }
    
    return extractFileInfo();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseMaterials() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析材料数据");
    return extractMaterials();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseBoundaries() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析边界条件数据");
    return extractBoundaries();
}

std::vector<nlohmann::json> MaxwellParserImpl::parseExcitations() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析激励源数据");
    return extractExcitations();
}

nlohmann::json MaxwellParserImpl::parseSolutionSetup() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析求解设置数据");
    return extractSolutionSetup();
}

nlohmann::json MaxwellParserImpl::parseGeometry() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析几何数据");
    return extractGeometry();
}

nlohmann::json MaxwellParserImpl::parseAllData() {
    if (!isParseTreeValid()) {
        FEEM_ERROR("解析树无效，请先调用parseFileInfo");
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, 
                              "解析树无效，请先调用parseFileInfo");
    }
    
    FEEM_DEBUG("开始解析所有数据");
    
    nlohmann::json result;
    result["file_info"] = extractFileInfo().toJson();
    result["materials"] = extractMaterials();
    result["boundaries"] = extractBoundaries();
    result["excitations"] = extractExcitations();
    result["solution_setup"] = extractSolutionSetup();
    result["geometry"] = extractGeometry();
    
    return result;
}

// 私有方法实现

MaxwellFileInfo MaxwellParserImpl::extractFileInfo() const {
    MaxwellFileInfo info;
    
    auto root = parser_.get_root();
    if (!root) {
        throw MaxwellException(MaxwellErrorCode::DATA_CORRUPTED, "根节点为空");
    }
    
    info.file_path = file_path_;
    info.file_format = "Maxwell .aedt";
    
    // 从根块提取基本信息
    if (auto created_prop = root->find_property("Created")) {
        info.created_date = std::get<std::string>(created_prop->value);
    }
    
    if (auto product_prop = root->find_property("Product")) {
        info.maxwell_version = std::get<std::string>(product_prop->value);
    }
    
    // 从Desktop块提取版本信息
    auto desktop_block = findBlock("Desktop", root);
    if (desktop_block) {
        // 可以进一步提取版本信息
    }
    
    // 设置文件大小
    std::filesystem::path path(file_path_);
    if (std::filesystem::exists(path)) {
        info.file_size = std::filesystem::file_size(path);
    }
    
    return info;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractMaterials() const {
    std::vector<nlohmann::json> materials;
    
    auto root = parser_.get_root();
    if (!root) return materials;
    
    // 查找Materials块
    auto materials_block = findBlock("Materials", root);
    if (!materials_block) {
        return materials;
    }
    
    // 提取所有材料子块
    for (const auto& material_block : materials_block->children) {
        materials.push_back(blockToJson(material_block));
    }
    
    return materials;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractBoundaries() const {
    std::vector<nlohmann::json> boundaries;
    
    // TODO: 实现边界条件数据提取
    // 需要根据实际Maxwell文件结构实现
    
    return boundaries;
}

std::vector<nlohmann::json> MaxwellParserImpl::extractExcitations() const {
    std::vector<nlohmann::json> excitations;
    
    // TODO: 实现激励源数据提取
    // 需要根据实际Maxwell文件结构实现
    
    return excitations;
}

nlohmann::json MaxwellParserImpl::extractSolutionSetup() const {
    nlohmann::json setup;
    
    // TODO: 实现求解设置数据提取
    // 需要根据实际Maxwell文件结构实现
    
    return setup;
}

nlohmann::json MaxwellParserImpl::extractGeometry() const {
    nlohmann::json geometry;
    
    // TODO: 实现几何数据提取
    // 需要根据实际Maxwell文件结构实现
    
    return geometry;
}

bool MaxwellParserImpl::isParseTreeValid() const {
    return parser_.get_root() != nullptr && parser_.validate();
}

std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode> 
MaxwellParserImpl::findBlock(const std::string& block_name,
                            const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& start_node) const {
    if (!start_node) return nullptr;
    
    // 广度优先搜索
    if (start_node->name == block_name) {
        return start_node;
    }
    
    for (const auto& child : start_node->children) {
        if (child->name == block_name) {
            return child;
        }
    }
    
    // 递归搜索子块
    for (const auto& child : start_node->children) {
        auto result = findBlock(block_name, child);
        if (result) {
            return result;
        }
    }
    
    return nullptr;
}

nlohmann::json MaxwellParserImpl::propertyToJson(const fe_em::tool::maxwell_parser::Property& prop) const {
    nlohmann::json result;
    
    result["name"] = prop.name;
    result["type"] = static_cast<int>(prop.type);
    result["line_number"] = prop.line_number;
    
    // 根据类型转换值
    std::visit([&](const auto& value) {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, std::string>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, double>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, bool>) {
            result["value"] = value;
        } else if constexpr (std::is_same_v<T, std::vector<fe_em::tool::maxwell_parser::Value>>) {
            nlohmann::json array;
            for (const auto& item : value) {
                // 递归处理数组元素
                if (std::holds_alternative<double>(item)) {
                    array.push_back(std::get<double>(item));
                } else if (std::holds_alternative<std::string>(item)) {
                    array.push_back(std::get<std::string>(item));
                } else if (std::holds_alternative<bool>(item)) {
                    array.push_back(std::get<bool>(item));
                }
            }
            result["value"] = array;
        } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
            result["value"] = value;
        }
    }, prop.value);
    
    return result;
}

nlohmann::json MaxwellParserImpl::blockToJson(const std::shared_ptr<fe_em::tool::maxwell_parser::BlockNode>& block) const {
    nlohmann::json result;
    
    if (!block) return result;
    
    result["name"] = block->name;
    result["start_line"] = block->start_line;
    result["end_line"] = block->end_line;
    
    // 转换属性
    nlohmann::json properties;
    for (const auto& prop : block->properties) {
        properties[prop.name] = propertyToJson(prop);
    }
    result["properties"] = properties;
    
    // 递归转换子块
    nlohmann::json children;
    for (const auto& child : block->children) {
        children.push_back(blockToJson(child));
    }
    result["children"] = children;
    
    return result;
}

} // namespace tool