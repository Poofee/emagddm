#include "tool/maxwell_parser.hpp"
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>
#include <iostream>
#include <functional>

namespace fe_em {
namespace tool {
namespace maxwell_parser {

// 正则表达式模式定义
const std::regex BLOCK_BEGIN_PATTERN(R"(\s*\$begin\s+'(.*?)'\s*)");
const std::regex BLOCK_END_PATTERN(R"(\s*\$end\s+'(.*?)'\s*)");
const std::regex PROPERTY_PATTERN(R"(\s*(\w+)\s*=\s*(.*)\s*)");
const std::regex FUNCTION_PATTERN(R"(\s*(\w+)\s*\((.*?)\)\s*)");
const std::regex ARRAY_PATTERN(R"(\s*\[(\d+):\s*(.*)\]\s*)");
const std::regex SET_PATTERN(R"(\s*set\((.*)\)\s*)");
const std::regex STRING_PATTERN(R"('([^']*)'|\"([^\"]*)\")");
const std::regex NUMBER_PATTERN(R"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)");
const std::regex BOOLEAN_PATTERN(R"(true|false)");

MaxwellParser::MaxwellParser() 
    : current_line_(0) {
}

bool MaxwellParser::parse_file(const std::string& file_path) {
    try {
        // 读取文件内容
        std::ifstream file(file_path);
        if (!file.is_open()) {
            throw ParseError("无法打开文件: " + file_path, 0);
        }
        
        std::stringstream buffer;
        buffer << file.rdbuf();
        file_content_ = buffer.str();
        
        return parse_content(file_content_);
    } catch (const ParseError& e) {
        std::cerr << "解析错误: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "文件读取错误: " << e.what() << std::endl;
        return false;
    }
}

bool MaxwellParser::parse_content(const std::string& content) {
    clear();
    file_content_ = content;
    
    try {
        // 预处理内容
        preprocess_content();
        
        // 开始解析
        current_line_ = 0;
        root_node_ = parse_block();
        
        if (!root_node_) {
            throw ParseError("解析失败，未找到有效的根块", 0);
        }
        
        return validate();
        
    } catch (const ParseError& e) {
        std::cerr << "解析错误: " << e.what() << std::endl;
        return false;
    }
}

void MaxwellParser::preprocess_content() {
    lines_.clear();
    
    std::istringstream stream(file_content_);
    std::string line;
    
    while (std::getline(stream, line)) {
        // 移除行尾的换行符
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        // 跳过空行
        if (line.empty()) {
            continue;
        }
        
        lines_.push_back(line);
    }
}

void MaxwellParser::skip_whitespace_and_comments() {
    while (current_line_ < lines_.size()) {
        std::string line = lines_[current_line_];
        
        // 跳过空白行
        if (line.empty() || std::all_of(line.begin(), line.end(), ::isspace)) {
            current_line_++;
            continue;
        }
        
        // 跳过注释行（以#开头）
        if (line[0] == '#') {
            current_line_++;
            continue;
        }
        
        break;
    }
}

std::shared_ptr<BlockNode> MaxwellParser::parse_block() {
    skip_whitespace_and_comments();
    
    if (current_line_ >= lines_.size()) {
        return nullptr;
    }
    
    std::string line = lines_[current_line_];
    std::smatch match;
    
    // 检查是否是块开始
    if (!std::regex_match(line, match, BLOCK_BEGIN_PATTERN)) {
        throw ParseError("期望块开始标记", current_line_ + 1);
    }
    
    std::string block_name = match[1];
    auto block = std::make_shared<BlockNode>(block_name, current_line_ + 1);
    current_line_++;
    
    // 解析块内容
    while (current_line_ < lines_.size()) {
        skip_whitespace_and_comments();
        
        if (current_line_ >= lines_.size()) {
            throw ParseError("块未正确结束: " + block_name, block->start_line);
        }
        
        line = lines_[current_line_];
        
        // 检查是否是块结束
        if (std::regex_match(line, match, BLOCK_END_PATTERN)) {
            std::string end_block_name = match[1].str();
            if (end_block_name != block_name) {
                throw ParseError("块结束标记不匹配: " + end_block_name + " != " + block_name, 
                                current_line_ + 1);
            }
            block->end_line = current_line_ + 1;
            current_line_++;
            return block;
        }
        
        // 检查是否是子块开始
        if (std::regex_match(line, match, BLOCK_BEGIN_PATTERN)) {
            auto child_block = parse_block();
            if (child_block) {
                block->add_child(child_block);
            }
            continue;
        }
        
        // 解析属性
        try {
            Property prop = parse_property();
            block->add_property(prop);
        } catch (const ParseError& e) {
            // 如果不是属性，可能是其他内容，跳过
            current_line_++;
        }
    }
    
    throw ParseError("块未正确结束: " + block_name, block->start_line);
}

Property MaxwellParser::parse_property() {
    if (current_line_ >= lines_.size()) {
        throw ParseError("期望属性定义", current_line_ + 1);
    }
    
    std::string line = lines_[current_line_];
    std::smatch match;
    
    if (!std::regex_match(line, match, PROPERTY_PATTERN)) {
        throw ParseError("无效的属性格式", current_line_ + 1);
    }
    
    std::string prop_name = match[1];
    std::string value_str = match[2];
    
    // 识别数据类型
    DataType data_type = identify_data_type(value_str);
    Value value;
    
    // 根据数据类型解析值
    switch (data_type) {
        case DataType::STRING:
            value = parse_string(value_str);
            break;
        case DataType::NUMBER:
            value = parse_number(value_str);
            break;
        case DataType::BOOLEAN:
            value = parse_boolean(value_str);
            break;
        case DataType::ARRAY:
            value = parse_array(value_str);
            break;
        case DataType::FUNCTION:
            value = parse_function(value_str);
            break;
        case DataType::SET:
            value = parse_set(value_str);
            break;
        default:
            // 默认为字符串
            value = value_str;
            break;
    }
    
    current_line_++;
    return Property(prop_name, value, data_type, current_line_);
}

DataType MaxwellParser::identify_data_type(const std::string& value_str) {
    std::smatch match;
    
    if (std::regex_match(value_str, match, STRING_PATTERN)) {
        return DataType::STRING;
    } else if (std::regex_match(value_str, match, NUMBER_PATTERN)) {
        return DataType::NUMBER;
    } else if (std::regex_match(value_str, match, BOOLEAN_PATTERN)) {
        return DataType::BOOLEAN;
    } else if (std::regex_match(value_str, match, ARRAY_PATTERN)) {
        return DataType::ARRAY;
    } else if (std::regex_match(value_str, match, FUNCTION_PATTERN)) {
        return DataType::FUNCTION;
    } else if (std::regex_match(value_str, match, SET_PATTERN)) {
        return DataType::SET;
    }
    
    return DataType::UNKNOWN;
}

Value MaxwellParser::parse_string(const std::string& value_str) {
    std::smatch match;
    if (std::regex_search(value_str, match, STRING_PATTERN)) {
        return match[1].matched ? match[1].str() : match[2].str();
    }
    return value_str; // 如果没有引号，返回原始字符串
}

Value MaxwellParser::parse_number(const std::string& value_str) {
    try {
        return std::stod(value_str);
    } catch (const std::exception&) {
        throw ParseError("无效的数值格式: " + value_str, current_line_ + 1);
    }
}

Value MaxwellParser::parse_boolean(const std::string& value_str) {
    return value_str == "true";
}

Value MaxwellParser::parse_array(const std::string& value_str) {
    std::smatch match;
    if (!std::regex_match(value_str, match, ARRAY_PATTERN)) {
        throw ParseError("无效的数组格式: " + value_str, current_line_ + 1);
    }
    
    int size = std::stoi(match[1]);
    std::string data_str = match[2];
    
    std::vector<double> numeric_array;
    std::vector<std::string> string_array;
    std::istringstream stream(data_str);
    std::string item;
    
    // 简单的逗号分隔解析
    while (std::getline(stream, item, ',')) {
        // 去除前后空白
        item.erase(0, item.find_first_not_of(" \t\n\r\f\v"));
        item.erase(item.find_last_not_of(" \t\n\r\f\v") + 1);
        
        if (!item.empty()) {
            DataType item_type = identify_data_type(item);
            switch (item_type) {
                case DataType::NUMBER:
                    numeric_array.push_back(std::get<double>(parse_number(item)));
                    break;
                case DataType::STRING:
                    string_array.push_back(std::get<std::string>(parse_string(item)));
                    break;
                default:
                    string_array.push_back(item);
                    break;
            }
        }
    }
    
    // 优先返回数值数组，如果没有数值则返回字符串数组
    if (!numeric_array.empty()) {
        if (static_cast<int>(numeric_array.size()) != size) {
            throw ParseError("数组大小不匹配: 声明=" + std::to_string(size) + 
                            ", 实际=" + std::to_string(numeric_array.size()), current_line_ + 1);
        }
        return numeric_array;
    } else {
        if (static_cast<int>(string_array.size()) != size) {
            throw ParseError("数组大小不匹配: 声明=" + std::to_string(size) + 
                            ", 实际=" + std::to_string(string_array.size()), current_line_ + 1);
        }
        return string_array;
    }
}

Value MaxwellParser::parse_function(const std::string& value_str) {
    std::smatch match;
    if (!std::regex_match(value_str, match, FUNCTION_PATTERN)) {
        throw ParseError("无效的函数格式: " + value_str, current_line_ + 1);
    }
    
    std::string func_name = match[1];
    std::string args_str = match[2];
    
    // 将函数调用转换为字符串表示
    return func_name + "(" + args_str + ")";
}

Value MaxwellParser::parse_set(const std::string& value_str) {
    std::smatch match;
    if (!std::regex_match(value_str, match, SET_PATTERN)) {
        throw ParseError("无效的集合格式: " + value_str, current_line_ + 1);
    }
    
    std::string items_str = match[1];
    std::vector<std::string> set_items;
    std::istringstream stream(items_str);
    std::string item;
    
    // 解析逗号分隔的集合项
    while (std::getline(stream, item, ',')) {
        // 去除前后空白和引号
        item.erase(0, item.find_first_not_of(" \t\n\r\f\v'\""));
        item.erase(item.find_last_not_of(" \t\n\r\f\v'\"") + 1);
        
        if (!item.empty()) {
            set_items.push_back(item);
        }
    }
    
    return set_items;
}

void MaxwellParser::clear() {
    file_content_.clear();
    root_node_.reset();
    lines_.clear();
    current_line_ = 0;
}

std::string MaxwellParser::get_error_info() const {
    if (!root_node_) {
        return "解析失败: 未生成有效解析树";
    }
    
    // 检查块是否完整结束
    if (root_node_->end_line == -1) {
        return "解析警告: 根块未正确结束";
    }
    
    return "解析成功: 生成有效解析树";
}

bool MaxwellParser::validate() const {
    if (!root_node_) {
        return false;
    }
    
    // 基本的验证逻辑
    // 可以扩展为更复杂的验证规则
    return root_node_->end_line != -1;
}

void MaxwellParser::print_tree(std::ostream& os, int indent) const {
    if (!root_node_) {
        os << "解析树为空" << std::endl;
        return;
    }
    
    std::function<void(const std::shared_ptr<BlockNode>&, int)> print_node;
    print_node = [&](const std::shared_ptr<BlockNode>& node, int level) {
        // 打印缩进
        std::string indent_str(level * 2, ' ');
        
        // 打印块信息
        os << indent_str << "Block: " << node->name 
           << " (lines " << node->start_line << "-" << node->end_line << ")" << std::endl;
        
        // 打印属性
        for (const auto& prop : node->properties) {
            os << indent_str << "  " << prop.name << " = ";
            
            // 根据类型打印值
            std::visit([&](const auto& value) {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, std::string>) {
                    os << "'" << value << "'";
                } else if constexpr (std::is_same_v<T, double>) {
                    os << value;
                } else if constexpr (std::is_same_v<T, bool>) {
                    os << (value ? "true" : "false");
                } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
                    os << "[字符串数组, 大小=" << value.size() << "]";
                } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                    os << "[数值数组, 大小=" << value.size() << "]";
                }
            }, prop.value);
            
            os << " (" << static_cast<int>(prop.type) << ")" << std::endl;
        }
        
        // 递归打印子块
        for (const auto& child : node->children) {
            print_node(child, level + 1);
        }
    };
    
    print_node(root_node_, indent);
}

} // namespace maxwell_parser
} // namespace tool
} // namespace fe_em