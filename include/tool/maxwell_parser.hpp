#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <variant>
#include <optional>

namespace fe_em {
namespace tool {

/**
 * @brief Maxwell文件解析器命名空间
 * 负责解析Maxwell .aedt项目文件的层次化结构
 */
namespace maxwell_parser {

/**
 * @brief 数据类型枚举
 * 定义Maxwell文件中支持的数据类型
 */
enum class DataType {
    STRING,      // 字符串类型
    NUMBER,      // 数值类型
    BOOLEAN,     // 布尔类型
    ARRAY,       // 数组类型
    FUNCTION,    // 函数调用
    SET,         // 集合类型
    UNKNOWN      // 未知类型
};

/**
 * @brief 数据值类型
 * 支持多种数据类型的存储
 */
using Value = std::variant<
    std::string,                    // 字符串值
    double,                         // 数值
    bool,                           // 布尔值
    std::vector<std::string>,       // 字符串数组
    std::vector<double>             // 数值数组
>;

/**
 * @brief 属性节点
 * 表示文件中的一个属性定义
 */
struct Property {
    std::string name;           // 属性名称
    Value value;                // 属性值
    DataType type;              // 数据类型
    int line_number;            // 行号（用于错误定位）
    
    Property(const std::string& n, const Value& v, DataType t, int ln)
        : name(n), value(v), type(t), line_number(ln) {}
};

/**
 * @brief 块节点
 * 表示文件中的一个$begin/$end块
 */
struct BlockNode : public std::enable_shared_from_this<BlockNode> {
    std::string name;                           // 块名称
    std::vector<Property> properties;          // 属性列表
    std::vector<std::shared_ptr<BlockNode>> children; // 子块列表
    std::weak_ptr<BlockNode> parent;           // 父块引用
    int start_line;                            // 起始行号
    int end_line;                              // 结束行号
    
    BlockNode(const std::string& n, int sl)
        : name(n), start_line(sl), end_line(-1) {}
    
    /**
     * @brief 添加属性到当前块
     */
    void add_property(const Property& prop) {
        properties.push_back(prop);
    }
    
    /**
     * @brief 添加子块到当前块
     */
    void add_child(std::shared_ptr<BlockNode> child) {
        child->parent = shared_from_this();
        children.push_back(child);
    }
    
    /**
     * @brief 查找指定名称的属性
     */
    std::optional<Property> find_property(const std::string& prop_name) const {
        for (const auto& prop : properties) {
            if (prop.name == prop_name) {
                return prop;
            }
        }
        return std::nullopt;
    }
    
    /**
     * @brief 查找指定名称的子块
     */
    std::vector<std::shared_ptr<BlockNode>> find_children(const std::string& child_name) const {
        std::vector<std::shared_ptr<BlockNode>> result;
        for (const auto& child : children) {
            if (child->name == child_name) {
                result.push_back(child);
            }
        }
        return result;
    }
};

/**
 * @brief 解析错误异常
 * 在解析过程中遇到错误时抛出
 */
class ParseError : public std::exception {
private:
    std::string message_;
    int line_number_;
    
public:
    ParseError(const std::string& msg, int line_num)
        : message_(msg + " at line " + std::to_string(line_num)), line_number_(line_num) {}
    
    const char* what() const noexcept override {
        return message_.c_str();
    }
    
    int line_number() const { return line_number_; }
};

/**
 * @brief Maxwell文件解析器主类
 * 负责解析.aedt文件的层次化结构
 */
class MaxwellParser {
private:
    std::string file_content_;                 // 文件内容
    std::shared_ptr<BlockNode> root_node_;     // 根节点
    std::vector<std::string> lines_;           // 按行分割的内容
    int current_line_;                         // 当前解析行号
    
    /**
     * @brief 预处理文件内容
     */
    void preprocess_content();
    
    /**
     * @brief 跳过空白行和注释
     */
    void skip_whitespace_and_comments();
    
    /**
     * @brief 解析单个块
     */
    std::shared_ptr<BlockNode> parse_block();
    
    /**
     * @brief 解析属性定义
     */
    Property parse_property();
    
    /**
     * @brief 识别数据类型
     */
    DataType identify_data_type(const std::string& value_str);
    
    /**
     * @brief 解析数值
     */
    Value parse_number(const std::string& value_str);
    
    /**
     * @brief 解析字符串
     */
    Value parse_string(const std::string& value_str);
    
    /**
     * @brief 解析布尔值
     */
    Value parse_boolean(const std::string& value_str);
    
    /**
     * @brief 解析数组
     */
    Value parse_array(const std::string& value_str);
    
    /**
     * @brief 解析函数调用
     */
    Value parse_function(const std::string& value_str);
    
    /**
     * @brief 解析集合
     */
    Value parse_set(const std::string& value_str);
    
public:
    /**
     * @brief 构造函数
     */
    MaxwellParser();
    
    /**
     * @brief 从文件路径解析
     */
    bool parse_file(const std::string& file_path);
    
    /**
     * @brief 从字符串内容解析
     */
    bool parse_content(const std::string& content);
    
    /**
     * @brief 获取解析结果根节点
     */
    std::shared_ptr<BlockNode> get_root() const { return root_node_; }
    
    /**
     * @brief 清空解析结果
     */
    void clear();
    
    /**
     * @brief 获取解析错误信息
     */
    std::string get_error_info() const;
    
    /**
     * @brief 验证解析结果
     */
    bool validate() const;
    
    /**
     * @brief 打印解析树结构（用于调试）
     */
    void print_tree(std::ostream& os = std::cout, int indent = 0) const;
};

} // namespace maxwell_parser
} // namespace tool
} // namespace fe_em