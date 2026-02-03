/**
 * @file xml_interface.hpp
 * @brief 基础工具层 - XML读写接口模块头文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#ifndef FE_EM_FETIDP_XML_INTERFACE_HPP
#define FE_EM_FETIDP_XML_INTERFACE_HPP

#include <string>
#include <vector>
#include <memory>
#include <optional>

// pugixml 前向声明
namespace pugi {
    class xml_document;
    class xml_node;
    class xml_attribute;
}

namespace tool {

/**
 * @brief XML节点类，封装pugixml的节点操作
 */
class XmlNode {
public:
    /**
     * @brief 默认构造函数
     */
    XmlNode();
    
    /**
     * @brief 构造函数
     * @param node pugixml节点指针
     */
    explicit XmlNode(pugi::xml_node node);
    
    /**
     * @brief 析构函数
     */
    ~XmlNode();
    
    /**
     * @brief 检查节点是否有效
     * @return 节点是否有效
     */
    bool isValid() const;
    
    /**
     * @brief 获取节点名称
     * @return 节点名称
     */
    std::string getName() const;
    
    /**
     * @brief 设置节点名称
     * @param name 节点名称
     */
    void setName(const std::string& name);
    
    /**
     * @brief 获取节点文本内容
     * @return 节点文本内容
     */
    std::string getText() const;
    
    /**
     * @brief 设置节点文本内容
     * @param text 文本内容
     */
    void setText(const std::string& text);
    
    /**
     * @brief 获取属性值
     * @param name 属性名称
     * @return 属性值，如果不存在返回std::nullopt
     */
    std::optional<std::string> getAttribute(const std::string& name) const;
    
    /**
     * @brief 设置属性值
     * @param name 属性名称
     * @param value 属性值
     */
    void setAttribute(const std::string& name, const std::string& value);
    
    /**
     * @brief 删除属性
     * @param name 属性名称
     * @return 是否删除成功
     */
    bool removeAttribute(const std::string& name);
    
    /**
     * @brief 获取所有属性名称
     * @return 属性名称列表
     */
    std::vector<std::string> getAttributeNames() const;
    
    /**
     * @brief 添加子节点
     * @param name 子节点名称
     * @return 新创建的节点
     */
    XmlNode appendChild(const std::string& name);
    
    /**
     * @brief 获取子节点
     * @param name 子节点名称
     * @return 子节点，如果不存在返回无效节点
     */
    XmlNode getChild(const std::string& name) const;
    
    /**
     * @brief 获取所有子节点
     * @return 子节点列表
     */
    std::vector<XmlNode> getChildren() const;
    
    /**
     * @brief 获取所有指定名称的子节点
     * @param name 子节点名称
     * @return 子节点列表
     */
    std::vector<XmlNode> getChildren(const std::string& name) const;
    
    /**
     * @brief 删除子节点
     * @param name 子节点名称
     * @return 是否删除成功
     */
    bool removeChild(const std::string& name);
    
    /**
     * @brief 获取父节点
     * @return 父节点
     */
    XmlNode getParent() const;
    
    /**
     * @brief 获取下一个兄弟节点
     * @return 下一个兄弟节点
     */
    XmlNode getNextSibling() const;
    
    /**
     * @brief 获取前一个兄弟节点
     * @return 前一个兄弟节点
     */
    XmlNode getPreviousSibling() const;
    
    /**
     * @brief 获取内部pugixml节点
     * @return pugixml节点
     */
    pugi::xml_node getInternalNode() const;

private:
    std::shared_ptr<pugi::xml_node> node_;
};

/**
 * @brief XML文档类，封装pugixml的文档操作
 */
class XmlDocument {
public:
    /**
     * @brief 默认构造函数
     */
    XmlDocument();
    
    /**
     * @brief 移动构造函数
     */
    XmlDocument(XmlDocument&& other) noexcept;
    
    /**
     * @brief 移动赋值运算符
     */
    XmlDocument& operator=(XmlDocument&& other) noexcept;
    
    /**
     * @brief 析构函数
     */
    ~XmlDocument();
    
    /**
     * @brief 从文件加载XML文档
     * @param filename 文件名
     * @return 是否加载成功
     */
    bool loadFromFile(const std::string& filename);
    
    /**
     * @brief 从字符串加载XML文档
     * @param xml_string XML字符串
     * @return 是否加载成功
     */
    bool loadFromString(const std::string& xml_string);
    
    /**
     * @brief 保存XML文档到文件
     * @param filename 文件名
     * @param indent 是否格式化输出
     * @return 是否保存成功
     */
    bool saveToFile(const std::string& filename, bool indent = true) const;
    
    /**
     * @brief 获取XML文档字符串
     * @param indent 是否格式化输出
     * @return XML字符串
     */
    std::string toString(bool indent = true) const;
    
    /**
     * @brief 获取根节点
     * @return 根节点
     */
    XmlNode getRootNode() const;
    
    /**
     * @brief 创建新的XML文档
     * @param root_name 根节点名称
     */
    void createNew(const std::string& root_name = "root");
    
    /**
     * @brief 检查文档是否有效
     * @return 文档是否有效
     */
    bool isValid() const;
    
    /**
     * @brief 重置文档
     */
    void reset();
    
    /**
     * @brief 获取内部pugixml文档
     * @return pugixml文档
     */
    pugi::xml_document* getInternalDocument() const;

private:
    std::unique_ptr<pugi::xml_document> doc_;
};

/**
 * @brief XML接口工厂类
 */
class XmlFactory {
public:
    /**
     * @brief 创建XML文档
     * @return XML文档对象
     */
    static XmlDocument createDocument();
    
    /**
     * @brief 从文件创建XML文档
     * @param filename 文件名
     * @return XML文档对象，如果失败返回无效文档
     */
    static XmlDocument createDocumentFromFile(const std::string& filename);
    
    /**
     * @brief 从字符串创建XML文档
     * @param xml_string XML字符串
     * @return XML文档对象，如果失败返回无效文档
     */
    static XmlDocument createDocumentFromString(const std::string& xml_string);
};

} // namespace tool

#endif // FE_EM_FETIDP_XML_INTERFACE_HPP