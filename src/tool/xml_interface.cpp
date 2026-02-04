/**
 * @file xml_interface.cpp
 * @brief 基础工具层 - XML读写接口模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "xml_interface.hpp"
#include "../lib/pugixml/src/pugixml.hpp"
#include <stdexcept>
#include <algorithm>
#include <sstream>

namespace tool {

// ============================================================================
// XmlNode 实现
// ============================================================================

XmlNode::XmlNode() : node_(nullptr) {
}

XmlNode::XmlNode(pugi::xml_node node) : node_(std::make_shared<pugi::xml_node>(node)) {
}

XmlNode::~XmlNode() {
}

bool XmlNode::isValid() const {
    return node_ && !node_->empty();
}

std::string XmlNode::getName() const {
    if (!isValid()) {
        return "";
    }
    return node_->name();
}

void XmlNode::setName(const std::string& name) {
    if (!isValid()) {
        throw std::runtime_error("Invalid XML node");
    }
    node_->set_name(name.c_str());
}

std::string XmlNode::getText() const {
    if (!isValid()) {
        return "";
    }
    return node_->text().get();
}

void XmlNode::setText(const std::string& text) {
    if (!isValid()) {
        throw std::runtime_error("Invalid XML node");
    }
    node_->text().set(text.c_str());
}

std::optional<std::string> XmlNode::getAttribute(const std::string& name) const {
    if (!isValid()) {
        return std::nullopt;
    }
    
    pugi::xml_attribute attr = node_->attribute(name.c_str());
    if (!attr) {
        return std::nullopt;
    }
    
    return std::string(attr.value());
}

void XmlNode::setAttribute(const std::string& name, const std::string& value) {
    if (!isValid()) {
        throw std::runtime_error("Invalid XML node");
    }
    
    pugi::xml_attribute attr = node_->attribute(name.c_str());
    if (attr) {
        attr.set_value(value.c_str());
    } else {
        node_->append_attribute(name.c_str()).set_value(value.c_str());
    }
}

bool XmlNode::removeAttribute(const std::string& name) {
    if (!isValid()) {
        return false;
    }
    
    pugi::xml_attribute attr = node_->attribute(name.c_str());
    if (!attr) {
        return false;
    }
    
    return node_->remove_attribute(attr);
}

std::vector<std::string> XmlNode::getAttributeNames() const {
    std::vector<std::string> names;
    
    if (!isValid()) {
        return names;
    }
    
    for (pugi::xml_attribute attr : node_->attributes()) {
        names.push_back(attr.name());
    }
    
    return names;
}

XmlNode XmlNode::appendChild(const std::string& name) {
    if (!isValid()) {
        throw std::runtime_error("Invalid XML node");
    }
    
    pugi::xml_node new_node = node_->append_child(name.c_str());
    return XmlNode(new_node);
}

XmlNode XmlNode::getChild(const std::string& name) const {
    if (!isValid()) {
        return XmlNode();
    }
    
    pugi::xml_node child = node_->child(name.c_str());
    return XmlNode(child);
}

std::vector<XmlNode> XmlNode::getChildren() const {
    std::vector<XmlNode> children;
    
    if (!isValid()) {
        return children;
    }
    
    for (pugi::xml_node child : node_->children()) {
        children.push_back(XmlNode(child));
    }
    
    return children;
}

std::vector<XmlNode> XmlNode::getChildren(const std::string& name) const {
    std::vector<XmlNode> children;
    
    if (!isValid()) {
        return children;
    }
    
    for (pugi::xml_node child : node_->children(name.c_str())) {
        children.push_back(XmlNode(child));
    }
    
    return children;
}

bool XmlNode::removeChild(const std::string& name) {
    if (!isValid()) {
        return false;
    }
    
    pugi::xml_node child = node_->child(name.c_str());
    if (!child) {
        return false;
    }
    
    return node_->remove_child(child);
}

XmlNode XmlNode::getParent() const {
    if (!isValid()) {
        return XmlNode();
    }
    
    pugi::xml_node parent = node_->parent();
    return XmlNode(parent);
}

XmlNode XmlNode::getNextSibling() const {
    if (!isValid()) {
        return XmlNode();
    }
    
    pugi::xml_node sibling = node_->next_sibling();
    return XmlNode(sibling);
}

XmlNode XmlNode::getPreviousSibling() const {
    if (!isValid()) {
        return XmlNode();
    }
    
    pugi::xml_node sibling = node_->previous_sibling();
    return XmlNode(sibling);
}

pugi::xml_node XmlNode::getInternalNode() const {
    if (!node_) {
        return pugi::xml_node();
    }
    return *node_;
}

// ============================================================================
// XmlDocument 实现
// ============================================================================

XmlDocument::XmlDocument() : doc_(std::make_unique<pugi::xml_document>()) {
}

XmlDocument::XmlDocument(XmlDocument&& other) noexcept : doc_(std::move(other.doc_)) {
}

XmlDocument& XmlDocument::operator=(XmlDocument&& other) noexcept {
    if (this != &other) {
        doc_ = std::move(other.doc_);
    }
    return *this;
}

XmlDocument::~XmlDocument() {
}

bool XmlDocument::loadFromFile(const std::string& filename) {
    if (!doc_) {
        return false;
    }
    
    pugi::xml_parse_result result = doc_->load_file(filename.c_str());
    return result;
}

bool XmlDocument::loadFromString(const std::string& xml_string) {
    if (!doc_) {
        return false;
    }
    
    pugi::xml_parse_result result = doc_->load_string(xml_string.c_str());
    return result;
}

bool XmlDocument::saveToFile(const std::string& filename, bool indent) const {
    if (!doc_) {
        return false;
    }
    
    return doc_->save_file(filename.c_str(), "  ", indent ? pugi::format_indent : pugi::format_raw);
}

std::string XmlDocument::toString(bool indent) const {
    if (!doc_) {
        return "";
    }
    
    std::ostringstream oss;
    doc_->save(oss, "  ", indent ? pugi::format_indent : pugi::format_raw);
    return oss.str();
}

XmlNode XmlDocument::getRootNode() const {
    if (!doc_) {
        return XmlNode();
    }
    
    pugi::xml_node root = doc_->document_element();
    return XmlNode(root);
}

void XmlDocument::createNew(const std::string& root_name) {
    if (!doc_) {
        return;
    }
    
    doc_->reset();
    pugi::xml_node decl = doc_->append_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";
    
    doc_->append_child(root_name.c_str());
}

bool XmlDocument::isValid() const {
    return doc_ && doc_->document_element();
}

void XmlDocument::reset() {
    if (doc_) {
        doc_->reset();
    }
}

pugi::xml_document* XmlDocument::getInternalDocument() const {
    return doc_.get();
}

// ============================================================================
// XmlFactory 实现
// ============================================================================

XmlDocument XmlFactory::createDocument() {
    XmlDocument doc;
    doc.createNew();
    return std::move(doc);
}

XmlDocument XmlFactory::createDocumentFromFile(const std::string& filename) {
    XmlDocument doc;
    if (doc.loadFromFile(filename)) {
        return std::move(doc);
    }
    return XmlDocument();
}

XmlDocument XmlFactory::createDocumentFromString(const std::string& xml_string) {
    XmlDocument doc;
    if (doc.loadFromString(xml_string)) {
        return std::move(doc);
    }
    return XmlDocument();
}

} // namespace tool