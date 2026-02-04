/**
 * @file test_xml.cpp
 * @brief XML接口测试程序
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "tool/xml_interface.hpp"
#include "tool/logger_factory.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>

/**
 * @brief 测试XML文档创建和保存
 */
void testXmlDocumentCreation() {
    FEEM_INFO("开始测试XML文档创建功能", "test_xml");
    
    try {
        // 创建新文档
        tool::XmlDocument doc = tool::XmlFactory::createDocument();
        
        if (!doc.isValid()) {
            FEEM_ERROR("文档创建失败", "test_xml");
            return;
        }
        
        FEEM_INFO("文档创建成功", "test_xml");
        
        // 获取根节点并设置属性
        tool::XmlNode root = doc.getRootNode();
        root.setAttribute("version", "1.0");
        root.setAttribute("created_by", "FE-EM-FETIDP");
        
        // 添加子节点
        tool::XmlNode config_node = root.appendChild("configuration");
        config_node.setAttribute("type", "electromagnetic");
        
        // 添加配置项
        tool::XmlNode solver_node = config_node.appendChild("solver");
        solver_node.setAttribute("name", "FETI-DP");
        solver_node.setAttribute("max_iterations", "1000");
        solver_node.setAttribute("tolerance", "1e-12");
        
        // 添加网格配置
        tool::XmlNode mesh_node = config_node.appendChild("mesh");
        mesh_node.setAttribute("type", "triangular");
        mesh_node.setAttribute("elements", "10000");
        
        // 添加材料属性
        tool::XmlNode materials_node = config_node.appendChild("materials");
        
        tool::XmlNode copper_node = materials_node.appendChild("material");
        copper_node.setAttribute("name", "copper");
        copper_node.setAttribute("conductivity", "5.96e7");
        copper_node.setAttribute("permeability", "1.0");
        
        tool::XmlNode air_node = materials_node.appendChild("material");
        air_node.setAttribute("name", "air");
        air_node.setAttribute("conductivity", "0.0");
        air_node.setAttribute("permeability", "1.0");
        
        // 保存到文件
        std::string filename = "test_config.xml";
        if (doc.saveToFile(filename)) {
            FEEM_INFO("XML文档保存成功: " + filename, "test_xml");
            
            // 读取并验证文件
            tool::XmlDocument loaded_doc = tool::XmlFactory::createDocumentFromFile(filename);
            if (loaded_doc.isValid()) {
                FEEM_INFO("XML文档加载验证成功", "test_xml");
                
                // 验证根节点属性
                tool::XmlNode loaded_root = loaded_doc.getRootNode();
                auto version = loaded_root.getAttribute("version");
                if (version && *version == "1.0") {
                    FEEM_INFO("根节点属性验证成功", "test_xml");
                } else {
                    FEEM_ERROR("根节点属性验证失败", "test_xml");
                }
                
                // 验证配置节点
                tool::XmlNode loaded_config = loaded_root.getChild("configuration");
                if (loaded_config.isValid()) {
                    FEEM_INFO("配置节点验证成功", "test_xml");
                }
                
                // 验证求解器配置
                tool::XmlNode loaded_solver = loaded_config.getChild("solver");
                auto solver_name = loaded_solver.getAttribute("name");
                if (solver_name && *solver_name == "FETI-DP") {
                    FEEM_INFO("求解器配置验证成功", "test_xml");
                }
                
                // 验证材料数量
                tool::XmlNode loaded_materials = loaded_config.getChild("materials");
                auto material_children = loaded_materials.getChildren("material");
                if (material_children.size() == 2) {
                    FEEM_INFO("材料数量验证成功: " + std::to_string(material_children.size()), "test_xml");
                }
                
                // 输出XML内容
                FEEM_DEBUG("生成的XML内容:", "test_xml");
                FEEM_DEBUG(doc.toString(), "test_xml");
                
            } else {
                FEEM_ERROR("XML文档加载验证失败", "test_xml");
            }
            
            // 清理测试文件
            std::filesystem::remove(filename);
            
        } else {
            FEEM_ERROR("XML文档保存失败", "test_xml");
        }
        
    } catch (const std::exception& e) {
        FEEM_ERROR("测试异常: " + std::string(e.what()), "test_xml");
    }
}

/**
 * @brief 测试XML字符串操作
 */
void testXmlStringOperations() {
    FEEM_INFO("开始测试XML字符串操作功能", "test_xml");
    
    try {
        // 从字符串创建文档
        std::string xml_string = R"(
            <simulation>
                <parameters>
                    <time_step>0.001</time_step>
                    <duration>1.0</duration>
                    <output_frequency>100</output_frequency>
                </parameters>
                <boundary_conditions>
                    <bc type="dirichlet" value="0.0"/>
                    <bc type="neumann" value="1.0"/>
                </boundary_conditions>
            </simulation>
        )";
        
        tool::XmlDocument doc = tool::XmlFactory::createDocumentFromString(xml_string);
        
        if (doc.isValid()) {
            FEEM_INFO("XML字符串解析成功", "test_xml");
            
            // 读取参数
            tool::XmlNode root = doc.getRootNode();
            tool::XmlNode params = root.getChild("parameters");
            
            tool::XmlNode time_step = params.getChild("time_step");
            FEEM_INFO("时间步长: " + time_step.getText(), "test_xml");
            
            tool::XmlNode duration = params.getChild("duration");
            FEEM_INFO("模拟时长: " + duration.getText(), "test_xml");
            
            // 修改参数
            time_step.setText("0.0005");
            FEEM_INFO("修改后的时间步长: " + time_step.getText(), "test_xml");
            
            // 读取边界条件
            tool::XmlNode bcs = root.getChild("boundary_conditions");
            auto bc_nodes = bcs.getChildren("bc");
            
            FEEM_INFO("边界条件数量: " + std::to_string(bc_nodes.size()), "test_xml");
            
            for (size_t i = 0; i < bc_nodes.size(); ++i) {
                auto type = bc_nodes[i].getAttribute("type");
                auto value = bc_nodes[i].getAttribute("value");
                
                if (type && value) {
                    FEEM_INFO("边界条件 " + std::to_string(i) + ": type=" + *type + ", value=" + *value, "test_xml");
                }
            }
            
            // 添加新的边界条件
            tool::XmlNode new_bc = bcs.appendChild("bc");
            new_bc.setAttribute("type", "mixed");
            new_bc.setAttribute("value", "0.5");
            
            FEEM_INFO("添加新边界条件后的数量: " + std::to_string(bcs.getChildren("bc").size()), "test_xml");
            
            // 输出修改后的XML
            FEEM_DEBUG("修改后的XML内容:", "test_xml");
            FEEM_DEBUG(doc.toString(), "test_xml");
            
        } else {
            FEEM_ERROR("XML字符串解析失败", "test_xml");
        }
        
    } catch (const std::exception& e) {
        FEEM_ERROR("测试异常: " + std::string(e.what()), "test_xml");
    }
}

/**
 * @brief 测试XML节点操作
 */
void testXmlNodeOperations() {
    FEEM_INFO("开始测试XML节点操作功能", "test_xml");
    
    try {
        // 创建测试文档
        tool::XmlDocument doc = tool::XmlFactory::createDocument();
        tool::XmlNode root = doc.getRootNode();
        root.setName("test_data");
        
        // 添加多个相同名称的节点
        for (int i = 0; i < 5; ++i) {
            tool::XmlNode data_node = root.appendChild("data_point");
            data_node.setAttribute("id", std::to_string(i));
            data_node.setAttribute("value", std::to_string(i * 10.0));
            data_node.setText("Point " + std::to_string(i));
        }
        
        // 测试获取所有子节点
        auto all_children = root.getChildren();
        FEEM_INFO("总子节点数量: " + std::to_string(all_children.size()), "test_xml");
        
        // 测试获取特定名称的子节点
        auto data_points = root.getChildren("data_point");
        FEEM_INFO("数据点数量: " + std::to_string(data_points.size()), "test_xml");
        
        // 验证每个数据点
        for (const auto& point : data_points) {
            auto id = point.getAttribute("id");
            auto value = point.getAttribute("value");
            
            if (id && value) {
                FEEM_INFO("数据点 ID=" + *id + ", 值=" + *value + ", 文本=" + point.getText(), "test_xml");
            }
        }
        
        // 测试删除节点
        if (root.removeChild("data_point")) {
            FEEM_INFO("删除第一个数据点成功", "test_xml");
            FEEM_INFO("删除后数据点数量: " + std::to_string(root.getChildren("data_point").size()), "test_xml");
        }
        
        // 测试属性操作
        tool::XmlNode first_point = root.getChild("data_point");
        if (first_point.isValid()) {
            // 获取所有属性名称
            auto attr_names = first_point.getAttributeNames();
            FEEM_INFO("第一个数据点的属性数量: " + std::to_string(attr_names.size()), "test_xml");
            
            for (const auto& name : attr_names) {
                auto value = first_point.getAttribute(name);
                if (value) {
                    FEEM_INFO("属性 " + name + " = " + *value, "test_xml");
                }
            }
            
            // 修改属性
            first_point.setAttribute("modified", "true");
            
            // 删除属性
            if (first_point.removeAttribute("id")) {
                FEEM_INFO("删除ID属性成功", "test_xml");
            }
            
            FEEM_INFO("修改后属性数量: " + std::to_string(first_point.getAttributeNames().size()), "test_xml");
        }
        
        FEEM_DEBUG("节点操作测试后的XML内容:", "test_xml");
        FEEM_DEBUG(doc.toString(), "test_xml");
        
    } catch (const std::exception& e) {
        FEEM_ERROR("测试异常: " + std::string(e.what()), "test_xml");
    }
}

/**
 * @brief 主测试函数
 */
int main() {
    FEEM_INFO("开始XML接口测试", "test_xml");
    
    try {
        // 测试XML文档创建和保存
        testXmlDocumentCreation();
        
        // 测试XML字符串操作
        testXmlStringOperations();
        
        // 测试XML节点操作
        testXmlNodeOperations();
        
        FEEM_INFO("XML接口测试完成", "test_xml");
        
    } catch (const std::exception& e) {
        FEEM_ERROR("测试程序异常: " + std::string(e.what()), "test_xml");
        return 1;
    }
    
    return 0;
}