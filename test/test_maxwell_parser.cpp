#include "tool/maxwell_parser.hpp"
#include <iostream>
#include <cassert>

using namespace fe_em::tool::maxwell_parser;

/**
 * @brief 测试基础解析功能
 */
void test_basic_parsing() {
    std::cout << "=== 测试基础解析功能 ===" << std::endl;
    
    MaxwellParser parser;
    
    // 测试简单块结构
    std::string test_content = R"(
$begin 'TestBlock'
    Name = 'TestName'
    Value = 123.45
    Flag = true
$end 'TestBlock'
)";
    
    bool result = parser.parse_content(test_content);
    assert(result && "基础解析失败");
    
    auto root = parser.get_root();
    assert(root && "根节点为空");
    assert(root->name == "TestBlock" && "块名称不匹配");
    
    // 验证属性
    auto name_prop = root->find_property("Name");
    assert(name_prop.has_value() && "Name属性未找到");
    assert(std::get<std::string>(name_prop->value) == "TestName" && "Name值不匹配");
    
    auto value_prop = root->find_property("Value");
    assert(value_prop.has_value() && "Value属性未找到");
    assert(std::get<double>(value_prop->value) == 123.45 && "Value值不匹配");
    
    auto flag_prop = root->find_property("Flag");
    assert(flag_prop.has_value() && "Flag属性未找到");
    assert(std::get<bool>(flag_prop->value) == true && "Flag值不匹配");
    
    std::cout << "基础解析测试通过" << std::endl;
}

/**
 * @brief 测试嵌套块结构
 */
void test_nested_blocks() {
    std::cout << "\n=== 测试嵌套块结构 ===" << std::endl;
    
    MaxwellParser parser;
    
    std::string test_content = R"(
$begin 'ParentBlock'
    ParentProp = 'ParentValue'
    $begin 'ChildBlock'
        ChildProp = 42
    $end 'ChildBlock'
$end 'ParentBlock'
)";
    
    bool result = parser.parse_content(test_content);
    assert(result && "嵌套块解析失败");
    
    auto root = parser.get_root();
    assert(root && "根节点为空");
    
    // 验证父块属性
    auto parent_prop = root->find_property("ParentProp");
    assert(parent_prop.has_value() && "ParentProp属性未找到");
    
    // 验证子块
    auto children = root->find_children("ChildBlock");
    assert(children.size() == 1 && "子块数量不匹配");
    
    auto child_prop = children[0]->find_property("ChildProp");
    assert(child_prop.has_value() && "ChildProp属性未找到");
    assert(std::get<double>(child_prop->value) == 42 && "ChildProp值不匹配");
    
    std::cout << "嵌套块测试通过" << std::endl;
}

/**
 * @brief 测试数组解析
 */
void test_array_parsing() {
    std::cout << "\n=== 测试数组解析 ===" << std::endl;
    
    MaxwellParser parser;
    
    std::string test_content = R"(
$begin 'ArrayTest'
    IntArray = [5: 1, 2, 3, 4, 5]
    DoubleArray = [3: 1.5, 2.5, 3.5]
$end 'ArrayTest'
)";
    
    bool result = parser.parse_content(test_content);
    assert(result && "数组解析失败");
    
    auto root = parser.get_root();
    
    // 测试整数数组
    auto int_array_prop = root->find_property("IntArray");
    assert(int_array_prop.has_value() && "IntArray属性未找到");
    
    auto int_array = std::get<std::vector<Value>>(int_array_prop->value);
    assert(int_array.size() == 5 && "整数数组大小不匹配");
    
    for (int i = 0; i < 5; ++i) {
        assert(std::get<double>(int_array[i]) == i + 1 && "整数数组值不匹配");
    }
    
    // 测试浮点数数组
    auto double_array_prop = root->find_property("DoubleArray");
    assert(double_array_prop.has_value() && "DoubleArray属性未找到");
    
    auto double_array = std::get<std::vector<Value>>(double_array_prop->value);
    assert(double_array.size() == 3 && "浮点数数组大小不匹配");
    
    for (int i = 0; i < 3; ++i) {
        assert(std::get<double>(double_array[i]) == i + 1.5 && "浮点数数组值不匹配");
    }
    
    std::cout << "数组解析测试通过" << std::endl;
}

/**
 * @brief 测试函数和集合解析
 */
void test_function_and_set() {
    std::cout << "\n=== 测试函数和集合解析 ===" << std::endl;
    
    MaxwellParser parser;
    
    std::string test_content = R"(
$begin 'FunctionTest'
    VersionFunc = Version(1, 0)
    StringSet = set('item1', 'item2', 'item3')
$end 'FunctionTest'
)";
    
    bool result = parser.parse_content(test_content);
    assert(result && "函数和集合解析失败");
    
    auto root = parser.get_root();
    
    // 测试函数调用
    auto func_prop = root->find_property("VersionFunc");
    assert(func_prop.has_value() && "VersionFunc属性未找到");
    assert(std::get<std::string>(func_prop->value) == "Version(1, 0)" && "函数值不匹配");
    
    // 测试集合
    auto set_prop = root->find_property("StringSet");
    assert(set_prop.has_value() && "StringSet属性未找到");
    
    auto set_items = std::get<std::vector<std::string>>(set_prop->value);
    assert(set_items.size() == 3 && "集合大小不匹配");
    assert(set_items[0] == "item1" && "集合项1不匹配");
    assert(set_items[1] == "item2" && "集合项2不匹配");
    assert(set_items[2] == "item3" && "集合项3不匹配");
    
    std::cout << "函数和集合解析测试通过" << std::endl;
}

/**
 * @brief 测试实际Maxwell文件片段
 */
void test_real_maxwell_snippet() {
    std::cout << "\n=== 测试实际Maxwell文件片段 ===" << std::endl;
    
    MaxwellParser parser;
    
    // 使用实际的Maxwell文件片段进行测试
    std::string test_content = R"($begin 'AnsoftProject'
	Created='Mon Oct 13 11:37:38 2025'
	Product='ElectronicsDesktop'
	FileOwnedByWorkbench=false
	$begin 'Desktop'
		Version(2024, 1)
		InfrastructureVersion(1, 0)
	$end 'Desktop'
	UsesAdvancedFeatures=false
	NextUniqueID=0
	MoveBackwards=false
$end 'AnsoftProject')";
    
    bool result = parser.parse_content(test_content);
    assert(result && "实际文件片段解析失败");
    
    auto root = parser.get_root();
    assert(root && "根节点为空");
    assert(root->name == "AnsoftProject" && "根块名称不匹配");
    
    // 验证根块属性
    auto created_prop = root->find_property("Created");
    assert(created_prop.has_value() && "Created属性未找到");
    
    auto product_prop = root->find_property("Product");
    assert(product_prop.has_value() && "Product属性未找到");
    assert(std::get<std::string>(product_prop->value) == "ElectronicsDesktop" && "Product值不匹配");
    
    // 验证子块
    auto desktop_blocks = root->find_children("Desktop");
    assert(desktop_blocks.size() == 1 && "Desktop块数量不匹配");
    
    std::cout << "实际Maxwell文件片段测试通过" << std::endl;
}

/**
 * @brief 测试错误处理
 */
void test_error_handling() {
    std::cout << "\n=== 测试错误处理 ===" << std::endl;
    
    MaxwellParser parser;
    
    // 测试不匹配的块结束标记
    std::string test_content = R"($begin 'TestBlock'
    Prop = 'Value'
$end 'WrongBlock')";
    
    bool result = parser.parse_content(test_content);
    assert(!result && "错误处理失败：应该检测到块结束标记不匹配");
    
    std::cout << "错误处理测试通过" << std::endl;
}

/**
 * @brief 测试文件解析
 */
void test_file_parsing() {
    std::cout << "\n=== 测试文件解析 ===" << std::endl;
    
    MaxwellParser parser;
    
    // 测试从文件解析
    bool result = parser.parse_file("docs/project/Temp.aedt");
    
    if (result) {
        auto root = parser.get_root();
        assert(root && "文件解析根节点为空");
        
        std::cout << "文件解析成功，根块名称: " << root->name << std::endl;
        
        // 打印解析树结构（用于调试）
        parser.print_tree();
    } else {
        std::cout << "文件解析失败，错误信息: " << parser.get_error_info() << std::endl;
    }
    
    std::cout << "文件解析测试完成" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始Maxwell解析器测试..." << std::endl;
    
    try {
        test_basic_parsing();
        test_nested_blocks();
        test_array_parsing();
        test_function_and_set();
        test_real_maxwell_snippet();
        test_error_handling();
        test_file_parsing();
        
        std::cout << "\n✅ 所有测试通过!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}