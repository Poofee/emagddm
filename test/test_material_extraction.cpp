/**
 * @file test_material_extraction.cpp
 * @brief 测试材料数据提取功能
 * @details 测试从Maxwell文件中提取材料数据的功能
 * @author AI Developer
 * @date 2026-04-03
 * @version 1.0
 */

#include "tool/maxwell_parser_impl.hpp"
#include "tool/logger_factory.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace tool;

/**
 * @brief 打印材料信息
 */
void print_material(const nlohmann::json& material) {
    std::cout << "\n=== 材料信息 ===" << std::endl;
    std::cout << "名称: " << material["name"].get<std::string>() << std::endl;
    
    // 线性属性
    if (material.contains("conductivity")) {
        std::cout << "电导率: " << material["conductivity"].get<double>() << " S/m" << std::endl;
    }
    if (material.contains("mass_density")) {
        std::cout << "质量密度: " << material["mass_density"].get<double>() << " kg/m³" << std::endl;
    }
    if (material.contains("thermal_conductivity")) {
        std::cout << "热导率: " << material["thermal_conductivity"].get<double>() << " W/(m·K)" << std::endl;
    }
    if (material.contains("specific_heat")) {
        std::cout << "比热容: " << material["specific_heat"].get<double>() << " J/(kg·K)" << std::endl;
    }
    
    // 磁导率
    if (material.contains("permeability_type")) {
        std::string perm_type = material["permeability_type"].get<std::string>();
        std::cout << "磁导率类型: " << perm_type << std::endl;
        
        if (perm_type == "linear" && material.contains("permeability")) {
            std::cout << "相对磁导率: " << material["permeability"].get<double>() << std::endl;
        } else if (perm_type == "nonlinear" && material.contains("bh_curve")) {
            auto bh_curve = material["bh_curve"];
            std::cout << "BH曲线点数: " << bh_curve.size() << std::endl;
            
            // 打印前5个点
            std::cout << "BH曲线（前5个点）:" << std::endl;
            for (size_t i = 0; i < std::min(bh_curve.size(), size_t(5)); ++i) {
                std::cout << "  [" << i << "] H=" << std::setw(10) << bh_curve[i][0].get<double>()
                          << " A/m, B=" << std::setw(10) << bh_curve[i][1].get<double>() << " T" << std::endl;
            }
        }
    }
    
    // 磁芯损耗系数
    if (material.contains("core_loss_kh")) {
        std::cout << "磁芯损耗系数 kh: " << material["core_loss_kh"].get<double>() << std::endl;
    }
    if (material.contains("core_loss_kc")) {
        std::cout << "磁芯损耗系数 kc: " << material["core_loss_kc"].get<double>() << std::endl;
    }
    if (material.contains("core_loss_ke")) {
        std::cout << "磁芯损耗系数 ke: " << material["core_loss_ke"].get<double>() << std::endl;
    }
    
    // 矫顽力
    if (material.contains("coercivity")) {
        auto coercivity = material["coercivity"];
        std::cout << "矫顽力: " << std::endl;
        if (coercivity.contains("magnitude")) {
            std::cout << "  幅值: " << coercivity["magnitude"].get<std::string>() << std::endl;
        }
        if (coercivity.contains("direction")) {
            auto dir = coercivity["direction"];
            std::cout << "  方向: [" << dir[0].get<double>() << ", " 
                      << dir[1].get<double>() << ", " 
                      << dir[2].get<double>() << "]" << std::endl;
        }
    }
}

/**
 * @brief 测试材料提取功能
 */
void test_material_extraction() {
    std::cout << "=== 测试材料数据提取功能 ===" << std::endl;
    
    // 初始化日志系统以便看到调试输出
    tool::LoggerFactory::initializeDefaultLogger("", true, tool::LoggerType::SPDLOG);
    tool::LoggerFactory::setDefaultLoggerLevel(tool::LogLevel::DEBUG);
    
    MaxwellParserImpl parser;
    
    // 测试文件路径 - 使用绝对路径
    std::string file_path = "d:/codes/emagddm/docs/project/Temp.aedt";
    
    std::cout << "测试文件路径: " << file_path << std::endl;
    
    // 检查文件是否可以解析
    if (!parser.canParse(file_path)) {
        std::cout << "文件格式不支持: " << file_path << std::endl;
        return;
    }
    
    std::cout << "文件格式支持，开始解析..." << std::endl;
    
    try {
        // 解析文件信息
        auto file_info = parser.parseFileInfo();
        std::cout << "\n文件信息:" << std::endl;
        std::cout << "  文件路径: " << file_info.file_path << std::endl;
        std::cout << "  文件格式: " << file_info.file_format << std::endl;
        std::cout << "  Maxwell版本: " << file_info.maxwell_version << std::endl;
        std::cout << "  文件大小: " << file_info.file_size << " 字节" << std::endl;
        
        // 解析材料数据
        auto materials = parser.parseMaterials();
        
        std::cout << "\n提取到 " << materials.size() << " 个材料" << std::endl;
        
        // 打印所有材料信息
        for (const auto& material : materials) {
            print_material(material);
        }
        
        // 验证特定材料
        bool found_copper = false;
        bool found_nonlinear = false;
        
        for (const auto& material : materials) {
            std::string name = material["name"].get<std::string>();
            
            // 检查线性材料（copper）
            if (name == "copper") {
                found_copper = true;
                std::cout << "\n验证 copper 材料:" << std::endl;
                assert(material.contains("conductivity"));
                assert(material["conductivity"].get<double>() == 58000000.0);
                std::cout << "  ✓ 电导率正确" << std::endl;
                
                assert(material.contains("permeability"));
                assert(material["permeability"].get<double>() == 0.999991);
                std::cout << "  ✓ 磁导率正确" << std::endl;
                
                assert(material.contains("mass_density"));
                assert(material["mass_density"].get<double>() == 8933.0);
                std::cout << "  ✓ 质量密度正确" << std::endl;
            }
            
            // 检查非线性材料
            if (material.contains("permeability_type") && 
                material["permeability_type"].get<std::string>() == "nonlinear") {
                found_nonlinear = true;
                std::cout << "\n验证非线性材料 " << name << ":" << std::endl;
                assert(material.contains("bh_curve"));
                assert(material["bh_curve"].size() > 0);
                std::cout << "  ✓ BH曲线提取成功，点数: " << material["bh_curve"].size() << std::endl;
            }
        }
        
        assert(found_copper && "未找到 copper 材料");
        assert(found_nonlinear && "未找到非线性材料");
        
        std::cout << "\n✅ 所有材料提取测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ 测试失败: " << e.what() << std::endl;
        throw;
    }
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始材料数据提取测试..." << std::endl;
    
    try {
        test_material_extraction();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}
