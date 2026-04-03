/**
 * @file test_material_simple.cpp
 * @brief 简单测试材料数据提取功能
 */

#include "tool/maxwell_parser_impl.hpp"
#include "tool/logger_factory.hpp"
#include <iostream>
#include <cassert>

using namespace tool;

int main() {
    std::cout << "=== 简单材料提取测试 ===" << std::endl;
    
    // 初始化日志
    LoggerFactory::initializeDefaultLogger("", true, LoggerType::SPDLOG);
    LoggerFactory::setDefaultLoggerLevel(LogLevel::INFO);
    
    MaxwellParserImpl parser;
    
    // 测试简单的材料数据
    std::string test_content = R"($begin 'AnsoftProject'
	$begin 'Definitions'
		$begin 'Materials'
			$begin 'copper'
				CoordinateSystemType='Cartesian'
				BulkOrSurfaceType=1
				permeability='0.999991'
				conductivity='58000000'
				mass_density='8933'
				thermal_conductivity='400'
				specific_heat='385'
			$end 'copper'
			$begin 'B27AV1400'
				CoordinateSystemType='Cartesian'
				BulkOrSurfaceType=1
				$begin 'permeability'
					property_type='nonlinear'
					BTypeForSingleCurve='normal'
					HUnit='A_per_meter'
					BUnit='tesla'
					$begin 'BHCoordinates'
						DimUnits[2: '', '']
						Points[6: 0, 0, 100, 0.5, 200, 1.0]
					$end 'BHCoordinates'
				$end 'permeability'
				conductivity='2000000'
				core_loss_kh='156.98695'
				core_loss_kc='0.27355194'
				core_loss_ke='0'
				mass_density='7420.5'
				$begin 'magnetic_coercivity'
					property_type='VectorProperty'
					Magnitude='0A_per_meter'
					DirComp1='0'
					DirComp2='0'
					DirComp3='0'
				$end 'magnetic_coercivity'
			$end 'B27AV1400'
		$end 'Materials'
	$end 'Definitions'
$end 'AnsoftProject')";
    
    // 使用内部解析器解析内容
    auto& maxwell_parser = parser.getParser();
    bool result = maxwell_parser.parse_content(test_content);
    
    if (!result) {
        std::cerr << "解析失败" << std::endl;
        return 1;
    }
    
    std::cout << "解析成功" << std::endl;
    
    // 提取材料
    auto materials = parser.parseMaterials();
    
    std::cout << "提取到 " << materials.size() << " 个材料" << std::endl;
    
    // 验证材料
    for (const auto& material : materials) {
        std::cout << "\n材料名称: " << material["name"].get<std::string>() << std::endl;
        
        if (material.contains("permeability")) {
            std::cout << "  磁导率: " << material["permeability"].get<double>() << std::endl;
        }
        
        if (material.contains("conductivity")) {
            std::cout << "  电导率: " << material["conductivity"].get<double>() << std::endl;
        }
        
        if (material.contains("bh_curve")) {
            std::cout << "  BH曲线点数: " << material["bh_curve"].size() << std::endl;
        }
    }
    
    // 验证 copper 材料
    bool found_copper = false;
    for (const auto& material : materials) {
        if (material["name"].get<std::string>() == "copper") {
            found_copper = true;
            assert(material["conductivity"].get<double>() == 58000000.0);
            assert(material["permeability"].get<double>() == 0.999991);
            assert(material["mass_density"].get<double>() == 8933.0);
            std::cout << "\n✓ copper 材料验证通过" << std::endl;
        }
    }
    
    assert(found_copper && "未找到 copper 材料");
    
    std::cout << "\n✅ 所有测试通过!" << std::endl;
    return 0;
}
