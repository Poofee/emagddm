/**
 * @file test_em_tools.cpp
 * @brief 电磁模块工具类单元测试程序（简化版）
 */

#include "tool/em_enums.hpp"
#include "tool/em_exception.hpp"
#include "tool/id_generator.hpp"
#include "tool/logger_factory.hpp"
#include <iostream>
#include <cassert>
#include <sstream>

int test_count = 0;
int pass_count = 0;
int fail_count = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        test_count++; \
        if (condition) { \
            pass_count++; \
            std::cout << "  [PASS] " << message << std::endl; \
        } else { \
            fail_count++; \
            std::cout << "  [FAIL] " << message << std::endl; \
        } \
    } while(0)

#define TEST_SECTION(name) \
    std::cout << "\n=== " << name << " ===" << std::endl

void testDimTypeEnum() {
    TEST_SECTION("维度类型枚举测试");
    
    TEST_ASSERT(tool::dimTypeToString(tool::DimType::D2) == "D2", 
                "D2转换为字符串正确");
    TEST_ASSERT(tool::dimTypeToString(tool::DimType::D3) == "D3", 
                "D3转换为字符串正确");
    TEST_ASSERT(tool::dimTypeToString(tool::DimType::AXIS) == "AXIS", 
                "AXIS转换为字符串正确");
    
    TEST_ASSERT(tool::stringToDimType("D2") == tool::DimType::D2, 
                "字符串D2解析正确");
    TEST_ASSERT(tool::stringToDimType("2D") == tool::DimType::D2, 
                "字符串2D解析正确");
    TEST_ASSERT(tool::stringToDimType("AXIS") == tool::DimType::AXIS, 
                "字符串AXIS解析正确");
}

void testFieldTypeEnum() {
    TEST_SECTION("场量类型枚举测试");
    
    TEST_ASSERT(tool::fieldTypeToString(tool::FieldType::SCALAR) == "SCALAR", 
                "SCALAR转换为字符串正确");
    TEST_ASSERT(tool::fieldTypeToString(tool::FieldType::VECTOR) == "VECTOR", 
                "VECTOR转换为字符串正确");
    
    TEST_ASSERT(tool::stringToFieldType("SCALAR") == tool::FieldType::SCALAR, 
                "字符串SCALAR解析正确");
    TEST_ASSERT(tool::stringToFieldType("Scalar") == tool::FieldType::SCALAR, 
                "字符串Scalar解析正确");
}

void testMatTypeEnum() {
    TEST_SECTION("材料类型枚举测试");
    
    TEST_ASSERT(tool::matTypeToString(tool::MatType::LINEAR_ISOTROPIC) == "LINEAR_ISOTROPIC", 
                "LINEAR_ISOTROPIC转换为字符串正确");
    TEST_ASSERT(tool::matTypeToString(tool::MatType::PERMANENT_MAGNET) == "PERMANENT_MAGNET", 
                "PERMANENT_MAGNET转换为字符串正确");
    
    TEST_ASSERT(tool::stringToMatType("LINEAR_ISOTROPIC") == tool::MatType::LINEAR_ISOTROPIC, 
                "字符串LINEAR_ISOTROPIC解析正确");
    TEST_ASSERT(tool::stringToMatType("Linear Isotropic") == tool::MatType::LINEAR_ISOTROPIC, 
                "字符串Linear Isotropic解析正确");
}

void testBndTypeEnum() {
    TEST_SECTION("边界类型枚举测试");
    
    TEST_ASSERT(tool::bndTypeToString(tool::BndType::DIRICHLET) == "DIRICHLET", 
                "DIRICHLET转换为字符串正确");
    TEST_ASSERT(tool::bndTypeToString(tool::BndType::PERFECT_E) == "PERFECT_E", 
                "PERFECT_E转换为字符串正确");
    TEST_ASSERT(tool::bndTypeToString(tool::BndType::BALLOON) == "BALLOON", 
                "BALLOON转换为字符串正确");
    
    TEST_ASSERT(tool::stringToBndType("DIRICHLET") == tool::BndType::DIRICHLET, 
                "字符串DIRICHLET解析正确");
    TEST_ASSERT(tool::stringToBndType("Fixed") == tool::BndType::DIRICHLET, 
                "字符串Fixed解析为DIRICHLET正确");
    TEST_ASSERT(tool::stringToBndType("PERFECT_E") == tool::BndType::PERFECT_E, 
                "字符串PERFECT_E解析正确");
}

void testExcitationTypeEnum() {
    TEST_SECTION("激励类型枚举测试");
    
    TEST_ASSERT(tool::excitationTypeToString(tool::ExcitationType::COIL) == "COIL", 
                "COIL转换为字符串正确");
    TEST_ASSERT(tool::excitationTypeToString(tool::ExcitationType::WINDING) == "WINDING", 
                "WINDING转换为字符串正确");
    
    TEST_ASSERT(tool::stringToExcitationType("COIL") == tool::ExcitationType::COIL, 
                "字符串COIL解析正确");
    TEST_ASSERT(tool::stringToExcitationType("Coil") == tool::ExcitationType::COIL, 
                "字符串Coil解析为COIL正确");
    TEST_ASSERT(tool::stringToExcitationType("Current Density") == tool::ExcitationType::CURRENT_DENSITY, 
                "字符串Current Density解析正确");
}

void testProjectFileTypeEnum() {
    TEST_SECTION("项目文件类型枚举测试");
    
    TEST_ASSERT(tool::projectFileTypeToString(tool::ProjectFileType::AEDT) == "AEDT", 
                "AEDT转换为字符串正确");
    TEST_ASSERT(tool::projectFileTypeToString(tool::ProjectFileType::EMF) == "EMF", 
                "EMF转换为字符串正确");
    
    TEST_ASSERT(tool::stringToProjectFileType("AEDT") == tool::ProjectFileType::AEDT, 
                "字符串AEDT解析正确");
    TEST_ASSERT(tool::stringToProjectFileType(".emf") == tool::ProjectFileType::EMF, 
                "字符串.emf解析正确");
}

void testMaxwellVersionEnum() {
    TEST_SECTION("Maxwell版本枚举测试");
    
    TEST_ASSERT(tool::maxwellVersionToString(tool::MaxwellVersion::R22) == "R22", 
                "R22转换为字符串正确");
    TEST_ASSERT(tool::maxwellVersionToString(tool::MaxwellVersion::R24) == "R24", 
                "R24转换为字符串正确");
    
    TEST_ASSERT(tool::stringToMaxwellVersion("R22") == tool::MaxwellVersion::R22, 
                "字符串R22解析正确");
    TEST_ASSERT(tool::stringToMaxwellVersion("2022") == tool::MaxwellVersion::R22, 
                "字符串2022解析为R22正确");
    TEST_ASSERT(tool::stringToMaxwellVersion("NEWER") == tool::MaxwellVersion::NEWER, 
                "字符串NEWER解析正确");
}

void testIDGenerator() {
    TEST_SECTION("ID生成器测试");
    
    tool::IDGenerator& id_gen = tool::IDGenerator::getInstance();
    
    uint64_t id1 = id_gen.generateID(tool::IDCategory::MATERIAL);
    TEST_ASSERT(id1 >= 1000, "材料ID在正确范围内");
    
    uint64_t id2 = id_gen.generateID(tool::IDCategory::GEOMETRY);
    TEST_ASSERT(id2 >= 10000, "几何ID在正确范围内");
    
    std::string id_str = id_gen.generateIDString(tool::IDCategory::MATERIAL, 1000);
    TEST_ASSERT(id_str == "Mat_1000", "生成ID字符串正确");
}

void testEntityIDGenerator() {
    TEST_SECTION("实体ID生成器测试");
    
    tool::EntityIDGenerator& entity_gen = tool::EntityIDGenerator::getInstance();
    
    uint64_t mat_id1 = entity_gen.generateMaterialID("Copper");
    TEST_ASSERT(mat_id1 >= 1000, "Copper材料ID在正确范围内");
    
    uint64_t mat_id2 = entity_gen.generateMaterialID("Aluminum");
    TEST_ASSERT(mat_id2 > mat_id1, "不同材料ID递增");
    
    uint64_t same_mat_id = entity_gen.generateMaterialID("Copper");
    TEST_ASSERT(same_mat_id == mat_id1, "相同材料返回相同ID");
}

void testExceptionClasses() {
    TEST_SECTION("异常类测试");
    
    try {
        throw tool::project::ProjectNotFoundException("test.aedt");
    } catch (const tool::project::ProjectNotFoundException& e) {
        TEST_ASSERT(std::string(e.what()).find("not found") != std::string::npos, 
                    "ProjectNotFoundException消息包含'not found'");
        TEST_ASSERT(e.getFilePath() == "test.aedt", 
                    "ProjectNotFoundException获取文件路径正确");
    }
    
    try {
        throw tool::format::XMLParseException(10, 5, "Invalid element");
    } catch (const tool::format::XMLParseException& e) {
        TEST_ASSERT(e.getLineNumber() == 10, 
                    "XMLParseException获取行号正确");
        TEST_ASSERT(e.getColumnNumber() == 5, 
                    "XMLParseException获取列号正确");
    }
    
    try {
        throw tool::boundary::BoundaryConflictException("Entity1", "Dirichlet", "Neumann");
    } catch (const tool::boundary::BoundaryConflictException& e) {
        TEST_ASSERT(std::string(e.what()).find("Boundary conflict") != std::string::npos, 
                    "BoundaryConflictException消息包含'Boundary conflict'");
    }
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   电磁模块工具类单元测试" << std::endl;
    std::cout << "========================================" << std::endl;
    
    try {
        testDimTypeEnum();
        testFieldTypeEnum();
        testMatTypeEnum();
        testBndTypeEnum();
        testExcitationTypeEnum();
        testProjectFileTypeEnum();
        testMaxwellVersionEnum();
        testIDGenerator();
        testEntityIDGenerator();
        testExceptionClasses();
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "   测试结果汇总" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "  总测试数: " << test_count << std::endl;
        std::cout << "  通过: " << pass_count << std::endl;
        std::cout << "  失败: " << fail_count << std::endl;
        std::cout << "========================================" << std::endl;
        
        if (fail_count > 0) {
            std::cout << "测试存在失败项，请检查上述输出" << std::endl;
            return 1;
        } else {
            std::cout << "所有测试通过！" << std::endl;
            return 0;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "测试过程中发生异常: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "测试过程中发生未知异常" << std::endl;
        return 1;
    }
}
