/**
 * @file test_maxwell_converter_main.cpp
 * @brief Maxwell数据转换器主测试程序
 * @details 专门测试Maxwell解析数据到内部数据模型的转换功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include "tool/maxwell_converter_impl.hpp"
#include "tool/maxwell_parser.hpp"
#include <memory>

using namespace tool;

class MaxwellConverterMainTest : public ::testing::Test {
protected:
    void SetUp() override {
        converter_ = std::make_unique<MaxwellConverterImpl>();
        parser_ = std::make_unique<fe_em::tool::maxwell_parser::MaxwellParser>();
    }

    void TearDown() override {
        converter_.reset();
        parser_.reset();
    }

    std::unique_ptr<MaxwellConverterImpl> converter_;
    std::unique_ptr<fe_em::tool::maxwell_parser::MaxwellParser> parser_;
};

/**
 * @brief 测试完整的边界条件转换流程
 */
TEST_F(MaxwellConverterMainTest, CompleteBoundaryConversionTest) {
    std::cout << "=== 测试边界条件转换功能 ===" << std::endl;
    
    // 测试边界条件数据
    std::string boundary_content = R"(
$begin 'Boundary1'
    Name = 'ZeroPotential'
    Type = 'Dirichlet'
    Value = 0.0
    Faces = 'Face1' 'Face2'
    Edges = 'Edge1'
    Objects = 'Object1'
    MasterName = 'MasterBoundary'
    SlaveName = 'SlaveBoundary'
    SubType = 'Impedance'
    RadiationDistance = 10.0
    PerfectESymmetry = true
    PerfectHSymmetry = false
    InfiniteSphereRadius = 100.0
$end 'Boundary1'
    )";

    // 解析边界条件数据
    bool parse_result = parser_->parse_content(boundary_content);
    ASSERT_TRUE(parse_result) << "边界条件数据解析失败";

    auto root = parser_->get_root();
    ASSERT_TRUE(root) << "根节点为空";

    std::cout << "解析成功，块名称: " << root->name << std::endl;

    // 转换边界条件数据
    auto boundary = converter_->convertBoundaryDirect(root);
    ASSERT_TRUE(boundary) << "边界条件转换失败";

    // 验证转换结果
    EXPECT_EQ(boundary->getName(), "ZeroPotential");
    EXPECT_EQ(boundary->getType(), BndType::DIRICHLET);
    EXPECT_EQ(boundary->getMaxwellBoundaryID(), "Boundary1");
    
    std::cout << "边界条件转换成功:" << std::endl;
    std::cout << "  名称: " << boundary->getName() << std::endl;
    std::cout << "  类型: Dirichlet" << std::endl;
    std::cout << "  Maxwell ID: " << boundary->getMaxwellBoundaryID() << std::endl;
    
    // 验证关联几何对象
    const auto& faces = boundary->getFaces();
    std::cout << "  关联面数: " << faces.size() << std::endl;
    
    const auto& edges = boundary->getEdges();
    std::cout << "  关联边数: " << edges.size() << std::endl;
    
    const auto& objects = boundary->getObjects();
    std::cout << "  关联对象数: " << objects.size() << std::endl;
    
    std::cout << "边界条件测试通过" << std::endl;
}

/**
 * @brief 测试完整的激励源转换流程
 */
TEST_F(MaxwellConverterMainTest, CompleteExcitationConversionTest) {
    std::cout << "=== 测试激励源转换功能 ===" << std::endl;
    
    // 测试激励源数据
    std::string excitation_content = R"(
$begin 'Excitation1'
    Name = 'CurrentSource'
    Type = 'Current'
    Value = 10.0
    Phase = 90.0
    Frequency = 50.0
    IsSolid = true
    CoilGroup = 'CoilGroup1'
    ConnectionType = 'Series'
    NumberOfTurns = 100
    WaveformType = 'AC'
    DutyCycle = 0.5
    WindingType = 'Solid'
    MotionType = 'Rotation'
    RotationSpeed = 1500.0
    Direction = 1
    PolygonPoints = 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0
$end 'Excitation1'
    )";

    // 解析激励源数据
    bool parse_result = parser_->parse_content(excitation_content);
    ASSERT_TRUE(parse_result) << "激励源数据解析失败";

    auto root = parser_->get_root();
    ASSERT_TRUE(root) << "根节点为空";

    std::cout << "解析成功，块名称: " << root->name << std::endl;

    // 转换激励源数据
    auto excitation = converter_->convertExcitationDirect(root);
    ASSERT_TRUE(excitation) << "激励源转换失败";

    // 验证转换结果
    EXPECT_EQ(excitation->getName(), "CurrentSource");
    EXPECT_EQ(excitation->getType(), ExcitationType::CURRENT_DENSITY);
    EXPECT_DOUBLE_EQ(excitation->getValue(), 10.0);
    EXPECT_DOUBLE_EQ(excitation->getPhase(), 90.0);
    EXPECT_DOUBLE_EQ(excitation->getFrequency(), 50.0);
    EXPECT_TRUE(excitation->isSolid());
    EXPECT_EQ(excitation->getCoilGroup(), "CoilGroup1");
    EXPECT_EQ(excitation->getConnectionType(), CoilConnectionType::SERIES);
    EXPECT_EQ(excitation->getNumberOfTurns(), 100);
    EXPECT_EQ(excitation->getWaveformType(), ExcitationWaveformType::SINUSOIDAL);
    EXPECT_DOUBLE_EQ(excitation->getDutyCycle(), 0.5);
    EXPECT_EQ(excitation->getWindingType(), WindingType::SOLID);
    EXPECT_EQ(excitation->getMotionType(), MotionType::ROTATION);
    EXPECT_DOUBLE_EQ(excitation->getRotationSpeed(), 1500.0);
    EXPECT_EQ(excitation->getDirection(), 1);
    EXPECT_EQ(excitation->getMaxwellExcitationID(), "Excitation1");
    
    std::cout << "激励源转换成功:" << std::endl;
    std::cout << "  名称: " << excitation->getName() << std::endl;
    std::cout << "  类型: Current Density" << std::endl;
    std::cout << "  值: " << excitation->getValue() << " A" << std::endl;
    std::cout << "  频率: " << excitation->getFrequency() << " Hz" << std::endl;
    std::cout << "  匝数: " << excitation->getNumberOfTurns() << std::endl;
    std::cout << "  转速: " << excitation->getRotationSpeed() << " RPM" << std::endl;
    std::cout << "  Maxwell ID: " << excitation->getMaxwellExcitationID() << std::endl;
    
    // 验证多边形点
    const auto& polygon_points = excitation->getPolygonPoints();
    std::cout << "  多边形点数: " << polygon_points.size() << std::endl;
    
    std::cout << "激励源测试通过" << std::endl;
}

/**
 * @brief 测试材料转换功能
 */
TEST_F(MaxwellConverterMainTest, CompleteMaterialConversionTest) {
    std::cout << "=== 测试材料转换功能 ===" << std::endl;
    
    // 测试材料数据
    std::string material_content = R"(
$begin 'Material1'
    Name = 'Copper'
    Type = 'LinearIsotropic'
    RelativePermeability = 1.0
    Conductivity = 5.8e7
    MassDensity = 8960.0
$end 'Material1'
    )";

    // 解析材料数据
    bool parse_result = parser_->parse_content(material_content);
    ASSERT_TRUE(parse_result) << "材料数据解析失败";

    auto root = parser_->get_root();
    ASSERT_TRUE(root) << "根节点为空";

    std::cout << "解析成功，块名称: " << root->name << std::endl;

    // 转换材料数据
    auto material = converter_->convertMaterialDirect(root);
    ASSERT_TRUE(material) << "材料转换失败";

    // 验证转换结果
    EXPECT_EQ(material->getName(), "Copper");
    EXPECT_EQ(material->getType(), MatType::LINEAR_ISOTROPIC);
    EXPECT_DOUBLE_EQ(material->getRelativePermeability(), 1.0);
    EXPECT_DOUBLE_EQ(material->getConductivity(), 5.8e7);
    EXPECT_DOUBLE_EQ(material->getMassDensity(), 8960.0);
    EXPECT_EQ(material->getMaxwellMaterialID(), "Material1");
    
    std::cout << "材料转换成功:" << std::endl;
    std::cout << "  名称: " << material->getName() << std::endl;
    std::cout << "  类型: Linear Isotropic" << std::endl;
    std::cout << "  相对磁导率: " << material->getRelativePermeability() << std::endl;
    std::cout << "  电导率: " << material->getConductivity() << " S/m" << std::endl;
    std::cout << "  质量密度: " << material->getMassDensity() << " kg/m³" << std::endl;
    std::cout << "  Maxwell ID: " << material->getMaxwellMaterialID() << std::endl;
    
    std::cout << "材料测试通过" << std::endl;
}

/**
 * @brief 测试空数据转换
 */
TEST_F(MaxwellConverterMainTest, EmptyConversionTest) {
    std::cout << "=== 测试空数据转换 ===" << std::endl;
    
    // 测试空边界条件
    auto boundary = converter_->convertBoundaryDirect(nullptr);
    EXPECT_FALSE(boundary) << "空边界条件应该返回nullptr";
    
    // 测试空激励源
    auto excitation = converter_->convertExcitationDirect(nullptr);
    EXPECT_FALSE(excitation) << "空激励源应该返回nullptr";
    
    // 测试空材料
    auto material = converter_->convertMaterialDirect(nullptr);
    EXPECT_FALSE(material) << "空材料应该返回nullptr";
    
    std::cout << "空数据转换测试通过" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "开始运行Maxwell数据转换器测试..." << std::endl;
    
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    
    if (result == 0) {
        std::cout << "所有测试通过!" << std::endl;
    } else {
        std::cout << "有测试失败!" << std::endl;
    }
    
    return result;
}