/**
 * @file test_maxwell_converter.cpp
 * @brief Maxwell数据转换器测试
 * @details 测试Maxwell解析数据到内部数据模型的转换功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include "tool/maxwell_converter_impl.hpp"
#include "tool/maxwell_parser.hpp"
#include <memory>

using namespace tool;

class MaxwellConverterTest : public ::testing::Test {
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
 * @brief 测试材料数据转换功能
 */
TEST_F(MaxwellConverterTest, MaterialConversionTest) {
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
    ASSERT_EQ(root->name, "Material1") << "块名称不匹配";

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
}

/**
 * @brief 测试非线性材料转换功能
 */
TEST_F(MaxwellConverterTest, NonlinearMaterialConversionTest) {
    // 测试非线性材料数据
    std::string material_content = R"(
$begin 'Material2'
    Name = 'SiliconSteel'
    Type = 'NonlinearIsotropic'
    RelativePermeability = 2000.0
    Conductivity = 2.0e6
    BHCurve = 0.0 0.0 100.0 0.5 200.0 1.0 300.0 1.2
    CoreLossEnabled = true
    CoreLossModel = 'Steinmetz'
    CoreLossCoefficients = 0.01 1.5 1.8 0.0
$end 'Material2'
    )";

    // 解析材料数据
    bool parse_result = parser_->parse_content(material_content);
    ASSERT_TRUE(parse_result) << "非线性材料数据解析失败";

    auto root = parser_->get_root();
    ASSERT_TRUE(root) << "根节点为空";

    // 转换材料数据
    auto material = converter_->convertMaterialDirect(root);
    ASSERT_TRUE(material) << "非线性材料转换失败";

    // 验证转换结果
    EXPECT_EQ(material->getName(), "SiliconSteel");
    EXPECT_EQ(material->getType(), MatType::NONLINEAR_ISOTROPIC);
    EXPECT_DOUBLE_EQ(material->getRelativePermeability(), 2000.0);
    EXPECT_DOUBLE_EQ(material->getConductivity(), 2.0e6);
    EXPECT_TRUE(material->isCoreLossEnabled());
    EXPECT_EQ(material->getCoreLossModel(), CoreLossModelType::STEINMETZ);

    // 验证B-H曲线数据
    const auto& bh_curve = material->getBHCurve();
    EXPECT_EQ(bh_curve.size(), 4);
    if (bh_curve.size() >= 4) {
        EXPECT_DOUBLE_EQ(bh_curve[0].h, 0.0);
        EXPECT_DOUBLE_EQ(bh_curve[0].b, 0.0);
        EXPECT_DOUBLE_EQ(bh_curve[1].h, 100.0);
        EXPECT_DOUBLE_EQ(bh_curve[1].b, 0.5);
        EXPECT_DOUBLE_EQ(bh_curve[2].h, 200.0);
        EXPECT_DOUBLE_EQ(bh_curve[2].b, 1.0);
        EXPECT_DOUBLE_EQ(bh_curve[3].h, 300.0);
        EXPECT_DOUBLE_EQ(bh_curve[3].b, 1.2);
    }
}

/**
 * @brief 测试各向异性材料转换功能
 */
TEST_F(MaxwellConverterTest, AnisotropicMaterialConversionTest) {
    // 测试各向异性材料数据
    std::string material_content = R"(
$begin 'Material3'
    Name = 'AnisotropicMaterial'
    Type = 'LinearAnisotropic'
    AnisotropicPermeability = 1000.0 500.0 1000.0
    AnisotropicConductivity = 1.0e6 2.0e6 1.0e6
    TemperatureCoefficient = 0.0039
$end 'Material3'
    )";

    // 解析材料数据
    bool parse_result = parser_->parse_content(material_content);
    ASSERT_TRUE(parse_result) << "各向异性材料数据解析失败";

    auto root = parser_->get_root();
    ASSERT_TRUE(root) << "根节点为空";

    // 转换材料数据
    auto material = converter_->convertMaterialDirect(root);
    ASSERT_TRUE(material) << "各向异性材料转换失败";

    // 验证转换结果
    EXPECT_EQ(material->getName(), "AnisotropicMaterial");
    EXPECT_EQ(material->getType(), MatType::LINEAR_ANISOTROPIC);
    EXPECT_DOUBLE_EQ(material->getTemperatureCoefficient(), 0.0039);

    // 验证各向异性参数
    const auto& perm_data = material->getAnisotropicPermeability();
    EXPECT_EQ(perm_data.size(), 3);
    if (perm_data.size() >= 3) {
        EXPECT_DOUBLE_EQ(perm_data[0], 1000.0);
        EXPECT_DOUBLE_EQ(perm_data[1], 500.0);
        EXPECT_DOUBLE_EQ(perm_data[2], 1000.0);
    }

    const auto& cond_data = material->getAnisotropicConductivity();
    EXPECT_EQ(cond_data.size(), 3);
    if (cond_data.size() >= 3) {
        EXPECT_DOUBLE_EQ(cond_data[0], 1.0e6);
        EXPECT_DOUBLE_EQ(cond_data[1], 2.0e6);
        EXPECT_DOUBLE_EQ(cond_data[2], 1.0e6);
    }
}

/**
 * @brief 测试空材料块转换
 */
TEST_F(MaxwellConverterTest, EmptyMaterialConversionTest) {
    // 测试空材料块
    auto material = converter_->convertMaterialDirect(nullptr);
    EXPECT_FALSE(material) << "空材料块应该返回nullptr";
}

/**
 * @brief 测试边界条件转换功能
 */
TEST_F(MaxwellConverterTest, BoundaryConversionTest) {
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
    ASSERT_EQ(root->name, "Boundary1") << "块名称不匹配";

    // 转换边界条件数据
    auto boundary = converter_->convertBoundaryDirect(root);
    ASSERT_TRUE(boundary) << "边界条件转换失败";

    // 验证转换结果
    EXPECT_EQ(boundary->getName(), "ZeroPotential");
    EXPECT_EQ(boundary->getType(), BndType::DIRICHLET);
    //EXPECT_DOUBLE_EQ(boundary->getVectorPotential(), 0.0);
    EXPECT_EQ(boundary->getMaxwellBoundaryID(), "Boundary1");
    
    // 验证关联几何对象
    const auto& faces = boundary->getFaces();
    EXPECT_EQ(faces.size(), 2);
    if (faces.size() >= 2) {
        EXPECT_EQ(faces[0], "Face1");
        EXPECT_EQ(faces[1], "Face2");
    }
    
    const auto& edges = boundary->getEdges();
    EXPECT_EQ(edges.size(), 1);
    if (edges.size() >= 1) {
        EXPECT_EQ(edges[0], "Edge1");
    }
    
    const auto& objects = boundary->getObjects();
    EXPECT_EQ(objects.size(), 1);
    if (objects.size() >= 1) {
        EXPECT_EQ(objects[0], "Object1");
    }
    
    EXPECT_EQ(boundary->getMasterName(), "MasterBoundary");
    EXPECT_EQ(boundary->getSlaveName(), "SlaveBoundary");
    //EXPECT_EQ(boundary->getBoundarySubType(), BoundarySubType::IMPEDANCE);
    EXPECT_DOUBLE_EQ(boundary->getRadiationDistance(), 10.0);
    EXPECT_TRUE(boundary->getPerfectESymmetry());
    EXPECT_FALSE(boundary->getPerfectHSymmetry());
    EXPECT_DOUBLE_EQ(boundary->getInfiniteSphereRadius(), 100.0);
}

/**
 * @brief 测试激励源转换功能
 */
TEST_F(MaxwellConverterTest, ExcitationConversionTest) {
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
    ASSERT_EQ(root->name, "Excitation1") << "块名称不匹配";

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
    //EXPECT_EQ(excitation->getWaveformType(), ExcitationWaveformType::AC);
    EXPECT_DOUBLE_EQ(excitation->getDutyCycle(), 0.5);
    EXPECT_EQ(excitation->getWindingType(), WindingType::SOLID);
    EXPECT_EQ(excitation->getMotionType(), MotionType::ROTATION);
    EXPECT_DOUBLE_EQ(excitation->getRotationSpeed(), 1500.0);
    EXPECT_EQ(excitation->getDirection(), 1);
    EXPECT_EQ(excitation->getMaxwellExcitationID(), "Excitation1");
    
    // 验证多边形点
    const auto& polygon_points = excitation->getPolygonPoints();
    EXPECT_EQ(polygon_points.size(), 4);
    if (polygon_points.size() >= 4) {
        EXPECT_DOUBLE_EQ(polygon_points[0].first, 0.0);
        EXPECT_DOUBLE_EQ(polygon_points[0].second, 0.0);
        EXPECT_DOUBLE_EQ(polygon_points[1].first, 1.0);
        EXPECT_DOUBLE_EQ(polygon_points[1].second, 0.0);
        EXPECT_DOUBLE_EQ(polygon_points[2].first, 1.0);
        EXPECT_DOUBLE_EQ(polygon_points[2].second, 1.0);
        EXPECT_DOUBLE_EQ(polygon_points[3].first, 0.0);
        EXPECT_DOUBLE_EQ(polygon_points[3].second, 1.0);
    }
}

/**
 * @brief 测试属性值解析功能
 */
TEST_F(MaxwellConverterTest, PropertyValueParsingTest) {
    // 测试数值解析
    EXPECT_DOUBLE_EQ(converter_->parseNumericValue("123.45"), 123.45);
    EXPECT_DOUBLE_EQ(converter_->parseNumericValue("-67.89"), -67.89);
    EXPECT_DOUBLE_EQ(converter_->parseNumericValue("1.23e-4"), 1.23e-4);
    
    // 测试布尔值解析
    EXPECT_TRUE(converter_->parseBooleanValue("true"));
    EXPECT_TRUE(converter_->parseBooleanValue("True"));
    EXPECT_TRUE(converter_->parseBooleanValue("1"));
    EXPECT_TRUE(converter_->parseBooleanValue("yes"));
    EXPECT_FALSE(converter_->parseBooleanValue("false"));
    EXPECT_FALSE(converter_->parseBooleanValue("False"));
    EXPECT_FALSE(converter_->parseBooleanValue("0"));
    EXPECT_FALSE(converter_->parseBooleanValue("no"));
    
    // 测试数值数组解析
    auto numeric_array = converter_->parseNumericArray("1.0 2.0 3.0 4.0");
    EXPECT_EQ(numeric_array.size(), 4);
    if (numeric_array.size() >= 4) {
        EXPECT_DOUBLE_EQ(numeric_array[0], 1.0);
        EXPECT_DOUBLE_EQ(numeric_array[1], 2.0);
        EXPECT_DOUBLE_EQ(numeric_array[2], 3.0);
        EXPECT_DOUBLE_EQ(numeric_array[3], 4.0);
    }
    
    // 测试字符串数组解析
    auto string_array = converter_->parseStringArray("'item1' 'item2' 'item3'");
    EXPECT_EQ(string_array.size(), 3);
    if (string_array.size() >= 3) {
        EXPECT_EQ(string_array[0], "item1");
        EXPECT_EQ(string_array[1], "item2");
        EXPECT_EQ(string_array[2], "item3");
    }
}

/**
 * @brief 测试枚举类型转换功能
 */
TEST_F(MaxwellConverterTest, EnumConversionTest) {
    // 测试材料类型转换
    EXPECT_EQ(converter_->convertMaterialType("LinearIsotropic"), MatType::LINEAR_ISOTROPIC);
    EXPECT_EQ(converter_->convertMaterialType("LinearAnisotropic"), MatType::LINEAR_ANISOTROPIC);
    EXPECT_EQ(converter_->convertMaterialType("NonlinearIsotropic"), MatType::NONLINEAR_ISOTROPIC);
    EXPECT_EQ(converter_->convertMaterialType("NonlinearAnisotropic"), MatType::NONLINEAR_ANISOTROPIC);
    
    // 测试B-H曲线类型转换
    EXPECT_EQ(converter_->convertBHCurveType("SingleCurve"), BHCurveType::SINGLE_CURVE);
    EXPECT_EQ(converter_->convertBHCurveType("TempDependent"), BHCurveType::TEMP_DEPENDENT);
    EXPECT_EQ(converter_->convertBHCurveType("FreqDependent"), BHCurveType::FREQ_DEPENDENT);
    EXPECT_EQ(converter_->convertBHCurveType("CustomCurve"), BHCurveType::CUSTOM_CURVE);
    
    // 测试磁芯损耗模型类型转换
    EXPECT_EQ(converter_->convertCoreLossModelType("Steinmetz"), CoreLossModelType::STEINMETZ);
    EXPECT_EQ(converter_->convertCoreLossModelType("Bertotti"), CoreLossModelType::Bertotti);
    EXPECT_EQ(converter_->convertCoreLossModelType("Custom"), CoreLossModelType::CUSTOM);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}