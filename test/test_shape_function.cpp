/**
 * @file test_shape_function.cpp
 * @brief 形函数模块功能测试（Google Test版本）
 * @details 覆盖Lagrange节点元、Nedelec棱边矢量元和工厂类的全部核心接口，
 *          包括形函数值、单位分解性、雅可比矩阵、棱边形函数、旋度计算等。
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "shape_function_base.hpp"
#include "lagrange_element.hpp"
#include "nedelec_element.hpp"
#include "shape_function_factory.hpp"

using namespace numeric;

// ==================== 测试固件 ====================

class ShapeFunctionTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

// ====================================================================
//  1.1 Lagrange 单元测试 - 形函数值与单位分解性
// ====================================================================

TEST_F(ShapeFunctionTest, TRI3_Centroid_ShapeFunctionValues) {
    auto tri3 = std::make_unique<LagrangeElement<2>>(ElementType::TRI3);

    LocalPoint xi(1.0 / 3.0, 1.0 / 3.0);
    Eigen::VectorXd N = tri3->evalN(xi);

    ASSERT_EQ(N.size(), 3);
    EXPECT_NEAR(N[0], 1.0 / 3.0, 1e-12) << "N1在形心处应为1/3";
    EXPECT_NEAR(N[1], 1.0 / 3.0, 1e-12) << "N2在形心处应为1/3";
    EXPECT_NEAR(N[2], 1.0 / 3.0, 1e-12) << "N3在形心处应为1/3";

    EXPECT_NEAR(N.sum(), 1.0, 1e-12) << "TRI3形函数应满足单位分解性 ΣNi=1";
}

TEST_F(ShapeFunctionTest, QUAD4_Center_ShapeFunctionValues) {
    auto quad4 = std::make_unique<LagrangeElement<2>>(ElementType::QUAD4);

    LocalPoint xi(0.0, 0.0);
    Eigen::VectorXd N = quad4->evalN(xi);

    ASSERT_EQ(N.size(), 4);
    EXPECT_NEAR(N[0], 0.25, 1e-12) << "N1在中心处应为0.25";
    EXPECT_NEAR(N[1], 0.25, 1e-12) << "N2在中心处应为0.25";
    EXPECT_NEAR(N[2], 0.25, 1e-12) << "N3在中心处应为0.25";
    EXPECT_NEAR(N[3], 0.25, 1e-12) << "N4在中心处应为0.25";
}

TEST_F(ShapeFunctionTest, TET4_Centroid_ShapeFunctionValues) {
    auto tet4 = std::make_unique<LagrangeElement<3>>(ElementType::TET4);

    LocalPoint xi(0.25, 0.25, 0.25);
    Eigen::VectorXd N = tet4->evalN(xi);

    ASSERT_EQ(N.size(), 4);
    EXPECT_NEAR(N[0], 0.25, 1e-12) << "N1在形心处应为0.25";
    EXPECT_NEAR(N[1], 0.25, 1e-12) << "N2在形心处应为0.25";
    EXPECT_NEAR(N[2], 0.25, 1e-12) << "N3在形心处应为0.25";
    EXPECT_NEAR(N[3], 0.25, 1e-12) << "N4在形心处应为0.25";
}

TEST_F(ShapeFunctionTest, HEX8_Center_ShapeFunctionValues) {
    auto hex8 = std::make_unique<LagrangeElement<3>>(ElementType::HEX8);

    LocalPoint xi(0.0, 0.0, 0.0);
    Eigen::VectorXd N = hex8->evalN(xi);

    ASSERT_EQ(N.size(), 8);
    for (int i = 0; i < 8; i++) {
        EXPECT_NEAR(N[i], 0.125, 1e-12)
            << "HEX8第" << (i + 1) << "个节点形函数在中心处应为0.125";
    }
}

TEST_F(ShapeFunctionTest, PRISM6_Centroid_ShapeFunctionValues) {
    auto prism6 = std::make_unique<LagrangeElement<3>>(ElementType::PRISM6);

    LocalPoint xi(1.0 / 3.0, 1.0 / 3.0, 0.0);
    Eigen::VectorXd N = prism6->evalN(xi);

    ASSERT_EQ(N.size(), 6);
    double expected_val = 1.0 / 6.0;
    for (int i = 0; i < 6; i++) {
        EXPECT_NEAR(N[i], expected_val, 1e-12)
            << "PRISM6第" << (i + 1) << "个节点形函数在几何中心处应为1/6";
    }

    EXPECT_NEAR(N.sum(), 1.0, 1e-12) << "PRISM6应满足单位分解性";
}

TEST_F(ShapeFunctionTest, PYRAMID5_NearBaseCenter_ShapeFunctionValues) {
    auto pyramid5 = std::make_unique<LagrangeElement<3>>(ElementType::PYRAMID5);

    LocalPoint xi(0.0, 0.0, 0.0);
    Eigen::VectorXd N = pyramid5->evalN(xi);

    ASSERT_EQ(N.size(), 5);

    // 在底面中心(ξ=0,η=0,ζ=0)：底面四个角点各贡献0.25，锥顶贡献0
    EXPECT_NEAR(N[0], 0.25, 1e-12) << "底面角点1在底面中心处应为0.25";
    EXPECT_NEAR(N[1], 0.25, 1e-12) << "底面角点2在底面中心处应为0.25";
    EXPECT_NEAR(N[2], 0.25, 1e-12) << "底面角点3在底面中心处应为0.25";
    EXPECT_NEAR(N[3], 0.25, 1e-12) << "底面角点4在底面中心处应为0.25";
    EXPECT_NEAR(N[4], 0.0, 1e-12) << "锥顶在底面中心处应为0";

    EXPECT_NEAR(N.sum(), 1.0, 1e-12) << "PYRAMID5应满足单位分解性";
}

TEST_F(ShapeFunctionTest, UnitPartitionProperty) {
    struct TestCase {
        ElementType type;
        int dim;
        std::vector<LocalPoint> points;
    };

    std::vector<TestCase> test_cases;

    // TRI3：面积坐标系下的多个积分点
    test_cases.push_back({ElementType::TRI3, 2,
                          {LocalPoint(0.1, 0.2),
                           LocalPoint(0.5, 0.3),
                           LocalPoint(0.33, 0.33)}});

    // QUAD4：参考域[-1,1]²内的多个点
    test_cases.push_back({ElementType::QUAD4, 2,
                          {LocalPoint(0.0, 0.0),
                           LocalPoint(0.5, -0.3),
                           LocalPoint(-0.7, 0.7)}});

    // TET4：体积坐标系内的多个点
    test_cases.push_back({ElementType::TET4, 3,
                          {LocalPoint(0.1, 0.2, 0.3),
                           LocalPoint(0.25, 0.25, 0.25),
                           LocalPoint(0.05, 0.15, 0.1)}});

    // HEX8：参考域[-1,1]³内的多个点
    test_cases.push_back({ElementType::HEX8, 3,
                          {LocalPoint(0.0, 0.0, 0.0),
                           LocalPoint(0.5, -0.3, 0.7),
                           LocalPoint(-0.5, 0.5, -0.5)}});

    // PRISM6
    test_cases.push_back({ElementType::PRISM6, 3,
                          {LocalPoint(0.2, 0.3, 0.0),
                           LocalPoint(0.33, 0.33, 0.5),
                           LocalPoint(0.1, 0.1, -0.7)}});

    // PYRAMID5（测试点均满足ζ<1，避开锥顶奇异性）
    test_cases.push_back({ElementType::PYRAMID5, 3,
                          {LocalPoint(0.0, 0.0, 0.0),
                           LocalPoint(0.3, -0.2, 0.3),
                           LocalPoint(-0.5, 0.5, 0.1)}});

    for (const auto& tc : test_cases) {
        std::unique_ptr<ShapeFunctionBase> elem;
        if (tc.dim == 2) {
            elem = std::make_unique<LagrangeElement<2>>(tc.type);
        } else {
            elem = std::make_unique<LagrangeElement<3>>(tc.type);
        }

        for (const auto& pt : tc.points) {
            Eigen::VectorXd N = elem->evalN(pt);
            EXPECT_NEAR(N.sum(), 1.0, 1e-10)
                << "单元类型=" << static_cast<int>(tc.type)
                << " 在点(" << pt.coords[0] << "," << pt.coords[1]
                << "," << pt.coords[2] << ")处不满足单位分解性, sum=" << N.sum();
        }
    }
}

TEST_F(ShapeFunctionTest, JacobianDeterminant_Positivity) {
    // QUAD4单元：使用单位正方形物理坐标
    auto quad4 = std::make_unique<LagrangeElement<2>>(ElementType::QUAD4);

    // 节点坐标矩阵（dim × node_count）：单位正方形
    Eigen::MatrixXd quad4_nodes(2, 4);
    quad4_nodes.col(0) = Eigen::Vector2d(0.0, 0.0);   // (-1,-1)→(0,0)
    quad4_nodes.col(1) = Eigen::Vector2d(1.0, 0.0);   // ( 1,-1)→(1,0)
    quad4_nodes.col(2) = Eigen::Vector2d(1.0, 1.0);   // ( 1, 1)→(1,1)
    quad4_nodes.col(3) = Eigen::Vector2d(0.0, 1.0);   // (-1, 1)→(0,1)

    LocalPoint xi_quad(0.0, 0.0);
    JacobianResult jaco_quad = quad4->calcJacobian(xi_quad, quad4_nodes);
    EXPECT_GT(jaco_quad.det_j, 0.0)
        << "QUAD4雅可比行列式应在正坐标下为正值, detJ=" << jaco_quad.det_j;

    // HEX8单元：使用单位立方体物理坐标
    auto hex8 = std::make_unique<LagrangeElement<3>>(ElementType::HEX8);

    Eigen::MatrixXd hex8_nodes(3, 8);
    hex8_nodes.col(0) = Eigen::Vector3d(0.0, 0.0, 0.0);
    hex8_nodes.col(1) = Eigen::Vector3d(1.0, 0.0, 0.0);
    hex8_nodes.col(2) = Eigen::Vector3d(1.0, 1.0, 0.0);
    hex8_nodes.col(3) = Eigen::Vector3d(0.0, 1.0, 0.0);
    hex8_nodes.col(4) = Eigen::Vector3d(0.0, 0.0, 1.0);
    hex8_nodes.col(5) = Eigen::Vector3d(1.0, 0.0, 1.0);
    hex8_nodes.col(6) = Eigen::Vector3d(1.0, 1.0, 1.0);
    hex8_nodes.col(7) = Eigen::Vector3d(0.0, 1.0, 1.0);

    LocalPoint xi_hex(0.0, 0.0, 0.0);
    JacobianResult jaco_hex = hex8->calcJacobian(xi_hex, hex8_nodes);
    EXPECT_GT(jaco_hex.det_j, 0.0)
        << "HEX8雅可比行列式应在正坐标下为正值, detJ=" << jaco_hex.det_j;

    // TET4单元：使用非退化四面体物理坐标（节点顺序遵循右手定则）
    auto tet4 = std::make_unique<LagrangeElement<3>>(ElementType::TET4);

    Eigen::MatrixXd tet4_nodes(3, 4);
    tet4_nodes.col(0) = Eigen::Vector3d(0.0, 0.0, 0.0);
    tet4_nodes.col(1) = Eigen::Vector3d(1.0, 0.0, 0.0);
    tet4_nodes.col(2) = Eigen::Vector3d(0.0, 0.0, 1.0);
    tet4_nodes.col(3) = Eigen::Vector3d(0.0, 1.0, 0.0);

    LocalPoint xi_tet(0.2, 0.2, 0.2);
    JacobianResult jaco_tet = tet4->calcJacobian(xi_tet, tet4_nodes);
    EXPECT_GT(jaco_tet.det_j, 0.0)
        << "TET4雅可比行列式应在正坐标下为正值, detJ=" << jaco_tet.det_j;
}

// ====================================================================
//  1.2 Nedelec 棱边元测试
// ====================================================================

TEST_F(ShapeFunctionTest, TRI3_EDGE_TangentialComponent) {
    auto tri3_edge = std::make_unique<NedelecElement<2>>(ElementType::TRI3_EDGE);

    // 在三角形内部某点计算第0条棱边的基函数
    // 棱边0对应节点1→2方向，基函数W₀=(η, -ξ)
    LocalPoint xi(0.4, 0.3);
    Eigen::Vector3d W = tri3_edge->evalEdgeFunction(0, xi);

    // W₀=(η, -ξ)=(0.3, -0.4)，沿棱边切向应有非零分量
    EXPECT_NE(W.norm(), 0.0) << "TRI3_EDGE棱边0基函数不应为零向量";

    // 验证基函数的具体分量值
    EXPECT_NEAR(W[0], 0.3, 1e-14) << "W₀_x分量应为η=0.3";
    EXPECT_NEAR(W[1], -0.4, 1e-14) << "W₀_y分量应为-ξ=-0.4";
    EXPECT_NEAR(W[2], 0.0, 1e-14) << "TRI3_EDGE为2D单元，z分量应为0";

    // 棱边0的方向大致为+x方向（从节点1到节点2），x分量应非零
    EXPECT_NE(std::abs(W[0]), 0.0) << "沿棱边0切向分量应非零";
}

TEST_F(ShapeFunctionTest, TET4_EDGE_DOFCount) {
    auto tet4_edge = std::make_unique<NedelecElement<3>>(ElementType::TET4_EDGE);

    EXPECT_EQ(tet4_edge->getNodeCount(), 6)
        << "TET4_EDGE应有6条棱边自由度";
    EXPECT_EQ(tet4_edge->getNodeType(), ElementType::TET4_EDGE)
        << "节点类型应为TET4_EDGE";
    EXPECT_EQ(tet4_edge->getDim(), 3)
        << "空间维度应为3";
}

TEST_F(ShapeFunctionTest, PRISM6_EDGE_DOFCount) {
    auto prism6_edge = std::make_unique<NedelecElement<3>>(ElementType::PRISM6_EDGE);

    EXPECT_EQ(prism6_edge->getNodeCount(), 9)
        << "PRISM6_EDGE应有9条棱边自由度";
    EXPECT_EQ(prism6_edge->getNodeType(), ElementType::PRISM6_EDGE)
        << "节点类型应为PRISM6_EDGE";
    EXPECT_EQ(prism6_edge->getDim(), 3)
        << "空间维度应为3";
}

TEST_F(ShapeFunctionTest, PYRAMID5_EDGE_DOFCount) {
    auto pyramid5_edge = std::make_unique<NedelecElement<3>>(ElementType::PYRAMID5_EDGE);

    EXPECT_EQ(pyramid5_edge->getNodeCount(), 8)
        << "PYRAMID5_EDGE应有8条棱边自由度";
    EXPECT_EQ(pyramid5_edge->getNodeType(), ElementType::PYRAMID5_EDGE)
        << "节点类型应为PYRAMID5_EDGE";
    EXPECT_EQ(pyramid5_edge->getDim(), 3)
        << "空间维度应为3";
}

TEST_F(ShapeFunctionTest, Nedelec_CurlCalculation) {
    // TRI3_EDGE旋度：一阶Nedelec元的旋度为常标量-2（存入z分量）
    auto tri3_edge = std::make_unique<NedelecElement<2>>(ElementType::TRI3_EDGE);
    LocalPoint xi_tri(0.3, 0.4);

    Eigen::Vector3d curl_tri0 = tri3_edge->evalCurlEdge(0, xi_tri);
    Eigen::Vector3d curl_tri1 = tri3_edge->evalCurlEdge(1, xi_tri);
    Eigen::Vector3d curl_tri2 = tri3_edge->evalCurlEdge(2, xi_tri);

    // TRI3_EDGE所有棱边的旋度均为常数(0, 0, -2)
    EXPECT_NEAR(curl_tri0[0], 0.0, 1e-14) << "TRI3_EDGE curl x分量应为0";
    EXPECT_NEAR(curl_tri0[1], 0.0, 1e-14) << "TRI3_EDGE curl y分量应为0";
    EXPECT_NEAR(curl_tri0[2], -2.0, 1e-14) << "TRI3_EDGE curl z分量应为-2";

    EXPECT_NEAR(curl_tri1[2], -2.0, 1e-14) << "TRI3_EDGE棱边1旋度z分量应为-2";
    EXPECT_NEAR(curl_tri2[2], -2.0, 1e-14) << "TRI3_EDGE棱边2旋度z分量应为-2";

    // TET4_EDGE旋度：为常矢量（不依赖位置），每条棱边不同
    auto tet4_edge = std::make_unique<NedelecElement<3>>(ElementType::TET4_EDGE);
    LocalPoint xi_tet(0.2, 0.2, 0.2);

    // 棱边0的旋度应为(0, 0, -2)
    Eigen::Vector3d curl_tet0 = tet4_edge->evalCurlEdge(0, xi_tet);
    EXPECT_NEAR(curl_tet0[0], 0.0, 1e-14) << "TET4_EDGE棱边0 curl x应为0";
    EXPECT_NEAR(curl_tet0[1], 0.0, 1e-14) << "TET4_EDGE棱边0 curl y应为0";
    EXPECT_NEAR(curl_tet0[2], -2.0, 1e-14) << "TET4_EDGE棱边0 curl z应为-2";

    // 棱边1的旋度应为(-2, 0, 0)
    Eigen::Vector3d curl_tet1 = tet4_edge->evalCurlEdge(1, xi_tet);
    EXPECT_NEAR(curl_tet1[0], -2.0, 1e-14) << "TET4_EDGE棱边1 curl x应为-2";
    EXPECT_NEAR(curl_tet1[1], 0.0, 1e-14) << "TET4_EDGE棱边1 curl y应为0";
    EXPECT_NEAR(curl_tet1[2], 0.0, 1e-14) << "TET4_EDGE棱边1 curl z应为0";

    // 验证旋度为常矢量（在不同点处结果相同）
    LocalPoint xi_tet2(0.1, 0.5, 0.1);
    Eigen::Vector3d curl_tet0_2 = tet4_edge->evalCurlEdge(0, xi_tet2);
    EXPECT_NEAR((curl_tet0 - curl_tet0_2).norm(), 0.0, 1e-14)
        << "TET4_EDGE旋度应为常矢量，不依赖空间位置";
}

// ====================================================================
//  1.3 工厂类测试
// ====================================================================

TEST_F(ShapeFunctionTest, Factory_CreateAllTypes) {
    std::vector<std::string> all_types = ShapeFunctionFactory::getSupportedTypes();

    EXPECT_EQ(static_cast<int>(all_types.size()), 22)
        << "工厂类应支持22种单元类型";

    std::vector<int> expected_node_counts = {
        2, 3,           // LINE2, LINE3
        3, 6, 4, 8, 9,  // TRI3, TRI6, QUAD4, QUAD8, QUAD9
        4, 10, 8, 20, 27, 6, 15, 5, 13,  // TET4..PYRAMID13
        3, 4,           // TRI3_EDGE, QUAD4_EDGE
        6, 12, 9, 8     // TET4_EDGE..PYRAMID5_EDGE
    };

    for (size_t i = 0; i < all_types.size(); i++) {
        const std::string& type_name = all_types[i];
        auto elem = ShapeFunctionFactory::create(type_name);

        ASSERT_NE(elem, nullptr)
            << "工厂创建[" << type_name << "]失败，返回了空指针";

        EXPECT_EQ(elem->getNodeCount(), expected_node_counts[i])
            << "[" << type_name << "]的节点数应为" << expected_node_counts[i]
            << ", 实际得到" << elem->getNodeCount();
    }
}

TEST_F(ShapeFunctionTest, Factory_UnknownType_ReturnsNullptr) {
    auto unknown = ShapeFunctionFactory::create("UNKNOWN_TYPE");
    EXPECT_EQ(unknown, nullptr)
        << "未知类型'UNKNOWN_TYPE'应返回nullptr";

    auto empty_str = ShapeFunctionFactory::create("");
    EXPECT_EQ(empty_str, nullptr)
        << "空字符串应返回nullptr";

    auto partial_match = ShapeFunctionFactory::create("TRI");
    EXPECT_EQ(partial_match, nullptr)
        << "部分匹配'TRI'应返回nullptr";
}

TEST_F(ShapeFunctionTest, Factory_IsSupported) {
    EXPECT_TRUE(ShapeFunctionFactory::isSupported("TRI3"))
        << "'TRI3'应被支持";
    EXPECT_TRUE(ShapeFunctionFactory::isSupported("HEX8"))
        << "'HEX8'应被支持";
    EXPECT_TRUE(ShapeFunctionFactory::isSupported("TET4_EDGE"))
        << "'TET4_EDGE'应被支持";
    EXPECT_TRUE(ShapeFunctionFactory::isSupported("PYRAMID5_EDGE"))
        << "'PYRAMID5_EDGE'应被支持";

    EXPECT_FALSE(ShapeFunctionFactory::isSupported("UNKNOWN"))
        << "'UNKNOWN'不应被支持";
    EXPECT_FALSE(ShapeFunctionFactory::isSupported("Tri3"))
        << "'Tri3'(大小写敏感)不应被支持";
    EXPECT_FALSE(ShapeFunctionFactory::isSupported(""))
        << "空字符串不应被支持";
}

TEST_F(ShapeFunctionTest, Factory_GetSupportedTypes) {
    std::vector<std::string> types = ShapeFunctionFactory::getSupportedTypes();

    EXPECT_EQ(static_cast<int>(types.size()), 22)
        << "支持的类型数量应为22";

    // 验证关键类型名称的存在性
    auto has_type = [&](const std::string& name) -> bool {
        return std::find(types.begin(), types.end(), name) != types.end();
    };

    EXPECT_TRUE(has_type("LINE2")) << "应包含LINE2";
    EXPECT_TRUE(has_type("LINE3")) << "应包含LINE3";
    EXPECT_TRUE(has_type("TRI3")) << "应包含TRI3";
    EXPECT_TRUE(has_type("TRI6")) << "应包含TRI6";
    EXPECT_TRUE(has_type("QUAD4")) << "应包含QUAD4";
    EXPECT_TRUE(has_type("QUAD8")) << "应包含QUAD8";
    EXPECT_TRUE(has_type("QUAD9")) << "应包含QUAD9";
    EXPECT_TRUE(has_type("TET4")) << "应包含TET4";
    EXPECT_TRUE(has_type("TET10")) << "应包含TET10";
    EXPECT_TRUE(has_type("HEX8")) << "应包含HEX8";
    EXPECT_TRUE(has_type("HEX20")) << "应包含HEX20";
    EXPECT_TRUE(has_type("HEX27")) << "应包含HEX27";
    EXPECT_TRUE(has_type("PRISM6")) << "应包含PRISM6";
    EXPECT_TRUE(has_type("PRISM15")) << "应包含PRISM15";
    EXPECT_TRUE(has_type("PYRAMID5")) << "应包含PYRAMID5";
    EXPECT_TRUE(has_type("PYRAMID13")) << "应包含PYRAMID13";
    EXPECT_TRUE(has_type("TRI3_EDGE")) << "应包含TRI3_EDGE";
    EXPECT_TRUE(has_type("QUAD4_EDGE")) << "应包含QUAD4_EDGE";
    EXPECT_TRUE(has_type("TET4_EDGE")) << "应包含TET4_EDGE";
    EXPECT_TRUE(has_type("HEX8_EDGE")) << "应包含HEX8_EDGE";
    EXPECT_TRUE(has_type("PRISM6_EDGE")) << "应包含PRISM6_EDGE";
    EXPECT_TRUE(has_type("PYRAMID5_EDGE")) << "应包含PYRAMID5_EDGE";
}

// ==================== 主函数 ====================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
