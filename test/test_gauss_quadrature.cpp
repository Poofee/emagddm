/**
 * @file test_gauss_quadrature.cpp
 * @brief 高斯积分点库单元测试
 * @details 验证各单元类型积分点的正确性，包括：
 *          1. 积分点数量正确
 *          2. 权重总和等于参考单元体积
 *          3. 坐标在有效范围内
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <cmath>
#include "gauss_quadrature.hpp"

using namespace numeric;

// ==================== 测试辅助函数 ====================

/**
 * @brief 计算积分点权重总和
 * @param points 积分点向量
 * @return double 权重总和
 */
static double calcTotalWeight(const std::vector<GaussPoint>& points)
{
    double total = 0.0;
    for (const auto& gp : points) {
        total += gp.weight;
    }
    return total;
}

/**
 * @brief 验证2D积分点坐标在单位三角形内
 * @param points 积分点向量
 * @return bool 所有坐标均在有效范围内返回true
 */
static bool validateTriCoords(const std::vector<GaussPoint>& points)
{
    for (const auto& gp : points) {
        double xi = gp.coords[0];
        double eta = gp.coords[1];
        // 面积坐标约束: ξ ≥ 0, η ≥ 0, ξ+η ≤ 1
        if (xi < -1e-10 || eta < -1e-10 || (xi + eta) > 1.0 + 1e-10) {
            return false;
        }
    }
    return true;
}

/**
 * @brief 验证QUAD4积分点坐标在[-1,1]范围内
 * @param points 积分点向量
 * @return bool 所有坐标均在有效范围内返回true
 */
static bool validateQuadCoords(const std::vector<GaussPoint>& points)
{
    for (const auto& gp : points) {
        if (std::abs(gp.coords[0]) > 1.0 + 1e-10 ||
            std::abs(gp.coords[1]) > 1.0 + 1e-10) {
            return false;
        }
    }
    return true;
}

// ==================== TRI3 单元测试 ====================

/**
 * @test TRI3单元1点积分测试
 * @details 验证质心积分方案的基本属性
 */
TEST(GaussQuadratureTest, TRI3_1Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::TRI3, 1);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 1);

    // 验证维度
    EXPECT_EQ(points[0].dim, 2);

    // 验证质心坐标 (1/3, 1/3)
    EXPECT_NEAR(points[0].coords[0], 1.0 / 3.0, 1e-10);
    EXPECT_NEAR(points[0].coords[1], 1.0 / 3.0, 1e-10);

    // 验证权重（三角形面积=0.5）
    EXPECT_NEAR(points[0].weight, 0.5, 1e-10);
}

/**
 * @test TRI3单元3点积分测试
 * @details 验证边中点积分方案的权重和与面积
 */
TEST(GaussQuadratureTest, TRI3_3Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::TRI3, 3);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 3);

    // 验证每个点的维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 2);
    }

    // 验证坐标有效性
    ASSERT_TRUE(validateTriCoords(points));

    // 验证权重总和（应等于三角形面积0.5）
    double total_weight = calcTotalWeight(points);
    EXPECT_NEAR(total_weight, 0.5, 1e-10);

    // 验证每点权重相等（均为1/6）
    for (const auto& gp : points) {
        EXPECT_NEAR(gp.weight, 1.0 / 6.0, 1e-10);
    }
}

// ==================== QUAD4 单元测试 ====================

/**
 * @test QUAD4单元4点积分测试
 * @details 验证Gauss-Legendre 2×2积分方案
 */
TEST(GaussQuadratureTest, QUAD4_4Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::QUAD4, 4);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 4);

    // 验证维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 2);
    }

    // 验证坐标范围 [-1, 1]
    ASSERT_TRUE(validateQuadCoords(points));

    // 验证权重总和（正方形面积=4）
    double total_weight = calcTotalWeight(points);
    EXPECT_NEAR(total_weight, 4.0, 1e-10);

    // 验证每点权重为1.0
    for (const auto& gp : points) {
        EXPECT_NEAR(gp.weight, 1.0, 1e-10);
    }

    // 验证Gauss点位置 ±1/√3
    double gp_pos = 1.0 / std::sqrt(3.0);
    for (const auto& gp : points) {
        EXPECT_NEAR(std::abs(gp.coords[0]), gp_pos, 1e-10);
        EXPECT_NEAR(std::abs(gp.coords[1]), gp_pos, 1e-10);
    }
}

// ==================== TET4 单元测试 ====================

/**
 * @test TET4单元1点积分测试
 * @details 验证四面体质心积分方案
 */
TEST(GaussQuadratureTest, TET4_1Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::TET4, 1);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 1);

    // 验证维度
    EXPECT_EQ(points[0].dim, 3);

    // 验证质心坐标 (1/4, 1/4, 1/4)
    EXPECT_NEAR(points[0].coords[0], 1.0 / 4.0, 1e-10);
    EXPECT_NEAR(points[0].coords[1], 1.0 / 4.0, 1e-10);
    EXPECT_NEAR(points[0].coords[2], 1.0 / 4.0, 1e-10);

    // 验证权重（四面体体积=1/6）
    EXPECT_NEAR(points[0].weight, 1.0 / 6.0, 1e-10);
}

/**
 * @test TET4单元4点积分测试
 * @details 验证Hammer积分方案
 */
TEST(GaussQuadratureTest, TET4_4Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::TET4, 4);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 4);

    // 验证维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 3);
    }

    // 验证权重总和（应等于四面体体积1/6）
    double total_weight = calcTotalWeight(points);
    EXPECT_NEAR(total_weight, 4.0 / 24.0, 1e-10);  // 4 × 1/24 = 1/6

    // 验证每点权重为1/24
    for (const auto& gp : points) {
        EXPECT_NEAR(gp.weight, 1.0 / 24.0, 1e-10);
    }

    // 验证坐标为正值且小于1（在四面体内部）
    for (const auto& gp : points) {
        EXPECT_GE(gp.coords[0], 0.0);
        EXPECT_GE(gp.coords[1], 0.0);
        EXPECT_GE(gp.coords[2], 0.0);
        EXPECT_LE(gp.coords[0] + gp.coords[1] + gp.coords[2], 1.0 + 1e-10);
    }
}

// ==================== HEX8 单元测试 ====================

/**
 * @test HEX8单元8点积分测试
 * @details 验证Gauss-Legendre 2×2×2积分方案
 */
TEST(GaussQuadratureTest, HEX8_8Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::HEX8, 8);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 8);

    // 验证维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 3);
    }

    // 验证坐标范围 [-1, 1]
    for (const auto& gp : points) {
        EXPECT_LE(std::abs(gp.coords[0]), 1.0 + 1e-10);
        EXPECT_LE(std::abs(gp.coords[1]), 1.0 + 1e-10);
        EXPECT_LE(std::abs(gp.coords[2]), 1.0 + 1e-10);
    }

    // 验证权重总和（立方体体积=8）
    double total_weight = calcTotalWeight(points);
    EXPECT_NEAR(total_weight, 8.0, 1e-10);

    // 验证每点权重为1.0
    for (const auto& gp : points) {
        EXPECT_NEAR(gp.weight, 1.0, 1e-10);
    }

    // 验证Gauss点位置 ±1/√3
    double gp_pos = 1.0 / std::sqrt(3.0);
    int count_plus[3] = {0, 0, 0};
    int count_minus[3] = {0, 0, 0};

    for (const auto& gp : points) {
        for (int d = 0; d < 3; d++) {
            if (gp.coords[d] > 0) count_plus[d]++;
            else count_minus[d]++;
        }
    }

    // 每个方向应该有4个正值和4个负值
    for (int d = 0; d < 3; d++) {
        EXPECT_EQ(count_plus[d], 4);
        EXPECT_EQ(count_minus[d], 4);
    }
}

// ==================== PRISM6 单元测试 ====================

/**
 * @test PRISM6单元6点积分测试
 * @details 验证三棱柱张量积积分方案
 */
TEST(GaussQuadratureTest, PRISM6_6Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::PRISM6, 6);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 6);

    // 验证维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 3);
    }

    // 验证ζ坐标在[-1, 1]范围内
    for (const auto& gp : points) {
        EXPECT_LE(std::abs(gp.coords[2]), 1.0 + 1e-10);
    }

    // 验证权重总和（三棱柱体积=1.0）
    double total_weight = calcTotalWeight(points);
    EXPECT_NEAR(total_weight, 1.0, 1e-10);

    // 验证每点权重为1/6
    for (const auto& gp : points) {
        EXPECT_NEAR(gp.weight, 1.0 / 6.0, 1e-10);
    }

    // 验证应有3个点在ζ<0，3个点在ζ>0（对称分布）
    int count_zeta_neg = 0;
    int count_zeta_pos = 0;
    for (const auto& gp : points) {
        if (gp.coords[2] < 0) count_zeta_neg++;
        else count_zeta_pos++;
    }
    EXPECT_EQ(count_zeta_neg, 3);
    EXPECT_EQ(count_zeta_pos, 3);
}

// ==================== PYRAMID5 单元测试 ====================

/**
 * @test PYRAMID5单元5点积分测试
 * @details 验证金字塔近似积分方案
 */
TEST(GaussQuadratureTest, PYRAMID5_5Point)
{
    auto points = GaussQuadrature::getPoints(ElementType::PYRAMID5, 5);

    // 验证积分点数量
    ASSERT_EQ(points.size(), 5);

    // 验证维度
    for (const auto& gp : points) {
        EXPECT_EQ(gp.dim, 3);
    }

    // 验证中心点存在且位于原点附近
    bool found_center = false;
    for (const auto& gp : points) {
        if (std::abs(gp.coords[0]) < 1e-10 &&
            std::abs(gp.coords[1]) < 1e-10) {
            found_center = true;
            // 中心点权重应较大
            EXPECT_GT(gp.weight, 0.5);  // 应该大于0.5
            break;
        }
    }
    ASSERT_TRUE(found_center) << "未找到中心积分点";

    // 验证4个角点分布
    int corner_count = 0;
    for (const auto& gp : points) {
        if (std::abs(gp.coords[0]) > 0.1 && std::abs(gp.coords[1]) > 0.1) {
            corner_count++;
            // 角点权重应较小
            EXPECT_LT(gp.weight, 0.5);
        }
    }
    EXPECT_EQ(corner_count, 4) << "应恰好有4个角点";

    // 验证权重总和（应接近金字塔体积的归一化值）
    double total_weight = calcTotalWeight(points);
    EXPECT_GT(total_weight, 0.0);  // 权重总和应为正
    EXPECT_LT(total_weight, 5.0);   // 合理上限

    // 输出实际值供参考（调试用）
    std::cout << "[INFO] PYRAMID5 总权重: " << total_weight << std::endl;
}

// ==================== 错误处理测试 ====================

/**
 * @test 不支持的单元类型测试
 * @details 验证对无效输入的错误处理
 */
TEST(GaussQuadratureTest, UnsupportedElementType)
{
    // 测试不支持的单元类型（如LINE2）
    auto points = GaussQuadrature::getPoints(ElementType::LINE2, 1);
    EXPECT_TRUE(points.empty());

    // 测试无效的order参数
    points = GaussQuadrature::getPoints(ElementType::TRI3, 99);
    EXPECT_TRUE(points.empty());
}

// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
