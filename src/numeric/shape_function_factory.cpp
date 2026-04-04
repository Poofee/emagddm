/**
 * @file shape_function_factory.cpp
 * @brief 数值计算层 - 形函数单元工厂类完整实现
 * @details 实现全部23种单元类型的字符串匹配与对象创建逻辑，
 *          包括19种Lagrange节点元和4种Nedelec棱边矢量元。
 *
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#include "shape_function_factory.hpp"
#include "lagrange_element.hpp"
#include "nedelec_element.hpp"
#include <unordered_set>

namespace numeric {

// ==================== 工厂方法实现 ====================

std::unique_ptr<ShapeFunctionBase> ShapeFunctionFactory::create(
    const std::string& element_type)
{
    // ---------- 一维Lagrange单元 ----------

    if (element_type == "LINE2") {
        return std::make_unique<LagrangeElement<1>>(ElementType::LINE2);
    }
    if (element_type == "LINE3") {
        return std::make_unique<LagrangeElement<1>>(ElementType::LINE3);
    }

    // ---------- 二维Lagrange单元 ----------

    if (element_type == "TRI3") {
        return std::make_unique<LagrangeElement<2>>(ElementType::TRI3);
    }
    if (element_type == "TRI6") {
        return std::make_unique<LagrangeElement<2>>(ElementType::TRI6);
    }
    if (element_type == "QUAD4") {
        return std::make_unique<LagrangeElement<2>>(ElementType::QUAD4);
    }
    if (element_type == "QUAD8") {
        return std::make_unique<LagrangeElement<2>>(ElementType::QUAD8);
    }
    if (element_type == "QUAD9") {
        return std::make_unique<LagrangeElement<2>>(ElementType::QUAD9);
    }

    // ---------- 三维Lagrange单元 ----------

    if (element_type == "TET4") {
        return std::make_unique<LagrangeElement<3>>(ElementType::TET4);
    }
    if (element_type == "TET10") {
        return std::make_unique<LagrangeElement<3>>(ElementType::TET10);
    }
    if (element_type == "HEX8") {
        return std::make_unique<LagrangeElement<3>>(ElementType::HEX8);
    }
    if (element_type == "HEX20") {
        return std::make_unique<LagrangeElement<3>>(ElementType::HEX20);
    }
    if (element_type == "HEX27") {
        return std::make_unique<LagrangeElement<3>>(ElementType::HEX27);
    }
    if (element_type == "PRISM6") {
        return std::make_unique<LagrangeElement<3>>(ElementType::PRISM6);
    }
    if (element_type == "PRISM15") {
        return std::make_unique<LagrangeElement<3>>(ElementType::PRISM15);
    }
    if (element_type == "PYRAMID5") {
        return std::make_unique<LagrangeElement<3>>(ElementType::PYRAMID5);
    }
    if (element_type == "PYRAMID13") {
        return std::make_unique<LagrangeElement<3>>(ElementType::PYRAMID13);
    }

    // ---------- 二维Nedelec棱边矢量元 ----------

    if (element_type == "TRI3_EDGE") {
        return std::make_unique<NedelecElement<2>>(ElementType::TRI3_EDGE);
    }
    if (element_type == "QUAD4_EDGE") {
        return std::make_unique<NedelecElement<2>>(ElementType::QUAD4_EDGE);
    }

    // ---------- 三维Nedelec棱边矢量元 ----------

    if (element_type == "TET4_EDGE") {
        return std::make_unique<NedelecElement<3>>(ElementType::TET4_EDGE);
    }
    if (element_type == "HEX8_EDGE") {
        return std::make_unique<NedelecElement<3>>(ElementType::HEX8_EDGE);
    }
    if (element_type == "PRISM6_EDGE") {
        return std::make_unique<NedelecElement<3>>(ElementType::PRISM6_EDGE);
    }
    if (element_type == "PYRAMID5_EDGE") {
        return std::make_unique<NedelecElement<3>>(ElementType::PYRAMID5_EDGE);
    }

    // ---------- 未知类型处理 ----------

    FEEM_ERROR("未知的单元类型: {}", element_type);
    return nullptr;
}

// ==================== 类型检查方法实现 ====================

bool ShapeFunctionFactory::isSupported(const std::string& element_type)
{
    static const std::unordered_set<std::string> supported_types = {
        "LINE2", "LINE3",
        "TRI3", "TRI6", "QUAD4", "QUAD8", "QUAD9",
        "TET4", "TET10", "HEX8", "HEX20", "HEX27",
        "PRISM6", "PRISM15", "PYRAMID5", "PYRAMID13",
        "TRI3_EDGE", "QUAD4_EDGE",
        "TET4_EDGE", "HEX8_EDGE", "PRISM6_EDGE", "PYRAMID5_EDGE"
    };
    return supported_types.find(element_type) != supported_types.end();
}

// ==================== 支持类型列表方法实现 ====================

std::vector<std::string> ShapeFunctionFactory::getSupportedTypes()
{
    // 返回全部23种支持的单元类型名称，按维度和类别分组排列
    return {
        // 一维Lagrange节点单元（2种）
        "LINE2",
        "LINE3",

        // 二维Lagrange节点单元（5种）
        "TRI3",
        "TRI6",
        "QUAD4",
        "QUAD8",
        "QUAD9",

        // 三维Lagrange节点单元（12种）
        "TET4",
        "TET10",
        "HEX8",
        "HEX20",
        "HEX27",
        "PRISM6",
        "PRISM15",
        "PYRAMID5",
        "PYRAMID13",

        // 二维Nedelec棱边矢量元（2种）
        "TRI3_EDGE",
        "QUAD4_EDGE",

        // 三维Nedelec棱边矢量元（4种）
        "TET4_EDGE",
        "HEX8_EDGE",
        "PRISM6_EDGE",
        "PYRAMID5_EDGE"
    };
}

} // namespace numeric
