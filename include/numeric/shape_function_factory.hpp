/**
 * @file shape_function_factory.hpp
 * @brief 数值计算层 - 形函数单元工厂类声明
 * @details 提供统一的形函数对象创建接口，根据字符串标识符实例化对应的
 *          Lagrange节点元或Nedelec棱边矢量元对象。
 *          支持全部23种单元类型的动态创建，封装了模板参数的显式指定逻辑。
 *
 * 支持的单元类型（共23种）：
 * - 1D Lagrange: LINE2, LINE3
 * - 2D Lagrange: TRI3, TRI6, QUAD4, QUAD8, QUAD9
 * - 3D Lagrange: TET4, TET10, HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13
 * - 2D Nedelec: TRI3_EDGE, QUAD4_EDGE
 * - 3D Nedelec: TET4_EDGE, HEX8_EDGE, PRISM6_EDGE, PYRAMID5_EDGE
 *
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include "shape_function_base.hpp"

namespace numeric {

/**
 * @class ShapeFunctionFactory
 * @brief 形函数单元静态工厂类
 * @details 纯静态工具类，禁止实例化。
 *          根据字符串标识符创建对应的ShapeFunctionBase派生类对象，
 *          将用户侧的字符串配置与内部模板实例化逻辑解耦。
 *
 * 使用示例：
 * @code
 *   // 通过字符串创建四边形单元
 *   auto elem = ShapeFunctionFactory::create("QUAD4");
 *   if (elem) {
 *       Eigen::VectorXd N = elem->evalN(LocalPoint(0.0, 0.0));
 *   }
 *
 *   // 检查类型是否支持
 *   if (ShapeFunctionFactory::isSupported("TET4_EDGE")) {
 *       auto nedelec = ShapeFunctionFactory::create("TET4_EDGE");
 *   }
 * @endcode
 *
 * @note 工厂方法返回std::unique_ptr<ShapeFunctionBase>，调用方持有所有权
 * @see ShapeFunctionBase 形函数抽象基类
 * @see ElementType 单元类型枚举
 */
class ShapeFunctionFactory {
public:
    /**
     * @brief 根据字符串标识符创建对应的形函数单元对象
     * @param element_type 单元类型的字符串标识符（如"TRI3"、"HEX8_EDGE"等）
     * @return std::unique_ptr<ShapeFunctionBase> 创建的单元对象指针，
     *         类型不支持或创建失败时返回nullptr
     *
     * @details 内部通过字符串匹配分派到对应的模板实例化：
     *          - Lagrange单元：LagrangeElement<Dim>(ElementType)
     *          - Nedelec单元：NedelecElement<Dim>(ElementType)
     *          其中Dim由单元类型隐式决定（1D/2D/3D）
     *
     * @warning 返回的智能指针可能为空，调用前需检查有效性
     */
    static std::unique_ptr<ShapeFunctionBase> create(const std::string& element_type);

    /**
     * @brief 检查指定的单元类型字符串是否被支持
     * @param element_type 单元类型字符串（如"TET10"、"PRISM6_EDGE"等）
     * @return bool 字符串匹配已知类型返回true，否则返回false
     *
     * @note 此方法不实际创建对象，仅做字符串存在性检查，开销极低
     */
    static bool isSupported(const std::string& element_type);

    /**
     * @brief 获取所有支持的单元类型名称列表
     * @return std::vector<std::string> 包含全部23种支持类型的字符串向量，
     *         顺序为：1D Lagrange → 2D Lagrange → 3D Lagrange → 2D Nedelec → 3D Nedelec
     */
    static std::vector<std::string> getSupportedTypes();

private:
    // 禁止实例化（纯静态工厂，无成员变量）
    ShapeFunctionFactory() = delete;
};

} // namespace numeric
