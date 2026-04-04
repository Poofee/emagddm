/**
 * @file lagrange_element.hpp
 * @brief 数值计算层 - Lagrange标量节点有限元单元模板类声明
 * @details 实现全部19种Lagrange单元类型的形函数与梯度计算，
 *          包括1D线单元、2D面单元、3D体单元。
 *          模板参数Dim指定空间维度（1/2/3），
 *          通过ElementType枚举分派到具体单元类型的计算逻辑。
 *
 * 支持的单元类型：
 * - 1D: LINE2(2节点), LINE3(3节点二次)
 * - 2D: TRI3, TRI6, QUAD4, QUAD8, QUAD9
 * - 3D: TET4, TET10, HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13
 *
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include "shape_function_base.hpp"

namespace numeric {

/**
 * @class LagrangeElement
 * @brief Lagrange标量节点有限元单元模板类
 * @tparam Dim 空间维度（1=一维, 2=二维, 3=三维）
 *
 * @details 继承ShapeFunctionBase抽象基类，实现Lagrange插值的标量形函数。
 *          所有形函数满足单位分解性 ΣN_i = 1，
 *          偏导数通过解析求导精确计算。
 *
 * 使用示例：
 * @code
 *   // 创建二维四边形单元
 *   auto quad = std::make_unique<LagrangeElement<2>>(ElementType::QUAD4);
 *
 *   // 计算形函数值
 *   LocalPoint xi(0.0, 0.0);  // 单元中心
 *   Eigen::VectorXd N = quad->evalN(xi);
 *
 *   // 计算形函数梯度
 *   Eigen::MatrixXd gradN = quad->evalGradN(xi);
 * @endcode
 */
template<int Dim>
class LagrangeElement : public ShapeFunctionBase {
public:
    /**
     * @brief 构造函数，初始化单元类型
     * @param element_type 有限元单元类型枚举值
     * @throws 无显式异常，非法类型通过日志警告并设置默认值
     */
    explicit LagrangeElement(ElementType element_type)
        : element_type_(element_type) {
        FEEM_INFO("初始化Lagrange单元: Dim={}, Type={}", Dim,
                  static_cast<int>(element_type_));
    }

    /**
     * @brief 析构函数（默认）
     */
    ~LagrangeElement() override = default;

    // ========== 纯虚方法实现 ==========

    /**
     * @brief 获取单元节点类型
     * @return ElementType 当前单元的类型枚举值
     */
    ElementType getNodeType() const override {
        return element_type_;
    }

    /**
     * @brief 获取单元节点总数
     * @return int 节点数量（如LINE2返回2，HEX8返回8）
     */
    int getNodeCount() const override {
        switch (element_type_) {
            case ElementType::LINE2:   return 2;
            case ElementType::LINE3:   return 3;
            case ElementType::TRI3:    return 3;
            case ElementType::TRI6:    return 6;
            case ElementType::QUAD4:   return 4;
            case ElementType::QUAD8:   return 8;
            case ElementType::QUAD9:   return 9;
            case ElementType::TET4:    return 4;
            case ElementType::TET10:   return 10;
            case ElementType::HEX8:    return 8;
            case ElementType::HEX20:   return 20;
            case ElementType::HEX27:   return 27;
            case ElementType::PRISM6:  return 6;
            case ElementType::PRISM15: return 15;
            case ElementType::PYRAMID5:  return 5;
            case ElementType::PYRAMID13: return 13;
            default:
                FEEM_WARN("未知的单元类型: {}, 返回0", static_cast<int>(element_type_));
                return 0;
        }
    }

    /**
     * @brief 获取单元空间维度
     * @return int 空间维度（由模板参数Dim决定）
     */
    int getDim() const override { return Dim; }

    /**
     * @brief 计算参考域形函数值向量
     * @param xi 参考域局部坐标点
     * @return Eigen::VectorXd 形函数值向量（node_count × 1）
     *
     * @note 内部根据element_type_分派到对应单元类型的计算函数
     * @warning 调用前需确保xi的dim与单元维度匹配
     */
    Eigen::VectorXd evalN(const LocalPoint& xi) const override {
        double x = xi.coords[0];
        double y = xi.coords[1];
        double z = xi.coords[2];
        if constexpr (Dim == 1) {
            switch (element_type_) {
                case ElementType::LINE2: return evalLine2N(x);
                case ElementType::LINE3: return evalLine3N(x);
                default:
                    FEEM_ERROR("evalN: 1D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::VectorXd();
            }
        } else if constexpr (Dim == 2) {
            switch (element_type_) {
                case ElementType::TRI3:  return evalTri3N(x, y);
                case ElementType::TRI6:  return evalTri6N(x, y);
                case ElementType::QUAD4: return evalQuad4N(x, y);
                case ElementType::QUAD8: return evalQuad8N(x, y);
                case ElementType::QUAD9: return evalQuad9N(x, y);
                default:
                    FEEM_ERROR("evalN: 2D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::VectorXd();
            }
        } else {
            switch (element_type_) {
                case ElementType::TET4:      return evalTet4N(x, y, z);
                case ElementType::TET10:     return evalTet10N(x, y, z);
                case ElementType::HEX8:      return evalHex8N(x, y, z);
                case ElementType::HEX20:     return evalHex20N(x, y, z);
                case ElementType::HEX27:     return evalHex27N(x, y, z);
                case ElementType::PRISM6:    return evalPrism6N(x, y, z);
                case ElementType::PRISM15:   return evalPrism15N(x, y, z);
                case ElementType::PYRAMID5:  return evalPyramid5N(x, y, z);
                case ElementType::PYRAMID13: return evalPyramid13N(x, y, z);
                default:
                    FEEM_ERROR("evalN: 3D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::VectorXd();
            }
        }
    }

    /**
     * @brief 计算参考域形函数梯度矩阵
     * @param xi 参考域局部坐标点
     * @return Eigen::MatrixXd 形函数梯度矩阵（node_count × dim）
     *
     * @note 对于三维单元，每行为[∂N_i/∂ξ, ∂N_i/∂η, ∂N_i/∂ζ]
     */
    Eigen::MatrixXd evalGradN(const LocalPoint& xi) const override {
        double x = xi.coords[0];
        double y = xi.coords[1];
        double z = xi.coords[2];
        if constexpr (Dim == 1) {
            switch (element_type_) {
                case ElementType::LINE2: return evalLine2GradN(x);
                case ElementType::LINE3: return evalLine3GradN(x);
                default:
                    FEEM_ERROR("evalGradN: 1D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::MatrixXd();
            }
        } else if constexpr (Dim == 2) {
            switch (element_type_) {
                case ElementType::TRI3:  return evalTri3GradN(x, y);
                case ElementType::TRI6:  return evalTri6GradN(x, y);
                case ElementType::QUAD4: return evalQuad4GradN(x, y);
                case ElementType::QUAD8: return evalQuad8GradN(x, y);
                case ElementType::QUAD9: return evalQuad9GradN(x, y);
                default:
                    FEEM_ERROR("evalGradN: 2D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::MatrixXd();
            }
        } else {
            switch (element_type_) {
                case ElementType::TET4:      return evalTet4GradN(x, y, z);
                case ElementType::TET10:     return evalTet10GradN(x, y, z);
                case ElementType::HEX8:      return evalHex8GradN(x, y, z);
                case ElementType::HEX20:     return evalHex20GradN(x, y, z);
                case ElementType::HEX27:     return evalHex27GradN(x, y, z);
                case ElementType::PRISM6:    return evalPrism6GradN(x, y, z);
                case ElementType::PRISM15:   return evalPrism15GradN(x, y, z);
                case ElementType::PYRAMID5:  return evalPyramid5GradN(x, y, z);
                case ElementType::PYRAMID13: return evalPyramid13GradN(x, y, z);
                default:
                    FEEM_ERROR("evalGradN: 3D单元不支持类型 {}", static_cast<int>(element_type_));
                    return Eigen::MatrixXd();
            }
        }
    }

private:
    ElementType element_type_;  ///< 单元类型枚举值

    // ========== 1D单元计算函数 ==========

    /**
     * @brief LINE2二节点线性单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 [N₁, N₂]
     *
     * @details N₁ = (1-ξ)/2, N₂ = (1+ξ)/2
     */
    static Eigen::VectorXd evalLine2N(double xi);

    /**
     * @brief LINE2单元形函数梯度计算
     * @param xi 局部坐标（未使用，梯度为常数）
     * @return Eigen::MatrixXd 梯度矩阵 [dN₁/dξ; dN₂/dξ] (2×1)
     */
    static Eigen::MatrixXd evalLine2GradN(double xi);

    /**
     * @brief LINE3三节点二次单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 [N₁, N₂, N₃]
     *
     * @details N₁ = ξ(ξ-1)/2, N₂ = 1-ξ², N₃ = ξ(ξ+1)/2
     */
    static Eigen::VectorXd evalLine3N(double xi);

    /**
     * @brief LINE3单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @return Eigen::MatrixXd 梯度矩阵 (3×1)
     */
    static Eigen::MatrixXd evalLine3GradN(double xi);

    // ========== 2D单元计算函数 ==========

    /**
     * @brief TRI3三节点三角形单元形函数计算
     * @param xi 第一面积坐标 ξ ≥ 0
     * @param eta 第二面积坐标 η ≥ 0
     * @return Eigen::VectorXd 形函数值 [N₁, N₂, N₃]
     *
     * @details 使用面积坐标：N₁=ξ, N₂=η, N₃=1-ξ-η
     *          约束条件：ξ≥0, η≥0, ξ+η≤1
     */
    static Eigen::VectorXd evalTri3N(double xi, double eta);

    /**
     * @brief TRI3单元形函数梯度计算
     * @param xi 面积坐标 ξ
     * @param eta 面积坐标 η
     * @return Eigen::MatrixXd 梯度矩阵 (3×2)，每行为[∂N/∂ξ, ∂N/∂η]
     */
    static Eigen::MatrixXd evalTri3GradN(double xi, double eta);

    /**
     * @brief TRI6六节点二次三角形单元形函数计算
     * @param xi 面积坐标 ξ
     * @param eta 面积坐标 η
     * @return Eigen::VectorXd 形函数值 (6×1)
     *
     * @details 角点：Nᵢ=Lᵢ(2Lᵢ-1), 边中：N₄=4ξη, N₅=4ηζ, N₆=4ζξ
     *          其中 ζ=1-ξ-η
     */
    static Eigen::VectorXd evalTri6N(double xi, double eta);

    /**
     * @brief TRI6单元形函数梯度计算
     * @param xi 面积坐标 ξ
     * @param eta 面积坐标 η
     * @return Eigen::MatrixXd 梯度矩阵 (6×2)
     */
    static Eigen::MatrixXd evalTri6GradN(double xi, double eta);

    /**
     * @brief QUAD4四节点双线性四边形单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (4×1)
     *
     * @details Nᵢ = ¼(1+ξᵢξ)(1+ηᵢη)，i=1..4
     *          节点顺序：(-1,-1)=1, (1,-1)=2, (1,1)=3, (-1,1)=4
     */
    static Eigen::VectorXd evalQuad4N(double xi, double eta);

    /**
     * @brief QUAD4单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @return Eigen::MatrixXd 梯度矩阵 (4×2)
     */
    static Eigen::MatrixXd evalQuad4GradN(double xi, double eta);

    /**
     * @brief QUAD8八节点Serendipity四边形单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (8×1)
     *
     * @details 角点和边中节点的Serendipity形函数
     */
    static Eigen::VectorXd evalQuad8N(double xi, double eta);

    /**
     * @brief QUAD8单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @return Eigen::MatrixXd 梯度矩阵 (8×2)
     */
    static Eigen::MatrixXd evalQuad8GradN(double xi, double eta);

    /**
     * @brief QUAD9九节点双二次四边形单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (9×1)
     *
     * @details 使用张量积形式：Nᵢⱼ = Lᵢ(ξ)Lⱼ(η)
     *          其中L为1D二次Lagrange多项式
     */
    static Eigen::VectorXd evalQuad9N(double xi, double eta);

    /**
     * @brief QUAD9单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @return Eigen::MatrixXd 梯度矩阵 (9×2)
     */
    static Eigen::MatrixXd evalQuad9GradN(double xi, double eta);

    // ========== 3D单元计算函数 ==========

    /**
     * @brief TET4四节点四面体单元形函数计算
     * @param xi 第一体积坐标 ξ
     * @param eta 第二体积坐标 η
     * @param zeta 第三体积坐标 ζ
     * @return Eigen::VectorXd 形函数值 (4×1)
     *
     * @details 使用体积坐标：N₁=ξ, N₂=η, N₃=ζ, N₄=1-ξ-η-ζ
     */
    static Eigen::VectorXd evalTet4N(double xi, double eta, double zeta);

    /**
     * @brief TET4单元形函数梯度计算
     * @param xi 体积坐标 ξ
     * @param eta 体积坐标 η
     * @param zeta 体积坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (4×3)
     */
    static Eigen::MatrixXd evalTet4GradN(double xi, double eta, double zeta);

    /**
     * @brief TET10十节点二次四面体单元形函数计算
     * @param xi 体积坐标 ξ
     * @param eta 体积坐标 η
     * @param zeta 体积坐标 ζ
     * @return Eigen::VectorXd 形函数值 (10×1)
     *
     * @details 角点：Nᵢ=Lᵢ(2Lᵢ-1)，边中：Nⱼ=4LₖLₗ
     */
    static Eigen::VectorXd evalTet10N(double xi, double eta, double zeta);

    /**
     * @brief TET10单元形函数梯度计算
     * @param xi 体积坐标 ξ
     * @param eta 体积坐标 η
     * @param zeta 体积坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (10×3)
     */
    static Eigen::MatrixXd evalTet10GradN(double xi, double eta, double zeta);

    /**
     * @brief HEX8八节点三线性六面体单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @param zeta 局部坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (8×1)
     *
     * @details Nᵢ = ⅛(1+ξᵢξ)(1+ηᵢη)(1+ζᵢζ)，i=1..8
     */
    static Eigen::VectorXd evalHex8N(double xi, double eta, double zeta);

    /**
     * @brief HEX8单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @param zeta 局部坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (8×3)
     */
    static Eigen::MatrixXd evalHex8GradN(double xi, double eta, double zeta);

    /**
     * @brief HEX20二十节点Serendipity六面体单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @param zeta 局部坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (20×1)
     */
    static Eigen::VectorXd evalHex20N(double xi, double eta, double zeta);

    /**
     * @brief HEX20单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @param zeta 局部坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (20×3)
     */
    static Eigen::MatrixXd evalHex20GradN(double xi, double eta, double zeta);

    /**
     * @brief HEX27二十七节点三二次六面体单元形函数计算
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @param zeta 局部坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (27×1)
     *
     * @details 张量积形式：Nᵢⱼₖ = Lᵢ(ξ)Lⱼ(η)Lₖ(ζ)
     */
    static Eigen::VectorXd evalHex27N(double xi, double eta, double zeta);

    /**
     * @brief HEX27单元形函数梯度计算
     * @param xi 局部坐标 ξ
     * @param eta 局部坐标 η
     * @param zeta 局部坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (27×3)
     */
    static Eigen::MatrixXd evalHex27GradN(double xi, double eta, double zeta);

    /**
     * @brief PRISM6六节点线性三棱柱单元形函数计算
     * @param xi 三角形面坐标 ξ ≥ 0
     * @param eta 三角形面坐标 η ≥ 0
     * @param zeta 棱柱高度坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (6×1)
     *
     * @details 底面：Nᵢ = Lᵢ·(1-ζ)/2，顶面：Nᵢ₊₃ = Lᵢ·(1+ζ)/2
     */
    static Eigen::VectorXd evalPrism6N(double xi, double eta, double zeta);

    /**
     * @brief PRISM6单元形函数梯度计算
     * @param xi 面坐标 ξ
     * @param eta 面坐标 η
     * @param zeta 高度坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (6×3)
     */
    static Eigen::MatrixXd evalPrism6GradN(double xi, double eta, double zeta);

    /**
     * @brief PRISM15十五节点二次三棱柱单元形函数计算
     * @param xi 面坐标 ξ
     * @param eta 面坐标 η
     * @param zeta 高度坐标 ζ
     * @return Eigen::VectorXd 形函数值 (15×1)
     */
    static Eigen::VectorXd evalPrism15N(double xi, double eta, double zeta);

    /**
     * @brief PRISM15单元形函数梯度计算
     * @param xi 面坐标 ξ
     * @param eta 面坐标 η
     * @param zeta 高度坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (15×3)
     */
    static Eigen::MatrixXd evalPrism15GradN(double xi, double eta, double zeta);

    /**
     * @brief PYRAMID5五节点线性金字塔单元形函数计算
     * @param xi 底面坐标 ξ ∈ [-1, 1]
     * @param eta 底面坐标 η ∈ [-1, 1]
     * @param zeta 高度坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (5×1)
     *
     * @details 使用有理多项式形式避免奇异性
     */
    static Eigen::VectorXd evalPyramid5N(double xi, double eta, double zeta);

    /**
     * @brief PYRAMID5单元形函数梯度计算
     * @param xi 底面坐标 ξ
     * @param eta 底面坐标 η
     * @param zeta 高度坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (5×3)
     */
    static Eigen::MatrixXd evalPyramid5GradN(double xi, double eta, double zeta);

    /**
     * @brief PYRAMID13十三节点二次金字塔单元形函数计算
     * @param xi 底面坐标 ξ ∈ [-1, 1]
     * @param eta 底面坐标 η ∈ [-1, 1]
     * @param zeta 高度坐标 ζ ∈ [-1, 1]
     * @return Eigen::VectorXd 形函数值 (13×1)
     */
    static Eigen::VectorXd evalPyramid13N(double xi, double eta, double zeta);

    /**
     * @brief PYRAMID13单元形函数梯度计算
     * @param xi 底面坐标 ξ
     * @param eta 底面坐标 η
     * @param zeta 高度坐标 ζ
     * @return Eigen::MatrixXd 梯度矩阵 (13×3)
     */
    static Eigen::MatrixXd evalPyramid13GradN(double xi, double eta, double zeta);
};

} // namespace numeric
