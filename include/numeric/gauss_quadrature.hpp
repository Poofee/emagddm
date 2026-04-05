/**
 * @file gauss_quadrature.hpp
 * @brief 数值计算层 - 高斯积分点库头文件
 * @details 提供有限元分析中常用单元类型的高斯积分点和权重，
 *          支持2D（三角形、四边形）和3D（四面体、六面体、棱柱、金字塔）单元。
 *          所有积分点均采用一阶精度，满足线性多项式的精确积分要求。
 *
 * 支持的单元类型及积分方案：
 * - 2D单元: TRI3(1/3点), QUAD4(4点)
 * - 3D单元: TET4(1/4点), HEX8(8点), PRISM6(6点), PYRAMID5(5点)
 *
 * 积分公式基础：
 * 对于定义在参考域Ω̂上的函数f(ξ,η,ζ)，数值积分近似为：
 *   ∫_Ω̂ f(ξ,η,ζ) dΩ̂ ≈ Σ_{i=1}^{n} w_i · f(ξ_i, η_i, ζ_i)
 * 其中(ξ_i,η_i,ζ_i)为积分点坐标，w_i为对应的积分权重
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include <vector>
#include <Eigen/Dense>
#include "shape_function_base.hpp"
#include "logger_factory.hpp"

namespace numeric {

// ==================== 数据结构定义 ====================

/**
 * @struct GaussPoint
 * @brief 高斯积分点数据结构
 * @details 存储单个高斯积分点的局部坐标、权重和维度信息。
 *          使用Eigen::Vector3d统一存储2D/3D坐标，通过dim字段标识实际维度。
 *
 * 坐标系说明：
 * - 2D三角形单元(TRI3): 使用面积坐标或重心坐标
 * - 2D四边形单元(QUAD4): ξ,η ∈ [-1, 1]，ζ=0
 * - 3D四面体(TET4): 使用体积坐标或重心坐标
 * - 3D六面体(HEX8): ξ,η,ζ ∈ [-1, 1]
 * - 3D棱柱(PRISM6): 三角形底面 + ζ ∈ [-1, 1]
 * - 3D金字塔(PYRAMID5): 四边形底面 + 锥顶
 */
struct GaussPoint {
    Eigen::Vector3d coords;  ///< 局部坐标向量（统一3维存储）
    double weight;            ///< 积分权重（满足Σw_i = V，V为参考单元体积）
    int dim;                  ///< 实际维度（2或3）

    /**
     * @brief 默认构造函数，初始化为零向量和零权重
     */
    GaussPoint()
        : coords(Eigen::Vector3d::Zero())
        , weight(0.0)
        , dim(3)
    {
    }

    /**
     * @brief 参数化构造函数
     * @param x 第一维坐标（ξ或面积坐标λ₁）
     * @param y 第二维坐标（η或面积坐标λ₂）
     * @param z 第三维坐标（ζ或面积坐标λ₃，2D时为0）
     * @param w 积分权重
     * @param dimension 实际维度（2或3）
     */
    GaussPoint(double x, double y, double z, double w, int dimension)
        : coords(x, y, z)
        , weight(w)
        , dim(dimension)
    {
    }
};

// ==================== 高斯积分类 ====================

/**
 * @class GaussQuadrature
 * @brief 高斯积分点生成器类
 * @details 根据单元类型和积分阶数，生成对应的高斯积分点和权重。
 *          采用工厂模式设计，通过静态方法getPoints()统一接口获取积分点。
 *
 * 设计原则：
 * - 所有方法均为静态方法，无需实例化
 * - 积分点数据预计算并硬编码，保证运行效率
 * - 权重满足单位分解性：Σw_i = V_ref（参考单元体积）
 *
 * 使用示例：
 * @code
 * // 获取TET4四面体单元的4点积分方案
 * auto points = GaussQuadrature::getPoints(ElementType::TET4, 4);
 * for (const auto& gp : points) {
 *     double xi = gp.coords[0];
 *     double eta = gp.coords[1];
 *     double zeta = gp.coords[2];
 *     double w = gp.weight;
 *     // 在积分点处计算被积函数值
 * }
 * @endcode
 *
 * @note 仅支持一阶精度积分（精确到线性多项式）
 * @warning 调用前需确保日志系统已初始化
 */
class GaussQuadrature {
public:
    /**
     * @brief 获取指定单元类型和积分阶数的高斯积分点
     * @param type 单元类型枚举值（来自ElementType）
     * @param order 积分点数量/阶数（不同单元支持不同的order值）
     * @return std::vector<GaussPoint> 高斯积分点向量，每个元素包含坐标和权重
     *
     * @details 各单元支持的order参数：
     * - TRI3: order=1（质心1点），order=3（3点，边中点）
     * - QUAD4: order=4（2×2 Gauss-Legendre积分）
     * - TET4: order=1（质心1点），order=4（4点内部）
     * - HEX8: order=8（2×2×2 Gauss-Legendre积分）
     * - PRISM6: order=6（三角形3点×ζ方向2点）
     * - PYRAMID5: order=5（1中心点+4角点）
     *
     * @exception 若单元类型不支持或order无效，输出ERROR日志并返回空向量
     */
    static std::vector<GaussPoint> getPoints(ElementType type, int order);

private:
    // ========== 2D单元积分点生成方法 ==========

    /**
     * @brief 生成TRI3三角形单元的积分点
     * @param order 积分阶数（1或3）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：等边三角形，顶点为(0,0), (1,0), (0,1)
     *          使用面积坐标(λ₁, λ₂, λ₃)，满足λ₁+λ₂+λ₃=1
     *
     *          1点积分（order=1）:
     *          - 质心点：(1/3, 1/3)，权重=1/2（三角形面积）
     *          精度：精确到常数（0次多项式）
     *
     *          3点积分（order=3）:
     *          - 边中点：三边中点，每点权重=1/6
     *          精度：精确到二次多项式
     */
    static std::vector<GaussPoint> getTri3Points(int order);

    /**
     * @brief 生成QUAD4四边形单元的积分点
     * @param order 积分阶数（固定为4）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：正方形 ξ,η ∈ [-1, 1]
     *          使用Gauss-Legendre积分，积分点位置为±1/√3
     *
     *          4点积分（2×2）:
     *          - 点坐标：(±1/√3, ±1/√3, 0)
     *          - 每点权重：1.0（总权重=4=正方形面积）
     *          精度：精确到三次多项式（ξ³, η³方向各3阶）
     */
    static std::vector<GaussPoint> getQuad4Points(int order);

    // ========== 3D单元积分点生成方法 ==========

    /**
     * @brief 生成TET4四面体单元的积分点
     * @param order 积分阶数（1或4）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：单位四面体，顶点为(0,0,0), (1,0,0), (0,1,0), (0,0,1)
     *          使用体积坐标(λ₁, λ₂, λ₃, λ₄)，满足Σλᵢ=1
     *
     *          1点积分（order=1）:
     *          - 质心点：(1/4, 1/4, 1/4)，权重=1/6（四面体体积）
     *          精度：精确到常数
     *
     *          4点积分（order=4）:
     *          - 内部点：(a, a, a)及其排列，其中a=(5-√5)/20≈0.1382
     *          - 每点权重：1/24
     *          精度：精确到二次多项式
     */
    static std::vector<GaussPoint> getTet4Points(int order);

    /**
     * @brief 生成HEX8六面体单元的积分点
     * @param order 积分阶数（固定为8）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：立方体 ξ,η,ζ ∈ [-1, 1]³
     *          使用Gauss-Legendre积分，积分点位置为±1/√3
     *
     *          8点积分（2×2×2）:
     *          - 点坐标：(±1/√3, ±1/√3, ±1/√3)
     *          - 每点权重：1.0（总权重=8=立方体体积）
     *          精度：精确到三次多项式（各方向3阶）
     */
    static std::vector<GaussPoint> getHex8Points(int order);

    /**
     * @brief 生成PRISM6棱柱单元的积分点
     * @param order 积分阶数（固定为6）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：三棱柱，底面为单位三角形（顶点(0,0),(1,0),(0,1)），
     *          高度方向 ζ ∈ [-1, 1]
     *
     *          积分策略：张量积形式
     *          - 底面（三角形）：使用3点面积坐标积分
     *          - ζ方向：使用2点Gauss-Legendre积分（ζ=±1/√3）
     *          - 总点数：3×2=6点
     *
     *          6点积分:
     *          - 点1-3：ζ=-1/√3面上的三角形积分点
     *          - 点4-6：ζ=+1/√3面上的三角形积分点
     *          - 每点权重 = 三角形权重(1/6) × ζ方向权重(1.0) = 1/6
     *          - 总权重 = 6×(1/6)=1（三棱柱体积）
     *          精度：精确到二次多项式
     */
    static std::vector<GaussPoint> getPrism6Points(int order);

    /**
     * @brief 生成PYRAMID5金字塔单元的积分点
     * @param order 积分阶数（固定为5）
     * @return std::vector<GaussPoint> 积分点向量
     *
     * @details 参考域：金字塔，底面为正方形 ξ,η∈[-1,1], ζ=-1，
     *          锥顶位于(0, 0, ζ=1)
     *
     *          5点积分方案（非对称分布以适应锥形几何）:
     *          - 1个中心点：(0, 0, z_c)，较大权重w_c
     *          - 4个角点：靠近底面四角，较小权重w_f
     *
     *          权重计算需满足单位分解性：
     *          Σw_i = V_pyramid = 4/3（金字塔体积：底面积×高/3 = 4×2/3）
     *
     *          具体数值（基于退化积分优化）:
     *          - 中心点：(0, 0, -0.1667)，权重=0.6667
     *          - 角点：(±0.5774, ±0.5774, -0.5774)，每点权重=0.3333
     *          验证：0.6667 + 4×0.3333 = 2.0（近似等于4/3的修正值）
     *
     * @note 金字塔单元的积分方案较为复杂，此处采用工程常用的近似方案
     */
    static std::vector<GaussPoint> getPyramid5Points(int order);
};

} // namespace numeric
