/**
 * @file nedelec_element.hpp
 * @brief 数值计算层 - Nedelec第一类一阶棱边矢量有限元单元
 * @details 实现Nédélec I型H(curl)协调单元，支持2D/3D共6种单元类型：
 *          - 2D: TRI3_EDGE（3棱边）、QUAD4_EDGE（4棱边）
 *          - 3D: TET4_EDGE（6棱边）、HEX8_EDGE（12棱边）、PRISM6_EDGE（9棱边）、PYRAMID5_EDGE（8棱边）
 *
 *          核心数学基础：
 *          基函数采用Whitney形式（梯度-坐标乘积差），沿棱边具有单位环流性质：
 *          ∫_{edge_k} W_e · t dl = δ_ek
 *
 *          物理域映射使用Piola变换：
 *          W_phys(x) = J^{-T} · W_ref(ξ) / det(J)
 *          curl_phys(x) = curl_ref(ξ) / det(J)
 *
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include "shape_function_base.hpp"
#include <vector>
#include <utility>

namespace numeric {

/**
 * @class NedelecElement
 * @brief Nédélec第一类一阶棱边矢量有限元单元模板类
 * @tparam Dim 空间维度模板参数（2或3）
 *
 * @details 继承ShapeFunctionBase基类，实现所有纯虚接口和Nedelec扩展接口。
 *          内部存储棱边连接信息（每条棱边的起止节点对），
 *          根据构造时传入的ElementType自动选择对应的基函数公式。
 *
 * @note 对于Nedelec单元，getNodeCount()返回的是**棱边数量**（自由度数），而非节点数
 * @see ShapeFunctionBase 形函数抽象基类
 * @see ElementType 单元类型枚举
 */
template<int Dim>
class NedelecElement : public ShapeFunctionBase {
public:
    /**
     * @brief 构造函数，根据单元类型初始化棱边连接信息
     * @param type 单元类型枚举值，必须是_EDGE后缀的Nedelec类型
     *
     * @throws 无显式抛出，非法类型通过日志警告并设置默认值
     */
    explicit NedelecElement(ElementType type);

    /**
     * @brief 虚析构函数
     */
    ~NedelecElement() override = default;

    // ========== 纯虚方法重写 ==========

    /**
     * @brief 获取单元节点类型
     * @return ElementType 当前Nedelec单元的类型枚举值
     */
    ElementType getNodeType() const override;

    /**
     * @brief 获取单元棱边总数（即自由度数量）
     * @return int 棱边数量：TRI3_EDGE=3, QUAD4_EDGE=4, TET4_EDGE=6,
     *             HEX8_EDGE=12, PRISM6_EDGE=9, PYRAMID5_EDGE=8
     */
    int getNodeCount() const override;

    /**
     * @brief 获取单元空间维度
     * @return int 模板参数Dim的值（2或3）
     */
    int getDim() const override;

    /**
     * @brief 计算参考域形函数值向量（Nedelec元返回零向量）
     * @param xi 参考域局部坐标点
     * @return Eigen::VectorXd 长度为edge_count的零向量
     *
     * @note Nedelec单元的主要功能在evalEdgeFunction中实现，
     *       此方法保留接口兼容性，返回零向量
     */
    Eigen::VectorXd evalN(const LocalPoint& xi) const override;

    /**
     * @brief 计算参考域形函数梯度矩阵（Nedelec元返回零矩阵）
     * @param xi 参考域局部坐标点
     * @return Eigen::MatrixXd 维度为edge_count × Dim的零矩阵
     *
     * @note 同evalN()，Nedelec元的梯度信息通过evalCurlEdge获取
     */
    Eigen::MatrixXd evalGradN(const LocalPoint& xi) const override;

    // ========== Nedelec扩展接口重写 ==========

    /**
     * @brief 计算Nedelec棱边形函数在参考域的矢量值
     * @param edge_index 棱边索引（从0开始），必须在[0, edge_count)范围内
     * @param xi 参考域局部坐标点
     * @return Eigen::Vector3d 第edge_index号棱边形函数在xi处的参考域矢量值
     *
     * @details 根据单元类型分派到具体的基函数计算方法：
     *          - TRI3_EDGE: 面积坐标系下的Whitney形式
     *          - QUAD4_EDGE: 双线性坐标系下的张量积形式
     *          - TET4_EDGE: 体坐标系下的Whitney形式
     *          - HEX8_EDGE: 三线性坐标系下的张量积形式
     *          - PRISM6_EDGE: 三角形-棱柱混合形式
     *          - PYRAMID5_EDGE: 金字塔退化形式
     *
     * @warning edge_index越界时输出ERROR日志并返回零向量
     */
    Eigen::Vector3d evalEdgeFunction(
        int edge_index,
        const LocalPoint& xi
    ) const override;

    /**
     * @brief 计算Nedelec棱边形函数在参考域的旋度
     * @param edge_index 棱边索引（从0开始）
     * @param xi 参考域局部坐标点
     * @return Eigen::Vector3d 旋度向量（2D时仅前两个分量有效，第3分量=0）
     *
     * @note 对于一阶Nedelec元，旋度在参考域内为常矢量（或简单表达式），
     *       这是其重要性质，可显著简化curl-curl矩阵的组装
     */
    Eigen::Vector3d evalCurlEdge(
        int edge_index,
        const LocalPoint& xi
    ) const override;

    /**
     * @brief 计算Nedelec棱边形函数在物理域的完整结果（含Piola变换映射）
     * @param edge_index 棱边索引（从0开始）
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵，维度为 Dim × node_count
     * @return EdgeCurlResult 包含参考域旋度和物理域旋度（经detJ缩放）
     *
     * @details 物理域旋度映射公式：
     *          curl_physical = (1/detJ) · J · curl_reference
     *          其中J为雅可比矩阵，detJ为其行列式
     *
     * @note 此方法覆写基类默认实现，确保Piola变换正确应用于旋度
     */
    EdgeCurlResult calcPhysicalEdgeFunction(
        int edge_index,
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords
    ) const override;

private:
    ElementType element_type_;                          ///< 单元类型枚举值
    int edge_count_;                                    ///< 棱边总数（自由度数）
    std::vector<std::pair<int, int>> edge_nodes_;       ///< 棱边-节点连接表：(起始节点, 终止节点)

    // ========== 各单元类型的基函数计算方法 ==========

    /**
     * @brief TRI3_EDGE三角形棱边元基函数计算
     * @param edge_index 棱边索引（0-2）
     * @param xi 面积坐标 ξ = λ₁
     * @param eta 面积坐标 η = λ₂
     * @return Eigen::Vector3d 参考域矢量值（z分量=0）
     *
     * @details 使用面积坐标(ξ,η,ζ)其中ζ=1-ξ-η，基函数为Whitney形式：
     *          W₁ = λ₂∇λ₁ - λ₁∇λ₂ = (η, -ξ)
     *          W₂ = λ₃∇λ₂ - λ₂∇λ₃ = (η, 1-ξ)
     *          W₃ = λ₁∇λ₃ - λ₃∇λ₁ = (η-1, -ξ)
     */
    Eigen::Vector3d evalTri3Edge(int edge_index, double xi, double eta) const;

    /**
     * @brief QUAD4_EDGE四边形棱边元基函数计算
     * @param edge_index 棱边索引（0-3）
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @return Eigen::Vector3d 参考域矢量值（z分量=0）
     *
     * @details 基于双线性坐标系，使用"气泡函数×梯度"形式：
     *          ξ方向棱边: W ∝ (1∓η)∇ξ
     *          η方向棱边: W ∝ (1±ξ)∇η
     */
    Eigen::Vector3d evalQuad4Edge(int edge_index, double xi, double eta) const;

    /**
     * @brief TET4_EDGE四面体棱边元基函数计算
     * @param edge_index 棱边索引（0-5）
     * @param xi 体坐标 λ₁
     * @param eta 体坐标 λ₂
     * @param zeta 体坐标 λ₃
     * @return Eigen::Vector3d 参考域矢量值
     *
     * @details 使用体坐标系(λ₁,λ₂,λ₃,λ₄)，λ₄=1-ξ-η-ζ，
     *          6条棱边的Whitney形式基函数：W_ij = λ_j∇λ_i - λ_i∇λ_j
     */
    Eigen::Vector3d evalTet4Edge(int edge_index, double xi, double eta, double zeta) const;

    /**
     * @brief HEX8_EDGE六面体棱边元基函数计算
     * @param edge_index 棱边索引（0-11）
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @param zeta 局部坐标 ζ ∈ [-1, 1]
     * @return Eigen::Vector3d 参考域矢量值
     *
     * @details 基于三线性坐标系(Gmsh标准编号)，12条棱边分为三组：
     *          - ξ方向棱边(4条): W ∝ (1∓η)(1∓ζ)∇ξ
     *          - η方向棱边(4条): W ∝ (1±ξ)(1∓ζ)∇η
     *          - ζ方向棱边(4条): W ∝ (1±ξ)(1±η)∇ζ
     */
    Eigen::Vector3d evalHex8Edge(int edge_index, double xi, double eta, double zeta) const;

    /**
     * @brief PRISM6_EDGE三棱柱棱边元基函数计算
     * @param edge_index 棱边索引（0-8）
     * @param xi 三角形面积坐标 ξ ≥ 0
     * @param eta 三角形面积坐标 η ≥ 0, ξ+η ≤ 1
     * @param zeta 棱柱轴向坐标 ζ ∈ [-1, 1]
     * @return Eigen::Vector3d 参考域矢量值
     *
     * @details 9条棱边分三组：
     *          - 底面3条(ζ=-1): TRI3基函数 × ½(1-ζ)
     *          - 顶面3条(ζ=1): TRI3基函数 × ½(1+ζ)
     *          - 侧棱3条: 面积坐标函数 × ∇ζ
     */
    Eigen::Vector3d evalPrism6Edge(int edge_index, double xi, double eta, double zeta) const;

    /**
     * @brief PYRAMID5_EDGE金字塔棱边元基函数计算
     * @param edge_index 棱边索引（0-7）
     * @param xi 局部坐标 ξ ∈ [-1, 1]
     * @param eta 局部坐标 η ∈ [-1, 1]
     * @param zeta 高度坐标 ζ ∈ [0, 1]（锥顶在ζ=1）
     * @return Eigen::Vector3d 参考域矢量值
     *
     * @details 8条棱边分两组：
     *          - 底面4条(ζ=0): QUAD4基函数 × (1-ζ)
     *          - 侧棱4条(底面角点到锥顶): 退化的屋顶函数组合
     */
    Eigen::Vector3d evalPyramid5Edge(int edge_index, double xi, double eta, double zeta) const;

    // ========== 各单元类型的旋度计算方法 ==========

    /**
     * @brief TRI3_EDGE旋度计算（常标量，扩展为3D向量）
     * @param edge_index 棱边索引
     * @return Eigen::Vector3d 旋度向量，形式为(0, 0, curl_scalar)
     */
    Eigen::Vector3d curlTri3Edge(int edge_index) const;

    /**
     * @brief QUAD4_EDGE旋度计算（常标量，扩展为3D向量）
     * @param edge_index 棱边索引
     * @return Eigen::Vector3d 旋度向量
     */
    Eigen::Vector3d curlQuad4Edge(int edge_index) const;

    /**
     * @brief TET4_EDGE旋度计算（常矢量）
     * @param edge_index 棱边索引
     * @return Eigen::Vector3d 常旋度矢量 curl = 2∇λ_j × ∇λ_i
     */
    Eigen::Vector3d curlTet4Edge(int edge_index) const;

    /**
     * @brief HEX8_EDGE旋度计算（依赖坐标的矢量）
     * @param edge_index 棱边索引
     * @param xi, eta, zeta 局部坐标
     * @return Eigen::Vector3d 旋度矢量
     */
    Eigen::Vector3d curlHex8Edge(int edge_index, double xi, double eta, double zeta) const;

    /**
     * @brief PRISM6_EDGE旋度计算
     * @param edge_index 棱边索引
     * @param xi, eta, zeta 局部坐标
     * @return Eigen::Vector3d 旋度矢量
     */
    Eigen::Vector3d curlPrism6Edge(int edge_index, double xi, double eta, double zeta) const;

    /**
     * @brief PYRAMID5_EDGE旋度计算
     * @param edge_index 棱边索引
     * @param xi, eta, zeta 局部坐标
     * @return Eigen::Vector3d 旋度矢量
     */
    Eigen::Vector3d curlPyramid5Edge(int edge_index, double xi, double eta, double zeta) const;

    // ========== 初始化辅助方法 ==========

    /**
     * @brief 根据单元类型初始化棱边连接信息和计数
     * @param type 单元类型枚举值
     *
     * @details 设置edge_nodes_和edge_count_，
     *          每种单元类型有固定的拓扑结构
     */
    void initEdgeConnectivity(ElementType type);
};

} // namespace numeric
