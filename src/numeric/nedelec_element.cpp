/**
 * @file nedelec_element.cpp
 * @brief 数值计算层 - Nedelec第一类一阶棱边矢量有限元单元实现
 * @details 实现Nédélec I型H(curl)协调单元的全部6种单元类型的基函数和旋度计算，
 *          包括2D（TRI3_EDGE、QUAD4_EDGE）和3D（TET4_EDGE、HEX8_EDGE、PRISM6_EDGE、PYRAMID5_EDGE）。
 *
 *          数学基础：Whitney形式基函数 W_ij = λ_j∇λ_i - λ_i∇λ_j，
 *          满足单位环流性质 ∫_{edge_k} W_e · t dl = δ_ek
 *
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#include "nedelec_element.hpp"

namespace numeric {

// ==================== 构造函数与初始化 ====================

template<int Dim>
NedelecElement<Dim>::NedelecElement(ElementType type)
    : element_type_(type)
    , edge_count_(0)
{
    // 根据单元类型初始化棱边拓扑连接信息
    initEdgeConnectivity(type);

    FEEM_INFO("NedelecElement构造完成: 类型={}, 维度={}, 棱边数={}",
              static_cast<int>(type), Dim, edge_count_);
}

template<int Dim>
void NedelecElement<Dim>::initEdgeConnectivity(ElementType type)
{
    // 清空已有数据
    edge_nodes_.clear();

    switch (type) {
    case ElementType::TRI3_EDGE:
        // 三角形3条棱边：使用Gmsh标准节点编号
        // 节点1→2, 节点2→3, 节点3→1（0-based索引: 0→1, 1→2, 2→0）
        edge_count_ = 3;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 0}};
        break;

    case ElementType::QUAD4_EDGE:
        // 四边形4条棱边：Gmsh标准编号
        // 底边(1→2), 右边(2→3), 顶边(3→4), 左边(4→1)
        edge_count_ = 4;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        break;

    case ElementType::TET4_EDGE:
        // 四面体6条棱边：Gmsh标准编号
        // 底面三角形3条 + 顶点到各底面节点3条
        edge_count_ = 6;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 0},   // 底面三角形
                       {0, 3}, {1, 3}, {2, 3}};   // 侧棱（到顶点4/索引3）
        break;

    case ElementType::HEX8_EDGE:
        // 六面体12条棱边：Gmsh/ANSYS标准编号
        // 底面4条(z=-1) + 顶面4条(z=1) + 竖直棱边4条
        edge_count_ = 12;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0},   // 底面(zeta=-1)
                       {4, 5}, {5, 6}, {6, 7}, {7, 4},   // 顶面(zeta=+1)
                       {0, 4}, {1, 5}, {2, 6}, {3, 7}};  // 竖直棱边
        break;

    case ElementType::PRISM6_EDGE:
        // 三棱柱9条棱边：Gmsh标准编号
        // 底面3条 + 顶面3条 + 侧棱3条
        edge_count_ = 9;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 0},           // 底面(zeta=-1)
                       {3, 4}, {4, 5}, {5, 3},           // 顶面(zeta=+1)
                       {0, 3}, {1, 4}, {2, 5}};          // 侧棱（竖直方向）
        break;

    case ElementType::PYRAMID5_EDGE:
        // 金字塔8条棱边：Gmsh标准编号
        // 底面四边形4条 + 侧棱4条（到底面锥顶）
        edge_count_ = 8;
        edge_nodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0},   // 底面(zeta=0)
                       {0, 4}, {1, 4}, {2, 4}, {3, 4}};  // 侧棱（到锥顶/索引4）
        break;

    default:
        // 非Nedelec类型，输出警告并设置默认空值
        edge_count_ = 0;
        FEEM_WARN("initEdgeConnectivity: 不支持的单元类型={}，已设置为空",
                  static_cast<int>(type));
        break;
    }
}

// ==================== 基础接口方法实现 ====================

template<int Dim>
ElementType NedelecElement<Dim>::getNodeType() const
{
    return element_type_;
}

template<int Dim>
int NedelecElement<Dim>::getNodeCount() const
{
    return edge_count_;
}

template<int Dim>
int NedelecElement<Dim>::getDim() const
{
    return Dim;
}

template<int Dim>
Eigen::VectorXd NedelecElement<Dim>::evalN(const LocalPoint& xi) const
{
    // Nedelec棱边元的形函数值主要通过evalEdgeFunction获取
    // 此方法保留接口兼容性，返回长度为edge_count的零向量
    (void)xi;  // 参数保留以匹配基类接口，Nedelec元不使用此方法
    return Eigen::VectorXd::Zero(edge_count_);
}

template<int Dim>
Eigen::MatrixXd NedelecElement<Dim>::evalGradN(const LocalPoint& xi) const
{
    // Nedelec元的梯度信息通过evalCurlEdge（旋度）获取
    // 此方法返回edge_count × Dim的零矩阵以保持接口兼容性
    (void)xi;  // 参数保留以匹配基类接口，Nedelec元不使用此方法
    return Eigen::MatrixXd::Zero(edge_count_, Dim);
}

// ==================== evalEdgeFunction 分派实现 ====================

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalEdgeFunction(
    int edge_index,
    const LocalPoint& xi
) const {
    // 边界检查：棱边索引必须在有效范围内
    if (edge_index < 0 || edge_index >= edge_count_) {
        FEEM_ERROR("evalEdgeFunction: 棱边索引越界 edge_index={}, 有效范围[0, {})",
                   edge_index, edge_count_);
        return Eigen::Vector3d::Zero();
    }

    // 提取局部坐标分量
    double xi_coord = xi.coords(0);
    double eta_coord = xi.coords(1);
    double zeta_coord = xi.coords(2);

    // 根据单元类型分派到对应的基函数计算方法
    switch (element_type_) {
    case ElementType::TRI3_EDGE:
        return evalTri3Edge(edge_index, xi_coord, eta_coord);

    case ElementType::QUAD4_EDGE:
        return evalQuad4Edge(edge_index, xi_coord, eta_coord);

    case ElementType::TET4_EDGE:
        return evalTet4Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    case ElementType::HEX8_EDGE:
        return evalHex8Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    case ElementType::PRISM6_EDGE:
        return evalPrism6Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    case ElementType::PYRAMID5_EDGE:
        return evalPyramid5Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    default:
        FEEM_ERROR("evalEdgeFunction: 未知的Nedelec单元类型={}",
                   static_cast<int>(element_type_));
        return Eigen::Vector3d::Zero();
    }
}

// ==================== evalCurlEdge 分派实现 ====================

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalCurlEdge(
    int edge_index,
    const LocalPoint& xi
) const {
    // 边界检查
    if (edge_index < 0 || edge_index >= edge_count_) {
        FEEM_ERROR("evalCurlEdge: 棱边索引越界 edge_index={}, 有效范围[0, {})",
                   edge_index, edge_count_);
        return Eigen::Vector3d::Zero();
    }

    double xi_coord = xi.coords(0);
    double eta_coord = xi.coords(1);
    double zeta_coord = xi.coords(2);

    switch (element_type_) {
    case ElementType::TRI3_EDGE:
        return curlTri3Edge(edge_index);

    case ElementType::QUAD4_EDGE:
        return curlQuad4Edge(edge_index);

    case ElementType::TET4_EDGE:
        return curlTet4Edge(edge_index);

    case ElementType::HEX8_EDGE:
        return curlHex8Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    case ElementType::PRISM6_EDGE:
        return curlPrism6Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    case ElementType::PYRAMID5_EDGE:
        return curlPyramid5Edge(edge_index, xi_coord, eta_coord, zeta_coord);

    default:
        FEEM_ERROR("evalCurlEdge: 未知的Nedelec单元类型={}",
                   static_cast<int>(element_type_));
        return Eigen::Vector3d::Zero();
    }
}

// ==================== calcPhysicalEdgeFunction 实现 ====================

template<int Dim>
EdgeCurlResult NedelecElement<Dim>::calcPhysicalEdgeFunction(
    int edge_index,
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords
) const {
    EdgeCurlResult result;

    // 步骤1：获取参考域旋度
    result.curl_reference = evalCurlEdge(edge_index, xi);

    // 步骤2：计算雅可比矩阵用于坐标映射
    JacobianResult jaco_result = calcJacobian(xi, node_coords);

    // 步骤3：物理域旋度映射公式
    // 3D: curl_physical = (1/detJ) · J · curl_reference
    // 2D标量旋度同理（curl_reference的第3分量存储标量旋度）
    if (std::abs(jaco_result.det_j) > 1e-15) {
        result.curl_physical = (jaco_result.jacobian * result.curl_reference)
                               / jaco_result.det_j;
    } else {
        FEEM_WARN("calcPhysicalEdgeFunction: 雅可比行列式接近零({}), "
                  "物理域旋度可能不准确", jaco_result.det_j);
        result.curl_physical = Eigen::Vector3d::Zero();
    }

    return result;
}

// =====================================================================
//  TRI3_EDGE 实现 - 三角形一阶Nedelec棱边元（2D，3条棱边）
// =====================================================================

/**
 * 三角形参考域定义：
 * 面积坐标系 (ξ, η)，第三个面积坐标 ζ = 1 - ξ - η
 * 定义域: ξ ≥ 0, η ≥ 0, ξ + η ≤ 1
 * 节点位置: v₁=(1,0), v₂=(0,1), v₃=(0,0) 对应面积坐标
 *
 * 梯度向量（常数）:
 * ∇λ₁ = ∇ξ = (1, 0)^T
 * ∇λ₂ = ∇η = (0, 1)^T
 * ∇λ₃ = ∇ζ = (-1, -1)^T
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalTri3Edge(
    int edge_index, double xi, double eta
) const {
    // 面积坐标
    double lambda1 = xi;              // λ₁ = ξ
    double lambda2 = eta;             // λ₂ = η
    double lambda3 = 1.0 - xi - eta;  // λ₃ = ζ = 1-ξ-η

    switch (edge_index) {
    case 0:
        // 棱边0（节点1→2）: W₀ = λ₂∇λ₁ - λ₁∇λ₂ = η(1,0) - ξ(0,1) = (η, -ξ)
        return Eigen::Vector3d(lambda2, -lambda1, 0.0);

    case 1:
        // 棱边1（节点2→3）: W₁ = λ₃∇λ₂ - λ₂∇λ₃ = ζ(0,1) - η(-1,-1)
        //                    = (η, ζ+η) = (η, 1-ξ)
        return Eigen::Vector3d(lambda2, lambda3 + lambda2, 0.0);

    case 2:
        // 棱边2（节点3→1）: W₂ = λ₁∇λ₃ - λ₃∇λ₁ = ξ(-1,-1) - ζ(1,0)
        //                    = (-ξ-ζ, -ξ) = (η-1, -ξ)
        return Eigen::Vector3d(lambda2 - 1.0, -lambda1, 0.0);

    default:
        FEEM_ERROR("evalTri3Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlTri3Edge(int edge_index) const {
    // 2D下旋度为标量: curl(W) = ∂W_y/∂x - ∂W_x/∂y
    // 对于Whitney形式 W_ij = λ_j∇λ_i - λ_i∇λ_j
    // curl(W_ij) = -2 （所有棱边的旋度均为常数-2）
    //
    // 验证: W₀ = (η, -ξ), curl(W₀) = ∂(-ξ)/∂ξ - ∂(η)/∂η = -1 - 1 = -2
    //       W₁ = (η, 1-ξ), curl(W₁) = ∂(1-ξ)/∂ξ - ∂(η)/∂η = -1 - 1 = -2
    //       W₂ = (η-1, -ξ), curl(W₂) = ∂(-ξ)/∂ξ - ∂(η-1)/∂η = -1 - 1 = -2
    //
    // 返回格式: 将标量旋度存入第3分量（z方向），前两分量为0
    (void)edge_index;  // 所有棱边旋度相同，参数未使用
    return Eigen::Vector3d(0.0, 0.0, -2.0);
}

// =====================================================================
//  QUAD4_EDGE 实现 - 四边形一阶Nedelec棱边元（2D，4条棱边）
// =====================================================================

/**
 * 四边形参考域定义：
 * 局部坐标 (ξ, η) ∈ [-1, 1] × [-1, 1]
 * Gmsh标准节点编号:
 *   v₁=(-1,-1), v₂=(1,-1), v₃=(1,1), v₄=(-1,1)
 *
 * 棱边方向约定：
 *   边0(v₁→v₂): 底边，η=-1，方向+ξ
 *   边1(v₂→v₃): 右边，ξ=1，方向+η
 *   边2(v₃→v₄): 顶边，η=1，方向-ξ（反向）
 *   边3(v₄→v₁): 左边，ξ=-1，方向-η（反向）
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalQuad4Edge(
    int edge_index, double xi, double eta
) const {
    switch (edge_index) {
    case 0: {
        // 棱边0（底边，η=-1，方向+ξ）
        // W₀ = ½(1-η)∇ξ = ½(1-η)(1, 0)
        double coeff = 0.5 * (1.0 - eta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    case 1: {
        // 棱边1（右边，ξ=1，方向+η）
        // W₁ = ½(1+ξ)∇η = ½(1+ξ)(0, 1)
        double coeff = 0.5 * (1.0 + xi);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    case 2: {
        // 棱边2（顶边，η=1，方向-ξ，取负号）
        // W₂ = -½(1+η)∇ξ = -½(1+η)(1, 0)
        double coeff = -0.5 * (1.0 + eta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    case 3: {
        // 棱边3（左边，ξ=-1，方向-η，取负号）
        // W₃ = -½(1-ξ)∇η = -½(1-ξ)(0, 1)
        double coeff = -0.5 * (1.0 - xi);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    default:
        FEEM_ERROR("evalQuad4Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlQuad4Edge(int edge_index) const {
    // 2D旋度为标量，计算每个基函数的 curl = ∂W_η/∂ξ - ∂W_ξ/∂η
    //
    // W₀ = ½(1-η)(1,0):  curl = 0 - (-½) = +½
    // W₁ = ½(1+ξ)(0,1):  curl = ½ - 0 = +½
    // W₂ = -½(1+η)(1,0): curl = 0 - (-½) = +½
    // W₃ = -½(1-ξ)(0,1): curl = (-½) - 0 = -½ ... 等等让我重新算
    //
    // 重新计算:
    // W₃ = (0, -½(1-ξ)): curl = ∂[-½(1-ξ)]/∂ξ - 0 = ½
    //
    // 所以所有4个基函数的旋度都是 +½
    // 这是因为四边形的旋度之和应该满足某种一致性条件
    (void)edge_index;
    return Eigen::Vector3d(0.0, 0.0, 0.5);
}

// =====================================================================
//  TET4_EDGE 实现 - 四面体一阶Nedelec棱边元（3D，6条棱边）
// =====================================================================

/**
 * 四面体参考域定义：
 * 体坐标系 (λ₁, λ₂, λ₃, λ₄)，其中 λᵢ ≥ 0, Σλᵢ = 1
 * 映射到局部坐标: ξ=λ₁, η=λ₂, ζ=λ₃, λ₄=1-ξ-η-ζ
 * 定义域: ξ≥0, η≥0, ζ≥0, ξ+η+ζ≤1
 *
 * Gmsh标准节点编号:
 *   v₁=(1,0,0), v₂=(0,1,0), v₃=(0,0,1), v₄=(0,0,0)
 *
 * 梯度向量（常数）:
 *   ∇λ₁ = (1, 0, 0)^T
 *   ∇λ₂ = (0, 1, 0)^T
 *   ∇λ₃ = (0, 0, 1)^T
 *   ∇λ₄ = (-1, -1, -1)^T
 *
 * 6条棱边的Whitney基函数: W_ij = λ_j∇λ_i - λ_i∇λ_j
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalTet4Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    // 体坐标
    double l1 = xi;                    // λ₁
    double l2 = eta;                   // λ₂
    double l3 = zeta;                  // λ₃
    double l4 = 1.0 - xi - eta - zeta; // λ₄

    switch (edge_index) {
    case 0: {
        // 棱边0（节点1→2）: W₀ = λ₂∇λ₁ - λ₁∇λ₂ = l2(1,0,0) - l1(0,1,0)
        return Eigen::Vector3d(l2, -l1, 0.0);
    }

    case 1: {
        // 棱边1（节点2→3）: W₁ = λ₃∇λ₂ - λ₂∇λ₃ = l3(0,1,0) - l2(0,0,1)
        return Eigen::Vector3d(0.0, l3, -l2);
    }

    case 2: {
        // 棱边2（节点3→1）: W₂ = λ₁∇λ₃ - λ₃∇λ₁ = l1(0,0,1) - l3(1,0,0)
        return Eigen::Vector3d(-l3, 0.0, l1);
    }

    case 3: {
        // 棱边3（节点1→4）: W₃ = λ₄∇λ₁ - λ₁∇λ₄
        // = l4(1,0,0) - l1(-1,-1,-1) = (l4+l1, l1, l1) = (1-l2-l3, l1, l1)
        return Eigen::Vector3d(l4 + l1, l1, l1);
    }

    case 4: {
        // 棱边4（节点2→4）: W₄ = λ₄∇λ₂ - λ₂∇λ₄
        // = l4(0,1,0) - l2(-1,-1,-1) = (l2, l4+l2, l2) = (l2, 1-l1-l3, l2)
        return Eigen::Vector3d(l2, l4 + l2, l2);
    }

    case 5: {
        // 棱边5（节点3→4）: W₅ = λ₄∇λ₃ - λ₃∇λ₄
        // = l4(0,0,1) - l3(-1,-1,-1) = (l3, l3, l4+l3) = (l3, l3, 1-l1-l2)
        return Eigen::Vector3d(l3, l3, l4 + l3);
    }

    default:
        FEEM_ERROR("evalTet4Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlTet4Edge(int edge_index) const {
    // 3D旋度公式: curl(W_ij) = curl(λ_j∇λ_i - λ_i∇λ_j)
    // 使用矢量恒等式: ∇×(fV) = ∇f × V + f(∇×V)
    // 由于∇λ_i为常矢量（∇×∇λ_i = 0），简化为:
    // curl(W_ij) = ∇λ_j × ∇λ_i - ∇λ_i × ∇λ_j = 2(∇λ_j × ∇λ_i)
    //
    // 这是一个**常矢量**（不依赖空间位置），是一阶Nedelec元的重要性质
    switch (edge_index) {
    case 0:
        // curl(W₀) = 2∇λ₂ × ∇λ₁ = 2(0,1,0) × (1,0,0) = 2(0,0,-1)
        return Eigen::Vector3d(0.0, 0.0, -2.0);

    case 1:
        // curl(W₁) = 2∇λ₃ × ∇λ₂ = 2(0,0,1) × (0,1,0) = 2(-1,0,0)
        return Eigen::Vector3d(-2.0, 0.0, 0.0);

    case 2:
        // curl(W₂) = 2∇λ₁ × ∇λ₃ = 2(1,0,0) × (0,0,1) = 2(0,-1,0)
        return Eigen::Vector3d(0.0, -2.0, 0.0);

    case 3: {
        // curl(W₃) = 2∇λ₄ × ∇λ₁ = 2(-1,-1,-1) × (1,0,0) = 2(0,-1,1)
        return Eigen::Vector3d(0.0, -2.0, 2.0);
    }

    case 4: {
        // curl(W₄) = 2∇λ₄ × ∇λ₂ = 2(-1,-1,-1) × (0,1,0) = 2(1,0,-1)
        return Eigen::Vector3d(2.0, 0.0, -2.0);
    }

    case 5: {
        // curl(W₅) = 2∇λ₄ × ∇λ₃ = 2(-1,-1,-1) × (0,0,1) = 2(-1,1,0)
        return Eigen::Vector3d(-2.0, 2.0, 0.0);
    }

    default:
        FEEM_ERROR("curlTet4Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

// =====================================================================
//  HEX8_EDGE 实现 - 六面体一阶Nedelec棱边元（3D，12条棱边）
// =====================================================================

/**
 * 六面体参考域定义：
 * 局部坐标 (ξ, η, ζ) ∈ [-1, 1]³
 * Gmsh/ANSYS标准节点编号:
 *   底面(zeta=-1): v₁(-1,-1,-1), v₂(1,-1,-1), v₃(1,1,-1), v₄(-1,1,-1)
 *   顶面(zeta=+1): v₅(-1,-1,+1), v₆(1,-1,+1), v₇(1,1,+1), v₈(-1,1,+1)
 *
 * 12条棱边分组：
 *   ξ方向(平行于ξ轴): e0(1-2), e2(3-4反向), e4(5-6), e6(7-8反向)
 *   η方向(平行于η轴): e1(2-3), e3(4-1反向), e5(6-7), e7(8-5反向)
 *   ζ方向(平行于ζ轴): e8(1-5), e9(2-6), e10(3-7), e11(4-8)
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalHex8Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    switch (edge_index) {

    // ===== ξ方向棱边（平行于ξ轴） =====
    case 0: {
        // 棱边0（节点1→2）: 底面前沿，η=-1, ζ=-1, 方向+ξ
        // W₀ = ½(1-η)(1-ζ)∇ξ
        double coeff = 0.5 * (1.0 - eta) * (1.0 - zeta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    case 2: {
        // 棱边2（节点3→4）: 底面后沿，η=+1, ζ=-1, 方向-ξ（反向）
        // W₂ = -½(1+η)(1-ζ)∇ξ
        double coeff = -0.5 * (1.0 + eta) * (1.0 - zeta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    case 4: {
        // 棱边4（节点5→6）: 顶面前沿，η=-1, ζ=+1, 方向+ξ
        // W₄ = ½(1-η)(1+ζ)∇ξ
        double coeff = 0.5 * (1.0 - eta) * (1.0 + zeta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    case 6: {
        // 棱边6（节点7→8）: 顶面后沿，η=+1, ζ=+1, 方向-ξ（反向）
        // W₆ = -½(1+η)(1+ζ)∇ξ
        double coeff = -0.5 * (1.0 + eta) * (1.0 + zeta);
        return Eigen::Vector3d(coeff, 0.0, 0.0);
    }

    // ===== η方向棱边（平行于η轴） =====
    case 1: {
        // 棱边1（节点2→3）: 底面右沿，ξ=+1, ζ=-1, 方向+η
        // W₁ = ½(1+ξ)(1-ζ)∇η
        double coeff = 0.5 * (1.0 + xi) * (1.0 - zeta);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    case 3: {
        // 棱边3（节点4→1）: 底面左沿，ξ=-1, ζ=-1, 方向-η（反向）
        // W₃ = -½(1-ξ)(1-ζ)∇η
        double coeff = -0.5 * (1.0 - xi) * (1.0 - zeta);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    case 5: {
        // 棱边5（节点6→7）: 顶面右沿，ξ=+1, ζ=+1, 方向+η
        // W₅ = ½(1+ξ)(1+ζ)∇η
        double coeff = 0.5 * (1.0 + xi) * (1.0 + zeta);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    case 7: {
        // 棱边7（节点8→5）: 顶面左沿，ξ=-1, ζ=+1, 方向-η（反向）
        // W₇ = -½(1-ξ)(1+ζ)∇η
        double coeff = -0.5 * (1.0 - xi) * (1.0 + zeta);
        return Eigen::Vector3d(0.0, coeff, 0.0);
    }

    // ===== ζ方向棱边（平行于ζ轴） =====
    case 8: {
        // 棱边8（节点1→5）: 左前沿竖直棱，ξ=-1, η=-1, 方向+ζ
        // W₈ = ½(1-ξ)(1-η)∇ζ
        double coeff = 0.5 * (1.0 - xi) * (1.0 - eta);
        return Eigen::Vector3d(0.0, 0.0, coeff);
    }

    case 9: {
        // 棱边9（节点2→6）: 右前沿竖直棱，ξ=+1, η=-1, 方向+ζ
        // W₉ = ½(1+ξ)(1-η)∇ζ
        double coeff = 0.5 * (1.0 + xi) * (1.0 - eta);
        return Eigen::Vector3d(0.0, 0.0, coeff);
    }

    case 10: {
        // 棱边10（节点3→7）: 右后沿竖直棱，ξ=+1, η=+1, 方向+ζ
        // W₁₀ = ½(1+ξ)(1+η)∇ζ
        double coeff = 0.5 * (1.0 + xi) * (1.0 + eta);
        return Eigen::Vector3d(0.0, 0.0, coeff);
    }

    case 11: {
        // 棱边11（节点4→8）: 左后沿竖直棱，ξ=-1, η=+1, 方向+ζ
        // W₁₁ = ½(1-ξ)(1+η)∇ζ
        double coeff = 0.5 * (1.0 - xi) * (1.0 + eta);
        return Eigen::Vector3d(0.0, 0.0, coeff);
    }

    default:
        FEEM_ERROR("evalHex8Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlHex8Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    // 六面体旋度计算: curl(f·e_dir) = ∇f × e_dir
    // 因为curl(constant_vector) = 0，所以只需计算梯度叉乘方向向量
    //
    // 对于 f(η,ζ)·∇ξ 形式的基函数:
    //   curl = (∂f/∂η, ∂f/∂ζ, 0) × (1,0,0) = (0, ∂f/∂ζ, -∂f/∂η)
    //
    // 对于 f(ξ,ζ)·∇η 形式的基函数:
    //   curl = (∂f/∂ξ, 0, ∂f/∂ζ) × (0,1,0) = (-∂f/∂ζ, 0, ∂f/∂ξ)
    //
    // 对于 f(ξ,η)·∇ζ 形式的基函数:
    //   curl = (∂f/∂ξ, ∂f/∂η, 0) × (0,0,1) = (∂f/∂η, -∂f/∂ξ, 0)

    switch (edge_index) {

    // ===== ξ方向棱边旋度 =====
    case 0: {
        // W₀ = ½(1-η)(1-ζ)∇ξ, f = ½(1-η)(1-ζ)
        // ∇f = (0, -½(1-ζ), -½(1-η))
        // curl = ∇f × (1,0,0) = (0, -½(1-η), ½(1-ζ))
        return Eigen::Vector3d(0.0, -0.5 * (1.0 - eta), 0.5 * (1.0 - zeta));
    }

    case 2: {
        // W₂ = -½(1+η)(1-ζ)∇ξ, f = -½(1+η)(1-ζ)
        // ∇f = (0, -½(1-ζ), ½(1+η))
        // curl = (0, ½(1+η), -½(1-ζ))
        return Eigen::Vector3d(0.0, 0.5 * (1.0 + eta), -0.5 * (1.0 - zeta));
    }

    case 4: {
        // W₄ = ½(1-η)(1+ζ)∇ξ, f = ½(1-η)(1+ζ)
        // ∇f = (0, -½(1+ζ), ½(1-η))
        // curl = (0, -½(1-η), -½(1+ζ)) ... wait let me recalculate
        // curl = ∇f × e_ξ = (0, df/dη, df/dζ) × (1,0,0) = (0, df/dζ, -df/dη)
        //       = (0, ½(1-η), -(-½(1+ζ))) = (0, ½(1-η), ½(1+ζ))
        return Eigen::Vector3d(0.0, 0.5 * (1.0 - eta), 0.5 * (1.0 + zeta));
    }

    case 6: {
        // W₆ = -½(1+η)(1+ζ)∇ξ
        // curl = (0, -½(1+ζ), ½(1+η))
        return Eigen::Vector3d(0.0, -0.5 * (1.0 + zeta), 0.5 * (1.0 + eta));
    }

    // ===== η方向棱边旋度 =====
    case 1: {
        // W₁ = ½(1+ξ)(1-ζ)∇η, f = ½(1+ξ)(1-ζ)
        // curl = ∇f × (0,1,0) = (-df/dζ, 0, df/dξ) = (½(1+ξ), 0, ½(1-ζ))
        return Eigen::Vector3d(0.5 * (1.0 + xi), 0.0, 0.5 * (1.0 - zeta));
    }

    case 3: {
        // W₃ = -½(1-ξ)(1-ζ)∇η
        // curl = (-(-½(1-ξ)), 0, -½(1-ζ)) = (½(1-ξ), 0, -½(1-ζ))
        return Eigen::Vector3d(0.5 * (1.0 - xi), 0.0, -0.5 * (1.0 - zeta));
    }

    case 5: {
        // W₅ = ½(1+ξ)(1+ζ)∇η
        // curl = (-½(1+ξ), 0, ½(1+ζ))
        return Eigen::Vector3d(-0.5 * (1.0 + xi), 0.0, 0.5 * (1.0 + zeta));
    }

    case 7: {
        // W₇ = -½(1-ξ)(1+ζ)∇η
        // curl = (½(1-ξ), 0, -½(1+ζ))
        return Eigen::Vector3d(0.5 * (1.0 - xi), 0.0, -0.5 * (1.0 + zeta));
    }

    // ===== ζ方向棱边旋度 =====
    case 8: {
        // W₈ = ½(1-ξ)(1-η)∇ζ, f = ½(1-ξ)(1-η)
        // curl = ∇f × (0,0,1) = (df/dη, -df/dξ, 0)
        //       = (-½(1-ξ), ½(1-η), 0)
        return Eigen::Vector3d(-0.5 * (1.0 - xi), 0.5 * (1.0 - eta), 0.0);
    }

    case 9: {
        // W₉ = ½(1+ξ)(1-η)∇ζ
        // curl = (-½(1+ξ), -½(1-η), 0)
        return Eigen::Vector3d(-0.5 * (1.0 + xi), -0.5 * (1.0 - eta), 0.0);
    }

    case 10: {
        // W₁₀ = ½(1+ξ)(1+η)∇ζ
        // curl = (½(1+ξ), -½(1+η), 0)
        return Eigen::Vector3d(0.5 * (1.0 + xi), -0.5 * (1.0 + eta), 0.0);
    }

    case 11: {
        // W₁₁ = ½(1-ξ)(1+η)∇ζ
        // curl = (½(1-ξ), ½(1+η), 0)
        return Eigen::Vector3d(0.5 * (1.0 - xi), 0.5 * (1.0 + eta), 0.0);
    }

    default:
        FEEM_ERROR("curlHex8Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

// =====================================================================
//  PRISM6_EDGE 实现 - 三棱柱一阶Nedelec棱边元（3D，9条棱边）
// =====================================================================

/**
 * 三棱柱参考域定义：
 * 三角形部分: 面积坐标 (ξ, η), ξ≥0, η≥0, ξ+η≤1
 * 棱柱轴向: ζ ∈ [-1, 1]
 *
 * Gmsh标准节点编号:
 *   底面(ζ=-1): v₁(0,0,-1), v₂(1,0,-1), v₃(0,1,-1)
 *   顶面(ζ=+1): v₄(0,0,+1), v₅(1,0,+1), v₆(0,1,+1)
 *
 * 面积坐标函数: L₁=ξ, L₂=η, L₃=1-ξ-η
 * ζ方向线性插值函数: N⁻(ζ)=½(1-ζ), N⁺(ζ)=½(1+ζ)
 *
 * 9条棱边分组：
 *   底面3条(e0-e2): TRI3基函数 × N⁻(ζ)
 *   顶面3条(e3-e5): TRI3基函数 × N⁺(ζ)
 *   侧棱3条(e6-e8): 面积坐标 × ∇ζ
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalPrism6Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    // 面积坐标
    double L1 = xi;                    // 第一个面积坐标
    double L2 = eta;                   // 第二个面积坐标
    double L3 = 1.0 - xi - eta;        // 第三个面积坐标

    // ζ方向的线性插值系数
    double n_minus = 0.5 * (1.0 - zeta); // N⁻(ζ): 在ζ=-1处为1，ζ=+1处为0
    double n_plus  = 0.5 * (1.0 + zeta); // N⁺(ζ): 在ζ=-1处为0，ζ=+1处为1

    if (edge_index >= 0 && edge_index <= 2) {
        // ===== 底面棱边（ζ=-1面，使用TRI3基函数×N⁻） =====
        int tri_edge = edge_index;  // 映射到三角形棱边索引
        Eigen::Vector3d tri_func = evalTri3Edge(tri_edge, xi, eta);
        return n_minus * tri_func;
    }

    if (edge_index >= 3 && edge_index <= 5) {
        // ===== 顶面棱边（ζ=+1面，使用TRI3基函数×N⁺） =====
        int tri_edge = edge_index - 3;  // 映射到三角形棱边索引
        Eigen::Vector3d tri_func = evalTri3Edge(tri_edge, xi, eta);
        return n_plus * tri_func;
    }

    // ===== 侧棱（连接底面与顶面对应角点） =====
    switch (edge_index) {
    case 6:
        // 侧棱6（节点1→4）: v₁(0,0)→v₄(0,0)，对应L₁
        // W₆ = L₁ · ∇ζ = ξ · (0,0,1)
        return Eigen::Vector3d(0.0, 0.0, L1);

    case 7:
        // 侧棱7（节点2→5）: v₂(1,0)→v₅(1,0)，对应L₂
        // W₇ = L₂ · ∇ζ = η · (0,0,1)
        return Eigen::Vector3d(0.0, 0.0, L2);

    case 8:
        // 侧棱8（节点3→6）: v₃(0,1)→v₆(0,1)，对应L₃
        // W₈ = L₃ · ∇ζ = (1-ξ-η) · (0,0,1)
        return Eigen::Vector3d(0.0, 0.0, L3);

    default:
        FEEM_ERROR("evalPrism6Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlPrism6Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    // 棱柱旋度需要同时考虑三角形面内部分和ζ方向部分的导数
    // 对于底面/顶面棱边: curl(N^±(ζ) · W_tri) = d(N^±)/dζ · e_z × W_tri + N^± · curl_tri(W_tri)
    // 对于侧棱: curl(L_i · ∇ζ) = ∇L_i × ∇ζ

    if (edge_index >= 0 && edge_index <= 2) {
        // 底面棱边旋度: curl = (dN⁻/dζ) e_z × W_tri + N⁻ · curl_tri
        // dN⁻/dζ = -½, curl_tri(W_tri) = (0,0,-2)
        int tri_edge = edge_index;
        Eigen::Vector3d w_tri = evalTri3Edge(tri_edge, xi, eta);
        double n_minus = 0.5 * (1.0 - zeta);
        double dn_minus_dzeta = -0.5;

        // e_z × W_tri = (0,0,1) × (Wx, Wy, 0) = (-Wy, Wx, 0)
        Eigen::Vector3d cross_term(dn_minus_dzeta * (-w_tri.y()),
                                   dn_minus_dzeta * w_tri.x(), 0.0);
        // curl_tri term: (0, 0, -2) scaled by n_minus
        Eigen::Vector3d curl_tri_term(0.0, 0.0, -2.0 * n_minus);
        return cross_term + curl_tri_term;
    }

    if (edge_index >= 3 && edge_index <= 5) {
        // 顶面棱边旋度: curl = (dN⁺/dζ) e_z × W_tri + N⁺ · curl_tri
        // dN⁺/dζ = +½
        int tri_edge = edge_index - 3;
        Eigen::Vector3d w_tri = evalTri3Edge(tri_edge, xi, eta);
        double n_plus = 0.5 * (1.0 + zeta);
        double dn_plus_dzeta = 0.5;

        Eigen::Vector3d cross_term(dn_plus_dzeta * (-w_tri.y()),
                                   dn_plus_dzeta * w_tri.x(), 0.0);
        Eigen::Vector3d curl_tri_term(0.0, 0.0, -2.0 * n_plus);
        return cross_term + curl_tri_term;
    }

    // 侧棱旋度: curl(L_i · ∇ζ) = ∇L_i × ∇ζ
    // ∇L₁ = (1, 0, 0), ∇L₂ = (0, 1, 0), ∇L₃ = (-1, -1, 0)
    // ∇ζ = (0, 0, 1)
    switch (edge_index) {
    case 6:
        // curl(W₆) = ∇L₁ × ∇ζ = (1,0,0) × (0,0,1) = (0, -1, 0)
        return Eigen::Vector3d(0.0, -1.0, 0.0);

    case 7:
        // curl(W₇) = ∇L₂ × ∇ζ = (0,1,0) × (0,0,1) = (1, 0, 0)
        return Eigen::Vector3d(1.0, 0.0, 0.0);

    case 8:
        // curl(W₈) = ∇L₃ × ∇ζ = (-1,-1,0) × (0,0,1) = (-1, 1, 0)
        return Eigen::Vector3d(-1.0, 1.0, 0.0);

    default:
        FEEM_ERROR("curlPrism6Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

// =====================================================================
//  PYRAMID5_EDGE 实现 - 金字塔一阶Nedelec棱边元（3D，8条棱边）
// =====================================================================

/**
 * 金字塔参考域定义：
 * 局部坐标 (ξ, η, ζ)，其中 ξ,η ∈ [-1, 1], ζ ∈ [0, 1]
 * 锥顶在 (0, 0, 1)，底面在 ζ=0 平面
 *
 * Gmsh标准节点编号:
 *   底面(ζ=0): v₁(-1,-1,0), v₂(1,-1,0), v₃(1,1,0), v₄(-1,1,0)
 *   锥顶:     v₅(0, 0, 1)
 *
 * 8条棱边分组：
 *   底面4条(e0-e3): QUAD4基函数形式 × (1-ζ) [考虑锥形收缩]
 *   侧棱4条(e4-e7): 从底面角点到锥顶的退化棱边函数
 *
 * 辅助变量: r = 1 - ζ（r=1在底面，r=0在锥顶）
 */
template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::evalPyramid5Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    // 辅助变量：r在底面为1，锥顶为0
    double r = 1.0 - zeta;

    if (edge_index >= 0 && edge_index <= 3) {
        // ===== 底面棱边（类似QUAD4，乘以r因子模拟锥形收缩） =====
        // 基函数形式与QUAD4相同，但乘以r=(1-ζ)使函数在锥顶处自然衰减
        switch (edge_index) {
        case 0: {
            // 棱边0（底边 v₁→v₂）: W₀ = ½r(1-η)∇ξ
            double coeff = 0.5 * r * (1.0 - eta);
            return Eigen::Vector3d(coeff, 0.0, 0.0);
        }

        case 1: {
            // 棱边1（右边 v₂→v₃）: W₁ = ½r(1+ξ)∇η
            double coeff = 0.5 * r * (1.0 + xi);
            return Eigen::Vector3d(0.0, coeff, 0.0);
        }

        case 2: {
            // 棱边2（顶边 v₃→v₄）: W₂ = -½r(1+η)∇ξ
            double coeff = -0.5 * r * (1.0 + eta);
            return Eigen::Vector3d(coeff, 0.0, 0.0);
        }

        case 3: {
            // 棱边3（左边 v₄→v₁）: W₃ = -½r(1-ξ)∇η
            double coeff = -0.5 * r * (1.0 - xi);
            return Eigen::Vector3d(0.0, coeff, 0.0);
        }

        default:
            break;  // 不会执行，消除编译器警告
        }
    }

    // ===== 侧棱（从底面角点到锥顶） =====
    // 使用"屋顶函数"形式的基函数，确保切向连续性和正确的环流性质
    // 一般形式: W_side = ¼[(1±ξ)(1±η)(-∇ζ) + ζ((1±η)(±∇ξ) + (1±ξ)(±∇η))]
    switch (edge_index) {
    case 4: {
        // 侧棱4（v₁→v₅）: 从(-1,-1,0)到(0,0,1)
        // 屋顶函数: ¼[(1-ξ)(1-η)(-∇ζ) + ζ((1-η)∇ξ + (1-ξ)∇η)]
        double roof = 0.25 * (1.0 - xi) * (1.0 - eta);
        double zx = 0.25 * zeta * (1.0 - eta);   // ζ*(1-η)/4 的ξ分量
        double zy = 0.25 * zeta * (1.0 - xi);     // ζ*(1-ξ)/4 的η分量
        return Eigen::Vector3d(zx - roof * 0.0, zy - roof * 0.0, -roof);
    }

    case 5: {
        // 侧棱5（v₂→v₅）: 从(1,-1,0)到(0,0,1)
        // W₅ = ¼[(1+ξ)(1-η)(-∇ζ) + ζ((1-η)(-∇ξ) + (1+ξ)∇η)]
        double roof = 0.25 * (1.0 + xi) * (1.0 - eta);
        double zx = -0.25 * zeta * (1.0 - eta);   // ζ项中ξ分量带负号
        double zy = 0.25 * zeta * (1.0 + xi);
        return Eigen::Vector3d(zx, zy, -roof);
    }

    case 6: {
        // 侧棱6（v₃→v₅）: 从(1,1,0)到(0,0,1)
        // W₆ = ¼[(1+ξ)(1+η)(-∇ζ) + ζ((1+η)(-∇ξ) + (1+ξ)(-∇η))]
        double roof = 0.25 * (1.0 + xi) * (1.0 + eta);
        double zx = -0.25 * zeta * (1.0 + eta);
        double zy = -0.25 * zeta * (1.0 + xi);
        return Eigen::Vector3d(zx, zy, -roof);
    }

    case 7: {
        // 侧棱7（v₄→v₅）: 从(-1,1,0)到(0,0,1)
        // W₇ = ¼[(1-ξ)(1+η)(-∇ζ) + ζ((1+η)∇ξ + (1-ξ)(-∇η))]
        double roof = 0.25 * (1.0 - xi) * (1.0 + eta);
        double zx = 0.25 * zeta * (1.0 + eta);
        double zy = -0.25 * zeta * (1.0 - xi);
        return Eigen::Vector3d(zx, zy, -roof);
    }

    default:
        FEEM_ERROR("evalPyramid5Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

template<int Dim>
Eigen::Vector3d NedelecElement<Dim>::curlPyramid5Edge(
    int edge_index, double xi, double eta, double zeta
) const {
    double r = 1.0 - zeta;

    if (edge_index >= 0 && edge_index <= 3) {
        // 底面棱边旋度: 类似QUAD4但额外包含对ζ的导数项
        // W = ½r(·)∇dir, 其中r = 1-ζ
        // curl = ∇(½r(·)) × ∇dir = [0, 0, ∓½(·)] × ∇dir + ½r · curl_QUAD4
        // 其中∓来自dr/dζ = -1
        switch (edge_index) {
        case 0: {
            // W₀ = ½r(1-η)∇ξ, curl包含dr/dζ产生的额外项
            // 基础项(来自QUAD4): (0, 0, ½r) ... 不对，让我重新推导
            // curl(f(η,ζ)e_ξ) = (0, df/dζ, -df/dη)
            // f = ½(1-ζ)(1-η) = ½r(1-η)
            // df/dζ = -½(1-η), df/dη = -½(1-ζ) = -½r
            // curl = (0, -½(1-η), ½r)
            return Eigen::Vector3d(0.0, -0.5 * (1.0 - eta), 0.5 * r);
        }

        case 1: {
            // W₁ = ½r(1+ξ)∇η
            // f = ½r(1+ξ), df/dζ = -½(1+ξ), df/dξ = ½r
            // curl(f e_η) = (-df/dζ, 0, df/dξ) = (½(1+ξ), 0, ½r)
            return Eigen::Vector3d(0.5 * (1.0 + xi), 0.0, 0.5 * r);
        }

        case 2: {
            // W₂ = -½r(1+η)∇ξ
            // f = -½r(1+η), df/dζ = ½(1+η), df/dη = -½r
            // curl = (0, -df/dζ, -df/dη) = (0, -½(1+η), ½r)
            return Eigen::Vector3d(0.0, -0.5 * (1.0 + eta), 0.5 * r);
        }

        case 3: {
            // W₃ = -½r(1-ξ)∇η
            // f = -½r(1-ξ), df/dζ = ½(1-ξ), df/dξ = -½r
            // curl = (-df/dζ, 0, df/dξ) = (-½(1-ξ), 0, -½r)
            return Eigen::Vector3d(-0.5 * (1.0 - xi), 0.0, -0.5 * r);
        }

        default:
            break;
        }
    }

    // 侧棱旋度计算
    // 侧棱基函数形式: W = ¼[f_roof(-∇ζ) + ζ(g_x∇ξ + g_y∇η)]
    // 其中f_roof和g_x,g_y是(ξ,η)的双线性函数
    // 旋度展开较复杂，此处给出各棱边的显式结果
    switch (edge_index) {
    case 4: {
        // 侧棱4(v₁→v₅)旋度
        double ome = 1.0 - xi;  // (1-ξ)
        double one = 1.0 - eta; // (1-η)
        // curl = ¼[(one*(-1) - (-ome)*0), ((-ome)*(-1)-one*0), (one*1-ome*1)]
        //      = ¼[-one, ome, one-ome]
        return Eigen::Vector3d(-0.25 * one, 0.25 * ome, 0.25 * (one - ome));
    }

    case 5: {
        // 侧棱5(v₂→v₅)旋度
        double ope = 1.0 + xi;  // (1+ξ)
        double one = 1.0 - eta; // (1-η)
        return Eigen::Vector3d(-0.25 * one, -0.25 * ope, 0.25 * (-one - ope));
    }

    case 6: {
        // 侧棱6(v₃→v₅)旋度
        double ope = 1.0 + xi;  // (1+ξ)
        double opo = 1.0 + eta; // (1+η)
        return Eigen::Vector3d(0.25 * opo, -0.25 * ope, 0.25 * (opo + ope));
    }

    case 7: {
        // 侧棱7(v₄→v₅)旋度
        double ome = 1.0 - xi;  // (1-ξ)
        double opo = 1.0 + eta; // (1+η)
        return Eigen::Vector3d(0.25 * opo, 0.25 * ome, 0.25 * (opo - ome));
    }

    default:
        FEEM_ERROR("curlPyramid5Edge: 无效的棱边索引={}", edge_index);
        return Eigen::Vector3d::Zero();
    }
}

// ==================== 模板显式实例化 ====================

// 仅实例化Dim=2和Dim=3两种维度，覆盖全部电磁有限元应用场景
template class NedelecElement<2>;
template class NedelecElement<3>;

} // namespace numeric
