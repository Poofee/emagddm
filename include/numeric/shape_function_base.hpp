/**
 * @file shape_function_base.hpp
 * @brief 数值计算层 - 有限元形函数抽象基类定义
 * @details 定义形函数计算的统一接口，支持Lagrange节点元和Nedelec棱边矢量元，
 *          提供雅可比矩阵计算、物理坐标梯度映射等通用实现。
 *          所有具体单元类型（线、三角、四边形、四面体、六面体等）均需继承此基类。
 * @author Poofee
 * @date 2026-04-04
 * @version 1.0
 */

#pragma once

#include <Eigen/Dense>
#include "logger_factory.hpp"

namespace numeric {

// ==================== 单元类型枚举 ====================

/**
 * @enum ElementType
 * @brief 有限元单元类型枚举
 * @details 覆盖全部常用单元类型，分为两大类：
 *          - Lagrange节点元：标量插值，用于电势、磁标势等标量场
 *          - Nedelec棱边矢量元：矢量插值（切向连续），用于磁场矢量位A等矢量场
 */
enum class ElementType {
    // ---------- 一维Lagrange单元 ----------
    LINE2 = 0,       ///< 二节点线单元（线性）
    LINE3,           ///< 三节点线单元（二次）

    // ---------- 二维Lagrange单元 ----------
    TRI3,            ///< 三节点三角形单元（线性）
    TRI6,            ///< 六节点三角形单元（二次）
    QUAD4,           ///< 四节点四边形单元（双线性）
    QUAD8,           ///< 八节点四边形单元（Serendipity二次）
    QUAD9,           ///< 九节点四边形单元（Lagrange二次）

    // ---------- 三维Lagrange单元 ----------
    TET4,            ///< 四节点四面体单元（线性）
    TET10,           ///< 十节点四面体单元（二次）
    HEX8,            ///< 八节点六面体单元（三线性）
    HEX20,           ///< 二十节点六面体单元（Serendipity二次）
    HEX27,           ///< 二十七节点六面体单元（Lagrange二次）
    PRISM6,          ///< 六节点棱柱单元（线性）
    PRISM15,         ///< 十五节点棱柱单元（二次）
    PYRAMID5,        ///< 五节点棱锥单元（线性）
    PYRAMID13,       ///< 十三节点棱锥单元（二次）

    // ---------- Nedelec棱边矢量元 ----------
    TRI3_EDGE,       ///< 三角形一阶Nedelec棱边元（3条棱边自由度）
    QUAD4_EDGE,      ///< 四边形一阶Nedelec棱边元（4条棱边自由度）
    TET4_EDGE,       ///< 四面体一阶Nedelec棱边元（6条棱边自由度）
    HEX8_EDGE,       ///< 六面体一阶Nedelec棱边元（12条棱边自由度）
    PRISM6_EDGE,     ///< 棱柱一阶Nedelec棱边元（9条棱边自由度）
    PYRAMID5_EDGE    ///< 棱锥一阶Nedelec棱边元（8条棱边自由度）
};

// ==================== 数据结构体 ====================

/**
 * @struct LocalPoint
 * @brief 局部坐标点（参考域坐标）
 * @details 使用Eigen::Vector3d统一存储1D/2D/3D局部坐标，
 *          通过dim成员标识实际维度，未使用的分量默认为零。
 *          参考域定义：
 *          - 1D: ξ ∈ [-1, 1]
 *          - 2D: (ξ, η) ∈ [-1, 1]²
 *          - 3D: (ξ, η, ζ) ∈ [-1, 1]³
 */
struct LocalPoint {
    Eigen::Vector3d coords;  ///< 局部坐标向量，统一3维存储
    int dim;                 ///< 实际维度（1/2/3）

    /**
     * @brief 默认构造函数，初始化为零向量，维度为3
     */
    LocalPoint()
        : coords(Eigen::Vector3d::Zero())
        , dim(3)
    {
    }

    /**
     * @brief 一维局部坐标构造函数
     * @param xi 一维坐标ξ ∈ [-1, 1]
     */
    explicit LocalPoint(double xi)
        : coords(xi, 0.0, 0.0)
        , dim(1)
    {
    }

    /**
     * @brief 二维局部坐标构造函数
     * @param xi 第一维坐标 ξ ∈ [-1, 1]
     * @param eta 第二维坐标 η ∈ [-1, 1]
     */
    LocalPoint(double xi, double eta)
        : coords(xi, eta, 0.0)
        , dim(2)
    {
    }

    /**
     * @brief 三维局部坐标构造函数
     * @param xi 第一维坐标 ξ ∈ [-1, 1]
     * @param eta 第二维坐标 η ∈ [-1, 1]
     * @param zeta 第三维坐标 ζ ∈ [-1, 1]
     */
    LocalPoint(double xi, double eta, double zeta)
        : coords(xi, eta, zeta)
        , dim(3)
    {
    }
};

/**
 * @struct JacobianResult
 * @brief 雅可比矩阵计算结果
 * @details 存储从参考域到物理域的坐标变换雅可比信息，
 *          用于形函数梯度的物理坐标映射和积分体积元的计算。
 *
 *          雅可比矩阵定义为：J = ∂x/∂ξ = dN/dξ · X_node
 *          其中X_node为节点物理坐标矩阵（dim × node_count）
 */
struct JacobianResult {
    Eigen::MatrixXd jacobian;   ///< 雅可比矩阵 J（dim × dim 或 dim × ref_dim）
    double det_j;               ///< 雅可比行列式值 det(J)，用于积分体积元 dV = |detJ| dξ
    Eigen::MatrixXd inv_jacobian;///< 雅可比逆矩阵 J⁻¹，用于梯度从参考域到物理域的映射

    /**
     * @brief 默认构造函数，初始化空矩阵和零值
     */
    JacobianResult()
        : det_j(0.0)
    {
    }
};

/**
 * @struct EdgeCurlResult
 * @brief Nedelec棱边元旋度计算结果
 * @details 存储Nedelec棱边矢量形函数的旋度信息，
 *          用于电磁场中curl-curl算子的矩阵组装（如∇×∇×A）。
 *          包含参考域旋度和物理域旋度两种表示。
 */
struct EdgeCurlResult {
    Eigen::Vector3d curl_reference;   ///< 参考域旋度向量 (∇_ξ × N_edge)
    Eigen::Vector3d curl_physical;    ///< 物理域旋度向量 (∇_x × N_edge)

    /**
     * @brief 默认构造函数，初始化零向量
     */
    EdgeCurlResult()
        : curl_reference(Eigen::Vector3d::Zero())
        , curl_physical(Eigen::Vector3d::Zero())
    {
    }
};

// ==================== 形函数抽象基类 ====================

/**
 * @class ShapeFunctionBase
 * @brief 有限元形函数抽象基类
 * @details 定义所有单元形函数的统一接口，是有限元数值积分的核心组件。
 *          支持两类形函数：
 *          - Lagrange节点元：提供evalN（形函数值）和evalGradN（形函数梯度）接口
 *          - Nedelec棱边矢量元：额外提供evalEdgeFunction（棱边形函数）、
 *            evalCurlEdge（旋度）和calcPhysicalEdgeFunction（物理域映射）接口
 *
 *          子类需实现的纯虚方法：
 *          - getNodeType(): 返回单元类型枚举
 *          - getNodeCount(): 返回单元节点/棱边数
 *          - getDim(): 返回空间维度
 *          - evalN(): 计算参考域形函数值向量
 *          - evalGradN(): 计算参考域形函数梯度矩阵
 *
 *          基类已提供默认实现的方法：
 *          - calcJacobian(): 基于evalGradN计算雅可比矩阵
 *          - calcPhysicalGradN(): 将参考域梯度映射到物理坐标
 *
 * @note 此类位于numeric命名空间下，属于数值计算层的核心基础模块
 * @see ElementType 单元类型枚举
 * @see LocalPoint 局部坐标点结构
 * @see JacobianResult 雅可比结果结构
 */
class ShapeFunctionBase {
public:
    /**
     * @brief 虚析构函数，确保派生类正确释放资源
     */
    virtual ~ShapeFunctionBase() = default;

    // ========== 纯虚方法：子类必须实现 ==========

    /**
     * @brief 获取单元节点类型
     * @return ElementType 当前单元的类型枚举值
     */
    virtual ElementType getNodeType() const = 0;

    /**
     * @brief 获取单元节点（或棱边）总数
     * @return int 节点/棱边数量，如TRI3返回3，HEX8返回8，TET4_EDGE返回6
     */
    virtual int getNodeCount() const = 0;

    /**
     * @brief 获取单元空间维度
     * @return int 空间维度（1/2/3），如LINE2返回1，QUAD4返回2，TET4返回3
     */
    virtual int getDim() const = 0;

    /**
     * @brief 计算参考域形函数值向量
     * @param xi 参考域局部坐标点
     * @return Eigen::VectorXd 形函数值向量，维度为node_count × 1，
     *         第i个元素为第i个节点在xi处的形函数值N_i(ξ)
     *
     * @note 形函数满足单位分解性：Σ N_i(ξ) = 1
     */
    virtual Eigen::VectorXd evalN(const LocalPoint& xi) const = 0;

    /**
     * @brief 计算参考域形函数梯度矩阵
     * @param xi 参考域局部坐标点
     * @return Eigen::MatrixXd 形函数梯度矩阵，维度为 node_count × dim_ref，
     *         第i行为第i个节点的梯度行向量 [∂N_i/∂ξ₁, ∂N_i/∂ξ₂, ...]
     *
     * @note 对于三维单元，每行为[∂N_i/∂ξ, ∂N_i/∂η, ∂N_i/∂ζ]；
     *       对于二维单元，每行为[∂N_i/∂ξ, ∂N_i/∂η]
     */
    virtual Eigen::MatrixXd evalGradN(const LocalPoint& xi) const = 0;

    // ========== 带默认实现的虚方法 ==========

    /**
     * @brief 计算雅可比矩阵及其相关量
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵，维度为 dim × node_count，
     *                    每列为一个节点的物理坐标向量
     * @return JacobianResult 雅可比计算结果，包含J、det(J)、J⁻¹
     *
     * @details 计算公式：J = ∂x/∂ξ = (∂N/∂ξ) · X_node
     *          其中∂N/∂ξ为evalGradN()返回的梯度矩阵(node_count × ref_dim)，
     *          X_node为节点坐标矩阵(dim × node_count)，
     *          最终J的维度为(dim × ref_dim)或方阵(dim × dim)
     *
     * @note 默认实现基于evalGradN()，子类可覆写以优化特殊单元的计算
     * @warning 当detJ接近零时（畸形单元），数值精度会显著下降
     */
    virtual JacobianResult calcJacobian(
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords
    ) const;

    /**
     * @brief 计算物理坐标系下的形函数梯度
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 物理域梯度矩阵，维度为 node_count × dim，
     *         第i行为第i个节点的物理梯度 [∂N_i/∂x, ∂N_i/∂y, ∂N_i/∂z]
     *
     * @details 计算公式：∇_x N = J⁻¹ · ∇_ξ N
     *          先调用calcJacobian获取J⁻¹，再与参考域梯度相乘完成映射
     *
     * @note 默认实现基于calcJacobian()和evalGradN()，子类可覆写以优化性能
     */
    virtual Eigen::MatrixXd calcPhysicalGradN(
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords
    ) const;

    // ========== Nedelec棱边元扩展接口（默认返回零/空）==========

    /**
     * @brief 计算Nedelec棱边形函数在参考域的矢量值
     * @param edge_index 棱边索引（从0开始编号），对应单元的第edge_index条棱边
     * @param xi 参考域局部坐标点
     * @return Eigen::Vector3d 棱边形函数在xi处的矢量值（参考域表示）
     *
     * @details Nedelec棱边元的形函数沿棱边切向具有单位积分特性：
     *          ∫_edge (N_edge · t) ds = δ_{edge_index, edge_j}
     *          其中t为棱边切向单位向量
     *
     * @note 默认实现返回零向量，仅Nedelec单元子类需要覆写此方法
     * @warning 对非Nedelec单元调用此方法将得到无意义的零结果
     */
    virtual Eigen::Vector3d evalEdgeFunction(
        int edge_index,
        const LocalPoint& xi
    ) const;

    /**
     * @brief 计算Nedelec棱边形函数在参考域的旋度
     * @param edge_index 棱边索引（从0开始编号）
     * @param xi 参考域局部坐标点
     * @return Eigen::Vector3d 棱边形函数旋度向量（参考域表示）
     *
     * @details 旋度用于组装电磁场的curl-curl项（如∇×∇×A），
     *          是涡流场和波动场有限元方程的核心组成部分
     *
     * @note 默认实现返回零向量，仅Nedelec单元子类需要覆写此方法
     */
    virtual Eigen::Vector3d evalCurlEdge(
        int edge_index,
        const LocalPoint& xi
    ) const;

    /**
     * @brief 计算Nedelec棱边形函数在物理域的完整结果（含旋度映射）
     * @param edge_index 棱边索引（从0开始编号）
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵，维度为 dim × node_count
     * @return EdgeCurlResult 包含参考域和物理域的旋度结果
     *
     * @details 物理域旋度的映射关系：∇_x × N = (1/detJ) · J · (∇_ξ × N)
     *          其中J为雅可比矩阵，detJ为其行列式
     *          此方法同时返回参考域旋度和映射后的物理域旋度
     *
     * @note 默认实现基于evalCurlEdge()和calcJacobian()，子类可覆写优化
     */
    virtual EdgeCurlResult calcPhysicalEdgeFunction(
        int edge_index,
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords
    ) const;
};

// ==================== 默认方法实现 ====================

inline JacobianResult ShapeFunctionBase::calcJacobian(
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords
) const {
    JacobianResult result;

    // 获取参考域形函数梯度矩阵（node_count × dim_ref）
    Eigen::MatrixXd grad_n = evalGradN(xi);

    // 计算雅可比矩阵 J = X_node · (∂N/∂ξ)
    // node_coords维度: dim × node_count
    // grad_n维度: node_count × ref_dim
    // 结果jacobian维度: dim × ref_dim（对于等参元为方阵 dim × dim）
    result.jacobian = node_coords * grad_n;

    // 雅可比矩阵为方阵时计算行列式和逆矩阵
    // 对于等参元，参考域维度等于物理域维度，J为方阵
    if (result.jacobian.rows() == result.jacobian.cols()) {
        result.det_j = result.jacobian.determinant();
        result.inv_jacobian = result.jacobian.inverse();
    } else {
        // 非方阵情况（如退化单元），使用伪逆近似
        result.det_j = 0.0;
        FEEM_WARN("雅可比矩阵非方阵（{}×{}），使用Moore-Penrose伪逆",
                  result.jacobian.rows(), result.jacobian.cols());
        result.inv_jacobian = result.jacobian.completeOrthogonalDecomposition().pseudoInverse();
    }

    return result;
}

inline Eigen::MatrixXd ShapeFunctionBase::calcPhysicalGradN(
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords
) const {
    // 计算雅可比矩阵获取逆映射矩阵 J⁻¹
    JacobianResult jaco_result = calcJacobian(xi, node_coords);

    // 获取参考域形函数梯度矩阵 ∇_ξ N（node_count × dim_ref）
    Eigen::MatrixXd grad_n_ref = evalGradN(xi);

    // 物理域梯度映射: ∇_x N = (∇_ξ N) · J⁻¹
    // grad_n_ref维度: node_count × dim_ref
    // inv_jacobian维度: ref_dim × dim（对于等参元 ref_dim = dim）
    // 结果维度: node_count × dim（物理空间维度）
    return grad_n_ref * jaco_result.inv_jacobian;
}

inline Eigen::Vector3d ShapeFunctionBase::evalEdgeFunction(
    int edge_index,
    const LocalPoint& xi
) const {
    (void)xi;
    // 默认返回零向量，仅Nedelec棱边元子类需覆写此方法
    FEEM_DEBUG("evalEdgeFunction: 非Nedelec单元调用，返回零向量, edge_index={}", edge_index);
    return Eigen::Vector3d::Zero();
}

inline Eigen::Vector3d ShapeFunctionBase::evalCurlEdge(
    int edge_index,
    const LocalPoint& xi
) const {
    (void)xi;
    // 默认返回零向量，仅Nedelec棱边元子类需覆写此方法
    FEEM_DEBUG("evalCurlEdge: 非Nedelec单元调用，返回零向量, edge_index={}", edge_index);
    return Eigen::Vector3d::Zero();
}

inline EdgeCurlResult ShapeFunctionBase::calcPhysicalEdgeFunction(
    int edge_index,
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords
) const {
    EdgeCurlResult result;

    // 获取参考域旋度（由子类Nedelec单元提供，默认为零）
    result.curl_reference = evalCurlEdge(edge_index, xi);

    // 计算雅可比矩阵用于坐标映射
    JacobianResult jaco_result = calcJacobian(xi, node_coords);

    // 物理域旋度映射公式: ∇_x × N = (1/detJ) · J · (∇_ξ × N)
    // 此公式适用于三维等参变换的Nedelec元旋度映射
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

} // namespace numeric
