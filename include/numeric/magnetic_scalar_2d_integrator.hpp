/**
 * @file magnetic_scalar_2d_integrator.hpp
 * @brief 数值计算层 - 二维静磁场/瞬态磁场标量位积分器定义
 * @details 实现基于标量磁势 φ_m 的二维有限元单元矩阵组装，
 *          支持静磁场（Poisson方程）和瞬态磁场（磁扩散方程）两种分析模式。
 *          仅支持2D Lagrange单元（TRI3三角形单元、QUAD4四边形单元）。
 *
 *          物理方程：
 *          - 静态磁场: ∇·(μ∇φ_m) = -ρ_m   （Poisson型方程）
 *          - 瞬态磁场: ∇·(μ∇φ_m) + σ∂φ_m/∂t = -ρ_m  （磁扩散/涡流方程）
 *
 *          弱形式离散：
 *          - 刚度矩阵 K_m = ∫_Ω (∇N)^T · μ · (∇N) dΩ  （磁阻项）
 *          - 质量矩阵 M_m = ∫_Ω N^T · σ · N dΩ            （涡流耗散项）
 *          - 源项向量 F_m = ∫_Ω N · ρ_m dΩ               （等效磁荷源）
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include <memory>
#include "em_element_integrator_base.hpp"
#include "lagrange_element.hpp"
#include "gauss_quadrature.hpp"
#include "logger_factory.hpp"

namespace numeric {

/**
 * @class MagneticScalar2DIntegrator
 * @brief 二维静磁场/瞬态磁场标量磁势积分器
 * @details 继承EMElementIntegratorBase抽象基类，实现二维标量磁势φ_m的单元级
 *          有限元矩阵组装。采用Lagrange节点形函数进行标量场插值，
 *          通过高斯数值积分完成刚度矩阵、质量矩阵和源项向量的计算。
 *
 *          适用场景：
 *          - 二维静磁场分析：电机气隙磁场、永磁体磁场分布等
 *          - 二维瞬态涡流分析：变压器涡流损耗、电机起动过程等
 *          - 磁屏蔽效果评估：非线性铁磁材料B-H曲线建模
 *
 *          支持的单元类型（仅2D）：
 *          - TRI3: 3节点线性三角形单元，使用1点或3点高斯积分
 *          - QUAD4: 4节点双线性四边形单元，使用2×2=4点高斯积分
 *
 *          材料模型支持：
 *          - 线性均匀材料：恒定磁导率μ和电导率σ
 *          - 非线性磁导率：通过覆写getEffectiveMu()实现B-H曲线插值
 *
 *          数值积分策略：
 *          - TRI3: order=1（质心1点，精确到常数）或order=3（边中3点，精确到二次）
 *          - QUAD4: order=4（2×2 Gauss-Legendre，精确到三次多项式）
 *
 * @note 与电场静电积分器的区别：
 *       1. 物理场不同：ε → μ（介电常数→磁导率）
 *       2. 维度限制：仅支持2D单元（TRI3/QUAD4）
 *       3. 矩阵命名：瞬态用质量矩阵M_m（非阻尼C_e），对应不同的物理机制
 *       4. 默认参数：使用真空磁导率μ₀=4π×10⁻⁷ H/m
 *
 * @see EMElementIntegratorBase 基类接口
 * @see LagrangeElement<2> 二维Lagrange形函数
 * @see GaussQuadrature 高斯积分类
 */
class MagneticScalar2DIntegrator : public EMElementIntegratorBase {
public:
    /**
     * @brief 构造函数，初始化二维标量磁势积分器
     * @param element_type 2D单元类型枚举值，仅支持TRI3或QUAD4
     *
     * @details 内部执行以下初始化操作：
     *          1. 校验单元类型的合法性（必须为2D Lagrange单元）
     *          2. 根据单元类型实例化对应的LagrangeElement<2>形函数对象
     *          3. 设置默认材料参数（真空磁导率μ₀、零电导率）
     *          4. 默认关闭瞬态模式（is_transient_=false）
     *
     * @exception 若传入非2D单元类型或非法枚举值，输出ERROR日志并抛出std::invalid_argument异常
     *
     * @code
     * // 创建三角形单元积分器
     * auto tri_integ = std::make_unique<MagneticScalar2DIntegrator>(ElementType::TRI3);
     *
     * // 创建四边形单元积分器
     * auto quad_integ = std::make_unique<MagneticScalar2DIntegrator>(ElementType::QUAD4);
     * @endcode
     */
    explicit MagneticScalar2DIntegrator(ElementType element_type);

    /**
     * @brief 虚析构函数，释放形函数智能指针资源
     */
    ~MagneticScalar2DIntegrator() override = default;

    // ========== 纯虚方法实现（来自EMElementIntegratorBase）==========

    /**
     * @brief 计算单元的全部矩阵和源项向量（一次性完整计算）
     * @param node_coords 节点坐标矩阵，维度为 2 × node_count，
     *                    每列对应一个节点的物理坐标 (x, y)^T
     * @return ElementMatrices 包含K_m（刚度）、M_m（质量）、C（阻尼，本模块为零）、F_m（源项）的结构体
     *
     * @details 核心入口方法，内部调用流程：
     *          1. 计算刚度矩阵 K_m = computeStiffnessMatrix(node_coords)
     *          2. 若瞬态模式开启，计算质量矩阵 M_m = computeMassMatrix(node_coords)
     *          3. 计算源项向量 F_m = computeSourceVector(node_coords)
     *          4. 组装返回ElementMatrices结构体
     *
     *          数学公式：
     *          K_m[i][j] = Σ_q w_q · |detJ_q| · μ_eff · (∇N_i·∇N_j)
     *          M_m[i][j] = Σ_q w_q · |detJ_q| · σ_eff · (N_i·N_j)   （仅瞬态模式）
     *          F_m[i] = 0                                           （当前零源项）
     *
     * @note 调用前需确保node_coords维度正确（2行×node_count列）
     * @warning 返回的C阻尼矩阵在本模块中始终为零矩阵（标量磁势无阻尼项）
     */
    ElementMatrices computeAllMatrices(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元刚度矩阵 K_m（磁阻项）
     * @param node_coords 节点坐标矩阵，维度为 2 × node_count
     * @return Eigen::MatrixXd 刚度矩阵 K_m，维度为 n_dof × n_dof
     *
     * @details 刚度矩阵对应磁扩散方程中的空间微分算子项 ∇·(μ∇φ_m)。
     *          采用Gauss-Legendre数值积分在参考域上计算：
     *
     *          K_m = ∫_Ω̂ (∇_x N)^T · μ · (∇_x N) · |detJ| dξ̂
     *              ≈ Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · μ_eff · G_q^T · G_q
     *
     *          其中：
     *          - G_q = ∇_x N(ξ_q) 为第q个积分点处的物理域形函数梯度矩阵（n_dof × 2）
     *          - detJ(ξ_q) 为雅可比行列式（2D时等于面积缩放因子）
     *          - μ_eff 为有效磁导率（可通过getEffectiveMu获取非线性值）
     *          - w_q 为第q个积分点的权重
     *
     *          矩阵特性：
     *          - 对称正定（当μ > 0且单元非退化时）
     *          - 条件数与单元形状质量相关（畸形单元条件数恶化）
     *          - TRI3: 3×3对称矩阵；QUAD4: 4×4对称矩阵
     *
     * @note 此方法同时用于静态和瞬态模式（K_m是共用项）
     */
    Eigen::MatrixXd computeStiffnessMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元质量矩阵 M_m（涡流耗散项，仅瞬态模式有意义）
     * @param node_coords 节点坐标矩阵，维度为 2 × node_count
     * @return Eigen::MatrixXd 质量矩阵 M_m，维度为 n_dof × n_dof
     *
     * @details 质量矩阵对应瞬态磁扩散方程中的时间一阶导数项 σ∂φ_m/∂t。
     *          在半离散化后形成 M_m · {dφ/dt} 项，用于时间推进求解。
     *
     *          积分公式：
     *          M_m = ∫_Ω̂ N^T · σ · N · |detJ| dξ̂
     *              ≈ Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · σ_eff · N_q^T · N_q
     *
     *          其中：
     *          - N_q = N(ξ_q) 为第q个积分点处的形函数值向量（n_dof × 1）
     *          - σ_eff 为有效电导率（导体区σ > 0，绝缘体/空气σ ≈ 0）
     *
     *          物理意义：
     *          - 导体区域（如铜绕组、铁芯）：σ >> 0，M_m显著影响瞬态响应
     *          - 绝缘区域（如空气、真空）：σ ≈ 0，M_m接近零矩阵
     *          - 涡流效应完全由M_m矩阵描述（与K_m共同构成磁扩散算子）
     *
     * @note 仅在瞬态模式（isTransient()=true）下此矩阵有物理意义
     * @note 当σ=0（默认值）时，返回零矩阵
     */
    Eigen::MatrixXd computeMassMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元阻尼矩阵 C（标量磁势模型中此项为零）
     * @param node_coords 节点坐标矩阵，维度为 2 × node_count
     * @return Eigen::MatrixXd 阻尼矩阵 C，始终返回零矩阵（n_dof × n_dof）
     *
     * @details 在标量磁势φ_m的公式体系中，不存在独立的零阶阻尼项。
     *          与矢量磁势A法不同（A法的C矩阵来自σA项），标量磁势法的
     *          时间导数项仅通过质量矩阵M_m描述（σ∂φ_m/∂t项）。
     *
     *          返回零矩阵以保持接口一致性，全局组装时可安全忽略此矩阵。
     *
     * @note 此方法的存在仅为满足基类接口契约，实际计算中不使用C矩阵
     */
    Eigen::MatrixXd computeDampingMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元源项向量 F_m（等效磁荷激励项）
     * @param node_coords 节点坐标矩阵，维度为 2 × node_count
     * @return Eigen::VectorXd 源项向量 F_m，当前实现返回零向量（n_dof × 1）
     *
     * @details 源项向量对应等效磁荷密度 ρ_m 的单元投影：
     *          F_m[i] = ∫_Ω N_i · ρ_m dΩ
     *
     *          当前版本假设无等效磁荷源（ρ_m = 0），适用于：
     *          - 永磁体磁场（通过边界条件或面电流模拟，而非体磁荷）
     *          - 励磁线圈磁场（通过矢量位耦合或边界激励施加）
     *          - 纯边值问题（狄利克雷/诺伊曼边界条件驱动）
     *
     *          后续扩展方向：
     *          - 支持永磁体的等效磁荷模型（ρ_m = -∇·M）
     *          - 支持外部磁场源的体投影
     *
     * @note 当前返回零向量，后续可扩展支持非零源项
     */
    Eigen::VectorXd computeSourceVector(
        const Eigen::MatrixXd& node_coords
    ) const override;

    // ========== 查询方法 ==========

    /**
     * @brief 获取当前积分器的单元类型
     * @return ElementType 当前使用的2D单元类型枚举值（TRI3或QUAD4）
     */
    ElementType getElementType() const;

    /**
     * @brief 获取当前单元的节点数量
     * @return int 节点数（TRI3返回3，QUAD4返回4）
     */
    int getNodeCount() const;

private:
    std::unique_ptr<LagrangeElement<2>> shape_func_;  ///< 二维Lagrange形函数智能指针
    ElementType element_type_;                         ///< 当前单元类型枚举值

    /**
     * @brief 根据单元类型确定高斯积分阶数
     * @return int 高斯积分点数量（TRI3返回3，QUAD4返回4）
     *
     * @details 选择保证精度要求的最低积分阶数：
     *          - TRI3（线性三角形）：3点积分，精确到二次多项式
     *          - QUAD4（双线性四边形）：4点（2×2）积分，精确到三次多项式
     */
    int getGaussOrder() const;
};

} // namespace numeric
