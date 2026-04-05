/**
 * @file electrostatic_integrator.hpp
 * @brief 数值计算层 - 静电场/瞬态电场单元积分器定义
 * @details 实现基于Lagrange节点元的静电场和瞬态电场有限元单元矩阵组装，
 *          继承EMElementIntegratorBase抽象基类，支持二维和三维常用单元类型。
 *
 *          物理方程背景：
 *          静电场控制方程（Poisson方程）：
 *              -∇·(ε∇φ) = ρ_e
 *          瞬态电场控制方程（电流连续性方程）：
 *              ∇·(ε∇φ) + σ∂φ/∂t = -ρ_e
 *          其中φ为电势，ε为介电常数，σ为电导率，ρ_e为电荷密度。
 *
 *          有限元离散后的矩阵形式：
 *          静态：K·φ = F
 *          瞬态：K·φ + C·∂φ/∂t = F   （或 K·φ + M·∂φ/∂t + C·φ = F）
 *
 * 支持的单元类型：
 * - 二维: TRI3（3节点三角）、QUAD4（4节点四边形）
 * - 三维: TET4（4节点四面体）、HEX8（8节点六面体）、PRISM6（6节点棱柱）、
 *         PRISM15（15节点棱柱）、PYRAMID5（5节点金字塔）、PYRAMID13（13节点金字塔）
 *
 * 核心矩阵公式：
 * - 刚度矩阵 K_e = ∫_Ω (∇N)^T · ε · (∇N) dΩ    （介电项，静态/瞬态共用）
 * - 阻尼矩阵 C_e = ∫_Ω N^T · σ · N dΩ           （电导耗散项，仅瞬态）
 * - 质量矩阵 M_e = ∫_Ω N^T · ε · N dΩ            （介电惯性项，瞬态可选）
 * - 源项向量 F_e = ∫_Ω N · ρ_e dΩ                （电荷激励项，当前为零源项）
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include "em_element_integrator_base.hpp"
#include "shape_function_factory.hpp"
#include <memory>
#include <string>

namespace numeric {

/**
 * @class ElectrostaticIntegrator
 * @brief 静电场/瞬态电场Lagrange节点元积分器
 * @details 继承EMElementIntegratorBase基类，实现标量电势场的单元级有限元矩阵计算。
 *          使用Lagrange形函数进行节点值插值，通过高斯数值积分完成参考域到物理域的矩阵组装。
 *
 *          设计特点：
 *          - 采用策略模式，通过ElementType枚举在构造时确定单元几何类型
 *          - 内部持有std::unique_ptr<ShapeFunctionBase>管理形函数对象生命周期
 *          - 支持非线性介电常数ε(x)和电导率σ(x)的空间变化
 *          - 瞬态模式下自动计算阻尼矩阵C和质量矩阵M
 *
 *          典型使用流程：
 *          1. 构造积分器对象，指定单元类型
 *          2. 调用setMaterialProperties()设置材料参数（ε, σ等）
 *          3. 如需瞬态分析，调用setTransientMode(true)
 *          4. 传入节点坐标，调用computeAllMatrices()或单独的computeXxxMatrix()
 *
 *          性能特征（单次computeAllMatrices调用）：
 *          - 时间复杂度: O(n_gp × n_dof²)，n_gp为高斯积分点数，n_dof为单元自由度数
 *          - 空间复杂度: O(n_dof²)，用于存储单元刚度/阻尼/质量矩阵
 *          - 对典型线性单元（TRI3: 3dof, QUAD4: 4dof, TET4: 4dof, HEX8: 8dof），
 *            单次计算耗时在微秒量级
 *
 * @code
 *   // 创建三维四面体单元的静电场积分器
 *   auto integrator = std::make_unique<ElectrostaticIntegrator>(ElementType::TET4);
 *
 *   // 设置材料参数（如硅: ε_r=11.7, σ≈0）
 *   MaterialProperties silicon(8.854e-12 * 11.7, 1.257e-6, 0.0);
 *   integrator->setMaterialProperties(silicon);
 *
 *   // 计算单元矩阵
 *   Eigen::MatrixXd coords(3, 4);  // 3D坐标 × 4节点
 *   coords << 0, 1, 0, 0.5,
 *             0, 0, 1, 0.5,
 *             0, 0, 0, 1.0;
 *   ElementMatrices mats = integrator->computeAllMatrices(coords);
 *   // mats.K 为 4×4 刚度矩阵
 *   // mats.F 为 4×1 源项向量（零向量）
 * @endcode
 *
 * @see EMElementIntegratorBase 基类接口定义
 * @see LagrangeElement 底层Lagrange形函数实现
 * @see GaussQuadrature 高斯积分点生成
 * @see MaterialProperties 材料参数结构体
 */
class ElectrostaticIntegrator : public EMElementIntegratorBase {
public:
    /**
     * @brief 构造函数，根据单元类型初始化形函数对象
     * @param element_type 有限元单元类型枚举值，
     *        支持2D: TRI3, QUAD4；3D: TET4, HEX8, PRISM6, PRISM15, PYRAMID5, PYRAMID13
     *
     * @details 内部执行以下初始化步骤：
     *          1. 存储单元类型枚举到element_type_
     *          2. 通过ShapeFunctionFactory创建对应类型的Lagrange形函数对象
     *          3. 根据单元类型自动选择默认高斯积分阶数gauss_order_
     *          4. 将is_transient_初始化为false（默认静态模式）
     *
     *          各单元类型的默认积分阶数：
     *          - TRI3: order=1（质心1点积分，精确到常数）
     *          - QUAD4: order=4（2×2 Gauss-Legendre积分）
     *          - TET4: order=1（质心1点积分）
     *          - HEX8: order=8（2×2×2 Gauss-Legendre积分）
     *          - PRISM6: order=6（三角形3点×ζ方向2点）
     *          - PRISM15: order=6（同PRISM6，二次单元需更高阶但此处用默认）
     *          - PYRAMID5: order=5（5点近似积分）
     *          - PYRAMID13: order=5（同PYRAMID5）
     *
     * @warning 若传入不支持的单元类型，输出ERROR日志并使用TRI3作为回退默认值
     */
    explicit ElectrostaticIntegrator(ElementType element_type);

    /**
     * @brief 虚析构函数，确保通过基类指针删除时正确释放形函数对象
     */
    ~ElectrostaticIntegrator() override = default;

    // ========== 纯虚方法实现（来自EMElementIntegratorBase）==========

    /**
     * @brief 计算单元的全部矩阵和源项向量
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count，
     *                    每列对应一个节点的物理坐标 (x, y[, z])^T
     * @return ElementMatrices 包含K、M、C、F的结构体
     *
     * @details 一次性完成所有单元矩阵的高斯积分计算。
     *          内部依次调用：
     *          1. computeStiffnessMatrix() → 填充 result.K
     *          2. computeMassMatrix()      → 填充 result.M（瞬态模式有效）
     *          3. computeDampingMatrix()   → 填充 result.C（瞬态模式有效）
     *          4. computeSourceVector()   → 填充 result.F
     *
     *          数学公式汇总：
     *          - K_ij = ∫_Ω ε(∇N_i)·(∇N_j) dV = Σ_q w_q|detJ_q| ε_q (∇N_i)_q·(∇N_j)_q
     *          - M_ij = ∫_Ω ε N_i N_j dV       = Σ_q w_q|detJ_q| ε_q N_iq N_jq
     *          - C_ij = ∫_Ω σ N_i N_j dV        = Σ_q w_q|detJ_q| σ_q N_iq N_jq
     *          - F_i  = ∫_Ω N_i ρ_e dV          = Σ_q w_q|detJ_q| ρ_eq N_iq = 0（零源项）
     *
     * @note 推荐优先使用此方法以避免重复计算公共中间量（如雅可比矩阵）
     */
    ElementMatrices computeAllMatrices(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元刚度矩阵（介电项）
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 刚度矩阵 K_e，维度为 n_nodes × n_nodes
     *
     * @details 刚度矩阵对应静电场Poisson方程中的∇·(ε∇φ)算子离散化结果。
     *
     *          数值积分公式：
     *          K_e = ∫_Ωe (∇N)^T · ε · (∇N) dΩ
     *             = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ε_eff(ξ_q)
     *               · [∇N(ξ_q)]^T · ∇N(ξ_q)
     *
     *          其中：
     *          - ∇N(ξ_q)：物理域形函数梯度矩阵（node_count × dim），通过
     *                     shape_func_->calcPhysicalGradN(xi, node_coords)获取
     *          - ε_eff(ξ_q)：积分点处的有效介电常数，支持非线性材料，
     *                     通过getEffectiveEpsilon(position)获取
     *          - detJ(ξ_q)：雅可比行列式，用于参考域→物理域体积元变换
     *          - w_q：第q个高斯积分点的权重
     *
     *          矩阵性质：
     *          - 对称正定（或半正定，当ε>0时严格正定）
     *          - 稀疏性取决于单元类型（带状结构）
     *          - 条件数与网格质量相关（畸形单元导致条件数恶化）
     *
     * @note 此方法是静电场求解的核心，其精度直接影响最终解的准确性
     */
    Eigen::MatrixXd computeStiffnessMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元质量矩阵（介电惯性项）
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 质量矩阵 M_e，维度为 n_nodes × n_nodes
     *
     * @details 质量矩阵对应瞬态电场方程中与时间二阶导数相关的介电惯性项。
     *          在标准静电场/准静态电场分析中，此项通常不出现；
     *          仅在高频电磁波传播（含位移电流 ∂D/∂t 的完整Maxwell方程组）
     *          或特定瞬态格式中需要。
     *
     *          数值积分公式：
     *          M_e = ∫_Ωe N^T · ε · N dΩ
     *             = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ε_eff(ξ_q)
     *               · N(ξ_q)^T · N(ξ_q)
     *
     *          其中：
     *          - N(ξ_q)：形函数值向量（node_count × 1），通过
     *                   shape_func_->evalN(xi)获取
     *          - ε_eff(ξ_q)：有效介电常数（同刚度矩阵）
     *
     *          矩阵性质：
     *          - 对称正定（当ε>0时）
     *          - 与刚度矩阵K相同的稀疏模式
     *          - 在集中质量法近似下可替换为对角矩阵
     *
     * @note 仅在瞬态模式下有完整物理意义；静态模式下返回的一致质量矩阵仍可用于
     *       特殊用途（如质量 lumping 技术）
     */
    Eigen::MatrixXd computeMassMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元阻尼矩阵（电导耗散项）
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 阻尼矩阵 C_e，维度为 n_nodes × n_nodes
     *
     * @details 阻尼矩阵对应瞬态电场方程中的电导率耗散项（欧姆损耗）。
     *          来源于电流连续性方程中的σ∂φ/∂t项（或更一般地，J=σE项）。
     *          对于理想绝缘体（σ=0），此矩阵为零矩阵。
     *
     *          数值积分公式：
     *          C_e = ∫_Ωe N^T · σ · N dΩ
     *             = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · σ_eff(ξ_q)
     *               · N(ξ_q)^T · N(ξ_q)
     *
     *          其中：
     *          - N(ξ_q)：形函数值向量（node_count × 1），通过
     *                   shape_func_->evalN(xi)获取
     *          - σ_eff(ξ_q)：积分点处的有效电导率，支持非线性材料，
     *                   通过getEffectiveSigma(position)获取
     *
     *          物理意义：
     *          - 描述导体中的焦耳热损耗功率 P = ∫ σ|E|² dV 的离散化
     *          - 在电路类比中相当于电阻矩阵 R = C⁻¹
     *          - 对于良导体（如铜 σ≈5.8×10⁷ S/m），C矩阵元素值很大
     *
     *          矩阵性质：
     *          - 对称半正定（当σ≥0时）
     *          - 当σ=0时退化为零矩阵（理想绝缘体无耗散）
     *
     * @warning 仅在瞬态模式（is_transient_=true）下此矩阵参与时域求解
     */
    Eigen::MatrixXd computeDampingMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 单独计算单元源项向量（电荷密度激励项）
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::VectorXd 源项向量 F_e，维度为 n_nodes × 1
     *
     * @details 源项向量对应Poisson方程右端项的电荷密度分布。
     *          当前实现采用零源项假设（ρ_e ≡ 0），适用于：
     *          - 拉普拉斯方程求解（无自由电荷区域）
     *          - 边界条件驱动的电场问题（电极施加电压）
     *          - 电荷中性区域的电场分布计算
     *
     *          一般形式的数值积分公式：
     *          F_e = ∫_Ωe N · ρ_e dΩ
     *             = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ρ_e(ξ_q) · N(ξ_q)
     *
     *          当前实现（ρ_e = 0）：
     *          F_e = 0 向量（所有分量均为零）
     *
     * @note 后续版本可扩展支持空间变化的电荷密度分布 ρ_e(x,y,z)，
     *       包括点电荷、线电荷、面电荷等特殊源项的处理
     */
    Eigen::VectorXd computeSourceVector(
        const Eigen::MatrixXd& node_coords
    ) const override;

    // ========== 查询方法 ==========

    /**
     * @brief 获取当前积分器的单元类型
     * @return ElementType 当前使用的单元类型枚举值
     */
    ElementType getElementType() const;

    /**
     * @brief 获取当前单元的节点数量
     * @return int 节点总数（如TRI3返回3，HEX8返回8）
     *
     * @details 直接委托给内部形函数对象的getNodeCount()方法
     */
    int getNodeCount() const;

    /**
     * @brief 获取当前单元的空间维度
     * @return int 空间维度（2或3）
     *
     * @details 直接委托给内部形函数对象的getDim()方法
     */
    int getDimension() const;

    /**
     * @brief 获取当前使用的高斯积分阶数
     * @return int 积分阶数（各单元类型的默认值见构造函数文档）
     */
    int getGaussOrder() const;

    /**
     * @brief 设置高斯积分阶数
     * @param order 新的积分阶数，必须为正整数
     *
     * @details 更改积分阶数会影响矩阵计算的精度和计算成本：
     *          - 较低阶数：计算快但可能精度不足（尤其对畸形单元）
     *          - 较高阶数：精度高但计算开销增大
     *
     * @note 调用前建议确认目标单元类型支持指定的order值
     * @warning 若order <= 0，输出WARN日志并保持原值不变
     */
    void setGaussOrder(int order);

private:
    /**
     * @brief 根据单元类型枚举值获取对应的字符串标识符
     * @param type 单元类型枚举值
     * @return std::string 对应的类型名称字符串（如"TRI3"、"HEX8"等），
     *         不支持的类型返回空字符串
     *
     * @details 用于将ElementType枚举转换为ShapeFunctionFactory所需的字符串参数。
     *          映射表覆盖全部支持的Lagrange单元类型。
     */
    static std::string elementTypeToString(ElementType type);

    /**
     * @brief 根据单元类型确定默认的高斯积分阶数
     * @param type 单元类型枚举值
     * @return int 该单元类型推荐的默认积分阶数
     *
     * @details 各单元类型的默认积分方案：
     *          - TRI3:   1 （质心单点积分）
     *          - QUAD4:  4 （2×2 Gauss-Legendre）
     *          - TET4:   1 （质心单点积分）
     *          - HEX8:   8 （2×2×2 Gauss-Legendre）
     *          - PRISM6: 6 （三角形3点×ζ方向2点）
     *          - PRISM15:6 （同PRISM6，二次单元可用更高阶）
     *          - PYRAMID5: 5 （5点近似积分）
     *          - PYRAMID13:5 （同PYRAMID5）
     */
    static int getDefaultGaussOrder(ElementType type);

    ElementType element_type_;                      ///< 单元类型枚举值
    int gauss_order_;                               ///< 高斯积分阶数
    std::unique_ptr<ShapeFunctionBase> shape_func_; ///< 形函数对象指针（RAII管理）
};

} // namespace numeric
