/**
 * @file em_element_integrator_base.hpp
 * @brief 数值计算层 - 电磁场单元积分器抽象基类定义
 * @details 定义电磁场有限元分析中单元级矩阵/向量计算的统一接口，
 *          支持静态场、瞬态涡流场和频域波动场的单元矩阵组装。
 *          所有具体积分器（Lagrange元、Nedelec元等）均需继承此基类。
 *
 *          核心数学基础：
 *          - 刚度矩阵 K_e = ∫_Ωe (∇×N_i)·(∇×N_j) dΩ  （curl-curl项，矢量位A）
 *          - 质量矩阵 M_e = ∫_Ωe N_i·N_j dΩ            （质量项，瞬态惯性）
 *          - 阻尼矩阵 C_e = ∫_Ωe σ N_i·N_j dΩ         （电导率耗散项）
 *          - 源项向量 F_e = ∫_Ωe N_i·J_src dΩ          （源电流激励项）
 *
 *          适用场景：
 *          - 低频电磁A-φ法：静磁场、涡流场（含K、C、F矩阵）
 *          - 高频电磁波：时谐场/瞬态波动场（含K、M、C、F全量矩阵）
 *          - 非线性材料：支持ε(x)、μ(x)、σ(x)的空间变化和非线性特性
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include <Eigen/Dense>
#include "shape_function_base.hpp"
#include "gauss_quadrature.hpp"
#include "logger_factory.hpp"

namespace numeric {

// ==================== 数据结构体定义 ====================

/**
 * @struct MaterialProperties
 * @brief 电磁材料参数结构体
 * @details 封装电磁场分析所需的全部材料物理参数，
 *          包含线性参数值和非线性标志位，用于驱动单元矩阵的材料相关计算。
 *
 *          物理量说明：
 *          - epsilon（介电常数ε）：描述介质极化能力，真空值 ε₀ ≈ 8.854×10⁻¹² F/m
 *          - mu（磁导率μ）：描述介质磁化能力，真空值 μ₀ = 4π×10⁻⁷ H/m
 *          - sigma（电导率σ）：描述导电性能，绝缘体≈0 S/m，铜≈5.8×10⁷ S/m
 */
struct MaterialProperties {
    double epsilon;               ///< 介电常数 ε (F/m)，默认值为真空介电常数 ε₀
    double mu;                    ///< 磁导率 μ (H/m)，默认值为真空磁导率 μ₀
    double sigma;                 ///< 电导率 σ (S/m)，绝缘体为0，导体为有限正值
    bool is_nonlinear_epsilon;    ///< 是否非线性介电常数（如铁电材料），false=线性常量
    bool is_nonlinear_mu;         ///< 是否非线性磁导率（如铁磁材料B-H曲线），false=线性常量
    bool is_nonlinear_sigma;      ///< 是否非线性电导率（如温度依赖），false=线性常量

    /**
     * @brief 默认构造函数，初始化为真空材料参数
     * @details 默认值：
     *          - ε = 8.854187817e-12 F/m（真空介电常数）
     *          - μ = 4π × 10⁻⁷ H/m（真空磁导率，约1.2566370614e-6）
     *          - σ = 0 S/m（理想绝缘体）
     *          - 非线性标志均为false（线性材料）
     */
    MaterialProperties()
        : epsilon(8.854187817e-12)
        , mu(1.2566370614359173e-6)
        , sigma(0.0)
        , is_nonlinear_epsilon(false)
        , is_nonlinear_mu(false)
        , is_nonlinear_sigma(false)
    {
    }

    /**
     * @brief 参数化构造函数
     * @param eps 介电常数 ε (F/m)，必须为正数
     * @param m 磁导率 μ (H/m)，必须为正数
     * @param sig 电导率 σ (S/m)，必须为非负数
     * @param nl_eps 是否非线性介电常数，默认false
     * @param nl_mu 是否非线性磁导率，默认false
     * @param nl_sig 是否非线性电导率，默认false
     */
    MaterialProperties(double eps, double m, double sig,
                       bool nl_eps = false, bool nl_mu = false, bool nl_sig = false)
        : epsilon(eps)
        , mu(m)
        , sigma(sig)
        , is_nonlinear_epsilon(nl_eps)
        , is_nonlinear_mu(nl_mu)
        , is_nonlinear_sigma(nl_sig)
    {
    }
};

/**
 * @struct ElementMatrices
 * @brief 单元级有限元矩阵结果集合
 * @details 存储单个单元的所有有限元矩阵和向量，
 *          用于后续的全局刚度矩阵组装（scatter操作）。
 *
 *          矩阵维度约定（以n_dof表示单元自由度数）：
 *          - K: n_dof × n_dof 刚度矩阵（对称正定或半正定）
 *          - M: n_dof × n_dof 质量矩阵（对称正定）
 *          - C: n_dof × n_dof 阻尼矩阵（对称半正定）
 *          - F: n_dof × 1 源项向量
 *
 *          矩阵的物理含义（以A-φ公式体系为例）：
 *          - K 来自 ∇×(ν∇×A) 项（磁阻率相关的curl-curl算子）
 *          - M 来自 σ∂A/∂t 项（瞬态时间导数的质量项）
 *          - C 来自 σA 项（涡流耗散阻尼项）
 *          - F 来自 J_src 源项（外加电流密度）
 */
struct ElementMatrices {
    Eigen::MatrixXd K;  ///< 刚度矩阵 K_e（静态/瞬态共用），维度: n_dof × n_dof
    Eigen::MatrixXd M;  ///< 质量矩阵 M_e（瞬态专用），维度: n_dof × n_dof
    Eigen::MatrixXd C;  ///< 阻尼矩阵 C_e（瞬态专用），维度: n_dof × n_dof
    Eigen::VectorXd F;  ///< 源项向量 F_e，维度: n_dof × 1

    /**
     * @brief 默认构造函数，初始化为未分配状态的空矩阵/向量
     */
    ElementMatrices()
    {
    }
};

// ==================== 抽象基类定义 ====================

/**
 * @class EMElementIntegratorBase
 * @brief 电磁场单元积分器抽象基类
 * @details 定义电磁场有限元分析中单元级矩阵/向量计算的统一接口框架，
 *          是数值计算层的核心抽象。所有具体积分器实现（基于Lagrange节点元、
 *          Nedelec棱边元等不同形函数类型）必须继承此基类并实现纯虚方法。
 *
 *          设计模式与架构定位：
 *          - 采用策略模式（Strategy Pattern）：通过虚函数分派到具体单元类型的积分算法
 *          - 位于形函数层（ShapeFunctionBase）之上、全局组装层之下
 *          - 与GaussQuadrature配合完成参考域上的高斯数值积分
 *
 *          继承关系：
 *          EMElementIntegratorBase（本类）
 *              ├── LagrangeEMIntegrator（Lagrange节点元积分器，用于标量场φ）
 *              └── NedelecEMIntegrator（Nedelec棱边元积分器，用于矢量场A）
 *
 *          协作组件：
 *          - ShapeFunctionBase：提供形函数值N_i、梯度∇N_i、旋度∇×N_i的计算
 *          - GaussQuadrature：提供高斯积分点和权重
 *          - MaterialProperties：提供材料参数（ε, μ, σ）及其非线性信息
 *
 *          时间复杂度说明（单次computeAllMatrices调用）：
 *          - 设n_gp为高斯积分点数，n_dof为单元自由度数
 *          - 刚度矩阵K: O(n_gp × n_dof²) —— 每个积分点计算n_dof²个元素
 *          - 质量矩阵M: O(n_gp × n_dof²)
 *          - 阻尼矩阵C: O(n_gp × n_dof²)
 *          - 源项向量F: O(n_gp × n_dof)
 *          - 总体复杂度: O(n_gp × n_dof²)
 *
 * @note 此类位于numeric命名空间下，属于数值计算层的核心模块
 * @see MaterialProperties 材料参数结构体
 * @see ElementMatrices 单元矩阵结果结构体
 * @see ShapeFunctionBase 形函数抽象基类
 * @see GaussQuadrature 高斯积分类
 */
class EMElementIntegratorBase {
public:
    /**
     * @brief 虚析构函数，确保派生类通过基类指针删除时正确释放资源
     */
    virtual ~EMElementIntegratorBase() = default;

    // ========== 纯虚方法：派生类必须实现 ==========

    /**
     * @brief 计算单元的全部矩阵和源项向量
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count，
     *                    每列对应一个节点的物理坐标 (x, y, z)^T
     * @return ElementMatrices 包含K（刚度）、M（质量）、C（阻尼）、F（源项）的结构体
     *
     * @details 这是积分器的核心入口方法，一次性完成所有单元矩阵的高斯积分计算。
     *          内部调用流程：
     *          1. 获取高斯积分点和权重（通过GaussQuadrature::getPoints）
     *          2. 对每个积分点ξ_q：
     *             a. 计算形函数值N_i(ξ_q)、梯度∇N_i(ξ_q)或旋度∇×N_i(ξ_q)
     *             b. 计算雅可比行列式detJ(ξ_q)用于体积元变换
     *             c. 获取有效材料参数ε_eff、μ_eff、σ_eff（考虑非线性）
     *             d. 计算被积函数并乘以权重w_q·|detJ|
     *          3. 累加所有积分点的贡献得到最终矩阵
     *
     *          数值积分公式（以刚度矩阵为例）：
     *          K_e = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · f_K(N, ∇N, material, ξ_q)
     *          其中f_K为刚度核函数，由派生类根据单元类型定义
     *
     * @note 调用前需确保已通过setMaterialProperties()设置材料参数
     * @warning 返回的矩阵尺寸由派生类的getNodeCount()决定，需与全局组装器协调
     */
    virtual ElementMatrices computeAllMatrices(
        const Eigen::MatrixXd& node_coords
    ) const = 0;

    /**
     * @brief 单独计算单元刚度矩阵
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 刚度矩阵 K_e，维度为 n_dof × n_dof
     *
     * @details 刚度矩阵对应电磁方程中的空间微分算子项。
     *          对于不同的物理场和形函数类型，刚度核函数不同：
     *          - 标量场（Lagrange元，电势φ）：K_ij = ∫ ε ∇N_i·∇N_j dΩ
     *          - 矢量场（Nedelec元，磁矢位A）：K_ij = ∫ (1/μ)(∇×N_i)·(∇×N_j) dΩ
     *
     * @note 当仅需刚度矩阵时（如静态场求解），调用此方法比computeAllMatrices更高效
     */
    virtual Eigen::MatrixXd computeStiffnessMatrix(
        const Eigen::MatrixXd& node_coords
    ) const = 0;

    /**
     * @brief 单独计算单元质量矩阵
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 质量矩阵 M_e，维度为 n_dof × n_dof
     *
     * @details 质量矩阵对应瞬态/时谐方程中的时间二阶导数项或一阶导数项。
     *          对于A-φ体系的瞬态涡流方程：
     *          M_ij = ∫ σ N_i·N_j dΩ  （电导率质量项，关联∂A/∂t）
     *          或对于波动方程：
     *          M_ij = ∫ ε N_i·N_j dΩ  （介电常数质量项，关联∂²A/∂t²）
     *
     *          具体的质量核函数由派生类根据所求解的物理方程决定。
     *
     * @note 仅在瞬态模式下（is_transient_=true）此矩阵有物理意义
     */
    virtual Eigen::MatrixXd computeMassMatrix(
        const Eigen::MatrixXd& node_coords
    ) const = 0;

    /**
     * @brief 单独计算单元阻尼矩阵
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::MatrixXd 阻尼矩阵 C_e，维度为 n_dof × n_dof
     *
     * @details 阻尼矩阵对应瞬态方程中的零阶（无时间导数）耗散项。
     *          对于A-φ体系的涡流方程：
     *          C_ij = ∫ σ N_i·N_j dΩ  （涡流损耗阻尼项）
     *
     *          在某些公式体系中，阻尼矩阵可能与质量矩阵使用相同的核函数
     *          但具有不同的物理系数（如频率因子jω）。
     *
     * @note 仅在瞬态模式下（is_transient_=true）此矩阵有物理意义
     */
    virtual Eigen::MatrixXd computeDampingMatrix(
        const Eigen::MatrixXd& node_coords
    ) const = 0;

    /**
     * @brief 单独计算单元源项向量
     * @param node_coords 节点坐标矩阵，维度为 dim × node_count
     * @return Eigen::VectorXd 源项向量 F_e，维度为 n_dof × 1
     *
     * @details 源项向量对应外部激励（如源电流密度J_src）在单元上的投影。
     *          一般形式：
     *          F_i = ∫ N_i · f_src dΩ
     *          其中f_src为源项密度函数（如电流密度J_exc、电荷密度ρ等）。
     *
     *          对于不同的激励类型：
     *          - 电流源激励：F_i = ∫ N_i · J_src dΩ
     *          - 电压源激励：通过边界条件施加，不直接进入F向量
     *
     * @note 若当前单元无外部激励，返回零向量
     */
    virtual Eigen::VectorXd computeSourceVector(
        const Eigen::MatrixXd& node_coords
    ) const = 0;

    // ========== 材料参数设置方法（带默认实现的虚方法）==========

    /**
     * @brief 设置单元的材料属性参数
     * @param props 材料属性结构体，包含ε、μ、σ及其非线性标志
     *
     * @details 将传入的材料参数存储到内部成员变量material_props_中，
     *          后续的矩阵计算将使用这些参数。
     *          典型调用场景：
     *          - 从网格文件读取每个单元的材料标签后设置
     *          - 多材料区域遍历时逐单元切换材料参数
     *
     * @code
     * MaterialProperties copper(8.854e-12, 1.257e-6, 5.8e7);
     * integrator->setMaterialProperties(copper);
     * auto mats = integrator->computeAllMatrices(node_coords);
     * @endcode
     *
     * @note 参数合法性检查：epsilon > 0, mu > 0, sigma >= 0
     * @warning 若参数非法，输出WARN日志并保持原值不变
     */
    virtual void setMaterialProperties(const MaterialProperties& props);

    /**
     * @brief 设置是否启用瞬态分析模式
     * @param is_transient true启用瞬态模式（计算M和C矩阵），false仅计算K和F
     *
     * @details 瞬态模式控制矩阵计算的范围：
     *          - 静态模式（is_transient=false）：仅计算刚度矩阵K和源项F
     *            适用于静磁场、静电场等稳态问题
     *          - 瞬态模式（is_transient=true）：计算全部矩阵K、M、C、F
     *            适用于涡流场、瞬态波动场等动态问题
     *
     *          性能影响：瞬态模式的计算量约为静态模式的2~3倍
     *          （需额外计算质量矩阵和阻尼矩阵的高斯积分）
     */
    virtual void setTransientMode(bool is_transient);

    // ========== 查询方法（带默认实现的虚方法）==========

    /**
     * @brief 查询当前是否处于瞬态分析模式
     * @return true 当前为瞬态模式（将计算M和C矩阵），false为静态模式
     */
    virtual bool isTransient() const;

    /**
     * @brief 获取当前设置的材料属性参数
     * @return MaterialProperties 当前生效的材料属性结构体副本
     */
    virtual MaterialProperties getMaterialProperties() const;

    // ========== 非线性材料接口（虚函数，默认返回线性值）==========

    /**
     * @brief 获取指定位置的有效介电常数（支持非线性覆写）
     * @param position 物理域位置向量 (x, y, z)^T，维度为3
     * @return double 该位置处的有效介电常数 ε_eff (F/m)
     *
     * @details 基类默认实现直接返回material_props_.epsilon的线性常量值。
     *          派生类可覆写此方法以支持以下非线性场景：
     *          - 空间非均匀材料：ε随位置变化（如分层介质）
     *          - 场依赖非线性：ε随电场强度E变化（铁电材料）
     *          - 温度依赖性：ε随温度T变化
     *
     *          非线性介电模型示例（派生类中实现）：
     *          ε_eff(E) = ε₀ · (1 + χ₁|E|² + χ₂|E|⁴ + ...)  （Kerr效应型）
     *          或 ε_eff(T) = ε_ref · [1 + α_ε(T - T_ref)]       （温度系数型）
     *
     * @note position参数的未使用分量（如2D问题的z分量）应忽略
     * @warning 仅当material_props_.is_nonlinear_epsilon=true时才应调用非线性版本
     */
    virtual double getEffectiveEpsilon(const Eigen::Vector3d& position) const;

    /**
     * @brief 获取指定位置的有效磁导率（支持非线性覆写）
     * @param position 物理域位置向量 (x, y, z)^T，维度为3
     * @return double 该位置处的有效磁导率 μ_eff (H/m)
     *
     * @details 基类默认实现直接返回material_props_.mu的线性常量值。
     *          派生类可覆写此方法以支持以下非线性场景：
     *          - 铁磁材料B-H曲线非线性：μ(B) = dB/dH，随磁感应强度B饱和
     *          - 磁滞效应：μ依赖于磁场历史状态
     *          - 各向异性：μ在不同方向取不同值（张量形式，需扩展接口）
     *
     *          常见非线性磁导率模型（派生类中实现）：
     *          - B-H曲线插值：给定离散(H_k, B_k)数据点，插值求μ = ΔB/ΔH
     *          - Langevin函数：μ_eff(H) = μ₀ · (1 + M_s/H · L(aH))
     *          - Frohlich-Kennely模型：μ_eff(H) = μ_max / (1 + βH²)
     *
     * @note 对于电机、变压器等含铁芯设备，非线性磁导率的准确建模至关重要
     * @warning 仅当material_props_.is_nonlinear_mu=true时才应调用非线性版本
     */
    virtual double getEffectiveMu(const Eigen::Vector3d& position) const;

    /**
     * @brief 获取指定位置的有效电导率（支持非线性覆写）
     * @param position 物理域位置向量 (x, y, z)^T，维度为3
     * @return double 该位置处的有效电导率 σ_eff (S/m)
     *
     * @details 基类默认实现直接返回material_props_.sigma的线性常量值。
     *          派生类可覆写此方法以支持以下非线性场景：
     *          - 温度依赖电导率：σ(T) = σ_ref / [1 + α(T - T_ref)]（金属电阻温升）
     *          - 场依赖电导率：半导体材料的σ(E)非线性（碰撞电离等效应）
     *          - 各向异性电导：石墨烯等材料的方向依赖电导率
     *
     *          常见非线性电导率模型：
     *          - 线性温度模型：σ(T) = σ₀ / (1 + γΔT)
     *          - Power Law模型：σ(J) = σ₀ · |J/J_c|^(n-1) （超导临界电流）
     *
     * @note 电导率为0时表示理想绝缘体，此时阻尼矩阵C为零矩阵
     * @warning 仅当material_props_.is_nonlinear_sigma=true时才应调用非线性版本
     */
    virtual double getEffectiveSigma(const Eigen::Vector3d& position) const;

protected:
    MaterialProperties material_props_;   ///< 当前单元的材料属性参数
    bool is_transient_;                  ///< 是否启用瞬态分析模式，默认false
};

// ==================== 默认方法实现 ====================

inline void EMElementIntegratorBase::setMaterialProperties(
    const MaterialProperties& props
) {
    // 参数合法性检查：介电常数和磁导率必须为正，电导率必须非负
    if (props.epsilon <= 0.0) {
        FEEM_WARN("setMaterialProperties: 介电常数epsilon={}必须为正数，拒绝设置", props.epsilon);
        return;
    }
    if (props.mu <= 0.0) {
        FEEM_WARN("setMaterialProperties: 磁导率mu={}必须为正数，拒绝设置", props.mu);
        return;
    }
    if (props.sigma < 0.0) {
        FEEM_WARN("setMaterialProperties: 电导率sigma={}不能为负数，拒绝设置", props.sigma);
        return;
    }

    material_props_ = props;
    FEEM_DEBUG("setMaterialProperties: 材料参数已更新 ε={:.4e}, μ={:.4e}, σ={:.4e}",
               props.epsilon, props.mu, props.sigma);
}

inline void EMElementIntegratorBase::setTransientMode(bool is_transient) {
    is_transient_ = is_transient;
    FEEM_DEBUG("setTransientMode: 瞬态模式设置为{}", is_transient ? "开启" : "关闭");
}

inline bool EMElementIntegratorBase::isTransient() const {
    return is_transient_;
}

inline MaterialProperties EMElementIntegratorBase::getMaterialProperties() const {
    return material_props_;
}

inline double EMElementIntegratorBase::getEffectiveEpsilon(
    const Eigen::Vector3d& position
) const {
    // 默认返回线性常量值，忽略位置参数
    (void)position;
    return material_props_.epsilon;
}

inline double EMElementIntegratorBase::getEffectiveMu(
    const Eigen::Vector3d& position
) const {
    // 默认返回线性常量值，忽略位置参数
    (void)position;
    return material_props_.mu;
}

inline double EMElementIntegratorBase::getEffectiveSigma(
    const Eigen::Vector3d& position
) const {
    // 默认返回线性常量值，忽略位置参数
    (void)position;
    return material_props_.sigma;
}

} // namespace numeric
