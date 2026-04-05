/**
 * @file electrostatic_integrator.cpp
 * @brief 数值计算层 - 静电场/瞬态电场单元积分器完整实现
 * @details 实现基于Lagrange节点元的静电场和瞬态电场有限元单元矩阵组装，
 *          包括刚度矩阵、质量矩阵、阻尼矩阵和源项向量的高斯数值积分计算。
 *
 * 数学公式参考：
 * - 刚度矩阵: K_e = ∫_Ω (∇N)^T · ε · (∇N) dΩ
 * - 质量矩阵: M_e = ∫_Ω N^T · ε · N dΩ
 * - 阻尼矩阵: C_e = ∫_Ω N^T · σ · N dΩ
 * - 源项向量: F_e = ∫_Ω N · ρ_e dΩ  （当前实现 ρ_e = 0）
 *
 * 支持的单元类型：
 * - 二维: TRI3, QUAD4
 * - 三维: TET4, HEX8, PRISM6, PRISM15, PYRAMID5, PYRAMID13
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#include "electrostatic_integrator.hpp"
#include <cmath>

namespace numeric {

// ==================== 构造函数与析构函数 ====================

/**
 * @brief 构造函数实现：初始化形函数对象并设置默认积分阶数
 * @param element_type 目标单元类型枚举值
 *
 * @details 执行步骤：
 * 1. 存储单元类型到 element_type_
 * 2. 通过 elementTypeToString() 将枚举转为字符串
 * 3. 调用 ShapeFunctionFactory::create() 创建 Lagrange 形函数对象
 * 4. 通过 getDefaultGaussOrder() 设置默认积分阶数
 * 5. 初始化瞬态模式标志为 false（静态分析默认）
 *
 * 错误处理：
 * - 若单元类型不受支持，回退使用 TRI3 并输出 ERROR 日志
 * - 若工厂创建失败（返回 nullptr），输出 ERROR 日志
 */
ElectrostaticIntegrator::ElectrostaticIntegrator(ElementType element_type)
    : element_type_(element_type)
    , gauss_order_(1)
    , shape_func_(nullptr)
{
    // 将 ElementType 枚举转换为字符串标识符，用于工厂方法调用
    std::string type_str = elementTypeToString(element_type_);

    if (type_str.empty()) {
        FEEM_ERROR("ElectrostaticIntegrator: 不支持的单元类型({}), "
                   "回退使用TRI3", static_cast<int>(element_type_));
        element_type_ = ElementType::TRI3;
        type_str = "TRI3";
    }

    // 通过形函数工厂创建对应的 Lagrange 单元对象
    shape_func_ = ShapeFunctionFactory::create(type_str);

    if (!shape_func_) {
        FEEM_ERROR("ElectrostaticIntegrator: 形函数创建失败(type={})", type_str);
        return;
    }

    // 根据单元类型自动选择默认高斯积分阶数
    gauss_order_ = getDefaultGaussOrder(element_type_);

    // 默认静态模式（基类 is_transient_ 已在默认构造中初始化为 false）

    FEEM_INFO("ElectrostaticIntegrator 初始化完成: type={}, nodes={}, "
              "dim={}, gauss_order={}",
              type_str, shape_func_->getNodeCount(),
              shape_func_->getDim(), gauss_order_);
}

// ==================== 核心矩阵计算方法 ====================

/**
 * @brief 计算全部单元矩阵（入口方法）
 * @param node_coords 节点坐标矩阵 (dim × node_count)
 * @return ElementMatrices 包含 K, M, C, F 的完整结果
 *
 * @details 按顺序计算四个矩阵/向量并组装到 ElementMatrices 结构体中。
 *          此方法的优点是避免重复计算公共中间量（雅可比矩阵等），
 *          但当前实现为清晰起见仍分别调用各独立计算方法。
 */
ElementMatrices ElectrostaticIntegrator::computeAllMatrices(
    const Eigen::MatrixXd& node_coords) const
{
    ElementMatrices result;

    // 刚度矩阵 K_e：介电项 ∫ ε ∇N·∇N dΩ（静态/瞬态共用）
    result.K = computeStiffnessMatrix(node_coords);

    // 质量矩阵 M_e：介电惯性项 ∫ ε N·N dΩ（瞬态/通用）
    result.M = computeMassMatrix(node_coords);

    // 阻尼矩阵 C_e：电导耗散项 ∫ σ N·N dΩ（仅瞬态有意义）
    result.C = computeDampingMatrix(node_coords);

    // 源项向量 F_e：电荷激励项 ∫ N·ρ_e dΩ（当前零源项）
    result.F = computeSourceVector(node_coords);

    FEEM_DEBUG("computeAllMatrices 完成: K={}×{}, M={}×{}, C={}×{}, F={}",
               result.K.rows(), result.K.cols(),
               result.M.rows(), result.M.cols(),
               result.C.rows(), result.C.cols(),
               result.F.size());

    return result;
}

/**
 * @brief 计算单元刚度矩阵 K_e
 * @param node_coords 节点坐标矩阵 (dim × node_count)
 * @return 刚度矩阵 (n_nodes × n_nodes)
 *
 * @details 数值积分公式：
 *   K_e = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ε_eff(x_q) · G_q^T · G_q
 *
 *   其中：
 *   - ξ_q: 第 q 个高斯积分点的局部坐标
 *   - w_q: 第 q 个积分点的权重
 *   - detJ(ξ_q): ξ_q 处的雅可比行列式
 *   - ε_eff(x_q): ξ_q 映射到物理域位置 x_q 处的有效介电常数
 *   - G_q = ∇N(ξ_q): ξ_q 处的物理域形函数梯度矩阵 (n_nodes × dim)
 *
 *   矩阵组装过程（对每个积分点 q）：
 *   1. 从 GaussQuadrature 获取积分点坐标 ξ_q 和权重 w_q
 *   2. 构造 LocalPoint 对象表示局部坐标
 *   3. 通过 shape_func_->calcPhysicalGradN() 计算 G_q
 *   4. 通过 shape_func_->calcJacobian() 获取 detJ
 *   5. 将 ξ_q 映射到物理坐标 x_q（通过节点坐标插值）
 *   6. 调用 getEffectiveEpsilon(x_q) 获取有效介电常数
 *   7. 累加贡献: K += (G_q^T · ε · G_q) · |detJ| · w_q
 */
Eigen::MatrixXd ElectrostaticIntegrator::computeStiffnessMatrix(
    const Eigen::MatrixXd& node_coords) const
{
    // 获取单元自由度数（节点数），初始化为零矩阵
    int n_nodes = shape_func_->getNodeCount();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_nodes, n_nodes);

    // 获取当前单元类型对应的高斯积分点和权重
    auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeStiffnessMatrix: 无法获取单元类型({})的积分点, "
                   "order={}", static_cast<int>(element_type_), gauss_order_);
        return K;
    }

    // 遍历所有高斯积分点，累加各点的被积函数贡献
    for (const auto& gp : gauss_points) {
        // 构造局部坐标点（从 GaussPoint 的 coords 向量提取分量）
        LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 计算物理域形函数梯度矩阵 G_q (n_nodes × dim)
        // 公式: G_q = (∂N/∂ξ) · J⁻¹，将参考域梯度映射到物理域
        Eigen::MatrixXd gradN = shape_func_->calcPhysicalGradN(xi, node_coords);

        // 计算雅可比矩阵及其行列式，用于体积元变换 dV = |detJ| dξ
        JacobianResult jaco = shape_func_->calcJacobian(xi, node_coords);

        // 检查雅可比行列式的合法性（畸形单元检测）
        if (std::abs(jaco.det_j) < 1e-15) {
            FEEM_WARN("computeStiffnessMatrix: 雅可比行列式接近零({:.6e}), "
                      "可能存在畸形单元, 积分点({}, {}, {})",
                      jaco.det_j, gp.coords[0], gp.coords[1], gp.coords[2]);
            continue;  // 跳过此积分点以避免数值不稳定
        }

        // 将积分点局部坐标映射到物理域坐标（用于非线性材料参数查询）
        // 物理坐标: x_q = N(ξ_q) · X_node （形函数值乘以节点坐标矩阵）
        Eigen::VectorXd N_vals = shape_func_->evalN(xi);
        Eigen::Vector3d physical_pos = (node_coords * N_vals).head<3>();

        // 获取积分点处的有效介电常数（支持非线性材料模型覆写）
        double eps_effective = getEffectiveEpsilon(physical_pos);

        // 组装刚度矩阵当前积分点的贡献:
        // ΔK_q = G_q^T · ε_eff · G_q · |detJ| · w_q
        //
        // 维度分析:
        //   G_q^T:      (dim × n_nodes)
        //   ε_eff:     标量（各向同性介质，若需各向异性可替换为张量）
        //   G_q:       (n_nodes × dim)
        //   |detJ|·w:  标量
        //   结果ΔK_q:  (n_nodes × n_nodes)
        K += (gradN.transpose() * eps_effective * gradN)
             * std::abs(jaco.det_j) * gp.weight;
    }

    FEEM_DEBUG("computeStiffnessMatrix: K矩阵范数={:.6e}, 条件数估计可用",
               K.norm());

    return K;
}

/**
 * @brief 计算单元质量矩阵 M_e
 * @param node_coords 节点坐标矩阵 (dim × node_count)
 * @return 质量矩阵 (n_nodes × n_nodes)
 *
 * @details 数值积分公式：
 *   M_e = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ε_eff(x_q) · N_q · N_q^T
 *
 *   其中 N_q 为第 q 个积分点处的形函数值向量 (n_nodes × 1)，
 *   外积 N_q · N_q^T 产生 (n_nodes × n_nodes) 矩阵。
 *
 *   物理意义：
 *   质量矩阵对应瞬态方程中的介电惯性项（位移电流相关），
 *   在静电场静态分析中不参与求解，但在时域波动问题中至关重要。
 *   对于准静态近似（忽略位移电流），M 矩阵通常不出现或设为零。
 *
 *   与阻尼矩阵的区别：
 *   M 使用介电常数 ε 作为系数（描述场的"惯性"）
 *   C 使用电导率 σ 作为系数（描述能量耗散）
 */
Eigen::MatrixXd ElectrostaticIntegrator::computeMassMatrix(
    const Eigen::MatrixXd& node_coords) const
{
    int n_nodes = shape_func_->getNodeCount();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_nodes, n_nodes);

    auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeMassMatrix: 无法获取积分点, type={}, order={}",
                   static_cast<int>(element_type_), gauss_order_);
        return M;
    }

    for (const auto& gp : gauss_points) {
        LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 计算形函数值向量 N_q (n_nodes × 1)
        // 每个 N_i 表示第 i 个节点形函数在该积分点处的值
        Eigen::VectorXd N = shape_func_->evalN(xi);

        // 计算雅可比行列式用于体积元变换
        JacobianResult jaco = shape_func_->calcJacobian(xi, node_coords);

        if (std::abs(jaco.det_j) < 1e-15) {
            FEEM_WARN("computeMassMatrix: 雅可比行列式接近零, 跳过积分点");
            continue;
        }

        // 映射到物理坐标以获取有效的材料参数
        Eigen::Vector3d physical_pos = (node_coords * N).head<3>();
        double eps_effective = getEffectiveEpsilon(physical_pos);

        // 组装质量矩阵贡献:
        // ΔM_q = N_q · ε_eff · N_q^T · |detJ| · w_q
        //
        // 外积形式: (n_nodes×1) · 标量 · (1×n_nodes) → (n_nodes×n_nodes)
        // 结果为对称正定矩阵（当 ε > 0 时）
        M += (N * eps_effective * N.transpose())
             * std::abs(jaco.det_j) * gp.weight;
    }

    FEEM_DEBUG("computeMassMatrix: M矩阵Frobenius范数={:.6e}", M.norm());

    return M;
}

/**
 * @brief 计算单元阻尼矩阵 C_e
 * @param node_coords 节点坐标矩阵 (dim × node_count)
 * @return 阻尼矩阵 (n_nodes × n_nodes)
 *
 * @details 数值积分公式：
 *   C_e = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · σ_eff(x_q) · N_q · N_q^T
 *
 *   其中 σ_eff(x_q) 为积分点处的有效电导率，通过 getEffectiveSigma() 获取，
 *   支持非线性电导率模型（如温度依赖、场依赖等）。
 *
 *   物理背景：
 *   阻尼矩阵来源于欧姆定律 J = σE 和电流连续性方程，
 *   描述导体中的焦耳热损耗。在电路类比中等效于电导矩阵 G = C。
 *
 *   特殊情况：
 *   - 理想绝缘体 (σ = 0): C 为零矩阵，无耗散
 *   - 良导体 (σ >> 0): C 矩阵元素值大，主导瞬态响应
 *   - 非均匀材料: σ 随空间变化，各积分点取不同值
 *
 *   注意：C 矩阵与 M 矩阵具有相同的数学结构（均为 N·N^T 型），
 *         区别仅在于系数不同（σ vs ε）。
 */
Eigen::MatrixXd ElectrostaticIntegrator::computeDampingMatrix(
    const Eigen::MatrixXd& node_coords) const
{
    int n_nodes = shape_func_->getNodeCount();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_nodes, n_nodes);

    auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeDampingMatrix: 无法获取积分点, type={}, order={}",
                   static_cast<int>(element_type_), gauss_order_);
        return C;
    }

    for (const auto& gp : gauss_points) {
        LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 计算形函数值向量 N_q (n_nodes × 1)
        Eigen::VectorXd N = shape_func_->evalN(xi);

        // 计算雅可比行列式
        JacobianResult jaco = shape_func_->calcJacobian(xi, node_coords);

        if (std::abs(jaco.det_j) < 1e-15) {
            FEEM_WARN("computeDampingMatrix: 雅可比行列式接近零, 跳过积分点");
            continue;
        }

        // 映射到物理坐标以获取有效电导率（支持非线性模型）
        Eigen::Vector3d physical_pos = (node_coords * N).head<3>();
        double sigma_effective = getEffectiveSigma(physical_pos);

        // 组装阻尼矩阵贡献:
        // ΔC_q = N_q · σ_eff · N_q^T · |detJ| · w_q
        //
        // 当 sigma_effective = 0 时（理想绝缘体），此项贡献为零
        // 当 sigma_effective > 0 时，产生对称半正定的耗散矩阵
        C += (N * sigma_effective * N.transpose())
             * std::abs(jaco.det_j) * gp.weight;
    }

    FEEM_DEBUG("computeDampingMatrix: C矩阵Frobenius范数={:.6e}", C.norm());

    return C;
}

/**
 * @brief 计算单元源项向量 F_e
 * @param node_coords 节点坐标矩阵 (dim × node_count)
 * @return 源项向量 (n_nodes × 1)
 *
 * @details 当前实现采用零源项假设 (ρ_e ≡ 0)，适用于以下场景：
 *   - 拉普拉斯方程求解区域（无自由电荷分布）
 *   - 电极边界条件驱动的电场问题
 *   - 电荷中性区域的电势分布计算
 *
 *   一般形式的源项积分公式（供后续扩展参考）：
 *   F_i = ∫_Ωe N_i(x) · ρ_e(x) dΩ
 *       = Σ_{q=1}^{n_gp} w_q · |detJ(ξ_q)| · ρ_e(x_q) · N_i(ξ_q)
 *
 *   后续可扩展支持的源项类型：
 *   - 均匀体电荷密度: ρ_e = const
 *   - 点电荷: ρ_e(x) = Q · δ(x - x₀) （需特殊处理）
 *   - 线电荷/面电荷: 通过降维积分处理
 *   - 空间变化的电荷分布: ρ_e(x,y,z) 用户自定义函数
 */
Eigen::VectorXd ElectrostaticIntegrator::computeSourceVector(
    const Eigen::MatrixXd& node_coords) const
{
    // 获取节点数量，初始化零向量
    int n_nodes = shape_func_->getNodeCount();
    Eigen::VectorXd F = Eigen::VectorXd::Zero(n_nodes);

    // 当前版本实现零源项（ρ_e = 0）
    // 无需执行高斯积分，直接返回零向量
    //
    // 如需启用非零源项，取消下方注释并实现 charge_density 函数:
    //
    // auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);
    // for (const auto& gp : gauss_points) {
    //     LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);
    //     Eigen::VectorXd N = shape_func_->evalN(xi);
    //     JacobianResult jaco = shape_func_->calcJacobian(xi, node_coords);
    //     double rho = getChargeDensity(/* position */);  // 用户定义的电荷密度
    //     F += N * rho * std::abs(jaco.det_j) * gp.weight;
    // }

    (void)node_coords;  // 抑制未使用参数警告

    FEEM_DEBUG("computeSourceVector: 返回零源项向量(size={})", n_nodes);

    return F;
}

// ==================== 查询与配置方法 ====================

/**
 * @brief 获取当前单元类型
 * @return ElementType 枚举值
 */
ElementType ElectrostaticIntegrator::getElementType() const {
    return element_type_;
}

/**
 * @brief 获取单元节点总数
 * @return int 节点数量
 */
int ElectrostaticIntegrator::getNodeCount() const {
    if (!shape_func_) {
        FEEM_WARN("getNodeCount: 形函数对象为空, 返回0");
        return 0;
    }
    return shape_func_->getNodeCount();
}

/**
 * @brief 获取单元空间维度
 * @return int 维度值（2 或 3）
 */
int ElectrostaticIntegrator::getDimension() const {
    if (!shape_func_) {
        FEEM_WARN("getDimension: 形函数对象为空, 返回0");
        return 0;
    }
    return shape_func_->getDim();
}

/**
 * @brief 获取当前高斯积分阶数
 * @return int 积分阶数
 */
int ElectrostaticIntegrator::getGaussOrder() const {
    return gauss_order_;
}

/**
 * @brief 设置新的高斯积分阶数
 * @param order 目标积分阶数（必须为正整数）
 *
 * @details 更改积分阶数的注意事项：
 *   - 较低阶数: 计算速度快，但对曲边单元或高精度需求可能不足
 *   - 较高阶数: 积分精度高，但计算开销随阶数呈多项式增长
 *   - 推荐保持默认值，除非有明确的精度或性能调优需求
 */
void ElectrostaticIntegrator::setGaussOrder(int order) {
    if (order <= 0) {
        FEEM_WARN("setGaussOrder: 积分阶数必须为正整数(收到order={}), "
                  "拒绝修改", order);
        return;
    }
    int old_order = gauss_order_;
    gauss_order_ = order;
    FEEM_DEBUG("setGaussOrder: 积分阶数从{}更改为{}", old_order, order);
}

// ==================== 私有辅助方法 ====================

/**
 * @brief 将 ElementType 枚举值转换为对应的字符串标识符
 * @param type 单元类型枚举值
 * @return std::string 类型名称字符串，不支持则返回空串
 *
 * @details 映射表覆盖全部 ElectrostaticIntegrator 支持的 Lagrange 单元类型。
 *          返回的字符串可直接传入 ShapeFunctionFactory::create() 方法。
 *
 *          支持的类型映射：
 *          二维: TRI3, QUAD4
 *          三维: TET4, HEX8, PRISM6, PRISM15, PYRAMID5, PYRAMID13
 */
std::string ElectrostaticIntegrator::elementTypeToString(ElementType type) {
    // 二维单元
    switch (type) {
        case ElementType::TRI3:  return "TRI3";
        case ElementType::QUAD4: return "QUAD4";

        // 三维单元
        case ElementType::TET4:    return "TET4";
        case ElementType::HEX8:    return "HEX8";
        case ElementType::PRISM6:  return "PRISM6";
        case ElementType::PRISM15: return "PRISM15";
        case ElementType::PYRAMID5:  return "PYRAMID5";
        case ElementType::PYRAMID13: return "PYRAMID13";

        default:
            return "";
    }
}

/**
 * @brief 根据单元类型确定默认的高斯积分阶数
 * @param type 单元类型枚举值
 * @return int 推荐的默认积分阶数
 *
 * @details 各单元类型的默认积分方案基于有限元理论中的精度要求：
 *
 *   精度说明（对于线性单元，刚度矩阵的被积函数为常数或线性）：
 *   - TRI3 (线性三角形): 1点质心积分精确到常数项
 *   - QUAD4 (双线性四边形): 2×2=4点 Gauss-Legendre 精确到双三次
 *   - TET4 (线性四面体): 1点质心积分精确到常数项
 *   - HEX8 (三线性六面体): 2×2×2=8点精确到三三次
 *   - PRISM6 (线性棱柱): 三角形3点×ζ方向2点=6点
 *   - PRISM15 (二次棱柱): 同PRISM6（更高精度可选但非必需）
 *   - PYRAMID5 (线性金字塔): 5点工程近似积分
 *   - PYRAMID13 (二次金字塔): 同PYRAMID5
 *
 * @note 这些默认值适用于大多数工程应用场景；
 *       特殊情况（如高度畸变网格）可能需要提高积分阶数
 */
int ElectrostaticIntegrator::getDefaultGaussOrder(ElementType type) {
    switch (type) {
        case ElementType::TRI3:     return 1;   // 质心单点积分
        case ElementType::QUAD4:    return 4;   // 2×2 Gauss-Legendre
        case ElementType::TET4:     return 1;   // 质心单点积分
        case ElementType::HEX8:     return 8;   // 2×2×2 Gauss-Legendre
        case ElementType::PRISM6:   return 6;   // 三角形3点 × ζ方向2点
        case ElementType::PRISM15:  return 6;   // 同 PRISM6（二次单元默认）
        case ElementType::PYRAMID5: return 5;   // 5点近似积分
        case ElementType::PYRAMID13:return 5;   // 同 PYRAMID5
        default:                    return 1;   // 回退默认值
    }
}

} // namespace numeric
