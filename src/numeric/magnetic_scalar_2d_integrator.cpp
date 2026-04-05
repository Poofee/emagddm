/**
 * @file magnetic_scalar_2d_integrator.cpp
 * @brief 数值计算层 - 二维静磁场/瞬态磁场标量位积分器实现
 * @details 实现MagneticScalar2DIntegrator类的全部方法，
 *          包括刚度矩阵、质量矩阵、阻尼矩阵和源项向量的高斯积分计算。
 *
 *          核心算法流程（以刚度矩阵为例）：
 *          1. 获取高斯积分点集 {ξ_q, w_q}, q=1..n_gp
 *          2. 对每个积分点：
 *             a. 计算物理域形函数梯度 G_q = calcPhysicalGradN(ξ_q, node_coords)
 *             b. 计算雅可比行列式 detJ_q 用于体积元变换
 *             c. 获取有效磁导率 μ_eff = getEffectiveMu(x_q)
 *             d. 累加贡献: K += w_q · |detJ_q| · μ_eff · G_q^T · G_q
 *          3. 返回组装完成的单元矩阵
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#include "magnetic_scalar_2d_integrator.hpp"

namespace numeric {

// ==================== 构造函数与析构函数 ====================

MagneticScalar2DIntegrator::MagneticScalar2DIntegrator(ElementType element_type)
    : element_type_(element_type)
{
    // 校验单元类型：仅支持2D Lagrange单元（TRI3和QUAD4）
    if (element_type != ElementType::TRI3 && element_type != ElementType::QUAD4) {
        FEEM_ERROR("MagneticScalar2DIntegrator: 不支持的单元类型={} "
                   "(仅支持TRI3和QUAD4两种2D单元)",
                   static_cast<int>(element_type));
        throw std::invalid_argument(
            "MagneticScalar2DIntegrator: 仅支持2D单元类型(TRI3, QUAD4)，"
            "传入类型: " + std::to_string(static_cast<int>(element_type)));
    }

    // 实例化二维Lagrange形函数对象（模板参数Dim=2）
    shape_func_ = std::make_unique<LagrangeElement<2>>(element_type);

    // 设置默认材料参数：真空磁导率μ₀、零电导率、线性材料
    MaterialProperties default_props;
    default_props.mu = 4.0 * M_PI * 1e-7;       // 真空磁导率 μ₀ = 4π×10⁻⁷ H/m
    default_props.sigma = 0.0;                    // 默认无耗散（绝缘体/空气）
    default_props.epsilon = 8.854187817e-12;     // 介电常数（本模块不使用，但需设置）
    default_props.is_nonlinear_mu = false;        // 默认线性磁导率
    default_props.is_nonlinear_sigma = false;     // 默认线性电导率
    material_props_ = default_props;

    // 默认关闭瞬态模式（静态磁场分析）
    is_transient_ = false;

    FEEM_INFO("MagneticScalar2DIntegrator 初始化完成: 单元类型={}, 节点数={}, "
              "默认μ₀={:.4e} H/m",
              static_cast<int>(element_type_),
              shape_func_->getNodeCount(),
              material_props_.mu);
}

// ==================== 核心计算方法实现 ====================

ElementMatrices MagneticScalar2DIntegrator::computeAllMatrices(
    const Eigen::MatrixXd& node_coords) const {
    ElementMatrices result;

    // 计算刚度矩阵 K_m（静态/瞬态共用项）
    result.K = computeStiffnessMatrix(node_coords);

    // 若启用瞬态模式，额外计算质量矩阵 M_m（涡流耗散项）
    if (is_transient_) {
        result.M = computeMassMatrix(node_coords);
    } else {
        int n_dof = shape_func_->getNodeCount();
        result.M = Eigen::MatrixXd::Zero(n_dof, n_dof);
    }

    // 阻尼矩阵 C 在标量磁势模型中始终为零（无独立零阶耗散项）
    int n_dof = shape_func_->getNodeCount();
    result.C = Eigen::MatrixXd::Zero(n_dof, n_dof);

    // 源项向量 F_m（当前为零源项，无等效磁荷）
    result.F = computeSourceVector(node_coords);

    FEEM_DEBUG("computeAllMatrices: 矩阵组装完成 K[{}×{}], M[{}×{}], C[{}×{}], F[{}]",
               result.K.rows(), result.K.cols(),
               result.M.rows(), result.M.cols(),
               result.C.rows(), result.C.cols(),
               result.F.size());

    return result;
}

Eigen::MatrixXd MagneticScalar2DIntegrator::computeStiffnessMatrix(
    const Eigen::MatrixXd& node_coords) const {
    // 获取单元自由度数（等于节点数，Lagrange元每节点1个DOF）
    int n_dof = shape_func_->getNodeCount();

    // 初始化刚度矩阵为零矩阵（n_dof × n_dof，对称结构）
    Eigen::MatrixXd K_m = Eigen::MatrixXd::Zero(n_dof, n_dof);

    // 获取当前单元类型对应的高斯积分点和权重
    int gauss_order = getGaussOrder();
    std::vector<GaussPoint> gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order);

    // 检查积分点获取是否成功
    if (gauss_points.empty()) {
        FEEM_ERROR("computeStiffnessMatrix: 无法获取单元类型={}的高斯积分点(order={})",
                   static_cast<int>(element_type_), gauss_order);
        return K_m;
    }

    // 遍历所有高斯积分点，累加被积函数的数值近似值
    for (const auto& gp : gauss_points) {
        // 构造参考域局部坐标点
        LocalPoint xi(gp.coords[0], gp.coords[1]);

        // 计算物理域形函数梯度矩阵 G_q = ∇_x N(ξ_q)，维度: n_dof × 2
        // 内部调用链: evalGradN → calcJacobian → 梯度映射 ∇_x N = (∇_ξ N)·J⁻¹
        Eigen::MatrixXd grad_n_physical = shape_func_->calcPhysicalGradN(xi, node_coords);

        // 计算雅可比行列式 |detJ|，用于参考域到物理域的体积元变换
        JacobianResult jaco_result = shape_func_->calcJacobian(xi, node_coords);
        double det_j = std::abs(jaco_result.det_j);

        // 安全检查：防止畸形单元导致detJ接近零（可能引起数值不稳定）
        if (det_j < 1e-15) {
            FEEM_WARN("computeStiffnessMatrix: 雅可比行列式接近零(detJ={:.2e})，"
                      "可能存在畸形单元，跳过此积分点", det_j);
            continue;
        }

        // 计算积分点的物理坐标位置（用于非线性材料参数查询）
        Eigen::Vector3d position = node_coords * shape_func_->evalN(xi);

        // 获取该位置的有效磁导率（支持非线性B-H曲线覆写）
        double mu_eff = getEffectiveMu(position);

        // 组装刚度矩阵的被积函数核:
        // K_m += w_q · |detJ| · μ_eff · G_q^T · G_q
        //
        // 数学推导说明：
        // - G_q 为 n_dof × 2 矩阵（每行是一个节点的物理梯度 [∂N_i/∂x, ∂N_i/∂y]）
        // - G_q^T · G_q 为 2 × 2 矩阵（梯度外积的压缩形式）
        // - 最终通过 outer product 展开为 n_dof × n_dof 的完整刚度矩阵
        //
        // 等价展开形式（逐元素）：
        // K_m[i][j] += w_q · |detJ| · μ_eff · Σ_{k=1}^{2} (∂N_i/∂x_k)(∂N_j/∂x_k)
        //              = w_q · |detJ| · μ_eff · (∇N_i · ∇N_j)   （内积形式）
        K_m += gp.weight * det_j * mu_eff * (grad_n_physical * grad_n_physical.transpose());
    }

    FEEM_DEBUG("computeStiffnessMatrix: K_m[{}×{}] 组装完成, ||K||_F={:.6e}",
               n_dof, n_dof, K_m.norm());

    return K_m;
}

Eigen::MatrixXd MagneticScalar2DIntegrator::computeMassMatrix(
    const Eigen::MatrixXd& node_coords) const {
    // 获取单元自由度数
    int n_dof = shape_func_->getNodeCount();

    // 初始化质量矩阵为零矩阵（n_dof × n_dof，对称正定结构）
    Eigen::MatrixXd M_m = Eigen::MatrixXd::Zero(n_dof, n_dof);

    // 快速路径：若电导率为零（绝缘体），直接返回零矩阵（无需积分计算）
    double sigma_eff_base = material_props_.sigma;
    if (sigma_eff_base == 0.0 && !material_props_.is_nonlinear_sigma) {
        FEEM_DEBUG("computeMassMatrix: 电导率σ=0且无非线性标记，返回零矩阵");
        return M_m;
    }

    // 获取高斯积分点和权重
    int gauss_order = getGaussOrder();
    std::vector<GaussPoint> gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeMassMatrix: 无法获取单元类型={}的高斯积分点(order={})",
                   static_cast<int>(element_type_), gauss_order);
        return M_m;
    }

    // 遍历所有高斯积分点，计算质量矩阵的数值积分
    for (const auto& gp : gauss_points) {
        LocalPoint xi(gp.coords[0], gp.coords[1]);

        // 计算形函数值向量 N_q = N(ξ_q)，维度: n_dof × 1
        // 第i个元素 N_i(ξ_q) 表示第i个形函数在积分点q处的插值权重
        Eigen::VectorXd N_q = shape_func_->evalN(xi);

        // 计算雅可比行列式用于体积元变换
        JacobianResult jaco_result = shape_func_->calcJacobian(xi, node_coords);
        double det_j = std::abs(jaco_result.det_j);

        if (det_j < 1e-15) {
            FEEM_WARN("computeMassMatrix: 雅可比行列式接近零(detJ={:.2e})，跳过此积分点", det_j);
            continue;
        }

        // 计算积分点物理坐标（用于非线性电导率查询）
        Eigen::Vector3d position = node_coords * N_q;

        // 获取有效电导率（支持温度依赖等非线性模型）
        double sigma_eff = getEffectiveSigma(position);

        // 组装质量矩阵的被积函数核:
        // M_m += w_q · |detJ| · σ_eff · N_q · N_q^T
        //
        // 物理意义说明：
        // - N_q · N_q^T 形成 n_dof × n_dof 矩阵（形函数值的outer product）
        // - 第(i,j)元素 = N_i(ξ_q) · N_j(ξ_q)，表示节点i和j在积分点q处的耦合强度
        // - 加权后累加到全局质量矩阵中
        //
        // 与时间离散化的关系：
        // 半离散方程: M_m · {dφ/dt} + K_m · {φ} = {F}
        // 其中 M_m 描述了涡流场的"惯性"特性（能量耗散的时间常数）
        M_m += gp.weight * det_j * sigma_eff * (N_q * N_q.transpose());
    }

    FEEM_DEBUG("computeMassMatrix: M_m[{}×{}] 组装完成, ||M||_F={:.6e}",
               n_dof, n_dof, M_m.norm());

    return M_m;
}

Eigen::MatrixXd MagneticScalar2DIntegrator::computeDampingMatrix(
    const Eigen::MatrixXd& node_coords) const {
    // 标量磁势φ_m的公式体系中不存在独立的零阶阻尼项C
    // （与矢量磁势A法的σA阻尼项不同）
    // 返回零矩阵以保持基类接口一致性
    (void)node_coords;

    int n_dof = shape_func_->getNodeCount();
    FEEM_DEBUG("computeDampingMatrix: 标量磁势模型无阻尼项，返回{}×{}零矩阵",
               n_dof, n_dof);

    return Eigen::MatrixXd::Zero(n_dof, n_dof);
}

Eigen::VectorXd MagneticScalar2DIntegrator::computeSourceVector(
    const Eigen::MatrixXd& node_coords) const {
    // 当前版本假设无等效磁荷源（ρ_m = 0），返回零向量
    // 后续可扩展支持永磁体等效磁荷模型: ρ_m = -∇·M（磁化强度的散度）
    (void)node_coords;

    int n_dof = shape_func_->getNodeCount();
    FEEM_DEBUG("computeSourceVector: 当前为零源项模式(ρ_m=0)，返回{}维零向量", n_dof);

    return Eigen::VectorXd::Zero(n_dof);
}

// ==================== 辅助查询方法实现 ====================

ElementType MagneticScalar2DIntegrator::getElementType() const {
    return element_type_;
}

int MagneticScalar2DIntegrator::getNodeCount() const {
    return shape_func_->getNodeCount();
}

int MagneticScalar2DIntegrator::getGaussOrder() const {
    // 根据单元类型选择保证精度要求的最低积分阶数
    switch (element_type_) {
        case ElementType::TRI3:
            // 三角形3点积分方案：精确到二次多项式
            // 积分点位于三边中点，每点权重=1/6
            // 对于线性TRI3单元，3点积分可精确计算刚度矩阵（被积函数为常数）
            return 3;

        case ElementType::QUAD4:
            // 四边形4点（2×2 Gauss-Legendre）积分方案
            // 积分点位置: (±1/√3, ±1/√3)，每点权重=1.0
            // 精确到三次多项式（双线性单元的刚度核函数为双线性，完全精确）
            return 4;

        default:
            // 不应到达此处（构造函数已校验合法性），安全回退
            FEEM_ERROR("getGaussOrder: 未知的单元类型={}, 返回默认值1",
                       static_cast<int>(element_type_));
            return 1;
    }
}

} // namespace numeric
