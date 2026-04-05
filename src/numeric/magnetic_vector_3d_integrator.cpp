/**
 * @file magnetic_vector_3d_integrator.cpp
 * @brief 数值计算层 - 三维静磁场/瞬态磁场矢量磁位(A)积分器实现
 * @details 基于Nedelec第一类一阶H(curl)协调棱边元，实现三维矢量磁位A的
 *          单元级有限元矩阵组装。包含旋度刚度矩阵、质量矩阵（涡流项）、
 *          阻尼矩阵和源项向量的高斯积分计算。
 *
 * 数学公式汇总：
 *
 * 【旋度刚度矩阵 K_curl】（静态/瞬态共用）：
 *   K_ij = ∫_Ω (∇×N_i)^T · ν · (∇×N_j) dΩ
 *        = Σ_q w_q |detJ_q| ν_q (∇_x × N_i)_q · (∇_x × N_j)_q
 *   其中 ν = 1/μ 为磁阻率，(∇_x × N) 通过 Piola 旋度映射获取
 *
 * 【质量矩阵 M_A】（仅瞬态模式，涡流项 σ∂A/∂t）：
 *   M_ij = ∫_Ω N_i^T · σ · N_j dΩ
 *        = Σ_q w_q |detJ_q| σ_q W_phys_i · W_phys_j
 *   其中 W_phys = J^{-T} W_ref / detJ （Piola 变换）
 *
 * 【阻尼矩阵 C】（等同于质量矩阵乘以电导率系数）：
 *   C_ij = ∫_Ω σ N_i · N_j dΩ = M_ij （当前实现）
 *
 * 【源项向量 F_J】（当前为零源项）：
 *   F_i = ∫_Ω N_i · J_s dΩ = 0 （J_s = 0）
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#include "magnetic_vector_3d_integrator.hpp"

#include "shape_function_factory.hpp"

#include <cmath>
#include <algorithm>

namespace numeric {

// ==================== 构造函数与初始化 ====================

MagneticVector3DIntegrator::MagneticVector3DIntegrator(ElementType element_type)
    : element_type_(element_type)
    , gauss_order_(0)
{
    // 步骤1：验证单元类型必须是合法的3D Nedelec棱边元
    if (!isValid3DNedelecType(element_type)) {
        FEEM_ERROR("MagneticVector3DIntegrator: 不支持的单元类型={}, "
                   "仅支持TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE",
                   static_cast<int>(element_type));
        return;
    }

    // 步骤2：创建NedelecElement<3>形函数实例（三维棱边元）
    nedelec_func_ = std::make_unique<NedelecElement<3>>(element_type);

    // 步骤3：根据单元类型确定高斯积分阶数
    gauss_order_ = determineGaussOrder(element_type);

    int edge_count = nedelec_func_->getNodeCount();

    FEEM_INFO("MagneticVector3DIntegrator构造完成: 类型={}, 棱边数(DOF)={}, "
              "高斯积分阶数={}",
              static_cast<int>(element_type), edge_count, gauss_order_);
}

// ==================== 静态工具方法 ====================

bool MagneticVector3DIntegrator::isValid3DNedelecType(ElementType type)
{
    // 仅接受四种三维Nedelec棱边元类型
    return type == ElementType::TET4_EDGE
        || type == ElementType::HEX8_EDGE
        || type == ElementType::PRISM6_EDGE
        || type == ElementType::PYRAMID5_EDGE;
}

int MagneticVector3DIntegrator::determineGaussOrder(ElementType type)
{
    // 各3D Nedelec单元类型对应的高斯积分阶数映射
    // 积分阶数的选择保证对线性/二次多项式的精确积分
    switch (type) {
    case ElementType::TET4_EDGE:
        return 4;      // 4点Hammer积分方案（精确到二次多项式）

    case ElementType::HEX8_EDGE:
        return 8;      // 2×2×2 Gauss-Legendre积分（精确到三次多项式）

    case ElementType::PRISM6_EDGE:
        return 6;      // 三角形3点×ζ方向2点张量积（精确到二次多项式）

    case ElementType::PYRAMID5_EDGE:
        return 5;      // 1中心+4角点工程近似方案

    default:
        FEEM_ERROR("determineGaussOrder: 未知的3D Nedelec单元类型={}, 返回默认阶数0",
                   static_cast<int>(type));
        return 0;
    }
}

// ==================== 核心计算方法：computeAllMatrices ====================

ElementMatrices MagneticVector3DIntegrator::computeAllMatrices(
    const Eigen::MatrixXd& node_coords
) const {
    ElementMatrices result;

    // 步骤1：始终计算旋度刚度矩阵K（静态和瞬态共用）
    result.K = computeStiffnessMatrix(node_coords);

    // 步骤2：根据模式决定是否计算额外矩阵
    if (is_transient_) {
        // 瞬态模式：计算质量矩阵M（涡流σ∂A/∂t项）
        result.M = computeMassMatrix(node_coords);

        // 阻尼矩阵C在当前A-φ公式中与M使用相同核函数
        result.C = computeDampingMatrix(node_coords);
    } else {
        // 静态模式：M和C置为空矩阵（无物理意义）
        int n_edges = nedelec_func_->getNodeCount();
        result.M = Eigen::MatrixXd::Zero(n_edges, n_edges);
        result.C = Eigen::MatrixXd::Zero(n_edges, n_edges);
    }

    // 步骤3：源项向量F（当前为零源项J_s=0）
    result.F = computeSourceVector(node_coords);

    FEEM_DEBUG("computeAllMatrices完成: K维度={}×{}, M维度={}×{}, F维度={}",
               result.K.rows(), result.K.cols(),
               result.M.rows(), result.M.cols(),
               result.F.size());

    return result;
}

// ==================== 核心计算方法：computeStiffnessMatrix（旋度刚度矩阵）====================

Eigen::MatrixXd MagneticVector3DIntegrator::computeStiffnessMatrix(
    const Eigen::MatrixXd& node_coords
) const {
    // 获取单元自由度数（即棱边数量）
    int n_edges = nedelec_func_->getNodeCount();

    // 初始化刚度矩阵为零矩阵（n_edges × n_edges）
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_edges, n_edges);

    // 获取高斯积分点和权重
    // 注意：GaussQuadrature::getPoints 使用的是Lagrange单元类型枚举，
    // 但Nedelec单元与同几何形状的Lagrange单元共享相同的参考域和积分方案
    auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeStiffnessMatrix: 获取高斯积分点失败，返回零矩阵");
        return K;
    }

    // 遍历所有高斯积分点进行数值积分
    for (const auto& gp : gauss_points) {
        // 构造参考域局部坐标点
        LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 计算雅可比矩阵及其行列式（用于物理域映射和体积元变换）
        JacobianResult jaco_result = nedelec_func_->calcJacobian(xi, node_coords);

        double abs_det_j = std::abs(jaco_result.det_j);

        // 安全检查：雅可比行列式不能接近零（否则单元严重畸变）
        if (abs_det_j < 1e-15) {
            FEEM_WARN("computeStiffnessMatrix: 积分点({:.4f},{:.4f},{:.4f})处"
                      "雅可比行列式接近零({}), 跳过该点",
                      gp.coords[0], gp.coords[1], gp.coords[2],
                      jaco_result.det_j);
            continue;
        }

        // 获取有效磁阻率 ν = 1/μ（支持非线性材料模型）
        // 通过getEffectiveMu()可覆写以支持B-H曲线等非线性特性
        // 使用积分点的实际物理坐标进行材料参数查询（支持空间变化材料）
        Eigen::Vector3d physical_pos = computePhysicalPosition(xi, node_coords);
        double mu_eff = getEffectiveMu(physical_pos);
        double nu_effective = 1.0 / mu_eff;  // ν = 1/μ，单位：m/H

        // 性能优化：预先计算所有棱边的物理域旋度向量
        // 避免在内层循环中重复调用calcPhysicalEdgeFunction
        std::vector<Eigen::Vector3d> curl_physical_list(n_edges);

        for (int i = 0; i < n_edges; ++i) {
            EdgeCurlResult curl_result =
                nedelec_func_->calcPhysicalEdgeFunction(i, xi, node_coords);
            curl_physical_list[i] = curl_result.curl_physical;
        }

        // 组装刚度矩阵：遍历所有棱边对 (i, j)
        // K[i][j] += w_q · |detJ| · ν · (∇×N_i)_phys · (∇×N_j)_phys
        //
        // 数学说明：
        // - curl_physical_list[i] 已包含 (1/detJ)·J·curl_ref 的完整Piola旋度映射
        // - 内积运算 curl_i^T · curl_j 给出两个旋度向量的点积（标量）
        // - 对于各向同性材料，ν为标量系数；各向异性时可扩展为张量
        double weight_factor = gp.weight * abs_det_j * nu_effective;

        for (int i = 0; i < n_edges; ++i) {
            const Eigen::Vector3d& curl_i = curl_physical_list[i];

            // 利用对称性：K[i][j] = K[j][i]，仅需计算上三角并复制到下三角
            for (int j = i; j < n_edges; ++j) {
                const Eigen::Vector3d& curl_j = curl_physical_list[j];

                // 旋度向量的内积：(∇×N_i)^T · (∇×N_j) = 标量
                double curl_dot_product = curl_i.dot(curl_j);

                // 累加当前积分点的贡献
                K(i, j) += weight_factor * curl_dot_product;

                // 利用对称性填充下三角（避免重复计算）
                if (j != i) {
                    K(j, i) = K(i, j);
                }
            }
        }
    }

    FEEM_DEBUG("computeStiffnessMatrix完成: 维度={}×{}, "
              "Frobenius范数={:.6e}",
              n_edges, n_edges, K.norm());

    return K;
}

// ==================== 核心计算方法：computeMassMatrix（质量矩阵/涡流项）====================

Eigen::MatrixXd MagneticVector3DIntegrator::computeMassMatrix(
    const Eigen::MatrixXd& node_coords
) const {
    // 获取单元自由度数（棱边数量）
    int n_edges = nedelec_func_->getNodeCount();

    // 初始化质量矩阵为零矩阵
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_edges, n_edges);

    // 获取高斯积分点
    auto gauss_points = GaussQuadrature::getPoints(element_type_, gauss_order_);

    if (gauss_points.empty()) {
        FEEM_ERROR("computeMassMatrix: 获取高斯积分点失败，返回零矩阵");
        return M;
    }

    // 遍历所有高斯积分点
    for (const auto& gp : gauss_points) {
        // 构造参考域局部坐标点
        LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);

        // 计算雅可比矩阵结果（供Piola变换使用）
        JacobianResult jaco_result = nedelec_func_->calcJacobian(xi, node_coords);

        double abs_det_j = std::abs(jaco_result.det_j);

        if (abs_det_j < 1e-15) {
            FEEM_WARN("computeMassMatrix: 积分点({:.4f},{:.4f},{:.4f})处"
                      "雅可比行列式接近零({}), 跳过该点",
                      gp.coords[0], gp.coords[1], gp.coords[2],
                      jaco_result.det_j);
            continue;
        }

        // 获取有效电导率 σ（支持非线性/温度依赖模型）
        // 使用积分点的实际物理坐标进行材料参数查询（支持空间变化材料）
        Eigen::Vector3d physical_pos_sigma = computePhysicalPosition(xi, node_coords);
        double sigma_eff = getEffectiveSigma(physical_pos_sigma);

        // 电导率为零时（绝缘体区域），该积分点对质量矩阵无贡献
        if (sigma_eff <= 0.0) {
            continue;
        }

        // 性能优化：预先计算所有棱边的物理域基函数值（经Piola变换）
        // Piola变换公式：W_phys = J^{-T} · W_ref / detJ
        std::vector<Eigen::Vector3d> edge_phys_list(n_edges);

        for (int i = 0; i < n_edges; ++i) {
            edge_phys_list[i] = mapEdgeFunctionToPhysical(
                i, xi, node_coords, jaco_result);
        }

        // 组装质量矩阵：遍历所有棱边对 (i, j)
        // M[i][j] += w_q · |detJ| · σ · W_phys_i · W_phys_j
        //
        // 完整展开形式（代入Piola变换）：
        // M_ij = Σ_q w_q · σ_q / |detJ_q| · (W_ref_i^T · J^{-1} · J^{-T} · W_ref_j)
        //
        // 物理意义：
        // - 质量矩阵关联时间导数项 σ∂A/∂t，描述涡流效应的能量耗散
        // - 对于导体区域（如铜σ≈5.8×10⁷ S/m），此项显著
        // - 对于绝缘体区域（σ=0），此项为零
        double weight_factor = gp.weight * abs_det_j * sigma_eff;

        for (int i = 0; i < n_edges; ++i) {
            const Eigen::Vector3d& W_phys_i = edge_phys_list[i];

            for (int j = i; j < n_edges; ++j) {
                const Eigen::Vector3d& W_phys_j = edge_phys_list[j];

                // 矢量基函数的内积：W_phys_i^T · W_phys_j = 标量
                double edge_dot_product = W_phys_i.dot(W_phys_j);

                // 累加当前积分点的贡献
                M(i, j) += weight_factor * edge_dot_product;

                // 利用对称性填充下三角
                if (j != i) {
                    M(j, i) = M(i, j);
                }
            }
        }
    }

    FEEM_DEBUG("computeMassMatrix完成: 维度={}×{} Frobenius范数={:.6e}",
              n_edges, n_edges, M.norm());

    return M;
}

// ==================== 核心计算方法：computeDampingMatrix（阻尼矩阵）====================

Eigen::MatrixXd MagneticVector3DIntegrator::computeDampingMatrix(
    const Eigen::MatrixXd& node_coords
) const {
    // 在A-φ公式的瞬态涡流方程中：
    // ∇×(ν∇×A) + σ∂A/∂t = J_s
    // 阻尼矩阵C和质量矩阵M使用相同的核函数（σ N_i · N_j）
    // 因为方程中只有一个含σ的时间导数项
    //
    // 注：在某些扩展公式体系（如A-V-Ω或含速度效应）中，
    // C可能与M有不同的物理系数，此时需分别实现

    return computeMassMatrix(node_coords);
}

// ==================== 核心计算方法：computeSourceVector（源项向量）====================

Eigen::VectorXd MagneticVector3DIntegrator::computeSourceVector(
    const Eigen::MatrixXd& node_coords
) const {
    // 获取单元自由度数
    int n_edges = nedelec_func_->getNodeCount();

    // 当前实现：返回零源项向量（J_s = 0）
    // 适用场景：无外加电流激励的磁场分析
    //
    // 理论公式（供后续扩展参考）：
    // F_J[i] = ∫_Ω N_edge_i(x) · J_s(x) dΩ
    //
    // 后续可支持的激励类型：
    // 1. 均匀电流密度：J_s = 常矢量 → 解析积分
    // 2. 线圈电流：线电流奇异性 → 分布式源近似
    // 3. 外电路耦合：从电路方程获取等效面电流密度
    // 4. 永磁体等效：通过表面磁荷密度建模
    //
    (void)node_coords;  // 参数保留供未来扩展使用

    return Eigen::VectorXd::Zero(n_edges);
}

// ==================== 查询方法 ====================

int MagneticVector3DIntegrator::getEdgeCount() const
{
    if (!nedelec_func_) {
        FEEM_WARN("getEdgeCount: Nedelec形函数对象未初始化");
        return 0;
    }
    return nedelec_func_->getNodeCount();
}

ElementType MagneticVector3DIntegrator::getElementType() const
{
    return element_type_;
}

// ==================== 内部工具方法：Piola变换映射 ====================

Eigen::Vector3d MagneticVector3DIntegrator::mapEdgeFunctionToPhysical(
    int edge_index,
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords,
    const JacobianResult& jaco_result
) const {
    // node_coords保留参数以保持接口一致性（jaco_result由调用方预计算）
    (void)node_coords;

    /**
     * Piola变换详解（H(curl)协调变换）：
     *
     * 对于矢量值形函数，从参考域到物理域的映射不能简单使用
     * 等参变换（那是针对标量值的），而必须使用Piola变换以保证
     * 切向分量在单元边界上的连续性（H(curl)协调性条件）。
     *
     * 变换公式：
     *   W_phys(x) = J^{-T}(ξ) · W_ref(ξ) / det(J)
     *
     * 各分量含义：
     *   - W_ref(ξ)：参考域中的Nedelec基函数矢量值（通过evalEdgeFunction获取）
     *   - J = ∂x/∂ξ：雅可比矩阵（3×3），将参考域微分映射到物理域
     *   - J^{-T} = (J^{-1})^T：雅可比逆矩阵的转置
     *   - det(J)：雅可比行列式
     *
     * 为什么需要Piola变换？
     *   标量形函数的等参变换保证函数值的连续性（C⁰连续），
     *   但对于矢量场，我们需要的不是函数值连续，而是切向分量连续
     *   （因为∇×A的边界条件涉及切向分量）。Piola变换恰好满足这一要求：
     *   t^T · W_phys 在跨越单元边界时保持连续，其中t为边界切向向量。
     */

    // 步骤1：获取参考域基函数矢量值
    Eigen::Vector3d W_ref = nedelec_func_->evalEdgeFunction(edge_index, xi);

    // 步骤2：安全检查雅可比行列式
    if (std::abs(jaco_result.det_j) < 1e-15) {
        FEEM_WARN("mapEdgeFunctionToPhysical: 雅可比行列式接近零({}), "
                  "返回零向量", jaco_result.det_j);
        return Eigen::Vector3d::Zero();
    }

    // 步骤3：应用Piola变换公式
    // W_phys = J^{-T} · W_ref / detJ
    // 等价于：W_phys = (J^{-1})^T · W_ref / detJ
    Eigen::Vector3d W_phys = (jaco_result.inv_jacobian.transpose() * W_ref)
                             / jaco_result.det_j;

    return W_phys;
}

// ==================== 辅助方法：computePhysicalPosition ====================

Eigen::Vector3d MagneticVector3DIntegrator::computePhysicalPosition(
    const LocalPoint& xi,
    const Eigen::MatrixXd& node_coords
) const {
    // 根据Nedelec单元的几何类型，确定对应的Lagrange标量单元类型字符串
    // 用于计算物理坐标的等参映射
    std::string lagrange_type_str;

    switch (element_type_) {
    case ElementType::TET4_EDGE:
        lagrange_type_str = "TET4";
        break;
    case ElementType::HEX8_EDGE:
        lagrange_type_str = "HEX8";
        break;
    case ElementType::PRISM6_EDGE:
        lagrange_type_str = "PRISM6";
        break;
    case ElementType::PYRAMID5_EDGE:
        lagrange_type_str = "PYRAMID5";
        break;
    default:
        FEEM_WARN("computePhysicalPosition: 未知的Nedelec单元类型={}, 使用质心近似",
                   static_cast<int>(element_type_));
        // 降级方案：返回所有节点的平均位置
        return node_coords.rowwise().mean().head<3>();
    }

    // 创建临时Lagrange形函数对象（仅用于坐标映射）
    auto shape_func = ShapeFunctionFactory::create(lagrange_type_str);
    if (!shape_func) {
        FEEM_ERROR("computePhysicalPosition: 无法创建Lagrange形函数对象, type={}",
                   lagrange_type_str);
        return Eigen::Vector3d::Zero();
    }

    // 计算参考域坐标点处的标量形函数值
    Eigen::VectorXd N_vals = shape_func->evalN(xi);

    // 物理坐标 = 节点坐标矩阵 × 形函数值向量
    // node_coords 维度: (dim+1) × n_nodes 或 dim × n_nodes
    // N_vals 维度: n_nodes × 1
    // 结果维度: (dim+1) × 1 或 dim × 1，取前3维作为3D坐标
    Eigen::Vector3d physical_pos = (node_coords * N_vals).head<3>();

    return physical_pos;
}

} // namespace numeric
