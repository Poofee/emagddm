/**
 * @file magnetic_vector_3d_integrator.hpp
 * @brief 数值计算层 - 三维静磁场/瞬态磁场矢量磁位(A)积分器
 * @details 基于Nedelec第一类一阶H(curl)协调棱边元，实现三维矢量磁位A的
 *          单元级有限元矩阵组装。支持静态磁场和瞬态涡流场两种分析模式。
 *
 *          物理方程与弱形式：
 *
 *          【静态磁场】（忽略位移电流，∇×H = J_s）：
 *          控制方程：∇×(ν∇×A) = J_s
 *          其中 ν = 1/μ 为磁阻率（m/H），J_s 为源电流密度（A/m²）
 *          弱形式刚度矩阵：
 *            K_curl[i][j] = ∫_Ω (∇×N_i)^T · ν · (∇×N_j) dΩ
 *
 *          【瞬态涡流场】（A-V公式化，仅A分量）：
 *          控制方程：∇×(ν∇×A) + σ∂A/∂t = J_s
 *          弱形式矩阵：
 *            K_curl[i][j] = ∫_Ω (∇×N_i)^T · ν · (∇×N_j) dΩ   （旋度刚度）
 *            M_A[i][j]    = ∫_Ω N_i^T · σ · N_j dΩ              （质量矩阵/涡流项）
 *            F_J[i]       = ∫_Ω N_i · J_s dΩ                    （源项向量，当前为零）
 *
 *          支持的单元类型（仅3D Nedelec棱边元）：
 *          - TET4_EDGE：四面体，6条棱边，6个自由度
 *          - HEX8_EDGE：六面体，12条棱边，12个自由度
 *          - PRISM6_EDGE：三棱柱，9条棱边，9个自由度
 *          - PYRAMID5_EDGE：金字塔，8条棱边，8个自由度
 *
 *          H(curl)协调性说明：
 *          Nedelec棱边元的基函数在共享棱边上具有切向连续性，
 *          保证∇×A在跨越单元边界时的正确性，由形函数库自动保证。
 *
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include "em_element_integrator_base.hpp"
#include "nedelec_element.hpp"
#include "gauss_quadrature.hpp"
#include "logger_factory.hpp"

#include <memory>
#include <vector>

namespace numeric {

/**
 * @class MagneticVector3DIntegrator
 * @brief 三维静磁场/瞬态磁场矢量磁位(A)积分器
 * @details 继承EMElementIntegratorBase基类，使用NedelecElement<3>作为形函数，
 *          实现基于H(curl)弱形式的单元矩阵高斯积分计算。
 *
 *          核心特性：
 *          - 形函数类型：Nedelec第一类一阶棱边元（矢量值），非Lagrange节点元（标量值）
 *          - 自由度含义：棱边环流（edge circulation），非节点值
 *          - 弱形式算子：curl-curl算子 ∇×(ν∇×)，非梯度散度算子 ∇·(μ∇)
 *          - 矩阵对称性：K_curl对于各向同性材料仍为对称正定（半正定）矩阵
 *
 *          物理域映射：
 *          - 旋度物理域映射：curl_phys = (1/detJ) · J · curl_ref
 *          - 基函数Piola变换：W_phys = J^{-T} · W_ref / detJ
 *          均通过calcPhysicalEdgeFunction()和雅可比矩阵自动处理
 *
 *          应用场景：
 *          - 三维静磁场分析（电机、变压器、电磁铁等）
 *          - 瞬态涡流场分析（感应加热、涡流制动、电磁屏蔽等）
 *          - 低频电磁兼容（EMC）仿真
 *
 * @note 与标量磁位积分器的本质区别在于使用矢量值Nedelec基函数和curl-curl算子
 * @see EMElementIntegratorBase 积分器抽象基类
 * @see NedelecElement<3> 三维Nedelec棱边元形函数
 * @see GaussQuadrature 高斯积分类
 */
class MagneticVector3DIntegrator : public EMElementIntegratorBase {
public:
    /**
     * @brief 构造函数，根据单元类型初始化Nedelec形函数对象
     * @param element_type 单元类型枚举值，必须是3D Nedelec棱边元类型
     *
     * @details 支持的类型及对应自由度数：
     *          - TET4_EDGE: 6条棱边（6 DOF）
     *          - HEX8_EDGE: 12条棱边（12 DOF）
     *          - PRISM6_EDGE: 9条棱边（9 DOF）
     *          - PYRAMID5_EDGE: 8条棱边（8 DOF）
     *
     *          内部操作：
     *          1. 验证单元类型的合法性（必须为3D _EDGE类型）
     *          2. 创建NedelecElement<3>形函数实例
     *          3. 根据单元类型设置对应的高斯积分阶数
     *
     * @warning 非法单元类型（如Lagrange节点元或2D Nedelec元）将输出ERROR日志
     */
    explicit MagneticVector3DIntegrator(ElementType element_type);

    /**
     * @brief 虚析构函数，释放Nedelec形函数对象资源
     */
    ~MagneticVector3DIntegrator() override = default;

    // ========== 纯虚方法实现（EMElementIntegratorBase接口）==========

    /**
     * @brief 计算单元的全部矩阵和源项向量（核心入口方法）
     * @param node_coords 节点坐标矩阵，维度为 3 × node_count，
     *                    每列对应一个节点的物理坐标 (x, y, z)^T
     * @return ElementMatrices 包含K、M、C、F的结构体
     *
     * @details 根据当前模式（静态/瞬态）计算相应的矩阵集合：
     *          - 静态模式（is_transient_=false）：计算K（旋度刚度）、F（零源项）
     *          - 瞬态模式（is_transient_=true）：计算全部K、M、C、F
     *
     *          内部调用链：
     *          computeAllMatrices → computeStiffnessMatrix (K)
     *                           → computeMassMatrix (M, 仅瞬态)
     *                           → computeDampingMatrix (C, 仅瞬态)
     *                           → computeSourceVector (F)
     *
     * @note K矩阵始终计算；M/C矩阵仅在瞬态模式下有意义
     */
    ElementMatrices computeAllMatrices(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 计算旋度刚度矩阵 K_curl（curl-curl项）
     * @param node_coords 节点坐标矩阵，维度为 3 × node_count
     * @return Eigen::MatrixXd 旋度刚度矩阵，维度 n_edges × n_edges
     *
     * @details 数值积分公式：
     *          K_curl[i][j] = Σ_{q=1}^{n_gp} w_q · |detJ_q| · ν_eff ·
     *                         (∇_x × N_i)_q^T · (∇_x × N_j)_q
     *
     *          其中：
     *          - (∇_x × N)_q：第q个积分点处基函数的物理域旋度向量，
     *            通过 nedelec_func_->calcPhysicalEdgeFunction() 获取
     *          - ν_eff = 1/μ_eff：有效磁阻率（m/H），
     *            通过 getEffectiveMu() 获取并取倒数
     *          - detJ_q：第q个积分点处的雅可比行列式
     *          - w_q：第q个积分点的权重
     *
     *          对于各向同性线性材料，K_curl为对称半正定矩阵。
     *          对于非线性材料（B-H曲线），ν随位置变化，需逐积分点求值。
     *
     * @note 此方法是静态和瞬态模式的共用核心计算
     */
    Eigen::MatrixXd computeStiffnessMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 计算质量矩阵 M_A（涡流项，对应σ∂A/∂t）
     * @param node_coords 节点坐标矩阵，维度为 3 × node_count
     * @return Eigen::MatrixXd 质量矩阵，维度 n_edges × n_edges
     *
     * @details 数值积分公式（矢量值基函数的内积）：
     *          M_A[i][j] = Σ_{q=1}^{n_gp} w_q · |detJ_q| · σ_eff ·
     *                      W_phys_i(ξ_q)^T · W_phys_j(ξ_q)
     *
     *          其中W_phys为经Piola变换后的物理域基函数值：
     *          W_phys(x) = J^{-T}(ξ) · W_ref(ξ) / detJ(ξ)
     *
     *          Piola变换保证H(curl)协调性：
     *          切向分量在单元边界上连续，确保∇×A的全局正确性
     *
     *          完整展开形式：
     *          M_ij = Σ_q w_q · σ_q / |detJ_q| ·
     *                 [W_ref_i^T · J^{-1} · J^{-T} · W_ref_j]
     *
     * @note 仅在瞬态模式下此矩阵有物理意义（电导率σ > 0的区域）
     * @warning 对于绝缘体区域（σ=0），此矩阵为零矩阵
     */
    Eigen::MatrixXd computeMassMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 计算阻尼矩阵 C（涡流耗散项）
     * @param node_coords 节点坐标矩阵，维度为 3 × node_count
     * @return Eigen::MatrixXd 阻尼矩阵，维度 n_edges × n_edges
     *
     * @details 在A-φ公式体系的瞬态涡流方程中，阻尼矩阵与质量矩阵
     *          使用相同的核函数（σ N_i · N_j），但可能关联不同的时间导数阶数
     *          或频率因子。当前实现中，阻尼矩阵等同于质量矩阵乘以电导率系数。
     *
     *          C_ij = ∫_Ω σ · N_edge_i · N_edge_j dΩ
     *
     * @note 当前实现返回与computeMassMatrix相同的结果（σ加权的一致质量矩阵）
     * @warning 仅在瞬态模式下调用
     */
    Eigen::MatrixXd computeDampingMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override;

    /**
     * @brief 计算源项向量 F_J（外部电流激励投影）
     * @param node_coords 节点坐标矩阵，维度为 3 × node_count
     * @return Eigen::VectorXd 源项向量，维度 n_edges × 1
     *
     * @details 源项向量的理论公式：
     *          F_J[i] = ∫_Ω N_edge_i · J_s dΩ
     *          其中 J_s 为外加源电流密度（A/m²）
     *
     *          当前实现返回零向量（J_s = 0），适用于无源激励的分析场景。
     *          后续扩展可支持：
     *          - 均匀电流密度：J_s = 常矢量
     *          - 线圈激励：通过线电流奇异性处理
     *          - 外电路耦合：从电路方程获取等效J_s
     *
     * @note 返回零向量不影响刚度矩阵的正确性，仅影响右端项
     */
    Eigen::VectorXd computeSourceVector(
        const Eigen::MatrixXd& node_coords
    ) const override;

    // ========== 查询方法 ==========

    /**
     * @brief 获取单元的棱边数量（即自由度数量）
     * @return int 棱边数量：TET4_EDGE=6, HEX8_EDGE=12,
     *             PRISM6_EDGE=9, PYRAMID5_EDGE=8
     *
     * @note 对于Nedelec单元，getNodeCount()返回的是棱边数而非节点数
     */
    int getEdgeCount() const;

    /**
     * @brief 获取当前使用的单元类型
     * @return ElementType 当前Nedelec单元的类型枚举值
     */
    ElementType getElementType() const;

private:
    std::unique_ptr<NedelecElement<3>> nedelec_func_;  ///< Nedelec三维棱边形函数对象
    ElementType element_type_;                          ///< 单元类型枚举值
    int gauss_order_;                                   ///< 高斯积分阶数

    // ========== 内部工具方法 ==========

    /**
     * @brief 根据单元类型确定对应的高斯积分阶数
     * @param type 单元类型枚举值
     * @return int 对应的高斯积分阶数
     *
     * @details 各单元类型与积分阶数的映射关系：
     *          - TET4_EDGE → order=4（4点Hammer积分，精确到二次多项式）
     *          - HEX8_EDGE → order=8（2×2×2 Gauss-Legendre积分）
     *          - PRISM6_EDGE → order=6（三角形3点×ζ方向2点张量积）
     *          - PYRAMID5_EDGE → order=5（1中心+4角点工程近似方案）
     *
     * @note 积分阶数的选择保证对线性多项式的精确积分
     */
    static int determineGaussOrder(ElementType type);

    /**
     * @brief 检查单元类型是否为合法的3D Nedelec棱边元
     * @param type 待检查的单元类型枚举值
     * @return true 类型合法（TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE）
     * @return false 类型非法
     */
    static bool isValid3DNedelecType(ElementType type);

    /**
     * @brief 计算参考域积分点对应的物理域坐标
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵（dim × node_count 或 (dim+1) × node_count）
     * @return Eigen::Vector3d 物理域三维坐标向量
     *
     * @details 对于Nedelec棱边元，由于evalN()返回零向量（无标量形函数概念），
     *          本方法通过构造等参映射的标量形函数来计算物理位置：
     *          - 使用与单元几何类型匹配的Lagrange标量形函数
     *          - 物理坐标 = Σ N_i(ξ) · X_i （节点坐标的加权和）
     *
     * @note 此方法主要用于非线性材料参数查询时的位置信息获取，
     *       保证 getEffectiveMu(position) 和 getEffectiveSigma(position)
     *       能够接收到正确的空间位置坐标
     */
    Eigen::Vector3d computePhysicalPosition(
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords
    ) const;

    /**
     * @brief 将参考域Nedelec基函数值通过Piola变换映射到物理域
     * @param edge_index 棱边索引
     * @param xi 参考域局部坐标点
     * @param node_coords 节点物理坐标矩阵
     * @param jaco_result 已计算的雅可比结果（避免重复计算）
     * @return Eigen::Vector3d 物理域基函数矢量值
     *
     * @details Piola变换公式（H(curl)协调变换）：
     *          W_phys(x) = J^{-T}(ξ) · W_ref(ξ) / det(J)
     *
     *          其中：
     *          - W_ref：参考域基函数值（通过evalEdgeFunction获取）
     *          - J：雅可比矩阵（3×3）
     *          - J^{-T}：雅可比逆矩阵的转置
     *          - det(J)：雅可比行列式
     *
     *          此变换保证基函数的切向分量在单元边界上连续，
     *          是H(curl)有限元方法的核心数学基础之一
     */
    Eigen::Vector3d mapEdgeFunctionToPhysical(
        int edge_index,
        const LocalPoint& xi,
        const Eigen::MatrixXd& node_coords,
        const JacobianResult& jaco_result
    ) const;
};

} // namespace numeric
