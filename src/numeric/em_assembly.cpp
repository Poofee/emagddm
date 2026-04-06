/**
 * @file em_assembly.cpp
 * @brief 数值计算层 - 电磁场有限元全局组装器实现
 * @details 实现电磁场有限元分析的全局系统矩阵/向量组装功能，
 *          负责将单元级矩阵（K_e, M_e, C_e）和源项向量（F_e）通过scatter操作
 *          组装为全局刚度矩阵（K）、质量矩阵（M）、阻尼矩阵（C）和右端项向量（F）。
 *
 *          核心组装流程：
 *          1. 预估非零元素数量并预分配COO矩阵内存
 *          2. 遍历所有单元，调用EMElementIntegratorBase计算单元矩阵
 *          3. 通过Local2Global映射将单元矩阵scatter到全局COO矩阵
 *          4. 将COO格式转换为CSR格式以支持高效求解
 *          5. 统计装配信息并进行对称性检查
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "em_assembly.hpp"
#include "logger_factory.hpp"

#include <chrono>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace numeric {

// ==================== 内部测试用积分器类 ====================

/**
 * @class DummyEMIntegrator
 * @brief 测试用占位符积分器（用于验证组装逻辑正确性）
 * @details 在正式的LagrangeEMIntegrator和NedelecEMIntegrator实现完成前，
 *          提供一个返回单位矩阵或简单测试数据的临时实现。
 *          此类仅用于开发和测试阶段，后续将被具体积分器替代。
 */
class DummyEMIntegrator : public EMElementIntegratorBase {
public:
    /**
     * @brief 构造函数
     * @param dof_count 单元自由度数（决定返回矩阵的维度）
     */
    explicit DummyEMIntegrator(int dof_count)
        : dof_count_(dof_count)
    {
    }

    /**
     * @brief 计算全部单元矩阵和源项向量
     * @param node_coords 节点坐标矩阵（未使用）
     * @return ElementMatrices 包含测试数据的单元矩阵集合
     */
    ElementMatrices computeAllMatrices(
        const Eigen::MatrixXd& node_coords
    ) const override {
        (void)node_coords;

        ElementMatrices mats;
        // 返回单位矩阵作为刚度矩阵（用于验证scatter逻辑）
        mats.K = Eigen::MatrixXd::Identity(dof_count_, dof_count_);
        // 返回对角质量矩阵（对角线元素=材料电导率σ）
        mats.M = Eigen::MatrixXd::Identity(dof_count_, dof_count_) * material_props_.sigma;
        // 返回零阻尼矩阵（简化模型）
        mats.C = Eigen::MatrixXd::Zero(dof_count_, dof_count_);
        // 返回零源项向量（无外部激励）
        mats.F = Eigen::VectorXd::Zero(dof_count_);

        FEEM_DEBUG("DummyEMIntegrator::computeAllMatrices: 返回{}×{}测试矩阵",
                   dof_count_, dof_count_);
        return mats;
    }

    /**
     * @brief 计算单元刚度矩阵
     * @param node_coords 节点坐标矩阵（未使用）
     * @return Eigen::MatrixXd 单位刚度矩阵
     */
    Eigen::MatrixXd computeStiffnessMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override {
        (void)node_coords;
        return Eigen::MatrixXd::Identity(dof_count_, dof_count_);
    }

    /**
     * @brief 计算单元质量矩阵
     * @param node_coords 节点坐标矩阵（未使用）
     * @return Eigen::MatrixXd 对角质量矩阵
     */
    Eigen::MatrixXd computeMassMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override {
        (void)node_coords;
        return Eigen::MatrixXd::Identity(dof_count_, dof_count_) * material_props_.sigma;
    }

    /**
     * @brief 计算单元阻尼矩阵
     * @param node_coords 节点坐标矩阵（未使用）
     * @return Eigen::VectorXd 零阻尼矩阵
     */
    Eigen::MatrixXd computeDampingMatrix(
        const Eigen::MatrixXd& node_coords
    ) const override {
        (void)node_coords;
        return Eigen::MatrixXd::Zero(dof_count_, dof_count_);
    }

    /**
     * @brief 计算单元源项向量
     * @param node_coords 节点坐标矩阵（未使用）
     * @return Eigen::VectorXd 零源项向量
     */
    Eigen::VectorXd computeSourceVector(
        const Eigen::MatrixXd& node_coords
    ) const override {
        (void)node_coords;
        return Eigen::VectorXd::Zero(dof_count_);
    }

private:
    int dof_count_;  ///< 单元自由度数
};

// ==================== EMAssembly 类实现 ====================

// ---------- 构造与析构 ----------

EMAssembly::EMAssembly()
    : K_csr_()
    , M_csr_()
    , C_csr_()
    , F_()
    , stats_()
    , num_threads_(1)
{
    clear();
}

// ---------- 核心组装接口 ----------

bool EMAssembly::assemble(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<fe_em::Local2Global>& elem_local_to_global,
    const std::vector<MaterialProperties>& materials,
    bool is_transient
) {
    // ========== Step 1: 输入参数校验 ==========

    FEEM_INFO("assemble: 开始全局矩阵组装流程");

    // 检查网格数据有效性
    if (mesh_data.nodes.empty()) {
        stats_.error_message = "错误：网格节点列表为空";
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    if (mesh_data.elements.empty()) {
        stats_.error_message = "错误：网格单元列表为空";
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    // 检查DOF映射表完整性
    if (elem_local_to_global.size() < mesh_data.elements.size()) {
        stats_.error_message = "错误：DOF映射表大小（" +
                              std::to_string(elem_local_to_global.size()) +
                              "）小于单元数量（" +
                              std::to_string(mesh_data.elements.size()) + "）";
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    // 检查材料参数有效性
    if (materials.empty()) {
        stats_.error_message = "错误：材料参数列表为空";
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    // 清空上次组装结果
    clear();

    // 计算（估算）总DOF数（用于COO矩阵维度）
    int total_dofs = 0;
    for (const auto& l2g : elem_local_to_global) {
        int max_idx = 0;
        for (int idx : l2g.indices) {
            if (idx > max_idx) {
                max_idx = idx;
            }
        }
        if (max_idx >= total_dofs) {
            total_dofs = max_idx + 1;
        }
    }

    if (total_dofs <= 0) {
        stats_.error_message = "错误：无效的总自由度数（<= 0）";
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    FEEM_INFO("assemble: 输入校验通过 - 单元数={}, 总DOF数={}, 瞬态模式={}",
              mesh_data.elements.size(), total_dofs,
              is_transient ? "开启" : "关闭");

    // ========== Step 2: 初始化COO矩阵和F向量 ==========

    // 估算非零元素数量（基于启发式算法）
    int estimated_nnz = estimateNNZ(mesh_data.elements, 15);

    FEEM_DEBUG("assemble: 预估非零元素数量NNZ={}", estimated_nnz);

    // 创建COO格式稀疏矩阵（预分配内存避免动态扩容）
    CooMatrix<double> K_coo(total_dofs, total_dofs, estimated_nnz);
    CooMatrix<double> M_coo;   // 瞬态场才使用
    CooMatrix<double> C_coo;   // 瞬态场才使用

    if (is_transient) {
        M_coo = CooMatrix<double>(total_dofs, total_dofs, estimated_nnz);
        C_coo = CooMatrix<double>(total_dofs, total_dofs, estimated_nnz);
    }

    // 初始化右端项源向量（全零）
    Eigen::VectorXd F = Eigen::VectorXd::Zero(total_dofs);

    // ========== Step 3: 单元遍历循环 ==========

    auto start_time = std::chrono::high_resolution_clock::now();

    int processed_elements = 0;
    int skipped_elements = 0;

    for (const auto& elem : mesh_data.elements) {
        // 3.1 获取该单元的Local2Global映射表
        const auto& l2g = elem_local_to_global[elem.id];

        // 3.2 提取单元节点坐标
        Eigen::MatrixXd node_coords;
        try {
            node_coords = buildNodeCoords(elem, mesh_data.nodes);
        } catch (const std::exception& e) {
            FEEM_WARN("assemble: 单元[{}]节点坐标提取失败: {} - 跳过此单元",
                     elem.id, e.what());
            skipped_elements++;
            continue;
        }

        // 3.3 获取材料参数
        if (elem.material_id < 0 || elem.material_id >= static_cast<int>(materials.size())) {
            FEEM_WARN("assemble: 单元[{}]的材料ID={}越界（有效范围[0, {}]）- 跳过此单元",
                     elem.id, elem.material_id, materials.size() - 1);
            skipped_elements++;
            continue;
        }
        const auto& material = materials[elem.material_id];

        // 3.4 创建积分器并计算单元矩阵
        auto integrator = createIntegrator(elem.elem_type, elem.dof_type);
        if (!integrator) {
            FEEM_WARN("assemble: 单元[{}]类型组合不支持（ElemType={}, DOFType={}）- 跳过",
                     elem.id,
                     static_cast<int>(elem.elem_type),
                     static_cast<int>(elem.dof_type));
            skipped_elements++;
            continue;
        }

        // 设置材料属性和瞬态模式
        integrator->setMaterialProperties(material);
        integrator->setTransientMode(is_transient);

        // 计算单元矩阵
        ElementMatrices elem_mats;
        try {
            if (is_transient) {
                elem_mats = integrator->computeAllMatrices(node_coords);
            } else {
                // 静态模式仅计算K和F（性能优化）
                elem_mats.K = integrator->computeStiffnessMatrix(node_coords);
                elem_mats.F = integrator->computeSourceVector(node_coords);
            }
        } catch (const std::exception& e) {
            FEEM_ERROR("assemble: 单元[{}]矩阵计算异常: {}", elem.id, e.what());
            skipped_elements++;
            continue;
        }

        // 3.5 Scatter单元矩阵到全局COO矩阵
        scatterElementMatrix(K_coo, elem_mats.K, l2g);

        if (is_transient) {
            scatterElementMatrix(M_coo, elem_mats.M, l2g);
            scatterElementMatrix(C_coo, elem_mats.C, l2g);
        }

        // 3.6 Scatter源项向量到全局F向量
        scatterElementVector(F, elem_mats.F, l2g);

        processed_elements++;

        // 定期输出进度日志（每1000个单元输出一次）
        if (processed_elements % 1000 == 0) {
            FEEM_DEBUG("assemble: 已处理 {}/{} 个单元",
                      processed_elements, mesh_data.elements.size());
        }
    }

    // ========== Step 4: COO→CSR转换 ==========

    FEEM_INFO("assemble: 单元处理完成 - 成功={}, 跳过={}",
              processed_elements, skipped_elements);

    try {
        // 初始化CSR矩阵尺寸（必须与COO矩阵匹配）
        K_csr_ = CsrMatrix<double>(total_dofs, total_dofs);

        // 将COO格式转换为CSR格式（排序+压缩）
        K_csr_.build_from_coo(K_coo);

        if (is_transient) {
            M_csr_ = CsrMatrix<double>(total_dofs, total_dofs);
            C_csr_ = CsrMatrix<double>(total_dofs, total_dofs);
            M_csr_.build_from_coo(M_coo);
            C_csr_.build_from_coo(C_coo);
        }

        // 保存右端项向量
        F_ = F;
    } catch (const std::exception& e) {
        stats_.error_message = "错误：COO→CSR转换失败 - " + std::string(e.what());
        FEEM_ERROR("assemble: {}", stats_.error_message);
        return false;
    }

    // ========== Step 5: 统计信息收集 ==========

    auto end_time = std::chrono::high_resolution_clock::now();
    double assembly_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time
    ).count();

    // 填充统计信息结构体
    stats_.assembly_time_ms = assembly_time_ms;
    stats_.total_elements = processed_elements;
    stats_.free_dofs = total_dofs;

    // 计算约束DOF数量（通过统计l2g中-1的数量估算）
    int constrained_count = 0;
    for (const auto& l2g : elem_local_to_global) {
        for (int idx : l2g.indices) {
            if (idx == -1) {
                constrained_count++;
            }
        }
    }
    // 去重处理（同一约束DOF可能被多个单元引用）
    stats_.constrained_dofs = std::min(constrained_count,
                                      static_cast<int>(mesh_data.nodes.size()) * 3);

    // 矩阵非零元素数量
    stats_.nnz_stiffness = K_csr_.nnz();
    stats_.nnz_mass = is_transient ? M_csr_.nnz() : 0;
    stats_.nnz_damping = is_transient ? C_csr_.nnz() : 0;

    // 对称性检查
    stats_.is_symmetric_k = checkSymmetry(K_csr_);
    stats_.is_symmetric_m = is_transient ? checkSymmetry(M_csr_) : false;
    stats_.is_symmetric_c = is_transient ? checkSymmetry(C_csr_) : false;

    // 输出汇总日志
    FEEM_INFO("assemble: 全局矩阵组装完成");
    FEEM_INFO("  - 装配耗时: {:.2f} ms", assembly_time_ms);
    FEEM_INFO("  - 处理单元: {} 个（跳过 {} 个）", processed_elements, skipped_elements);
    FEEM_INFO("  - 自由DOF: {}, 约束DOF: {}", total_dofs, stats_.constrained_dofs);
    FEEM_INFO("  - 刚度矩阵K: {}×{}, NNZ={}, 对称={}",
              K_csr_.rows(), K_csr_.cols(), stats_.nnz_stiffness,
              stats_.is_symmetric_k ? "是" : "否");

    if (is_transient) {
        FEEM_INFO("  - 质量矩阵M: {}×{}, NNZ={}, 对称={}",
                  M_csr_.rows(), M_csr_.cols(), stats_.nnz_mass,
                  stats_.is_symmetric_m ? "是" : "否");
        FEEM_INFO("  - 阻尼矩阵C: {}×{}, NNZ={}, 对称={}",
                  C_csr_.rows(), C_csr_.cols(), stats_.nnz_damping,
                  stats_.is_symmetric_c ? "是" : "否");
    }

    return true;
}

// ---------- 结果访问接口 ----------

const CsrMatrix<double>& EMAssembly::getStiffnessMatrix() const {
    return K_csr_;
}

const CsrMatrix<double>& EMAssembly::getMassMatrix() const {
    return M_csr_;
}

const CsrMatrix<double>& EMAssembly::getDampingMatrix() const {
    return C_csr_;
}

const Eigen::VectorXd& EMAssembly::getSourceVector() const {
    return F_;
}

AssemblyStats EMAssembly::getStats() const {
    return stats_;
}

// ---------- 资源管理 ----------

void EMAssembly::clear() {
    // 清空所有CSR矩阵
    K_csr_.clear();
    M_csr_.clear();
    C_csr_.clear();

    // 清空右端项向量（释放Eigen向量内存）
    F_ = Eigen::VectorXd();

    // 重置统计信息为默认状态
    stats_ = AssemblyStats();

    FEEM_DEBUG("clear: 组装器状态已重置");
}

// ---------- 并行控制接口（OpenMP预留）----------

void EMAssembly::setNumThreads(int num_threads) {
    if (num_threads <= 0) {
        FEEM_WARN("setNumThreads: 无效线程数={}，强制设置为1", num_threads);
        num_threads = 1;
    }
    num_threads_ = num_threads;
    FEEM_DEBUG("setNumThreads: 线程数设置为{}", num_threads_);
}

int EMAssembly::getNumThreads() const {
    return num_threads_;
}

// ---------- 私有辅助方法 ----------

std::unique_ptr<EMElementIntegratorBase> EMAssembly::createIntegrator(
    fe_em::ElemType elem_type,
    fe_em::DOFType dof_type
) {
    // ========================================================================
    // ⚠️  重要提示：当前使用 DummyEMIntegrator 占位符实现
    // ========================================================================
    //
    // 【状态】待实现（PENDING IMPLEMENTATION）
    // 【优先级】高（HIGH PRIORITY）- 影响装配结果的物理正确性
    // 【计划版本】v1.1.0 - 单元积分器集成
    //
    // 【需要实现的类】:
    //   1. LagrangeEMIntegrator : EMElementIntegratorBase
    //      - 用于标量场单元（SCALAR_ONLY / MIXED_AV的标量部分）
    //      - 支持单元类型: TRI3, TRI6, QUAD4, QUAD8, QUAD9, TET4, TET10,
    //                      HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13
    //      - 物理方程: ∫ ε∇N_i·∇N_j dΩ (静电/静磁)
    //
    //   2. NedelecEMIntegrator : EMElementIntegratorBase
    //      - 用于矢量场单元（VECTOR_EDGE_ONLY / MIXED_AV的矢量部分）
    //      - 支持单元类型: TET4_EDGE, HEX8_EDGE, PRISM6_EDGE, PYRAMID5_EDGE
    //      - 物理方程: ∫ (1/μ)(∇×N_i)·(∇×N_j) dΩ (磁场矢量位A)
    //
    //   3. MixedAVIntegrator : EMElementIntegratorBase （可选，或组合使用上述两个）
    //      - 用于A-V混合格式单元（MIXED_AV）
    //      - 同时处理标量电位V和矢量磁位A
    //
    // 【当前占位符说明】:
    //   DummyEMIntegrator 返回单位矩阵（identity matrix），
    //   仅用于验证组装逻辑的正确性（矩阵维度、对称性、约束处理等）。
    //   **不能用于实际物理仿真！**
    //
    // 【替换步骤】:
    //   1. 实现上述积分器类（参考 em_element_integrator_base.hpp 接口定义）
    //   2. 在此处根据 elem_type 和 dof_type 实例化正确的积分器
    //   3. 删除 DummyEMIntegrator 相关代码
    //   4. 更新本注释为"已完成"
    //
    // ========================================================================

    // 运行时提醒（每次调用都会输出，提醒开发者注意）
    static bool warning_shown = false;
    if (!warning_shown) {
        FEEM_WARN("createIntegrator: ⚠️ 当前使用DummyEMIntegrator占位符！");
        FEEM_WARN("  请尽快实现LagrangeEMIntegrator和NedelecEMIntegrator以获得物理正确结果");
        FEEM_WARN("  详细信息请查看源代码中的TODO注释（em_assembly.cpp:createIntegrator）");
        warning_shown = true;
    }

    int dof_count = 0;

    // 根据单元类型和DOF类型估算自由度数
    // 注意：实际DOF数应由形函数模块提供，此处为简化估计
    switch (elem_type) {
        case numeric::ElementType::TRI3:
        case numeric::ElementType::TRI6:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 3 : 6;
            break;

        case numeric::ElementType::QUAD4:
        case numeric::ElementType::QUAD8:
        case numeric::ElementType::QUAD9:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 4 : 8;
            break;

        case numeric::ElementType::TET4:
        case numeric::ElementType::TET10:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 4 :
                        (dof_type == fe_em::DOFType::VECTOR_EDGE_ONLY) ? 6 : 10;
            break;

        case numeric::ElementType::HEX8:
        case numeric::ElementType::HEX20:
        case numeric::ElementType::HEX27:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 8 :
                        (dof_type == fe_em::DOFType::VECTOR_EDGE_ONLY) ? 12 : 20;
            break;

        case numeric::ElementType::PRISM6:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 6 :
                        (dof_type == fe_em::DOFType::MIXED_AV) ? 15 : 9;
            break;

        case numeric::ElementType::PRISM15:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 15 :
                        (dof_type == fe_em::DOFType::MIXED_AV) ? 30 : 18;
            break;

        case numeric::ElementType::PYRAMID5:
        case numeric::ElementType::PYRAMID13:
            dof_count = (dof_type == fe_em::DOFType::SCALAR_ONLY) ? 5 :
                        (dof_type == fe_em::DOFType::VECTOR_EDGE_ONLY) ? 8 : 13;
            break;

        case numeric::ElementType::TET4_EDGE:
        case numeric::ElementType::HEX8_EDGE:
        case numeric::ElementType::PRISM6_EDGE:
        case numeric::ElementType::PYRAMID5_EDGE:
            dof_count = (dof_type == fe_em::DOFType::VECTOR_EDGE_ONLY) ? 
                        ((elem_type == numeric::ElementType::TET4_EDGE) ? 6 :
                         (elem_type == numeric::ElementType::HEX8_EDGE) ? 12 :
                         (elem_type == numeric::ElementType::PRISM6_EDGE) ? 9 : 8) : 0;
            break;

        default:
            FEEM_ERROR("createIntegrator: 不支持的单元类型（ElemType={})",
                      static_cast<int>(elem_type));
            return nullptr;
    }

    FEEM_DEBUG("createIntegrator: 创建DummyEMIntegrator（DOF={}, ElemType={}, DOFType={})",
               dof_count, static_cast<int>(elem_type), static_cast<int>(dof_type));

    // 返回测试用占位符积分器
    return std::make_unique<DummyEMIntegrator>(dof_count);
}

Eigen::MatrixXd EMAssembly::buildNodeCoords(
    const fe_em::Element& element,
    const std::vector<fe_em::Node>& nodes
) {
    // 提取单元节点数量
    int num_nodes = static_cast<int>(element.node_ids.size());

    if (num_nodes <= 0) {
        throw std::runtime_error("buildNodeCoords: 单元节点列表为空");
    }

    // 构建坐标矩阵（3D坐标，行=维度，列=节点编号）
    Eigen::MatrixXd coords(3, num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
        int node_id = element.node_ids[i];

        // 边界检查：确保节点ID在有效范围内
        if (node_id < 0 || node_id >= static_cast<int>(nodes.size())) {
            throw std::out_of_range(
                "buildNodeCoords: 节点ID=" + std::to_string(node_id) +
                " 越界（有效范围[0, " + std::to_string(nodes.size() - 1) + "]）"
            );
        }

        // 提取节点坐标（假设nodes按ID顺序存储）
        const auto& node = nodes[node_id];
        coords(0, i) = node.x;  // X坐标
        coords(1, i) = node.y;  // Y坐标
        coords(2, i) = node.z;  // Z坐标（2D问题时z=0）
    }

    FEEM_DEBUG("buildNodeCoords: 成功构建{}×{}坐标矩阵", coords.rows(), coords.cols());

    return coords;
}

int EMAssembly::estimateNNZ(
    const std::vector<fe_em::Element>& elements,
    int avg_dof_per_elem
) {
    // 保守估计公式：
    // NNZ_estimate = n_elements × avg_dof_per_elem²
    //
    // 物理依据：
    // 每个单元贡献一个 local_dofs × local_dofs 的稠密子矩阵到全局矩阵，
    // 最坏情况下每个元素都是非零的（不考虑重叠合并）。
    // 实际NNZ会因共享节点的重复坐标合并而小于此估计值。

    if (elements.empty() || avg_dof_per_elem <= 0) {
        return 0;
    }

    int64_t estimate = static_cast<int64_t>(elements.size()) *
                       static_cast<int64_t>(avg_dof_per_elem) *
                       static_cast<int64_t>(avg_dof_per_elem);

    // 防止整数溢出（限制最大值为2^31-1）
    const int MAX_NNZ = 2147483647;
    if (estimate > MAX_NNZ) {
        FEEM_WARN("estimateNNZ: 估算值溢出（{}），截断到{}", estimate, MAX_NNZ);
        estimate = MAX_NNZ;
    }

    return static_cast<int>(estimate);
}

void EMAssembly::scatterElementMatrix(
    CooMatrix<double>& coo_matrix,
    const Eigen::MatrixXd& elem_matrix,
    const fe_em::Local2Global& l2g
) {
    // 获取单元局部DOF数量
    int n_dof = static_cast<int>(l2g.indices.size());

    // 检查维度匹配
    if (elem_matrix.rows() != n_dof || elem_matrix.cols() != n_dof) {
        FEEM_ERROR("scatterElementMatrix: 维度不匹配（矩阵{}×{} vs DOF映射{}）",
                  elem_matrix.rows(), elem_matrix.cols(), n_dof);
        return;
    }

    // 获取全局矩阵尺寸用于边界检查
    int global_size = coo_matrix.rows();
    int error_count = 0;

    // 双重循环遍历单元矩阵的所有元素
    for (int i = 0; i < n_dof; ++i) {
        // 获取第i行对应的全局行索引
        int global_row = l2g.indices[i];

        // 约束跳过规则：行索引为-1表示该DOF受Dirichlet边界条件约束
        if (global_row == -1) {
            continue;  // 跳过整行（消去法边界条件处理）
        }

        // 边界检查（防御性编程）：全局行索引必须在有效范围内
        if (global_row < 0 || global_row >= global_size) {
            FEEM_ERROR("scatterElementMatrix: 全局行索引{}越界（有效范围[0, {}]，局部索引i={})",
                      global_row, global_size - 1, i);
            error_count++;
            continue;
        }

        for (int j = 0; j < n_dof; ++j) {
            // 获取第j列对应的全局列索引
            int global_col = l2g.indices[j];

            // 约束跳过规则：列索引为-1表示该DOF受约束
            if (global_col == -1) {
                continue;  // 跳过该元素
            }

            // 边界检查：全局列索引必须在有效范围内
            if (global_col < 0 || global_col >= global_size) {
                FEEM_ERROR("scatterElementMatrix: 全局列索引{}越界（有效范围[0, {}]，局部索引j={})",
                          global_col, global_size - 1, j);
                error_count++;
                continue;
            }

            // 将局部矩阵元素累加到全局COO矩阵的对应位置
            // 注：COO格式允许重复坐标，将在CSR转换时自动合并
            coo_matrix.add_value(global_row, global_col, elem_matrix(i, j));
        }
    }

    // 如果存在边界错误，输出汇总警告
    if (error_count > 0) {
        FEEM_WARN("scatterElementMatrix: 单元装配完成但存在{}个索引越界错误（已跳过）", error_count);
    }
}

void EMAssembly::scatterElementVector(
    Eigen::VectorXd& global_vector,
    const Eigen::VectorXd& elem_vector,
    const fe_em::Local2Global& l2g
) {
    // 获取单元局部DOF数量
    int n_dof = static_cast<int>(l2g.indices.size());

    // 检查维度匹配
    if (elem_vector.size() != n_dof) {
        FEEM_ERROR("scatterElementVector: 维度不匹配（向量长度{} vs DOF映射{}）",
                  elem_vector.size(), n_dof);
        return;
    }

    // 获取全局向量尺寸用于边界检查
    int global_size = static_cast<int>(global_vector.size());
    int error_count = 0;

    // 单层循环遍历单元向量的所有分量
    for (int i = 0; i < n_dof; ++i) {
        // 获取第i个分量对应的全局行索引
        int global_row = l2g.indices[i];

        // 约束跳过规则：索引为-1表示该DOF受约束
        if (global_row == -1) {
            continue;  // 跳过该分量
        }

        // 边界检查（防御性编程）：全局索引必须在有效范围内
        if (global_row < 0 || global_row >= global_size) {
            FEEM_ERROR("scatterElementVector: 全局索引{}越界（有效范围[0, {}]，局部索引i={})",
                      global_row, global_size - 1, i);
            error_count++;
            continue;
        }

        // 将局部向量分量累加到全局向量的对应位置
        global_vector(global_row) += elem_vector(i);
    }

    // 如果存在边界错误，输出汇总警告
    if (error_count > 0) {
        FEEM_WARN("scatterElementVector: 向量装配完成但存在{}个索引越界错误（已跳过）", error_count);
    }
}

bool EMAssembly::checkSymmetry(const CsrMatrix<double>& matrix) const {
    // 对称性检查算法：
    // 遍历CSR矩阵的所有非零元素A(i,j)，
    // 对于每个非对角线元素（i≠j），查找转置位置(j,i)的元素A(j,i)，
    // 验证 |A(i,j) - A(j,i)| < tolerance。

    // 快速路径：空矩阵或单行/单列矩阵视为对称
    if (matrix.nnz() == 0 || matrix.rows() <= 1) {
        return true;
    }

    // 方阵检查（非方阵必不对称）
    if (matrix.rows() != matrix.cols()) {
        return false;
    }

    // 获取CSR矩阵的内部数据
    const auto& row_ptr = matrix.get_row_ptr();
    const auto& col_indices = matrix.get_col_indices();
    const auto& values = matrix.get_values();

    int n_rows = matrix.rows();

    // 定义数值容差
    // 采用混合容差策略：绝对容差 + 相对容差
    const double TOL_ABS = 1e-10;   // 绝对容差（适用于接近零的元素）
    const double TOL_REL = 1e-8;    // 相对容差（适用于大幅值元素）

    // 遍历每一行
    for (int i = 0; i < n_rows; ++i) {
        // 获取当前行的起始和结束位置
        int row_start = row_ptr[i];
        int row_end = row_ptr[i + 1];

        // 遍历当前行的所有非零元素
        for (int k = row_start; k < row_end; ++k) {
            int j = col_indices[k];      // 列索引
            double a_ij = values[k];      // 元素值 A(i,j)

            // 跳过对角线元素（自动满足对称性）
            if (i == j) {
                continue;
            }

            // 在第j行中查找列索引为i的元素（即A(j,i)）
            bool found_transpose = false;
            double a_ji = 0.0;

            int j_row_start = row_ptr[j];
            int j_row_end = row_ptr[j + 1];

            // 二分查找（利用CSR存储的行列有序性）
            // 注：若col_indices未排序，需退化为线性搜索
            auto it_lower = col_indices.begin() + j_row_start;
            auto it_upper = col_indices.begin() + j_row_end;
            auto it_find = std::lower_bound(it_lower, it_upper, i);

            if (it_find != it_upper && *it_find == i) {
                // 找到转置位置元素
                found_transpose = true;
                int pos = static_cast<int>(std::distance(col_indices.begin(), it_find));
                a_ji = values[pos];
            }

            // 若转置元素不存在或超出容差范围，则矩阵非对称
            if (!found_transpose) {
                FEEM_DEBUG("checkSymmetry: 非对称 - A({},{})={} 存在但 A({},{}) 不存在",
                          i, j, a_ij, j, i);
                return false;
            }

            // 计算偏差（采用混合容差）
            double diff = std::abs(a_ij - a_ji);
            double max_val = std::max(std::abs(a_ij), std::abs(a_ji));
            double tolerance = TOL_ABS + TOL_REL * max_val;

            if (diff > tolerance) {
                FEEM_DEBUG("checkSymmetry: 非对称 - |A({},{})-A({},{})|={} > 容差{}",
                          i, j, j, i, diff, tolerance);
                return false;
            }
        }
    }

    // 所有非零元素均通过对称性检验
    return true;
}

} // namespace numeric
