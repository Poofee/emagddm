/**
 * @file em_dof_manager.hpp
 * @brief 电磁物理层 - 自由度(DOF)管理器核心类头文件
 * @details 实现有限元离散系统的自由度编号、约束处理和映射表构建，
 *          是电磁场求解器的最核心模块。支持标量节点DOF和矢量棱边DOF的混合管理。
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 *
 * @note 核心算法说明：
 *       采用四步编号流程实现DOF管理：
 *       1. 预编号（assignPrenumbering）：为所有物理实体分配临时连续编号
 *       2. 标记约束（markConstrainedDOFs）：识别Dirichlet边界条件对应的受约束DOF
 *       3. 分块重编号（renumberFreeDOFs）：消去法重编号，自由DOF在前从0开始连续编号
 *       4. 构建映射表（buildMappingTables）：生成每个单元的Local2Global映射关系
 *
 * @note 与后续模块对接：
 *       - 单元组装模块：使用getElemLocalToGlobal()获取映射表进行矩阵scatter
 *       - 求解器模块：使用getNumFreeDOFs()确定系统维度
 *       - 边界条件施加：使用getConstrainedDOFValues()获取约束值向量
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <set>
#include <memory>
#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"
#include "element_geometry.hpp"
#include "global_edge_id_generator.hpp"
#include "tree_gauge.hpp"
#include "logger_factory.hpp"

namespace fe_em {

/**
 * @class EMDOFManager
 * @brief 电磁场有限元自由度管理器（核心类）
 * @details 管理整个离散系统的自由度分配、约束识别和局部-全局映射构建。
 *          支持三种DOF类型的混合网格：
 *          - SCALAR_ONLY: 纯标量单元（电位φ、磁位ψ）
 *          - VECTOR_EDGE_ONLY: 纯矢量单元（磁矢势A的切向分量）
 *          - MIXED_AV: A-V混合格式单元（涡流场导体区）
 *
 * @note 设计原则：
 *       - 构造时仅保存引用，不拷贝网格数据（零开销）
 *       - 必须显式调用build()执行完整编号流程
 *       - 内部采用四步流程确保逻辑清晰、可维护性强
 *       - 所有中间数据在四步间通过私有成员传递
 *
 * @note 典型使用流程：
 *       @code
 *       // 1. 创建DOF管理器（传入网格数据和全局棱边映射）
 *       EMDOFManager dof_manager(mesh_data, edge_mapping);
 *
 *       // 2. 执行完整DOF编号流程
 *       dof_manager.build();
 *
 *       // 3. 获取结果
 *       int num_free = dof_manager.getNumFreeDOFs();
 *       const auto& mappings = dof_manager.getElemLocalToGlobal();
 *       const auto& constrained_vals = dof_manager.getConstrainedDOFValues();
 *       @endcode
 */
class EMDOFManager {
public:
    /**
     * @brief 构造函数
     * @param mesh_data 有限元网格拓扑数据的常量引用（仅读取，不修改）
     * @param elem_local_to_global_edge 全局棱边ID映射表（来自GlobalEdgeIDGenerator）
     *        可为空向量（当网格中无矢量单元时），默认为空
     * @param edge_gen 全局棱边ID生成器的常量指针（可选，用于O(1)棱边查找）
     *        传入后可显著提升矢量边界约束标记的效率。
     *        若为nullptr，则回退到O(n)遍历方式
     *
     * @note 构造函数仅保存引用和数据副本，不执行任何编号操作。
     *       实际的DOF编号需显式调用 build() 方法。
     * @warning 调用build()前确保mesh_data已正确填充节点、单元和边界标记数据
     */
    explicit EMDOFManager(const EMMeshData& mesh_data,
                          const std::vector<std::vector<int>>& elem_local_to_global_edge = {},
                          const GlobalEdgeIDGenerator* edge_gen = nullptr);

    /**
     * @brief 执行完整DOF编号流程
     * @details 按顺序调用四个内部步骤：
     *          1. assignPrenumbering(): 为节点和棱边分配预编号
     *          2. markConstrainedDOFs(): 识别Dirichlet边界条件的受约束DOF
     *          3. renumberFreeDOFs(): 消去法重编号，自由DOF在前连续编号
     *          4. buildMappingTables(): 生成每个单元的Local2Global映射表
     *
     * @note 完成后输出统计日志（总DOF数、自由DOF数、约束DOF数、单元映射表大小）
     * @warning 此方法可重复调用，每次调用会重置之前的结果
     * @exception 若遇到不支持的单元类型或数据不一致，抛出std::invalid_argument
     */
    void build();

    /**
     * @brief 获取自由DOF总数
     * @return int 自由DOF的数量（经约束处理后重新连续编号的总数）
     *
     * @note 必须在build()调用后使用，否则返回0
     * @note 返回值即为全局系统矩阵的维度
     */
    int getNumFreeDOFs() const;

    /**
     * @brief 获取每个单元的Local2Global映射表
     * @return const std::vector<Local2Global>& 映射表向量的常量引用
     *         - 外层维度：按单元索引（size == elements.size()）
     *         - 每个元素：该单元的局部DOF到全局DOF的映射
     *
     * @details Local2Global.indices的含义：
     *          - indices[i] == -1: 第i个局部DOF为约束DOF（Dirichlet边界）
     *          - indices[i] >= 0: 第i个局部DOF为自由DOF，值为全局系统中的编号
     *
     * @note 对于MIXED_AV单元：
     *       - 前num_scalar_dofs个：标量节点DOF（按节点顺序）
     *       - 后num_vector_dofs个：矢量棱边DOF（按局部棱边顺序）
     *
     * @note 必须在build()调用后使用
     * @note 返回值为const引用，外部不可修改
     */
    const std::vector<Local2Global>& getElemLocalToGlobal() const;

    /**
     * @brief 获取约束DOF值向量（按全局预编号顺序）
     * @return const std::vector<double>& 约束值向量的常量引用
     *         - 向量长度等于预编号总数（包含自由和约束DOF）
     *         - 约束DOF位置存储对应的Dirichlet边界值
     *         - 自由DOF位置存储0.0（占位符）
     *
     * @note 此向量用于求解后恢复完整解向量（含约束DOF的值）
     * @note 必须在build()调用后使用
     */
    const std::vector<double>& getConstrainedDOFValues() const;

    /**
     * @brief 应用Tree Gauge规范（可选调用）
     * @details 在build()之后调用，进一步压缩矢量棱边DOF：
     *          - 树边DOF从系统中剔除（A=0）
     *          - 仅保留余树边DOF参与求解
     *          - 显著减少系统规模（通常减少30-50%）
     *
     * @note 调用后，getNumFreeDOFs()返回的是约化后的DOF数量
     * @note getElemLocalToGlobal()中的树边DOF标记为-1
     * @warning 此方法仅对VECTOR_EDGE_ONLY和MIXED_AV类型的矢量DOF有效
     *
     * @throws std::runtime_error 当未先调用build()时抛出异常
     *
     * @note 典型使用流程：
     *       @code
     *       dof_manager.build();
     *       int original_dofs = dof_manager.getNumFreeDOFs();
     *       dof_manager.applyTreeGauge();
     *       int reduced_dofs = dof_manager.getNumFreeDOFs();
     *       // reduced_dofs < original_dofs (压缩成功)
     *       @endcode
     */
    void applyTreeGauge();

    /**
     * @brief 获取Tree Gauge对象（用于查询树边/余树边信息）
     * @return const TreeGauge* Tree Gauge对象的常量指针（未调用applyTreeGauge时返回nullptr）
     *
     * @note 返回的指针生命周期与EMDOFManager实例相同
     * @note 可用于查询树边/余树边分类、环量自由度等信息
     * @note 也可用于解向量恢复（mapReducedSolutionToFull）
     */
    const TreeGauge* getTreeGauge() const;

private:
    // ==================== 输入数据 ====================
    const EMMeshData& mesh_data_;                              ///< 网格数据常量引用
    const std::vector<std::vector<int>> elem_local_to_global_edge_; ///< 全局棱边ID映射表
    const GlobalEdgeIDGenerator* edge_gen_ = nullptr;            ///< 全局棱边ID生成器指针（可选，用于O(1)棱边查找）

    // ==================== 输出结果 ====================
    int num_free_dofs_ = 0;                                    ///< 自由DOF总数
    std::vector<Local2Global> elem_local_to_global_;           ///< 单元局部→全局映射表
    std::vector<double> constrained_dof_values_;               ///< 约束DOF值向量（预编号顺序）

    // ==================== Tree Gauge集成数据 ====================
    std::unique_ptr<TreeGauge> tree_gauge_;                    ///< Tree Gauge对象（延迟初始化）
    int original_num_free_dofs_ = 0;                           ///< applyTreeGauge前的原始自由DOF数
    std::unordered_map<int, int> free_to_reduced_;             ///< 原始自由编号→约化编号映射

    // ==================== 中间数据（四步间传递）====================

    /**
     * @struct PrenumberingData
     * @brief 预编号阶段的中间数据结构
     * @details 存储第一步assignPrenumbering()的结果，供后续步骤使用
     */
    struct PrenumberingData {
        std::unordered_map<int, int> node_to_prenum;   ///< 节点全局ID → 预编号（标量DOF）
        std::unordered_map<int, int> edge_to_prenum;   ///< 全局棱边ID → 预编号（矢量DOF）
        int total_prenum = 0;                          ///< 预编号总数（标量+矢量）
        int num_scalar_prenum = 0;                     ///< 标量预编号数量
        int num_vector_prenum = 0;                     ///< 矢量预编号数量
    };

    /**
     * @struct ConstraintData
     * @brief 约束标记阶段的中间数据结构
     * @details 存储第二步markConstrainedDOFs()的结果
     */
    struct ConstraintData {
        std::set<int> constrained_scalar_nodes;        ///< 受约束的标量节点ID集合
        std::set<int> constrained_edges;               ///< 受约束的全局棱边ID集合
    };

    PrenumberingData prenum_data_;                      ///< 预编号数据
    ConstraintData constraint_data_;                    ///< 约束标记数据
    std::unordered_map<int, int> prenum_to_free_;       ///< 预编号 → 自由编号映射（第三步结果）

    // ==================== 内部方法：四步编号流程 ====================

    /**
     * @brief 第一步：预编号
     * @details 为所有物理实体分配临时连续编号：
     *          - 标量节点DOF从0开始连续编号（遍历被SCALAR_ONLY/MIXED_AV单元引用的节点）
     *          - 矢量棱边DOF紧接标量后续编号（遍历VECTOR_EDGE_ONLY/MIXED_AV单元的棱边）
     *
     * @note 结果存储在prenum_data_中，供后续步骤使用
     * @note 时间复杂度：O(N_nodes + N_edges)，其中N_edges为矢量单元总棱边数
     */
    void assignPrenumbering();

    /**
     * @brief 第二步：标记约束DOF
     * @details 遍历所有边界标记（EMBoundaryMarker），识别Dirichlet边界条件对应的受约束DOF：
     *          - 标量DOF约束：当 bnd_type==DIRICHLET 且 target_ids 持有 vector<int>（节点列表）时，
     *            将每个节点ID加入 constraint_data_.constrained_scalar_nodes 集合
     *          - 矢量DOF约束：当 bnd_type==DIRICHLET 且 target_ids 持有 vector<vector<int>>（棱边节点对列表）时，
     *            通过 edge_gen_->getGlobalEdgeID() O(1)查找全局棱边ID，加入 constraint_data_.constrained_edges 集合
     *            同时将 Dirichlet 值记录到 constrained_dof_values_ 的对应预编号位置
     *
     * @note 查找策略：
     *       - 若 edge_gen_ 非空，使用 O(1) 哈希查找（推荐）
     *       - 若 edge_gen_ 为空，回退到 O(N_elements × N_edges) 遍历查找
     *
     * @note 结果存储在 constraint_data_ 和 constrained_dof_values_ 中
     * @note constrained_dof_values_ 按预编号顺序排列：约束位置存实际Dirichlet值，自由位置存0.0
     *
     * @note 时间复杂度：
     *       - 有 edge_gen_: O(N_boundary_markers × N_target_entities) 平均O(1)每次查找
     *       - 无 edge_gen_: O(N_boundary_markers × N_target_entities × N_elements × N_edges_per_elem)
     *
     * @warning 异常场景处理：
     *       - target_ids 不匹配任何已知 variant 类型：跳过该边界标记并输出 WARN 日志
     *       - 棱边节点对中的节点数<2：跳过该条目并输出 WARN 日志
     *       - 棱边全局ID未找到（可能属于标量单元）：输出 WARN 日志，不标记为约束
     */
    void markConstrainedDOFs();

    /**
     * @brief 第三步：分块重编号（消去法）
     * @details 将自由DOF重新连续编号，消除约束DOF造成的编号间隙：
     *          - 自由DOF在前，从0开始连续编号
     *          - 约束DOF标记为-1（在最终映射表中体现）
     *          - 构建prenum_to_free_映射：原始预编号 → 新的自由编号
     *
     * @note 结果存储在prenum_to_free_和num_free_dofs_中
     * @note 这是"消去法"的核心步骤，将约束DOF从系统中剔除
     */
    void renumberFreeDOFs();

    /**
     * @brief 第四步：构建单元映射表
     * @details 遍历所有单元，根据dof_type分别生成Local2Global映射：
     *          - SCALAR_ONLY单元：indices.size() = 节点数，按节点顺序查映射
     *          - VECTOR_EDGE_ONLY单元：indices.size() = 棱边数，按棱边顺序查映射
     *          - MIXED_AV单元：indices.size() = 节点数 + 棱边数，前半段标量+后半段矢量
     *
     * @note 对于被约束的DOF，indices中填入-1
     * @note 结果存储在elem_local_to_global_中
     */
    void buildMappingTables();

    // ==================== 内部辅助方法 ====================

    /**
     * @brief 判断DOFType是否涉及标量分量
     * @param dof_type 单元DOF类型枚举
     * @return true 如果包含标量DOF（SCALAR_ONLY或MIXED_AV）
     * @return false 如果纯矢量（VECTOR_EDGE_ONLY）
     */
    static bool hasScalarDOF(DOFType dof_type);

    /**
     * @brief 判断DOFType是否涉及矢量分量
     * @param dof_type 单元DOF类型枚举
     * @return true 如果包含矢量DOF（VECTOR_EDGE_ONLY或MIXED_AV）
     * @return false 如果纯标量（SCALAR_ONLY）
     */
    static bool hasVectorDOF(DOFType dof_type);

    /**
     * @brief 构建单个SCALAR_ONLY单元的Local2Global映射
     * @param elem 单元数据的常量引用
     * @return Local2Global 该单元的局部-全局映射
     */
    Local2Global buildScalarElementMapping(const Element& elem) const;

    /**
     * @brief 构建单个VECTOR_EDGE_ONLY单元的Local2Global映射
     * @param elem_id 单元全局ID（用于索引elem_local_to_global_edge_）
     * @param elem_type 单元类型枚举
     * @return Local2Global 该单元的局部-全局映射
     */
    Local2Global buildVectorElementMapping(int elem_id, ElemType elem_type) const;

    /**
     * @brief 构建单个MIXED_AV单元的Local2Global映射
     * @param elem 单元数据的常量引用
     * @return Local2Global 该单元的局部-全局映射（前半段标量+后半段矢量）
     */
    Local2Global buildMixedAVElementMapping(const Element& elem) const;

    /**
     * @brief 使用Tree Gauge结果重新编号矢量DOF
     * @details 在applyTreeGauge()中调用，根据Tree Gauge的分类结果重新编号：
     *          - 树边DOF标记为-1（从系统中剔除）
     *          - 余树边DOF分配新的约化编号
     *          - 标量DOF保持不变
     *
     * @note 此方法会修改elem_local_to_global_和num_free_dofs_
     * @note 标量DOF编号保持不变，仅对矢量棱边DOF进行重编号
     */
    void renumberVectorDOFsWithTreeGauge();
};

} // namespace fe_em
