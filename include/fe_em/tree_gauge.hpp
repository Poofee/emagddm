/**
 * @file tree_gauge.hpp
 * @brief 电磁物理层 - Tree Gauge核心应用类头文件
 * @details 整合MeshGraphExtractor、SpanningTreeBuilder、MultiConnectedDomainDetector三个模块，
 *          实现完整的Tree Gauge规范构建流程。作为整个Tree Gauge模块的核心入口点。
 *
 * @note 核心功能：
 *       - 5步流程：图结构提取 → 生成树构建 → 多连通域检测 → 边分类 → 约化映射
 *       - 树边/余树边分类与Dirichlet边界条件自动兼容
 *       - 多连通域环量自由度的保留与管理
 *       - 约化解到完整解的映射功能
 *
 * @note 数据流设计：
 *       EMMeshData + edge_mapping
 *           ↓ MeshGraphExtractor
 *       GraphResult (邻接表 + 棱边映射)
 *           ↓ SpanningTreeBuilder (带边界优先)
 *       SpanningTreeResult (树边 + 余树边)
 *           ↓ MultiConnectedDomainDetector
 *       DetectionResult (孔洞数 + 环量自由度)
 *           ↓ TreeGauge::classifyEdges()
 *       最终分类结果 (final_tree_edges + final_cotree_edges)
 *           ↓ TreeGauge::buildReducedMapping()
 *       约化DOF映射 (cotree_edge_to_reduced_dof_)
 *
 * @note 与前后模块对接：
 *       - 输入：EMMeshData（网格拓扑）、elem_local_to_global_edge（棱边ID映射）、constrained_edges（约束边）
 *       - 输出：树边/余树边分类、约化DOF映射、解向量转换能力
 *       - 下游：EMDOFManager使用约化DOF信息压缩自由度空间
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include <stdexcept>

#include "em_mesh_data.hpp"
#include "mesh_graph_extractor.hpp"
#include "spanning_tree.hpp"
#include "multi_connected_domain.hpp"
#include "logger_factory.hpp"

namespace fe_em {

/**
 * @class TreeGauge
 * @brief Tree Gauge规范核心应用类
 * @details 整合图提取、生成树构建、多连通域检测三个底层模块，
 *          提供一键式build()接口完成完整的Tree Gauge规范构建。
 *
 * @note 算法原理：
 *       Tree Gauge是电磁场Nedelec棱边元求解中的规范条件选择策略。
 *       其核心思想是将自由度空间划分为：
 *       - **树边（Tree Edges）**：A值被置为0（由生成树定义）
 *       - **余树边（Cotree Edges）**：需要求解的自由度
 *       - **环量自由度（Circulation DOFs）**：多连通域中不可置零的特殊自由度
 *
 *       对于有V个节点、E条棱边、g个孔洞的网格：
 *       - 树边数 = V - C + g （C为连通分量数）
 *       - 余树边数 = E - V + C - g
 *       - 约化DOF数 = 余树边数 + 环量自由度数
 *
 * @note 设计原则：
 *       - 构造时仅保存const引用（零拷贝，无内存开销）
 *       - 延迟计算模式：构造时不执行操作，需显式调用build()
 *       - 结果缓存：中间结果存储在成员变量中供查询
 *       - 异常安全：每个步骤独立捕获异常，提供详细错误信息
 *       - 线程安全：实例级别非线程安全（但可并行创建多个实例）
 *
 * @note 典型使用流程：
 *       @code
 *       // 第一步：准备输入数据
 *       EMMeshData mesh_data = ...;  // 从MaxwellParser加载
 *       GlobalEdgeIDGenerator edge_gen(mesh_data);
 *       edge_gen.generate();
 *       const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();
 *
 *       // 第二步：创建TreeGauge并构建
 *       std::set<int> constrained_edges = {0, 5, 10};  // Dirichlet边界边
 *       TreeGauge gauge(mesh_data, edge_mapping, constrained_edges);
 *       gauge.build();
 *
 *       // 第三步：查询结果
 *       int reduced_dofs = gauge.getReducedNumDofs();
 *       const auto& tree_edges = gauge.getTreeEdges();
 *       const auto& cotree_edges = gauge.getCotreeEdges();
 *
 *       // 第四步：解向量转换（求解后）
 *       std::vector<double> reduced_solution(reduced_dofs);  // 求解器输出
 *       std::vector<double> full_solution;
 *       gauge.mapReducedSolutionToFull(reduced_solution, full_solution);
 *       // full_solution[i] = 第i条全局棱边的A值（树边为0）
 *       @endcode
 *
 * @note 性能特征：
 *       - 时间复杂度：O(V + E) （V=节点数，E=棱边数）
 *       - 空间复杂度：O(V + E) （存储邻接表和各种映射表）
 *       - 构建时间：<50ms（对于中型模型，10000节点级别）
 *       - 内存优化：预分配容器容量，避免动态扩容
 */
class TreeGauge {
public:
    /**
     * @brief 构造函数
     * @param mesh_data 有限元网格数据的常量引用（仅读取，不修改）
     * @param elem_local_to_global_edge 全局棱边ID映射表（来自GlobalEdgeIDGenerator）
     * @param constrained_edges 被约束的棱边集合（默认为空集，可选）
     *
     * @details 构造函数仅保存引用，不执行任何构建操作。
     *          实际的Tree Gauge构建需显式调用 build() 方法。
     *
     * @param constrained_edges 说明：
     *        - 这些边的A值将被强制置为0（归入树边集合）
     *        - 典型用途：Dirichlet边界条件 A×n = 0 的棱边
     *        - 即使这些边在生成树中被归类为余树边，也会被重新分类为树边
     *        - 若不需要特殊约束，可传空集（使用默认参数即可）
     *
     * @warning 调用 build() 前，mesh_data 和 elem_local_to_global_edge 必须已正确填充
     * @throws 无（构造函数不执行实质性操作）
     *
     * @note 内存管理：
     *       - 所有输入通过const引用传入，零拷贝
     *       - 调用者需确保引用的对象在TreeGauge生命周期内有效
     */
    explicit TreeGauge(
        const EMMeshData& mesh_data,
        const std::vector<std::vector<int>>& elem_local_to_global_edge,
        const std::set<int>& constrained_edges = {}
    );

    /**
     * @brief 执行完整的Tree Gauge构建流程（5步算法）
     * @throws std::runtime_error 当任何步骤失败时抛出（包含详细错误信息）
     *
     * @details 核心算法步骤：
     *          **步骤1：提取图结构**
     *          - 使用MeshGraphExtractor从网格数据中提取无向图
     *          - 生成邻接表（adjacency_list_）和棱边映射（edge_to_global_id_）
     *          - 统计节点数（num_nodes_）和棱边数（num_edges_）
     *
     *          **步骤2：构建生成树（带边界优先）**
     *          - 使用SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority()
     *          - 自动将constrained_edges_传入以实现边界优先策略
     *          - 结果存储在tree_result_中（包含树边、余树边、父节点映射）
     *
     *          **步骤3：检测多连通域**
     *          - 使用MultiConnectedDomainDetector::detect()识别孔洞
     *          - 检测环量自由度（不可置零的磁通回路）
     *          - 结果存储在domain_result_中（孔洞数、回路棱边、环量DOFs）
     *
     *          **步骤4：边分类（最终分类）**
     *          - 合并生成树结果和约束边信息
     *          - 生成最终树边集合（final_tree_edges_）
     *          - 生成最终余树边集合（final_cotree_edges_）
     *
     *          **步骤5：构建约化映射**
     *          - 为每个余树边分配约化DOF索引
     *          - 为每个环量自由度分配约化DOF索引
     *          - 建立全局棱边ID → 约化DOF索引的双向映射
     *
     * @warning 此方法只能调用一次（重复调用会抛出异常）
     * @note 每个步骤都有详细的日志输出，便于调试和性能分析
     * @note 日志输出示例：
     *       [INFO] ========== 开始Tree Gauge构建 ==========
     *       [INFO] 步骤1[图结构提取]完成 - 节点数: 500, 棱边数: 1200
     *       [INFO] 步骤2[生成树构建]完成 - 树边数: 499, 余树边数: 701
     *       [INFO] 步骤3[多连通域检测]完成 - 孔洞数: 1, 环量自由度数: 1
     *       [INFO] 步骤4[边分类]完成 - 最终树边: 500, 最终余树边: 699
     *       [INFO] 步骤5[约化映射]完成 - 约化DOF数: 700
     *       [INFO] ========== Tree Gauge构建完成 ==========
     *       [INFO] 压缩率: 41.7% (原1200→700)
     */
    void build();

    // ==================== 查询接口 ====================

    /**
     * @brief 获取约化后的总自由度数
     * @return int 约化DOF数量（余树边数 + 环量自由度数）
     *
     * @note 必须在 build() 后调用，否则返回0
     * @note 返回值含义：
     *       - 这是求解器实际需要求解的自由度数量
     *       - 相比原始的总棱边数，通常能减少30-60%
     *       - 减少程度取决于网格拓扑（越接近树状，压缩率越高）
     */
    int getReducedNumDofs() const;

    /**
     * @brief 获取最终树边集合（A=0的边）
     * @return const std::set<int>& 树边的全局棱边ID集合（常量引用）
     *
     * @note 树边包括：
     *       - 生成树中的所有边
     *       - 被约束的余树边（强制A=0）
     * @note 必须在 build() 后调用
     */
    const std::set<int>& getTreeEdges() const;

    /**
     * @brief 获取最终余树边集合（待求解的边）
     * @return const std::set<int>& 余树边的全局棱边ID集合（常量引用）
     *
     * @note 余树边不包括：
     *       - 树边（已排除）
     *       - 被约束的边（已归入树边）
     *       - 环量自由度（单独处理，但也属于约化DOF的一部分）
     * @note 必须在 build() 后调用
     */
    const std::set<int>& getCotreeEdges() const;

    /**
     * @brief 获取环量自由度集合（多连通域中的特殊自由度）
     * @return const std::set<int>& 环量自由度的全局棱边ID集合（常量引用）
     *
     * @note 环量自由度特点：
     *       - 属于余树边的子集
     *       - 代表穿过孔洞的独立磁通回路
     *       - 不可被置为零（否则违反物理守恒定律）
     *       - 在单连通域中此集合为空
     * @note 必须在 build() 后调用
     */
    const std::set<int>& getCirculationDOFs() const;

    /**
     * @brief 获取余树边到约化DOF索引的映射
     * @return const std::unordered_map<int, int>& 映射表（键：全局棱边ID，值：约化DOF索引）
     *
     * @note 映射关系说明：
     *       - 键（key）：全局棱边ID（来自GlobalEdgeIDGenerator）
     *       - 值（value）：在reduced_solution向量中的索引位置
     *       - 索引范围：[0, getReducedNumDofs()-1]
     *       - 映射包含：普通余树边 + 环量自由度
     * @note 典型用法：
     *       @code
     *       const auto& mapping = gauge.getCotreeEdgeToReducedDOF();
     *       for (auto [edge_id, reduced_idx] : mapping) {
     *           reduced_solution[reduced_idx] = ...;  // 组装RHS或提取解
     *       }
     *       @endcode
     * @note 必须在 build() 后调用
     */
    const std::unordered_map<int, int>& getCotreeEdgeToReducedDOF() const;

    // ==================== 解向量操作 ====================

    /**
     * @brief 将约化解向量映射回完整解向量
     * @param reduced_solution 约化后的解向量（长度=getReducedNumDofs()）
     * @param[out] full_solution 输出的完整解向量（长度=num_edges_，树边位置为0）
     *
     * @details 映射规则：
     *          1. 初始化full_solution为全0向量（长度=总棱边数）
     *          2. 遍历cotree_edge_to_reduced_dof_映射表
     *          3. 对每个余树边/环量自由度：
     *             - 从reduced_solution中读取对应位置的值
     *             - 写入full_solution的全局棱边ID对应的位置
     *          4. 树边位置保持0.0（Tree Gauge规范要求）
     *
     * @note 输入输出说明：
     *       - reduced_solution：求解器的输出（仅包含约化DOF的值）
     *       - full_solution：完整解（可用于后处理可视化）
     *       - full_solution会被resize和填充，无需预先分配
     *
     * @warning reduced_solution的长度必须等于getReducedNumDofs()
     * @throws std::invalid_argument 当输入维度不匹配时抛出
     * @note 必须在 build() 后调用
     *
     * @note 典型用法：
     *       @code
     *       // 求解器得到约化解
     *       std::vector<double> reduced_sol(gauge.getReducedNumDofs());
     *       solver.solve(reduced_sol);  // 假设的求解器接口
     *
     *       // 转换为完整解用于后处理
     *       std::vector<double> full_sol;
     *       gauge.mapReducedSolutionToFull(reduced_sol, full_sol);
     *
     *       // full_sol[i] = 第i条全局棱边的A值
     *       for (int i = 0; i < full_sol.size(); ++i) {
     *           if (gauge.isTreeEdge(i)) {
     *               assert(full_sol[i] == 0.0);  // 树边必须为0
     *           }
     *       }
     *       @endcode
     */
    void mapReducedSolutionToFull(
        const std::vector<double>& reduced_solution,
        std::vector<double>& full_solution
    ) const;

    // ==================== 判断方法 ====================

    /**
     * @brief 判断指定棱边是否为树边
     * @param global_edge_id 全局棱边ID
     * @return bool 如果是树边返回true，否则返回false
     *
     * @note 时间复杂度：O(log N)，N为树边数量（使用std::set查找）
     * @note 必须在 build() 后调用，否则返回false
     */
    bool isTreeEdge(int global_edge_id) const;

    /**
     * @brief 判断指定棱边是否为余树边
     * @param global_edge_id 全局棱边ID
     * @return bool 如果是余树边返回true，否则返回false
     *
     * @note 注意：环量自由度虽然在技术上是余树边，
     *       但此方法返回的是"纯余树边"（不包括环量自由度）
     *       若需判断是否为环量自由度，请使用getCirculationDOFs().count()
     * @note 时间复杂度：O(log N)，N为余树边数量
     * @note 必须在 build() 后调用，否则返回false
     */
    bool isCotreeEdge(int global_edge_id) const;

private:
    // ==================== 输入数据（const引用，零拷贝）====================

    const EMMeshData& mesh_data_;                              ///< 网格数据常量引用
    const std::vector<std::vector<int>>& elem_local_to_global_edge_;  ///< 全局棱边ID映射表引用
    const std::set<int> constrained_edges_;                     ///< 被约束的棱边集合

    // ==================== 图结构数据（来自MeshGraphExtractor）====================

    std::vector<std::vector<int>> adjacency_list_;              ///< 节点邻接表
    std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> edge_to_global_id_;  ///< 棱边→全局ID映射
    int num_nodes_ = 0;                                         ///< 图中总节点数
    int num_edges_ = 0;                                         ///< 图中总棱边数

    // ==================== 生成树结果（来自SpanningTreeBuilder）====================

    SpanningTreeResult tree_result_;                            ///< 生成树构建结果

    // ==================== 多连通域检测结果（来自MultiConnectedDomainDetector）====================

    DetectionResult domain_result_;                             ///< 多连通域检测结果

    // ==================== 最终分类结果 ====================

    std::set<int> final_tree_edges_;                            ///< 最终树边集合（A=0）
    std::set<int> final_cotree_edges_;                          ///< 最终余树边集合（待求解）
    std::unordered_map<int, int> cotree_edge_to_reduced_dof_;   ///< 余树边→约化DOF索引映射
    int reduced_num_dofs_ = 0;                                  ///< 约化总自由度数

    // ==================== 内部状态标记 ====================

    bool is_built_ = false;                                     ///< 是否已完成build()

    // ==================== 内部方法（5步流程）====================

    /**
     * @brief 步骤1：提取图结构
     * @throws std::runtime_error 当图提取失败时抛出
     * @details 使用MeshGraphExtractor从网格数据中提取无向图结构，
     *          并将结果保存到adjacency_list_、edge_to_global_id_、num_nodes_、num_edges_
     */
    void extractGraphStructure();

    /**
     * @brief 步骤2：构建生成树（带边界优先策略）
     * @throws std::runtime_error 当生成树构建失败时抛出
     * @details 使用SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority()构建生成树，
     *          将constrained_edges_作为Dirichlet边界传入以确保边界兼容性
     */
    void buildSpanningTree();

    /**
     * @brief 步骤3：检测多连通域
     * @throws std::runtime_error 当检测失败时抛出
     * @details 使用MultiConnectedDomainDetector::detect()识别孔洞和环量自由度，
     *          结果保存在domain_result_中
     */
    void detectMultiConnectedDomain();

    /**
     * @brief 步骤4：边分类（合并生成树结果和约束信息）
     * @throws 无（内部逻辑不会失败）
     * @details 分类规则：
     *          - final_tree_edges_ = tree_result_.tree_edges ∪ (被约束的余树边)
     *          - final_cotree_edges_ = tree_result_.cotree_edges - 被约束的边 - 环量自由度
     */
    void classifyEdges();

    /**
     * @brief 步骤5：构建约化DOF映射
     * @throws 无（内部逻辑不会失败）
     * @details 为每个最终余树边和环量自由度分配连续的约化DOF索引，
     *          建立cotree_edge_to_reduced_dof_映射表
     */
    void buildReducedMapping();
};

} // namespace fe_em
