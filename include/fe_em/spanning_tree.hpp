/**
 * @file spanning_tree.hpp
 * @brief 电磁物理层 - 生成树构建器头文件
 * @details 提供BFS/DFS算法从无向图中构建生成树（Spanning Tree），
 *          用于Tree Gauge规范的树边/余树边分类。
 *
 * @note 核心功能：
 *       - BFS生成树构建（广度优先，缓存友好，推荐使用）
 *       - DFS生成树构建（深度优先，显式栈实现避免递归溢出）
 *       - Dirichlet边界优先策略（确保边界条件与Tree Gauge一致）
 *       - 智能根节点选择（支持多种策略枚举）
 *
 * @note 与前后模块对接：
 *       - MeshGraphExtractor：提供邻接表和棱边映射作为输入
 *       - TreeGauge：使用树边/余树边分类结果执行规范条件
 *       - EMDOFManager：基于分类结果压缩自由度空间
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include "global_edge_id_generator.hpp"
#include "logger_factory.hpp"

namespace fe_em {

// ==================== 枚举类型定义 ====================

/**
 * @enum RootNodeStrategy
 * @brief 根节点选择策略枚举
 * @details 控制BFS/DFS算法起始节点的选择方式，
 *          不同策略影响生成树的形状和边界兼容性。
 */
enum class RootNodeStrategy {
    FIRST_NODE,           ///< 选择索引0的节点（默认，简单快速）
    MAX_DEGREE,           ///< 选择度数最大的节点（推荐，生成更平衡的树）
    DIRICHLET_BOUNDARY,   ///< 从Dirichlet边界棱边中选择一个节点（保证边界兼容性）
    USER_SPECIFIED        ///< 使用用户指定的节点ID
};

// ==================== 输出数据结构 ====================

/**
 * @struct SpanningTreeResult
 * @brief 生成树构建结果数据结构
 * @details 包含完整的生成树信息，用于后续Tree Gauge处理：
 *          - 树边集合：这些边的A值将被置为0
 *          - 余树边集合：需要求解的自由度
 *          - 父节点映射：用于路径重建和调试
 *          - 根节点ID：树的根
 *
 * @note 数据结构设计原则：
 *       - tree_edges和cotree_edges使用std::set（自动去重、有序遍历）
 *       - node_parent使用unordered_map（O(1)查找性能）
 *       - 所有棱边ID均为全局棱边ID（来自edge_to_global_id映射）
 */
struct SpanningTreeResult {
    std::set<int> tree_edges;                                          ///< 树边的全局棱边ID集合（A=0的边）
    std::set<int> cotree_edges;                                        ///< 余树边的全局棱边ID集合（待求解的边）
    std::unordered_map<int, int> node_parent;                          ///< 节点索引 → 父节点索引的映射（根节点的父节点为-1）
    int root_node = 0;                                                 ///< 生成树的根节点索引

    /**
     * @brief 检查结果是否有效
     * @return bool 如果包含有效数据返回true
     */
    bool isValid() const {
        return root_node >= 0 && !node_parent.empty();
    }
};

// ==================== 核心工具类 ====================

/**
 * @class SpanningTreeBuilder
 * @brief 生成树构建器
 * @details 基于图论算法从无向图中构建生成树，支持BFS和DFS两种策略，
 *          并提供Dirichlet边界优先选项以确保边界条件一致性。
 *
 * @note 算法流程（以BFS为例）：
 *       1. **根节点选择**：根据RootNodeStrategy选择起始节点
 *       2. **初始化**：标记根节点为已访问，加入队列/栈
 *       3. **遍历**：按BFS/DFS顺序访问所有可达节点
 *       4. **分类**：首次访问的边归入tree_edges，其余归入cotree_edges
 *       5. **多连通分量**：对每个未访问分量重复上述过程
 *
 * @note 设计原则：
 *       - 静态方法设计（无状态，线程安全）
 *       - 输入来自MeshGraphExtractor::GraphResult
 *       - 时间复杂度O(V+E)，空间复杂度O(V+E)
 *       - 支持百万级节点的大规模网格
 *       - 处理非连通图（多个connected components）
 *
 * @note 典型使用流程：
 *       @code
 *       // 第一步：提取图结构（已在前面完成）
 *       GraphResult graph = extractor.extract();
 *
 *       // 第二步：构建BFS生成树
 *       auto result = SpanningTreeBuilder::buildBFSTree(
 *           graph.adjacency_list,
 *           graph.edge_to_global_id,
 *           0  // 从节点0开始
 *       );
 *
 *       // 第三步：使用结果
 *       int tree_edge_count = result.tree_edges.size();      // V - C (C为连通分量数)
 *       int cotree_edge_count = result.cotree_edges.size();  // E - V + C
 *       @endcode
 *
 * @note 性能特征：
 *       - BFS：队列操作，缓存友好，适合大规模稀疏图
 *       - DFS：栈操作，可能产生更深但更窄的树
 *       - 边界优先：略微增加开销，但保证Dirichlet边在树中
 *
 * @note 日志输出示例：
 *       [INFO] ========== SpanningTreeBuilder: BFS生成树构建完成 ==========
 *       [INFO] 输入统计: 节点数=10000, 棱边数=55000, 连通分量数=1
 *       [INFO] 根节点选择: 节点42(策略=MAX_DEGREE, 度数=12)
 *       [INFO] 生成树: 树边数=9999, 余树边数=45001
 *       [INFO] =======================================================
 */
class SpanningTreeBuilder {
public:
    /**
     * @brief 使用BFS算法构建生成树（推荐方法）
     * @param adjacency_list 节点邻接表（adjacency_list[i] = 节点i的邻居列表）
     * @param edge_to_global_id 棱边键→全局棱边ID的映射（键为排序后的节点对）
     * @param root_node 起始根节点索引（默认为0）
     * @return SpanningTreeResult 包含树边、余树边、父节点映射的完整结果
     *
     * @details BFS算法步骤：
     *          1. 初始化visited数组，标记root_node为已访问
     *          2. 将root_node加入std::queue
     *          3. 循环直到队列为空：
     *             a) 出队当前节点u
     *             b) 遍历u的所有邻居v
     *             c) 若v未被访问：
     *                - 标记v为已访问
     *                - 通过edge_to_global_id查找边(u,v)的全局ID
     *                - 将该全局ID加入tree_edges
     *                - 设置node_parent[v] = u
     *                - 将v入队
     *          4. 所有不在tree_edges中的边即为cotree_edges
     *
     * @warning adjacency_list和edge_to_global_id必须来自同一个GraphResult实例
     * @note 此方法可正确处理非连通图（会为每个连通分量构建子树）
     *
     * @throws 无（内部错误通过日志输出警告）
     */
    static SpanningTreeResult buildBFSTree(
        const std::vector<std::vector<int>>& adjacency_list,
        const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
        int root_node = 0
    );

    /**
     * @brief 使用DFS算法构建生成树（备选方法）
     * @param adjacency_list 节点邻接表
     * @param edge_to_global_id 棱边键→全局棱边ID的映射
     * @param root_node 起始根节点索引（默认为0）
     * @return SpanningTreeResult 完整的生成树结果
     *
     * @details DFS算法特点：
     *          - 使用显式栈（std::stack）而非递归（避免栈溢出）
     *          - 产生深度优先的遍历顺序
     *          - 可能生成比BFS更高但更窄的树
     *          - 对于某些拓扑结构可能有更好的局部性
     *
     * @note 算法逻辑与buildBFSTree相同，仅遍历顺序不同
     * @note 推荐优先使用buildBFSTree（更好的缓存行为）
     */
    static SpanningTreeResult buildDFSTree(
        const std::vector<std::vector<int>>& adjacency_list,
        const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
        int root_node = 0
    );

    /**
     * @brief 使用带Dirichlet边界优先策略的BFS构建生成树
     * @param adjacency_list 节点邻接表
     * @param edge_to_global_id 棱边键→全局棱边ID的映射
     * @param dirichlet_edges Dirichlet边界棱边的全局ID集合（这些边应优先包含在树中）
     * @param strategy 根节点选择策略（默认为DIRICHLET_BOUNDARY）
     * @param user_root_node 用户指定的根节点（仅当strategy=USER_SPECIFIED时有效）
     * @return SpanningTreeResult 完整的生成树结果
     *
     * @details 边界优先策略原理：
     *          1. **根节点选择**：若strategy=DIRICHLET_BOUNDARY且dirichlet_edges非空，
     *             则从dirichlet_edges中选择一条边的一个端点作为根节点
     *          2. **优先遍历**：在BFS遍历每个节点的邻居时，
     *             先检查邻居连接是否为Dirichlet边，若是则优先处理
     *          3. **效果**：确保所有Dirichlet边界棱边被包含在tree_edges中，
     *             使其A值被置为0，与Dirichlet边界条件A=0一致
     *
     * @note 典型应用场景：
     *       - 静磁场求解器：强制磁矢位在边界上为0
     *       - 保证数值解满足物理边界条件
     *       - 避免因规范选择导致的边界不一致问题
     *
     * @warning dirichlet_edges中的ID必须是有效的全局棱边ID
     * @note 若dirichlet_edges为空，则退化为普通BFS
     */
    static SpanningTreeResult buildBFSTreeWithBoundaryPriority(
        const std::vector<std::vector<int>>& adjacency_list,
        const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
        const std::set<int>& dirichlet_edges,
        RootNodeStrategy strategy = RootNodeStrategy::DIRICHLET_BOUNDARY,
        int user_root_node = 0
    );

private:
    /**
     * @brief 根据策略选择根节点
     * @param adjacency_list 节点邻接表（用于计算度数）
     * @param edge_to_global_id 棱边映射（用于DIRICHLET_BOUNDARY策略）
     * @param dirichlet_edges Dirichlet边界棱边集合
     * @param strategy 根节点选择策略
     * @param user_root_node 用户指定的根节点
     * @param num_nodes 总节点数
     * @return int 选定的根节点索引
     *
     * @details 各策略的实现逻辑：
     *          - FIRST_NODE：直接返回0
     *          - MAX_DEGREE：遍历adjacency_list找最大size的节点
     *          - DIRICHLET_BOUNDARY：从dirichlet_edges反向查找到端点节点
     *          - USER_SPECIFIED：返回user_root_node（需验证范围）
     */
    static int selectRootNode(
        const std::vector<std::vector<int>>& adjacency_list,
        const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
        const std::set<int>& dirichlet_edges,
        RootNodeStrategy strategy,
        int user_root_node,
        int num_nodes
    );

    /**
     * @brief 根据两个节点索引查找对应的全局棱边ID
     * @param node1 第一个节点索引
     * @param node2 第二个节点索引
     * @param edge_to_global_id 棱边映射表
     * @return int 全局棱边ID，若未找到返回-1
     *
     * @details 内部将节点对排序为(min, max)后在哈希表中查找。
     *          这是无向图的标准化表示方式。
     */
    static int findGlobalEdgeID(int node1, int node2,
                                const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id);

    /**
     * @brief 构建余树边集合（所有不在树中的边）
     * @param tree_edges 树边的全局棱边ID集合
     * @param edge_to_global_id 完整的棱边映射表
     * @return std::set<int> 余树边的全局棱边ID集合
     *
     * @details 遍历edge_to_global_id的所有条目，
     *          将不在tree_edges中的值收集到cotree_edges集合中。
     */
    static std::set<int> buildCotreeEdges(
        const std::set<int>& tree_edges,
        const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id
    );
};

} // namespace fe_em
