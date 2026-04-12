/**
 * @file mesh_graph_extractor.hpp
 * @brief 电磁物理层 - 网格图结构提取器头文件
 * @details 从EMMeshData中提取无向图结构，用于Tree Gauge生成树构建。
 *          仅处理矢量单元（VECTOR_EDGE_ONLY / MIXED_AV），忽略纯标量单元。
 *
 * @note 核心功能：
 *       - 构建节点邻接表（adjacency_list）
 *       - 建立棱边到全局ID的映射（edge_to_global_id）
 *       - 支持非连续节点ID的紧凑索引映射
 *
 * @note 与后续模块对接：
 *       - SpanningTreeBuilder：使用邻接表和棱边映射构建生成树
 *       - TreeGauge：基于图结构执行树-余树分类
 *       - MultiConnectedDomainDetector：基于图拓扑检测多连通域
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "em_mesh_data.hpp"
#include "global_edge_id_generator.hpp"
#include "logger_factory.hpp"

namespace fe_em {

// ==================== 输出数据结构 ====================

/**
 * @struct GraphResult
 * @brief 图结构提取结果数据结构
 * @details 包含Tree Gauge算法所需的完整图论信息：
 *          - 节点邻接表：用于BFS/DFS遍历
 *          - 棱边映射表：用于树边/余树边分类
 *          - 统计信息：节点数、棱边数
 *
 * @note 数据结构设计原则：
 *       - adjacency_list使用vector<vector<int>>（缓存友好，优于set）
 *       - edge_to_global_id使用unordered_map（O(1)查找性能）
 *       - 所有ID均为全局ID（与EMMeshData中的原始ID一致）
 */
struct GraphResult {
    std::vector<std::vector<int>> adjacency_list;                                     ///< 节点邻接表（adjacency_list[node_index] = 相邻节点的node_index列表）
    std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> edge_to_global_id;      ///< 无序节点对 → 全局棱边ID的映射（键为排序后的节点对）
    std::unordered_map<int, int> global_node_to_index;                                ///< 全局节点ID → 紧凑索引的映射（处理非连续ID）
    std::vector<int> index_to_global_node;                                            ///< 紧引索引 → 全局节点ID的反向映射
    int num_nodes = 0;                                                                ///< 图中总节点数（参与矢量单元的节点）
    int num_edges = 0;                                                                ///< 图中总棱边数（去重后的唯一物理棱边）

    /**
     * @brief 检查结果是否有效
     * @return bool 如果包含有效数据返回true
     */
    bool isValid() const {
        return num_nodes > 0 && num_edges >= 0 && !adjacency_list.empty();
    }
};

// ==================== 核心工具类 ====================

/**
 * @class MeshGraphExtractor
 * @brief 网格图结构提取器
 * @details 从有限元网格数据中提取用于Tree Gauge计算的无向图结构。
 *
 * @note 算法流程：
 *       1. 遍历所有矢量单元（VECTOR_EDGE_ONLY / MIXED_AV）
 *       2. 对每个单元的每条局部棱边，获取两个端点的全局节点ID
 *       3. 将全局节点ID映射为紧凑连续索引（处理非连续ID情况）
 *       4. 构建邻接表（避免重复添加邻居节点）
 *       5. 记录每条唯一棱边对应的全局棱边ID（来自GlobalEdgeIDGenerator）
 *
 * @note 设计原则：
 *       - 仅处理矢量DOF类型单元，标量单元自动跳过
 *       - 构造时仅保存引用，不拷贝网格数据（零开销）
 *       - 使用哈希表实现O(1)棱边去重和查找
 *       - 支持百万级节点的大规模网格
 *
 * @note 典型使用流程：
 *       @code
 *       // 第一步：生成全局棱边ID
 *       GlobalEdgeIDGenerator edge_gen(mesh_data);
 *       edge_gen.generate();
 *       const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();
 *
 *       // 第二步：提取图结构
 *       MeshGraphExtractor extractor(mesh_data, edge_mapping);
 *       GraphResult result = extractor.extract();
 *
 *       // 第三步：使用结果
 *       int n = result.num_nodes;
 *       int e = result.num_edges;
 *       const auto& adj = result.adjacency_list;
 *       @endcode
 *
 * @note 性能特征：
 *       - 时间复杂度：O(E_local)，E_local为所有矢量单元的局部棱边总数
 *       - 空间复杂度：O(V + E)，V为节点数，E为唯一棱边数
 *       - 内存优化：预分配容器容量，避免动态扩容
 */
class MeshGraphExtractor {
public:
    /**
     * @brief 构造函数
     * @param mesh_data 有限元网格数据的常量引用（仅读取，不修改）
     * @param elem_local_to_global_edge 全局棱边ID映射表（来自GlobalEdgeIDGenerator）
     *
     * @note 构造函数仅保存引用，不执行任何提取操作。
     *       实际的图提取需显式调用 extract() 方法。
     *
     * @warning 调用 extract() 前，elem_local_to_global_edge 必须已由
     *          GlobalEdgeIDGenerator::generate() 正确填充
     *
     * @throws 无（构造函数不执行实质性操作）
     */
    explicit MeshGraphExtractor(
        const EMMeshData& mesh_data,
        const std::vector<std::vector<int>>& elem_local_to_global_edge
    );

    /**
     * @brief 执行图结构提取
     * @return GraphResult 包含邻接表、棱边映射、统计信息的完整结果
     * @throws std::invalid_argument 当网格数据为空或遇到不支持的单元类型时抛出
     *
     * @details 核心算法步骤：
     *          1. **第一阶段：收集节点ID**
     *             - 遍历所有矢量单元，收集出现的全局节点ID
     *             - 建立全局ID → 紧凑索引的双向映射
     *             - 处理节点ID可能不连续的情况
     *
     *          2. **第二阶段：构建邻接表和棱边映射**
     *             - 初始化邻接表（num_nodes × 空vector）
     *             - 对每个矢量单元的每条局部棱边：
     *               a) 获取局部节点ID → 全局节点ID → 索引
     *               b) 构建排序后的节点对 (min, max)
     *               c) 若该棱边首次出现：
     *                  - 记录到 edge_to_global_id（值为全局棱边ID）
     *                  - 双向添加到邻接表（无向图）
     *               d) 若该棱边已存在：跳过（避免重复）
     *
     *          3. **第三阶段：后处理**
     *             - 对每个节点的邻居列表排序（便于调试和确定性输出）
     *             - 统计并记录日志信息
     *
     * @warning 调用前确保 mesh_data_ 和 elem_local_to_global_edge_ 已正确填充
     * @note 此方法可重复调用，每次调用会重新计算并返回新结果
     *
     * @note 日志输出示例：
     *       [INFO] ========== MeshGraphExtractor: 图结构提取完成 ==========
     *       [INFO] 输入统计: 总单元数=100, 矢量单元数=80
     *       [INFO] 图结构:   节点数=500, 唯一棱边数=1200
     *       [INFO] 映射信息: 全局ID范围=[0, 520], 紧凑索引范围=[0, 499]
     *       [INFO] =======================================================
     */
    GraphResult extract();

private:
    const EMMeshData& mesh_data_;                              ///< 网格数据常量引用
    const std::vector<std::vector<int>>& elem_local_to_global_edge_;  ///< 全局棱边ID映射表引用

    /**
     * @brief 收集所有矢量单元中出现的全局节点ID
     * @return std::vector<int> 去重后的全局节点ID列表（未排序）
     * @details 遍历所有VECTOR_EDGE_ONLY和MIXED_AV类型单元，
     *          收集其node_ids字段中的所有全局节点ID，
     *          使用临时set去重后转为vector返回
     */
    std::vector<int> collectVectorNodes() const;

    /**
     * @brief 构建全局节点ID到紧凑索引的双向映射
     * @param global_node_ids 去重后的全局节点ID列表
     * @param out_global_to_index 输出：全局ID → 索引的映射
     * @param out_index_to_global 输出：索引 → 全局ID的映射
     * @details 将可能不连续的全局节点ID映射为[0, N-1]范围的紧凑索引，
     *          便于数组索引访问邻接表
     */
    void buildNodeIndexMapping(
        const std::vector<int>& global_node_ids,
        std::unordered_map<int, int>& out_global_to_index,
        std::vector<int>& out_index_to_global
    ) const;

    /**
     * @brief 检查单元是否为矢量DOF类型
     * @param dof_type 单元的DOFType枚举值
     * @return bool 如果是VECTOR_EDGE_ONLY或MIXED_AV返回true
     * @details 辅助方法，用于过滤需要处理的单元类型
     */
    static bool isVectorElementType(DOFType dof_type);
};

} // namespace fe_em
