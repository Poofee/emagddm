/**
 * @file global_edge_id_generator.hpp
 * @brief 电磁物理层 - 全局棱边ID生成器头文件
 * @details 基于有限元网格拓扑数据，为所有矢量单元（VECTOR_EDGE_ONLY / MIXED_AV）
 *          的棱边分配全局唯一ID，构建「局部棱边ID → 全局棱边ID」映射表。
 *          是Nedelec棱边元DOF管理和FETI-DP棱边自由度编号的基础工具。
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 *
 * @note 算法原理：
 *       - 每条物理棱边由一对全局节点ID（升序）唯一标识
 *       - 遍历所有矢量单元的局部棱边，将局部节点对映射为全局节点对
 *       - 使用哈希表去重，首次出现的棱边分配新ID（从0连续编号）
 *       - 共享棱边的不同单元自动获得相同的全局棱边ID
 *
 * @note 与后续模块对接：
 *       - DOF管理模块：基于全局棱边ID分配Nedelec自由度编号
 *       - FETI-DP求解器：基于全局棱边ID构建拉格朗日乘子空间
 *       - 后处理模块：基于全局棱边ID输出切向场分量
 */

#pragma once

#include <vector>
#include <unordered_map>
#include "em_mesh_data.hpp"
#include "logger_factory.hpp"

namespace fe_em {

// ==================== 内部辅助：跨平台兼容的棱边键哈希函数 ====================

/**
 * @struct EdgeKeyHash
 * @brief 棱边键哈希函数对象（解决libc++未提供std::hash<std::tuple>特化的问题）
 * @details macOS的libc++标准库未提供std::tuple的哈希特化，
 *          使用std::pair配合自定义哈希函数实现跨平台兼容。
 *          哈希策略：将两个int组合为64位值，避免碰撞
 */
struct EdgeKeyHash {
    size_t operator()(const std::pair<int, int>& p) const {
        static constexpr size_t prime = 0x9e3779b97f4a7c15ULL;
        return (static_cast<size_t>(p.first) * prime) ^ static_cast<size_t>(p.second);
    }
};

/**
 * @class GlobalEdgeIDGenerator
 * @brief 全局棱边ID生成器
 * @details 遍历网格中所有矢量单元（VECTOR_EDGE_ONLY / MIXED_AV），
 *          为每条物理棱边分配唯一的全局ID，并记录每个单元的
 *          局部棱边到全局棱边的映射关系。
 *
 * @note 设计原则：
 *       - 仅处理矢量DOF类型的单元，标量单元对应位置留空
 *       - 全局棱边ID从0开始连续编号，便于数组索引
 *       - 使用unordered_map实现O(1)棱边查找和去重
 *       - 构造时仅保存引用，不拷贝网格数据（零开销）
 *
 * @note 典型使用流程：
 *       @code
 *       GlobalEdgeIDGenerator generator(mesh_data);
 *       generator.generate();
 *       int num_edges = generator.getNumGlobalEdges();
 *       const auto& mapping = generator.getElemLocalToGlobalEdge();
 *       // mapping[elem_idx][local_edge_idx] = global_edge_id
 *       @endcode
 */
class GlobalEdgeIDGenerator {
public:
    /**
     * @brief 构造函数
     * @param mesh_data 有限元网格数据的常量引用（仅读取，不修改）
     *
     * @note 构造函数仅保存引用，不执行任何生成操作。
     *       实际的ID生成需显式调用 generate() 方法。
     */
    explicit GlobalEdgeIDGenerator(const EMMeshData& mesh_data);

    /**
     * @brief 执行全局棱边ID生成
     * @throws std::invalid_argument 当遇到不支持的单元类型时抛出
     *
     * @details 核心算法步骤：
     *          1. 创建 unordered_map 存储唯一棱边（键：排序后的全局节点对）
     *          2. 遍历所有单元，仅处理 VECTOR_EDGE_ONLY / MIXED_AV 类型
     *          3. 对每条局部棱边：局部节点ID → 全局节点ID → 排序为(min, max)
     *          4. 哈希表中不存在则分配新全局ID（从0递增）
     *          5. 记录 局部棱边ID → 全局棱边ID 映射到内部表
     *          6. 输出统计日志（总单元数、矢量单元数、全局棱边数）
     *
     * @warning 调用前确保 mesh_data_ 已正确填充节点和单元数据
     * @note 此方法可重复调用，每次调用会重置之前的结果
     */
    void generate();

    /**
     * @brief 获取全局棱边总数
     * @return int 全局唯一棱边的数量（从0到num_global_edges_-1编号）
     *
     * @note 必须在 generate() 调用后使用，否则返回0
     */
    int getNumGlobalEdges() const;

    /**
     * @brief 获取每个单元的局部棱边到全局棱边的映射表
     * @return const std::vector<std::vector<int>>& 二维向量引用
     *         - 外层维度：按单元索引（size == elements.size()）
     *         - 内层维度：按局部棱边索引（size == 该单元的局部棱边数）
     *         - 值：对应的全局棱边ID（从0开始）
     *
     * @details 返回值语义说明：
     *          - 矢量单元（VECTOR_EDGE_ONLY / MIXED_AV）：内层vector非空，
     *            每个元素为该单元第i条局部棱边对应的全局棱边ID
     *          - 标量单元（SCALAR_ONLY）：内层vector为空（size == 0）
     *
     * @note 必须在 generate() 调用后使用
     * @note 返回值为const引用，外部不可修改
     */
    const std::vector<std::vector<int>>& getElemLocalToGlobalEdge() const;

    /**
     * @brief 根据两个节点ID查找对应的全局棱边ID（O(1)哈希查找）
     * @param node1 棱边的第一个全局节点ID
     * @param node2 棱边的第二个全局节点ID
     * @return int 对应的全局棱边ID（从0开始），若未找到返回-1
     *
     * @details 内部将节点对排序为(min, max)后在持久化哈希表中查找。
     *          此方法为EMDOFManager的矢量边界约束标记提供高效查找能力，
     *          避免了O(N_elements × N_edges)的遍历开销。
     *
     * @note 必须在 generate() 调用后使用，否则返回-1
     * @note 节点顺序无关：(node1=3, node2=5) 与 (node1=5, node2=3) 返回相同结果
     *
     * @code
     * // 典型用法：在DOF管理器中快速定位Dirichlet棱边约束
     * int global_edge_id = generator.getGlobalEdgeID(constrained_node_a, constrained_node_b);
     * if (global_edge_id >= 0) {
     *     // 找到对应的自由度，执行约束标记
     * }
     * @endcode
     */
    int getGlobalEdgeID(int node1, int node2) const;

private:
    const EMMeshData& mesh_data_;                              ///< 网格数据常量引用
    int num_global_edges_ = 0;                                 ///< 全局棱边总数
    std::vector<std::vector<int>> elem_local_to_global_edge_;  ///< 单元局部→全局棱边映射表
    std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> edge_map_;  ///< 持久化棱边哈希表（支持O(1)反向查找）
};

} // namespace fe_em
