/**
 * @file global_edge_id_generator.cpp
 * @brief 电磁物理层 - 全局棱边ID生成器实现
 * @details 实现基于网格拓扑的全局棱边唯一ID分配算法，
 *          为Nedelec棱边元的自由度管理提供基础编号服务。
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "global_edge_id_generator.hpp"
#include "element_geometry.hpp"
#include <unordered_map>
#include <algorithm>

namespace fe_em {

// ==================== 构造与析构 ====================

GlobalEdgeIDGenerator::GlobalEdgeIDGenerator(const EMMeshData& mesh_data)
    : mesh_data_(mesh_data)
    , num_global_edges_(0)
{
    FEEM_DEBUG("GlobalEdgeIDGenerator: 构造完成, 网格包含{}个单元, {}个节点",
               mesh_data_.elements.size(), mesh_data_.nodes.size());
}

// ==================== 核心算法：全局棱边ID生成 ====================

void GlobalEdgeIDGenerator::generate()
{
    // 重置内部状态（支持重复调用）
    num_global_edges_ = 0;
    elem_local_to_global_edge_.clear();
    elem_local_to_global_edge_.resize(mesh_data_.elements.size());

    // 哈希表：键为排序后的全局节点对(min_node_id, max_node_id)，值为全局棱边ID
    // 使用成员变量edge_map_持久化存储，支持后续O(1)反向查找（getGlobalEdgeID接口）
    edge_map_.clear();
    edge_map_.reserve(mesh_data_.elements.size() * 6);

    int vector_elem_count = 0;  // 矢量DOF单元计数（用于日志统计）

    // 遍历所有单元，构建局部→全局棱边映射
    for (size_t elem_idx = 0; elem_idx < mesh_data_.elements.size(); ++elem_idx) {
        const auto& elem = mesh_data_.elements[elem_idx];

        // 仅处理矢量DOF类型的单元（Nedelec棱边元需要棱边自由度）
        if (elem.dof_type != DOFType::VECTOR_EDGE_ONLY &&
            elem.dof_type != DOFType::MIXED_AV) {
            // 标量单元不涉及棱边自由度，对应位置保留空vector
            continue;
        }

        ++vector_elem_count;

        // 获取当前单元类型的局部棱边定义（每条棱边为局部节点ID对）
        std::vector<std::tuple<int, int>> local_edges;
        try {
            local_edges = ElementGeometry::get_local_edges(elem.elem_type);
        } catch (const std::invalid_argument& e) {
            FEEM_ERROR("GlobalEdgeIDGenerator: 单元[{}]类型[{}]不支持棱边查询: {}",
                       elem.id, static_cast<int>(elem.elem_type), e.what());
            throw std::invalid_argument(
                "GlobalEdgeIDGenerator: 不支持的单元类型 (elem_id=" +
                std::to_string(elem.id) + ")");
        }

        // 为当前单元预分配局部→全局映射向量
        elem_local_to_global_edge_[elem_idx].resize(local_edges.size());

        // 遍历当前单元的每条局部棱边，建立局部→全局映射
        for (size_t local_edge_idx = 0; local_edge_idx < local_edges.size(); ++local_edge_idx) {
            const auto& [local_node_a, local_node_b] = local_edges[local_edge_idx];

            // 将局部节点ID转换为全局节点ID
            // node_ids存储的是全局节点ID列表，local_node_a/b是索引位置
            int global_node_a = elem.node_ids[local_node_a];
            int global_node_b = elem.node_ids[local_node_b];

            // 将全局节点对排序为升序pair作为唯一键
            // 排序确保(1,5)和(5,1)映射到同一条物理棱边
            auto edge_key = std::make_pair(
                std::min(global_node_a, global_node_b),
                std::max(global_node_a, global_node_b));

            // 在哈希表中查找该棱边是否已分配全局ID
            auto it = edge_map_.find(edge_key);
            if (it == edge_map_.end()) {
                // 首次出现：分配新的全局棱边ID（从0连续递增）
                int new_global_id = static_cast<int>(edge_map_.size());
                edge_map_[edge_key] = new_global_id;
                elem_local_to_global_edge_[elem_idx][local_edge_idx] = new_global_id;
            } else {
                // 已存在：复用已有的全局棱边ID（共享棱边场景）
                elem_local_to_global_edge_[elem_idx][local_edge_idx] = it->second;
            }
        }
    }

    // 记录最终统计结果
    num_global_edges_ = static_cast<int>(edge_map_.size());

    // 输出统计信息日志
    FEEM_INFO("GlobalEdgeIDGenerator: 棱边ID生成完成");
    FEEM_INFO("  - 总单元数: {}", mesh_data_.elements.size());
    FEEM_INFO("  - 矢量单元数(VECTOR_EDGE_ONLY/MIXED_AV): {}", vector_elem_count);
    FEEM_INFO("  - 全局唯一棱边数: {}", num_global_edges_);
    FEEM_INFO("  - 哈希表负载因子: {:.2f}", edge_map_.load_factor());
}

// ==================== 查询接口 ====================

int GlobalEdgeIDGenerator::getNumGlobalEdges() const
{
    return num_global_edges_;
}

const std::vector<std::vector<int>>& GlobalEdgeIDGenerator::getElemLocalToGlobalEdge() const
{
    return elem_local_to_global_edge_;
}

int GlobalEdgeIDGenerator::getGlobalEdgeID(int node1, int node2) const
{
    // 将节点对排序为升序（与generate()中的键格式一致）
    auto edge_key = std::make_pair(
        std::min(node1, node2),
        std::max(node1, node2));

    // O(1)哈希查找
    auto it = edge_map_.find(edge_key);
    if (it != edge_map_.end()) {
        return it->second;
    }

    // 未找到：该棱边可能不属于任何矢量单元（如纯标量网格中的虚拟棱边）
    FEEM_DEBUG("getGlobalEdgeID: 未找到棱边({}, {})对应的全局ID（可能属于标量单元）", node1, node2);
    return -1;
}

} // namespace fe_em
