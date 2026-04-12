/**
 * @file mesh_graph_extractor.cpp
 * @brief 电磁物理层 - 网格图结构提取器实现文件
 * @details 实现从EMMeshData中提取无向图结构的核心算法。
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "mesh_graph_extractor.hpp"
#include "element_geometry.hpp"
#include <set>
#include <stdexcept>

namespace fe_em {

// ==================== 构造函数 ====================

MeshGraphExtractor::MeshGraphExtractor(
    const EMMeshData& mesh_data,
    const std::vector<std::vector<int>>& elem_local_to_global_edge
) : mesh_data_(mesh_data),
    elem_local_to_global_edge_(elem_local_to_global_edge) {
}

// ==================== 公共接口：extract() ====================

GraphResult MeshGraphExtractor::extract() {
    FEEM_INFO("========== MeshGraphExtractor: 开始提取图结构 ==========");

    // ---------- 输入验证 ----------
    if (mesh_data_.elements.empty()) {
        FEEM_WARN("网格数据为空（0个单元），返回空图结构");
        return GraphResult{};
    }

    if (elem_local_to_global_edge_.size() != mesh_data_.elements.size()) {
        throw std::invalid_argument(
            "MeshGraphExtractor::extract(): elem_local_to_global_edge尺寸(" +
            std::to_string(elem_local_to_global_edge_.size()) +
            ")与单元数(" + std::to_string(mesh_data_.elements.size()) + ")不匹配"
        );
    }

    // ========== 第一阶段：收集节点并建立索引映射 ==========
    FEEM_DEBUG("第一阶段：收集矢量单元节点ID");
    std::vector<int> global_node_ids = collectVectorNodes();

    if (global_node_ids.empty()) {
        FEEM_WARN("未找到任何矢量单元（VECTOR_EDGE_ONLY/MIXED_AV），返回空图结构");
        return GraphResult{};
    }

    GraphResult result;
    buildNodeIndexMapping(global_node_ids, result.global_node_to_index, result.index_to_global_node);
    result.num_nodes = static_cast<int>(result.index_to_global_node.size());

    FEEM_INFO("映射信息: 全局节点数={}, 紧引索引范围=[0, {}]",
              global_node_ids.size(), result.num_nodes - 1);

    // ========== 第二阶段：构建邻接表和棱边映射 ==========
    FEEM_DEBUG("第二阶段：构建邻接表和棱边映射");

    // 预分配邻接表空间（每个节点一个空vector）
    result.adjacency_list.resize(result.num_nodes);

    // 统计信息
    int total_elements = static_cast<int>(mesh_data_.elements.size());
    int vector_element_count = 0;
    int total_local_edges_processed = 0;

    // 遍历所有单元
    for (int elem_idx = 0; elem_idx < total_elements; ++elem_idx) {
        const auto& element = mesh_data_.elements[elem_idx];

        // 跳过非矢量单元
        if (!isVectorElementType(element.dof_type)) {
            continue;
        }
        vector_element_count++;

        // 检查该单元是否有全局棱边映射（标量单元的映射为空）
        if (elem_local_to_global_edge_[elem_idx].empty()) {
            FEEM_DEBUG("单元{}({})的局部→全局棱边映射为空，跳过",
                      elem_idx, dofTypeToString(element.dof_type));
            continue;
        }

        // 获取单元的局部棱边定义
        std::vector<std::tuple<int, int>> local_edges;
        try {
            local_edges = ElementGeometry::get_local_edges(element.elem_type);
        } catch (const std::invalid_argument& e) {
            FEEM_ERROR("单元{}: 不支持的单元类型{} - {}",
                      elem_idx, static_cast<int>(element.elem_type), e.what());
            throw;
        }

        // 遍历每条局部棱边
        const size_t num_local_edges = local_edges.size();
        for (size_t local_edge_idx = 0; local_edge_idx < num_local_edges; ++local_edge_idx) {
            const auto& [local_node1, local_node2] = local_edges[local_edge_idx];

            // 局部节点ID → 全局节点ID
            int global_node1 = element.node_ids[local_node1];
            int global_node2 = element.node_ids[local_node2];

            // 全局节点ID → 紧凑索引
            auto it1 = result.global_node_to_index.find(global_node1);
            auto it2 = result.global_node_to_index.find(global_node2);

            // 安全检查：节点应在映射中（因为是从同一批单元收集的）
            if (it1 == result.global_node_to_index.end() || it2 == result.global_node_to_index.end()) {
                FEEM_WARN("单元{}的局部棱边{}包含未映射的节点({}, {})",
                         elem_idx, local_edge_idx, global_node1, global_node2);
                continue;
            }

            int index1 = it1->second;
            int index2 = it2->second;

            // 构建排序后的节点对（作为无向图的键）
            std::pair<int, int> edge_key = (index1 < index2)
                                          ? std::make_pair(index1, index2)
                                          : std::make_pair(index2, index1);

            // 获取全局棱边ID
            int global_edge_id = elem_local_to_global_edge_[elem_idx][local_edge_idx];

            // 检查是否已存在该棱边（去重）
            if (result.edge_to_global_id.find(edge_key) == result.edge_to_global_id.end()) {
                // 首次出现：记录棱边映射
                result.edge_to_global_id[edge_key] = global_edge_id;

                // 双向添加到邻接表（无向图）
                result.adjacency_list[index1].push_back(index2);
                result.adjacency_list[index2].push_back(index1);
            }

            total_local_edges_processed++;
        }
    }

    result.num_edges = static_cast<int>(result.edge_to_global_id.size());

    // ========== 第三阶段：后处理 ==========
    FEEM_DEBUG("第三阶段：后处理（排序邻居列表）");

    // 对每个节点的邻居列表排序（便于调试和确定性输出）
    for (auto& neighbors : result.adjacency_list) {
        std::sort(neighbors.begin(), neighbors.end());
    }

    // ========== 日志输出 ==========
    FEEM_INFO("========== MeshGraphExtractor: 图结构提取完成 ==========");
    FEEM_INFO("输入统计: 总单元数={}, 矢量单元数={}", total_elements, vector_element_count);
    FEEM_INFO("图结构:   节点数={}, 唯一棱边数={}", result.num_nodes, result.num_edges);
    FEEM_INFO("处理统计: 局部棱边处理总数={}", total_local_edges_processed);

    // 计算平均度数（用于性能分析）
    if (result.num_nodes > 0) {
        double avg_degree = 0.0;
        for (const auto& neighbors : result.adjacency_list) {
            avg_degree += neighbors.size();
        }
        avg_degree /= result.num_nodes;
        FEEM_INFO("拓扑特征: 平均节点度数={:.2f}", avg_degree);
    }

    FEEM_INFO("======================================================");

    return result;
}

// ==================== 私有辅助方法 ====================

std::vector<int> MeshGraphExtractor::collectVectorNodes() const {
    std::set<int> node_set;

    for (const auto& element : mesh_data_.elements) {
        if (isVectorElementType(element.dof_type)) {
            for (int node_id : element.node_ids) {
                node_set.insert(node_id);
            }
        }
    }

    return std::vector<int>(node_set.begin(), node_set.end());
}

void MeshGraphExtractor::buildNodeIndexMapping(
    const std::vector<int>& global_node_ids,
    std::unordered_map<int, int>& out_global_to_index,
    std::vector<int>& out_index_to_global
) const {
    out_index_to_global.resize(global_node_ids.size());

    for (size_t i = 0; i < global_node_ids.size(); ++i) {
        int global_id = global_node_ids[i];
        int index = static_cast<int>(i);

        out_global_to_index[global_id] = index;
        out_index_to_global[index] = global_id;
    }
}

bool MeshGraphExtractor::isVectorElementType(DOFType dof_type) {
    return dof_type == DOFType::VECTOR_EDGE_ONLY || dof_type == DOFType::MIXED_AV;
}

} // namespace fe_em
