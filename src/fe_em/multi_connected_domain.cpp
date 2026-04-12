/**
 * @file multi_connected_domain.cpp
 * @brief 电磁物理层 - 多连通域检测器实现文件
 * @details 实现基于网格拓扑的多连通域检测算法，识别孔洞并选择环量自由度。
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "multi_connected_domain.hpp"
#include "spanning_tree.hpp"
#include <algorithm>
#include <stdexcept>
#include <set>

namespace fe_em {

// ==================== 公共接口：detect() ====================

DetectionResult MultiConnectedDomainDetector::detect(
    const EMMeshData& mesh_data,
    const SpanningTreeResult& tree_result
) {
    FEEM_INFO("========== MultiConnectedDomainDetector: 开始多连通域检测 ==========");

    // ---------- 第一步：输入验证 ----------
    validateInput(mesh_data, tree_result);

    // ---------- 第二步：边界信息分析 ----------
    std::set<int> boundary_types;
    int inner_boundary_count = 0;
    analyzeBoundaryMarkers(mesh_data, boundary_types, inner_boundary_count);

    FEEM_DEBUG("边界分析完成: 边界类型数={}, 潜在内部边界数={}",
               boundary_types.size(), inner_boundary_count);

    // ---------- 第三步：环量自由度选择 ----------
    std::set<int> circulation_dofs = selectCirculationDOFs(
        tree_result.cotree_edges,
        tree_result.tree_edges,
        tree_result.node_parent,
        inner_boundary_count
    );

    // ---------- 第四步：构建孔洞回路 ----------
    std::vector<std::vector<int>> hole_loop_edges = buildHoleLoopEdges(
        circulation_dofs,
        tree_result.tree_edges,
        tree_result.node_parent
    );

    // ---------- 第五步：组装结果 ----------
    DetectionResult result;
    result.num_holes = static_cast<int>(circulation_dofs.size());
    result.hole_loop_edges = std::move(hole_loop_edges);
    result.circulation_dofs = std::move(circulation_dofs);

    // ---------- 第六步：结果验证 ----------
    if (!validateResult(result, tree_result.cotree_edges)) {
        FEEM_WARN("多连通域检测结果一致性验证失败，返回空结果");
        return DetectionResult{};
    }

    // ---------- 日志输出 ----------
    FEEM_INFO("========== MultiConnectedDomainDetector: 多连通域检测完成 ==========");
    FEEM_INFO("输入统计: 节点数={}, 单元数={}, 边界标记数={}",
              mesh_data.getNodeCount(), mesh_data.getElementCount(),
              mesh_data.getBoundaryMarkerCount());
    FEEM_INFO("生成树信息: 树边数={}, 余树边数={}",
              tree_result.tree_edges.size(), tree_result.cotree_edges.size());
    FEEM_INFO("检测结果: 孔洞数={}, 环量自由度数={}",
              result.num_holes, result.circulation_dofs.size());

    if (result.num_holes > 0) {
        std::string circulation_str;
        for (int edge_id : result.circulation_dofs) {
            circulation_str += std::to_string(edge_id) + " ";
        }
        FEEM_INFO("环量边: [{}]", circulation_str);
    }

    FEEM_INFO("============================================================");

    return result;
}

// ==================== 私有方法：输入验证 ====================

void MultiConnectedDomainDetector::validateInput(
    const EMMeshData& mesh_data,
    const SpanningTreeResult& tree_result
) {
    // 验证网格数据不为空
    if (mesh_data.nodes.empty()) {
        throw std::invalid_argument(
            "MultiConnectedDomainDetector::detect(): 网格数据为空（无节点）"
        );
    }

    if (mesh_data.elements.empty()) {
        throw std::invalid_argument(
            "MultiConnectedDomainDetector::detect(): 网格数据为空（无单元）"
        );
    }

    // 验证生成树结果有效
    if (!tree_result.isValid()) {
        throw std::invalid_argument(
            "MultiConnectedDomainDetector::detect(): 生成树结果无效"
        );
    }

    // 验证余树边集合非空（否则无法选择环量自由度）
    if (tree_result.cotree_edges.empty()) {
        FEEM_WARN("余树边集合为空，无法选择环量自由度（可能为单连通域或退化情况）");
    }

    FEEM_DEBUG("输入验证通过: 节点数={}, 单元数={}, 树边数={}, 余树边数={}",
               mesh_data.getNodeCount(), mesh_data.getElementCount(),
               tree_result.tree_edges.size(), tree_result.cotree_edges.size());
}

// ==================== 私有方法：边界标记分析 ====================

void MultiConnectedDomainDetector::analyzeBoundaryMarkers(
    const EMMeshData& mesh_data,
    std::set<int>& out_boundary_types,
    int& out_inner_boundary_count
) {
    out_boundary_types.clear();
    out_inner_boundary_count = 0;

    // 如果没有边界标记，直接返回（单连通域）
    if (mesh_data.boundary_markers.empty()) {
        FEEM_DEBUG("无边界标记信息，假设为单连通域");
        return;
    }

    // 遍历所有边界标记，收集类型和潜在孔洞信息
    for (const auto& marker : mesh_data.boundary_markers) {
        out_boundary_types.insert(static_cast<int>(marker.bnd_type));

        // 启发式规则1：检查边界名称是否包含"hole"、"inner"、"cavity"等关键词
        std::string name_lower = marker.name;
        std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                      [](unsigned char c) { return std::tolower(c); });

        if (name_lower.find("hole") != std::string::npos ||
            name_lower.find("inner") != std::string::npos ||
            name_lower.find("cavity") != std::string::npos ||
            name_lower.find("孔") != std::string::npos) {
            out_inner_boundary_count++;
            FEEM_DEBUG("发现潜在内部边界: ID={}, 名称={}, 类型={}",
                       marker.id, marker.name,
                       static_cast<int>(marker.bnd_type));
        }

        // 启发式规则2：PERFECT_E/PERFECT_H类型的边界可能是内部导体边界
        if (marker.bnd_type == BndType::PERFECT_E ||
            marker.bnd_type == BndType::PERFECT_H) {
            // 仅当目标实体数量较少时才认为是内部边界（避免误判外部边界）
            bool is_small_target = false;
            if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                const auto& target_vec = std::get<std::vector<int>>(marker.target_ids);
                if (static_cast<size_t>(target_vec.size()) < static_cast<size_t>(mesh_data.getNodeCount() / 2)) {
                    is_small_target = true;
                }
            } else if (std::holds_alternative<std::vector<std::vector<int>>>(marker.target_ids)) {
                const auto& target_vec = std::get<std::vector<std::vector<int>>>(marker.target_ids);
                if (static_cast<size_t>(target_vec.size()) < static_cast<size_t>(mesh_data.getElementCount() / 2)) {
                    is_small_target = true;
                }
            }

            if (is_small_target && name_lower.find("hole") == std::string::npos) {
                out_inner_boundary_count++;
                FEEM_DEBUG("发现潜在内部导体边界: ID={}, 名称={}", marker.id, marker.name);
            }
        }
    }

    FEEM_INFO("边界分析统计: 总边界标记数={}, 不同类型数={}, 潜在孔洞数={}",
              mesh_data.boundary_markers.size(), out_boundary_types.size(),
              out_inner_boundary_count);
}

// ==================== 私有方法：环量自由度选择 ====================

std::set<int> MultiConnectedDomainDetector::selectCirculationDOFs(
    const std::set<int>& cotree_edges,
    const std::set<int>& /*tree_edges*/,
    const std::unordered_map<int, int>& /*node_parent*/,
    int estimated_hole_count
) {
    std::set<int> selected_dofs;

    // 情况1：估计的孔洞数为0 → 返回空集（单连通域）
    if (estimated_hole_count <= 0) {
        FEEM_INFO("估计孔洞数=0，不选择环量自由度（单连通域）");
        return selected_dofs;
    }

    // 情况2：有余树边但无法精确识别孔洞 → 使用简化策略
    if (cotree_edges.empty()) {
        FEEM_WARN("余树边为空，无法选择环量自由度");
        return selected_dofs;
    }

    // 实际可选数量不超过余树边数量和估计孔洞数的较小值
    int num_to_select = std::min(estimated_hole_count, static_cast<int>(cotree_edges.size()));

    FEEM_INFO("开始选择环量自由度: 目标数量={}, 可用余树边数={}",
              num_to_select, cotree_edges.size());

    // 将cotree_edges转为向量以便索引访问
    std::vector<int> cotree_vector(cotree_edges.begin(), cotree_edges.end());

    // Phase 1简化策略：贪心选择前N条余树边
    // 选择标准：优先选择ID较大或分布较均匀的边（简单启发式）
    for (int i = 0; i < num_to_select; ++i) {
        // 使用间隔采样策略，使选中的边分布更均匀
        int index = (i * static_cast<int>(cotree_vector.size())) / num_to_select;
        index = std::min(index, static_cast<int>(cotree_vector.size()) - 1);

        selected_dofs.insert(cotree_vector[index]);
    }

    // TODO Phase 2优化方向：
    // 1. 基于基本圈基（fundamental cycles）的精确选择：
    //    - 对每条余树边(u,v)，计算其在生成树中的基本圈
    //    - 分析基本圈是否包围某个孔洞（需要几何信息辅助）
    //    - 选择真正与孔洞相关的独立回路
    //
    // 2. 基于对偶图的方法：
    //    - 构建网格的对偶图（面→节点，邻接面→边）
    //    - 在对偶图中寻找独立回路对应原问题的孔洞
    //    - 对偶图中的生成树余树边即为原问题的环量自由度
    //
    // 3. 代数拓扑方法（Betti数）：
    //    - 计算网格的第一Betti数 β₁ = E - V + C（C=连通分量数）
    //    - β₁ 即为独立的孔洞数（亏格genus）
    //    - 构建同调群的基向量作为环量自由度

    FEEM_DEBUG("环量自由度选择完成: 选中{}条边", selected_dofs.size());

    return selected_dofs;
}

// ==================== 私有方法：构建孔洞回路 ====================

std::vector<std::vector<int>> MultiConnectedDomainDetector::buildHoleLoopEdges(
    const std::set<int>& circulation_dofs,
    const std::set<int>& /*tree_edges*/,
    const std::unordered_map<int, int>& /*node_parent*/,
    const std::vector<std::vector<int>>& /*adjacency_list*/
) {
    std::vector<std::vector<int>> hole_loops;

    // Phase 1简化实现：每个环量自由度单独构成一个"回路"
    // 实际上这只是占位符，完整的回路应该包含基本圈的所有棱边
    for (int edge_id : circulation_dofs) {
        hole_loops.push_back({edge_id});
    }

    // TODO Phase 2完整实现：
    // 对于每个环量边 edge_id = (u, v)：
    // 1. 在生成树中查找从u到v的唯一路径 P（使用node_parent反向追踪）
    // 2. 基本圈 = P ∪ {(u, v)} （路径加上这条余树边构成闭合回路）
    // 3. 将基本圈中的所有全局棱边ID收集到hole_loops[i]中
    //
    // 示例代码框架：
    // for (int edge_id : circulation_dofs) {
    //     // 假设我们能获取edge_id的两个端点节点
    //     int node_u = getEdgeEndpoint(edge_id, 0);
    //     int node_v = getEdgeEndpoint(edge_id, 1);
    //
    //     // 在生成树中追踪从node_v到node_u的路径
    //     std::vector<int> path;
    //     int current = node_v;
    //     while (current != node_u && current != -1) {
    //         int parent = node_parent.at(current);
    //         if (parent == -1) break;
    //         // 查找边(current, parent)的全局ID并加入path
    //         int path_edge = findGlobalEdgeID(current, parent, edge_to_global_id);
    //         path.push_back(path_edge);
    //         current = parent;
    //     }
    //
    //     // 构建完整的基本圈
    //     std::vector<int> cycle = path;
    //     cycle.push_back(edge_id);  // 加上余树边本身
    //     hole_loops.push_back(cycle);
    // }

    return hole_loops;
}

// ==================== 私有方法：结果验证 ====================

bool MultiConnectedDomainDetector::validateResult(
    const DetectionResult& result,
    const std::set<int>& cotree_edges
) {
    // 验证1：孔洞数非负
    if (result.num_holes < 0) {
        FEEM_ERROR("验证失败: 孔洞数为负 ({})", result.num_holes);
        return false;
    }

    // 验证2：回路列表大小与孔洞数一致
    if (static_cast<int>(result.hole_loop_edges.size()) != result.num_holes) {
        FEEM_ERROR("验证失败: 回路列表大小({})与孔洞数({})不一致",
                   result.hole_loop_edges.size(), result.num_holes);
        return false;
    }

    // 验证3：每个回路至少包含一条棱边
    for (int i = 0; i < result.num_holes; ++i) {
        if (result.hole_loop_edges[i].empty()) {
            FEEM_ERROR("验证失败: 孔洞{}的回路为空", i);
            return false;
        }
    }

    // 验证4：所有环量自由度必须属于余树边集合
    for (int dof : result.circulation_dofs) {
        if (cotree_edges.find(dof) == cotree_edges.end()) {
            FEEM_ERROR("验证失败: 环量自由度{}不在余树边集合中", dof);
            return false;
        }
    }

    // 验证5：环量自由度数应>= 孔洞数（每个孔洞至少一个自由度）
    if (static_cast<int>(result.circulation_dofs.size()) < result.num_holes) {
        FEEM_WARN("警告: 环量自由度数({})小于孔洞数({}), 可能存在未充分约束的孔洞",
                  result.circulation_dofs.size(), result.num_holes);
        // 这不是致命错误，仅发出警告
    }

    FEEM_DEBUG("结果验证通过: 孔洞数={}, 回路数={}, 环量自由度数={}",
               result.num_holes, result.hole_loop_edges.size(),
               result.circulation_dofs.size());

    return true;
}

} // namespace fe_em
