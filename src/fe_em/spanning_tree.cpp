/**
 * @file spanning_tree.cpp
 * @brief 电磁物理层 - 生成树构建器实现文件
 * @details 实现BFS/DFS生成树构建算法的核心逻辑，
 *          包括边界优先策略和多连通分量处理。
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "spanning_tree.hpp"
#include <algorithm>
#include <stdexcept>
#include <deque>
#include <vector>

namespace fe_em {

// ==================== 公共接口：buildBFSTree() ====================

SpanningTreeResult SpanningTreeBuilder::buildBFSTree(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
    int root_node
) {
    FEEM_INFO("========== SpanningTreeBuilder: 开始构建BFS生成树 ==========");

    // ---------- 输入验证 ----------
    int num_nodes = static_cast<int>(adjacency_list.size());

    if (num_nodes == 0) {
        FEEM_WARN("邻接表为空（0个节点），返回空结果");
        return SpanningTreeResult{};
    }

    if (edge_to_global_id.empty()) {
        FEEM_WARN("棱边映射表为空，返回空结果");
        return SpanningTreeResult{};
    }

    // 验证根节点范围
    if (root_node < 0 || root_node >= num_nodes) {
        FEEM_WARN("根节点{}超出范围[0, {})，调整为0", root_node, num_nodes);
        root_node = 0;
    }

    FEEM_DEBUG("输入统计: 节点数={}, 棱边数={}, 根节点={}",
              num_nodes, edge_to_global_id.size(), root_node);

    // ========== 初始化数据结构 ==========
    SpanningTreeResult result;
    result.root_node = root_node;

    std::vector<bool> visited(num_nodes, false);
    std::deque<int> bfs_queue;

    // ---------- 从指定根节点开始BFS ----------
    visited[root_node] = true;
    result.node_parent[root_node] = -1;  // 根节点的父节点为-1
    bfs_queue.push_back(root_node);

    int tree_edge_count = 0;

    FEEM_DEBUG("开始BFS遍历（起始节点={}）", root_node);

    while (!bfs_queue.empty()) {
        int u = bfs_queue.front();
        bfs_queue.pop_front();

        FEEM_DEBUG("处理节点{}（邻居数={})", u, adjacency_list[u].size());

        // 遍历当前节点的所有邻居
        for (int v : adjacency_list[u]) {
            if (!visited[v]) {
                // 标记为已访问
                visited[v] = true;

                // 查找边(u,v)对应的全局棱边ID
                int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);

                if (global_edge_id >= 0) {
                    // 将该边加入树边集合
                    result.tree_edges.insert(global_edge_id);
                    tree_edge_count++;

                    FEEM_DEBUG("  添加树边: 节点{}→节点{} (全局棱边ID={})",
                              u, v, global_edge_id);
                } else {
                    FEEM_WARN("  未找到边({},{})的全局棱边ID映射", u, v);
                }

                // 设置父节点关系
                result.node_parent[v] = u;

                // 将v入队
                bfs_queue.push_back(v);
            }
        }
    }

    // ========== 处理其他连通分量 ==========
    for (int i = 0; i < num_nodes; ++i) {
        if (!visited[i]) {
            FEEM_INFO("发现未访问的连通分量，从节点{}开始新的BFS", i);

            visited[i] = true;
            result.node_parent[i] = -1;  // 新连通分量的根节点
            bfs_queue.push_back(i);

            while (!bfs_queue.empty()) {
                int u = bfs_queue.front();
                bfs_queue.pop_front();

                for (int v : adjacency_list[u]) {
                    if (!visited[v]) {
                        visited[v] = true;

                        int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);
                        if (global_edge_id >= 0) {
                            result.tree_edges.insert(global_edge_id);
                            tree_edge_count++;
                        }

                        result.node_parent[v] = u;
                        bfs_queue.push_back(v);
                    }
                }
            }
        }
    }

    // ========== 构建余树边集合 ==========
    result.cotree_edges = buildCotreeEdges(result.tree_edges, edge_to_global_id);

    // ========== 日志输出 ==========
    int connected_components = 0;
    for (const auto& [node, parent] : result.node_parent) {
        if (parent == -1) {
            connected_components++;
        }
    }

    FEEM_INFO("========== SpanningTreeBuilder: BFS生成树构建完成 ==========");
    FEEM_INFO("输入统计: 节点数={}, 棱边数={}, 连通分量数={}",
              num_nodes, edge_to_global_id.size(), connected_components);
    FEEM_INFO("根节点选择: 节点{}(策略=FIRST_NODE)", root_node);
    FEEM_INFO("生成树: 树边数={}, 余树边数={}",
              result.tree_edges.size(), result.cotree_edges.size());
    FEEM_INFO("理论验证: 树边数应为(V-C)={}, 余树边数应为(E-V+C)={}",
              num_nodes - connected_components,
              edge_to_global_id.size() - num_nodes + connected_components);
    FEEM_INFO("======================================================");

    return result;
}


// ==================== 公共接口：buildDFSTree() ====================

SpanningTreeResult SpanningTreeBuilder::buildDFSTree(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
    int root_node
) {
    FEEM_INFO("========== SpanningTreeBuilder: 开始构建DFS生成树 ==========");

    // ---------- 输入验证 ----------
    int num_nodes = static_cast<int>(adjacency_list.size());

    if (num_nodes == 0) {
        FEEM_WARN("邻接表为空（0个节点），返回空结果");
        return SpanningTreeResult{};
    }

    if (edge_to_global_id.empty()) {
        FEEM_WARN("棱边映射表为空，返回空结果");
        return SpanningTreeResult{};
    }

    if (root_node < 0 || root_node >= num_nodes) {
        FEEM_WARN("根节点{}超出范围[0, {})，调整为0", root_node, num_nodes);
        root_node = 0;
    }

    FEEM_DEBUG("输入统计: 节点数={}, 棱边数={}, 根节点={}",
              num_nodes, edge_to_global_id.size(), root_node);

    // ========== 初始化数据结构 ==========
    SpanningTreeResult result;
    result.root_node = root_node;

    std::vector<bool> visited(num_nodes, false);
    std::vector<int> dfs_stack;

    // ---------- 从指定根节点开始DFS ----------
    visited[root_node] = true;
    result.node_parent[root_node] = -1;
    dfs_stack.push_back(root_node);

    int tree_edge_count = 0;

    FEEM_DEBUG("开始DFS遍历（起始节点={}）", root_node);

    while (!dfs_stack.empty()) {
        int u = dfs_stack.back();
        dfs_stack.pop_back();

        FEEM_DEBUG("处理节点{}（邻居数={})", u, adjacency_list[u].size());

        // 遍历当前节点的所有邻居（使用显式栈的非递归DFS）
        for (int v : adjacency_list[u]) {
            if (!visited[v]) {
                // 标记为已访问
                visited[v] = true;

                // 查找边(u,v)对应的全局棱边ID
                int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);

                if (global_edge_id >= 0) {
                    result.tree_edges.insert(global_edge_id);
                    tree_edge_count++;

                    FEEM_DEBUG("  添加树边: 节点{}→节点{} (全局棱边ID={})",
                              u, v, global_edge_id);
                } else {
                    FEEM_WARN("  未找到边({},{})的全局棱边ID映射", u, v);
                }

                // 设置父节点关系
                result.node_parent[v] = u;

                // 将v压栈（深度优先）
                dfs_stack.push_back(v);
            }
        }
    }

    // ========== 处理其他连通分量 ==========
    for (int i = 0; i < num_nodes; ++i) {
        if (!visited[i]) {
            FEEM_INFO("发现未访问的连通分量，从节点{}开始新的DFS", i);

            visited[i] = true;
            result.node_parent[i] = -1;
            dfs_stack.push_back(i);

            while (!dfs_stack.empty()) {
                int u = dfs_stack.back();
                dfs_stack.pop_back();

                for (int v : adjacency_list[u]) {
                    if (!visited[v]) {
                        visited[v] = true;

                        int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);
                        if (global_edge_id >= 0) {
                            result.tree_edges.insert(global_edge_id);
                            tree_edge_count++;
                        }

                        result.node_parent[v] = u;
                        dfs_stack.push_back(v);
                    }
                }
            }
        }
    }

    // ========== 构建余树边集合 ==========
    result.cotree_edges = buildCotreeEdges(result.tree_edges, edge_to_global_id);

    // ========== 日志输出 ==========
    int connected_components = 0;
    for (const auto& [node, parent] : result.node_parent) {
        if (parent == -1) {
            connected_components++;
        }
    }

    FEEM_INFO("========== SpanningTreeBuilder: DFS生成树构建完成 ==========");
    FEEM_INFO("输入统计: 节点数={}, 棱边数={}, 连通分量数={}",
              num_nodes, edge_to_global_id.size(), connected_components);
    FEEM_INFO("根节点选择: 节点{}(策略=FIRST_NODE)", root_node);
    FEEM_INFO("生成树: 树边数={}, 余树边数={}",
              result.tree_edges.size(), result.cotree_edges.size());
    FEEM_INFO("======================================================");

    return result;
}


// ==================== 公共接口：buildBFSTreeWithBoundaryPriority() ====================

SpanningTreeResult SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
    const std::set<int>& dirichlet_edges,
    RootNodeStrategy strategy,
    int user_root_node
) {
    FEEM_INFO("========== SpanningTreeBuilder: 开始构建带边界优先的BFS生成树 ==========");

    // ---------- 输入验证 ----------
    int num_nodes = static_cast<int>(adjacency_list.size());

    if (num_nodes == 0) {
        FEEM_WARN("邻接表为空（0个节点），返回空结果");
        return SpanningTreeResult{};
    }

    if (edge_to_global_id.empty()) {
        FEEM_WARN("棱边映射表为空，返回空结果");
        return SpanningTreeResult{};
    }

    FEEM_DEBUG("输入统计: 节点数={}, 棱边数={}, Dirichlet边数={}",
              num_nodes, edge_to_global_id.size(), dirichlet_edges.size());

    // ========== 智能根节点选择 ==========
    int root_node = selectRootNode(adjacency_list, edge_to_global_id,
                                   dirichlet_edges, strategy, user_root_node, num_nodes);

    FEEM_INFO("根节点选择: 节点{}(策略={})", root_node, static_cast<int>(strategy));

    // ========== 初始化数据结构 ==========
    SpanningTreeResult result;
    result.root_node = root_node;

    std::vector<bool> visited(num_nodes, false);
    std::deque<int> bfs_queue;

    // ---------- 构建Dirichlet边的反向映射（全局ID → 节点对） ----------
    std::unordered_map<int, std::pair<int, int>> dirichlet_edge_to_nodes;
    for (const auto& [edge_key, global_id] : edge_to_global_id) {
        if (dirichlet_edges.count(global_id) > 0) {
            dirichlet_edge_to_nodes[global_id] = edge_key;
        }
    }

    FEEM_DEBUG("Dirichlet边反向映射构建完成（共{}条）", dirichlet_edge_to_nodes.size());

    // ---------- 从根节点开始带优先级的BFS ----------
    visited[root_node] = true;
    result.node_parent[root_node] = -1;
    bfs_queue.push_back(root_node);

    int tree_edge_count = 0;
    int dirichlet_tree_edge_count = 0;

    FEEM_DEBUG("开始带边界优先的BFS遍历（起始节点={}）", root_node);

    while (!bfs_queue.empty()) {
        int u = bfs_queue.front();
        bfs_queue.pop_front();

        FEEM_DEBUG("处理节点{}（邻居数={})", u, adjacency_list[u].size());

        // 分离Dirichlet边邻居和普通邻居
        std::vector<std::pair<int, int>> dirichlet_neighbors;  // (neighbor, global_edge_id)
        std::vector<int> normal_neighbors;

        for (int v : adjacency_list[u]) {
            if (visited[v]) {
                continue;  // 已访问的跳过
            }

            // 查找边(u,v)的全局棱边ID
            int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);

            if (global_edge_id < 0) {
                FEEM_WARN("  未找到边({},{})的全局棱边ID映射", u, v);
                continue;
            }

            // 判断是否为Dirichlet边
            if (dirichlet_edges.count(global_edge_id) > 0) {
                dirichlet_neighbors.emplace_back(v, global_edge_id);
                FEEM_DEBUG("  发现Dirichlet边: 节点{}→节点{} (全局棱边ID={})",
                          u, v, global_edge_id);
            } else {
                normal_neighbors.push_back(v);
            }
        }

        // 优先处理Dirichlet边邻居（确保它们先加入队列）
        for (const auto& [v, global_edge_id] : dirichlet_neighbors) {
            if (!visited[v]) {  // 再次检查（防止重复）
                visited[v] = true;
                result.tree_edges.insert(global_edge_id);
                tree_edge_count++;
                dirichlet_tree_edge_count++;
                result.node_parent[v] = u;
                bfs_queue.push_back(v);

                FEEM_DEBUG("  [优先] 添加Dirichlet树边: 节点{}→节点{} (全局棱边ID={})",
                          u, v, global_edge_id);
            }
        }

        // 然后处理普通邻居
        for (int v : normal_neighbors) {
            if (!visited[v]) {  // 再次检查（防止重复）
                visited[v] = true;

                int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);
                if (global_edge_id >= 0) {
                    result.tree_edges.insert(global_edge_id);
                    tree_edge_count++;
                }

                result.node_parent[v] = u;
                bfs_queue.push_back(v);
            }
        }
    }

    // ========== 处理其他连通分量 ==========
    for (int i = 0; i < num_nodes; ++i) {
        if (!visited[i]) {
            FEEM_INFO("发现未访问的连通分量，从节点{}开始新的BFS", i);

            visited[i] = true;
            result.node_parent[i] = -1;
            bfs_queue.push_back(i);

            while (!bfs_queue.empty()) {
                int u = bfs_queue.front();
                bfs_queue.pop_front();

                for (int v : adjacency_list[u]) {
                    if (!visited[v]) {
                        visited[v] = true;

                        int global_edge_id = findGlobalEdgeID(u, v, edge_to_global_id);
                        if (global_edge_id >= 0) {
                            result.tree_edges.insert(global_edge_id);
                            tree_edge_count++;
                        }

                        result.node_parent[v] = u;
                        bfs_queue.push_back(v);
                    }
                }
            }
        }
    }

    // ========== 构建余树边集合 ==========
    result.cotree_edges = buildCotreeEdges(result.tree_edges, edge_to_global_id);

    // ========== 验证Dirichlet边是否都在树中 ==========
    int dirichlet_in_tree = 0;
    int dirichlet_in_cotree = 0;
    for (int dirichlet_edge_id : dirichlet_edges) {
        if (result.tree_edges.count(dirichlet_edge_id) > 0) {
            dirichlet_in_tree++;
        } else if (result.cotree_edges.count(dirichlet_edge_id) > 0) {
            dirichlet_in_cotree++;
            FEEM_WARN("Dirichlet边(全局ID={})在余树中（可能影响边界条件一致性）",
                     dirichlet_edge_id);
        }
    }

    // ========== 日志输出 ==========
    int connected_components = 0;
    for (const auto& [node, parent] : result.node_parent) {
        if (parent == -1) {
            connected_components++;
        }
    }

    FEEM_INFO("========== SpanningTreeBuilder: 带边界优先BFS生成树构建完成 ==========");
    FEEM_INFO("输入统计: 节点数={}, 棱边数={}, 连通分量数={}",
              num_nodes, edge_to_global_id.size(), connected_components);
    FEEM_INFO("根节点选择: 节点{}(策略={})", root_node, static_cast<int>(strategy));
    FEEM_INFO("生成树: 树边数={}, 余树边数={}", result.tree_edges.size(), result.cotree_edges.size());
    FEEM_INFO("Dirichlet边: 总数={}, 在树中={}, 在余树中={}",
              dirichlet_edges.size(), dirichlet_in_tree, dirichlet_in_cotree);
    FEEM_INFO("======================================================");

    return result;
}


// ==================== 私有辅助方法：selectRootNode() ====================

int SpanningTreeBuilder::selectRootNode(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id,
    const std::set<int>& dirichlet_edges,
    RootNodeStrategy strategy,
    int user_root_node,
    int num_nodes
) {
    switch (strategy) {
        case RootNodeStrategy::FIRST_NODE:
            FEEM_DEBUG("根节点策略: FIRST_NODE → 选择节点0");
            return 0;

        case RootNodeStrategy::MAX_DEGREE: {
            int max_degree = -1;
            int best_node = 0;

            for (int i = 0; i < num_nodes; ++i) {
                int degree = static_cast<int>(adjacency_list[i].size());
                if (degree > max_degree) {
                    max_degree = degree;
                    best_node = i;
                }
            }

            FEEM_DEBUG("根节点策略: MAX_DEGREE → 选择节点{}(度数={})", best_node, max_degree);
            return best_node;
        }

        case RootNodeStrategy::DIRICHLET_BOUNDARY: {
            if (!dirichlet_edges.empty()) {
                // 从dirichlet_edges中选择第一条边的一个端点
                int first_dirichlet_edge = *dirichlet_edges.begin();

                // 反向查找该全局棱边ID对应的节点对
                for (const auto& [edge_key, global_id] : edge_to_global_id) {
                    if (global_id == first_dirichlet_edge) {
                        int selected_node = edge_key.first;  // 选择较小的节点索引
                        FEEM_DEBUG("根节点策略: DIRICHLET_BOUNDARY → 选择节点{}(来自Dirichlet边全局ID={})",
                                  selected_node, first_dirichlet_edge);
                        return selected_node;
                    }
                }

                FEEM_WARN("未找到Dirichlet边全局ID={}对应的节点对，退化为FIRST_NODE",
                         first_dirichlet_edge);
                return 0;
            } else {
                FEEM_WARN("Dirichlet边集合为空，DIRICHLET_BOUNDARY策略退化为FIRST_NODE");
                return 0;
            }
        }

        case RootNodeStrategy::USER_SPECIFIED:
            if (user_root_node >= 0 && user_root_node < num_nodes) {
                FEEM_DEBUG("根节点策略: USER_SPECIFIED → 选择节点{}", user_root_node);
                return user_root_node;
            } else {
                FEEM_WARN("用户指定的根节点{}无效（应在[0, {})范围内），退化为FIRST_NODE",
                         user_root_node, num_nodes);
                return 0;
            }

        default:
            FEEM_WARN("未知的根节点策略，退化为FIRST_NODE");
            return 0;
    }
}


// ==================== 私有辅助方法：findGlobalEdgeID() ====================

int SpanningTreeBuilder::findGlobalEdgeID(
    int node1,
    int node2,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id
) {
    // 构建排序后的节点对（无向图的标准表示）
    std::pair<int, int> edge_key = (node1 < node2)
                                   ? std::make_pair(node1, node2)
                                   : std::make_pair(node2, node1);

    // 在哈希表中查找
    auto it = edge_to_global_id.find(edge_key);
    if (it != edge_to_global_id.end()) {
        return it->second;  // 返回全局棱边ID
    }

    return -1;  // 未找到
}


// ==================== 私有辅助方法：buildCotreeEdges() ====================

std::set<int> SpanningTreeBuilder::buildCotreeEdges(
    const std::set<int>& tree_edges,
    const std::unordered_map<std::pair<int, int>, int, EdgeKeyHash>& edge_to_global_id
) {
    std::set<int> cotree_edges;

    // 遍历所有棱边，将不在tree_edges中的归入cotree_edges
    for (const auto& [edge_key, global_id] : edge_to_global_id) {
        if (tree_edges.count(global_id) == 0) {
            cotree_edges.insert(global_id);
        }
    }

    return cotree_edges;
}

} // namespace fe_em
