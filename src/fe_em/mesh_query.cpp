/**
 * @file mesh_query.cpp
 * @brief 电磁物理层 - 网格几何查询工具类实现
 * @details 实现基于EMMeshData拓扑数据的查询操作。
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#include "mesh_query.hpp"
#include <algorithm>
#include <iostream>

namespace fe_em {

// ==================== 公共接口实现 ====================

std::vector<Node> MeshQuery::get_element_nodes(
    const Element& element, 
    const std::vector<Node>& all_nodes
) {
    std::vector<Node> elem_nodes;
    elem_nodes.reserve(element.node_ids.size());
    
    // 统一使用id字段进行精确匹配（避免下标访问与ID字段混淆）
    for (int node_id : element.node_ids) {
        auto it = std::find_if(all_nodes.begin(), all_nodes.end(),
                              [node_id](const Node& n) { return n.id == node_id; });
        
        if (it != all_nodes.end()) {
            elem_nodes.push_back(*it);
        } else {
            // 节点未找到时记录警告（使用cerr避免依赖日志系统）
            std::cerr << "[MeshQuery] Warning: Node id=" << node_id 
                      << " not found in mesh\n";
        }
    }
    
    return elem_nodes;
}

bool MeshQuery::is_node_in_region(
    int node_id, 
    int region_id, 
    const std::vector<Node>& all_nodes
) {
    // 通过id字段查找节点并检查region_id
    for (const auto& node : all_nodes) {
        if (node.id == node_id) {
            return node.region_id == region_id;
        }
    }
    return false;
}

bool MeshQuery::is_element_in_region(
    const Element& element, 
    int region_id
) {
    return element.region_id == region_id;
}

std::vector<int> MeshQuery::get_elements_in_region(
    int region_id, 
    const std::vector<Element>& all_elements
) {
    std::vector<int> result;
    
    for (const auto& elem : all_elements) {
        if (elem.region_id == region_id) {
            result.push_back(elem.id);
        }
    }
    
    return result;
}

std::vector<int> MeshQuery::get_nodes_in_region(
    int region_id, 
    const std::vector<Node>& all_nodes
) {
    std::vector<int> result;
    
    for (const auto& node : all_nodes) {
        if (node.region_id == region_id) {
            result.push_back(node.id);
        }
    }
    
    return result;
}

} // namespace fe_em
