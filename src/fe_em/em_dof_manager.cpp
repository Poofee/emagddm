/**
 * @file em_dof_manager.cpp
 * @brief 电磁物理层 - 自由度(DOF)管理器核心类实现
 * @details 实现四步DOF编号流程：预编号→标记约束→分块重编号→构建映射表
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include "em_dof_manager.hpp"
#include <algorithm>
#include <tuple>

namespace fe_em {

// ==================== 构造函数与公共接口 ====================

EMDOFManager::EMDOFManager(const EMMeshData& mesh_data,
                           const std::vector<std::vector<int>>& elem_local_to_global_edge,
                           const GlobalEdgeIDGenerator* edge_gen)
    : mesh_data_(mesh_data)
    , elem_local_to_global_edge_(elem_local_to_global_edge)
    , edge_gen_(edge_gen) {
    FEEM_DEBUG("EMDOFManager构造完成，网格节点数: {}, 单元数: {}, edge_gen: {}",
               mesh_data_.getNodeCount(), mesh_data_.getElementCount(),
               edge_gen_ ? "已提供(O(1)查找)" : "未提供(O(n)回退)");
}

void EMDOFManager::build() {
    // 重置所有中间数据和结果数据
    prenum_data_ = PrenumberingData{};
    constraint_data_ = ConstraintData{};
    prenum_to_free_.clear();
    num_free_dofs_ = 0;
    elem_local_to_global_.clear();
    constrained_dof_values_.clear();

    FEEM_INFO("========== 开始DOF编号流程 ==========");

    // 第一步：预编号（为节点和棱边分配临时连续编号）
    assignPrenumbering();
    FEEM_INFO("第一步[预编号]完成 - 标量DOF: {}, 矢量DOF: {}, 总计: {}",
              prenum_data_.num_scalar_prenum, prenum_data_.num_vector_prenum,
              prenum_data_.total_prenum);

    // 第二步：标记约束（识别Dirichlet边界条件对应的受约束DOF）
    markConstrainedDOFs();
    FEEM_INFO("第二步[标记约束]完成 - 约束标量节点: {}, 约束棱边: {}",
              constraint_data_.constrained_scalar_nodes.size(),
              constraint_data_.constrained_edges.size());

    // 第三步：分块重编号（消去法，自由DOF在前从0开始连续编号）
    renumberFreeDOFs();
    FEEM_INFO("第三步[分块重编号]完成 - 自由DOF数: {}", num_free_dofs_);

    // 第四步：构建单元映射表（生成每个单元的Local2Global映射）
    buildMappingTables();
    FEEM_INFO("第四步[构建映射表]完成 - 单元映射表大小: {}", elem_local_to_global_.size());

    FEEM_INFO("========== DOF编号流程结束 ==========");
    FEEM_INFO("统计摘要: 总DOF={}, 自由DOF={}, 约束DOF={}, 单元数={}",
              prenum_data_.total_prenum, num_free_dofs_,
              prenum_data_.total_prenum - num_free_dofs_,
              mesh_data_.getElementCount());
}

int EMDOFManager::getNumFreeDOFs() const {
    return num_free_dofs_;
}

const std::vector<Local2Global>& EMDOFManager::getElemLocalToGlobal() const {
    return elem_local_to_global_;
}

const std::vector<double>& EMDOFManager::getConstrainedDOFValues() const {
    return constrained_dof_values_;
}

// ==================== 静态辅助方法 ====================

bool EMDOFManager::hasScalarDOF(DOFType dof_type) {
    return dof_type == DOFType::SCALAR_ONLY || dof_type == DOFType::MIXED_AV;
}

bool EMDOFManager::hasVectorDOF(DOFType dof_type) {
    return dof_type == DOFType::VECTOR_EDGE_ONLY || dof_type == DOFType::MIXED_AV;
}

// ==================== 第一步：预编号 ====================

void EMDOFManager::assignPrenumbering() {
    int scalar_counter = 0;

    // 收集被标量单元（SCALAR_ONLY或MIXED_AV）引用的节点ID
    std::set<int> scalar_nodes_referenced;
    for (const auto& elem : mesh_data_.elements) {
        if (hasScalarDOF(elem.dof_type)) {
            for (int node_id : elem.node_ids) {
                scalar_nodes_referenced.insert(node_id);
            }
        }
    }

    // 为被引用的标量节点分配预编号（从0开始连续编号）
    for (int node_id : scalar_nodes_referenced) {
        prenum_data_.node_to_prenum[node_id] = scalar_counter++;
    }
    prenum_data_.num_scalar_prenum = scalar_counter;

    // 矢量棱边DOF紧接标量后续编号
    int vector_counter = scalar_counter;
    std::set<int> edges_numbered;  // 已编号的全局棱边ID集合（避免重复）

    // 遍历所有矢量单元（VECTOR_EDGE_ONLY和MIXED_AV），为棱边分配预编号
    for (size_t elem_idx = 0; elem_idx < mesh_data_.elements.size(); ++elem_idx) {
        const auto& elem = mesh_data_.elements[elem_idx];

        if (!hasVectorDOF(elem.dof_type)) {
            continue;
        }

        // 检查全局棱边映射表是否有效
        if (elem_idx >= elem_local_to_global_edge_.size()) {
            throw std::invalid_argument(
                "单元" + std::to_string(elem.id) + "的全局棱边映射缺失");
        }

        const auto& local_to_global_edge = elem_local_to_global_edge_[elem_idx];
        int num_edges = ElementGeometry::get_num_edges(elem.elem_type);

        if (static_cast<int>(local_to_global_edge.size()) != num_edges) {
            throw std::invalid_argument(
                "单元" + std::to_string(elem.id) + "的棱边映射维度不匹配");
        }

        // 按局部棱边顺序遍历，为每条全局棱边分配预编号（如尚未编号）
        for (int local_edge_idx = 0; local_edge_idx < num_edges; ++local_edge_idx) {
            int global_edge_id = local_to_global_edge[local_edge_idx];

            if (edges_numbered.find(global_edge_id) == edges_numbered.end()) {
                prenum_data_.edge_to_prenum[global_edge_id] = vector_counter++;
                edges_numbered.insert(global_edge_id);
            }
        }
    }
    prenum_data_.num_vector_prenum = vector_counter - scalar_counter;
    prenum_data_.total_prenum = vector_counter;
}

// ==================== 第二步：标记约束 ====================

void EMDOFManager::markConstrainedDOFs() {
    // 初始化约束值向量（按预编号顺序，初始值全为0.0）
    constrained_dof_values_.resize(prenum_data_.total_prenum, 0.0);

    // 遍历所有边界标记，识别Dirichlet边界条件
    for (const auto& marker : mesh_data_.boundary_markers) {
        if (marker.bnd_type != BndType::DIRICHLET) {
            continue;  // 仅处理Dirichlet边界条件
        }

        // 判断边界条件的DOF类型并分别处理
        if (hasScalarDOF(marker.dof_type)) {
            // 标量DOF约束：target_ids应持有vector<int>（节点ID列表）
            if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                const auto& node_ids = std::get<std::vector<int>>(marker.target_ids);
                for (int node_id : node_ids) {
                    constraint_data_.constrained_scalar_nodes.insert(node_id);

                    // 记录约束值到对应预编号位置
                    auto it = prenum_data_.node_to_prenum.find(node_id);
                    if (it != prenum_data_.node_to_prenum.end()) {
                        constrained_dof_values_[it->second] = marker.value;
                    }
                }
            }
        }

        if (hasVectorDOF(marker.dof_type)) {
            // 矢量DOF约束：target_ids应持有vector<vector<int>>（棱边的节点对列表）
            if (std::holds_alternative<std::vector<std::vector<int>>>(marker.target_ids)) {
                const auto& edge_node_pairs = std::get<std::vector<std::vector<int>>>(marker.target_ids);

                for (const auto& node_pair : edge_node_pairs) {
                    if (node_pair.size() < 2) {
                        FEEM_WARN("跳过无效的棱边节点对（节点数<2）");
                        continue;
                    }

                    // 将节点对排序，用于匹配全局棱边ID
                    int n1 = std::min(node_pair[0], node_pair[1]);
                    int n2 = std::max(node_pair[0], node_pair[1]);

                    // 优先使用 edge_gen_ 的 O(1) 哈希查找（若已提供）
                    int global_edge_id = -1;
                    if (edge_gen_) {
                        global_edge_id = edge_gen_->getGlobalEdgeID(n1, n2);
                    } else {
                        // 回退方案：遍历 elem_local_to_global_edge_ 进行O(n)匹配
                        bool found = false;
                        for (size_t elem_idx = 0; elem_idx < mesh_data_.elements.size() && !found; ++elem_idx) {
                            const auto& elem = mesh_data_.elements[elem_idx];
                            if (!hasVectorDOF(elem.dof_type)) continue;
                            if (elem_idx >= elem_local_to_global_edge_.size()) continue;

                            const auto& local_edges = ElementGeometry::get_local_edges(elem.elem_type);
                            const auto& local_to_global_edge = elem_local_to_global_edge_[elem_idx];

                            for (size_t local_edge_idx = 0; local_edge_idx < local_edges.size(); ++local_edge_idx) {
                                const auto& [local_n1, local_n2] = local_edges[local_edge_idx];
                                int global_n1 = elem.node_ids[local_n1];
                                int global_n2 = elem.node_ids[local_n2];
                                int gn1 = std::min(global_n1, global_n2);
                                int gn2 = std::max(global_n1, global_n2);
                                if (gn1 == n1 && gn2 == n2) {
                                    global_edge_id = local_to_global_edge[local_edge_idx];
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            FEEM_WARN("未找到棱边({}, {})对应的全局棱边ID（回退模式）", n1, n2);
                        }
                    }

                    if (global_edge_id >= 0) {
                        constraint_data_.constrained_edges.insert(global_edge_id);

                        auto it = prenum_data_.edge_to_prenum.find(global_edge_id);
                        if (it != prenum_data_.edge_to_prenum.end()) {
                            constrained_dof_values_[it->second] = marker.value;
                        }
                    } else {
                        FEEM_WARN("未找到棱边({}, {})对应的全局棱边ID，跳过该约束", n1, n2);
                    }
                }
            }
        }
    }
}

// ==================== 第三步：分块重编号（消去法）====================

void EMDOFManager::renumberFreeDOFs() {
    int free_counter = 0;
    prenum_to_free_.clear();

    // 遍历所有预编号，将未被约束的DOF重新连续编号
    for (int prenum = 0; prenum < prenum_data_.total_prenum; ++prenum) {
        bool is_constrained = false;

        // 检查该预编号是否属于被约束的标量节点
        for (const auto& [node_id, node_prenum] : prenum_data_.node_to_prenum) {
            if (node_prenum == prenum) {
                if (constraint_data_.constrained_scalar_nodes.count(node_id) > 0) {
                    is_constrained = true;
                }
                break;
            }
        }

        // 检查该预编号是否属于被约束的矢量棱边
        if (!is_constrained) {
            for (const auto& [edge_id, edge_prenum] : prenum_data_.edge_to_prenum) {
                if (edge_prenum == prenum) {
                    if (constraint_data_.constrained_edges.count(edge_id) > 0) {
                        is_constrained = true;
                    }
                    break;
                }
            }
        }

        // 自由DOF分配新的连续编号，约束DOF标记为-1
        if (is_constrained) {
            prenum_to_free_[prenum] = -1;
        } else {
            prenum_to_free_[prenum] = free_counter++;
        }
    }

    num_free_dofs_ = free_counter;
}

// ==================== 第四步：构建单元映射表 ====================

void EMDOFManager::buildMappingTables() {
    // 预分配内存（避免频繁扩容）
    elem_local_to_global_.reserve(mesh_data_.elements.size());

    // 遍历所有单元，根据dof_type分别生成映射
    for (const auto& elem : mesh_data_.elements) {
        Local2Global mapping;

        switch (elem.dof_type) {
            case DOFType::SCALAR_ONLY:
                mapping = buildScalarElementMapping(elem);
                break;
            case DOFType::VECTOR_EDGE_ONLY:
                mapping = buildVectorElementMapping(elem.id, elem.elem_type);
                break;
            case DOFType::MIXED_AV:
                mapping = buildMixedAVElementMapping(elem);
                break;
            default:
                throw std::invalid_argument(
                    "不支持的DOF类型: " + dofTypeToString(elem.dof_type));
        }

        mapping.element_id = elem.id;
        mapping.elem_type = elem.elem_type;
        elem_local_to_global_.push_back(std::move(mapping));
    }
}

// ==================== 内部辅助方法：构建单类单元映射 ====================

Local2Global EMDOFManager::buildScalarElementMapping(const Element& elem) const {
    int num_nodes = ElementGeometry::get_num_nodes(elem.elem_type);
    Local2Global mapping(static_cast<size_t>(num_nodes));
    mapping.num_scalar_dofs = num_nodes;
    mapping.num_vector_dofs = 0;

    // 按节点顺序：对每个node_ids[i]，查node_to_prenum → 查prenum_to_free → 填入indices[i]
    for (int i = 0; i < num_nodes; ++i) {
        int node_id = elem.node_ids[i];

        // 查找节点的预编号
        auto prenum_it = prenum_data_.node_to_prenum.find(node_id);
        if (prenum_it == prenum_data_.node_to_prenum.end()) {
            // 该节点未被任何标量单元引用，标记为-1
            mapping.indices[i] = -1;
            continue;
        }

        int prenum = prenum_it->second;

        // 查找预编号对应的自由编号
        auto free_it = prenum_to_free_.find(prenum);
        if (free_it == prenum_to_free_.end()) {
            mapping.indices[i] = -1;
            continue;
        }

        // 填入自由编号（-1表示约束DOF，>=0表示自由DOF）
        mapping.indices[i] = free_it->second;
    }

    return mapping;
}

Local2Global EMDOFManager::buildVectorElementMapping(int elem_id, ElemType elem_type) const {
    int num_edges = ElementGeometry::get_num_edges(elem_type);
    Local2Global mapping(static_cast<size_t>(num_edges));
    mapping.num_scalar_dofs = 0;
    mapping.num_vector_dofs = num_edges;

    // 获取该单元的全局棱边映射
    size_t elem_idx = static_cast<size_t>(elem_id);
    if (elem_idx >= elem_local_to_global_edge_.size()) {
        throw std::invalid_argument(
            "单元" + std::to_string(elem_id) + "的全局棱边映射缺失");
    }

    const auto& local_to_global_edge = elem_local_to_global_edge_[elem_idx];

    // 按局部棱边顺序：查全局棱边ID → 查edge_to_prenum → 查prenum_to_free
    for (int local_edge_idx = 0; local_edge_idx < num_edges; ++local_edge_idx) {
        int global_edge_id = local_to_global_edge[local_edge_idx];

        // 查找棱边的预编号
        auto prenum_it = prenum_data_.edge_to_prenum.find(global_edge_id);
        if (prenum_it == prenum_data_.edge_to_prenum.end()) {
            mapping.indices[local_edge_idx] = -1;
            continue;
        }

        int prenum = prenum_it->second;

        // 查找预编号对应的自由编号
        auto free_it = prenum_to_free_.find(prenum);
        if (free_it == prenum_to_free_.end()) {
            mapping.indices[local_edge_idx] = -1;
            continue;
        }

        // 填入自由编号（-1表示约束DOF，>=0表示自由DOF）
        mapping.indices[local_edge_idx] = free_it->second;
    }

    return mapping;
}

Local2Global EMDOFManager::buildMixedAVElementMapping(const Element& elem) const {
    int num_nodes = ElementGeometry::get_num_nodes(elem.elem_type);
    int num_edges = ElementGeometry::get_num_edges(elem.elem_type);
    int total_dofs = num_nodes + num_edges;

    Local2Global mapping(static_cast<size_t>(total_dofs));
    mapping.num_scalar_dofs = num_nodes;
    mapping.num_vector_dofs = num_edges;

    // ===== 前半部分 [0, num_scalar_dofs)：标量节点DOF（按节点顺序）=====
    for (int i = 0; i < num_nodes; ++i) {
        int node_id = elem.node_ids[i];

        auto prenum_it = prenum_data_.node_to_prenum.find(node_id);
        if (prenum_it == prenum_data_.node_to_prenum.end()) {
            mapping.indices[i] = -1;
            continue;
        }

        int prenum = prenum_it->second;
        auto free_it = prenum_to_free_.find(prenum);
        mapping.indices[i] = (free_it != prenum_to_free_.end()) ? free_it->second : -1;
    }

    // ===== 后半部分 [num_scalar_dofs, end)：矢量棱边DOF（按局部棱边顺序）=====
    size_t elem_idx = static_cast<size_t>(elem.id);
    if (elem_idx >= elem_local_to_global_edge_.size()) {
        throw std::invalid_argument(
            "单元" + std::to_string(elem.id) + "的全局棱边映射缺失");
    }

    const auto& local_to_global_edge = elem_local_to_global_edge_[elem_idx];

    for (int local_edge_idx = 0; local_edge_idx < num_edges; ++local_edge_idx) {
        int global_idx = num_nodes + local_edge_idx;  // 在indices中的位置
        int global_edge_id = local_to_global_edge[local_edge_idx];

        auto prenum_it = prenum_data_.edge_to_prenum.find(global_edge_id);
        if (prenum_it == prenum_data_.edge_to_prenum.end()) {
            mapping.indices[global_idx] = -1;
            continue;
        }

        int prenum = prenum_it->second;
        auto free_it = prenum_to_free_.find(prenum);
        mapping.indices[global_idx] = (free_it != prenum_to_free_.end()) ? free_it->second : -1;
    }

    return mapping;
}

} // namespace fe_em
