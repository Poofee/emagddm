/**
 * @file tree_gauge.cpp
 * @brief 电磁物理层 - Tree Gauge核心应用类实现
 * @details 实现TreeGauge类的所有方法，包括5步构建流程、查询接口和解向量操作。
 *
 * @note 实现要点：
 *       - 严格遵循异常安全原则（每个步骤独立try-catch）
 *       - 详细的日志输出（便于调试和性能分析）
 *       - 输入验证（防止无效数据导致未定义行为）
 *       - 内存管理（使用STL容器自动管理资源）
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "tree_gauge.hpp"

namespace fe_em {

// ==================== 构造函数 ====================

TreeGauge::TreeGauge(
    const EMMeshData& mesh_data,
    const std::vector<std::vector<int>>& elem_local_to_global_edge,
    const std::set<int>& constrained_edges
) : mesh_data_(mesh_data),
    elem_local_to_global_edge_(elem_local_to_global_edge),
    constrained_edges_(constrained_edges),
    is_built_(false)
{
    FEEM_DEBUG("TreeGauge构造完成 - 约束边数: {}", constrained_edges_.size());
}

// ==================== 核心构建流程 ====================

void TreeGauge::build() {
    if (is_built_) {
        FEEM_WARN("TreeGauge::build()已被调用过，重复调用将被忽略");
        return;
    }

    FEEM_INFO("========== 开始Tree Gauge构建 ==========");

    try {
        // 步骤1：提取图结构
        extractGraphStructure();
        FEEM_INFO("步骤1[图结构提取]完成 - 节点数: {}, 棱边数: {}", num_nodes_, num_edges_);

        // 步骤2：构建生成树（带边界优先策略）
        buildSpanningTree();
        FEEM_INFO("步骤2[生成树构建]完成 - 树边数: {}, 余树边数: {}",
                  tree_result_.tree_edges.size(), tree_result_.cotree_edges.size());

        // 步骤3：检测多连通域
        detectMultiConnectedDomain();
        FEEM_INFO("步骤3[多连通域检测]完成 - 孔洞数: {}, 环量自由度数: {}",
                  domain_result_.num_holes, domain_result_.circulation_dofs.size());

        // 步骤4：分类边
        classifyEdges();
        FEEM_INFO("步骤4[边分类]完成 - 最终树边: {}, 最终余树边: {}",
                  final_tree_edges_.size(), final_cotree_edges_.size());

        // 步骤5：构建约化映射
        buildReducedMapping();
        FEEM_INFO("步骤5[约化映射]完成 - 约化DOF数: {}", reduced_num_dofs_);

        is_built_ = true;

        // 输出统计摘要
        double compression_rate = (num_edges_ > 0) ?
            (100.0 * (num_edges_ - reduced_num_dofs_) / num_edges_) : 0.0;
        FEEM_INFO("========== Tree Gauge构建完成 ==========");
        FEEM_INFO("压缩率: {:.1f}% (原{}→{})", compression_rate, num_edges_, reduced_num_dofs_);

    } catch (const std::exception& e) {
        FEEM_ERROR("Tree Gauge构建失败: {}", e.what());
        throw std::runtime_error(std::string("TreeGauge::build() 失败: ") + e.what());
    }
}

// ==================== 步骤1：提取图结构 ====================

void TreeGauge::extractGraphStructure() {
    FEEM_INFO("执行步骤1: 提取图结构...");

    try {
        MeshGraphExtractor extractor(mesh_data_, elem_local_to_global_edge_);
        GraphResult result = extractor.extract();

        if (!result.isValid()) {
            throw std::runtime_error("MeshGraphExtractor返回无效结果");
        }

        adjacency_list_ = std::move(result.adjacency_list);
        edge_to_global_id_ = std::move(result.edge_to_global_id);
        num_nodes_ = result.num_nodes;
        num_edges_ = result.num_edges;

        FEEM_DEBUG("图结构提取成功 - 邻接表大小: {}, 棱边映射大小: {}",
                   adjacency_list_.size(), edge_to_global_id_.size());

    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("图结构提取失败: ") + e.what());
    }
}

// ==================== 步骤2：构建生成树 ====================

void TreeGauge::buildSpanningTree() {
    FEEM_INFO("执行步骤2: 构建生成树（带边界优先策略）...");

    try {
        tree_result_ = SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority(
            adjacency_list_,
            edge_to_global_id_,
            constrained_edges_,
            RootNodeStrategy::DIRICHLET_BOUNDARY
        );

        if (!tree_result_.isValid()) {
            throw std::runtime_error("SpanningTreeBuilder返回无效结果");
        }

        FEEM_DEBUG("生成树构建成功 - 根节点: {}, 父节点映射大小: {}",
                   tree_result_.root_node, tree_result_.node_parent.size());

    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("生成树构建失败: ") + e.what());
    }
}

// ==================== 步骤3：检测多连通域 ====================

void TreeGauge::detectMultiConnectedDomain() {
    FEEM_INFO("执行步骤3: 检测多连通域...");

    try {
        domain_result_ = MultiConnectedDomainDetector::detect(
            mesh_data_,
            tree_result_
        );

        if (!domain_result_.isValid()) {
            throw std::runtime_error("MultiConnectedDomainDetector返回无效结果");
        }

        FEEM_DEBUG("多连通域检测成功 - 孔洞数: {}, 回路棱边组数: {}",
                   domain_result_.num_holes, domain_result_.hole_loop_edges.size());

    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("多连通域检测失败: ") + e.what());
    }
}

// ==================== 步骤4：边分类 ====================

void TreeGauge::classifyEdges() {
    FEEM_DEBUG("执行步骤4: 边分类...");

    final_tree_edges_.clear();
    final_cotree_edges_.clear();

    final_tree_edges_ = tree_result_.tree_edges;

    int promoted_count = 0;
    for (int edge_id : tree_result_.cotree_edges) {
        if (constrained_edges_.count(edge_id) > 0) {
            final_tree_edges_.insert(edge_id);
            ++promoted_count;
        }
    }

    if (promoted_count > 0) {
        FEEM_INFO("边分类: {}条被约束的余树边提升为树边", promoted_count);
    }

    for (int edge_id : tree_result_.cotree_edges) {
        if (final_tree_edges_.count(edge_id) == 0 &&
            domain_result_.circulation_dofs.count(edge_id) == 0) {
            final_cotree_edges_.insert(edge_id);
        }
    }

    FEEM_DEBUG("边分类完成 - 树边: {}, 余树边: {}, 环量DOF: {}",
               final_tree_edges_.size(),
               final_cotree_edges_.size(),
               domain_result_.circulation_dofs.size());
}

// ==================== 步骤5：构建约化映射 ====================

void TreeGauge::buildReducedMapping() {
    FEEM_DEBUG("执行步骤5: 构建约化映射...");

    cotree_edge_to_reduced_dof_.clear();
    reduced_num_dofs_ = 0;

    for (int edge_id : final_cotree_edges_) {
        cotree_edge_to_reduced_dof_[edge_id] = reduced_num_dofs_++;
    }

    for (int edge_id : domain_result_.circulation_dofs) {
        if (cotree_edge_to_reduced_dof_.count(edge_id) == 0) {
            cotree_edge_to_reduced_dof_[edge_id] = reduced_num_dofs_++;
        }
    }

    FEEM_DEBUG("约化映射构建完成 - 总约化DOF: {}", reduced_num_dofs_);
}

// ==================== 查询接口实现 ====================

int TreeGauge::getReducedNumDofs() const {
    if (!is_built_) {
        FEEM_WARN("getReducedNumDofs()在build()前被调用，返回0");
        return 0;
    }
    return reduced_num_dofs_;
}

const std::set<int>& TreeGauge::getTreeEdges() const {
    static const std::set<int> empty_set;
    if (!is_built_) {
        FEEM_WARN("getTreeEdges()在build()前被调用，返回空集");
        return empty_set;
    }
    return final_tree_edges_;
}

const std::set<int>& TreeGauge::getCotreeEdges() const {
    static const std::set<int> empty_set;
    if (!is_built_) {
        FEEM_WARN("getCotreeEdges()在build()前被调用，返回空集");
        return empty_set;
    }
    return final_cotree_edges_;
}

const std::set<int>& TreeGauge::getCirculationDOFs() const {
    static const std::set<int> empty_set;
    if (!is_built_) {
        FEEM_WARN("getCirculationDOFs()在build()前被调用，返回空集");
        return empty_set;
    }
    return domain_result_.circulation_dofs;
}

const std::unordered_map<int, int>& TreeGauge::getCotreeEdgeToReducedDOF() const {
    static const std::unordered_map<int, int> empty_map;
    if (!is_built_) {
        FEEM_WARN("getCotreeEdgeToReducedDOF()在build()前被调用，返回空映射");
        return empty_map;
    }
    return cotree_edge_to_reduced_dof_;
}

// ==================== 解向量操作实现 ====================

void TreeGauge::mapReducedSolutionToFull(
    const std::vector<double>& reduced_solution,
    std::vector<double>& full_solution
) const {
    if (!is_built_) {
        throw std::invalid_argument("TreeGauge::mapReducedSolutionToFull(): 必须先调用build()");
    }

    if (static_cast<int>(reduced_solution.size()) != reduced_num_dofs_) {
        throw std::invalid_argument(
            "TreeGauge::mapReducedSolutionToFull(): reduced_solution长度不匹配 "
            "(期望" + std::to_string(reduced_num_dofs_) + ", 实际" +
            std::to_string(reduced_solution.size()) + ")"
        );
    }

    full_solution.assign(static_cast<size_t>(num_edges_), 0.0);

    for (const auto& [edge_id, reduced_idx] : cotree_edge_to_reduced_dof_) {
        if (reduced_idx >= 0 && reduced_idx < static_cast<int>(reduced_solution.size())) {
            if (static_cast<size_t>(edge_id) < full_solution.size()) {
                full_solution[static_cast<size_t>(edge_id)] = reduced_solution[static_cast<size_t>(reduced_idx)];
            } else {
                full_solution.resize(static_cast<size_t>(edge_id) + 1, 0.0);
                full_solution[static_cast<size_t>(edge_id)] = reduced_solution[static_cast<size_t>(reduced_idx)];
            }
        }
    }

    FEEM_DEBUG("解向量映射完成 - reduced_size: {}, full_size: {}",
               reduced_solution.size(), full_solution.size());
}

// ==================== 判断方法实现 ====================

bool TreeGauge::isTreeEdge(int global_edge_id) const {
    if (!is_built_) {
        return false;
    }
    return final_tree_edges_.count(global_edge_id) > 0;
}

bool TreeGauge::isCotreeEdge(int global_edge_id) const {
    if (!is_built_) {
        return false;
    }
    return final_cotree_edges_.count(global_edge_id) > 0;
}

} // namespace fe_em
