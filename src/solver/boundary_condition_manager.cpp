/**
 * @file boundary_condition_manager.cpp
 * @brief 求解器层 - 边界条件管理器实现文件
 * @details 实现BoundaryConditionManager类的全部方法，包括：
 *          - 从EMBoundaryMarker提取Dirichlet约束DOF
 *          - 消去法修正右端项向量（重新计算单元矩阵K_e）
 *          - 回扩完整解向量
 *          - 边界条件冲突检测
 *          - 静电场边界类型映射
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "boundary_condition_manager.hpp"
#include "electrostatic_integrator.hpp"
#include "logger_factory.hpp"

#include <algorithm>
#include <set>
#include <unordered_map>
#include <stdexcept>

namespace solver {

// ==================== 构造函数 ====================

BoundaryConditionManager::BoundaryConditionManager()
    : total_dofs_(0) {
}

// ==================== 公共方法 ====================

bool BoundaryConditionManager::extractDirichletDOFs(
    const std::vector<tool::Boundary>& boundaries,
    const fe_em::EMMeshData& mesh_data,
    const fe_em::EMDOFManager& dof_manager
) {
    clear();

    /* 构建边界名称→tool::Boundary映射，用于查找非DIRICHLET类型的约束值 */
    std::unordered_map<std::string, const tool::Boundary*> name_to_boundary;
    for (const auto& bnd : boundaries) {
        name_to_boundary[bnd.getName()] = &bnd;
    }

    /* 遍历所有边界标记，提取Dirichlet等效约束 */
    for (const auto& marker : mesh_data.boundary_markers) {
        tool::BndType mapped_type = mapToElectrostaticBndType(marker.bnd_type);

        if (mapped_type != tool::BndType::DIRICHLET) {
            continue;
        }

        /* 确定约束值 */
        double dirichlet_value = marker.value;

        /* ODD_SYMMETRY等效为φ=0 */
        if (marker.bnd_type == tool::BndType::ODD_SYMMETRY) {
            dirichlet_value = 0.0;
        }

        /* PERFECT_E：若marker.value为0但Boundary定义了非零电压，使用Boundary的电压值 */
        if (marker.bnd_type == tool::BndType::PERFECT_E && dirichlet_value == 0.0) {
            auto it = name_to_boundary.find(marker.name);
            if (it != name_to_boundary.end()) {
                double voltage = it->second->getVoltage();
                if (voltage != 0.0) {
                    dirichlet_value = voltage;
                }
            }
        }

        /* 从target_ids提取约束节点ID */
        if (!std::holds_alternative<std::vector<int>>(marker.target_ids)) {
            FEEM_WARN("边界标记'{}'的target_ids不是节点ID列表，跳过（BndType={}）",
                      marker.name, tool::bndTypeToString(marker.bnd_type));
            continue;
        }

        const auto& target_node_ids = std::get<std::vector<int>>(marker.target_ids);

        for (int node_id : target_node_ids) {
            DirichletDOF dof;
            dof.global_dof_id = -1;
            dof.value = dirichlet_value;
            dof.boundary_id = marker.id;
            dof.node_id = node_id;
            dirichlet_dofs_.push_back(dof);

            node_dirichlet_values_[node_id] = dirichlet_value;
        }

        FEEM_DEBUG("提取Dirichlet边界: name='{}', type={}, value={}, nodes={}",
                   marker.name, tool::bndTypeToString(marker.bnd_type),
                   dirichlet_value, target_node_ids.size());
    }

    /* 构建节点→自由DOF编号映射 */
    buildNodeToFreeDOFMapping(mesh_data, dof_manager);

    /* 构建自由DOF编号→预编号映射 */
    buildFreeDOFToPrenumMapping(mesh_data);

    /* 设置总DOF数 */
    total_dofs_ = static_cast<int>(dof_manager.getConstrainedDOFValues().size());
    if (total_dofs_ == 0 && !mesh_data.elements.empty()) {
        /* 如果DOF管理器尚未构建，从节点数估算 */
        std::set<int> scalar_nodes;
        for (const auto& elem : mesh_data.elements) {
            if (elem.dof_type == fe_em::DOFType::SCALAR_ONLY ||
                elem.dof_type == fe_em::DOFType::MIXED_AV) {
                for (int nid : elem.node_ids) {
                    scalar_nodes.insert(nid);
                }
            }
        }
        total_dofs_ = static_cast<int>(scalar_nodes.size());
    }

    /* 更新DirichletDOF的global_dof_id（使用预编号） */
    std::set<int> scalar_nodes_sorted;
    for (const auto& elem : mesh_data.elements) {
        if (elem.dof_type == fe_em::DOFType::SCALAR_ONLY ||
            elem.dof_type == fe_em::DOFType::MIXED_AV) {
            for (int nid : elem.node_ids) {
                scalar_nodes_sorted.insert(nid);
            }
        }
    }
    std::map<int, int> node_to_prenum;
    int prenum_counter = 0;
    for (int nid : scalar_nodes_sorted) {
        node_to_prenum[nid] = prenum_counter++;
    }
    for (auto& dof : dirichlet_dofs_) {
        auto it = node_to_prenum.find(dof.node_id);
        if (it != node_to_prenum.end()) {
            dof.global_dof_id = it->second;
        }
    }

    /* 构建constrained_values_映射（global_dof_id → value） */
    for (const auto& dof : dirichlet_dofs_) {
        constrained_values_[dof.global_dof_id] = dof.value;
    }

    /* 冲突检测 */
    if (!detectConflicts()) {
        FEEM_ERROR("Dirichlet边界条件冲突检测失败，同一节点施加了不同约束值");
        return false;
    }

    FEEM_INFO("提取Dirichlet约束完成: 约束DOF数={}, 约束节点数={}",
              dirichlet_dofs_.size(), node_dirichlet_values_.size());

    return true;
}

void BoundaryConditionManager::applyDirichletToRHS(
    Eigen::VectorXd& rhs,
    const fe_em::EMMeshData& mesh_data,
    const std::vector<fe_em::Local2Global>& elem_l2g,
    const std::vector<numeric::MaterialProperties>& materials
) {
    if (node_dirichlet_values_.empty()) {
        FEEM_DEBUG("无Dirichlet约束，跳过RHS修正");
        return;
    }

    FEEM_INFO("开始Dirichlet RHS修正（消去法）...");

    int correction_count = 0;

    /* 遍历所有单元，重新计算K_e并修正RHS */
    for (size_t elem_idx = 0; elem_idx < mesh_data.elements.size(); ++elem_idx) {
        const auto& elem = mesh_data.elements[elem_idx];

        /* 仅处理标量单元（静电场求解器） */
        if (elem.dof_type != fe_em::DOFType::SCALAR_ONLY) {
            continue;
        }

        /* 检查该单元是否包含约束节点（快速跳过无约束单元） */
        bool has_constrained = false;
        for (int node_id : elem.node_ids) {
            if (node_dirichlet_values_.find(node_id) != node_dirichlet_values_.end()) {
                has_constrained = true;
                break;
            }
        }
        if (!has_constrained) {
            continue;
        }

        /* 获取该单元的Local2Global映射 */
        if (elem_idx >= elem_l2g.size()) {
            FEEM_WARN("单元{}的Local2Global映射缺失，跳过RHS修正", elem.id);
            continue;
        }
        const auto& l2g = elem_l2g[elem_idx];

        /* 获取材料属性 */
        if (elem.material_id < 0 ||
            elem.material_id >= static_cast<int>(materials.size())) {
            FEEM_WARN("单元{}的material_id={}越界，跳过RHS修正",
                      elem.id, elem.material_id);
            continue;
        }

        /* 创建积分器并计算单元刚度矩阵K_e */
        auto integrator = std::make_unique<numeric::ElectrostaticIntegrator>(elem.elem_type);
        integrator->setMaterialProperties(materials[elem.material_id]);

        /* 构建单元节点坐标矩阵（3×num_nodes） */
        int num_nodes = static_cast<int>(elem.node_ids.size());
        Eigen::MatrixXd node_coords(3, num_nodes);
        for (int n = 0; n < num_nodes; ++n) {
            int nid = elem.node_ids[n];
            if (nid < 0 || nid >= static_cast<int>(mesh_data.nodes.size())) {
                FEEM_WARN("单元{}的节点ID={}越界，跳过RHS修正", elem.id, nid);
                continue;
            }
            const auto& node = mesh_data.nodes[nid];
            node_coords(0, n) = node.x;
            node_coords(1, n) = node.y;
            node_coords(2, n) = node.z;
        }

        Eigen::MatrixXd K_e = integrator->computeStiffnessMatrix(node_coords);

        /* 对每个(自由DOF行i, 约束DOF列j)对，修正RHS */
        int num_scalar_dofs = l2g.num_scalar_dofs;
        for (int j = 0; j < num_scalar_dofs; ++j) {
            /* 跳过自由DOF列 */
            if (l2g.indices[j] >= 0) {
                continue;
            }

            /* 获取约束节点的Dirichlet值 */
            int constrained_node_id = elem.node_ids[j];
            auto val_it = node_dirichlet_values_.find(constrained_node_id);
            if (val_it == node_dirichlet_values_.end()) {
                continue;
            }
            double phi_c = val_it->second;

            /* 对该约束列的所有自由DOF行进行修正 */
            for (int i = 0; i < num_scalar_dofs; ++i) {
                /* 跳过约束DOF行 */
                if (l2g.indices[i] < 0) {
                    continue;
                }

                int free_dof_index = l2g.indices[i];
                if (free_dof_index >= rhs.size()) {
                    FEEM_WARN("自由DOF索引{}超出RHS向量范围{}，跳过",
                              free_dof_index, rhs.size());
                    continue;
                }

                rhs(free_dof_index) -= K_e(i, j) * phi_c;
                correction_count++;
            }
        }
    }

    FEEM_INFO("Dirichlet RHS修正完成: 修正项数={}", correction_count);
}

Eigen::VectorXd BoundaryConditionManager::expandSolution(
    const Eigen::VectorXd& phi_free,
    const fe_em::EMDOFManager& dof_manager
) const {
    if (total_dofs_ <= 0) {
        FEEM_WARN("总DOF数未初始化，返回空解向量");
        return Eigen::VectorXd();
    }

    Eigen::VectorXd full_solution(total_dofs_);

    /* 初始化：先用DOF管理器的约束值向量填充（约束位置有值，自由位置为0） */
    const auto& constrained_values = dof_manager.getConstrainedDOFValues();
    for (int i = 0; i < total_dofs_ && i < static_cast<int>(constrained_values.size()); ++i) {
        full_solution(i) = constrained_values[i];
    }

    /* 用自由DOF解值填充对应位置 */
    for (const auto& [free_dof_idx, prenum] : free_dof_to_prenum_) {
        if (free_dof_idx >= 0 && free_dof_idx < phi_free.size() &&
            prenum >= 0 && prenum < total_dofs_) {
            full_solution(prenum) = phi_free(free_dof_idx);
        }
    }

    /* 确保约束DOF位置使用node_dirichlet_values_中的值
     * （覆盖constrained_values中可能为0的ODD_SYMMETRY约束） */
    for (const auto& [node_id, value] : node_dirichlet_values_) {
        auto it = node_to_free_dof_.find(node_id);
        if (it != node_to_free_dof_.end() && it->second < 0) {
            /* 该节点是约束节点，需要找到其预编号位置 */
            /* 通过遍历free_dof_to_prenum_的反向映射查找 */
            /* 更直接的方式：重构node_to_prenum映射 */
            /* 由于我们在extractDirichletDOFs中已经设置了global_dof_id，
             * 可以从dirichlet_dofs_中查找 */
            for (const auto& dof : dirichlet_dofs_) {
                if (dof.node_id == node_id && dof.global_dof_id >= 0 &&
                    dof.global_dof_id < total_dofs_) {
                    full_solution(dof.global_dof_id) = value;
                    break;
                }
            }
        }
    }

    FEEM_DEBUG("解向量回扩完成: 自由DOF={}, 总DOF={}",
               phi_free.size(), total_dofs_);

    return full_solution;
}

const std::vector<DirichletDOF>& BoundaryConditionManager::getDirichletDOFs() const {
    return dirichlet_dofs_;
}

const std::map<int, double>& BoundaryConditionManager::getNodeDirichletValues() const {
    return node_dirichlet_values_;
}

bool BoundaryConditionManager::hasDirichletBC() const {
    return !dirichlet_dofs_.empty();
}

void BoundaryConditionManager::clear() {
    dirichlet_dofs_.clear();
    constrained_values_.clear();
    node_dirichlet_values_.clear();
    node_to_free_dof_.clear();
    free_dof_to_prenum_.clear();
    total_dofs_ = 0;
}

// ==================== 私有方法 ====================

bool BoundaryConditionManager::detectConflicts() const {
    /* 检测同一节点是否施加了不同值的Dirichlet条件 */
    std::map<int, double> node_to_value;
    for (const auto& dof : dirichlet_dofs_) {
        auto it = node_to_value.find(dof.node_id);
        if (it != node_to_value.end()) {
            /* 同一节点已有约束，检查值是否一致 */
            double diff = std::abs(it->second - dof.value);
            double max_val = std::max(std::abs(it->second), std::abs(dof.value));
            double tol = 1e-12 + 1e-10 * max_val;
            if (diff > tol) {
                FEEM_ERROR("Dirichlet冲突: 节点{}被施加了不同约束值 {} 和 {}",
                           dof.node_id, it->second, dof.value);
                return false;
            }
        } else {
            node_to_value[dof.node_id] = dof.value;
        }
    }
    return true;
}

tool::BndType BoundaryConditionManager::mapToElectrostaticBndType(tool::BndType bnd_type) {
    switch (bnd_type) {
        /* 等效为Dirichlet的边界类型 */
        case tool::BndType::DIRICHLET:
            return tool::BndType::DIRICHLET;

        case tool::BndType::PERFECT_E:
            /* 理想电边界：切向电场为0，等效于φ=const */
            return tool::BndType::DIRICHLET;

        case tool::BndType::ODD_SYMMETRY:
            /* 奇对称：φ=0 */
            return tool::BndType::DIRICHLET;

        /* 自然边界条件（Neumann等效，无需DOF约束） */
        case tool::BndType::NEUMANN:
            return tool::BndType::NEUMANN;

        case tool::BndType::INSULATION:
            /* 绝缘边界：法向电流为0，等效Neumann ∂φ/∂n=0 */
            return tool::BndType::NEUMANN;

        case tool::BndType::BALLOON:
            /* 气球边界：远场近似，等效Neumann */
            return tool::BndType::NEUMANN;

        case tool::BndType::EVEN_SYMMETRY:
            /* 偶对称：法向导数为0，等效Neumann ∂φ/∂n=0 */
            return tool::BndType::NEUMANN;

        /* 需要特殊处理的边界类型（保留原类型） */
        case tool::BndType::ROBIN:
            return tool::BndType::ROBIN;

        case tool::BndType::PERIODIC:
            return tool::BndType::PERIODIC;

        case tool::BndType::ANTIPERIODIC:
            return tool::BndType::ANTIPERIODIC;

        case tool::BndType::MASTER_SLAVE:
            return tool::BndType::MASTER_SLAVE;

        /* 其他类型默认为自然边界条件 */
        default:
            FEEM_DEBUG("边界类型{}在静电场中默认映射为Neumann（自然边界条件）",
                       tool::bndTypeToString(bnd_type));
            return tool::BndType::NEUMANN;
    }
}

void BoundaryConditionManager::buildNodeToFreeDOFMapping(
    const fe_em::EMMeshData& mesh_data,
    const fe_em::EMDOFManager& dof_manager
) {
    node_to_free_dof_.clear();

    const auto& elem_l2g = dof_manager.getElemLocalToGlobal();

    for (size_t elem_idx = 0; elem_idx < mesh_data.elements.size(); ++elem_idx) {
        const auto& elem = mesh_data.elements[elem_idx];

        /* 仅处理标量单元 */
        if (elem.dof_type != fe_em::DOFType::SCALAR_ONLY &&
            elem.dof_type != fe_em::DOFType::MIXED_AV) {
            continue;
        }

        if (elem_idx >= elem_l2g.size()) {
            continue;
        }

        const auto& l2g = elem_l2g[elem_idx];
        int num_scalar = l2g.num_scalar_dofs;

        for (int i = 0; i < num_scalar; ++i) {
            int node_id = elem.node_ids[i];
            int free_dof = l2g.indices[i];

            /* 如果该节点已有映射，保持一致性 */
            auto it = node_to_free_dof_.find(node_id);
            if (it != node_to_free_dof_.end()) {
                /* 同一节点在不同单元中应映射到相同的自由DOF编号 */
                if (it->second != free_dof && it->second >= 0 && free_dof >= 0) {
                    FEEM_WARN("节点{}的自由DOF映射不一致: {} vs {}",
                              node_id, it->second, free_dof);
                }
                /* 优先保留非-1的映射 */
                if (free_dof >= 0) {
                    node_to_free_dof_[node_id] = free_dof;
                }
            } else {
                node_to_free_dof_[node_id] = free_dof;
            }
        }
    }
}

void BoundaryConditionManager::buildFreeDOFToPrenumMapping(
    const fe_em::EMMeshData& mesh_data
) {
    free_dof_to_prenum_.clear();

    /* 收集所有被标量单元引用的节点ID，按升序排列
     * 这与EMDOFManager::assignPrenumbering()中的std::set遍历顺序一致 */
    std::set<int> scalar_nodes;
    for (const auto& elem : mesh_data.elements) {
        if (elem.dof_type == fe_em::DOFType::SCALAR_ONLY ||
            elem.dof_type == fe_em::DOFType::MIXED_AV) {
            for (int node_id : elem.node_ids) {
                scalar_nodes.insert(node_id);
            }
        }
    }

    /* 按排序顺序分配预编号（与DOF管理器的预编号一致） */
    std::map<int, int> node_to_prenum;
    int prenum = 0;
    for (int node_id : scalar_nodes) {
        node_to_prenum[node_id] = prenum++;
    }

    /* 构建自由DOF编号→预编号映射 */
    for (const auto& [node_id, free_dof] : node_to_free_dof_) {
        if (free_dof >= 0) {
            auto prenum_it = node_to_prenum.find(node_id);
            if (prenum_it != node_to_prenum.end()) {
                free_dof_to_prenum_[free_dof] = prenum_it->second;
            }
        }
    }
}

} // namespace solver
