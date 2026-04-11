/**
 * @file electrostatic_force_calculator.cpp
 * @brief 求解器层 - 静电力计算器实现
 * @details 采用虚功法（Virtual Work Method）计算导体上的静电力和力矩。
 *          核心算法：对导体节点施加微小位移δ，通过能量差分计算力：
 *          F_i = -(W(s+δ) - W(s-δ)) / (2δ)
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "electrostatic_force_calculator.hpp"
#include "electrostatic_solver.hpp"
#include "field_data.hpp"
#include "logger_factory.hpp"
#include "math_constants.hpp"

#include <cmath>
#include <algorithm>
#include <unordered_set>

namespace solver {

// ============================================================
// 构造/析构
// ============================================================

ElectrostaticForceCalculator::ElectrostaticForceCalculator() {
    FEEM_DEBUG("ElectrostaticForceCalculator 构造完成");
}

// ============================================================
// 公共接口：力计算
// ============================================================

Eigen::Vector3d ElectrostaticForceCalculator::computeForce(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    const std::string& conductor_name,
    tool::DimType dim_type,
    double delta
) {
    FEEM_INFO("========================================");
    FEEM_INFO("开始计算导体 [{}] 上的静电力", conductor_name);
    FEEM_INFO("  虚位移量: {} m", delta);
    FEEM_INFO("  维度类型: {}", static_cast<int>(dim_type));

    // 步骤1：获取导体节点ID列表
    auto conductor_node_ids = getConductorNodeIds(boundaries, mesh_data, conductor_name);
    if (conductor_node_ids.empty()) {
        FEEM_ERROR("未找到导体 [{}] 的节点信息，无法计算力", conductor_name);
        return Eigen::Vector3d::Zero();
    }
    FEEM_INFO("  导体包含 {} 个节点", conductor_node_ids.size());

    // 步骤2：构建节点ID快速查找集合
    std::unordered_set<int> node_set(conductor_node_ids.begin(), conductor_node_ids.end());

    // 步骤3：对每个坐标方向进行中心差分求力
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> directions = {
        Eigen::Vector3d(1.0, 0.0, 0.0),   // X方向
        Eigen::Vector3d(0.0, 1.0, 0.0),   // Y方向
        Eigen::Vector3d(0.0, 0.0, 1.0)    // Z方向
    };

    // 2D问题只计算X和Y方向的力
    int num_directions = (dim_type == tool::DimType::D2) ? 2 : 3;

    for (int i = 0; i < num_directions; ++i) {
        const auto& dir = directions[i];

        // 正向位移
        auto mesh_plus = displaceConductorNodes(mesh_data, conductor_node_ids, dir * delta);
        double W_plus = solveForEnergy(mesh_plus, materials, boundaries, dim_type);

        if (std::isnan(W_plus)) {
            FEEM_ERROR("正向位移求解失败（方向{}），跳过该分量", i);
            continue;
        }

        // 负向位移
        auto mesh_minus = displaceConductorNodes(mesh_data, conductor_node_ids, -dir * delta);
        double W_minus = solveForEnergy(mesh_minus, materials, boundaries, dim_type);

        if (std::isnan(W_minus)) {
            FEEM_ERROR("负向位移求解失败（方向{}），跳过该分量", i);
            continue;
        }

        // 中心差分：F_i = -(W_+ - W_-) / (2*delta)
        force(i) = -(W_plus - W_minus) / (2.0 * delta);

        FEEM_DEBUG("  方向 {}: W_+={:.6e} J, W_-={:.6e} J, F_{}={:.6e} N",
                   i, W_plus, W_minus, i, force(i));
    }

    FEEM_INFO("  静电力结果: Fx={:.6e} N, Fy={:.6e} N, Fz={:.6e} N",
              force.x(), force.y(), force.z());
    FEEM_INFO("========================================");

    return force;
}

// ============================================================
// 公共接口：力矩计算
// ============================================================

double ElectrostaticForceCalculator::computeTorque(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    const std::string& conductor_name,
    const Eigen::Vector3d& axis_point,
    const Eigen::Vector3d& axis_direction,
    tool::DimType dim_type,
    double delta
) {
    FEEM_INFO("========================================");
    FEEM_INFO("开始计算导体 [{}] 上的静电力矩", conductor_name);
    FEEM_INFO("  旋转轴点: ({}, {}, {})", axis_point.x(), axis_point.y(), axis_point.z());
    FEEM_INFO("  旋转轴方向: ({}, {}, {})",
              axis_direction.x(), axis_direction.y(), axis_direction.z());
    FEEM_INFO("  虚角位移: {} rad", delta);

    // 步骤1：获取导体节点ID列表
    auto conductor_node_ids = getConductorNodeIds(boundaries, mesh_data, conductor_name);
    if (conductor_node_ids.empty()) {
        FEEM_ERROR("未找到导体 [{}] 的节点信息，无法计算力矩", conductor_name);
        return std::nan("");
    }
    FEEM_INFO("  导体包含 {} 个节点", conductor_node_ids.size());

    // 步骤2：归一化旋转轴方向向量
    Eigen::Vector3d k = axis_direction.normalized();
    if (k.norm() < 1e-12) {
        FEEM_ERROR("旋转轴方向向量为零，无法归一化");
        return std::nan("");
    }

    // 步骤3：正向旋转 +δ
    auto mesh_plus = rotateConductorNodes(
        mesh_data, conductor_node_ids, axis_point, k, delta
    );
    double W_plus = solveForEnergy(mesh_plus, materials, boundaries, dim_type);

    if (std::isnan(W_plus)) {
        FEEM_ERROR("正向旋转求解失败，无法计算力矩");
        return std::nan("");
    }

    // 步骤4：负向旋转 -δ
    auto mesh_minus = rotateConductorNodes(
        mesh_data, conductor_node_ids, axis_point, k, -delta
    );
    double W_minus = solveForEnergy(mesh_minus, materials, boundaries, dim_type);

    if (std::isnan(W_minus)) {
        FEEM_ERROR("负向旋转求解失败，无法计算力矩");
        return std::nan("");
    }

    // 步骤5：中心差分计算力矩 T = -(W_+ - W_-) / (2*delta)
    double torque = -(W_plus - W_minus) / (2.0 * delta);

    FEEM_INFO("  力矩结果: T={:.6e} N·m", torque);
    FEEM_INFO("========================================");

    return torque;
}

// ============================================================
// 公共接口：清空
// ============================================================

void ElectrostaticForceCalculator::clear() {
    FEEM_DEBUG("ElectrostaticForceCalculator 数据已清空");
}

// ============================================================
// 私有方法：获取导体节点ID列表
// ============================================================

std::vector<int> ElectrostaticForceCalculator::getConductorNodeIds(
    const std::vector<tool::Boundary>& boundaries,
    const fe_em::EMMeshData& mesh_data,
    const std::string& conductor_name
) {
    // 在边界条件列表中查找名称匹配的边界
    for (const auto& bc : boundaries) {
        if (bc.getName() == conductor_name) {
            // 找到匹配的边界条件，从其关联的几何对象中提取节点
            // 方法1：通过Boundary的objects查找对应的EMBoundaryMarker
            for (const auto& marker : mesh_data.boundary_markers) {
                if (marker.name == conductor_name || marker.name == bc.getName()) {
                    // 从target_ids中提取节点ID
                    if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                        return std::get<std::vector<int>>(marker.target_ids);
                    }
                }
            }

            // 方法2：如果Boundary有关联的对象名，尝试在boundary_markers中查找
            const auto& objects = bc.getObjects();
            for (const auto& obj_name : objects) {
                for (const auto& marker : mesh_data.boundary_markers) {
                    if (marker.name == obj_name) {
                        if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                            return std::get<std::vector<int>>(marker.target_ids);
                        }
                    }
                }
            }
        }
    }

    // 如果直接名称匹配未找到，尝试模糊匹配
    for (const auto& bc : boundaries) {
        const auto& objects = bc.getObjects();
        for (const auto& obj : objects) {
            if (obj.find(conductor_name) != std::string::npos ||
                conductor_name.find(obj) != std::string::npos) {
                for (const auto& marker : mesh_data.boundary_markers) {
                    if (marker.name == obj || marker.name == bc.getName()) {
                        if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                            FEEM_DEBUG("通过模糊匹配找到导体 [{}] -> 边界标记 [{}]",
                                       conductor_name, marker.name);
                            return std::get<std::vector<int>>(marker.target_ids);
                        }
                    }
                }
            }
        }
    }

    FEEM_WARN("未找到导体 [{}] 对应的节点信息", conductor_name);
    return {};
}

// ============================================================
// 私有方法：对导体节点施加位移
// ============================================================

fe_em::EMMeshData ElectrostaticForceCalculator::displaceConductorNodes(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<int>& conductor_node_ids,
    const Eigen::Vector3d& displacement
) {
    // 深拷贝网格数据
    fe_em::EMMeshData modified_mesh = mesh_data;

    // 构建节点ID快速查找集合
    std::unordered_set<int> node_set(conductor_node_ids.begin(),
                                      conductor_node_ids.end());

    // 对每个属于导体的节点施加位移
    for (auto& node : modified_mesh.nodes) {
        if (node_set.count(node.id) > 0) {
            node.x += displacement.x();
            node.y += displacement.y();
            node.z += displacement.z();
        }
    }

    FEEM_DEBUG("已对 {} 个导体节点施加位移 ({}, {}, {})",
               conductor_node_ids.size(),
               displacement.x(), displacement.y(), displacement.z());

    return modified_mesh;
}

// ============================================================
// 私有方法：对导体节点施加旋转变换（Rodrigues公式）
// ============================================================

fe_em::EMMeshData ElectrostaticForceCalculator::rotateConductorNodes(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<int>& conductor_node_ids,
    const Eigen::Vector3d& axis_point,
    const Eigen::Vector3d& axis_unit_dir,
    double angle
) {
    // 深拷贝网格数据
    fe_em::EMMeshData modified_mesh = mesh_data;

    // 构建节点ID快速查找集合
    std::unordered_set<int> node_set(conductor_node_ids.begin(),
                                      conductor_node_ids.end());

    // Rodrigues旋转公式的预计算量
    double cos_a = std::cos(angle);
    double sin_a = std::sin(angle);
    double one_minus_cos_a = 1.0 - cos_a;

    // 对每个属于导体的节点应用旋转变换
    for (auto& node : modified_mesh.nodes) {
        if (node_set.count(node.id) > 0) {
            // 将节点平移到以axis_point为原点的坐标系
            Eigen::Vector3d p(node.x - axis_point.x(),
                              node.y - axis_point.y(),
                              node.z - axis_point.z());

            // Rodrigues旋转公式:
            // p_rot = p*cos(θ) + (k×p)*sin(θ) + k*(k·p)*(1-cos(θ))
            Eigen::Vector3d k_cross_p = axis_unit_dir.cross(p);
            double k_dot_p = axis_unit_dir.dot(p);

            Eigen::Vector3d p_rot = p * cos_a
                                  + k_cross_p * sin_a
                                  + axis_unit_dir * (k_dot_p * one_minus_cos_a);

            // 平移回原始坐标系
            node.x = p_rot.x() + axis_point.x();
            node.y = p_rot.y() + axis_point.y();
            node.z = p_rot.z() + axis_point.z();
        }
    }

    FEEM_DEBUG("已对 {} 个导体节点绕轴 ({},{},{}) 旋转 {:.6f} rad",
               conductor_node_ids.size(),
               axis_unit_dir.x(), axis_unit_dir.y(), axis_unit_dir.z(), angle);

    return modified_mesh;
}

// ============================================================
// 私有方法：求解静电场并返回能量
// ============================================================

double ElectrostaticForceCalculator::solveForEnergy(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    tool::DimType dim_type
) {
    // 创建静电场求解器实例
    solver::ElectrostaticSolver solver(dim_type);

    // 执行完整的setup→solve→postProcess流程
    if (!solver.setup(mesh_data, materials, boundaries, {})) {
        FEEM_WARN("静电场求解器初始化失败，返回NaN能量值");
        return std::nan("");
    }

    if (!solver.solve()) {
        FEEM_WARN("静电场求解失败，返回NaN能量值");
        return std::nan("");
    }

    if (!solver.postProcess()) {
        FEEM_WARN("静电场后处理失败，返回NaN能量值");
        return std::nan("");
    }

    // 从FieldData中获取静电能量
    const auto& field_data = solver.getFieldData();
    double energy = field_data.computeElectrostaticEnergy(
        solver.getStiffnessMatrix(),
        solver.getFullSolution()
    );

    FEEM_DEBUG("静电场求解完成，能量 W={:.6e} J", energy);

    return energy;
}

} // namespace solver
