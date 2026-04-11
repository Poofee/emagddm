/**
 * @file electrostatic_force_calculator.hpp
 * @brief 求解器层 - 静电力计算器头文件
 * @details 采用虚功法（Virtual Work Method）计算导体上的静电力和力矩，
 *          对标ANSYS Maxwell的力计算功能。
 *
 *          虚功法原理：
 *          F = -∂W/∂s，数值实现：F ≈ -(W(s+δ) - W(s-δ)) / (2δ)
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include "em_element_integrator_base.hpp"
#include "project_data.hpp"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>

namespace solver {

/**
 * @class ElectrostaticForceCalculator
 * @brief 静电力计算器
 * @details 采用虚功法计算导体上的静电力和力矩。
 *          核心思路：对导体表面节点施加微小位移，重新求解，
 *          通过能量差分计算力/力矩。
 */
class ElectrostaticForceCalculator {
public:
    ElectrostaticForceCalculator();
    ~ElectrostaticForceCalculator() = default;

    /**
     * @brief 计算导体上的静电力（虚功法）
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     * @param boundaries 边界条件列表
     * @param conductor_name 导体名称（匹配Boundary::getName()）
     * @param dim_type 维度类型
     * @param delta 虚位移量（默认1e-6 m）
     * @return Eigen::Vector3d 力向量（牛顿），2D时z=0
     */
    Eigen::Vector3d computeForce(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::string& conductor_name,
        tool::DimType dim_type = tool::DimType::D2,
        double delta = 1e-6
    );

    /**
     * @brief 计算导体上的力矩（虚功法）
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     * @param boundaries 边界条件列表
     * @param conductor_name 导体名称
     * @param axis_point 旋转轴上一点
     * @param axis_direction 旋转轴方向（单位向量）
     * @param dim_type 维度类型
     * @param delta 虚角位移量（弧度，默认1e-6）
     * @return double 力矩值（牛顿·米）
     */
    double computeTorque(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::string& conductor_name,
        const Eigen::Vector3d& axis_point,
        const Eigen::Vector3d& axis_direction,
        tool::DimType dim_type = tool::DimType::D2,
        double delta = 1e-6
    );

    /**
     * @brief 清空内部数据
     */
    void clear();

private:
    /**
     * @brief 从边界条件中获取导体节点ID列表
     * @param boundaries 边界条件列表
     * @param mesh_data 网格数据（用于通过边界标记查找节点）
     * @param conductor_name 导体名称
     * @return std::vector<int> 节点ID列表
     */
    static std::vector<int> getConductorNodeIds(
        const std::vector<tool::Boundary>& boundaries,
        const fe_em::EMMeshData& mesh_data,
        const std::string& conductor_name
    );

    /**
     * @brief 对网格中指定导体的节点施加位移（深拷贝+修改）
     * @param mesh_data 原始网格数据
     * @param conductor_node_ids 导体节点ID列表
     * @param displacement 位移向量
     * @return fe_em::EMMeshData 修改后的网格数据（深拷贝）
     */
    static fe_em::EMMeshData displaceConductorNodes(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<int>& conductor_node_ids,
        const Eigen::Vector3d& displacement
    );

    /**
     * @brief 对网格中指定导体的节点施加旋转变换（Rodrigues旋转公式）
     * @param mesh_data 原始网格数据
     * @param conductor_node_ids 导体节点ID列表
     * @param axis_point 旋转轴上一点
     * @param axis_unit_dir 旋转轴单位方向向量
     * @param angle 旋转角度（弧度）
     * @return fe_em::EMMeshData 修改后的网格数据（深拷贝）
     */
    static fe_em::EMMeshData rotateConductorNodes(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<int>& conductor_node_ids,
        const Eigen::Vector3d& axis_point,
        const Eigen::Vector3d& axis_unit_dir,
        double angle
    );

    /**
     * @brief 求解静电场并返回静电能量
     * @param mesh_data 网格数据
     * @param materials 材料参数
     * @param boundaries 边界条件
     * @param dim_type 维度类型
     * @return double 静电能量（焦耳），求解失败返回NaN
     */
    static double solveForEnergy(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        tool::DimType dim_type
    );
};

} // namespace solver
