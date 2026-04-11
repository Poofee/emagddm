/**
 * @file capacitance_calculator.hpp
 * @brief 求解器层 - 电容矩阵计算器头文件
 * @details 计算多导体系统的电容矩阵，对标ANSYS Maxwell的电容矩阵输出功能。
 *          支持两种计算方法：
 *          1. 能量法：依次对每个导体施加单位电压，求解N次，提取能量
 *          2. 电荷法：通过高斯面积分计算导体表面电荷
 *
 *          电容矩阵定义（Maxwell电容矩阵）：
 *          - C_ii：导体i的自电容（对地电容）
 *          - C_ij：导体i与导体j之间的互电容（负值）
 *          - 矩阵对称：C_ij = C_ji
 *          - 每行之和：Σ_j C_ij = 导体i对地的总电容
 *          - Q_i = Σ_j C_ij * V_j
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include "em_element_integrator_base.hpp"
#include "em_linear_solver.h"
#include "project_data.hpp"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>

namespace solver {

/**
 * @struct ConductorInfo
 * @brief 导体信息结构体
 */
struct ConductorInfo {
    int id = -1;                     ///< 导体ID（对应boundary_id或region_id）
    std::string name;                ///< 导体名称
    std::vector<int> node_ids;       ///< 导体表面节点ID列表
    double voltage = 0.0;            ///< 导体电压值
    tool::BndType bnd_type = tool::BndType::DIRICHLET; ///< 边界类型
};

/**
 * @class CapacitanceCalculator
 * @brief 电容矩阵计算器
 * @details 采用电荷法计算多导体系统的Maxwell电容矩阵：
 *          1. 对N个导体系统，依次对第i个导体施加单位电压V_i=1V，其余导体接地V_j=0
 *          2. 求解N次静电场问题
 *          3. 每次求解后，通过高斯面积分计算各导体表面电荷
 *          4. 由电荷构建电容矩阵：C_ji = Q_j / V_i = Q_j（V_i=1V时）
 *          5. 强制矩阵对称：C = (C + C^T) / 2
 *
 *          电荷计算公式：
 *          Q_j = -∫_Γj D·n_element dΓ
 *          其中D为电位移矢量，n_element为相邻单元面的外法向（指向导体），
 *          负号将法向从单元外侧转换为导体外侧。
 */
class CapacitanceCalculator {
public:
    CapacitanceCalculator();
    ~CapacitanceCalculator() = default;

    /**
     * @brief 计算多导体系统的电容矩阵（电荷法）
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     * @param boundaries 边界条件列表
     * @param conductor_ids 导体ID列表（对应boundary_id）
     * @param dim_type 维度类型
     * @return Eigen::MatrixXd N×N电容矩阵（法拉），2D返回C/m，3D返回C
     */
    Eigen::MatrixXd computeCapacitanceMatrix(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::vector<int>& conductor_ids,
        tool::DimType dim_type = tool::DimType::D2
    );

    /**
     * @brief 计算两个导体之间的电容（互电容）
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     * @param boundaries 边界条件列表
     * @param conductor1_id 导体1的ID
     * @param conductor2_id 导体2的ID
     * @param dim_type 维度类型
     * @return double 互电容值（法拉），正值
     */
    double computeCapacitance(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        int conductor1_id,
        int conductor2_id,
        tool::DimType dim_type = tool::DimType::D2
    );

    /**
     * @brief 获取最近一次计算的导体信息列表
     * @return const std::vector<ConductorInfo>& 导体信息列表
     */
    const std::vector<ConductorInfo>& getConductorInfos() const;

    /**
     * @brief 清空内部数据
     */
    void clear();

private:
    std::vector<ConductorInfo> conductor_infos_;   ///< 导体信息列表
    Eigen::MatrixXd capacitance_matrix_;            ///< 电容矩阵缓存
    std::vector<tool::Boundary> original_boundaries_; ///< 原始边界条件列表缓存

    /**
     * @brief 从边界条件中提取导体信息
     * @param boundaries 边界条件列表
     * @param mesh_data 网格拓扑数据
     * @param conductor_ids 导体ID列表
     * @return bool 提取成功返回true
     */
    bool extractConductorInfos(
        const std::vector<tool::Boundary>& boundaries,
        const fe_em::EMMeshData& mesh_data,
        const std::vector<int>& conductor_ids
    );

    /**
     * @brief 构建第i个导体加单位电压时的边界条件列表
     * @param conductor_index 当前加电压的导体索引
     * @param num_conductors 导体总数
     * @return std::vector<tool::Boundary> 修改后的边界条件列表
     */
    std::vector<tool::Boundary> buildBoundaryConditions(
        int conductor_index,
        int num_conductors
    ) const;
};

} // namespace solver
