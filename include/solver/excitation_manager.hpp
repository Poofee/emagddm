/**
 * @file excitation_manager.hpp
 * @brief 求解器层 - 激励管理器头文件
 * @details 将tool::Excitation转化为右端项向量贡献或等效边界条件。
 *          静电场支持的激励类型：
 *          - VOLTAGE_SOURCE → 转化为Dirichlet边界条件
 *          - CURRENT_DENSITY → 转化为右端项体积分贡献（体电荷密度源）
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "em_enums.hpp"
#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"
#include "project_data.hpp"
#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>

namespace solver {

/**
 * @class ExcitationManager
 * @brief 激励管理器
 * @details 负责将tool::Excitation转化为有限元系统的激励施加：
 *          1. 电压源激励 → 转化为Dirichlet边界条件
 *          2. 体电荷密度源 → 转化为右端项体积分贡献
 *          3. 激励合法性校验
 */
class ExcitationManager {
public:
    ExcitationManager();
    ~ExcitationManager() = default;

    /**
     * @brief 处理激励列表，生成等效边界条件
     * @param excitations 激励源列表
     * @param mesh_data 网格拓扑数据
     * @return std::vector<tool::Boundary> 等效边界条件列表
     *         （电压源激励转化为Dirichlet边界条件）
     */
    std::vector<tool::Boundary> processExcitations(
        const std::vector<tool::Excitation>& excitations,
        const fe_em::EMMeshData& mesh_data
    );

    /**
     * @brief 计算体电荷密度源对右端项的贡献
     * @details F_i += ∫ N_i * ρ dΩ
     * @param rhs 右端项向量（会被修改）
     * @param mesh_data 网格拓扑数据
     * @param elem_l2g 单元局部-全局DOF映射表
     * @param charge_density 体电荷密度映射表（region_id → ρ值）
     */
    void applyVolumeChargeSource(
        Eigen::VectorXd& rhs,
        const fe_em::EMMeshData& mesh_data,
        const std::vector<fe_em::Local2Global>& elem_l2g,
        const std::map<int, double>& charge_density
    );

    /**
     * @brief 校验激励合法性
     * @param excitations 激励源列表
     * @return bool 全部合法返回true
     */
    bool validateExcitations(
        const std::vector<tool::Excitation>& excitations
    ) const;

    /**
     * @brief 清空所有数据
     */
    void clear();

    /**
     * @brief 获取内部体电荷密度映射表
     * @return const std::map<int, double>& region_id → ρ值的只读引用
     */
    const std::map<int, double>& getVolumeChargeDensity() const;

private:
    std::map<int, double> volume_charge_density_;  ///< 体电荷密度映射（region_id → ρ）

    /**
     * @brief 根据激励名称在网格边界标记中查找对应的region_id
     * @param name 激励名称
     * @param mesh_data 网格拓扑数据
     * @return int 匹配的region_id，未找到返回-1
     */
    int findRegionIdByName(
        const std::string& name,
        const fe_em::EMMeshData& mesh_data
    ) const;

    /**
     * @brief 将ElementType枚举转换为字符串标识符
     * @param type 单元类型枚举值
     * @return std::string 对应的字符串标识符，不支持返回空字符串
     */
    std::string elementTypeToString(numeric::ElementType type) const;

    /**
     * @brief 根据单元类型获取默认高斯积分阶数
     * @param type 单元类型枚举值
     * @return int 默认积分阶数
     */
    int getDefaultGaussOrder(numeric::ElementType type) const;

    /**
     * @brief 从网格数据中提取单元节点坐标矩阵
     * @param elem 单元数据
     * @param mesh_data 网格拓扑数据
     * @return Eigen::MatrixXd 节点坐标矩阵，维度为 dim × node_count
     */
    Eigen::MatrixXd extractNodeCoords(
        const fe_em::Element& elem,
        const fe_em::EMMeshData& mesh_data
    ) const;
};

} // namespace solver
