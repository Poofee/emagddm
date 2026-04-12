/**
 * @file field_data.hpp
 * @brief 求解器层 - 场数据与后处理模块头文件
 * @details 存储有限元求解结果数据，并提供常用物理量的后处理计算功能，
 *          包括电场强度、电位移矢量、静电能量、电流密度、焦耳热等。
 *          支持VTK和CSV格式的结果导出。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.1
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>

namespace fe_em { class EMMeshData; }
namespace numeric {
    struct MaterialProperties;
    template<typename T> class CsrMatrix;
}

namespace solver {

/**
 * @class FieldData
 * @brief 场数据存储与后处理计算类
 * @details 封装有限元求解结果的存储和派生物理量的计算功能。
 *          数据组织方式：
 *          - 节点电位：长度为num_total_nodes的向量
 *          - 单元电场强度：每个单元一个Vector3d
 *          - 单元电位移矢量：每个单元一个Vector3d
 *          - 静电能量：标量值
 */
class FieldData {
public:
    FieldData();
    ~FieldData() = default;

    // ========== 节点电位 ==========

    /**
     * @brief 设置节点电位值
     * @param phi 节点电位向量（包含约束DOF的完整向量）
     * @param num_total_dofs 总DOF数量
     */
    void setNodalPotential(const Eigen::VectorXd& phi, int num_total_dofs);

    /**
     * @brief 获取节点电位（只读）
     * @return const Eigen::VectorXd& 节点电位向量
     */
    const Eigen::VectorXd& getNodalPotential() const;

    // ========== 电场强度计算 ==========

    /**
     * @brief 计算单元电场强度 E = -∇φ
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     */
    void computeElectricField(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials
    );

    /**
     * @brief 获取单元电场强度（只读）
     * @return const std::vector<Eigen::Vector3d>& 每个单元的电场强度
     */
    const std::vector<Eigen::Vector3d>& getElementElectricField() const;

    // ========== 电位移矢量计算 ==========

    /**
     * @brief 计算单元电位移矢量 D = εE
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表
     */
    void computeElectricDisplacement(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials
    );

    /**
     * @brief 获取单元电位移矢量（只读）
     * @return const std::vector<Eigen::Vector3d>& 每个单元的电位移矢量
     */
    const std::vector<Eigen::Vector3d>& getElementElectricDisplacement() const;

    // ========== 静电能量计算 ==========

    /**
     * @brief 计算静电能量 W = 0.5 * φ^T * K * φ
     * @param K 刚度矩阵（CSR格式）
     * @param phi 自由DOF电位向量
     * @return double 静电能量（焦耳）
     */
    double computeElectrostaticEnergy(
        const numeric::CsrMatrix<double>& K,
        const Eigen::VectorXd& phi
    ) const;

    // ========== 电流密度计算 ==========

    /**
     * @brief 计算单元电流密度 J = σE
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表（需包含 sigma 字段）
     */
    void computeCurrentDensity(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials
    );

    /**
     * @brief 获取单元电流密度（只读）
     * @return const std::vector<Eigen::Vector3d>& 每个单元的电流密度向量 J (A/m²)
     */
    const std::vector<Eigen::Vector3d>& getElementCurrentDensity() const;

    // ========== 焦耳热计算 ==========

    /**
     * @brief 计算总焦耳热功率 P = φ^T K φ
     * @param K 刚度矩阵（CSR格式）
     * @param phi 自由DOF电位向量
     * @return double 总焦耳热功率（瓦特）
     */
    double computeJouleHeating(
        const numeric::CsrMatrix<double>& K,
        const Eigen::VectorXd& phi
    ) const;

    /**
     * @brief 计算每个单元的焦耳热功率密度 p = σ|E|²
     * @details 结果存入 elem_joule_power_density_，需先调用 computeCurrentDensity()
     */
    void computeElementPowerDensity();

    /**
     * @brief 获取单元焦耳热功率密度（只读）
     * @return const std::vector<double>& 每个单元的焦耳热功率密度 p (W/m³)
     */
    const std::vector<double>& getElementPowerDensity() const;

    // ========== 总电流计算 ==========

    /**
     * @brief 计算通过指定区域边界的总电流 I = ∫ J·n dS
     * @param mesh_data 网格拓扑数据
     * @param elem_J 单元电流密度向量列表
     * @return double 总电流值（安培），初版返回所有单元电流密度的统计近似值
     */
    double computeTotalCurrent(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<Eigen::Vector3d>& elem_J
    ) const;

    // ========== 电压降计算 ==========

    /**
     * @brief 计算两节点间的电压降
     * @param node_id_1 节点1 ID
     * @param node_id_2 节点2 ID
     * @return double 电压降 V(node₁) - V(node₂)（伏特）
     */
    double computeVoltageDrop(int node_id_1, int node_id_2) const;

    // ========== 结果导出 ==========

    /**
     * @brief 导出结果为VTK格式（Legacy VTK）
     * @param filename 输出文件路径
     * @param mesh_data 网格拓扑数据
     * @return bool 导出成功返回true
     */
    bool exportVTK(const std::string& filename,
                   const fe_em::EMMeshData& mesh_data) const;

    /**
     * @brief 导出结果为CSV格式
     * @param filename 输出文件路径
     * @param mesh_data 网格拓扑数据
     * @return bool 导出成功返回true
     */
    bool exportCSV(const std::string& filename,
                   const fe_em::EMMeshData& mesh_data) const;

    // ========== 资源管理 ==========

    /**
     * @brief 清空所有数据
     */
    void clear();

    /**
     * @brief 检查是否有有效的电位数据
     * @return bool 有有效数据返回true
     */
    bool hasData() const;

private:
    Eigen::VectorXd nodal_potential_;                          ///< 节点电位向量 φ
    std::vector<Eigen::Vector3d> elem_electric_field_;         ///< 单元电场强度 E = -∇φ
    std::vector<Eigen::Vector3d> elem_electric_displacement_;  ///< 单元电位移矢量 D = εE
    std::vector<Eigen::Vector3d> elem_current_density_;        ///< 单元电流密度 J=σE (A/m²)
    std::vector<double> elem_joule_power_density_;             ///< 单元焦耳热功率密度 p=σ|E|² (W/m³)
    double electrostatic_energy_;                              ///< 静电能量 W (J)
    double total_current_;                                     ///< 总电流 I (A)
    double joule_heating_total_;                               ///< 总焦耳热功率 P (W)
    bool has_potential_;                                       ///< 是否已设置节点电位
    bool has_electric_field_;                                  ///< 是否已计算电场强度
    bool has_electric_displacement_;                           ///< 是否已计算电位移矢量
    bool has_current_density_;                                 ///< 是否已计算电流密度
    bool has_joule_heating_;                                   ///< 是否已计算焦耳热
};

} // namespace solver
