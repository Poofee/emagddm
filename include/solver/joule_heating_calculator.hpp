/**
 * @file joule_heating_calculator.hpp
 * @brief 求解器层 - 焦耳热计算器头文件
 * @details 计算导电区域内的焦耳热损耗功率和功率密度分布。
 *          核心公式：P = φ^T K φ（总功率），p = σ|E|²（功率密度）
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include "em_element_integrator_base.hpp"
#include "csr_matrix.hpp"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <utility>

namespace solver {

/**
 * @class JouleHeatingCalculator
 * @brief 焦耳热计算器
 * @details 封装焦耳热损耗的全部计算功能，
 *          包括总功率、单元功率密度分布、区域功率统计等。
 *
 * 物理公式：
 * - 总焦耳热功率：P_total = φ^T · K · φ  （W）
 * - 单元功率密度：p_e = σ |E_e|² = J_e · E_e  （W/m³）
 * - 区域总功率：P_region = Σ_{e∈region} p_e · V_e  （W）
 */
class JouleHeatingCalculator {
public:
    JouleHeatingCalculator();
    ~JouleHeatingCalculator() = default;

    /**
     * @brief 计算总焦耳热功率（使用刚度矩阵和电位向量）
     * @param K 刚度矩阵（CSR格式，已构建）
     * @param phi 自由DOF电位向量
     * @return double 总焦耳热功率（瓦特），矩阵未构建返回0.0
     *
     * @details 公式：P = φ^T · K · φ
     *          与静电能量 W = 0.5 φ^T K φ 的区别在于没有 0.5 系数
     */
    double computeTotalPower(
        const numeric::CsrMatrix<double>& K,
        const Eigen::VectorXd& phi
    );

    /**
     * @brief 计算每个单元的焦耳热功率密度
     * @param elem_E_field 单元电场强度列表
     * @param materials 材料参数映射表（需包含 sigma）
     * @param mesh_data 网格拓扑数据
     * @return std::vector<double> 每个单元的功率密度（W/m³），长度=单元数
     */
    std::vector<double> computeElementPowerDensity(
        const std::vector<Eigen::Vector3d>& elem_E_field,
        const std::map<int, numeric::MaterialProperties>& materials,
        const fe_em::EMMeshData& mesh_data
    );

    /**
     * @brief 计算指定材料区域的总焦耳热功率
     * @param element_power_density 单元功率密度列表
     * @param mesh_data 网格数据
     * @param material_id 目标材料ID
     * @return double 该材料区域的总功率（W）
     */
    double computeRegionPower(
        const std::vector<double>& element_power_density,
        const fe_em::EMMeshData& mesh_data,
        int material_id
    );

    /**
     * @brief 计算整个模型的平均功率密度
     * @param element_power_density 单元功率密度列表
     * @param mesh_data 网格数据
     * @return double 平均功率密度（W/m³）
     */
    double computeAveragePowerDensity(
        const std::vector<double>& element_power_density,
        const fe_em::EMMeshData& mesh_data
    );

    /**
     * @brief 获取最大功率密度的单元索引和值
     * @param element_power_density 单元功率密度列表
     * @return std::pair<int, double> (单元索引, 最大功率密度值)
     */
    std::pair<int, double> findMaxPowerDensity(
        const std::vector<double>& element_power_density
    ) const;

    /**
     * @brief 验证焦耳热一致性：P = V²/R = I²R
     * @param total_power 计算得到的总功率
     * @param voltage 施加的电压差（V）
     * @param resistance 电阻值（Ω）
     * @return double 相对误差 |P_calc - V²/R| / (V²/R)
     */
    double verifyConsistency(
        double total_power,
        double voltage,
        double resistance
    ) const;

    /**
     * @brief 清空内部缓存
     */
    void clear();

private:
    double cached_total_power_;              ///< 缓存的总功率
    std::vector<double> cached_power_density_; ///< 缓存的单元功率密度
    bool has_cached_result_;                 ///< 是否有有效缓存
};

} // namespace solver
