/**
 * @file joule_heating_calculator.cpp
 * @brief 求解器层 - 焦耳热计算器实现文件
 * @details 实现焦耳热功率和功率密度的全部计算功能。
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "joule_heating_calculator.hpp"
#include "logger_factory.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace solver {

JouleHeatingCalculator::JouleHeatingCalculator()
    : cached_total_power_(0.0)
    , has_cached_result_(false)
{
}

double JouleHeatingCalculator::computeTotalPower(
    const numeric::CsrMatrix<double>& K,
    const Eigen::VectorXd& phi
) {
    if (!K.is_built()) {
        FEEM_WARN("computeTotalPower: 刚度矩阵未构建，返回0.0");
        return 0.0;
    }

    if (phi.size() == 0) {
        FEEM_WARN("computeTotalPower: 电位向量为空，返回0.0");
        return 0.0;
    }

    if (static_cast<int>(phi.size()) != K.rows()) {
        FEEM_ERROR("computeTotalPower: 维度不匹配 phi.size()={} != K.rows()={}",
                   phi.size(), K.rows());
        return 0.0;
    }

    try {
        const auto& eigen_K = K.get_eigen_real();
        Eigen::VectorXd y = eigen_K * phi;
        double total_power = phi.dot(y);

        cached_total_power_ = total_power;
        has_cached_result_ = true;

        FEEM_INFO("computeTotalPower: 总焦耳热功率 P = {:.6e} W", total_power);
        return total_power;
    } catch (const std::exception& e) {
        FEEM_ERROR("computeTotalPower: 矩阵向量运算异常 - {}", e.what());
        return 0.0;
    }
}

std::vector<double> JouleHeatingCalculator::computeElementPowerDensity(
    const std::vector<Eigen::Vector3d>& elem_E_field,
    const std::map<int, numeric::MaterialProperties>& materials,
    const fe_em::EMMeshData& mesh_data
) {
    size_t num_elements = mesh_data.getElementCount();

    if (elem_E_field.size() != num_elements) {
        FEEM_ERROR("computeElementPowerDensity: 电场列表长度{}与单元数{}不匹配",
                   elem_E_field.size(), num_elements);
        return {};
    }

    std::vector<double> power_density(num_elements, 0.0);

    for (size_t e = 0; e < num_elements; ++e) {
        int mat_id = mesh_data.elements[e].material_id;

        auto it = materials.find(mat_id);
        if (it == materials.end()) {
            FEEM_DEBUG("computeElementPowerDensity: 单元{} material_id={}未找到对应材料",
                       e, mat_id);
            continue;
        }

        double sigma = it->second.sigma;
        double E_squared = elem_E_field[e].squaredNorm();
        power_density[e] = sigma * E_squared;
    }

    cached_power_density_ = power_density;
    has_cached_result_ = true;

    double max_pd = *std::max_element(power_density.begin(), power_density.end());
    double sum_pd = std::accumulate(power_density.begin(), power_density.end(), 0.0);
    FEEM_INFO("computeElementPowerDensity: 单元功率密度计算完成, "
              "单元数={}, Σp={:.6e} W/m³, max(p)={:.6e} W/m³",
              num_elements, sum_pd, max_pd);

    return power_density;
}

double JouleHeatingCalculator::computeRegionPower(
    const std::vector<double>& element_power_density,
    const fe_em::EMMeshData& mesh_data,
    int material_id
) {
    size_t num_elements = mesh_data.getElementCount();

    if (element_power_density.size() != num_elements) {
        FEEM_ERROR("computeRegionPower: 功率密度列表长度{}与单元数{}不匹配",
                   element_power_density.size(), num_elements);
        return 0.0;
    }

    double region_power = 0.0;
    int element_count = 0;

    for (size_t e = 0; e < num_elements; ++e) {
        if (mesh_data.elements[e].material_id == material_id) {
            region_power += element_power_density[e];
            ++element_count;
        }
    }

    FEEM_INFO("computeRegionPower: material_id={} 区域总功率 P={:.6e} W ({}个单元)",
              material_id, region_power, element_count);

    return region_power;
}

double JouleHeatingCalculator::computeAveragePowerDensity(
    const std::vector<double>& element_power_density,
    const fe_em::EMMeshData& mesh_data
) {
    size_t num_elements = mesh_data.getElementCount();

    if (element_power_density.empty() || num_elements == 0) {
        FEEM_WARN("computeAveragePowerDensity: 输入数据为空，返回0.0");
        return 0.0;
    }

    double sum = std::accumulate(element_power_density.begin(),
                                 element_power_density.end(), 0.0);
    double avg = sum / static_cast<double>(num_elements);

    FEEM_INFO("computeAveragePowerDensity: 平均功率密度 p_avg={:.6e} W/m³ ({}个单元)",
              avg, num_elements);

    return avg;
}

std::pair<int, double> JouleHeatingCalculator::findMaxPowerDensity(
    const std::vector<double>& element_power_density
) const {
    if (element_power_density.empty()) {
        FEEM_WARN("findMaxPowerDensity: 功率密度列表为空");
        return {-1, 0.0};
    }

    auto max_it = std::max_element(element_power_density.begin(),
                                   element_power_density.end());
    int index = static_cast<int>(std::distance(element_power_density.begin(), max_it));

    FEEM_DEBUG("findMaxPowerDensity: 最大功率密度位于单元{}, 值={:.6e} W/m³",
               index, *max_it);

    return {index, *max_it};
}

double JouleHeatingCalculator::verifyConsistency(
    double total_power,
    double voltage,
    double resistance
) const {
    if (resistance <= 0.0) {
        FEEM_ERROR("verifyConsistency: 电阻值必须为正数, R={}", resistance);
        return -1.0;
    }

    double expected_power = voltage * voltage / resistance;
    double relative_error = std::abs(total_power - expected_power) / expected_power;

    FEEM_INFO("verifyConsistency: P_calc={:.6e}W, P_expected=V²/R={:.6e}W "
              "(V={:.2f}V, R={:.6e}Ω), 相对误差={:.4f}%",
              total_power, expected_power, voltage, resistance,
              relative_error * 100.0);

    return relative_error;
}

void JouleHeatingCalculator::clear() {
    cached_total_power_ = 0.0;
    cached_power_density_.clear();
    has_cached_result_ = false;
    FEEM_DEBUG("JouleHeatingCalculator: 内部缓存已清空");
}

} // namespace solver
