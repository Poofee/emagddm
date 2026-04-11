/**
 * @file excitation_manager.cpp
 * @brief 求解器层 - 激励管理器实现
 * @details 实现ExcitationManager类，将tool::Excitation转化为有限元系统的
 *          右端项向量贡献或等效边界条件。
 *
 *          静电场激励类型映射：
 *          - VOLTAGE_SOURCE → tool::Boundary(DIRICHLET, voltage=value)
 *          - CURRENT_DENSITY → volume_charge_density_[region_id] = ρ
 *            （在静电场中，电流密度激励等效为体电荷密度源）
 *
 *          体电荷密度源右端项贡献公式：
 *          F_i += ∫_Ω N_i * ρ dΩ = Σ_q w_q * |detJ_q| * ρ * N_i(ξ_q)
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "excitation_manager.hpp"
#include "logger_factory.hpp"
#include "shape_function_factory.hpp"
#include "gauss_quadrature.hpp"
#include <cmath>
#include <algorithm>

namespace solver {

// ==================== 构造函数 ====================

/**
 * @brief 默认构造函数
 * @details 初始化空的体电荷密度映射表
 */
ExcitationManager::ExcitationManager()
    : volume_charge_density_()
{
}

// ==================== 核心方法实现 ====================

/**
 * @brief 处理激励列表，生成等效边界条件
 * @details 遍历所有激励，根据激励类型进行转化：
 *          - VOLTAGE_SOURCE: 创建Dirichlet边界条件，设置电压值
 *          - CURRENT_DENSITY: 在静电场中映射为体电荷密度源，
 *            通过匹配网格边界标记的名称确定region_id
 *          - 其他类型: 记录警告日志，静电场不支持
 *
 * @param excitations 激励源列表
 * @param mesh_data 网格拓扑数据（用于查找region_id映射）
 * @return std::vector<tool::Boundary> 生成的等效边界条件列表
 */
std::vector<tool::Boundary> ExcitationManager::processExcitations(
    const std::vector<tool::Excitation>& excitations,
    const fe_em::EMMeshData& mesh_data
) {
    std::vector<tool::Boundary> boundaries;
    volume_charge_density_.clear();

    for (const auto& excitation : excitations) {
        switch (excitation.getType()) {
            case tool::ExcitationType::VOLTAGE_SOURCE: {
                tool::Boundary bnd(excitation.getName());
                bnd.setType(tool::BndType::DIRICHLET);
                bnd.setVoltage(excitation.getValue());
                boundaries.push_back(bnd);
                FEEM_INFO("ExcitationManager: 电压源激励 '{}' 转化为Dirichlet边界条件, "
                          "电压值={}", excitation.getName(), excitation.getValue());
                break;
            }

            case tool::ExcitationType::CURRENT_DENSITY: {
                double rho = excitation.getValue();
                int region_id = findRegionIdByName(excitation.getName(), mesh_data);
                if (region_id >= 0) {
                    volume_charge_density_[region_id] = rho;
                    FEEM_INFO("ExcitationManager: 电流密度激励 '{}' 映射为体电荷密度源, "
                              "region_id={}, ρ={}", excitation.getName(), region_id, rho);
                } else {
                    FEEM_WARN("ExcitationManager: 电流密度激励 '{}' 无法映射到网格区域, "
                              "请检查激励名称与网格边界标记的对应关系",
                              excitation.getName());
                }
                break;
            }

            default: {
                FEEM_WARN("ExcitationManager: 激励 '{}' 类型为 {}, "
                          "静电场分析不支持该激励类型，已跳过",
                          excitation.getName(),
                          tool::excitationTypeToString(excitation.getType()));
                break;
            }
        }
    }

    FEEM_INFO("ExcitationManager: 激励处理完成, 生成{}个边界条件, "
              "{}个体电荷密度源",
              boundaries.size(), volume_charge_density_.size());

    return boundaries;
}

/**
 * @brief 计算体电荷密度源对右端项的贡献
 * @details 对每个单元，检查其region_id是否在charge_density映射表中。
 *          若存在，则计算源项贡献：F_i += ∫ N_i * ρ dΩ
 *          使用高斯数值积分完成参考域到物理域的积分计算。
 *
 *          数值积分公式：
 *          F_i = Σ_{q=1}^{n_gp} w_q * |detJ(ξ_q)| * ρ * N_i(ξ_q)
 *
 *          其中：
 *          - w_q: 第q个高斯积分点的权重
 *          - detJ(ξ_q): 雅可比行列式，用于体积元变换
 *          - ρ: 体电荷密度（常数，由charge_density映射表提供）
 *          - N_i(ξ_q): 第i个形函数在积分点处的值
 *
 * @param rhs 右端项向量（会被原地修改，累加源项贡献）
 * @param mesh_data 网格拓扑数据（提供节点坐标和单元连接关系）
 * @param elem_l2g 单元局部-全局DOF映射表
 * @param charge_density 体电荷密度映射表（region_id → ρ值）
 */
void ExcitationManager::applyVolumeChargeSource(
    Eigen::VectorXd& rhs,
    const fe_em::EMMeshData& mesh_data,
    const std::vector<fe_em::Local2Global>& elem_l2g,
    const std::map<int, double>& charge_density
) {
    if (charge_density.empty()) {
        FEEM_DEBUG("ExcitationManager: 体电荷密度映射表为空, 跳过源项计算");
        return;
    }

    if (mesh_data.getElementCount() != elem_l2g.size()) {
        FEEM_ERROR("ExcitationManager: 单元数量({})与L2G映射表大小({})不匹配",
                   mesh_data.getElementCount(), elem_l2g.size());
        return;
    }

    size_t applied_count = 0;

    for (size_t e = 0; e < mesh_data.getElementCount(); ++e) {
        const auto& elem = mesh_data.elements[e];
        const auto& l2g = elem_l2g[e];

        auto it = charge_density.find(elem.region_id);
        if (it == charge_density.end()) {
            continue;
        }

        double rho = it->second;
        if (std::abs(rho) < 1e-30) {
            continue;
        }

        std::string type_str = elementTypeToString(elem.elem_type);
        if (type_str.empty()) {
            FEEM_WARN("ExcitationManager: 单元{}的类型({})不支持, 跳过源项计算",
                      elem.id, static_cast<int>(elem.elem_type));
            continue;
        }

        auto shape_func = numeric::ShapeFunctionFactory::create(type_str);
        if (!shape_func) {
            FEEM_WARN("ExcitationManager: 单元{}的形函数创建失败(type={}), 跳过",
                      elem.id, type_str);
            continue;
        }

        int gauss_order = getDefaultGaussOrder(elem.elem_type);
        auto gauss_points = numeric::GaussQuadrature::getPoints(elem.elem_type, gauss_order);
        if (gauss_points.empty()) {
            FEEM_WARN("ExcitationManager: 单元{}的高斯积分点获取失败, 跳过", elem.id);
            continue;
        }

        Eigen::MatrixXd node_coords = extractNodeCoords(elem, mesh_data);

        int n_nodes = shape_func->getNodeCount();
        Eigen::VectorXd local_rhs = Eigen::VectorXd::Zero(n_nodes);

        for (const auto& gp : gauss_points) {
            numeric::LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);
            if (gp.dim == 2) {
                xi.dim = 2;
            }

            Eigen::VectorXd N = shape_func->evalN(xi);

            numeric::JacobianResult jaco = shape_func->calcJacobian(xi, node_coords);
            double det_j = jaco.det_j;

            double dv = gp.weight * std::abs(det_j);
            local_rhs += dv * rho * N;
        }

        for (int i = 0; i < n_nodes && i < static_cast<int>(l2g.indices.size()); ++i) {
            int global_dof = l2g.indices[i];
            if (global_dof >= 0 && global_dof < rhs.size()) {
                rhs(global_dof) += local_rhs(i);
            }
        }

        ++applied_count;
    }

    FEEM_INFO("ExcitationManager: 体电荷密度源项施加完成, "
              "影响{}个单元", applied_count);
}

/**
 * @brief 校验激励合法性
 * @details 检查内容：
 *          1. 电压源值不能为NaN
 *          2. 电荷密度值必须有限（非NaN、非Inf）
 *          3. 目标实体引用有效性检查（通过名称匹配网格边界标记）
 *          4. 不支持的激励类型记录警告
 *
 * @param excitations 激励源列表
 * @return bool 全部合法返回true，存在非法项返回false
 */
bool ExcitationManager::validateExcitations(
    const std::vector<tool::Excitation>& excitations
) const {
    bool all_valid = true;

    for (const auto& excitation : excitations) {
        double value = excitation.getValue();

        switch (excitation.getType()) {
            case tool::ExcitationType::VOLTAGE_SOURCE: {
                if (std::isnan(value)) {
                    FEEM_ERROR("ExcitationManager: 电压源激励 '{}' 的值为NaN, 非法",
                               excitation.getName());
                    all_valid = false;
                }
                break;
            }

            case tool::ExcitationType::CURRENT_DENSITY: {
                if (!std::isfinite(value)) {
                    FEEM_ERROR("ExcitationManager: 电流密度激励 '{}' 的值非有限"
                               "（NaN或Inf）, 非法", excitation.getName());
                    all_valid = false;
                }
                break;
            }

            default: {
                FEEM_WARN("ExcitationManager: 激励 '{}' 类型为 {}, "
                          "静电场分析不支持该类型",
                          excitation.getName(),
                          tool::excitationTypeToString(excitation.getType()));
                break;
            }
        }
    }

    if (all_valid) {
        FEEM_INFO("ExcitationManager: 激励合法性校验通过, 共{}项",
                  excitations.size());
    } else {
        FEEM_ERROR("ExcitationManager: 激励合法性校验失败, 存在非法项");
    }

    return all_valid;
}

/**
 * @brief 清空所有数据
 * @details 清空内部体电荷密度映射表
 */
void ExcitationManager::clear() {
    volume_charge_density_.clear();
    FEEM_DEBUG("ExcitationManager: 数据已清空");
}

/**
 * @brief 获取内部体电荷密度映射表
 * @return const std::map<int, double>& region_id → ρ值的只读引用
 */
const std::map<int, double>& ExcitationManager::getVolumeChargeDensity() const {
    return volume_charge_density_;
}

// ==================== 私有辅助方法 ====================

/**
 * @brief 根据激励名称在网格边界标记中查找对应的region_id
 * @details 遍历mesh_data.boundary_markers，查找名称匹配的标记，
 *          从标记的target_ids中提取节点ID，再从节点中获取region_id。
 *          若无法确定region_id，则遍历单元查找包含该标记节点的单元region_id。
 *
 * @param name 激励名称
 * @param mesh_data 网格拓扑数据
 * @return int 匹配的region_id，未找到返回-1
 */
int ExcitationManager::findRegionIdByName(
    const std::string& name,
    const fe_em::EMMeshData& mesh_data
) const {
    for (const auto& marker : mesh_data.boundary_markers) {
        if (marker.name == name) {
            if (std::holds_alternative<std::vector<int>>(marker.target_ids)) {
                const auto& node_ids = std::get<std::vector<int>>(marker.target_ids);
                if (!node_ids.empty()) {
                    int first_node_id = node_ids[0];
                    for (const auto& node : mesh_data.nodes) {
                        if (node.id == first_node_id && node.region_id >= 0) {
                            return node.region_id;
                        }
                    }
                    for (const auto& elem : mesh_data.elements) {
                        auto it = std::find(elem.node_ids.begin(),
                                            elem.node_ids.end(), first_node_id);
                        if (it != elem.node_ids.end()) {
                            return elem.region_id;
                        }
                    }
                }
            }
        }
    }

    return -1;
}

/**
 * @brief 将ElementType枚举转换为字符串标识符
 * @details 用于ShapeFunctionFactory::create()的参数。
 *          仅支持Lagrange节点元类型（静电场标量位问题）。
 *
 * @param type 单元类型枚举值
 * @return std::string 对应的字符串标识符，不支持返回空字符串
 */
std::string ExcitationManager::elementTypeToString(
    numeric::ElementType type
) const {
    switch (type) {
        case numeric::ElementType::LINE2:  return "LINE2";
        case numeric::ElementType::LINE3:  return "LINE3";
        case numeric::ElementType::TRI3:   return "TRI3";
        case numeric::ElementType::TRI6:   return "TRI6";
        case numeric::ElementType::QUAD4:  return "QUAD4";
        case numeric::ElementType::QUAD8:  return "QUAD8";
        case numeric::ElementType::QUAD9:  return "QUAD9";
        case numeric::ElementType::TET4:   return "TET4";
        case numeric::ElementType::TET10:  return "TET10";
        case numeric::ElementType::HEX8:   return "HEX8";
        case numeric::ElementType::HEX20:  return "HEX20";
        case numeric::ElementType::HEX27:  return "HEX27";
        case numeric::ElementType::PRISM6: return "PRISM6";
        case numeric::ElementType::PRISM15:return "PRISM15";
        case numeric::ElementType::PYRAMID5:  return "PYRAMID5";
        case numeric::ElementType::PYRAMID13: return "PYRAMID13";
        default: return "";
    }
}

/**
 * @brief 根据单元类型获取默认高斯积分阶数
 * @details 各单元类型的默认积分方案与ElectrostaticIntegrator保持一致：
 *          - TRI3: 1点（质心积分）
 *          - QUAD4: 4点（2×2 Gauss-Legendre）
 *          - TET4: 1点（质心积分）
 *          - HEX8: 8点（2×2×2 Gauss-Legendre）
 *          - PRISM6: 6点（三角形3点×ζ方向2点）
 *          - PYRAMID5: 5点
 *
 * @param type 单元类型枚举值
 * @return int 默认积分阶数，不支持返回1
 */
int ExcitationManager::getDefaultGaussOrder(
    numeric::ElementType type
) const {
    switch (type) {
        case numeric::ElementType::TRI3:   return 1;
        case numeric::ElementType::TRI6:   return 3;
        case numeric::ElementType::QUAD4:  return 4;
        case numeric::ElementType::QUAD8:  return 4;
        case numeric::ElementType::QUAD9:  return 4;
        case numeric::ElementType::TET4:   return 1;
        case numeric::ElementType::TET10:  return 4;
        case numeric::ElementType::HEX8:   return 8;
        case numeric::ElementType::HEX20:  return 8;
        case numeric::ElementType::HEX27:  return 8;
        case numeric::ElementType::PRISM6: return 6;
        case numeric::ElementType::PRISM15:return 6;
        case numeric::ElementType::PYRAMID5:  return 5;
        case numeric::ElementType::PYRAMID13: return 5;
        case numeric::ElementType::LINE2:  return 2;
        case numeric::ElementType::LINE3:  return 3;
        default: return 1;
    }
}

/**
 * @brief 从网格数据中提取单元节点坐标矩阵
 * @details 构建dim × node_count的坐标矩阵，每列对应一个节点的物理坐标。
 *          二维节点z坐标为0，三维节点使用实际z坐标。
 *
 * @param elem 单元数据（包含node_ids列表）
 * @param mesh_data 网格拓扑数据（包含节点坐标信息）
 * @return Eigen::MatrixXd 节点坐标矩阵，维度为 dim × node_count
 */
Eigen::MatrixXd ExcitationManager::extractNodeCoords(
    const fe_em::Element& elem,
    const fe_em::EMMeshData& mesh_data
) const {
    int n_nodes = static_cast<int>(elem.node_ids.size());
    int dim = 3;

    Eigen::MatrixXd coords(dim, n_nodes);

    for (int i = 0; i < n_nodes; ++i) {
        int node_id = elem.node_ids[i];
        bool found = false;
        for (const auto& node : mesh_data.nodes) {
            if (node.id == node_id) {
                coords(0, i) = node.x;
                coords(1, i) = node.y;
                coords(2, i) = node.z;
                found = true;
                break;
            }
        }
        if (!found) {
            FEEM_WARN("ExcitationManager: 节点ID={}未找到, 坐标设为零", node_id);
            coords(0, i) = 0.0;
            coords(1, i) = 0.0;
            coords(2, i) = 0.0;
        }
    }

    return coords;
}

} // namespace solver
