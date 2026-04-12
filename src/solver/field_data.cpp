/**
 * @file field_data.cpp
 * @brief 求解器层 - 场数据与后处理模块实现
 * @details 实现有限元求解结果的存储和派生物理量的计算功能，
 *          包括电场强度E=-∇φ、电位移矢量D=εE、静电能量W=0.5φᵀKφ、
 *          电流密度J=σE、焦耳热功率P=φᵀKφ、总电流I、电压降ΔV等，
 *          以及VTK/CSV格式的结果导出。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.1
 */

#include "field_data.hpp"
#include "em_mesh_data.hpp"
#include "shape_function_factory.hpp"
#include "gauss_quadrature.hpp"
#include "mesh_query.hpp"
#include "element_geometry.hpp"
#include "math_constants.hpp"
#include "logger_factory.hpp"
#include "csr_matrix.hpp"
#include "em_element_integrator_base.hpp"

#include <fstream>
#include <map>
#include <cmath>
#include <unordered_set>

namespace solver {

// ==================== 内部辅助函数 ====================

/**
 * @brief 将ElementType枚举转换为ShapeFunctionFactory所需的字符串标识符
 * @param type 单元类型枚举值
 * @return std::string 对应的类型名称字符串，不支持时返回空字符串
 */
static std::string elementTypeToString(numeric::ElementType type) {
    switch (type) {
        case numeric::ElementType::LINE2:    return "LINE2";
        case numeric::ElementType::LINE3:    return "LINE3";
        case numeric::ElementType::TRI3:     return "TRI3";
        case numeric::ElementType::TRI6:     return "TRI6";
        case numeric::ElementType::QUAD4:    return "QUAD4";
        case numeric::ElementType::QUAD8:    return "QUAD8";
        case numeric::ElementType::QUAD9:    return "QUAD9";
        case numeric::ElementType::TET4:     return "TET4";
        case numeric::ElementType::TET10:    return "TET10";
        case numeric::ElementType::HEX8:     return "HEX8";
        case numeric::ElementType::HEX20:    return "HEX20";
        case numeric::ElementType::HEX27:    return "HEX27";
        case numeric::ElementType::PRISM6:   return "PRISM6";
        case numeric::ElementType::PRISM15:  return "PRISM15";
        case numeric::ElementType::PYRAMID5: return "PYRAMID5";
        case numeric::ElementType::PYRAMID13:return "PYRAMID13";
        default: return "";
    }
}

/**
 * @brief 获取单元的空间维度
 * @param type 单元类型枚举值
 * @return int 空间维度（1/2/3）
 */
static int getElementDimension(numeric::ElementType type) {
    switch (type) {
        case numeric::ElementType::LINE2:
        case numeric::ElementType::LINE3:
            return 1;
        case numeric::ElementType::TRI3:
        case numeric::ElementType::TRI6:
        case numeric::ElementType::QUAD4:
        case numeric::ElementType::QUAD8:
        case numeric::ElementType::QUAD9:
            return 2;
        default:
            return 3;
    }
}

/**
 * @brief 获取单元类型的默认高斯积分阶数
 * @param type 单元类型枚举值
 * @return int 默认积分阶数
 */
static int getDefaultGaussOrder(numeric::ElementType type) {
    switch (type) {
        case numeric::ElementType::TRI3:     return 1;
        case numeric::ElementType::TRI6:     return 3;
        case numeric::ElementType::QUAD4:    return 4;
        case numeric::ElementType::QUAD8:    return 4;
        case numeric::ElementType::QUAD9:    return 9;
        case numeric::ElementType::TET4:     return 1;
        case numeric::ElementType::TET10:    return 4;
        case numeric::ElementType::HEX8:     return 8;
        case numeric::ElementType::HEX20:    return 8;
        case numeric::ElementType::HEX27:    return 8;
        case numeric::ElementType::PRISM6:   return 6;
        case numeric::ElementType::PRISM15:  return 6;
        case numeric::ElementType::PYRAMID5: return 5;
        case numeric::ElementType::PYRAMID13:return 5;
        default:                    return 1;
    }
}

/**
 * @brief 获取单元类型对应的VTK单元类型编号
 * @param type 单元类型枚举值
 * @return int VTK单元类型编号，不支持时返回0
 */
static int getVtkCellType(numeric::ElementType type) {
    switch (type) {
        case numeric::ElementType::TRI3:     return 5;
        case numeric::ElementType::TRI6:     return 22;
        case numeric::ElementType::QUAD4:    return 9;
        case numeric::ElementType::QUAD8:    return 23;
        case numeric::ElementType::QUAD9:    return 28;
        case numeric::ElementType::TET4:     return 10;
        case numeric::ElementType::TET10:    return 24;
        case numeric::ElementType::HEX8:     return 12;
        case numeric::ElementType::HEX20:    return 25;
        case numeric::ElementType::HEX27:    return 29;
        case numeric::ElementType::PRISM6:   return 13;
        case numeric::ElementType::PRISM15:  return 26;
        case numeric::ElementType::PYRAMID5: return 14;
        case numeric::ElementType::PYRAMID13:return 27;
        default:                    return 0;
    }
}

/**
 * @brief 判断单元类型是否为Lagrange节点元（非Nedelec棱边元）
 * @param type 单元类型枚举值
 * @return bool 是Lagrange元返回true
 */
static bool isLagrangeElement(numeric::ElementType type) {
    switch (type) {
        case numeric::ElementType::LINE2:
        case numeric::ElementType::LINE3:
        case numeric::ElementType::TRI3:
        case numeric::ElementType::TRI6:
        case numeric::ElementType::QUAD4:
        case numeric::ElementType::QUAD8:
        case numeric::ElementType::QUAD9:
        case numeric::ElementType::TET4:
        case numeric::ElementType::TET10:
        case numeric::ElementType::HEX8:
        case numeric::ElementType::HEX20:
        case numeric::ElementType::HEX27:
        case numeric::ElementType::PRISM6:
        case numeric::ElementType::PRISM15:
        case numeric::ElementType::PYRAMID5:
        case numeric::ElementType::PYRAMID13:
            return true;
        default:
            return false;
    }
}

// ==================== 构造函数 ====================

FieldData::FieldData()
    : electrostatic_energy_(0.0)
    , total_current_(0.0)
    , joule_heating_total_(0.0)
    , has_potential_(false)
    , has_electric_field_(false)
    , has_electric_displacement_(false)
    , has_current_density_(false)
    , has_joule_heating_(false)
{
}

// ==================== 节点电位 ====================

void FieldData::setNodalPotential(const Eigen::VectorXd& phi, int num_total_dofs) {
    if (phi.size() != num_total_dofs) {
        FEEM_ERROR("setNodalPotential: 电位向量大小({})与总DOF数({})不匹配",
                   phi.size(), num_total_dofs);
        return;
    }
    nodal_potential_ = phi;
    has_potential_ = true;
    FEEM_INFO("setNodalPotential: 节点电位已设置, DOF数={}", num_total_dofs);
}

const Eigen::VectorXd& FieldData::getNodalPotential() const {
    return nodal_potential_;
}

// ==================== 电场强度计算 ====================

void FieldData::computeElectricField(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials)
{
    (void)materials;

    if (!has_potential_) {
        FEEM_ERROR("computeElectricField: 节点电位数据未设置，无法计算电场强度");
        return;
    }

    size_t num_elements = mesh_data.getElementCount();
    elem_electric_field_.assign(num_elements, Eigen::Vector3d::Zero());

    for (size_t e = 0; e < num_elements; ++e) {
        const auto& elem = mesh_data.elements[e];

        if (!isLagrangeElement(elem.elem_type)) {
            FEEM_WARN("computeElectricField: 跳过非Lagrange单元, element_id={}, type={}",
                      elem.id, static_cast<int>(elem.elem_type));
            continue;
        }

        auto elem_nodes = fe_em::MeshQuery::get_element_nodes(elem, mesh_data.nodes);
        if (elem_nodes.empty()) {
            FEEM_WARN("computeElectricField: 单元节点获取失败, element_id={}", elem.id);
            continue;
        }

        int node_count = static_cast<int>(elem_nodes.size());
        int dim = getElementDimension(elem.elem_type);

        Eigen::MatrixXd node_coords(dim, node_count);
        for (int i = 0; i < node_count; ++i) {
            node_coords(0, i) = elem_nodes[i].x;
            if (dim >= 2) node_coords(1, i) = elem_nodes[i].y;
            if (dim >= 3) node_coords(2, i) = elem_nodes[i].z;
        }

        std::string type_str = elementTypeToString(elem.elem_type);
        if (type_str.empty()) {
            FEEM_WARN("computeElectricField: 不支持的单元类型, element_id={}", elem.id);
            continue;
        }

        auto shape_func = numeric::ShapeFunctionFactory::create(type_str);
        if (!shape_func) {
            FEEM_WARN("computeElectricField: 形函数创建失败, element_id={}, type={}",
                      elem.id, type_str);
            continue;
        }

        int gauss_order = getDefaultGaussOrder(elem.elem_type);
        auto gauss_points = numeric::GaussQuadrature::getPoints(elem.elem_type, gauss_order);
        if (gauss_points.empty()) {
            FEEM_WARN("computeElectricField: 高斯积分点获取失败, element_id={}", elem.id);
            continue;
        }

        Eigen::VectorXd phi_elem(node_count);
        for (int i = 0; i < node_count; ++i) {
            int node_id = elem_nodes[i].id;
            if (node_id < 0 || node_id >= nodal_potential_.size()) {
                FEEM_ERROR("computeElectricField: 节点ID({})超出电位向量范围({}), element_id={}",
                           node_id, nodal_potential_.size(), elem.id);
                phi_elem(i) = 0.0;
            } else {
                phi_elem(i) = nodal_potential_(node_id);
            }
        }

        Eigen::Vector3d E_sum = Eigen::Vector3d::Zero();
        for (const auto& gp : gauss_points) {
            numeric::LocalPoint xi(gp.coords[0], gp.coords[1], gp.coords[2]);
            xi.dim = gp.dim;

            Eigen::MatrixXd grad_N = shape_func->calcPhysicalGradN(xi, node_coords);

            Eigen::VectorXd E_local = -grad_N.transpose() * phi_elem;

            E_sum[0] += E_local[0];
            if (E_local.size() > 1) E_sum[1] += E_local[1];
            if (E_local.size() > 2) E_sum[2] += E_local[2];
        }

        elem_electric_field_[e] = E_sum / static_cast<double>(gauss_points.size());
    }

    has_electric_field_ = true;
    FEEM_INFO("computeElectricField: 电场强度计算完成, 单元数={}", num_elements);
}

const std::vector<Eigen::Vector3d>& FieldData::getElementElectricField() const {
    return elem_electric_field_;
}

// ==================== 电位移矢量计算 ====================

void FieldData::computeElectricDisplacement(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials)
{
    if (!has_electric_field_) {
        FEEM_INFO("computeElectricDisplacement: 电场未计算，自动调用computeElectricField");
        computeElectricField(mesh_data, materials);
    }

    if (!has_electric_field_) {
        FEEM_ERROR("computeElectricDisplacement: 电场强度计算失败，无法计算电位移矢量");
        return;
    }

    size_t num_elements = mesh_data.getElementCount();
    elem_electric_displacement_.assign(num_elements, Eigen::Vector3d::Zero());

    for (size_t e = 0; e < num_elements; ++e) {
        const auto& elem = mesh_data.elements[e];

        double epsilon = numeric::EPSILON0;
        auto mat_it = materials.find(elem.material_id);
        if (mat_it != materials.end()) {
            epsilon = mat_it->second.epsilon;
        } else {
            FEEM_DEBUG("computeElectricDisplacement: 材料ID={}未找到，使用真空介电常数, element_id={}",
                       elem.material_id, elem.id);
        }

        elem_electric_displacement_[e] = epsilon * elem_electric_field_[e];
    }

    has_electric_displacement_ = true;
    FEEM_INFO("computeElectricDisplacement: 电位移矢量计算完成, 单元数={}", num_elements);
}

const std::vector<Eigen::Vector3d>& FieldData::getElementElectricDisplacement() const {
    return elem_electric_displacement_;
}

// ==================== 静电能量计算 ====================

double FieldData::computeElectrostaticEnergy(
    const numeric::CsrMatrix<double>& K,
    const Eigen::VectorXd& phi) const
{
    if (!K.is_built()) {
        FEEM_ERROR("computeElectrostaticEnergy: 刚度矩阵未构建");
        return 0.0;
    }

    int n = static_cast<int>(phi.size());

    std::vector<double> phi_vec(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        phi_vec[static_cast<size_t>(i)] = phi(i);
    }

    std::vector<double> k_phi_vec(static_cast<size_t>(n), 0.0);
    K.mat_vec(phi_vec, k_phi_vec);

    double energy = 0.0;
    for (int i = 0; i < n; ++i) {
        energy += phi(i) * k_phi_vec[static_cast<size_t>(i)];
    }
    energy *= 0.5;

    FEEM_INFO("computeElectrostaticEnergy: 静电能量 = {:.6e} J", energy);
    return energy;
}

// ==================== 电流密度计算 ====================

void FieldData::computeCurrentDensity(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials)
{
    if (!has_electric_field_) {
        FEEM_INFO("computeCurrentDensity: 电场未计算，自动调用computeElectricField");
        computeElectricField(mesh_data, materials);
    }

    if (!has_electric_field_) {
        FEEM_ERROR("computeCurrentDensity: 电场强度计算失败，无法计算电流密度");
        return;
    }

    size_t num_elements = mesh_data.getElementCount();
    elem_current_density_.assign(num_elements, Eigen::Vector3d::Zero());

    for (size_t e = 0; e < num_elements; ++e) {
        const auto& elem = mesh_data.elements[e];

        double sigma = 0.0;
        auto mat_it = materials.find(elem.material_id);
        if (mat_it != materials.end()) {
            sigma = mat_it->second.sigma;
        } else {
            FEEM_DEBUG("computeCurrentDensity: 材料ID={}未找到，使用默认电导率0, element_id={}",
                       elem.material_id, elem.id);
        }

        elem_current_density_[e] = sigma * elem_electric_field_[e];
    }

    has_current_density_ = true;
    FEEM_INFO("computeCurrentDensity: 电流密度计算完成, 单元数={}", num_elements);
}

const std::vector<Eigen::Vector3d>& FieldData::getElementCurrentDensity() const {
    return elem_current_density_;
}

// ==================== 焦耳热计算 ====================

double FieldData::computeJouleHeating(
    const numeric::CsrMatrix<double>& K,
    const Eigen::VectorXd& phi) const
{
    if (!K.is_built()) {
        FEEM_ERROR("computeJouleHeating: 刚度矩阵未构建");
        return 0.0;
    }

    int n = static_cast<int>(phi.size());

    std::vector<double> phi_vec(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        phi_vec[static_cast<size_t>(i)] = phi(i);
    }

    std::vector<double> k_phi_vec(static_cast<size_t>(n), 0.0);
    K.mat_vec(phi_vec, k_phi_vec);

    double power = 0.0;
    for (int i = 0; i < n; ++i) {
        power += phi(i) * k_phi_vec[static_cast<size_t>(i)];
    }

    FEEM_INFO("computeJouleHeating: 总焦耳热功率 = {:.6e} W", power);
    return power;
}

void FieldData::computeElementPowerDensity() {
    if (!has_current_density_ || !has_electric_field_) {
        FEEM_ERROR("computeElementPowerDensity: 需要先计算电流密度和电场强度");
        return;
    }

    size_t num_elements = elem_electric_field_.size();
    elem_joule_power_density_.resize(num_elements);

    for (size_t e = 0; e < num_elements; ++e) {
        double p_e = elem_current_density_[e].dot(elem_electric_field_[e]);
        elem_joule_power_density_[e] = p_e;
    }

    has_joule_heating_ = true;
    FEEM_INFO("computeElementPowerDensity: 单元焦耳热功率密度计算完成, 单元数={}", num_elements);
}

const std::vector<double>& FieldData::getElementPowerDensity() const {
    return elem_joule_power_density_;
}

// ==================== 总电流计算辅助函数 ====================

/**
 * @brief 判断是否为二维单元类型
 * @param elem_type 单元类型枚举值
 * @return true 如果是二维单元（TRI3/TRI6/QUAD4/QUAD8/QUAD9）
 */
static bool is2DElement(fe_em::ElemType elem_type) {
    switch (elem_type) {
        case fe_em::ElemType::TRI3:
        case fe_em::ElemType::TRI6:
        case fe_em::ElemType::QUAD4:
        case fe_em::ElemType::QUAD8:
        case fe_em::ElemType::QUAD9:
            return true;
        default:
            return false;
    }
}

/**
 * @brief 判断是否为三维单元类型
 * @param elem_type 单元类型枚举值
 * @return true 如果是三维单元（TET4/TET10/HEX8/HEX20/HEX27/PRISM6/PRISM15/PYRAMID5/PYRAMID13）
 */
static bool is3DElement(fe_em::ElemType elem_type) {
    switch (elem_type) {
        case fe_em::ElemType::TET4:
        case fe_em::ElemType::TET10:
        case fe_em::ElemType::HEX8:
        case fe_em::ElemType::HEX20:
        case fe_em::ElemType::HEX27:
        case fe_em::ElemType::PRISM6:
        case fe_em::ElemType::PRISM15:
        case fe_em::ElemType::PYRAMID5:
        case fe_em::ElemType::PYRAMID13:
            return true;
        default:
            return false;
    }
}

/**
 * @brief 获取二维单元的局部边定义（节点局部索引列表）
 * @param elem_type 二维单元类型枚举值
 * @return std::vector<std::vector<int>> 每条边的局部节点索引列表
 */
static std::vector<std::vector<int>> get2DEdgeDefinitions(fe_em::ElemType elem_type) {
    switch (elem_type) {
        case fe_em::ElemType::TRI3:
        case fe_em::ElemType::TRI6:
            return {{0, 1}, {1, 2}, {2, 0}};
        case fe_em::ElemType::QUAD4:
        case fe_em::ElemType::QUAD8:
        case fe_em::ElemType::QUAD9:
            return {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        default:
            return {};
    }
}

/**
 * @brief 根据面的节点数返回角点数量（处理退化情况）
 * @param face_node_count 面定义中的节点数
 * @return int 实际角点数量
 */
static int getCornerNodeCount(int face_node_count) {
    if (face_node_count >= 4) return 4;      /* 四边形面 */
    if (face_node_count == 3) return 3;       /* 三角形面 */
    if (face_node_count == 2) return 2;       /* 边 */
    return face_node_count;
}

/**
 * @brief O(1) 获取节点的物理坐标（通过节点ID直接索引）
 * @param node_id 节点ID
 * @param mesh_data 网格数据
 * @return Eigen::Vector3d 节点坐标 (x, y, z)，未找到返回零向量
 */
static Eigen::Vector3d getNodeCoords(int node_id, const fe_em::EMMeshData& mesh_data) {
    /* 节点ID通常等于其在nodes向量中的索引（0-based），优先直接索引 */
    if (node_id >= 0 && node_id < static_cast<int>(mesh_data.nodes.size())) {
        const auto& nd = mesh_data.nodes[node_id];
        if (nd.id == node_id) {
            return Eigen::Vector3d(nd.x, nd.y, nd.z);
        }
    }
    /* 回退：线性查找（应对ID不连续的情况） */
    for (const auto& nd : mesh_data.nodes) {
        if (nd.id == node_id) {
            return Eigen::Vector3d(nd.x, nd.y, nd.z);
        }
    }
    return Eigen::Vector3d::Zero();
}

/**
 * @brief 计算单元的质心坐标
 * @param elem 单元数据
 * @param mesh_data 网格拓扑数据
 * @return Eigen::Vector3d 质心坐标
 */
static Eigen::Vector3d computeElementCentroid(
    const fe_em::Element& elem,
    const fe_em::EMMeshData& mesh_data) {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    int count = 0;
    for (int node_id : elem.node_ids) {
        centroid += getNodeCoords(node_id, mesh_data);
        ++count;
    }
    if (count > 0) {
        centroid /= static_cast<double>(count);
    }
    return centroid;
}

/**
 * @brief 计算二维边的外法向量和长度（确保指向单元外部）
 * @param p0 边的第一个角节点坐标
 * @param p1 边的第二个角节点坐标
 * @param centroid 所属单元的质心坐标
 * @param normal [out] 边的外法向量（从单元内部指向外部）
 * @param length [out] 边的长度
 */
static void computeEdgeNormalAndLength(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& centroid,
    Eigen::Vector3d& normal,
    double& length) {
    Eigen::Vector3d edge = p1 - p0;
    length = std::sqrt(edge[0] * edge[0] + edge[1] * edge[1]);

    if (length > 1e-15) {
        normal[0] = edge[1] / length;
        normal[1] = -edge[0] / length;
        normal[2] = 0.0;

        Eigen::Vector3d edge_mid = 0.5 * (p0 + p1);
        Eigen::Vector3d to_outside = edge_mid - centroid;
        if (normal.dot(to_outside) < 0.0) {
            normal = -normal;
        }
    } else {
        normal.setZero();
        length = 0.0;
    }
}

/**
 * @brief 计算三维三角形面的外法向量和面积（确保指向单元外部）
 * @param p0/p1/p2 三个角节点坐标
 * @param centroid 所属单元质心
 * @param normal [out] 外法向量
 * @param area [out] 面积
 */
static void computeTriangleNormalAndArea(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& p2,
    const Eigen::Vector3d& centroid,
    Eigen::Vector3d& normal,
    double& area) {
    Eigen::Vector3d v1 = p1 - p0;
    Eigen::Vector3d v2 = p2 - p0;
    normal = v1.cross(v2);
    double cross_norm = normal.norm();
    area = 0.5 * cross_norm;

    if (cross_norm > 1e-15) {
        normal /= cross_norm;

        Eigen::Vector3d face_centroid = (p0 + p1 + p2) / 3.0;
        Eigen::Vector3d to_outside = face_centroid - centroid;
        if (normal.dot(to_outside) < 0.0) {
            normal = -normal;
        }
    } else {
        normal.setZero();
        area = 0.0;
    }
}

/**
 * @brief 计算三维四边形面的外法向量和面积（确保指向单元外部）
 * @param p0/p1/p2/p3 四个角节点坐标
 * @param centroid 所属单元质心
 * @param normal [out] 外法向量
 * @param area [out] 面积
 */
static void computeQuadNormalAndArea(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& p2,
    const Eigen::Vector3d& p3,
    const Eigen::Vector3d& centroid,
    Eigen::Vector3d& normal,
    double& area) {
    /* 将四边形拆分为两个三角形分别计算后平均 */
    Eigen::Vector3d n1, n2;
    double a1, a2;
    computeTriangleNormalAndArea(p0, p1, p2, centroid, n1, a1);
    computeTriangleNormalAndArea(p0, p2, p3, centroid, n2, a2);

    double total_area = a1 + a2;
    if (total_area > 1e-15) {
        normal = (n1 * a1 + n2 * a2) / total_area;
        normal.normalize();
        area = total_area;
    } else {
        normal.setZero();
        area = 0.0;
    }
}

// ==================== 总电流计算 ====================

double FieldData::computeTotalCurrent(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<Eigen::Vector3d>& elem_J) const
{
    if (elem_J.empty()) {
        FEEM_WARN("computeTotalCurrent: 电流密度向量为空，返回0");
        return 0.0;
    }

    size_t num_elements = mesh_data.getElementCount();
    if (elem_J.size() != num_elements) {
        FEEM_ERROR("computeTotalCurrent: 电流密度向量大小({})与单元数({})不匹配",
                   elem_J.size(), num_elements);
        return 0.0;
    }

    /* ========== 边界面积分法：I = ∫_Γ J·n dS ==========
     *
     * 遍历所有单元的面/边，找出与Dirichlet端子共享边界的面，
     * 对每个外边界执行 I_face = J·n × Area 并累加。
     *
     * 实现参考 ResistanceCalculator::computeTotalCurrent 的成熟模式：
     *   1. 单元循环 → 局部面定义 → 端子归属检查
     *   2. O(1) 节点坐标查找（getNodeCoords 直接索引）
     *   3. 成熟法向/面积计算函数（computeEdgeNormalAndLength 等）
     *   4. 多端子遍历：统计每个Dirichlet端子的电流，返回绝对值最大者 */

    /* 提取所有 Dirichlet 端子的节点集合 */
    std::vector<std::unordered_set<int>> all_terminals;
    std::vector<double> terminal_currents;

    for (const auto& marker : mesh_data.boundary_markers) {
        if (marker.bnd_type != tool::BndType::DIRICHLET) continue;
        if (!std::holds_alternative<std::vector<int>>(marker.target_ids)) continue;
        const auto& node_ids = std::get<std::vector<int>>(marker.target_ids);
        if (node_ids.empty()) continue;

        std::unordered_set<int> terminal_nodes(node_ids.begin(), node_ids.end());
        all_terminals.push_back(std::move(terminal_nodes));
    }

    if (all_terminals.empty()) {
        FEEM_WARN("computeTotalCurrent: 未找到 Dirichlet 端子，返回0");
        return 0.0;
    }

    /* 对每个端子分别计算电流 */
    terminal_currents.resize(all_terminals.size(), 0.0);

    for (size_t t_idx = 0; t_idx < all_terminals.size(); ++t_idx) {
        const auto& terminal_nodes = all_terminals[t_idx];
        int contributing_faces = 0;

    /* 遍历所有单元，检查每个面是否属于端子 */
    for (size_t e = 0; e < num_elements; ++e) {
        const auto& elem = mesh_data.elements[e];
        if (elem.dof_type != fe_em::DOFType::SCALAR_ONLY) continue;
        if (e >= elem_J.size()) continue;

        Eigen::Vector3d J_e = elem_J[e];

        /* 获取单元局部面/边定义 */
        std::vector<std::vector<int>> local_faces;
        if (is2DElement(elem.elem_type)) {
            local_faces = get2DEdgeDefinitions(elem.elem_type);
        } else if (is3DElement(elem.elem_type)) {
            try {
                local_faces = fe_em::ElementGeometry::get_local_faces(elem.elem_type);
            } catch (const std::invalid_argument&) {
                continue;
            }
        } else {
            continue;
        }

        /* 对每个局部面检查是否完全属于目标端子 */
        for (int f_idx = 0; f_idx < static_cast<int>(local_faces.size()); ++f_idx) {
            const auto& face_def = local_faces[f_idx];

            /* 检查面的所有角节点是否都属于该端子 */
            bool all_on_terminal = true;
            for (int ln : face_def) {
                if (ln < 0 || ln >= static_cast<int>(elem.node_ids.size())) {
                    all_on_terminal = false;
                    break;
                }
                int nid = elem.node_ids[ln];
                if (terminal_nodes.find(nid) == terminal_nodes.end()) {
                    all_on_terminal = false;
                    break;
                }
            }
            if (!all_on_terminal) continue;

            /* 获取面角节点坐标（O(1) 直接索引） */
            int corner_count = getCornerNodeCount(static_cast<int>(face_def.size()));
            std::vector<Eigen::Vector3d> corner_coords;
            corner_coords.reserve(corner_count);
            for (int k = 0; k < corner_count && k < static_cast<int>(face_def.size()); ++k) {
                int local_node = face_def[k];
                if (local_node >= 0 && local_node < static_cast<int>(elem.node_ids.size())) {
                    corner_coords.push_back(
                        getNodeCoords(elem.node_ids[local_node], mesh_data));
                }
            }

            if (corner_coords.size() < 2) continue;

            /* 计算外法向和面积（使用成熟辅助函数） */
            Eigen::Vector3d normal = Eigen::Vector3d::Zero();
            double area = 0.0;
            Eigen::Vector3d centroid = computeElementCentroid(elem, mesh_data);

            if (corner_coords.size() == 2) {
                computeEdgeNormalAndLength(
                    corner_coords[0], corner_coords[1], centroid, normal, area);
            } else if (corner_coords.size() == 3) {
                computeTriangleNormalAndArea(
                    corner_coords[0], corner_coords[1], corner_coords[2],
                    centroid, normal, area);
            } else if (corner_coords.size() >= 4) {
                computeQuadNormalAndArea(
                    corner_coords[0], corner_coords[1],
                    corner_coords[2], corner_coords[3],
                    centroid, normal, area);
            }

            /* 电流贡献: I_face = J·n × Area */
            double current_contribution = J_e.dot(normal) * area;
            terminal_currents[t_idx] += current_contribution;
            ++contributing_faces;
        }
    }  /* end for (e): 单元循环 */
    }  /* end for (t_idx): 端子循环 */

    /* 选择绝对值最大的端子电流作为总电流（避免正负端子相互抵消） */
    double total_current = 0.0;
    int best_terminal = -1;
    for (size_t t = 0; t < terminal_currents.size(); ++t) {
        if (std::abs(terminal_currents[t]) > std::abs(total_current)) {
            total_current = terminal_currents[t];
            best_terminal = static_cast<int>(t);
        }
    }

    FEEM_INFO("computeTotalCurrent: 总电流（边界面积分）= {:.6e} A (端子#{}, 共{}个端子)",
              total_current, best_terminal + 1, static_cast<int>(all_terminals.size()));
    return total_current;
}

// ==================== 电压降计算 ====================

double FieldData::computeVoltageDrop(int node_id_1, int node_id_2) const {
    if (!has_potential_) {
        FEEM_ERROR("computeVoltageDrop: 节点电位数据未设置");
        return 0.0;
    }

    if (node_id_1 < 0 || node_id_1 >= nodal_potential_.size()) {
        FEEM_ERROR("computeVoltageDrop: 节点ID({})超出电位向量范围({})",
                   node_id_1, nodal_potential_.size());
        return 0.0;
    }

    if (node_id_2 < 0 || node_id_2 >= nodal_potential_.size()) {
        FEEM_ERROR("computeVoltageDrop: 节点ID({})超出电位向量范围({})",
                   node_id_2, nodal_potential_.size());
        return 0.0;
    }

    double delta_V = nodal_potential_(node_id_1) - nodal_potential_(node_id_2);
    FEEM_DEBUG("computeVoltageDrop: V({}) - V({}) = {:.6e} V", node_id_1, node_id_2, delta_V);
    return delta_V;
}

// ==================== 结果导出 ====================

bool FieldData::exportVTK(const std::string& filename,
                          const fe_em::EMMeshData& mesh_data) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        FEEM_ERROR("exportVTK: 无法打开输出文件: {}", filename);
        return false;
    }

    std::map<int, int> node_id_to_vtk_index;
    for (size_t i = 0; i < mesh_data.nodes.size(); ++i) {
        node_id_to_vtk_index[mesh_data.nodes[i].id] = static_cast<int>(i);
    }

    size_t num_nodes = mesh_data.getNodeCount();
    size_t num_elements = mesh_data.getElementCount();

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Electrostatic Field Output\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    ofs << "POINTS " << num_nodes << " double\n";
    ofs.precision(15);
    for (size_t i = 0; i < num_nodes; ++i) {
        const auto& node = mesh_data.nodes[i];
        ofs << node.x << " " << node.y << " " << node.z << "\n";
    }

    size_t cell_list_size = 0;
    for (size_t e = 0; e < num_elements; ++e) {
        cell_list_size += 1 + mesh_data.elements[e].node_ids.size();
    }

    ofs << "CELLS " << num_elements << " " << cell_list_size << "\n";
    for (size_t e = 0; e < num_elements; ++e) {
        const auto& elem = mesh_data.elements[e];
        ofs << elem.node_ids.size();
        for (int node_id : elem.node_ids) {
            auto it = node_id_to_vtk_index.find(node_id);
            if (it != node_id_to_vtk_index.end()) {
                ofs << " " << it->second;
            } else {
                ofs << " 0";
                FEEM_WARN("exportVTK: 节点ID={}未找到映射, element_id={}", node_id, elem.id);
            }
        }
        ofs << "\n";
    }

    ofs << "CELL_TYPES " << num_elements << "\n";
    for (size_t e = 0; e < num_elements; ++e) {
        ofs << getVtkCellType(mesh_data.elements[e].elem_type) << "\n";
    }

    if (has_potential_) {
        ofs << "POINT_DATA " << num_nodes << "\n";
        ofs << "SCALARS Potential double 1\n";
        ofs << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < num_nodes; ++i) {
            int node_id = mesh_data.nodes[i].id;
            if (node_id >= 0 && node_id < nodal_potential_.size()) {
                ofs << nodal_potential_(node_id) << "\n";
            } else {
                ofs << "0.0\n";
            }
        }
    }

    if (has_electric_field_ && !elem_electric_field_.empty()) {
        ofs << "CELL_DATA " << num_elements << "\n";
        ofs << "VECTORS ElectricField double\n";
        for (size_t e = 0; e < num_elements; ++e) {
            const auto& E = elem_electric_field_[e];
            ofs << E[0] << " " << E[1] << " " << E[2] << "\n";
        }
    }

    ofs.close();
    FEEM_INFO("exportVTK: VTK文件导出成功: {}", filename);
    return true;
}

bool FieldData::exportCSV(const std::string& filename,
                          const fe_em::EMMeshData& mesh_data) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        FEEM_ERROR("exportCSV: 无法打开输出文件: {}", filename);
        return false;
    }

    ofs << "node_id,x,y,z,potential\n";
    ofs.precision(15);

    for (const auto& node : mesh_data.nodes) {
        ofs << node.id << "," << node.x << "," << node.y << "," << node.z << ",";
        if (has_potential_ && node.id >= 0 && node.id < nodal_potential_.size()) {
            ofs << nodal_potential_(node.id);
        } else {
            ofs << "0.0";
        }
        ofs << "\n";
    }

    ofs.close();
    FEEM_INFO("exportCSV: CSV文件导出成功: {}", filename);
    return true;
}

// ==================== 资源管理 ====================

void FieldData::clear() {
    nodal_potential_.resize(0);
    elem_electric_field_.clear();
    elem_electric_displacement_.clear();
    electrostatic_energy_ = 0.0;
    has_potential_ = false;
    has_electric_field_ = false;
    has_electric_displacement_ = false;

    elem_current_density_.clear();
    elem_joule_power_density_.clear();
    total_current_ = 0.0;
    joule_heating_total_ = 0.0;
    has_current_density_ = false;
    has_joule_heating_ = false;
}

bool FieldData::hasData() const {
    return has_potential_;
}

} // namespace solver
