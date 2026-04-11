/**
 * @file field_data.cpp
 * @brief 求解器层 - 场数据与后处理模块实现
 * @details 实现有限元求解结果的存储和派生物理量的计算功能，
 *          包括电场强度E=-∇φ、电位移矢量D=εE、静电能量W=0.5φᵀKφ，
 *          以及VTK/CSV格式的结果导出。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
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
    , has_potential_(false)
    , has_electric_field_(false)
    , has_electric_displacement_(false)
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
}

bool FieldData::hasData() const {
    return has_potential_;
}

} // namespace solver
