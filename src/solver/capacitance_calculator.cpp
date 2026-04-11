/**
 * @file capacitance_calculator.cpp
 * @brief 求解器层 - 电容矩阵计算器实现文件
 * @details 基于电荷法实现多导体系统Maxwell电容矩阵的计算。
 *          核心流程：对每个导体依次施加单位电压求解静电场，
 *          通过高斯面积分提取各导体表面电荷，构建电容矩阵。
 *
 *          电荷法原理：
 *          Q_i = Σ_j C_ij * V_j
 *          当V_k = 1V, V_j = 0V (j≠k)时：
 *          Q_i = C_ik，即C_ik = Q_i（V_k = 1V时导体i上的电荷）
 *
 *          表面电荷计算：
 *          Q = -∫_Γ D·n_elem dΓ
 *          其中n_elem为单元面外法向（指向导体），负号转换为导体外法向。
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#include "capacitance_calculator.hpp"
#include "electrostatic_solver.hpp"
#include "field_data.hpp"
#include "logger_factory.hpp"
#include "em_enums.hpp"
#include "mesh_query.hpp"
#include "element_geometry.hpp"
#include "math_constants.hpp"

#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <unordered_set>

namespace solver {

/* ==================== 内部辅助结构体 ==================== */

/**
 * @brief 面信息结构体，用于边界面到单元的映射
 */
struct FaceInfo {
    int element_index = -1;     ///< 相邻单元在mesh_data.elements中的索引
    int face_local_index = -1;  ///< 面在单元局部定义中的索引
};

/**
 * @brief 面键类型，使用排序后的面角节点ID列表作为唯一标识
 */
using FaceKey = std::vector<int>;

/**
 * @brief 边界面到单元映射类型
 */
using FaceToElementMap = std::map<FaceKey, FaceInfo>;

/* ==================== 内部辅助函数 ==================== */

/**
 * @brief 判断单元类型是否为二维单元
 * @param elem_type 单元类型枚举
 * @return bool 二维单元返回true
 */
static bool is2DElement(fe_em::ElemType elem_type) {
    return elem_type == numeric::ElementType::TRI3
        || elem_type == numeric::ElementType::TRI6
        || elem_type == numeric::ElementType::QUAD4
        || elem_type == numeric::ElementType::QUAD8
        || elem_type == numeric::ElementType::QUAD9;
}

/**
 * @brief 判断单元类型是否为三维单元
 * @param elem_type 单元类型枚举
 * @return bool 三维单元返回true
 */
static bool is3DElement(fe_em::ElemType elem_type) {
    return elem_type == numeric::ElementType::TET4
        || elem_type == numeric::ElementType::TET10
        || elem_type == numeric::ElementType::HEX8
        || elem_type == numeric::ElementType::HEX20
        || elem_type == numeric::ElementType::HEX27
        || elem_type == numeric::ElementType::PRISM6
        || elem_type == numeric::ElementType::PRISM15
        || elem_type == numeric::ElementType::PYRAMID5
        || elem_type == numeric::ElementType::PYRAMID13;
}

/**
 * @brief 获取二维单元的边定义（保留方向信息用于法向计算）
 * @param elem_type 单元类型枚举
 * @return std::vector<std::vector<int>> 边定义列表，每条边为有序局部节点ID
 *
 * @note 边的节点按逆时针单元的外法向约定排列，
 *       即对于边(n0,n1)，外法向为(dy,-dx)方向（右手定则）
 */
static std::vector<std::vector<int>> get2DEdgeDefinitions(fe_em::ElemType elem_type) {
    switch (elem_type) {
        case numeric::ElementType::TRI3:
            return {{0, 1}, {1, 2}, {2, 0}};
        case numeric::ElementType::TRI6:
            return {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
        case numeric::ElementType::QUAD4:
            return {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        case numeric::ElementType::QUAD8:
            return {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
        case numeric::ElementType::QUAD9:
            return {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
        default:
            return {};
    }
}

/**
 * @brief 从面定义中提取角节点数量
 * @param face_node_count 面的总节点数
 * @return int 角节点数量（2/3/4）
 */
static int getCornerNodeCount(int face_node_count) {
    if (face_node_count <= 2) return 2;
    if (face_node_count <= 3) return 3;
    if (face_node_count <= 6) return 3;
    return 4;
}

/**
 * @brief 从网格数据中获取节点坐标
 * @param node_id 节点ID（作为nodes向量的索引）
 * @param mesh_data 网格拓扑数据
 * @return Eigen::Vector3d 节点坐标
 */
static Eigen::Vector3d getNodeCoords(int node_id, const fe_em::EMMeshData& mesh_data) {
    if (node_id < 0 || node_id >= static_cast<int>(mesh_data.nodes.size())) {
        return Eigen::Vector3d::Zero();
    }
    const auto& node = mesh_data.nodes[node_id];
    return Eigen::Vector3d(node.x, node.y, node.z);
}

/**
 * @brief 构建边界面到单元的映射
 * @param mesh_data 网格拓扑数据
 * @return FaceToElementMap 边界面映射（仅包含属于一个单元的面）
 *
 * @details 遍历所有单元，提取每个单元的面/边定义，
 *          通过排序后的角节点ID作为键。若同一键出现两次，
 *          则该面为内部面（被两个单元共享），从映射中移除；
 *          最终映射中仅保留边界面（仅属于一个单元的面）。
 */
static FaceToElementMap buildBoundaryFaceMap(const fe_em::EMMeshData& mesh_data) {
    FaceToElementMap all_faces;

    for (size_t e = 0; e < mesh_data.elements.size(); ++e) {
        const auto& elem = mesh_data.elements[e];

        if (elem.dof_type != fe_em::DOFType::SCALAR_ONLY) {
            continue;
        }

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

        for (size_t f = 0; f < local_faces.size(); ++f) {
            const auto& face_def = local_faces[f];
            int corner_count = getCornerNodeCount(static_cast<int>(face_def.size()));

            FaceKey key;
            key.reserve(corner_count);
            for (int k = 0; k < corner_count && k < static_cast<int>(face_def.size()); ++k) {
                int local_node = face_def[k];
                if (local_node >= 0 && local_node < static_cast<int>(elem.node_ids.size())) {
                    key.push_back(elem.node_ids[local_node]);
                }
            }
            std::sort(key.begin(), key.end());

            if (key.empty()) continue;

            auto it = all_faces.find(key);
            if (it != all_faces.end()) {
                all_faces.erase(it);
            } else {
                all_faces[key] = {static_cast<int>(e), static_cast<int>(f)};
            }
        }
    }

    return all_faces;
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
 * @param p0 第一个角节点坐标
 * @param p1 第二个角节点坐标
 * @param p2 第三个角节点坐标
 * @param centroid 所属单元的质心坐标
 * @param normal [out] 面的外法向量
 * @param area [out] 面的面积
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
    Eigen::Vector3d cross = v1.cross(v2);
    double cross_norm = cross.norm();
    area = 0.5 * cross_norm;

    if (cross_norm > 1e-15) {
        normal = cross / cross_norm;

        Eigen::Vector3d face_center = (p0 + p1 + p2) / 3.0;
        Eigen::Vector3d to_outside = face_center - centroid;
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
 * @param p0 第一个角节点坐标
 * @param p1 第二个角节点坐标
 * @param p2 第三个角节点坐标
 * @param p3 第四个角节点坐标
 * @param centroid 所属单元的质心坐标
 * @param normal [out] 面的外法向量
 * @param area [out] 面的面积
 */
static void computeQuadNormalAndArea(
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& p2,
    const Eigen::Vector3d& p3,
    const Eigen::Vector3d& centroid,
    Eigen::Vector3d& normal,
    double& area) {
    Eigen::Vector3d v1 = p1 - p0;
    Eigen::Vector3d v2 = p3 - p0;
    Eigen::Vector3d cross = v1.cross(v2);
    double cross_norm = cross.norm();

    if (cross_norm > 1e-15) {
        normal = cross / cross_norm;

        Eigen::Vector3d face_center = 0.25 * (p0 + p1 + p2 + p3);
        Eigen::Vector3d to_outside = face_center - centroid;
        if (normal.dot(to_outside) < 0.0) {
            normal = -normal;
        }
    } else {
        normal.setZero();
    }

    double area1 = 0.5 * (p1 - p0).cross(p2 - p0).norm();
    double area2 = 0.5 * (p2 - p0).cross(p3 - p0).norm();
    area = area1 + area2;
}

/**
 * @brief 计算导体表面电荷（高斯面积分法）
 * @param field_data 场数据（包含D场）
 * @param mesh_data 网格拓扑数据
 * @param conductor_info 导体信息
 * @param boundary_face_map 边界面映射
 * @param dim_type 维度类型
 * @return double 导体表面电荷（2D: C/m, 3D: C, 轴对称: C）
 *
 * @details 对导体表面所有边界面/边进行积分：
 *          Q = -Σ_f (D_f · n_elem_f) * A_f
 *          其中n_elem_f为单元面外法向（指向导体外部），
 *          负号将法向从单元外侧转换为导体外侧（进入周围介质方向）。
 *
 *          轴对称情况下附加2πr因子：
 *          Q = -Σ_f (D_f · n_elem_f) * A_f * 2π * r_f
 */
static double computeConductorCharge(
    const FieldData& field_data,
    const fe_em::EMMeshData& mesh_data,
    const ConductorInfo& conductor_info,
    const FaceToElementMap& boundary_face_map,
    tool::DimType dim_type) {
    std::unordered_set<int> conductor_nodes(
        conductor_info.node_ids.begin(), conductor_info.node_ids.end());

    const auto& D_field = field_data.getElementElectricDisplacement();

    double total_charge = 0.0;

    for (const auto& [key, face_info] : boundary_face_map) {
        bool on_conductor = true;
        for (int node_id : key) {
            if (conductor_nodes.find(node_id) == conductor_nodes.end()) {
                on_conductor = false;
                break;
            }
        }
        if (!on_conductor) continue;

        int elem_idx = face_info.element_index;
        if (elem_idx < 0 || elem_idx >= static_cast<int>(D_field.size())) {
            continue;
        }

        Eigen::Vector3d D = D_field[elem_idx];

        const auto& elem = mesh_data.elements[elem_idx];

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

        if (face_info.face_local_index < 0
            || face_info.face_local_index >= static_cast<int>(local_faces.size())) {
            continue;
        }

        const auto& face_def = local_faces[face_info.face_local_index];
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

        double D_dot_n = D.dot(normal);
        double charge_contribution = D_dot_n * area;

        if (dim_type == tool::DimType::AXIS) {
            double r_center = 0.0;
            for (const auto& coord : corner_coords) {
                r_center += coord[0];
            }
            r_center /= static_cast<double>(corner_coords.size());
            charge_contribution *= 2.0 * numeric::PI * r_center;
        }

        total_charge -= charge_contribution;
    }

    return total_charge;
}

/* ==================== 构造函数 ==================== */

CapacitanceCalculator::CapacitanceCalculator()
    : conductor_infos_()
    , capacitance_matrix_()
    , original_boundaries_() {
}

/* ==================== 公共方法 ==================== */

Eigen::MatrixXd CapacitanceCalculator::computeCapacitanceMatrix(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    const std::vector<int>& conductor_ids,
    tool::DimType dim_type) {
    clear();

    int N = static_cast<int>(conductor_ids.size());
    if (N == 0) {
        FEEM_ERROR("CapacitanceCalculator::computeCapacitanceMatrix 导体ID列表为空");
        return Eigen::MatrixXd();
    }

    if (mesh_data.elements.empty()) {
        FEEM_ERROR("CapacitanceCalculator::computeCapacitanceMatrix 网格单元为空");
        return Eigen::MatrixXd();
    }

    if (!extractConductorInfos(boundaries, mesh_data, conductor_ids)) {
        FEEM_ERROR("CapacitanceCalculator::computeCapacitanceMatrix 导体信息提取失败");
        return Eigen::MatrixXd();
    }

    original_boundaries_ = boundaries;

    FEEM_INFO("CapacitanceCalculator: 开始构建边界面映射...");
    FaceToElementMap boundary_face_map = buildBoundaryFaceMap(mesh_data);
    FEEM_INFO("CapacitanceCalculator: 边界面映射构建完成, 边界面数量={}",
              boundary_face_map.size());

    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

    for (int i = 0; i < N; ++i) {
        FEEM_INFO("CapacitanceCalculator: 求解第{}/{}个导体施加单位电压 (导体ID={}, 名称={})",
                  i + 1, N, conductor_infos_[i].id, conductor_infos_[i].name);

        std::vector<tool::Boundary> modified_bcs = buildBoundaryConditions(i, N);

        ElectrostaticSolver solver(dim_type);

        std::vector<tool::Excitation> no_excitations;
        if (!solver.setup(mesh_data, materials, modified_bcs, no_excitations)) {
            FEEM_ERROR("CapacitanceCalculator: 求解器setup失败, 导体索引={}", i);
            return Eigen::MatrixXd();
        }

        if (!solver.solve()) {
            FEEM_ERROR("CapacitanceCalculator: 求解失败, 导体索引={}", i);
            return Eigen::MatrixXd();
        }

        if (!solver.postProcess()) {
            FEEM_ERROR("CapacitanceCalculator: 后处理失败, 导体索引={}", i);
            return Eigen::MatrixXd();
        }

        const FieldData& field_data = solver.getFieldData();

        for (int j = 0; j < N; ++j) {
            double charge = computeConductorCharge(
                field_data, mesh_data, conductor_infos_[j],
                boundary_face_map, dim_type);
            C(j, i) = charge;
        }

        FEEM_INFO("CapacitanceCalculator: 导体{}求解完成, 电荷分布:", i);
        for (int j = 0; j < N; ++j) {
            FEEM_INFO("  Q[{}] = {:.6e}", j, C(j, i));
        }
    }

    C = (C + C.transpose()) / 2.0;

    capacitance_matrix_ = C;

    FEEM_INFO("CapacitanceCalculator: 电容矩阵计算完成 ({}x{})", N, N);
    for (int i = 0; i < N; ++i) {
        std::string row_str;
        for (int j = 0; j < N; ++j) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%+.6e ", C(i, j));
            row_str += buf;
        }
        FEEM_INFO("  C[{}] = {}", i, row_str);
    }

    return C;
}

double CapacitanceCalculator::computeCapacitance(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    int conductor1_id,
    int conductor2_id,
    tool::DimType dim_type) {
    std::vector<int> conductor_ids = {conductor1_id, conductor2_id};

    Eigen::MatrixXd C = computeCapacitanceMatrix(
        mesh_data, materials, boundaries, conductor_ids, dim_type);

    if (C.size() == 0) {
        FEEM_ERROR("CapacitanceCalculator::computeCapacitance 电容矩阵计算失败");
        return 0.0;
    }

    double mutual_capacitance = -C(0, 1);

    if (mutual_capacitance < 0.0) {
        FEEM_WARN("CapacitanceCalculator::computeCapacitance 互电容为负值({:.6e}), "
                  "可能存在网格或边界条件问题", mutual_capacitance);
    }

    FEEM_INFO("CapacitanceCalculator::computeCapacitance 导体{}与导体{}之间的互电容 = {:.6e}",
              conductor1_id, conductor2_id, mutual_capacitance);

    return mutual_capacitance;
}

const std::vector<ConductorInfo>& CapacitanceCalculator::getConductorInfos() const {
    return conductor_infos_;
}

void CapacitanceCalculator::clear() {
    conductor_infos_.clear();
    capacitance_matrix_.resize(0, 0);
    original_boundaries_.clear();
}

/* ==================== 私有方法 ==================== */

bool CapacitanceCalculator::extractConductorInfos(
    const std::vector<tool::Boundary>& boundaries,
    const fe_em::EMMeshData& mesh_data,
    const std::vector<int>& conductor_ids) {
    conductor_infos_.clear();

    std::map<std::string, const tool::Boundary*> name_to_boundary;
    for (const auto& bnd : boundaries) {
        name_to_boundary[bnd.getName()] = &bnd;
    }

    for (int cid : conductor_ids) {
        ConductorInfo info;
        info.id = cid;

        const fe_em::EMBoundaryMarker* matched_marker = nullptr;
        for (const auto& marker : mesh_data.boundary_markers) {
            if (marker.id == cid) {
                matched_marker = &marker;
                break;
            }
        }

        if (!matched_marker) {
            FEEM_ERROR("CapacitanceCalculator::extractConductorInfos "
                       "未找到ID={}的边界标记", cid);
            return false;
        }

        info.name = matched_marker->name;
        info.bnd_type = matched_marker->bnd_type;

        if (matched_marker->bnd_type != tool::BndType::DIRICHLET
            && matched_marker->bnd_type != tool::BndType::PERFECT_E) {
            FEEM_WARN("CapacitanceCalculator::extractConductorInfos "
                       "边界标记ID={}的类型不是DIRICHLET或PERFECT_E (类型={}), "
                       "仍尝试作为导体处理", cid,
                       static_cast<int>(matched_marker->bnd_type));
        }

        if (std::holds_alternative<std::vector<int>>(matched_marker->target_ids)) {
            info.node_ids = std::get<std::vector<int>>(matched_marker->target_ids);
        } else if (std::holds_alternative<std::vector<std::vector<int>>>(
                       matched_marker->target_ids)) {
            const auto& face_ids =
                std::get<std::vector<std::vector<int>>>(matched_marker->target_ids);
            std::set<int> unique_nodes;
            for (const auto& face : face_ids) {
                for (int nid : face) {
                    unique_nodes.insert(nid);
                }
            }
            info.node_ids.assign(unique_nodes.begin(), unique_nodes.end());
        } else {
            FEEM_ERROR("CapacitanceCalculator::extractConductorInfos "
                       "边界标记ID={}的target_ids为空", cid);
            return false;
        }

        auto it = name_to_boundary.find(matched_marker->name);
        if (it != name_to_boundary.end()) {
            info.voltage = it->second->getVoltage();
        } else {
            info.voltage = matched_marker->value;
        }

        FEEM_INFO("CapacitanceCalculator: 提取导体信息 ID={}, 名称={}, "
                  "边界节点数={}, 电压={:.6e}",
                  info.id, info.name, info.node_ids.size(), info.voltage);

        conductor_infos_.push_back(info);
    }

    return true;
}

std::vector<tool::Boundary> CapacitanceCalculator::buildBoundaryConditions(
    int conductor_index,
    int num_conductors) const {
    std::vector<tool::Boundary> result;

    std::set<std::string> conductor_names;
    for (int i = 0; i < num_conductors; ++i) {
        conductor_names.insert(conductor_infos_[i].name);
    }

    for (const auto& orig_bnd : original_boundaries_) {
        if (conductor_names.find(orig_bnd.getName()) != conductor_names.end()) {
            continue;
        }
        result.push_back(orig_bnd);
    }

    for (int i = 0; i < num_conductors; ++i) {
        tool::Boundary bnd(conductor_infos_[i].name);
        bnd.setType(conductor_infos_[i].bnd_type);

        double voltage = (i == conductor_index) ? 1.0 : 0.0;
        bnd.setVoltage(voltage);

        result.push_back(bnd);
    }

    return result;
}

} // namespace solver
