/**
 * @file resistance_calculator.cpp
 * @brief 求解器层 - 电阻计算器实现文件
 * @details 基于欧姆定律法实现任意两端子间直流电阻的计算。
 *          核心流程：对高压端施加单位电压求解直流电流场，
 *          通过高斯面积分提取端面总电流，由 R = V/I 计算电阻。
 *
 *          欧姆定律法原理：
 *          I = ∫_Γ J·n dS = ∫_Γ σE·n dS
 *          当V_1 = 1V, V_2 = 0V时：
 *          R = V / I = 1.0 / I （I为端子1的总流出电流）
 *
 *          与CapacitanceCalculator的对比：
 *          - 电容：Q = ∫ D·dS, C = Q/V （电荷法）
 *          - 电阻：I = ∫ J·dS, R = V/I （欧姆定律法）
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include "resistance_calculator.hpp"
#include "dc_conduction_solver.hpp"
#include "field_data.hpp"
#include "logger_factory.hpp"
#include "em_enums.hpp"
#include "mesh_query.hpp"
#include "element_geometry.hpp"

#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <unordered_set>

namespace solver {

/* ==================== 内部辅助函数（匿名命名空间，避免-Wunused-function） ==================== */

namespace {

/**
 * @brief 判断单元类型是否为二维单元
 * @param elem_type 单元类型枚举
 * @return bool 二维单元返回true
 */
bool is2DElement(fe_em::ElemType elem_type) {
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
bool is3DElement(fe_em::ElemType elem_type) {
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
std::vector<std::vector<int>> get2DEdgeDefinitions(fe_em::ElemType elem_type) {
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
int getCornerNodeCount(int face_node_count) {
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
Eigen::Vector3d getNodeCoords(int node_id, const fe_em::EMMeshData& mesh_data) {
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
FaceToElementMap buildBoundaryFaceMap(const fe_em::EMMeshData& mesh_data) {
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
Eigen::Vector3d computeElementCentroid(
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
void computeEdgeNormalAndLength(
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
void computeTriangleNormalAndArea(
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
void computeQuadNormalAndArea(
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

} // anonymous namespace

/* ==================== 构造函数 ==================== */

ResistanceCalculator::ResistanceCalculator()
    : terminal_infos_()
    , last_resistance_(0.0)
    , original_boundaries_() {
}

/* ==================== 公共方法 ==================== */

double ResistanceCalculator::computeResistance(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials,
    const std::vector<tool::Boundary>& boundaries,
    int terminal1_id,
    int terminal2_id,
    tool::DimType dim_type) {
    clear();

    FEEM_INFO("ResistanceCalculator::computeResistance 开始计算电阻 "
              "(端子1_ID={}, 端子2_ID={}, 维度={})",
              terminal1_id, terminal2_id,
              static_cast<int>(dim_type));

    if (mesh_data.elements.empty()) {
        FEEM_ERROR("ResistanceCalculator::computeResistance 网格单元为空");
        return -1.0;
    }

    if (!extractTerminalInfos(boundaries, mesh_data, terminal1_id, terminal2_id)) {
        FEEM_ERROR("ResistanceCalculator::computeResistance 端子信息提取失败");
        return -1.0;
    }

    original_boundaries_ = boundaries;

    FEEM_INFO("ResistanceCalculator: 构建边界面映射...");
    FaceToElementMap boundary_face_map = buildBoundaryFaceMap(mesh_data);
    FEEM_INFO("ResistanceCalculator: 边界面映射构建完成, 边界面数量={}",
              boundary_face_map.size());

    std::vector<tool::Boundary> resistance_bcs = buildResistanceBCs(terminal1_id, terminal2_id, 1.0);

    DCConductionSolver solver(dim_type);

    std::vector<tool::Excitation> no_excitations;
    if (!solver.setup(mesh_data, materials, resistance_bcs, no_excitations)) {
        FEEM_ERROR("ResistanceCalculator: 求解器setup失败");
        return -1.0;
    }

    if (!solver.solve()) {
        FEEM_ERROR("ResistanceCalculator: 求解失败");
        return -1.0;
    }

    if (!solver.postProcess()) {
        FEEM_ERROR("ResistanceCalculator: 后处理失败");
        return -1.0;
    }

    const FieldData& field_data = solver.getFieldData();

    const auto& elem_J = field_data.getElementCurrentDensity();
    if (elem_J.empty()) {
        FEEM_ERROR("ResistanceCalculator: 电流密度数据为空");
        return -1.0;
    }

    double total_current = computeTerminalCurrentByGaussIntegral(
        field_data, mesh_data, terminal_infos_[0],
        boundary_face_map, dim_type);

    if (std::abs(total_current) < 1e-15) {
        FEEM_ERROR("ResistanceCalculator: 端子1的总电流接近零({:.6e}), "
                   "可能存在绝缘或边界条件错误", total_current);
        return -1.0;
    }

    double voltage_applied = 1.0;
    double resistance = voltage_applied / total_current;

    terminal_infos_[0].current = total_current;
    terminal_infos_[0].voltage = voltage_applied;
    terminal_infos_[1].current = -total_current;
    terminal_infos_[1].voltage = 0.0;

    last_resistance_ = resistance;

    FEEM_INFO("ResistanceCalculator::computeResistance 电阻计算完成:");
    FEEM_INFO("  端子1 (ID={}, 名称={}): V={:.6f} V, I={:.6e} A",
              terminal_infos_[0].id, terminal_infos_[0].name,
              terminal_infos_[0].voltage, terminal_infos_[0].current);
    FEEM_INFO("  端子2 (ID={}, 名称={}): V={:.6f} V, I={:.6e} A",
              terminal_infos_[1].id, terminal_infos_[1].name,
              terminal_infos_[1].voltage, terminal_infos_[1].current);
    FEEM_INFO("  电阻 R = V/I = {:.6e} Ω", resistance);

    return resistance;
}

double ResistanceCalculator::computeTotalCurrent(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<Eigen::Vector3d>& elem_current_density,
    int terminal_region_id) {
    const fe_em::EMBoundaryMarker* matched_marker = nullptr;
    for (const auto& marker : mesh_data.boundary_markers) {
        if (marker.id == terminal_region_id) {
            matched_marker = &marker;
            break;
        }
    }

    if (!matched_marker) {
        FEEM_ERROR("ResistanceCalculator::computeTotalCurrent "
                   "未找到ID={}的边界标记", terminal_region_id);
        return 0.0;
    }

    std::unordered_set<int> terminal_nodes;
    if (std::holds_alternative<std::vector<int>>(matched_marker->target_ids)) {
        const auto& node_ids = std::get<std::vector<int>>(matched_marker->target_ids);
        terminal_nodes.insert(node_ids.begin(), node_ids.end());
    } else if (std::holds_alternative<std::vector<std::vector<int>>>(
                   matched_marker->target_ids)) {
        const auto& face_ids =
            std::get<std::vector<std::vector<int>>>(matched_marker->target_ids);
        for (const auto& face : face_ids) {
            for (int nid : face) {
                terminal_nodes.insert(nid);
            }
        }
    }

    if (terminal_nodes.empty()) {
        FEEM_WARN("ResistanceCalculator::computeTotalCurrent "
                  "端子区域ID={}的节点集合为空", terminal_region_id);
        return 0.0;
    }

    /* ========== 边界面积分法：I = ∫_Γ J·n dS ==========
     * 遍历所有单元，找出与端子共享边/面的单元，
     * 对每个外边界面计算 I_face = J·n × Area 并累加 */

    double total_current = 0.0;
    int contributing_faces = 0;

    for (size_t e = 0; e < mesh_data.elements.size(); ++e) {
        const auto& elem = mesh_data.elements[e];
        if (elem.dof_type != fe_em::DOFType::SCALAR_ONLY) continue;
        if (e >= elem_current_density.size()) continue;

        Eigen::Vector3d J_e = elem_current_density[e];

        /* 获取单元局部面定义 */
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

        /* 对每个面检查是否属于端子 */
        for (int f_idx = 0; f_idx < static_cast<int>(local_faces.size()); ++f_idx) {
            const auto& face_def = local_faces[f_idx];

            /* 检查面的所有角节点是否都属于端子 */
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

            /* 获取面角节点坐标 */
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

            /* 计算外法向和面积 */
            Eigen::Vector3d normal(0, 0, 0);
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
            total_current += current_contribution;
            ++contributing_faces;
        }
    }

    FEEM_DEBUG("ResistanceCalculator::computeTotalCurrent 边界面积分结果: "
               "贡献面数={}, 总电流={:.6e} A",
               contributing_faces, total_current);

    return total_current;
}

const std::vector<TerminalInfo>& ResistanceCalculator::getTerminalInfos() const {
    return terminal_infos_;
}

void ResistanceCalculator::clear() {
    terminal_infos_.clear();
    last_resistance_ = 0.0;
    original_boundaries_.clear();
}

/* ==================== 私有方法 ==================== */

bool ResistanceCalculator::extractTerminalInfos(
    const std::vector<tool::Boundary>& boundaries,
    const fe_em::EMMeshData& mesh_data,
    int terminal1_id,
    int terminal2_id) {
    terminal_infos_.clear();

    std::map<std::string, const tool::Boundary*> name_to_boundary;
    for (const auto& bnd : boundaries) {
        name_to_boundary[bnd.getName()] = &bnd;
    }

    std::vector<int> terminal_ids = {terminal1_id, terminal2_id};

    for (int tid : terminal_ids) {
        TerminalInfo info;
        info.id = tid;

        const fe_em::EMBoundaryMarker* matched_marker = nullptr;
        for (const auto& marker : mesh_data.boundary_markers) {
            if (marker.id == tid) {
                matched_marker = &marker;
                break;
            }
        }

        if (!matched_marker) {
            FEEM_ERROR("ResistanceCalculator::extractTerminalInfos "
                       "未找到ID={}的边界标记", tid);
            return false;
        }

        info.name = matched_marker->name;

        if (matched_marker->bnd_type != tool::BndType::DIRICHLET
            && matched_marker->bnd_type != tool::BndType::PERFECT_E) {
            FEEM_WARN("ResistanceCalculator::extractTerminalInfos "
                       "边界标记ID={}的类型不是DIRICHLET或PERFECT_E (类型={}), "
                       "仍尝试作为端子处理", tid,
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
            FEEM_ERROR("ResistanceCalculator::extractTerminalInfos "
                       "边界标记ID={}的target_ids为空", tid);
            return false;
        }

        auto it = name_to_boundary.find(matched_marker->name);
        if (it != name_to_boundary.end()) {
            info.voltage = it->second->getVoltage();
        } else {
            info.voltage = matched_marker->value;
        }

        FEEM_INFO("ResistanceCalculator: 提取端子信息 ID={}, 名称={}, "
                  "边界节点数={}, 电压={:.6e}",
                  info.id, info.name, info.node_ids.size(), info.voltage);

        terminal_infos_.push_back(info);
    }

    return true;
}

std::vector<tool::Boundary> ResistanceCalculator::buildResistanceBCs(
    int terminal1_id,
    int terminal2_id,
    double voltage) const {
    std::vector<tool::Boundary> result;

    std::set<std::string> terminal_names;
    for (const auto& term_info : terminal_infos_) {
        terminal_names.insert(term_info.name);
    }

    for (const auto& orig_bnd : original_boundaries_) {
        if (terminal_names.find(orig_bnd.getName()) != terminal_names.end()) {
            continue;
        }
        result.push_back(orig_bnd);
    }

    for (size_t i = 0; i < terminal_infos_.size(); ++i) {
        tool::Boundary bnd(terminal_infos_[i].name);
        bnd.setType(tool::BndType::DIRICHLET);

        double terminal_voltage = 0.0;
        if (terminal_infos_[i].id == terminal1_id) {
            terminal_voltage = voltage;
        } else if (terminal_infos_[i].id == terminal2_id) {
            terminal_voltage = 0.0;
        }

        bnd.setVoltage(terminal_voltage);

        result.push_back(bnd);
    }

    return result;
}

double ResistanceCalculator::computeTerminalCurrentByGaussIntegral(
    const FieldData& field_data,
    const fe_em::EMMeshData& mesh_data,
    const TerminalInfo& terminal_info,
    const FaceToElementMap& boundary_face_map,
    tool::DimType dim_type) {
    std::unordered_set<int> terminal_nodes(
        terminal_info.node_ids.begin(), terminal_info.node_ids.end());

    const auto& J_field = field_data.getElementCurrentDensity();

    double total_current = 0.0;

    for (const auto& [key, face_info] : boundary_face_map) {
        bool on_terminal = true;
        for (int node_id : key) {
            if (terminal_nodes.find(node_id) == terminal_nodes.end()) {
                on_terminal = false;
                break;
            }
        }
        if (!on_terminal) continue;

        int elem_idx = face_info.element_index;
        if (elem_idx < 0 || elem_idx >= static_cast<int>(J_field.size())) {
            continue;
        }

        Eigen::Vector3d J = J_field[elem_idx];

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

        double J_dot_n = J.dot(normal);
        double current_contribution = J_dot_n * area;

        if (dim_type == tool::DimType::AXIS) {
            double r_center = 0.0;
            for (const auto& coord : corner_coords) {
                r_center += coord[0];
            }
            r_center /= static_cast<double>(corner_coords.size());
            current_contribution *= 2.0 * M_PI * r_center;
        }

        total_current += current_contribution;
    }

    return total_current;
}

} // namespace solver
