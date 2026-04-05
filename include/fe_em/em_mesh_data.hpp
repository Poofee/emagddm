/**
 * @file em_mesh_data.hpp
 * @brief 电磁物理层 - 有限元网格核心拓扑数据结构定义
 * @details 定义电磁场有限元求解所需的网格数据结构，包括节点、单元、边界标记等。
 *          本模块只存储**拓扑数据**（节点坐标、单元连接关系），不存储材料属性和完整边界定义。
 *          材料属性通过 Element::material_id 索引 ProjectManager 中的 tool::Material 对象获取。
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 * 
 * @note 架构设计说明：
 *       - tool::Mesh: 存储**网格剖分配置参数**（max_size, boundary_layer等）
 *       - fe_em::EMMeshData: 存储**网格拓扑数据**（本文件定义）
 *       - tool::Material: 存储**材料属性**（μ_r, ε_r, σ, B-H曲线等）
 *       - 三者职责分离，通过ProjectManager统一管理
 */

#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <variant>
#include "em_enums.hpp"  // 复用tool::BndType枚举（消除重复定义）
#include "shape_function_base.hpp"  // 复用numeric::ElementType枚举（消除字符串类型）

namespace fe_em {

// 类型别名：复用tool层的边界类型枚举（消除重复定义）
using BndType = ::tool::BndType;  // 23种通用边界类型的别名（使用全局路径避免歧义）

// 类型别名：复用numeric层的单元类型枚举（消除std::string，节省内存90%+）
using ElemType = ::numeric::ElementType;  // 18种单元类型的强类型枚举

// ==================== 单元DOF类型枚举 ====================

/**
 * @enum DOFType
 * @brief 单元自由度类型枚举，直接对接后续DOF管理模块
 * @details 决定该单元在求解时使用的未知量类型：
 *          - SCALAR_ONLY: 标量场（电位φ、磁位ψ），用于静电/静磁/非导电区3D磁场
 *          - VECTOR_EDGE_ONLY: 矢量场（磁矢势A的切向分量），用于Nedelec棱单元
 *          - MIXED_AV: A-V混合格式（矢量A + 标量V），用于涡流场导体区域
 */
enum class DOFType {
    SCALAR_ONLY,      ///< 纯标量：2D/3D电场、2D磁场、3D非导电区（空气）
    VECTOR_EDGE_ONLY, ///< 纯矢量：3D磁场（非导电区，简化模型）
    MIXED_AV          ///< A-V混合：3D导电区（涡流场）
};

/**
 * @brief 将DOFType枚举转换为字符串表示
 * @param dof_type DOFType枚举值
 * @return std::string 字符串表示
 */
inline std::string dofTypeToString(DOFType dof_type) {
    switch (dof_type) {
        case DOFType::SCALAR_ONLY:      return "SCALAR_ONLY";
        case DOFType::VECTOR_EDGE_ONLY: return "VECTOR_EDGE_ONLY";
        case DOFType::MIXED_AV:          return "MIXED_AV";
        default:                        return "UNKNOWN";
    }
}

/**
 * @brief 将字符串转换为DOFType枚举
 * @param str 字符串表示
 * @return DOFType 枚举值
 */
inline DOFType stringToDOFType(const std::string& str) {
    if (str == "SCALAR_ONLY")      return DOFType::SCALAR_ONLY;
    if (str == "VECTOR_EDGE_ONLY") return DOFType::VECTOR_EDGE_ONLY;
    if (str == "MIXED_AV")          return DOFType::MIXED_AV;
    throw std::invalid_argument("未知的DOFType字符串: " + str);
}

// ==================== 节点结构体 ====================

/**
 * @struct Node
 * @brief 有限元网格节点数据结构
 * @details 存储节点的三维坐标和区域信息。
 *          二维问题中z坐标设为0。
 * 
 * @note 与后续模块对接：
 *       - 形函数模块：通过node_ids获取节点坐标计算形函数值
 *       - 单元矩阵模块：节点坐标用于雅可比矩阵计算
 *       - 后处理模块：节点坐标用于结果可视化
 */
struct Node {
    int id;           ///< 节点全局ID（从1开始或从0开始，需统一）
    double x;         ///< X坐标 (m)
    double y;         ///< Y坐标 (m)
    double z;      
     //< Z坐标 (m)，二维问题时z=0
    int region_id;    ///< 节点所属区域ID（可选，为多区域标记预留）
};

// ==================== 单元结构体 ====================

/**
 * @struct Element
 * @brief 有限元网格单元数据结构（核心对接点）
 * @details 存储单元的拓扑连接关系和属性标记。
 *          这是与后续形函数模块、单元矩阵模块、DOF管理模块的主要接口。
 * 
 * @note 关键字段说明：
 *       - node_ids: 按标准有限元顺序排列的局部节点ID列表
 *                    （顺序遵循ANSYS Fluent/VTK标准）
 *       - type: 单元类型字符串，用于ElementGeometry查询几何定义
 *       - dof_type: **关键字段**，直接告诉后续模块这是什么DOF类型
 *       - material_id: 索引ProjectManager.materials_中的tool::Material对象
 *                      （不在此处存储材料参数，避免数据冗余）
 * 
 * @note 与后续模块对接：
 *       - 形函数模块：根据type选择形函数，根据node_ids获取坐标
 *       - DOF管理模块：根据dof_type确定自由度分配策略
 *       - 单元矩阵模块：根据material_id获取μ_r, ε_r, σ计算系数矩阵
 */
struct Element {
    int id;                          ///< 单元全局ID
    std::vector<int> node_ids;       ///< 单元包含的节点ID列表（按标准有限元顺序排列）
    ElemType elem_type;              ///< 单元类型枚举（复用numeric::ElementType，强类型检查+内存优化）
    DOFType dof_type;                ///< 单元DOF类型（**关键字段**，直接对接DOF管理模块）
    int material_id;                 ///< 单元所属材料ID（索引ProjectManager的materials_）
    int region_id;                   ///< 单元所属区域ID（用于多区域问题）
};

// ==================== 边界条件标记结构体（组合设计）====================

/**
 * @struct EMBoundaryMarker
 * @brief 电磁场边界条件轻量级标记结构体（组合设计）
 * @details 采用 **BndType + DOFType** 组合替代独立的EMBoundaryType枚举，
 *          消除与tool::BndType的重复定义。
 * 
 * @note 设计理念（组合优于继承/独立枚举）：
 *       - BndType (tool层): 定义"这是什么边界"（23种通用类型）
 *         * DIRICHLET: 狄利克雷（固定值）
 *         * NEUMANN: 诺伊曼（指定通量/导数）
 *         * PERIODIC: 周期性边界
 *         * PERFECT_E/H: 理想电/磁壁
 *         ... 等19种其他类型
 *       
 *       - DOFType (fe_em层): 定义"施加到什么DOF"
 *         * SCALAR_ONLY: 标量场（φ, ψ, V）
 *         * VECTOR_EDGE_ONLY: 矢量场（A的切向分量）
 *         * MIXED_AV: A-V混合格式
 *       
 *       - target_ids (variant): 隐含定义"几何实体类型"
 *         * vector<int>: 节点列表（节点标记）
 *         * vector<vector<int>>: 面/棱边列表（面或棱边标记）
 * 
 * @note 典型组合示例（电磁场专用场景）：
 *       | 场景                  | BndType        | DOFType           | target_ids类型      |
 *       |----------------------|---------------|-------------------|---------------------|
 *       | 固定电位节点          | DIRICHLET      | SCALAR_ONLY       | vector<int>        |
 *       | 磁通量面              | NEUMANN        | SCALAR_ONLY       | vector<vector<int>>|
 *       | 理想导体棱边(A×n=0)  | PERFECT_E      | VECTOR_EDGE_ONLY | vector<vector<int>>|
 *       | 周期性边界面          | PERIODIC       | 任意              | vector<vector<int>>|
 * 
 * @note 与后续模块对接：
 *       - DOF管理模块：根据bnd_type和dof_type识别受约束的自由度
 *       - 单元组装模块：将边界条件施加到系统矩阵中
 */
struct EMBoundaryMarker {
    int id;                                              ///< 边界标记ID
    BndType bnd_type;                                     ///< 边界条件类型（复用tool::BndType，23种通用类型）
    DOFType dof_type;                                      ///< 目标DOF类型（SCALAR_ONLY/VECTOR_EDGE_ONLY/MIXED_AV）
    std::variant<
        std::vector<int>,                                ///< 节点ID列表（用于节点标记）
        std::vector<std::vector<int>>                    ///< 面/棱边的节点ID列表（用于面/棱边标记）
    > target_ids;                                        ///< 目标几何实体ID（variant隐含几何实体类型）
    double value;                                         ///< 边界条件值（标量）
    std::string name;                                     ///< 边界名称（如 "Ground", "Ideal_Conductor"）
};

// ==================== 整体网格模型聚合结构体 ====================

/**
 * @struct EMMeshData
 * @brief 电磁场有限元网格拓扑数据聚合结构体
 * @details 聚合所有网格拓扑数据，作为数值计算层的最底层数据基石。
 *          **重要**：本结构体只存储拓扑数据，不存储以下内容：
 *          - ❌ 材料属性（通过Element::material_id索引外部tool::Material）
 *          - ❌ 完整边界定义（通过EMBoundaryMarker::id索引外部tool::Boundary）
 *          - ❌ 网格配置参数（由tool::Mesh负责）
 * 
 * @note 数据访问模式：
 *       @code
 *       // 获取单元的材料相对磁导率
 *       const auto& elem = em_mesh_data.elements[0];
 *       const auto* mat = project_manager.getMaterial(elem.material_id);
 *       double mu_r = mat->getRelativePermeability();
 *       @endcode
 * 
 * @note 与后续模块对接：
 *       - MeshQuery工具类：基于nodes和elements进行几何查询
 *       - ElementGeometry工具类：基于Element.type查询局部棱边/面定义
 *       - 求解器初始化：直接读取EMMeshData构建离散系统
 *       - 项目I/O：MaxwellParser填充此结构体，序列化到项目文件
 */
struct EMMeshData {
    std::vector<Node> nodes;                              ///< 网格节点列表
    std::vector<Element> elements;                         ///< 网格单元列表
    std::vector<EMBoundaryMarker> boundary_markers;        ///< 边界条件标记列表
    
    /**
     * @brief 获取节点数量
     * @return size_t 节点总数
     */
    size_t getNodeCount() const { return nodes.size(); }
    
    /**
     * @brief 获取单元数量
     * @return size_t 单元总数
     */
    size_t getElementCount() const { return elements.size(); }
    
    /**
     * @brief 获取边界标记数量
     * @return size_t 边界标记总数
     */
    size_t getBoundaryMarkerCount() const { return boundary_markers.size(); }
    
    /**
     * @brief 清空所有数据
     */
    void clear() {
        nodes.clear();
        elements.clear();
        boundary_markers.clear();
    }
    
    /**
     * @brief 检查是否为空
     * @return bool 如果没有任何数据返回true
     */
    bool isEmpty() const {
        return nodes.empty() && elements.empty();
    }
};

} // namespace fe_em
