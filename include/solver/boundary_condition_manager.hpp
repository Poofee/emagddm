/**
 * @file boundary_condition_manager.hpp
 * @brief 求解器层 - 边界条件管理器头文件
 * @details 将tool::Boundary转化为对系统矩阵和右端向量的修改，
 *          支持Dirichlet/Neumann/Robin/周期性等多种边界条件类型。
 *          采用消去法（Elimination）处理Dirichlet边界条件。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "em_enums.hpp"
#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"
#include "em_dof_manager.hpp"
#include "em_element_integrator_base.hpp"
#include "project_data.hpp"
#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>

namespace numeric { 
    template<typename T> class CsrMatrix;
    template<typename T> class CooMatrix;
}

namespace solver {

/**
 * @struct DirichletDOF
 * @brief Dirichlet约束DOF信息
 */
struct DirichletDOF {
    int global_dof_id;      ///< 全局DOF编号（预编号，-1表示未知）
    double value;           ///< 约束值（电位值）
    int boundary_id;        ///< 来源边界条件ID
    int node_id;            ///< 来源节点ID（用于RHS修正时查找约束值）
};

/**
 * @class BoundaryConditionManager
 * @brief 边界条件管理器
 * @details 负责将tool::Boundary转化为对有限元系统的修改：
 *          1. 解析Boundary列表，提取Dirichlet/Neumann/Robin等条件
 *          2. 将Dirichlet条件映射到DOF约束（配合EMDOFManager）
 *          3. 修正右端项向量（消去法：F_free -= K_fc * phi_c）
 *          4. 处理边界条件冲突检测
 *
 * @note 消去法RHS修正原理：
 *       组装阶段EMAssembly跳过了约束DOF（indices==-1）的行和列，
 *       因此组装后的K_csr仅包含自由DOF行和列。
 *       消去法要求对右端项进行修正：F_free = F_free - K_fc * phi_c，
 *       其中K_fc是刚度矩阵中自由DOF行、约束DOF列的子矩阵。
 *       由于K_csr中不包含约束DOF列的信息，需要重新计算单元矩阵K_e
 *       来获取K_fc的贡献，逐单元修正右端项向量。
 */
class BoundaryConditionManager {
public:
    /**
     * @brief 默认构造函数
     */
    BoundaryConditionManager();

    /**
     * @brief 默认析构函数
     */
    ~BoundaryConditionManager() = default;

    /**
     * @brief 从Boundary列表中提取所有Dirichlet约束DOF
     * @param boundaries 边界条件列表（tool层定义）
     * @param mesh_data 网格拓扑数据（包含EMBoundaryMarker，提供target_ids）
     * @param dof_manager DOF管理器（提供Local2Global映射）
     * @return bool 提取成功返回true
     *
     * @details 处理流程：
     *          1. 遍历mesh_data.boundary_markers，通过mapToElectrostaticBndType
     *             将各类边界类型映射为静电场等效类型
     *          2. 对Dirichlet等效类型（DIRICHLET/PERFECT_E/ODD_SYMMETRY），
     *             从marker.target_ids提取约束节点ID和约束值
     *          3. 构建node_dirichlet_values_映射（节点ID→约束值）
     *          4. 构建node_to_free_dof_映射（节点ID→自由DOF编号）
     *          5. 构建free_dof_to_prenum_映射（自由DOF编号→预编号）
     *          6. 执行冲突检测（同一节点不允许施加不同值的Dirichlet条件）
     */
    bool extractDirichletDOFs(
        const std::vector<tool::Boundary>& boundaries,
        const fe_em::EMMeshData& mesh_data,
        const fe_em::EMDOFManager& dof_manager
    );

    /**
     * @brief 修正右端项向量（消去法）
     * @details 对每个单元重新计算K_e，对自由DOF行i和约束DOF列j，
     *          执行 rhs[indices[i]] -= K_e(i,j) * phi_constrained_j
     *          其中phi_constrained_j从node_dirichlet_values_中查找。
     * @param rhs 右端项向量（仅自由DOF部分，会被修改）
     * @param mesh_data 网格拓扑数据（提供单元节点ID和坐标）
     * @param elem_l2g 所有单元的Local2Global映射表
     * @param materials 材料属性列表（用于创建积分器）
     */
    void applyDirichletToRHS(
        Eigen::VectorXd& rhs,
        const fe_em::EMMeshData& mesh_data,
        const std::vector<fe_em::Local2Global>& elem_l2g,
        const std::vector<numeric::MaterialProperties>& materials
    );

    /**
     * @brief 将约束DOF值回扩到完整解向量
     * @param phi_free 自由DOF解向量
     * @param dof_manager DOF管理器
     * @return Eigen::VectorXd 完整解向量（包含约束DOF值，按预编号顺序）
     *
     * @details 利用free_dof_to_prenum_映射将自由DOF解值填入正确位置，
     *          约束DOF位置填入node_dirichlet_values_中的约束值。
     */
    Eigen::VectorXd expandSolution(
        const Eigen::VectorXd& phi_free,
        const fe_em::EMDOFManager& dof_manager
    ) const;

    /**
     * @brief 获取所有Dirichlet约束DOF列表
     * @return const std::vector<DirichletDOF>& 约束DOF列表
     */
    const std::vector<DirichletDOF>& getDirichletDOFs() const;

    /**
     * @brief 获取节点Dirichlet值映射
     * @return const std::map<int, double>& 节点ID→约束值映射
     */
    const std::map<int, double>& getNodeDirichletValues() const;

    /**
     * @brief 检查是否存在至少一个Dirichlet边界条件
     * @return bool 存在返回true
     */
    bool hasDirichletBC() const;

    /**
     * @brief 清空所有数据
     */
    void clear();

private:
    std::vector<DirichletDOF> dirichlet_dofs_;       ///< Dirichlet约束DOF列表
    std::map<int, double> constrained_values_;        ///< 约束DOF值映射（global_dof_id → value）
    std::map<int, double> node_dirichlet_values_;     ///< 节点Dirichlet值映射（node_id → value）
    std::map<int, int> node_to_free_dof_;             ///< 节点到自由DOF编号映射（node_id → free_dof_index，-1表示约束）
    std::map<int, int> free_dof_to_prenum_;           ///< 自由DOF编号到预编号映射（free_dof_index → prenum）
    int total_dofs_;                                  ///< 总DOF数（预编号总数）

    /**
     * @brief 检测Dirichlet边界条件冲突
     * @details 同一节点不能施加不同值的Dirichlet条件
     * @return bool 无冲突返回true
     */
    bool detectConflicts() const;

    /**
     * @brief 将BndType映射为静电场等效边界条件
     * @param bnd_type 原始边界类型
     * @return BndType 等效边界类型
     *
     * @details 映射规则：
     *          - DIRICHLET → DIRICHLET（直接使用）
     *          - PERFECT_E → DIRICHLET（理想电边界，φ=const）
     *          - ODD_SYMMETRY → DIRICHLET（奇对称，φ=0）
     *          - NEUMANN, INSULATION, BALLOON, EVEN_SYMMETRY → NEUMANN（自然边界条件）
     *          - ROBIN → ROBIN（需矩阵修正，暂存）
     *          - PERIODIC, ANTIPERIODIC, MASTER_SLAVE → 原样保留（高级特性）
     *          - 其他类型 → NEUMANN（默认自然边界条件）
     */
    static tool::BndType mapToElectrostaticBndType(tool::BndType bnd_type);

    /**
     * @brief 构建节点到自由DOF编号的映射
     * @param mesh_data 网格拓扑数据
     * @param dof_manager DOF管理器
     *
     * @details 遍历所有SCALAR_ONLY单元，从Local2Global映射中提取
     *          每个节点的自由DOF编号（indices[i]），构建node_to_free_dof_映射。
     *          约束节点的自由DOF编号为-1。
     */
    void buildNodeToFreeDOFMapping(
        const fe_em::EMMeshData& mesh_data,
        const fe_em::EMDOFManager& dof_manager
    );

    /**
     * @brief 构建自由DOF编号到预编号的映射
     * @param mesh_data 网格拓扑数据
     *
     * @details 重构EMDOFManager内部的预编号顺序：
     *          收集所有被SCALAR_ONLY单元引用的节点ID，按升序排列
     *          （与EMDOFManager::assignPrenumbering中的std::set遍历顺序一致），
     *          依次分配预编号，然后结合node_to_free_dof_构建free_dof_to_prenum_映射。
     */
    void buildFreeDOFToPrenumMapping(const fe_em::EMMeshData& mesh_data);
};

} // namespace solver
