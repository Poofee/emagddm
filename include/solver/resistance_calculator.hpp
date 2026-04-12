/**
 * @file resistance_calculator.hpp
 * @brief 求解器层 - 电阻计算器头文件
 * @details 采用欧姆定律法计算任意两端子之间的直流电阻 R = V/I。
 *          对标ANSYS Maxwell的电阻参数提取功能。
 *
 *          计算方法（欧姆定律法）：
 *          1. 对端子i施加单位电压 V_i = 1V，端子j接地 V_j = 0V，其余绝缘
 *          2. 使用DCConductionSolver求解直流电流场问题
 *          3. 计算端子i表面的总流出电流 I_i = ∫_Γi J·dS = ∫_Γi σE·ndS
 *          4. 电阻 R_ij = V_i / I_i = 1 / I_i
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include "em_element_integrator_base.hpp"
#include "em_linear_solver.h"
#include "project_data.hpp"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>

namespace solver {

class FieldData;  ///< 前向声明

/**
 * @struct TerminalInfo
 * @brief 端子信息结构体
 * @details 存储端子的几何位置、电压和电流信息，
 *          用于电阻计算过程中的数据传递和结果存储。
 */
struct TerminalInfo {
    int id = -1;                     ///< 端子ID（对应boundary_id或region_id）
    std::string name;                ///< 端子名称
    std::vector<int> node_ids;       ///< 端子表面/区域节点ID列表
    double voltage = 0.0;            ///< 端子电压值（V）
    double current = 0.0;            ///< 端子总电流（A），计算后填充
};

/**
 * @brief 边界面到单元映射类型（使用排序后的面角节点ID列表作为键）
 */
using FaceKey = std::vector<int>;

/**
 * @brief 面信息结构体，用于边界面到单元的映射
 */
struct FaceInfo {
    int element_index = -1;     ///< 相邻单元在mesh_data.elements中的索引
    int face_local_index = -1;  ///< 面在单元局部定义中的索引
};

using FaceToElementMap = std::map<FaceKey, FaceInfo>;

/**
 * @class ResistanceCalculator
 * @brief 电阻计算器
 * @details 采用欧姆定律法 R = V / I 计算端子间直流电阻。
 *
 *          核心流程：
 *          1. extractTerminalInfos(): 从边界条件中提取两个端子的节点信息
 *          2. buildResistanceBCs(): 构建特殊边界条件（V1=1V, V2=0V）
 *          3. DCConductionSolver::setup(): 初始化求解器
 *          4. DCConductionSolver::solve(): 求解电位分布
 *          5. DCConductionSolver::postProcess(): 计算电流密度J
 *          6. computeTotalCurrent(): 计算端子1的总流出电流I
 *          7. R = V / I = 1.0 / I: 返回电阻值
 *
 *          与CapacitanceCalculator的关系：
 *          - CapacitanceCalculator采用电荷法 Q = ∫ D·dS 计算电容 C = Q/V
 *          - ResistanceCalculator采用欧姆定律法 I = ∫ J·dS 计算电阻 R = V/I
 *          - 两者结构相似，但物理量和求解器不同
 */
class ResistanceCalculator {
public:
    /**
     * @brief 默认构造函数
     * @details 初始化内部状态变量
     */
    ResistanceCalculator();

    /**
     * @brief 析构函数（默认）
     */
    ~ResistanceCalculator() = default;

    /**
     * @brief 计算两个端子之间的直流电阻（欧姆定律法）
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表（需包含电导率sigma字段）
     * @param boundaries 边界条件列表
     * @param terminal1_id 高压端端子ID（施加1V电压）
     * @param terminal2_id 接地端端子ID（接地0V）
     * @param dim_type 维度类型（2D/3D/轴对称）
     * @return double 电阻值（欧姆），失败返回-1.0
     *
     * @note 计算步骤：
     *       1. 提取两个端子的节点信息
     *       2. 构建边界条件（端子1加1V，端子2接地，其余绝缘）
     *       3. 求解直流电流场问题
     *       4. 计算端子1表面的总流出电流
     *       5. 返回 R = V/I = 1.0/I
     *
     * @warning 调用前需确保materials中包含有效的电导率值（sigma > 0）
     */
    double computeResistance(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        int terminal1_id,
        int terminal2_id,
        tool::DimType dim_type = tool::DimType::D3
    );

    /**
     * @brief 计算通过指定端面的总电流
     * @param mesh_data 网格数据
     * @param elem_current_density 单元电流密度向量列表（每个单元一个Vector3d）
     * @param terminal_region_id 端面所在区域的boundary marker ID
     * @return double 总电流（安培），正值表示流出端面
     *
     * @details 通过高斯面积分计算总电流：
     *          I = ∫_Γ J·n dS ≈ Σ_f (J_f · n_f) * A_f
     *          其中n_f为面外法向，A_f为面面积
     *
     * @note 初版实现使用简化的面积加权平均近似
     *       完整版需要构建边界面到单元的映射（参考CapacitanceCalculator）
     */
    double computeTotalCurrent(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<Eigen::Vector3d>& elem_current_density,
        int terminal_region_id
    );

    /**
     * @brief 获取最近一次计算的端子信息列表
     * @return const std::vector<TerminalInfo>& 端子信息列表的常量引用
     *
     * @note 在computeResistance()调用后可获取详细的端子信息
     */
    const std::vector<TerminalInfo>& getTerminalInfos() const;

    /**
     * @brief 清空内部数据
     * @details 重置所有内部状态变量，释放内存
     */
    void clear();

private:
    std::vector<TerminalInfo> terminal_infos_;      ///< 端子信息列表（最近一次计算）
    double last_resistance_;                         ///< 最近一次计算的电阻值（Ω）
    std::vector<tool::Boundary> original_boundaries_; ///< 原始边界条件缓存

    /**
     * @brief 从边界条件中提取端子信息
     * @param boundaries 边界条件列表
     * @param mesh_data 网格拓扑数据
     * @param terminal1_id 第一个端子ID
     * @param terminal2_id 第二个端子ID
     * @return bool 提取成功返回true，失败返回false
     *
     * @details 查找boundary_markers中匹配terminal_id的标记：
     *          1. 遍历mesh_data.boundary_markers查找匹配ID
     *          2. 提取端子名称、边界类型、关联节点
     *          3. 从boundaries中查找对应的电压值
     *          4. 填充TerminalInfo结构体
     */
    bool extractTerminalInfos(
        const std::vector<tool::Boundary>& boundaries,
        const fe_em::EMMeshData& mesh_data,
        int terminal1_id,
        int terminal2_id
    );

    /**
     * @brief 构建用于电阻计算的边界条件列表
     * @param terminal1_id 高压端端子ID
     * @param terminal2_id 接地端端子ID
     * @param voltage 施加的高压端电压值（默认1.0V）
     * @return std::vector<tool::Boundary> 修改后的边界条件列表
     *
     * @details 边界条件构建规则：
     *          1. 复制原始边界条件（排除两个端子的原有设置）
     *          2. 端子1设置为Dirichlet边界，电压=voltage（默认1.0V）
     *          3. 端子2设置为Dirichlet边界，电压=0.0V（接地）
     *          4. 其余边界保持不变（绝缘或原始设置）
     *
     * @note 此方法确保每次计算使用一致的边界条件配置
     */
    std::vector<tool::Boundary> buildResistanceBCs(
        int terminal1_id,
        int terminal2_id,
        double voltage = 1.0
    ) const;

    /**
     * @brief 通过高斯面积分计算端子表面的总电流
     * @param field_data 场数据（包含电流密度J）
     * @param mesh_data 网格拓扑数据
     * @param terminal_info 端子信息
     * @param boundary_face_map 边界面到单元的映射
     * @param dim_type 维度类型
     * @return double 总电流值（安培），正值表示流出端面
     *
     * @details 对端子表面所有边界面/边进行积分：
     *          I = Σ_f (J_f · n_f) * A_f
     *          其中n_f为单元面外法向（指向外部），
     *          A_f为面面积。
     *
     *          轴对称情况下附加2πr因子：
     *          I = Σ_f (J_f · n_f) * A_f * 2π * r_f
     */
    double computeTerminalCurrentByGaussIntegral(
        const FieldData& field_data,
        const fe_em::EMMeshData& mesh_data,
        const TerminalInfo& terminal_info,
        const FaceToElementMap& boundary_face_map,
        tool::DimType dim_type
    );
};

} // namespace solver
