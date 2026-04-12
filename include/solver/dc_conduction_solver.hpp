/**
 * @file dc_conduction_solver.hpp
 * @brief 求解器层 - 直流电流场求解器头文件
 * @details 实现完整的直流电流场有限元求解流程，对标ANSYS Maxwell DC Conduction Solver。
 *          支持二维/三维/轴对称模型，Dirichlet/Neumann/混合边界条件，
 *          电压源/电流源激励，电流密度计算，焦耳热功率计算等功能。
 *
 *          控制方程：
 *          - 二维/三维：∇·(σ∇V) = 0
 *          - 轴对称（r-z坐标系）：1/r · ∂/∂r(r·σ·∂V/∂r) + ∂/∂z(σ·∂V/∂z) = 0
 *
 *          弱形式：
 *          ∫Ω ∇N^T · σ · ∇V dΩ = ∫Γ N · (σ·∂V/∂n) dΓ
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include "physics_field.hpp"
#include "field_data.hpp"
#include "boundary_condition_manager.hpp"
#include "excitation_manager.hpp"
#include "em_dof_manager.hpp"
#include "em_assembly.hpp"
#include "em_linear_solver.h"
#include "em_solver_factory.h"
#include "dc_conduction_integrator.hpp"
#include "coo_matrix.hpp"

#include <memory>
#include <string>

namespace solver {

/**
 * @class DCConductionSolver
 * @brief 直流电流场求解器
 * @details 实现完整的直流电流场有限元求解流程，继承自PhysicsField抽象基类。
 *          核心求解流程：
 *          1. setup(): 前处理检查→材料映射→DOF分配→矩阵装配→边界条件→激励
 *          2. solve(): 选择线性求解器→求解KV=F→解向量回扩
 *          3. postProcess(): 计算电场/电流密度/焦耳热→存储到FieldData
 */
class DCConductionSolver : public PhysicsField {
public:
    DCConductionSolver();
    explicit DCConductionSolver(tool::DimType dim_type);
    ~DCConductionSolver() override = default;

    // ========== PhysicsField接口实现 ==========

    bool setup(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::vector<tool::Excitation>& excitations
    ) override;

    bool solve() override;

    bool postProcess() override;

    const FieldData& getFieldData() const override;

    tool::SimulationType getSimulationType() const override;

    std::string getSolverName() const override;

    void clear() override;

    // ========== 直流电流场专用接口 ==========

    /**
     * @brief 设置维度类型（2D/3D/轴对称）
     * @param dim_type 维度类型枚举
     */
    void setDimType(tool::DimType dim_type);

    /**
     * @brief 获取维度类型
     * @return DimType 维度类型
     */
    tool::DimType getDimType() const;

    /**
     * @brief 获取DOF管理器（只读）
     * @return const EMDOFManager& DOF管理器引用
     */
    const fe_em::EMDOFManager& getDOFManager() const;

    /**
     * @brief 获取刚度矩阵（只读）
     * @return const CsrMatrix<double>& CSR格式刚度矩阵
     */
    const numeric::CsrMatrix<double>& getStiffnessMatrix() const;

    /**
     * @brief 获取右端项向量（只读）
     * @return const Eigen::VectorXd& 右端项向量
     */
    const Eigen::VectorXd& getSourceVector() const;

    /**
     * @brief 获取完整解向量（包含约束DOF）
     * @return const Eigen::VectorXd& 完整电位向量
     */
    const Eigen::VectorXd& getFullSolution() const;

private:
    // ========== 内部状态 ==========
    tool::DimType dim_type_;                    ///< 维度类型
    bool setup_done_;                           ///< setup是否已完成
    bool solve_done_;                           ///< solve是否已完成

    // ========== 核心模块 ==========
    std::unique_ptr<fe_em::EMDOFManager> dof_manager_;       ///< DOF管理器
    std::unique_ptr<numeric::EMAssembly> assembly_;           ///< 全局组装器
    std::unique_ptr<numeric::EMLinearSolverBase> linear_solver_; ///< 线性求解器
    BoundaryConditionManager bc_manager_;                     ///< 边界条件管理器
    ExcitationManager excitation_manager_;                    ///< 激励管理器
    FieldData field_data_;                                    ///< 场数据

    // ========== 求解结果 ==========
    Eigen::VectorXd phi_free_;                 ///< 自由DOF解向量
    Eigen::VectorXd phi_full_;                 ///< 完整解向量（含约束DOF）

    // ========== 修正后的右端项 ==========
    Eigen::VectorXd f_free_;                   ///< 修正后的自由DOF右端项向量（施加边界条件后）

    // ========== 输入数据缓存 ==========
    const fe_em::EMMeshData* mesh_data_ptr_;               ///< 网格数据指针（不拥有）
    std::map<int, numeric::MaterialProperties> materials_;  ///< 材料参数映射
    std::vector<fe_em::Local2Global> elem_l2g_;             ///< DOF映射表
    std::vector<numeric::MaterialProperties> material_list_; ///< 按material_id索引的材料列表（保留原始ε和σ值，不做hack）

    // ========== 直接装配数据（绕过EMAssembly的DummyEMIntegrator） ==========
    std::unique_ptr<numeric::CsrMatrix<double>> direct_K_csr_;  ///< DCConductionIntegrator直接装配的刚度矩阵
    Eigen::VectorXd direct_F_;                                ///< 直接装配的源项向量（通常为零）

    // ========== 私有方法 ==========

    /**
     * @brief 前处理检查（含电导率合法性校验）
     * @return bool 检查通过返回true
     */
    bool preCheck();

    /**
     * @brief 从MaterialProperties构建材料映射（保留原始ε和σ值）
     * @param materials 项目材料列表
     * @note DCConductionIntegrator通过getEffectiveSigma()直接读取sigma字段，
     *       无需任何字段替换操作。
     */
    void buildMaterialMapping(
        const std::map<int, numeric::MaterialProperties>& materials
    );

    /**
     * @brief 执行DOF分配
     * @return bool 分配成功返回true
     */
    bool allocateDOFs();

    /**
     * @brief 执行全局矩阵装配（使用DCConductionIntegrator直接装配）
     * @return bool 装配成功返回true
     */
    bool assembleSystem();

    /**
     * @brief 施加边界条件和激励
     * @param boundaries 边界条件列表
     * @param excitations 激励列表
     * @return bool 施加成功返回true
     */
    bool applyBoundaryAndExcitation(
        const std::vector<tool::Boundary>& boundaries,
        const std::vector<tool::Excitation>& excitations
    );

    /**
     * @brief 选择并创建线性求解器
     * @return std::unique_ptr<EMLinearSolverBase> 求解器实例
     */
    std::unique_ptr<numeric::EMLinearSolverBase> createLinearSolver();
};

} // namespace solver
