/**
 * @file solver_scheduler.hpp
 * @brief 求解器层 - 求解器调度器头文件
 * @details 负责求解流程的统一调度与管理，是连接应用层(SolverApp)与物理场求解器的桥梁。
 *          核心职责：
 *          1. 从ProjectManager提取求解所需全部数据
 *          2. 根据SimulationType创建对应的PhysicsField实例
 *          3. 驱动setup→solve→postProcess完整流程
 *          4. 协调结果输出
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "physics_field.hpp"
#include "field_data.hpp"
#include <memory>
#include <string>

namespace tool { class ProjectManager; }

namespace solver {

/**
 * @class SolverScheduler
 * @brief 求解器调度器
 * @details 统一调度各类物理场求解器的执行流程。
 *          设计模式：策略模式（Strategy Pattern）—— 通过多态调用不同物理场求解器。
 */
class SolverScheduler {
public:
    SolverScheduler();
    ~SolverScheduler();

    /**
     * @brief 初始化调度器，从ProjectManager提取数据
     * @param pm 项目管理器引用
     * @return bool 初始化成功返回true
     */
    bool initialize(tool::ProjectManager& pm);

    /**
     * @brief 执行完整求解流程
     * @return bool 求解成功返回true
     */
    bool run();

    /**
     * @brief 设置自定义物理场求解器（覆盖自动创建逻辑）
     * @param field 物理场求解器智能指针
     */
    void setPhysicsField(std::unique_ptr<PhysicsField> field);

    /**
     * @brief 获取场数据结果（只读）
     * @return const FieldData& 场数据引用
     */
    const FieldData& getFieldData() const;

    /**
     * @brief 获取当前物理场求解器（只读）
     * @return const PhysicsField* 求解器指针（可能为nullptr）
     */
    const PhysicsField* getPhysicsField() const;

    /**
     * @brief 清空所有数据和资源
     */
    void clear();

private:
    std::unique_ptr<PhysicsField> physics_field_;   ///< 当前物理场求解器实例
    FieldData field_data_;                          ///< 场数据缓存
    tool::DimType dim_type_;                        ///< 维度类型
    tool::SimulationType sim_type_;                 ///< 求解类型

    /**
     * @brief 根据SimulationType和DimType创建对应的物理场求解器
     * @param sim_type 求解类型枚举
     * @param dim_type 维度类型枚举
     * @return std::unique_ptr<PhysicsField> 求解器实例，不支持的类型返回nullptr
     */
    std::unique_ptr<PhysicsField> createPhysicsField(
        tool::SimulationType sim_type,
        tool::DimType dim_type
    ) const;

    /**
     * @brief 从ProjectManager的Material列表构建MaterialProperties映射表
     * @details 将tool::Material对象转换为numeric::MaterialProperties，
     *          以material_id为键构建映射表供求解器使用。
     * @param pm 项目管理器引用
     * @return std::map<int, MaterialProperties> material_id → 材料属性
     */
    static std::map<int, numeric::MaterialProperties> buildMaterialMap(
        tool::ProjectManager& pm
    );

    /**
     * @brief 从ProjectManager的Boundary列表构建边界条件向量
     * @param pm 项目管理器引用
     * @return std::vector<tool::Boundary> 边界条件对象副本列表
     */
    static std::vector<tool::Boundary> buildBoundaryList(
        tool::ProjectManager& pm
    );

    /**
     * @brief 从ProjectManager的Excitation列表构建激励源向量
     * @param pm 项目管理器引用
     * @return std::vector<tool::Excitation> 激励源对象副本列表
     */
    static std::vector<tool::Excitation> buildExcitationList(
        tool::ProjectManager& pm
    );
};

} // namespace solver
