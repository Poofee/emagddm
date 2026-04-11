/**
 * @file physics_field.hpp
 * @brief 求解器层 - 物理场抽象基类定义
 * @details 定义所有物理场求解器的统一接口契约，
 *          支持setup/solve/postProcess完整求解流程。
 *          所有具体物理场求解器（静电场、静磁场、涡流场等）均需继承此基类。
 * @author Poofee
 * @date 2026-04-11
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include "em_element_integrator_base.hpp"
#include "em_enums.hpp"
#include "project_data.hpp"
#include <map>
#include <vector>
#include <string>
#include <memory>

namespace solver { class FieldData; }

namespace solver {

/**
 * @class PhysicsField
 * @brief 物理场抽象基类
 * @details 定义所有物理场求解器的通用接口契约，采用模板方法模式（Template Method）：
 *          - setup(): 前处理+DOF分配+矩阵装配+边界条件+激励
 *          - solve(): 线性求解+解向量回扩
 *          - postProcess(): 场量计算+结果存储
 *
 *          继承关系：
 *          PhysicsField（本类）
 *              ├── ElectrostaticSolver（静电场求解器）
 *              ├── MagnetostaticSolver（静磁场求解器，待开发）
 *              └── EddyCurrentSolver（涡流场求解器，待开发）
 */
class PhysicsField {
public:
    virtual ~PhysicsField() = default;

    /**
     * @brief 求解器初始化与矩阵装配
     * @param mesh_data 网格拓扑数据
     * @param materials 材料参数映射表（material_id → MaterialProperties）
     * @param boundaries 边界条件列表
     * @param excitations 激励源列表
     * @return bool 初始化成功返回true
     */
    virtual bool setup(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::vector<tool::Excitation>& excitations
    ) = 0;

    /**
     * @brief 执行线性求解
     * @return bool 求解成功返回true
     */
    virtual bool solve() = 0;

    /**
     * @brief 后处理计算（场量计算、能量计算等）
     * @return bool 后处理成功返回true
     */
    virtual bool postProcess() = 0;

    /**
     * @brief 获取场数据（只读）
     * @return const FieldData& 场数据常量引用
     */
    virtual const FieldData& getFieldData() const = 0;

    /**
     * @brief 获取物理场类型
     * @return SimulationType 物理场类型枚举
     */
    virtual tool::SimulationType getSimulationType() const = 0;

    /**
     * @brief 获取求解器名称
     * @return std::string 求解器名称字符串
     */
    virtual std::string getSolverName() const = 0;

    /**
     * @brief 清空所有内部数据，释放资源
     */
    virtual void clear() = 0;
};

} // namespace solver
