/**
 * @file element_geometry.hpp
 * @brief 电磁物理层 - 有限元单元几何定义工具类头文件
 * @details 提供标准有限元单元的局部棱边和面定义，硬编码遵循ANSYS Fluent/VTK标准顺序。
 *          为Nedelec棱单元、DOF管理模块、边界条件施加提供几何基础。
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0 (v2.0: 重构为枚举类型参数)
 * 
 * @note 标准来源说明：
 *       - 棱边定义：遵循VTK/ANSYS Fluent标准，局部节点ID按升序排列
 *       - 面定义：遵循右手定则，外法向朝外
 *       - 节点顺序：与Gmsh/ANSYS/VTK通用格式兼容
 */

#pragma once

#include <vector>
#include <tuple>
#include <stdexcept>
#include "em_mesh_data.hpp"  // 引入ElemType类型别名

namespace fe_em {

/**
 * @class ElementGeometry
 * @brief 有限元单元局部几何定义工具类（v2.0: 强类型枚举版本）
 * @details 提供静态方法查询各种单元类型的：
 *          - 局部棱边定义（每条棱边的两个端点局部节点ID）
 *          - 局部面定义（每个面的边界节点局部ID列表）
 *          - 基本拓扑信息（节点数、棱边数、面数）
 * 
 * @note 设计原则：
 *       - 所有方法为static，无状态持有，线程安全
 *       - 硬编码标准定义，运行时零开销查询（O(1)枚举索引）
 *       - 使用ElemType强类型枚举替代std::string，节省90%内存
 *       - 支持一阶和二阶单元共12种Lagrange类型 + Nedelec扩展
 * 
 * @note 性能对比（v1.0 vs v2.0）：
 *       - 内存: string(32-56B) → enum(4B)，节省90%+
 *       - 比较: O(n)字符串匹配 → O(1)整数比较
 *       - 类型安全: 运行时拼写错误 → 编译期强制检查
 * 
 * @note 与后续模块对接：
 *       - 形函数模块：直接传递ElemType给ShapeFunctionBase，无需转换
 *       - Nedelec棱单元：使用get_local_edges()确定棱边自由度位置
 *       - DOF管理模块：根据单元类型分配全局DOF编号
 *       - 边界条件施加：使用get_local_faces()识别边界面单元
 */
class ElementGeometry {
public:
    /**
     * @brief 获取单元的局部棱边定义
     * @param elem_type 单元类型枚举（如 ElemType::TET4, ElemType::HEX8）
     * @return std::vector<std::tuple<int, int>> 棱边列表，每条棱边为(局部节点ID1, 局部节点ID2)升序元组
     * @throws std::invalid_argument 如果单元类型不支持
     * 
     * @note 棱边定义标准：
     *       - 每条棱边的两个节点ID按升序排列（小ID在前）
     *       - 棱边顺序遵循VTK/ANSYS Fluent标准
     *       - 用于Nedelec棱单元的自由度定义
     */
    static std::vector<std::tuple<int, int>> get_local_edges(ElemType elem_type);
    
    /**
     * @brief 获取单元的局部面定义
     * @param elem_type 单元类型枚举
     * @return std::vector<std::vector<int>> 面列表，每个面为一组有序的局部节点ID
     * @throws std::invalid_argument 如果单元类型不支持
     * 
     * @note 面定义标准：
     *       - 面的节点按逆时针顺序排列（从外部看）
     *       - 遵循右手定则，外法向朝外
     *       - 二维单元只有1个面（即单元本身）
     */
    static std::vector<std::vector<int>> get_local_faces(ElemType elem_type);
    
    /**
     * @brief 获取单元的节点数
     * @param elem_type 单元类型枚举
     * @return int 节点数量
     * @throws std::invalid_argument 如果单元类型不支持
     */
    static int get_num_nodes(ElemType elem_type);
    
    /**
     * @brief 获取单元的棱边数
     * @param elem_type 单元类型枚举
     * @return int 棱边数量
     * @throws std::invalid_argument 如果单元类型不支持
     */
    static int get_num_edges(ElemType elem_type);
    
    /**
     * @brief 获取单元的面数
     * @param elem_type 单元类型枚举
     * @return int 面数量（二维单元返回1）
     * @throws std::invalid_argument 如果单元类型不支持
     */
    static int get_num_faces(ElemType elem_type);
    
    /**
     * @brief 检查单元类型是否受支持
     * @param elem_type 单元类型枚举
     * @return bool 如果支持返回true
     */
    static bool is_supported(ElemType elem_type);
    
    /**
     * @brief 获取所有支持的单元类型列表
     * @return std::vector<ElemType> 支持的单元类型枚举列表
     */
    static std::vector<ElemType> get_supported_types();

private:
    // 禁止实例化
    ElementGeometry() = delete;
};

} // namespace fe_em
