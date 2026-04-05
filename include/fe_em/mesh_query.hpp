/**
 * @file mesh_query.hpp
 * @brief 电磁物理层 - 网格几何查询工具类头文件
 * @details 提供基于网格拓扑数据的几何查询功能，用于辅助数值计算和后处理。
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 */

#pragma once

#include "em_mesh_data.hpp"
#include <vector>

namespace fe_em {

/**
 * @class MeshQuery
 * @brief 网格几何查询工具类
 * @details 提供静态方法对EMMeshData中的拓扑数据进行查询，
 *          包括获取单元节点坐标、区域判断等操作。
 * 
 * @note 设计原则：
 *       - 所有方法为static，无状态持有，线程安全
 *       - 基于传入的数据引用进行查询，不拷贝数据
 *       - 轻量级设计，避免复杂计算（体积/面积计算留空）
 * 
 * @note 与后续模块对接：
 *       - 形函数模块：使用get_element_nodes()获取节点坐标
 *       - 单元矩阵模块：使用区域判断进行材料属性查找
 *       - 后处理模块：使用查询结果进行可视化数据处理
 */
class MeshQuery {
public:
    /**
     * @brief 获取单元的节点坐标列表
     * @param element 单元对象
     * @param all_nodes 所有节点的列表
     * @return std::vector<Node> 该单元包含的所有节点（按node_ids顺序）
     * 
     * @note 使用示例：
     *       @code
     *       auto elem_nodes = MeshQuery::get_element_nodes(elem, mesh_data.nodes);
     *       for (const auto& node : elem_nodes) {
     *           // 访问 node.x, node.y, node.z
     *       }
     *       @endcode
     */
    static std::vector<Node> get_element_nodes(
        const Element& element, 
        const std::vector<Node>& all_nodes
    );
    
    /**
     * @brief 判断节点是否在某个区域内
     * @param node_id 节点全局ID
     * @param region_id 区域ID
     * @param all_nodes 所有节点的列表
     * @return bool 如果节点属于该区域返回true
     */
    static bool is_node_in_region(
        int node_id, 
        int region_id, 
        const std::vector<Node>& all_nodes
    );
    
    /**
     * @brief 判断单元是否在某个区域内
     * @param element 单元对象
     * @param region_id 区域ID
     * @return bool 如果单元属于该区域返回true（通过element.region_id判断）
     */
    static bool is_element_in_region(
        const Element& element, 
        int region_id
    );
    
    /**
     * @brief 获取某个区域内的所有单元ID
     * @param region_id 区域ID
     * @param all_elements 所有单元的列表
     * @return std::vector<int> 属于该区域的单元ID列表
     */
    static std::vector<int> get_elements_in_region(
        int region_id, 
        const std::vector<Element>& all_elements
    );
    
    /**
     * @brief 获取某个区域内的所有节点ID
     * @param region_id 区域ID
     * @param all_nodes 所有节点的列表
     * @return std::vector<int> 属于该区域的节点ID列表
     */
    static std::vector<int> get_nodes_in_region(
        int region_id, 
        const std::vector<Node>& all_nodes
    );

private:
    // 禁止实例化
    MeshQuery() = delete;
};

} // namespace fe_em
