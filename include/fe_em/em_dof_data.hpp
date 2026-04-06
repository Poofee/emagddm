/**
 * @file em_dof_data.hpp
 * @brief 电磁物理层 - 自由度(DOF)映射数据结构定义
 * @details 定义单元局部自由度到全局自由度的映射关系，是有限元组装与求解的核心数据结构。
 *          本模块存储每个单元的DOF映射信息，包括标量节点DOF和矢量棱边DOF的分配。
 *          通过Local2Global映射，单元矩阵可正确组装到全局系统矩阵的对应位置。
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 *
 * @note 核心概念说明：
 *       - 局部DOF: 单元内部的自由度编号（从0开始连续编号）
 *       - 全局DOF: 整个离散系统的自由度编号（经约束处理后重新编号）
 *       - 约束DOF: 受边界条件（如Dirichlet）约束的自由度，映射值为-1
 *       - 自由DOF: 未受约束的自由度，映射值>=0表示其在全局系统中的编号
 *
 * @note 与后续模块对接：
 *       - DOF管理模块: 构建并填充Local2Global映射表
 *       - 单元组装模块: 根据indices将单元矩阵 scatter 到全局矩阵
 *       - 求解器模块: 利用映射关系提取子向量、施加边界条件
 */

#pragma once

#include <string>
#include <vector>
#include <tuple>
#include "em_mesh_data.hpp"  // 获取ElemType类型别名（fe_em::ElemType = numeric::ElementType）

namespace fe_em {

// ==================== 局部-全局DOF映射结构体 ====================

/**
 * @struct Local2Global
 * @brief 单元局部自由度到全局自由度的映射数据结构
 * @details 存储单个单元的所有局部DOF到全局DOF的映射关系，
 *          是有限元离散系统组装与求解的核心桥梁。
 *
 * @note 映射规则：
 *       - indices[i] == -1: 第i个局部DOF为约束DOF（如Dirichlet边界），不参与求解
 *       - indices[i] >= 0:  第i个局部DOF为自由DOF，值为该DOF在全局系统中的编号
 *
 * @note DOF排列顺序（MIXED_AV场景）：
 *       - 前num_scalar_dofs个: 标量节点DOF（电位V或磁位ψ）
 *       - 后num_vector_dofs个: 矢量棱边DOF（磁矢势A的切向分量）
 *
 * @note 典型使用流程：
 *       @code
 *       // 1. DOF管理器构建映射
 *       Local2Global l2g(elem.num_dofs);
 *       l2g.element_id = elem.id;
 *       l2g.elem_type = elem.elem_type;
 *       // 填充indices...
 *
 *       // 2. 组装时使用映射scatter单元矩阵
 *       for (int i = 0; i < local_size; ++i) {
 *           if (l2g.indices[i] >= 0) {  // 仅处理自由DOF
 *               int global_i = l2g.indices[i];
 *               // 将local_matrix[i][j]累加到global_matrix[global_i][global_j]
 *           }
 *       }
 *       @endcode
 */
struct Local2Global {
    std::vector<int> indices;        ///< 核心映射向量（-1=约束DOF, >=0=自由DOF全局编号）
    int element_id = -1;             ///< 所属单元ID（未设置时为-1）
    int num_scalar_dofs = 0;         ///< 标量节点DOF数量（用于SCALAR_ONLY或MIXED_AV的前半段）
    int num_vector_dofs = 0;         ///< 矢量棱边DOF数量（用于VECTOR_EDGE_ONLY或MIXED_AV的后半段）
    ElemType elem_type;              ///< 单元类型枚举（复用numeric::ElementType）

    /**
     * @brief 默认构造函数
     * @details 初始化所有成员为默认值，indices为空向量
     */
    Local2Global() = default;

    /**
     * @brief 指定大小的构造函数
     * @param size 映射向量的初始大小（即该单元的总局部DOF数）
     * @details 将indices初始化为size个-1，表示所有DOF初始状态均为约束DOF
     *          （后续由DOF管理器根据实际约束情况更新为>=0的自由DOF编号）
     */
    explicit Local2Global(size_t size)
        : indices(size, -1) {}

    /**
     * @brief 判断是否为混合单元（同时包含标量和矢量DOF）
     * @return true 如果标量DOF数量和矢量DOF数量均大于0
     * @return false 如果仅含一种类型的DOF或两种都为0
     * @note MIXED_AV格式（A-V涡流场）的单元返回true，
     *       SCALAR_ONLY和VECTOR_EDGE_ONLY返回false
     */
    bool is_mixed() const {
        return num_scalar_dofs > 0 && num_vector_dofs > 0;
    }

    /**
     * @brief 提取标量DOF的局部-全局映射
     * @return std::vector<int> 前num_scalar_dofs个映射值的副本
     * @details 从indices向量的起始位置提取num_scalar_dofs个元素，
     *          对应单元中所有标量节点DOF的全局映射
     * @note 若num_scalar_dofs为0，返回空向量；若num_scalar_dofs超过indices大小则截断
     */
    std::vector<int> get_scalar_indices() const {
        if (num_scalar_dofs <= 0) {
            return {};
        }
        size_t count = static_cast<size_t>(num_scalar_dofs);
        count = std::min(count, indices.size());
        return std::vector<int>(indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(count));
    }

    /**
     * @brief 提取矢量DOF的局部-全局映射
     * @return std::vector<int> 后num_vector_dofs个映射值的副本
     * @details 对于MIXED_AV单元，从indices的第num_scalar_dofs位置开始提取；
     *          对于纯VECTOR_EDGE_ONLY单元，从索引0开始提取
     * @note 若num_vector_dofs为0，返回空向量；若超出范围则截断
     */
    std::vector<int> get_vector_indices() const {
        if (num_vector_dofs <= 0) {
            return {};
        }
        size_t start = is_mixed()
            ? static_cast<size_t>(num_scalar_dofs)
            : 0;
        size_t count = static_cast<size_t>(num_vector_dofs);
        start = std::min(start, indices.size());
        count = std::min(count, indices.size() - start);
        return std::vector<int>(
            indices.begin() + static_cast<std::ptrdiff_t>(start),
            indices.begin() + static_cast<std::ptrdiff_t>(start + count)
        );
    }
};

} // namespace fe_em
