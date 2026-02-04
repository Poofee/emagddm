/**
 * @file fe_em_mesh.hpp
 * @brief 电磁物理层 - 电磁网格管理模块头文件
 * @details 基于Gmsh API实现2D/3D网格的读取、解析、子域拆分、界面DOF提取，适配低频电磁A-φ法的网格特性
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include <vector>
#include <string>

namespace fe_em {

/**
 * @class EMMesh
 * @brief 电磁网格管理类，负责网格的读取、解析和子域划分
 */
class EMMesh {
public:
    /**
     * @brief 构造函数
     */
    EMMesh();

    /**
     * @brief 析构函数
     */
    ~EMMesh();

    /**
     * @brief 从Gmsh文件加载网格
     * @param mesh_file 网格文件路径
     * @return bool 加载成功返回true，失败返回false
     */
    bool loadFromGmsh(const std::string& mesh_file);

    /**
     * @brief 获取网格节点数
     * @return int 节点数量
     */
    int getNodeCount() const;

    /**
     * @brief 获取单元数
     * @return int 单元数量
     */
    int getElementCount() const;

private:
    int node_count_;
    int element_count_;
};

} // namespace fe_em