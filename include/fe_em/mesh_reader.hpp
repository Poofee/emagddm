/**
 * @file mesh_reader.hpp
 * @brief 电磁物理层 - 网格输入虚基类头文件
 * @details 定义网格文件读取器的抽象接口，为支持多种网格格式（MESH, VTK, ANSYS CDB等）预留扩展点。
 * @author Poofee
 * @date 2026-04-05
 * @version 1.0
 * 
 * @note 使用示例：
 *       @code
 *       // 从Gmsh文件读取
 *       MeshReaderGmsh reader;
 *       auto mesh_data = reader.read("model.msh");
 *       project_manager.setEMMeshData(std::move(mesh_data));
 *       @endcode
 */

#pragma once

#include "em_mesh_data.hpp"
#include <memory>
#include <string>

namespace fe_em {

/**
 * @class MeshReader
 * @brief 网格文件读取器抽象基类
 * @details 定义统一的网格文件读取接口。
 *          派生类实现特定格式的解析逻辑（如Gmsh .msh, VTK .vtk等）。
 */
class MeshReader {
public:
    /**
     * @brief 虚析构函数
     */
    virtual ~MeshReader() = default;
    
    /**
     * @brief 从文件读取网格数据
     * @param file_path 网格文件路径
     * @return std::unique_ptr<EMMeshData> 读取到的网格拓扑数据（调用者获得所有权）
     * 
     * @throws std::runtime_error 如果文件读取失败或格式错误
     * 
     * @note 返回的EMMeshData只包含拓扑数据（nodes, elements, boundary_markers），
     *       不包含材料属性。材料需另外设置或从默认库获取。
     */
    virtual std::unique_ptr<EMMeshData> read(const std::string& file_path) = 0;
    
    /**
     * @brief 检查是否支持指定格式的文件
     * @param file_path 文件路径（通过扩展名判断）
     * @return bool 如果支持该格式返回true
     */
    virtual bool supports_format(const std::string& file_path) const = 0;

protected:
    // 允许派生类构造
    MeshReader() = default;
    MeshReader(const MeshReader&) = delete;
    MeshReader& operator=(const MeshReader&) = delete;
};

} // namespace fe_em
