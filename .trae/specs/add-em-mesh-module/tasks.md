# Tasks

## 阶段一：核心拓扑数据结构定义
- [x] Task 1: 创建核心数据结构头文件 `include/fe_em/em_mesh_data.hpp`
  - [x] 1.1 定义Node结构体（节点数据：id, x, y, z, region_id）
  - [x] 1.2 定义DOFType枚举（单元DOF类型：SCALAR_ONLY, VECTOR_EDGE_ONLY, MIXED_AV）
  - [x] 1.3 定义Element结构体（单元数据：id, node_ids, type, dof_type, material_id, region_id）
  - [x] 1.4 定义EMBoundaryType枚举和EMBoundaryMarker结构体（轻量级边界标记）
  - [x] 1.5 定义EMMeshData聚合结构体（只存拓扑数据：nodes, elements, boundary_markers）
  - [x] 1.6 添加完整注释，标注与后续模块的对接点
  - [x] 1.7 添加材料参数访问辅助函数说明（通过ID索引外部Material）

## 阶段二：几何定义工具类
- [x] Task 2: 实现ElementGeometry工具类
  - [x] 2.1 创建头文件 `include/fe_em/element_geometry.hpp`，声明静态方法接口
  - [x] 2.2 实现文件 `src/fe_em/element_geometry.cpp`，硬编码所有单元类型的几何定义
  - [x] 2.3 实现一阶单元定义：TRI3(3节点,3棱,1面), QUAD4(4节点,4棱,1面), TET4(4节点,6棱,4面), HEX8(8节点,12棱,6面)
  - [x] 2.4 实现二阶单元定义：TRI6, QUAD8, TET10, HEX20
  - [x] 2.5 实现其他单元：PRISM6, PRISM15, PYRAMID5, PYRAMID13
  - [x] 2.6 标注标准来源（ANSYS Fluent/VTK顺序）

## 阶段三：几何查询工具类
- [x] Task 3: 实现MeshQuery工具类
  - [x] 3.1 创建头文件 `include/fe_em/mesh_query.hpp`
  - [x] 3.2 实现文件 `src/fe_em/mesh_query.cpp`
  - [x] 3.3 实现get_element_nodes()方法
  - [x] 3.4 实现is_node_in_region()方法
  - [x] 3.5 实现is_element_in_region()方法

## 阶段四：网格输入接口
- [x] Task 4: 创建MeshReader虚基类
  - [x] 4.1 创建头文件 `include/fe_em/mesh_reader.hpp`
  - [x] 4.2 定义纯虚接口read()返回EMMeshData

## 阶段五：ProjectManager集成
- [x] Task 5: 扩展ProjectManager支持EMMeshData管理
  - [x] 5.1 在`include/tool/project_manager.hpp`中添加em_mesh_data_成员（unique_ptr<EMMeshData>）
  - [x] 5.2 添加setEMMeshData()接口（转移所有权）
  - [x] 5.3 添加getEMMeshData()接口（返回原始指针）
  - [x] 5.4 添加hasEMMeshData()接口
  - [x] 5.5 在`src/tool/project_manager.cpp`中实现上述接口
  - [x] 5.6 确保与现有mesh_成员共存且职责清晰
  - [x] 5.7 添加注释说明两个Mesh成员的区别（配置 vs 拓扑）

## 阶段六：测试用例
- [x] Task 6: 编写二维静磁标量位测试用例 `tests/test_2d_mesh.cpp`
  - [x] 6.1 构建EMMeshData包含TRI3单元模型（3个节点）
  - [x] 6.2 使用tool::Material创建空气材料并设置到ProjectManager
  - [x] 6.3 将EMMeshData设置到ProjectManager
  - [x] 6.4 标记Dirichlet节点边界条件
  - [x] 6.5 验证从ProjectManager可正确获取拓扑数据和材料数据
  - [x] 6.6 验证ElementGeometry返回TRI3的棱边/面定义正确性

- [x] Task 7: 编写三维涡流场A-V混合测试用例 `tests/test_3d_mesh.cpp`
  - [x] 7.1 构建EMMeshData包含TET4单元模型（4个节点）
  - [x] 7.2 使用tool::Material创建导体材料并设置到ProjectManager（μ_r=1, σ=1e6）
  - [x] 7.3 将EMMeshData设置到ProjectManager
  - [x] 7.4 设置MIXED_AV DOF类型
  - [x] 7.5 标记理想导体棱边边界条件（Dirichlet Edge）
  - [x] 7.6 验证从ProjectManager可正确获取拓扑数据和材料参数
  - [x] 7.7 验证ElementGeometry返回TET4的棱边/面定义正确性

## 阶段七：构建系统集成
- [x] Task 8: 更新CMakeLists.txt并验证编译
  - [x] 8.1 在CMakeLists.txt中添加新的源文件（em_mesh_data无需cpp，element_geometry.cpp, mesh_query.cpp）
  - [x] 8.2 确保project_manager.cpp能找到fe_em头文件（添加include路径）
  - [x] 8.3 添加测试目标配置
  - [x] 8.4 验证主代码编译通过（包括ProjectManager修改）
  - [x] 8.5 运行test_2d_mesh测试用例确保通过 ✅ **全部通过**
  - [x] 8.6 运行test_3d_mesh测试用例确保通过 ✅ **全部通过**

# Task Dependencies
- [Task 2] depends on [Task 1] (ElementGeometry依赖em_mesh_data.h中的Node/Element/DOFType)
- [Task 3] depends on [Task 1] (MeshQuery依赖em_mesh_data.h中的数据结构)
- [Task 4] depends on [Task 1] (MeshReader依赖EMMeshData结构体)
- [Task 5] depends on [Task 1] (ProjectManager集成依赖EMMeshData定义)
- [Task 6] depends on [Task 1], [Task 2], [Task 3], [Task 5] (测试需核心数据+工具类+PM集成)
- [Task 7] depends on [Task 1], [Task 2], [Task 3], [Task 5] (测试需核心数据+工具类+PM集成)
- [Task 8] depends on [Task 2], [Task 3], [Task 4], [Task 5], [Task 6], [Task 7] (最后构建验证)

# ==================== 执行结果汇总 ====================
**状态**: ✅ 全部完成  
**编译**: ✅ 通过 (0 error, 0 warning)  
**测试**: 
- ✅ test_2d_mesh: 5/5 测试通过
- ✅ test_3d_mesh: 5/5 测试通过
**架构方案**: 方案C（混合分层）- EMMeshData + tool::Material 复用
