# 电磁场有限元自由度（DOF）管理模块 Spec

## Why
电磁场有限元求解需要统一的自由度（DOF）编号和局部-全局映射管理，作为单元矩阵装配和全局求解的核心寻址中枢。当前项目已有完整的网格数据结构（EMMeshData、Element、Node、EMBoundaryMarker）和几何工具类（ElementGeometry），但缺少DOF编号与映射功能，无法支撑后续的单元矩阵模块和全局装配模块开发。

## What Changes
- 新增核心公共数据结构头文件 `em_dof_data.hpp`（Local2Global映射表）
- 新增全局棱边ID生成器 `global_edge_id_generator.hpp/.cpp`（GlobalEdgeIDGenerator类）
- 新增DOF管理核心类 `em_dof_manager.hpp/.cpp`（EMDOFManager类）
- 新增3个测试用例验证各DOF类型场景的正确性
- 更新CMake构建配置以包含新模块

## Impact
- Affected specs: 网格模块（em_mesh_data.hpp）、几何工具（element_geometry.hpp）、数值层（shape_function_base.hpp ElementType枚举）
- Affected code: 
  - `include/fe_em/` — 新增3个头文件
  - `src/fe_em/` — 新增2个源文件
  - `test/` — 新增3个测试文件
  - `test/CMakeLists.txt` — 新增测试目标配置

## ADDED Requirements

### Requirement: 核心公共数据结构 Local2Global
系统 SHALL 提供 `Local2Global` 结构体，作为单元局部DOF索引到全局方程索引的唯一映射标准。

#### 数据成员：
- `std::vector<int> indices`：核心映射向量，-1表示约束DOF，>=0表示自由DOF全局编号
- `int element_id`：所属单元ID
- `int num_scalar_dofs`：标量节点DOF数量
- `int num_vector_dofs`：矢量棱边DOF数量
- 辅助方法：`is_mixed()`、`get_scalar_indices()`、`get_vector_indices()`

#### Scenario: MIXED_AV单元映射表拼接
- **WHEN** 创建MIXED_AV类型PRISM6单元的Local2Global映射
- **THEN** indices前半部分为标量节点DOF（按节点顺序），后半部分为矢量棱边DOF（按局部棱边顺序），num_scalar_dofs=6，num_vector_dofs=9

### Requirement: 全局棱边ID生成器 GlobalEdgeIDGenerator
系统 SHALL 提供 `GlobalEdgeIDGenerator` 类，为Nedelec棱边元生成全局唯一的棱边ID。

#### 功能：
- 输入：`const EMMeshData& mesh_data`
- 输出：
  - `int getNumGlobalEdges()`：全局棱边总数
  - `std::vector<std::vector<int>> getElemLocalToGlobalEdge()`：每个单元的「局部棱边ID→全局棱边ID」映射表
- 实现逻辑：使用`std::unordered_map<std::tuple<int,int>, int>`存储唯一棱边（升序节点ID对）并分配连续ID
- 仅处理 VECTOR_EDGE_ONLY 和 MIXED_AV 类型单元

#### Scenario: PRISM6单元全局棱边生成
- **WHEN** 1个PRISM6单元（6节点9棱边），DOFType=MIXED_AV
- **THEN** 生成9条全局棱边，每条棱边的局部→全局映射正确

### Requirement: DOF管理核心类 EMDOFManager
系统 SHALL 提供 `EMDOFManager` 类，执行完整的DOF编号流程。

#### 输入：
- `const EMMeshData& mesh_data`：网格拓扑数据
- 可选：`const std::vector<std::vector<int>>& elem_local_to_global_edge`（全局棱边ID映射表）

#### 输出：
- `int getNumFreeDOFs()`：全局自由DOF总数（最终求解规模）
- `std::vector<Local2Global> getElemLocalToGlobal()`：每个单元的Local2Global映射表
- `std::vector<double> getConstrainedDOFValues()`：约束DOF值向量

#### 编号流程（四步法）：
1. **预编号**：标量节点DOF从0开始 → 矢量棱边DOF紧接后续
2. **标记约束**：遍历边界条件，标记DIRICHLET类型的DOF为-1
3. **分块重编号**（消去法）：自由DOF在前连续编号，约束DOF标记为-1
4. **生成映射表**：根据单元DOFType分别生成SCALAR_ONLY / VECTOR_EDGE_ONLY / MIXED_AV三类映射表

#### Scenario: PRISM6 MIXED_AV完整DOF编号
- **WHEN** 1个PRISM6单元（6节点9棱边），MIXED_AV类型，1个节点Dirichlet边界，1条棱边Dirichlet边界
- **THEN** 自由DOF=13（5标量+8矢量），约束DOF=2，Local2Global映射表中2个位置为-1

### Requirement: 单元类型完备覆盖
系统 SHALL 支持以下所有单元类型的DOF管理：

#### 标量拉格朗日节点单元（SCALAR_ONLY / MIXED_AV）：
TRI3, TRI6, QUAD4, QUAD8, QUAD9, TET4, TET10, HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13

#### Nedelec一阶矢量棱边单元（VECTOR_EDGE_ONLY / MIXED_AV）：
TET4_EDGE(6棱边), HEX8_EDGE(12棱边), PRISM6_EDGE(9棱边), PYRAMID5_EDGE(8棱边)

## MODIFIED Requirements
无。本模块为全新功能，不修改现有代码逻辑。

## REMOVED Requirements
无。
