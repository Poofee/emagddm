# 工业级电磁场有限元网格模块 Spec

## Why
当前项目缺少完整的有限元网格**拓扑数据**结构和几何定义工具。现有的`tool::Mesh`类仅存储**网格剖分配置参数**（如max_element_size, boundary_layer等），不包含节点坐标、单元连接关系等数值计算所需的**拓扑数据**。需要构建工业级的网格数据模块作为求解器最底层数据基石。

## What Changes
- 新增 `include/fe_em/em_mesh_data.hpp` - 核心数据结构定义（Node, Element, DOFType, EMBoundaryType, EMBoundaryMarker, EMMeshData）
- 新增 `include/fe_em/element_geometry.hpp` - 三维单元局部棱边/面预定义工具类
- 新增 `src/fe_em/element_geometry.cpp` - ElementGeometry实现
- 新增 `include/fe_em/mesh_query.hpp` - 几何查询工具类
- 新增 `src/fe_em/mesh_query.cpp` - MeshQuery实现
- 新增 `include/fe_em/mesh_reader.hpp` - 网格输入虚基类
- **修改** `include/tool/project_manager.hpp` - 添加EMMeshData成员和管理接口
- **修改** `src/tool/project_manager.cpp` - 实现EMMeshData管理逻辑
- 新增 `tests/test_2d_mesh.cpp` - 二维静磁标量位测试用例
- 新增 `tests/test_3d_mesh.cpp` - 三维涡流场A-V混合测试用例
- 更新 `CMakeLists.txt` 添加新源文件编译支持

## Impact
- Affected specs: add-shape-function-module, add-em-element-integrator（为后续模块提供数据接口）
- Affected code: include/fe_em/, src/fe_em/, include/tool/project_manager.hpp, tests/

## 架构决策：混合分层方案（方案C）

### 核心理念：职责分离 + 统一管理

```
┌─────────────────────────────────────────────────────┐
│              ProjectManager (项目管理层)              │
│                                                      │
│  ┌─ 配置数据 ─┐    ┌─ 拓扑数据（新增）─┐            │
│  │ mesh_      │    │ em_mesh_data_     │            │
│  │ (MeshPtr)  │    │ (EMMeshData*)     │            │
│  │            │    │                   │            │
│  │ • max_size │    │ • nodes[]         │            │
│  │ • min_size │    │ • elements[]      │            │
│  │ • boundary │    │ • boundary_markers│            │
│  │ • layer    │    └───────┬───────────┘            │
│  └────────────┘            │                         │
│                             ↓ const pointer          │
│  ┌─ 业务数据 ─────────────────────────┐              │
│  │ materials_   → tool::Material[]    │              │
│  │ boundaries_  → tool::Boundary[]    │              │
│  └────────────────────────────────────┘              │
└─────────────────────────────────────────────────────┘
                      │
                      ↓ 数据流
┌─────────────────────────────────────────────────────┐
│           fe_em namespace (数值计算层)                │
│                                                      │
│  Node, Element        ← 核心拓扑结构                  │
│  DOFType              ← 单元DOF类型枚举              │
│  EMBoundaryMarker     ← 边界标记（轻量）             │
│  ElementGeometry      ← 几何定义工具                  │
│  MeshQuery            ← 几何查询工具                  │
│  MeshReader           ← 网格输入接口                  │
└─────────────────────────────────────────────────────┘
```

### 设计原则

#### 1️⃣ 职责分离（单一职责原则）

| 类 | 职责 | 存储内容 |
|----|------|----------|
| **tool::Mesh** | 网格**配置** | 剖分参数、加密设置、边界层 |
| **fe_em::EMMeshData** | 网格**拓扑** | 节点坐标、单元连接、DOF标记 |
| **tool::Material** | 材料**属性** | μ_r, ε_r, σ, B-H曲线 |
| **tool::Boundary** | 边界**定义** | 类型、值、几何对象关联 |

#### 2️⃣ 零拷贝引用（避免数据冗余）

```cpp
struct EMMeshData {
    std::vector<Node> nodes;                    // ✅ 存储（拓扑核心）
    std::vector<Element> elements;              // ✅ 存储（拓扑核心）
    
    // ❌ 不存储materials! 通过material_id索引ProjectManager
    // ❌ 不存储完整Boundary! 只存轻量标记
    
    std::vector<EMBoundaryMarker> boundary_markers; // ✅ 存储轻量标记
};

// 材料访问方式：
// Element.material_id → ProjectManager.getMaterial(id) → μ_r, ε_r, σ
```

#### 3️⃣ 统一管理入口（ProjectManager扩展）

```cpp
class ProjectManager {
    // ... 现有接口 ...
    
private:
    // ... 现有成员 ...
    std::optional<MeshPtr> mesh_;                    // 网格配置（已有）
    std::unique_ptr<fe_em::EMMeshData> em_mesh_data_; // 网格拓扑（新增）
    
public:
    // 新增接口：
    void setEMMeshData(std::unique_ptr<fe_em::EMMeshData> data);
    fe_em::EMMeshData* getEMMeshData() const;
    bool hasEMMeshData() const;
};
```

#### 4️⃣ 数据流场景

**场景A：从Maxwell文件加载**
```
Maxwell .aedt文件 → MaxwellParser解析
    ↓
ProjectManager.materials_ ← 材料数据
ProjectManager.boundaries_ ← 边界条件
ProjectManager.mesh_ ← 网格配置参数
ProjectManager.em_mesh_data_ ← 网格拓扑数据（节点+单元）
```

**场景B：从Gmsh/VTK文件加载**
```
Gmsh .msh文件 → MeshReaderGmsh::read()
    ↓
返回 EMMeshData（包含nodes + elements）
    ↓
ProjectManager.setEMMeshData(em_data)
材料需另外设置或从默认库获取
```

**场景C：数值计算使用**
```
求解器初始化:
    em_data = ProjectManager.getEMMeshData();
    
单元矩阵循环:
    for (elem : em_data->elements) {
        nodes = MeshQuery::get_element_nodes(elem, em_data->nodes);
        mat = ProjectManager.getMaterial(elem.material_id);
        mu = mat->getRelativePermeability();  // 通过ID索引
        // ... 计算单元矩阵 ...
    }
```

## ADDED Requirements

### Requirement: 核心拓扑数据结构（EMMeshData）
系统 SHALL 提供轻量级的有限元网格拓扑数据结构：
- **Node结构体**：三维坐标统一存储，支持区域标记
- **Element结构体**：包含节点ID列表、单元类型、DOF类型、材料ID（索引）、区域ID
- **DOFType枚举**：SCALAR_ONLY, VECTOR_EDGE_ONLY, MIXED_AV
- **EMBoundaryMarker结构体**：轻量级边界标记（节点/面/棱边ID列表 + 值）
- **EMMeshData聚合结构体**：只存拓扑数据，通过ID索引外部Material/Boundary

#### Scenario: 拓扑数据与业务数据解耦
- **WHEN** 需要访问单元的材料属性
- **THEN** 通过Element.material_id索引ProjectManager的materials_获取tool::Material对象
- **WHEN** 需要访问边界条件详细定义
- **THEN** 通过EMBoundaryMarker.id索引ProjectManager的boundaries_
- **WHEN** 仅需遍历拓扑进行数值计算
- **THEN** 直接使用EMMeshData的nodes和elements，无需访问ProjectManager

### Requirement: ProjectManager集成
系统 SHALL 扩展ProjectManager以支持EMMeshData的管理：
- 添加em_mesh_data_成员变量（unique_ptr持有所有权）
- 提供set/get/has接口
- 在项目保存/加载时序列化和反序列化EMMeshData
- 保持与现有mesh_（配置）成员的共存关系

#### Scenario: 从项目文件加载完整网格数据
- **WHEN** 调用ProjectManager::openProject()加载.aedt文件
- **THEN** 同时填充mesh_（剖分参数）和em_mesh_data_（拓扑数据）
- **WHEN** 访问网格配置信息
- **THEN** 使用getMesh()返回的Mesh对象
- **WHEN** 访问网格拓扑数据进行计算
- **THEN** 使用getEMMeshData()返回的EMMeshData对象

### Requirement: 三维单元几何定义工具
系统 SHALL 提供ElementGeometry工具类，硬编码标准有限元定义：
- 支持的单元类型：TRI3, TRI6, QUAD4, QUAD8, TET4, TET10, HEX8, HEX20, PRISM6, PRISM15, PYRAMID5, PYRAMID13
- 提供局部棱边定义（每个棱边存局部节点ID升序元组）
- 提供局部面定义（每个面存局部节点ID顺序列表）
- 提供节点数/棱边数/面数查询
- 遵循ANSYS Fluent/VTK标准有限元顺序

#### Scenario: 获取TET4单元几何信息
- **WHEN** 调用ElementGeometry::get_local_edges("TET4")
- **THEN** 返回6条棱边的局部节点ID对：(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
- **WHEN** 调用ElementGeometry::get_local_faces("TET4")
- **THEN** 返回4个面的局部节点ID列表

### Requirement: 几何查询工具
系统 SHALL 提供MeshQuery工具类，支持：
- 获取单元的节点坐标列表
- 判断节点是否在某个区域内
- 判断单元是否在某个区域内
- 所有查询操作基于传入的数据引用，不持有状态（线程安全）

#### Scenario: 查询单元信息
- **WHEN** 调用MeshQuery::get_element_nodes(element, all_nodes)
- **THEN** 返回该单元所有节点的完整坐标信息

### Requirement: 网格输入扩展接口
系统 SHALL 提供MeshReader虚基类，为未来支持MESH/VTK/ANSYS CDB等格式预留接口：
- 定义纯虚函数read()用于读取网格文件返回EMMeshData
- 允许派生类实现特定格式的解析器

#### Scenario: 扩展新的网格格式支持
- **WHEN** 需要支持新的网格格式（如VTK）
- **THEN** 继承MeshReader类并实现read()方法即可

### Requirement: 测试验证
系统 SHALL 提供两个完整的测试用例：
- 二维测试模型：构建TRI3单元EMMeshData，通过ProjectManager管理，验证数据流
- 三维测试模型：构建TET4单元EMMeshData，通过ProjectManager管理，验证A-V混合DOF类型

#### Scenario: 运行二维测试
- **WHEN** 执行test_2d_mesh.cpp
- **THEN** 创建EMMeshData包含TRI3单元，设置到ProjectManager，验证可正确获取材料和拓扑数据

#### Scenario: 运行三维测试
- **WHEN** 执行test_3d_mesh.cpp
- **THEN** 创建EMMeshData包含TET4单元，设置MIXED_AV DOF类型，验证理想导体边界标记

## MODIFIED Requirements
无

## REMOVED Requirements
无
