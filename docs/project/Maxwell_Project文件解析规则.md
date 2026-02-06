# Maxwell Project文件解析规则

## 1. 概述

本文档详细说明Maxwell Project文件（.aedt/.aedtz/.xml/.amat）的解析规则，为电磁模块数据读取模块开发提供技术指导。

## 2. 文件格式分析

### 2.1 .aedt/.aedtz文件格式

**文件结构**：
- 项目根目录：包含项目配置文件、几何文件、材料文件等
- 主要文件：
  - `Project.aedt`：项目主配置文件
  - `Project.aedt.auto`：自动保存文件
  - `Project.aedtz`：压缩项目包

**核心节点结构**：
```xml
<!-- 示例结构，具体节点名称需根据实际文件分析 -->
<Project>
    <Metadata>
        <Version>2024.1</Version>
        <SimulationType>Magnetostatic</SimulationType>
        <Dimension>2D</Dimension>
    </Metadata>
    <Geometry>
        <!-- 几何模型定义 -->
    </Geometry>
    <Materials>
        <!-- 材料定义 -->
    </Materials>
    <Boundaries>
        <!-- 边界条件定义 -->
    </Boundaries>
    <Excitations>
        <!-- 激励源定义 -->
    </Excitations>
    <Mesh>
        <!-- 网格设置 -->
    </Mesh>
    <Solution>
        <!-- 求解设置 -->
    </Solution>
</Project>
```

### 2.2 XML导出格式

Maxwell支持将项目导出为XML格式，包含完整项目数据：

**主要节点**：
- `ProjectInfo`：项目基本信息
- `Model`：几何模型数据
- `MaterialLibrary`：材料库数据
- `BoundaryConditions`：边界条件数据
- `Excitations`：激励源数据
- `MeshOperations`：网格操作数据
- `AnalysisSetup`：分析设置数据

### 2.3 .amat材料库文件格式

**文件结构**：
```xml
<MaterialLibrary>
    <Material Name="Copper">
        <Property Name="conductivity" Value="5.8e7" Unit="S/m"/>
        <Property Name="relative_permeability" Value="1"/>
        <BH_Curve>
            <Point H="0" B="0"/>
            <Point H="100" B="1.2"/>
            <!-- ... -->
        </BH_Curve>
    </Material>
    <!-- 更多材料定义 -->
</MaterialLibrary>
```

## 3. 核心数据解析规则

### 3.1 项目元数据解析

**解析目标**：提取项目基本信息
- 软件版本信息
- 仿真类型（静磁/瞬态/涡流等）
- 维度类型（2D/3D/轴对称）
- 项目名称和描述

**关键字段**：
```cpp
struct ProjectMeta {
    std::string maxwell_version;    // Maxwell软件版本
    std::string simulation_type;    // 仿真类型
    std::string dimension;          // 维度
    std::string project_name;       // 项目名称
    std::string description;        // 项目描述
    std::string created_date;       // 创建日期
    std::string modified_date;      // 修改日期
};
```

### 3.2 几何数据解析

**解析目标**：提取几何模型拓扑关系
- 顶点坐标数据
- 边、面、体拓扑关系
- 几何变换矩阵
- 用户自定义原语（UDP）数据

**关键字段**：
```cpp
struct GeometryData {
    std::vector<Vertex> vertices;           // 顶点列表
    std::vector<Edge> edges;               // 边列表
    std::vector<Face> faces;               // 面列表
    std::vector<Solid> solids;             // 体列表
    std::vector<Transform> transforms;     // 变换矩阵
    std::vector<UDPData> udp_geometries;   // UDP几何数据
};
```

### 3.3 材料数据解析

**解析目标**：提取材料属性参数
- 线性/非线性材料参数
- 各向同性/各向异性参数
- B-H曲线数据
- 磁芯损耗模型参数
- 温度相关参数

**关键字段**：
```cpp
struct MaterialData {
    std::string name;                      // 材料名称
    MaterialType type;                     // 材料类型
    double conductivity;                   // 电导率
    double relative_permeability;          // 相对磁导率
    std::vector<BHPoint> bh_curve;        // B-H曲线数据
    CoreLossModel coreloss_model;         // 磁芯损耗模型
    TemperatureDependence temp_dependence; // 温度相关参数
};
```

### 3.4 边界条件解析

**解析目标**：提取边界条件设置
- 边界类型（Dirichlet/Neumann/Periodic等）
- 边界参数值
- 关联几何实体
- 边界细分参数

**关键字段**：
```cpp
struct BoundaryData {
    std::string name;                      // 边界名称
    BoundaryType type;                     // 边界类型
    std::vector<double> parameters;        // 边界参数
    std::vector<int> geometry_ids;        // 关联几何ID
    BoundarySubdivision subdivision;      // 边界细分参数
};
```

### 3.5 激励源数据解析

**解析目标**：提取激励源设置
- 激励类型（电流源/电压源/磁通源等）
- 线圈绕组定义
- 时变激励波形参数
- 运动相关激励参数

**关键字段**：
```cpp
struct ExcitationData {
    std::string name;                      // 激励名称
    ExcitationType type;                   // 激励类型
    double amplitude;                      // 幅值
    double frequency;                      // 频率（时变激励）
    WaveformType waveform;                 // 波形类型
    CoilDefinition coil;                   // 线圈定义
    MotionParameters motion;               // 运动参数
};
```

## 4. 数据映射规则

### 4.1 Maxwell数据到自研数据模型映射

| Maxwell数据字段 | 自研数据模型字段 | 转换规则 |
|----------------|------------------|----------|
| `Material.conductivity` | `MaterialData.conductivity` | 直接映射，单位转换 |
| `Material.BH_Curve` | `MaterialData.bh_curve` | 点对点映射，插值处理 |
| `Boundary.value` | `BoundaryData.parameters` | 参数解析，单位转换 |
| `Excitation.amplitude` | `ExcitationData.amplitude` | 直接映射，单位转换 |
| `Geometry.vertex` | `GeometryData.vertices` | 坐标转换，单位统一 |

### 4.2 单位系统转换

Maxwell使用国际单位制（SI），需要确保单位转换准确：
- 长度：m → m（保持不变）
- 电导率：S/m → S/m（保持不变）
- 磁导率：H/m → H/m（保持不变）
- 频率：Hz → Hz（保持不变）

## 5. 异常处理规则

### 5.1 文件格式异常

**处理规则**：
- 文件不存在：抛出`FileNotFoundException`
- 文件格式错误：抛出`InvalidFormatException`
- 文件损坏：抛出`CorruptedFileException`

### 5.2 数据解析异常

**处理规则**：
- 数据格式不匹配：抛出`DataFormatException`
- 必填字段缺失：抛出`MissingFieldException`
- 数据范围异常：抛出`DataRangeException`

### 5.3 版本兼容性异常

**处理规则**：
- 版本不兼容：抛出`VersionCompatibilityException`
- 提供降级处理建议

## 6. 测试用例设计

### 6.1 标准测试用例

**测试文件**：
- Maxwell标准2D静磁算例
- Maxwell标准3D涡流算例
- Maxwell标准轴对称算例

**验证指标**：
- 数据读取准确率 ≥ 99%
- 核心数据无偏差
- 单位转换准确

### 6.2 边界测试用例

**测试场景**：
- 空项目文件
- 损坏的项目文件
- 超大项目文件（100万+节点）
- 多版本兼容性测试

## 7. 实施计划

### 7.1 第一阶段：基础解析模块
- 实现XML格式基础解析
- 完成项目元数据读取
- 建立基本异常处理机制

### 7.2 第二阶段：核心数据解析
- 实现几何数据解析
- 实现材料数据解析
- 实现边界条件解析

### 7.3 第三阶段：高级功能
- 实现激励源数据解析
- 支持.aedt/.aedtz格式
- 优化解析性能

## 8. 参考资料

1. Maxwell官方文档：文件格式说明
2. ANSYS Maxwell用户手册
3. XML解析技术文档

---

**文档版本**：1.0  
**创建日期**：2026-02-06  
**更新记录**：初始版本创建