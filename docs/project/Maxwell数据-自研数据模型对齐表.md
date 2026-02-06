# Maxwell数据-自研数据模型对齐表

## 1. 概述

本文档详细说明Maxwell Project文件中的数据与自研数据模型的映射关系，确保数据读取的准确性和完整性。

## 2. 项目元数据映射

### 2.1 Maxwell项目元数据 → 自研ProjectMeta

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `ProjectInfo.Version` | `maxwell_version` | `MaxwellVersion` | 版本号映射 | 支持R15-R24及更新版本 |
| `ProjectInfo.SimulationType` | `simulation_type` | `SimulationType` | 枚举映射 | 静磁/瞬态/涡流等 |
| `ProjectInfo.Dimension` | `dimension` | `DimType` | 维度映射 | 2D/3D/轴对称 |
| `ProjectInfo.Name` | `project_name` | `std::string` | 直接映射 | 项目名称 |
| `ProjectInfo.Description` | `description` | `std::string` | 直接映射 | 项目描述 |
| `ProjectInfo.CreatedDate` | `created_date` | `std::string` | 日期格式转换 | ISO 8601格式 |
| `ProjectInfo.ModifiedDate` | `modified_date` | `std::string` | 日期格式转换 | ISO 8601格式 |

## 3. 几何数据映射

### 3.1 Maxwell几何数据 → 自研GeometryData

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `Model.Vertices` | `vertices` | `std::vector<Vertex>` | 坐标转换 | 单位：m |
| `Model.Edges` | `edges` | `std::vector<Edge>` | 拓扑关系映射 | 边-顶点关联 |
| `Model.Faces` | `faces` | `std::vector<Face>` | 拓扑关系映射 | 面-边关联 |
| `Model.Solids` | `solids` | `std::vector<Solid>` | 拓扑关系映射 | 体-面关联 |
| `Model.Transform` | `transforms` | `std::vector<Transform>` | 矩阵转换 | 变换矩阵 |
| `Model.UDPGeometries` | `udp_geometries` | `std::vector<UDPData>` | 自定义几何映射 | 用户定义原语 |

## 4. 材料数据映射

### 4.1 Maxwell材料数据 → 自研MaterialData

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `MaterialLibrary.Material.Name` | `name` | `std::string` | 直接映射 | 材料名称 |
| `MaterialLibrary.Material.Type` | `type` | `MatType` | 枚举映射 | 材料类型 |
| `MaterialLibrary.Material.Conductivity` | `conductivity` | `double` | 单位转换 | S/m |
| `MaterialLibrary.Material.RelativePermeability` | `relative_permeability` | `double` | 直接映射 | 相对磁导率 |
| `MaterialLibrary.Material.BHCurve` | `bh_curve` | `std::vector<BHPoint>` | 点对点映射 | B-H曲线数据 |
| `MaterialLibrary.Material.CoreLossModel` | `coreloss_model` | `CoreLossModelType` | 枚举映射 | 铁损模型类型 |
| `MaterialLibrary.Material.CoreLossCoefficients` | `coreloss_coefficients` | `CoreLossCoefficients` | 系数映射 | Steinmetz系数等 |
| `MaterialLibrary.Material.TemperatureDependence` | `temp_dependence` | `TemperatureDependence` | 温度相关参数 | 温度系数 |

### 4.2 磁芯损耗模型映射

| Maxwell铁损模型 | 自研铁损模型 | 参数映射 | 备注 |
|----------------|--------------|----------|------|
| `Steinmetz` | `STEINMETZ` | ks, alpha, beta | 经典Steinmetz方程 |
| `GeneralizedSteinmetz` | `GENERALIZED_STEINMETZ` | ks, alpha, beta, kn | 广义Steinmetz方程 |
| `IorioVaraldi` | `IORIO_VARALDI` | 特定参数集 | Iorio-Varaldi模型 |
| `Preisach` | `PREISACH` | 磁滞参数 | Preisach磁滞模型 |
| `Bertotti` | `BERTOTTI` | 分离损耗参数 | Bertotti分离铁损模型 |

## 5. 边界条件映射

### 5.1 Maxwell边界条件 → 自研BoundaryData

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `BoundaryConditions.Boundary.Name` | `name` | `std::string` | 直接映射 | 边界名称 |
| `BoundaryConditions.Boundary.Type` | `type` | `BndType` | 枚举映射 | 边界类型 |
| `BoundaryConditions.Boundary.Value` | `parameters` | `std::vector<double>` | 参数解析 | 边界参数值 |
| `BoundaryConditions.Boundary.GeometryIDs` | `geometry_ids` | `std::vector<int>` | ID映射 | 关联几何实体 |
| `BoundaryConditions.Boundary.SubType` | `subdivision` | `BoundarySubdivision` | 子类型映射 | 边界细分参数 |

### 5.2 Maxwell特殊边界映射

| Maxwell边界类型 | 自研边界类型 | 特殊参数 | 备注 |
|----------------|--------------|----------|------|
| `PerfectE` | `PERFECT_E` | 电场边界 | 理想电边界 |
| `PerfectH` | `PERFECT_H` | 磁场边界 | 理想磁边界 |
| `Radiation` | `RADIATION` | 辐射参数 | 辐射边界 |
| `InfiniteSphere` | `INFINITE_BOX` | 无限远参数 | 无限远边界 |
| `MasterSlave` | `MASTER_SLAVE` | 主从映射 | 周期性边界 |

## 6. 激励源数据映射

### 6.1 Maxwell激励源 → 自研ExcitationData

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `Excitations.Excitation.Name` | `name` | `std::string` | 直接映射 | 激励名称 |
| `Excitations.Excitation.Type` | `type` | `ExcitationType` | 枚举映射 | 激励类型 |
| `Excitations.Excitation.Amplitude` | `amplitude` | `double` | 直接映射 | 激励幅值 |
| `Excitations.Excitation.Frequency` | `frequency` | `double` | 直接映射 | 激励频率 |
| `Excitations.Excitation.Waveform` | `waveform` | `WaveformType` | 枚举映射 | 波形类型 |
| `Excitations.Excitation.Coil` | `coil` | `CoilDefinition` | 线圈定义映射 | 线圈参数 |
| `Excitations.Excitation.Motion` | `motion` | `MotionParameters` | 运动参数映射 | 运动激励 |

### 6.2 时变激励波形映射

| Maxwell波形类型 | 自研波形类型 | 参数映射 | 备注 |
|----------------|--------------|----------|------|
| `DC` | `DC` | 幅值 | 直流波形 |
| `Sinusoidal` | `SINUSOIDAL` | 幅值、频率、相位 | 正弦波形 |
| `Square` | `SQUARE` | 幅值、频率、占空比 | 方波波形 |
| `Triangular` | `TRIANGULAR` | 幅值、频率 | 三角波形 |
| `Custom` | `CUSTOM` | 自定义采样点 | 自定义波形 |

## 7. 网格数据映射

### 7.1 Maxwell网格数据 → 自研MeshData

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `MeshOperations.Mesh.Type` | `generation_type` | `MeshGenerationType` | 枚举映射 | 网格生成类型 |
| `MeshOperations.Mesh.MaxElementSize` | `max_element_size` | `double` | 直接映射 | 最大单元尺寸 |
| `MeshOperations.Mesh.MinElementSize` | `min_element_size` | `double` | 直接映射 | 最小单元尺寸 |
| `MeshOperations.Mesh.BoundaryLayer` | `boundary_layer` | `bool` | 直接映射 | 边界层加密 |
| `MeshOperations.Mesh.SkinDepth` | `skin_depth_params` | `SkinDepthParams` | 趋肤深度参数 | 涡流网格加密 |
| `MeshOperations.Mesh.CoreLoss` | `coreloss_refinement` | `CoreLossRefinement` | 铁损加密参数 | 铁损区域加密 |

## 8. 求解设置映射

### 8.1 Maxwell求解设置 → 自研SolverSettings

| Maxwell字段 | 自研字段 | 数据类型 | 转换规则 | 备注 |
|------------|----------|----------|----------|------|
| `AnalysisSetup.Solver.Type` | `solver_type` | `SolverType` | 枚举映射 | 求解器类型 |
| `AnalysisSetup.Solver.Convergence` | `convergence_type` | `ConvergenceType` | 枚举映射 | 收敛准则 |
| `AnalysisSetup.Solver.MaxIterations` | `maximum_iterations` | `int` | 直接映射 | 最大迭代次数 |
| `AnalysisSetup.HPC.Mode` | `hpc_mode` | `HPCParallelMode` | 枚举映射 | HPC并行模式 |
| `AnalysisSetup.HPC.NumCores` | `num_cores` | `int` | 直接映射 | 核心数量 |
| `AnalysisSetup.HPC.SolverMode` | `hpc_solver_mode` | `HPCSolverMode` | 枚举映射 | HPC求解模式 |

## 9. 单位系统映射

### 9.1 Maxwell单位 → 自研单位系统

| 物理量 | Maxwell单位 | 自研单位 | 转换系数 | 备注 |
|--------|-------------|----------|----------|------|
| 长度 | m | m | 1.0 | 米 |
| 电导率 | S/m | S/m | 1.0 | 西门子/米 |
| 磁导率 | H/m | H/m | 1.0 | 亨利/米 |
| 频率 | Hz | Hz | 1.0 | 赫兹 |
| 电流 | A | A | 1.0 | 安培 |
| 电压 | V | V | 1.0 | 伏特 |
| 磁通密度 | T | T | 1.0 | 特斯拉 |
| 磁场强度 | A/m | A/m | 1.0 | 安培/米 |

## 10. 数据验证规则

### 10.1 必填字段验证

**项目元数据必填字段**：
- 软件版本
- 仿真类型
- 维度类型
- 项目名称

**材料数据必填字段**：
- 材料名称
- 材料类型
- 电导率（导体材料）
- 相对磁导率（磁性材料）

### 10.2 数据范围验证

**数值范围验证**：
- 电导率：≥ 0 S/m
- 相对磁导率：≥ 1
- 频率：≥ 0 Hz
- 迭代次数：≥ 1

### 10.3 数据一致性验证

**一致性规则**：
- 几何实体ID必须唯一
- 材料名称必须唯一
- 边界条件不能冲突
- 激励源参数必须合理

## 11. 异常处理映射

### 11.1 解析异常类型映射

| Maxwell解析异常 | 自研异常类型 | 处理策略 |
|----------------|--------------|----------|
| 文件不存在 | `FileNotFoundException` | 提示用户检查文件路径 |
| 文件格式错误 | `InvalidFormatException` | 提供格式修复建议 |
| 版本不兼容 | `VersionCompatibilityException` | 提供版本转换工具 |
| 数据缺失 | `MissingFieldException` | 使用默认值或提示用户 |
| 数据范围异常 | `DataRangeException` | 数据修正或拒绝导入 |

## 12. 测试用例设计

### 12.1 标准测试用例

**测试文件**：
- Maxwell标准2D静磁算例
- Maxwell标准3D涡流算例
- Maxwell标准轴对称算例

**验证指标**：
- 数据读取准确率 ≥ 99%
- 核心数据无偏差
- 单位转换准确
- 映射关系正确

### 12.2 边界测试用例

**测试场景**：
- 空项目文件
- 超大项目文件（100万+节点）
- 多版本兼容性测试
- 异常数据格式测试

---

**文档版本**：1.0  
**创建日期**：2026-02-06  
**更新记录**：初始版本创建