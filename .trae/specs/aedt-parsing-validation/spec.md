# AEDT文件解析与数据验证 Spec

## Why

当前项目已具备基础的Maxwell AEDT文件解析框架（`MaxwellParser` + `MaxwellParserImpl` + `MaxwellConverterImpl`），但`MaxwellParserImpl`中的边界条件、激励源、求解设置、几何数据提取方法均为TODO空实现。需要对`docs/project/Temp.aedt`（42029行，Maxwell 2D瞬态电机仿真文件）实现完整的数据解析，并验证读取数据的正确性，为后续FETI-DP求解器提供可靠的数据输入。

## What Changes

- **完善`MaxwellParserImpl`的数据提取方法**：实现`extractBoundaries()`、`extractExcitations()`、`extractSolutionSetup()`、`extractGeometry()`的完整逻辑
- **增强材料数据提取**：完善`extractMaterials()`，正确处理非线性B-H曲线、磁芯损耗系数、温度相关属性等
- **修复`MaxwellParser`的解析健壮性**：处理AEDT文件中的特殊格式（如带空格的属性名、多行属性值、特殊函数调用格式等）
- **编写AEDT文件解析验证测试**：创建`test_aedt_parsing.cpp`，使用Google Test对Temp.aedt的解析结果进行系统性验证
- **更新CMakeLists.txt**：添加新测试目标的构建配置

## Impact

- Affected specs: 无（新功能开发）
- Affected code:
  - `include/tool/maxwell_parser.hpp` - 可能需要增强Value类型定义
  - `src/tool/maxwell_parser.cpp` - 增强解析器对特殊格式的处理
  - `src/tool/maxwell_parser_impl.cpp` - 核心修改：实现所有数据提取方法
  - `test/test_aedt_parsing.cpp` - 新增：AEDT解析验证测试
  - `test/CMakeLists.txt` - 新增测试目标

## ADDED Requirements

### Requirement: 完整AEDT文件解析

系统SHALL能够完整解析`docs/project/Temp.aedt`文件（Maxwell 2D瞬态仿真），提取以下数据并转换为内部数据模型：

1. **项目元数据**：文件版本、产品信息、设计名称、求解类型（Transient）、几何模式（XY）
2. **材料数据**：所有定义的材料（B27AV1400硅钢片、N45UH永磁体、copper铜、vacuum真空、Arnold_Magnetics_N42M永磁体、JFE_35JN200硅钢片等），包括：
   - 线性/非线性磁导率
   - B-H曲线数据点（含H-B坐标对）
   - 电导率
   - 磁芯损耗参数（kh, kc, ke）
   - 矫顽力矢量
   - 质量密度
3. **边界条件数据**：所有边界定义（VectorPotential1矢量势边界、各相线圈边界PhA/PhB/PhC及其正负极性），包括：
   - 边界类型（Vector Potential / Coil）
   - 关联的边/对象ID列表
   - 线圈参数（导体数、绕组号、极性）
4. **激励源/线圈数据**：从BoundarySetup中提取的线圈配置
5. **求解设置数据**：Setup1瞬态求解配置，包括：
   - 求解类型（Transient）
   - 时间步长（0.0001s）
   - 停止时间（0.025s）
   - 非线性残差（0.0001）
   - 自适应时间步设置
6. **网格设置数据**：网格操作配置（长度-based网格、曲面近似等）

#### Scenario: 解析Temp.aedt文件成功

- **WHEN** 调用`MaxwellParserImpl::parseAllData()`解析`docs/project/Temp.aedt`
- **THEN** 返回完整的JSON数据结构，包含非空的materials、boundaries、excitations、solution_setup字段
- **AND** 材料数量 >= 8（至少包含上述列出的主要材料）
- **AND** 边界条件数量 >= 20（包含矢量势边界和各相线圈边界）
- **AND** 求解设置中SetupType为"Transient"

### Requirement: 材料数据正确性验证

系统SHALL能够验证从AEDT文件提取的材料数据与文件原始值一致：

#### Scenario: 验证关键材料属性

- **WHEN** 提取材料`copper`的数据
- **THEN** conductivity = 58000000 (S/m)
- **AND** permeability ≈ 0.999991 (相对磁导率)
- **AND** mass_density = 8933 (kg/m³)

- **WHEN** 提取材料`vacuum`的数据
- **THEN** permittivity = 1
- **AND** LibLocation = "SysLibrary"

- **WHEN** 提取材料`B27AV1400_2DSF0.970`的非线性B-H曲线
- **THEN** B-H曲线包含108个数据点
- **AND** 第一个数据点为 H=0, B=0
- **AND** 最后一个数据点的B值约等于2.086 T

- **WHEN** 提取材料`N45UH_100C_1.251_2DSF1.000_X`的永磁体参数
- **THEN** permeability = 1.0739096
- **AND** 磁性矫顽力Magnitude = -927000 A/m
- **AND** 方向为DirComp1=1（X方向）

### Requirement: 边界条件数据正确性验证

系统SHALL能够验证边界条件数据正确提取：

#### Scenario: 验证矢量势边界

- **WHEN** 提取名为`VectorPotential1`的边界条件
- **THEN** BoundType = "Vector Potential"
- **AND** Value = "0"（狄利克雷边界条件，Az=0）
- **AND** 关联18条边（Edges(18)）

#### Scenario: 验证线圈边界

- **WHEN** 提取名为`PhA_0`的边界条件
- **THEN** BoundType = "Coil"
- **AND** Winding = 1（A相绕组）
- **AND** PolarityType = "Positive"
- **AND** ParentBndID = 1（父边界ID）

### Requirement: 求解设置数据正确性验证

系统SHALL能够验证求解设置参数正确：

#### Scenario: 验证瞬态求解设置

- **WHEN** 提取Setup1的求解配置
- **THEN** SetupType = "Transient"
- **AND** StopTime = "0.025s"
- **AND** TimeStep = "0.0001s"
- **AND** NonlinearSolverResidual = "0.0001"
- **AND** FrequencyOfAddedVoltageSource = "308.06667Hz"

### Requirement: 项目元数据正确性验证

#### Scenario: 验证项目基本信息

- **WHEN** 提取项目元数据
- **THEN** Product = "ElectronicsDesktop"
- **AND** Desktop Version = (2024, 1)
- **AND** 设计名称 = "quan_ID0x45deg_20C3_temp"
- **AND** SolutionType = "Transient"
- **AND** GeometryMode = "XY"
- **AND** ModelDepth = "140mm"

## MODIFIED Requirements

### Requirement: MaxwellParser解析健壮性增强

现有的`MaxwellParser`需要增强以处理AEDT文件中的特殊格式：

1. **带空格/特殊字符的属性名**：如`'Conductor number'`、`'Perform Minimal validation'`、`'N Steps'`、`'Plane Background'`等
2. **函数调用形式属性**：如`Version(2024, 1)`、`Temperatures(6, '22cel', 11, ...)`、`Edges(18)`、`Objects(2331)`等
3. **Name()形式属性**：如`Name('Bt1')`、`ID=0`等
4. **长数组属性**：如Temperatures包含大量键值对
5. **属性值中嵌套引号**：如`Magnitude='-927000A_per_meter'`

现有PROPERTY_PATTERN正则`(\w+)\s*=\s*(.*)`需扩展为支持属性名中包含空格和特殊字符的情况。
