# AEDT数据结构增强与Maxwell数据完整性 Spec

## Why

当前AEDT项目加载器（`AedtProjectLoader` + `MaxwellParserImpl`）仅能提取Maxwell AEDT文件中的基础数据子集。以`docs/project/Temp.aedt`（42000+行，Maxwell 2D瞬态电机仿真）为对标样例，现有解析器在以下方面存在严重的数据缺失：

1. **材料数据不完整**：缺少温度相关B-H曲线、多频率铁损曲线（CoreLossMultiCurveData）、热修正器（ThermalModifier）、材料外观属性、堆叠类型、介电常数、力学参数等
2. **绕组数据完全缺失**：PhaseA/B/C绕组组的电流表达式、电阻、电感、并联支路数等关键电机参数未被提取
3. **运动设置完全缺失**：Band旋转运动设置（角速度、初始位置、运动轴等）未被提取
4. **网格操作不完整**：全局网格设置和局部对象级网格操作（长度-based、曲面近似、柱面间隙等）未完整映射到数据结构
5. **设计变量和输出变量缺失**：项目参数（fractions, delta, conds, R1, Le1等）和输出表达式（Lad/Laq, Flux_d/q, I_d/q等20个变量）未提取
6. **温度设置缺失**：物体-温度映射关系未提取
7. **求解设置不完整**：瞬态求解的30+个参数中仅提取了约10个

## What Changes

- **扩展Material数据结构**：增加温度相关B-H曲线、多频率铁损曲线、热修正器、材料元信息（库来源、修改时间）、外观数据、堆叠类型、附加物理属性
- **新增Winding数据结构**：专门存储绕组组信息（电流/电压表达式、电阻、电感、并联支路数、实体/绞线标志）
- **新增MotionSetup数据结构**：存储运动设置（运动类型、轴、初始位置、角速度、关联对象列表）
- **扩展Boundary数据结构**：增加ParentBndID层级关系、ConductorNumber变量引用、坐标系ID、组件标志
- **扩展SolutionSetup数据结构**：补全瞬态求解全部参数（时间积分方法、稳态检测、场保存策略、自适应时间步完整配置等）
- **新增MeshOperation数据结构**：存储网格操作详情（类型、启用状态、对象列表、长度限制、曲面近似参数）
- **新增DesignVariable数据结构**：存储设计变量（名称、值、单位、表达式）
- **新增OutputVariable数据结构**：存储输出变量表达式（名称、表达式字符串、单位）
- **新增TemperatureSettings数据结构**：存储物体-温度映射关系
- **增强MaxwellParserImpl解析能力**：实现上述所有新数据字段的提取逻辑
- **更新AedtProjectLoader转换逻辑**：将新提取的JSON数据转换为对应的C++数据结构对象

## Impact

- Affected specs: `aedt-parsing-validation`（需同步更新验证场景）
- Affected code:
  - `include/tool/project_data.hpp` - 核心修改：扩展/新增所有数据结构类
  - `include/tool/em_enums.hpp` - 新增枚举类型（运动类型细分、网格操作类型、绕组类型等）
  - `src/tool/maxwell_parser_impl.cpp` - 核心修改：实现所有新数据提取方法
  - `include/tool/maxwell_parser_impl.hpp` - 新增方法声明
  - `src/tool/aedt_project_loader.cpp` - 更新转换逻辑以处理新数据结构
  - `include/tool/aedt_project_loader.hpp` - 可能需要调整接口
  - `include/tool/project_manager.hpp` - 新增容器成员（windings_, motion_setups_, mesh_operations_, design_variables_, output_variables_, temperature_settings_）

## ADDED Requirements

### Requirement: 完整材料数据模型

系统SHALL能够保存AEDT文件中定义的全部材料属性：

#### 材料基础属性（已有，需确认完整性）
- 名称、磁导率（线性/非线性）、电导率、质量密度
- B-H曲线数据点（含normal/intrinsic类型标识）
- 磁芯损耗系数（kh, kc, ke）
- 矫顽力矢量（Magnitude + DirComp1/2/3方向分量）

#### 材料扩展属性（新增）
- 物理类型集合（PhysicsTypes: Electromagnetic, Thermal, Structural）
- 坐标系类型（CoordinateSystemType: Cartesian）
- 体/表面类型标志（BulkOrSurfaceType）
- 介电常数（permittivity）
- 热导率（thermal_conductivity）
- 比热容（specific_heat）
- 杨氏模量（youngs_modulus）
- 泊松比（poisons_ratio）
- 热膨胀系数（thermal_expansion_coefficient）
- 铁损类型选择（core_loss_type: Electrical Steel / Power Ferrite / Custom）
- 堆叠类型（stacking_type: Solid / Laminated）
- 堆叠因子（stacking_factor，如有）
- B-H曲线类型标识（BTypeForSingleCurve: normal / intrinsic）
- B-H曲线单位（HUnit, BUnit）
- 温度相关标志（IsTemperatureDependent）

#### 温度相关B-H曲线（新增）
- 多温度B-H曲线数据：`std::map<double /*温度*/, std::vector<BHDataPoint>>`
- 支持从`$begin 'Temperatures' → $begin 'Temperature'(TempRef=...) → BHCoordinates`结构提取

#### 多频率铁损曲线（新增，来自CoreLossMultiCurveData）
- 铁损单位（coreloss_unit: w_per_kg 或 w_per_m3）
- 多频率B-P曲线：`std::map<double /*频率Hz*/, std::vector<std::pair<double/*B*/, double/*P*/>>>`
- 支持从`$begin 'CoreLossMultiCurveData' → $begin 'AllCurves' → $begin 'OneCurve'(Frequency=...) → Coordinates`结构提取

#### 铁损系数设置数据（新增，来自CoefficientSetupData）
- 系数设置模式（coefficient_setup: w_per_kg / w_per_m3）
- 参考频率（Frequency）
- 叠片厚度（Thickness）
- 参考电导率（Conductivity）
- 单频率B-P曲线

#### 热修正器数据（新增，来自ModifierData > ThermalModifierData）
- 属性名称（Property: permeability / magnetic_coercivity）
- 修正公式（free_form_value: 如 "1-0.0010876*(Temp-20)"）
- 支持多属性修正器列表

#### 材料元信息（新增）
- 库来源（Library: Materials / ArnoldMagnetics / Granta Producers... / Project / PersonalLibrary / ''）
- 库位置（LibLocation: SysLibrary / Project / PersonalLibrary）
- 是否自库修改（ModSinceLib: true/false）
- 修改时间戳（ModTime: Unix时间戳）

#### 材料外观数据（新增，来自AttachedData > MatAppearanceData）
- RGB颜色值（Red, Green, Blue: 0-255）
- 透明度（Transparency: 0.0-1.0）

#### 材料备注和关键词（新增）
- 备注文本（Notes: 多行文本）
- 关键词列表（Keywords: 字符串）

#### Scenario: 完整提取B27AV1400硅钢片材料
- **WHEN** 解析包含B27AV1400材料的AEDT文件
- **THEN** Material对象应包含：
  - 非线性B-H曲线（106个数据点，BTypeForSingleCurve='normal'）
  - 铁损系数 kh=152.277, kc=0.265, ke=0
  - 铁损类型='Electrical Steel'
  - 堆叠类型='Solid'
  - 质量密度=7650 kg/m³
  - 电导率=2000000 S/m
  - **且** CoreLossMultiCurveData包含10条频率曲线（50Hz~2000Hz），每条约62-94个B-P数据点
  - **且** CoefficientSetupData包含参考频率340Hz、厚度0.35mm的78点B-P曲线

#### Scenario: 完整提取Arnold_Magnetics_N42M_20C永磁体带热修正器
- **WHEN** 解析包含Arnold_Magnetics_N42M_20C的材料
- **THEN** Material对象应包含：
  - 非线性B-H曲线（174个数据点，BTypeForSingleCurve='intrinsic'）
  - 矫顽力Magnitude=-1169650.466 A/m，方向=(1,0,0)
  - **且** ThermalModifierData包含2个修正器：
    - permeability修正: "1-0.0010876*(Temp-20)"
    - magnetic_coercivity修正: "1-0.006474*(Temp-20)"
  - **且** 包含热导率=7.6, 质量密度=7500, 比热容=460

#### Scenario: 提取温度相关B-H曲线材料（Arnold_Magnetics_N42M_20C_temptable）
- **WHEN** 解析包含温度相关B-H曲线的材料
- **THEN** Material对象应标记IsTemperatureDependent=true
  - **且** 包含基准B-H曲线（174点，20°C参考或默认）
  - **且** Temperatures map包含：
    - 20°C → 148点B-H曲线
    - 120°C → 172点B-H曲线

### Requirement: 绕组组（Winding Group）数据结构

系统SHALL新增Winding类来存储Maxwell绕组组定义：

#### 绕组组核心字段
- 名称（name: "PhaseA" / "PhaseB" / "PhaseC"）
- 边界类型（bound_type: "Winding Group"）
- 激励类型（type: "Current" / "Voltage"）
- 实体标志（is_solid: false = 绞线, true = 实心导体）
- 电流表达式（current: 可含变量的公式字符串，如 "Irms*sqrt(2)*sin(2*pi*few*time+ele_angle)"）
- 电压表达式（voltage: 如 "0"）
- 电阻（resistance: "1ohm"）
- 电感（inductance: "0H"）
- 并联支路数（parallel_branches_num: "4"）
- 关联的线圈边界列表（coils_: 引用Coil类型边界的ID/名称列表）

#### Scenario: 提取三相绕组组
- **WHEN** 解析Temp.aedt的边界条件
- **THEN** 系统应识别并创建3个Winding对象：
  - PhaseA: Current类型, 表达式="Irms*sqrt(2)*sin(2*pi*few*time+ele_angle)", 4并联支路, 关联14个Coil边界(PhA_0~13 + PhARe_0~13)
  - PhaseB: Current类型, 表达式含"-2/3*pi"相位偏移, 4并联支路, 关联14个Coil边界
  - PhaseC: Current类型, 表达式含"-4/3*pi"相位偏移, 4并联支路, 关联14个Coil边界

### Requirement: 运动设置（MotionSetup）数据结构

系统SHALL新增MotionSetup类来存储Maxwell运动设置：

#### 运动设置核心字段
- 名称（name: "MotionSetup1"）
- 运动类型（motion_type: "Band" / "Rotation" / "Translation"）
- 移动类型（move_type: "Rotate" / "Linear"）
- 坐标系ID（coordinate_system: 1）
- 运动轴（axis: "Z" / "X" / "Y"）
- 正方向标志（is_positive: true）
- 初始位置（initial_position: "7.5deg"）
- 是否有旋转限制（has_rotate_limit: false）
- 是否非柱面（non_cylindrical: false）
- 是否考虑机械瞬态（consider_mechanical_transient: false）
- 角速度/线速度（angular_velocity: "3000rpm" 或 linear_velocity）
- 关联对象列表（objects_: 对象ID列表）

#### Band/Moving区域字段（Moving子项）
- Band名称引用（band_name: 0，指向MotionSetup1）
- 运动部件对象列表（moving_objects_: 转子上的所有对象ID）

#### Scenario: 提取Band旋转运动设置
- **WHEN** 解析Temp.aedt的运动设置
- **THEN** MotionSetup对象应包含：
  - MotionType='Band', MoveType='Rotate', Axis='Z'
  - InitPos='7.5deg', AngularVelocity='3000rpm'
  - 关联对象=[6]（Band面对象）
  - **且** Moving子项包含21个转子对象ID（3406, 11, 2386, 2680, 2693, ... 2924）

### Requirement: 边界条件数据结构增强

系统SHALL扩展Boundary数据结构以支持完整的Maxwell边界定义：

#### 新增边界字段
- 父边界ID（parent_bnd_id: int, 用于Coil→Winding Group的层级关系，如PhA_0的ParentBndID=1指向PhaseA）
- 组件标志（is_component: bool）
- 坐标系ID（coordinate_system: int, 默认-1）
- 导体数（conductor_number: string, 可为变量引用如"conds"而非纯数字）

#### Scenario: 提取完整线圈边界层级
- **WHEN** 解析Temp.aedt的VectorPotential1边界
- **THEN** Boundary对象应包含：
  - BoundType='Vector Potential', Value='0'
  - Edges=[18条边]
  - ParentBndID=-1（无父边界，顶层边界）
  - CoordinateSystem=-1

- **WHEN** 解析PhA_0线圈边界
- **THEN** Boundary对象应包含：
  - BoundType='Coil', PolarityType='Positive'
  - Objects=[1857], Winding=1, ParentBndID=1（指向PhaseA）
  - ConductorNumber='conds'（变量引用）

### Requirement: 求解设置数据结构增强

系统SHALL扩展SolutionSetup以完整存储Maxwell瞬态求解配置：

#### 新增瞬态求解字段
- 启用标志（enabled: bool）
- 网格导入（mesh_link_import: bool）
- 时间积分方法（time_integration_method: int, 0=后向欧拉）
- B-H曲线平滑（smooth_bh_curve: bool）
- 输出误差（output_error: bool）
- 按对象输出铁损（output_per_object_core_loss: bool）
- 按对象输出实心损耗（output_per_object_solid_loss: bool）
- 使用控制程序（use_control_program: bool）
- 控制程序名/参数（control_program_name/arg: string）
- 快速到达稳态（fast_reach_steady_state: bool）
- 自动检测稳态（auto_detect_steady_state: bool）
- 添加电压源频率（frequency_of_added_voltage_source: string）
- 停止准则（stop_criterion: double, 如0.005）
- 是否一般瞬态（is_general_transient: bool）
- 是否半周期瞬态（is_half_periodic_transient: bool）
- 场保存类型（save_fields_type: "Every N Steps" / "Every N Seconds" / ...）
- N步数（n_steps: string, 如"1"）
- 步数范围（steps_from/to: string, 如"0s"/"0.025s"）
- 使用非线性迭代次数（use_nonlinear_iter_num: bool）
- 缓存保存类型（cache_save_kind: "Count" / "Time"）
- 求解步数（number_solve_steps: int）
- 求解范围（range_start/end: string）
- 使用自适应时间步（use_adaptive_time_step: bool）
- 初始时间步（initial_time_step: string, 如"0.002s"）
- 最小时间步（min_time_step: string, 如"0.001s"）
- 最大时间步（max_time_step: string, 如"0.003s"）
- 时间步误差容差（time_step_err_tolerance: double, 如0.0001）

#### Scenario: 完整提取瞬态求解设置
- **WHEN** 解析Temp.aedt的Setup1
- **THEN** SolutionSetup对象应包含全部40+个瞬态参数：
  - SetupType='Transient', StopTime='0.025s', TimeStep='0.0001s'
  - NonlinearSolverResidual='0.0001'
  - FrequencyOfAddedVoltageSource='308.06667Hz'
  - **且** 自适应时间步配置：InitialTimeStep='0.002s', MinTimeStep='0.001s', MaxTimeStep='0.003s'
  - **且** 场保存：SaveFieldsType='Every N Steps', NSteps='1', StepsFrom='0s', StepsTo='0.025s'

### Requirement: 网格操作（MeshOperation）数据结构

系统SHALL新增MeshOperation类和增强Mesh类：

#### 全局网格设置（Mesh类增强）
- 曲面近似模式（surf_approx_choice: "UseSlider" / "ManualSettings"）
- 曲面近似滑块值（slider_mesh_settings: int, 1-10）
- 是否使用自动长度（use_auto_length: bool）

#### 网格操作（MeshOperation类，新增）
- 名称（name: "M" / "C" / "R" / "S" / "SurfApprox1" / "CylindricalGap1"）
- ID（id: int）
- 类型（type: "LengthBased" / "SurfApproxBased" / "CylindricalGap"）
- 启用标志（enabled: bool）
- 是否组件（is_component: bool）
- 是否全局（is_global: bool）
- 关联对象列表（objects_: 对象ID列表）
- 是否内部细化（refine_inside: bool）
- 是否限制元素数（restrict_elem: bool）
- 最大元素数（num_max_elem: string, 如"1000"）
- 是否限制长度（restrict_length: bool）
- 最大长度（max_length: string, 如"1mm"）
- 是否应用于初始网格（apply_to_initial_mesh: bool）
- 曲面近似详细设置（手动模式下）：
  - 曲面偏差选择（surf_dev_choice: int）
  - 曲面偏差值（surf_dev: string, 如"0.15mm"）
  - 法向偏差选择（normal_dev_choice: int）
  - 法向偏差值（normal_dev: string, 如"30deg"）
  - 长宽比选择（aspect_ratio_choice: int）
- Band间隙设置（CylindricalGap类型）：
  - 是否使用Band映射角度（use_band_mapping_angle: bool）
  - Band映射角度（band_mapping_angle: string, 如"1deg"）

#### Scenario: 提取完整网格设置
- **WHEN** 解析Temp.aedt的网格设置
- **THEN** Mesh对象应包含全局设置：
  - SurfApproxChoice='UseSlider', SliderMeshSettings=5
  - UseAutoLength=true
  - **且** 应创建6个MeshOperation对象：
    - SurfApprox1: 手动曲面近似, 63个对象, SurfDev='0.15mm', NormalDev='30deg'
    - M: 长度based, MaxLength='1mm', 17个永磁体/转子对象
    - C: 长度based, MaxLength='1mm', 42个定子线圈对象
    - R: 长度based, MaxLength='4mm', 1个转子铁心对象
    - S: 长度based, MaxLength='5mm', 1个定子铁心对象
    - CylindricalGap1: 柱面间隙, BandMappingAngle='1deg', 1个Band面对象

### Requirement: 设计变量（DesignVariable）数据结构

系统SHALL新增DesignVariable类来存储项目设计变量：

#### 设计变量字段
- 名称（name: "fractions" / "delta" / "conds" / "R1" / "Le1" / "Irms" / "few" / "ele_angle"）
- 值（value: string, 如 "1", "109.987399716271deg", "8", "0.00524327280994245ohm"）
- 单位（unit: string, 如 "", "deg", "", "ohm", "H", "A", "Hz", "deg"）
- 是否索引变量（is_indexed: bool, 来自VariableOrders）
- 是否非索引变量（is_non_indexed: bool, 来自NonIndexedVariables）
- 对于非索引变量：表达式（expression: string, 如 "if(Time<0.005s,20, if(Time<0.015s,200,20))"）

#### Scenario: 提取设计变量
- **WHEN** 解析Temp.aedt的设计变量
- **THEN** 系统应创建9个DesignVariable对象：
  - 8个索引变量: fractions=1, delta=109.99°, conds=8, R1=0.00524Ω, Le1=3.71μH, Irms=500A, few=200Hz, ele_angle=45°
  - 1个非索引变量: Temp1=条件表达式（时间分段温度函数）

### Requirement: 输出变量（OutputVariable）数据结构

系统SHALL新增OutputVariable类来存储场计算器输出表达式：

#### 输出变量字段
- 名称（name: "pos" / "cos0"~"cos2" / "sin0"~"sin2" / "Lad"~"Lcq" / "L_d"/"L_q" / "Flux_d"/"Flux_q" / "I_d"/"I_q" / "Irms" / "Pcu"）
- ID（id: int, 0~20）
- 表达式（expression: string, 如 "(Moving1.Position -7.5 * PI/180) * 4 + PI"）
- 结果单位（result_unit: string, 如 "deg", "", "nH", "Wb", "mA", "ohm"）
- 显示单位（display_unit: string）

#### Scenario: 提取输出变量
- **WHEN** 解析Temp.aedt的输出变量
- **THEN** 系统应创建21个OutputVariable对象，包括：
  - 位置变换: pos, cos0~2, sin0~2
  - 电感矩阵变换: Lad, Laq, Lbd, Lbq, Lcd, Lcq, L_d, L_q
  - 磁链变换: Flux_d, Flux_q
  - 电流变换: I_d, I_q
  - 衍生量: Irms（有效值）, Pcu（铜耗）

### Requirement: 温度设置（TemperatureSettings）数据结构

系统SHALL新增TemperatureSettings类：

#### 温度设置字段
- 是否启用温度依赖（include_temperature_dependence: bool）
- 是否启用反馈（enable_feedback: bool）
- 物体-温度映射（object_temperature_map: `std::unordered_map<int/*对象ID*/, std::string/*温度引用*/>`）

#### Scenario: 提取温度设置
- **WHEN** 解析Temp.aedt的温度设置
- **THEN** TemperatureSettings对象应包含：
  - IncludeTemperatureDependence=true, EnableFeedback=false
  - **且** 映射表包含50+个条目：
    - 大部分对象 → '22cel'（22°C环境温度）
    - 转子对象（2680~2924）→ 'Temp1'（时变温度）

## MODIFIED Requirements

### Requirement: MaxwellParserImpl解析能力增强

现有的`MaxwellParserImpl`需要新增以下提取方法：

1. **`extractWindings()`** - 从BoundarySetup > Boundaries中提取BoundType='Winding Group'的边界
2. **`extractMotionSetups()`** - 从MotionSetupList中提取运动设置
3. **`extractMeshOperations()`** - 从MeshSetup > MeshOperations中提取网格操作
4. **`extractDesignVariables()`** - 从ModelSetup > Properties中提取设计变量
5. **`extractOutputVariables()`** - 从OutputVariable > OutputVariables中提取输出变量
6. **`extractTemperatureSettings()`** - 从TemperatureSettings中提取温度设置
7. **`extractGlobalBoundData()`** - 从GlobalBoundData中提取全局边界数据
8. **增强`extractSingleMaterial()`** - 补充所有新材料字段的提取
9. **增强`extractSolutionSetup()`** - 补充所有瞬态求解参数

### Requirement: AedtProjectLoader转换逻辑增强

现有的`AedtProjectLoader`转换方法需要更新：

1. **`convertAndAddMaterials()`** - 映射新的Material字段
2. **新增`convertAndAddWindings()`** - 将JSON绕组数据转为Winding对象
3. **新增`convertAndAddMotionSetups()`** - 将JSON运动数据转为MotionSetup对象
4. **新增`convertAndAddMeshOperations()`** - 将JSON网格操作数据转为MeshOperation对象
5. **新增`convertAndAddDesignVariables()`** - 将JSON设计变量转为DesignVariable对象
6. **新增`convertAndAddOutputVariables()`** - 将JSON输出变量转为OutputVariable对象
7. **`convertAndAddBoundaries()`** - 映射新的Boundary字段（parent_bnd_id等）
8. **`convertAndAddSolutionSetup()`** - 映射新的SolutionSetup字段

### Requirement: ProjectManager容器扩展

ProjectManager需要新增以下容器成员和访问接口：

```cpp
// 新增容器
std::unordered_map<std::string, WindingPtr> windings_;
std::vector<MotionSetupPtr> motion_setups_;           // 通常只有1个运动设置
std::unordered_map<std::string, MeshOperationPtr> mesh_operations_;
std::unordered_map<std::string, DesignVariablePtr> design_variables_;
std::unordered_map<std::string, OutputVariablePtr> output_variables_;
std::optional<TemperatureSettingsPtr> temperature_settings_;

// 新增访问接口
void addWinding(WindingPtr winding);
void addMotionSetup(MotionSetupPtr setup);
void addMeshOperation(MeshOperationPtr operation);
void addDesignVariable(DesignVariablePtr variable);
void addOutputVariable(OutputVariablePtr variable);
void setTemperatureSettings(TemperatureSettingsPtr settings);
```

## REMOVED Requirements

无删除需求。所有现有功能保持兼容，新功能为增量添加。
