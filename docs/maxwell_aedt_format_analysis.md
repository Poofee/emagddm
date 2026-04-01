# Maxwell AEDT 文件数据格式分析

## 概述

本文档基于对 `Project49.aedt` 文件的分析，详细描述了Ansys Maxwell软件使用的AEDT（Ansys Electronics Desktop）文件格式。AEDT文件采用层次化的文本格式存储电磁仿真项目的完整配置信息。

## 文件结构

### 1. 基本结构
AEDT文件采用类似XML的层次结构，但使用`$begin`和`$end`标签来定义块级元素：

```
$begin 'AnsoftProject'
    Created='Tue Feb 10 09:02:30 2026'
    Product='ElectronicsDesktop'
    FileOwnedByWorkbench=false
    $begin 'Desktop'
        Version(2022, 1)
        InfrastructureVersion(1, 0)
    $end 'Desktop'
    ...
$end 'AnsoftProject'
```

### 2. 项目元数据
- **创建时间**：`Created='Tue Feb 10 09:02:30 2026'`
- **产品类型**：`Product='ElectronicsDesktop'`
- **版本信息**：`Version(2022, 1)`
- **环境配置**：包含HFSS、Maxwell、Q3D等多种仿真环境

## 材料属性数据格式

### 材料定义结构
```
$begin 'vacuum'
    CoordinateSystemType='Cartesian'
    BulkOrSurfaceType=1
    $begin 'PhysicsTypes'
        set('Electromagnetic')
    $end 'PhysicsTypes'
    $begin 'AttachedData'
        $begin 'MatAppearanceData'
            property_data='appearance_data'
            Red=230
            Green=230
            Blue=230
            Transparency=0.949999988079071
        $end 'MatAppearanceData'
    $end 'AttachedData'
    permittivity='1'
    ModTime=1499970477
    Library='Materials'
    LibLocation='SysLibrary'
    ModSinceLib=false
$end 'vacuum'
```

### 材料属性类型
1. **电磁属性**：
   - `permittivity`：介电常数
   - `permeability`：磁导率
   - `conductivity`：电导率
   - `magnetic_coercivity`：矫顽力

2. **热学属性**：
   - `thermal_conductivity`：热导率
   - `specific_heat`：比热容
   - `thermal_expansion_coefficient`：热膨胀系数

3. **结构属性**：
   - `youngs_modulus`：杨氏模量
   - `poissons_ratio`：泊松比
   - `mass_density`：质量密度

4. **外观属性**：
   - RGB颜色值
   - 透明度设置

## 几何模型数据格式

### 几何体定义
```
$begin 'GeometryPart'
    $begin 'Attributes'
        Name='Box1'
        Color='(143 175 143)'
        MaterialValue='"TDK_FB6H_60cel"'
        SolveInside=true
    $end 'Attributes'
    $begin 'Operations'
        $begin 'Operation'
            OperationType='Box'
            $begin 'BoxParameters'
                XPosition='-0.6mm'
                YPosition='-2.2mm'
                ZPosition='0mm'
                XSize='1.4mm'
                YSize='1mm'
                ZSize='1.8mm'
            $end 'BoxParameters'
        $end 'Operation'
    $end 'Operations'
$end 'GeometryPart'
```

### 几何操作类型
- **Box**：长方体
- **拓扑信息**：包含面、边、顶点数量统计
- **坐标系**：基于笛卡尔坐标系
- **单位系统**：`Units='mm'`

## 边界条件和激励源数据格式

### 边界条件类型
1. **绝缘边界**：`BoundType='Insulating'`
2. **磁场边界**：`BoundType='Tangential H Field'`
3. **对称边界**：`BoundType='Symmetry'`
4. **独立/依赖边界对**：用于周期性边界条件

### 激励源类型
1. **电流激励**：
   ```
   $begin 'Current1'
       BoundType='Current'
       Current='10mA'
       IsSolid=true
   $end 'Current1'
   ```

2. **电流密度激励**：
   ```
   $begin 'CurrentDensity1'
       BoundType='Current Density'
       CurrentDensityX='10'
       CurrentDensityY='0'
       CurrentDensityZ='0'
   $end 'CurrentDensity1'
   ```

3. **电压激励**：
   ```
   $begin 'Voltage1'
       BoundType='Voltage'
       Voltage='10mV'
   $end 'Voltage1'
   ```

## 求解器设置数据格式

### 求解器类型
- **静磁求解**：`SolutionType='Magnetostatic'`
- **涡流求解**：`SolutionType='EddyCurrent'`

### 求解参数
```
$begin 'Setup1'
    ID=0
    SetupType='Magnetostatic'
    MaximumPasses=10
    MinimumPasses=2
    MinimumConvergedPasses=1
    PercentRefinement=30
    PercentError=1
    UseIterativeSolver=false
    RelativeResidual=1e-06
    NonLinearResidual=0.001
$end 'Setup1'
```

### 网格设置
```
$begin 'MeshSettings'
    $begin 'GlobalSurfApproximation'
        CurvedSurfaceApproxChoice='UseSlider'
        SliderMeshSettings=5
    $end 'GlobalSurfApproximation'
    MeshMethod='Auto'
$end 'MeshSettings'
```

## 结果和数据集数据格式

### 求解结果管理
```
$begin 'SolutionManager'
    $begin 'SimSetup'
        TypeName='BaseSetup'
        ID=18
        Name='Setup1'
        $begin 'Solution'
            ID=19
            Name='AdaptivePass'
        $end 'Solution'
        $begin 'Solution'
            ID=20
            Name='LastAdaptive'
        $end 'Solution'
    $end 'SimSetup'
$end 'SolutionManager'
```

### 后处理扫描
```
$begin 'PostprocessSweep'
    Variable='NormalizedDistance'
    RegularSweep=1
    Units=''
    Minimum=0
    Maximum=1
    Increment=0.01
$end 'PostprocessSweep'
```

## 数据格式特点总结

### 1. 层次化结构
- 使用`$begin`/`$end`标签定义嵌套层次
- 每个块包含特定的属性和子块
- 支持复杂的数据关系描述

### 2. 类型化属性
- 数值属性：`Current='10mA'`
- 字符串属性：`Name='Box1'`
- 布尔属性：`IsSolid=true`
- 数组属性：`Light0[4: 6710886, 0, -1, -0.150000005960464]`

### 3. 单位系统
- 明确指定物理量单位
- 支持国际单位制
- 自动单位转换

### 4. 版本控制
- 包含版本标识符
- 支持向后兼容
- 模块化版本管理

### 5. 库管理
- 系统库和项目库分离
- 材料库引用机制
- 符号和组件库管理

## 与FETI-DP求解器的相关性

### 数据映射关系
1. **材料属性** → 子域材料参数矩阵
2. **几何分区** → 子域划分依据
3. **边界条件** → 界面连续性条件
4. **求解设置** → 迭代收敛准则

### 转换建议
1. 提取材料电磁参数用于构建子域刚度矩阵
2. 根据几何分区自动生成子域划分
3. 将边界条件映射为FETI-DP的拉格朗日乘子
4. 继承求解精度和收敛性设置

## 结论

Maxwell AEDT文件格式提供了完整的电磁仿真项目配置信息，采用层次化的文本格式存储，便于解析和转换。该格式包含了材料属性、几何模型、边界条件、求解设置等关键信息，为开发FETI-DP求解器提供了丰富的数据源。

建议在开发过程中建立AEDT到FETI-DP数据格式的转换工具，充分利用现有的Maxwell项目配置，提高开发效率和兼容性。