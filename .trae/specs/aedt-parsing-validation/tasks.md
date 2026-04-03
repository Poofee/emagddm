# Tasks

- [x] Task 1: 增强MaxwellParser解析器健壮性
  - [x] SubTask 1.1: 扩展PROPERTY_PATTERN正则表达式，支持带空格/特殊字符的属性名（如`'Conductor number'`、`'N Steps'`）
  - [x] SubTask 1.2: 增强函数调用形式属性的解析（如`Version(2024, 1)`、`Edges(18)`、`Objects(2331)`）
  - [x] SubTask 1.3: 处理Name()形式特殊属性（如`Name('Bt1')`）
  - [x] SubTask 1.4: 处理长数组/键值对形式属性（如Temperatures多参数列表）

- [x] Task 2: 实现MaxwellParserImpl材料数据完整提取
  - [x] SubTask 2.1: 完善`extractMaterials()`方法，从Definitions/Materials块提取所有材料
  - [x] SubTask 2.2: 正确解析线性磁导率（permeability='1.0739096'）
  - [x] SubTask 2.3: 正确解析非线性B-H曲线（permeability块内的BHCoordinates/Points数组）
  - [x] SubTask 2.4: 解析磁芯损耗系数（core_loss_kh, core_loss_kc, core_loss_ke）
  - [x] SubTask 2.5: 解析矫顽力矢量（magnetic_coercivity块的Magnitude和DirComp）
  - [x] SubTask 2.6: 解析电导率、质量密度等基本电磁属性

- [x] Task 3: 实现MaxwellParserImpl边界条件数据提取
  - [x] SubTask 3.1: 实现`extractBoundaries()`方法，从Maxwell2DModel/BoundarySetup/Boundaries块提取
  - [x] SubTask 3.2: 解析边界类型（BoundType: Vector Potential / Coil等）
  - [x] SubTask 3.3: 解析关联边/对象列表（Edges(n)/Objects(n)）
  - [x] SubTask 3.4: 解析线圈参数（Winding, PolarityType, Conductor number, ParentBndID）

- [x] Task 4: 实现MaxwellParserImpl激励源数据提取
  - [x] SubTask 4.1: 实现`extractExcitations()`方法，从边界条件中识别并提取线圈激励源
  - [x] SubTask 4.2: 关联线圈边界与绕组信息

- [x] Task 5: 实现MaxwellParserImpl求解设置数据提取
  - [x] SubTask 5.1: 实现`extractSolutionSetup()`方法，从AnalysisSetup/SolveSetups/Setup1提取
  - [x] SubTask 5.2: 解析瞬态求解参数（StopTime, TimeStep, NonlinearSolverResidual等）
  - [x] SubTask 5.3: 解析自适应时间步设置

- [x] Task 6: 实现MaxwellParserImpl几何与网格设置数据提取
  - [x] SubTask 6.1: 实现`extractGeometry()`方法，提取基本几何信息（设计名称、模型深度、几何模式）
  - [x] SubTask 6.2: 实现网格设置摘要提取

- [x] Task 7: 编写AEDT解析验证测试（test_aedt_parsing.cpp）
  - [x] SubTask 7.1: 创建测试文件，使用Google Test框架
  - [x] SubTask 7.2: 编写项目元数据验证测试用例
  - [x] SubTask 7.3: 编写材料数据正确性验证测试用例（copper, vacuum, B27AV1400, N45UH永磁体）
  - [x] SubTask 7.4: 编写B-H曲线数据验证测试用例
  - [x] SubTask 7.5: 编写边界条件验证测试用例（VectorPotential1, PhA_0线圈）
  - [x] SubTask 7.6: 编写求解设置验证测试用例（Setup1瞬态配置）
  - [x] SubTask 7.7: 编材料/边界/求解设置数量断言测试

- [x] Task 8: 更新CMakeLists.txt并编译验证
  - [x] SubTask 8.1: 在test/CMakeLists.txt中添加test_aedt_parsing目标
  - [x] SubTask 8.2: 编译并通过所有测试

# Task Dependencies
- [Task 2, 3, 4, 5, 6] depends on [Task 1]
- [Task 7] depends on [Task 1, 2, 3, 4, 5, 6]
- [Task 8] depends on [Task 7]
