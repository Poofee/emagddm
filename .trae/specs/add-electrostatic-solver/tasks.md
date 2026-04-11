# Tasks

- [x] Task 1: 创建物理场抽象基类 PhysicsField
  - [x] SubTask 1.1: 创建 include/solver/physics_field.hpp，定义PhysicsField抽象基类接口（setup/solve/postProcess/getFieldData/getSimulationType/getSolverName/clear）
  - [x] SubTask 1.2: 创建 src/solver/ 目录结构
  - [x] SubTask 1.3: 更新 src/CMakeLists.txt 添加 solver 子目录

- [x] Task 2: 创建场数据与后处理模块 FieldData
  - [x] SubTask 2.1: 创建 include/solver/field_data.hpp，定义FieldData类（节点电位存储、电场强度/电位移矢量计算、静电能量计算、VTK/CSV导出接口）
  - [x] SubTask 2.2: 创建 src/solver/field_data.cpp，实现核心后处理计算逻辑（电场强度E=-∇φ、电位移D=εE、能量W=0.5φᵀKφ）
  - [x] SubTask 2.3: 实现VTK格式导出（节点电位标量数据 + 单元电场矢量数据）
  - [x] SubTask 2.4: 实现CSV格式导出

- [x] Task 3: 创建边界条件管理器 BoundaryConditionManager
  - [x] SubTask 3.1: 创建 include/solver/boundary_condition_manager.hpp，定义BoundaryConditionManager类
  - [x] SubTask 3.2: 创建 src/solver/boundary_condition_manager.cpp，实现Dirichlet边界条件施加（消去法：DOF约束标记 + 右端项修正）
  - [x] SubTask 3.3: 实现Neumann边界条件处理（自然边界条件 + 非零法向导数面积分）
  - [x] SubTask 3.4: 实现Robin/阻抗边界条件（矩阵对角线修改 + 右端项修改）
  - [x] SubTask 3.5: 实现边界条件冲突检测与优先级管理
  - [x] SubTask 3.6: 实现PERIODIC/ANTIPERIODIC周期性边界条件（主从DOF约束）

- [x] Task 4: 创建激励管理器 ExcitationManager
  - [x] SubTask 4.1: 创建 include/solver/excitation_manager.hpp，定义ExcitationManager类
  - [x] SubTask 4.2: 创建 src/solver/excitation_manager.cpp，实现电压源激励→Dirichlet边界条件转化
  - [x] SubTask 4.3: 实现体电荷密度源→右端项体积分贡献
  - [x] SubTask 4.4: 实现激励合法性校验

- [x] Task 5: 创建静电场求解器 ElectrostaticSolver
  - [x] SubTask 5.1: 创建 include/solver/electrostatic_solver.hpp，定义ElectrostaticSolver类（继承PhysicsField）
  - [x] SubTask 5.2: 创建 src/solver/electrostatic_solver.cpp，实现setup()方法（前处理检查→材料映射→DOF分配→矩阵装配→边界条件施加→激励施加）
  - [x] SubTask 5.3: 实现solve()方法（选择线性求解器→求解Kφ=F→解向量回扩）
  - [x] SubTask 5.4: 实现postProcess()方法（调用FieldData计算电场/电位移/能量）
  - [x] SubTask 5.5: 实现轴对称模式支持（积分权重乘以2πr因子，r=0奇异性处理）
  - [x] SubTask 5.6: 实现前处理检查（网格有效性、材料合理性、边界条件完整性、冲突检测）

- [x] Task 6: 创建电容矩阵计算器 CapacitanceCalculator
  - [x] SubTask 6.1: 创建 include/solver/capacitance_calculator.hpp，定义CapacitanceCalculator类
  - [x] SubTask 6.2: 创建 src/solver/capacitance_calculator.cpp，实现能量法电容矩阵计算（多次求解+能量提取）
  - [x] SubTask 6.3: 实现电荷法电容计算（高斯面积分计算导体电荷）

- [x] Task 7: 创建静电力计算器 ElectrostaticForceCalculator
  - [x] SubTask 7.1: 创建 include/solver/electrostatic_force_calculator.hpp，定义ElectrostaticForceCalculator类
  - [x] SubTask 7.2: 创建 src/solver/electrostatic_force_calculator.cpp，实现虚功法静电力计算
  - [x] SubTask 7.3: 实现虚功法扭矩计算

- [x] Task 8: 创建求解器调度器 SolverScheduler
  - [x] SubTask 8.1: 创建 include/solver/solver_scheduler.hpp，定义SolverScheduler类
  - [x] SubTask 8.2: 创建 src/solver/solver_scheduler.cpp，实现initialize()（从ProjectManager提取数据）
  - [x] SubTask 8.3: 实现run()（根据SimulationType创建PhysicsField→setup→solve→postProcess）
  - [x] SubTask 8.4: 实现结果输出调度（调用FieldData导出接口）

- [x] Task 9: 集成到SolverApp
  - [x] SubTask 9.1: 修改 src/app/solver_app.cpp 的 run() 方法，集成SolverScheduler调用流程
  - [x] SubTask 9.2: 更新JSON配置文件格式，添加静电场求解器配置项（SimulationType、导体ID列表等）

- [x] Task 10: 编写测试用例
  - [x] SubTask 10.1: 创建 test/test_electrostatic_solver.cpp
  - [x] SubTask 10.2: 实现2D平行板电容器测试（电位线性分布、电场均匀、电容精度<5%）
  - [x] SubTask 10.3: 实现2D同轴电缆测试（圆柱坐标解析解对比，误差<2%）
  - [x] SubTask 10.4: 实现3D立方体电容器测试（3D求解流程验证）
  - [x] SubTask 10.5: 实现轴对称同轴电缆测试（轴对称积分权重验证）
  - [x] SubTask 10.6: 实现多导体电容矩阵测试（3导体系统，对称性验证）
  - [x] SubTask 10.7: 实现边界条件施加测试（Dirichlet/Neumann/冲突检测）
  - [x] SubTask 10.8: 实现错误处理测试（缺少边界条件、材料异常）
  - [x] SubTask 10.9: 更新 test/CMakeLists.txt 添加测试目标

# Task Dependencies

- [Task 2] depends on [Task 1] (FieldData被PhysicsField接口引用)
- [Task 3] depends on [Task 1] (BoundaryConditionManager被ElectrostaticSolver使用)
- [Task 4] depends on [Task 1] (ExcitationManager被ElectrostaticSolver使用)
- [Task 5] depends on [Task 1, Task 2, Task 3, Task 4] (ElectrostaticSolver组合所有模块)
- [Task 6] depends on [Task 5] (CapacitanceCalculator需要调用ElectrostaticSolver)
- [Task 7] depends on [Task 5] (ElectrostaticForceCalculator需要调用ElectrostaticSolver)
- [Task 8] depends on [Task 5] (SolverScheduler调度ElectrostaticSolver)
- [Task 9] depends on [Task 8] (SolverApp集成SolverScheduler)
- [Task 10] depends on [Task 5, Task 6, Task 7] (测试需要完整求解器功能)

# Parallelizable Work

- Task 2 (FieldData) 和 Task 3 (BoundaryConditionManager) 和 Task 4 (ExcitationManager) 可并行开发
- Task 6 (CapacitanceCalculator) 和 Task 7 (ElectrostaticForceCalculator) 可并行开发
