# Tasks

- [x] Task 1: 扩展 SimulationType 枚举 + 更新 SolverScheduler
  - [x] SubTask 1.1: 在 `include/tool/em_enums.hpp` 的 `SimulationType` 枚举中添加 `DC_CONDUCTION` 值（位于 ELECTROSTATIC 之后）
  - [x] SubTask 1.2: 在 `include/solver/solver_scheduler.cpp` 的 `createPhysicsField()` 方法中添加 `DC_CONDITION` 分支，创建 DCConductionSolver 实例

- [x] Task 2: 创建直流电流场单元积分器 DCConductionIntegrator
  - [x] SubTask 2.1: 创建 `include/solver/dc_conduction_integrator.hpp`，定义 DCConductionIntegrator 类（继承 EMElementIntegratorBase，刚度矩阵使用 σ 作为材料系数）
  - [x] SubTask 2.2: 创建 `src/solver/dc_conduction_integrator.cpp`，实现 computeStiffnessMatrix()（K_e = ∫∇N^T·σ·∇NdΩ）、computeSourceVector()、computeAllMatrices()
  - [x] SubTask 2.3: 支持与 ElectrostaticIntegrator 相同的单元类型列表（TRI3/QUAD4/TET4/HEX8/PRISM6/PYRAMID5 等）

- [x] Task 3: 创建 FieldData 电流场后处理扩展
  - [x] SubTask 3.1: 在 `include/solver/field_data.hpp` 中新增电流场专用接口：computeCurrentDensity()、getElementCurrentDensity()、computeJouleHeating()、computeTotalCurrent()、computeVoltageDrop()
  - [x] SubTask 3.2: 在 `src/solver/field_data.cpp` 中实现：J = σE 计算、焦耳热功率密度 p = σ|E|² 计算、总电流 I = ∫J·ndS 积分、总功率 P = φ^T K φ

- [x] Task 4: 创建电阻计算器 ResistanceCalculator
  - [x] SubTask 4.1: 创建 `include/solver/resistance_calculator.hpp`，定义 ResistanceCalculator 类（computeResistance 端子间电阻、computeTotalCurrent 端面总电流）
  - [x] SubTask 4.2: 创建 `src/solver/resistance_calculator.cpp`，实现欧姆定律法电阻计算（施加单位电压→求解→积分端面电流→R=V/I）

- [x] Task 5: 创建焦耳热计算器 JouleHeatingCalculator
  - [x] SubTask 5.1: 创建 `include/solver/joule_heating_calculator.hpp`，定义 JouleHeatingCalculator 类（computeTotalPower、computeElementPowerDensity）
  - [x] SubTask 5.2: 创建 `src/solver/joule_heating_calculator.cpp`，实现 P_total = φ^T K φ 和 p_element = σ|E|²

- [x] Task 6: 创建直流电流场求解器 DCConductionSolver
  - [x] SubTask 6.1: 创建 `include/solver/dc_conduction_solver.hpp`，定义 DCConductionSolver 类（继承 PhysicsField，复用 ElectrostaticSolver 的完整架构模式）
  - [x] SubTask 6.2: 创建 `src/solver/dc_conduction_solver.cpp`，实现 setup()（前检查σ>0 → 材料映射(σ作为主参数) → DOF分配 → 装配(DCConductionIntegrator) → BC施加 → 激励施加）、solve()（SPD求解）、postProcess()（E/J/P 计算）
  - [x] SubTask 6.3: 实现 preCheck() 中电导率合法性校验（σ ≥ 0，导体区域必须有 σ > 0）

- [x] Task 7: 集成到 SolverScheduler + 测试用例
  - [x] SubTask 7.1: 验证 SolverScheduler 对 DC_CONDUCTION 类型的自动创建和调度
  - [x] SubTask 7.2: 创建 `test/test_dc_conduction_solver.cpp`
  - [x] SubTask 7.3: 实现均匀电阻条测试（电位线性分布、电流均匀、电阻 R=L/(σA) 精度<5%）
  - [x] SubTask 7.4: 实现焦耳热验证测试（P=V²/R=I²R）
  - [x] SubTask 7.5: 实现边界条件施加测试（Dirichlet电压/Neumann电流）
  - [x] SubTask 7.6: 实现错误处理测试（零电导率检测）
  - [x] SubTask 7.7: 更新 `test/CMakeLists.txt` 添加 test_dc_conduction_solver 目标

# Task Dependencies

- [Task 1] depends on [] (无依赖，可最先开始)
- [Task 2] depends on [Task 1] (需要 SimulationType 枚举已定义)
- [Task 3] depends on [] (FieldData 扩展独立于其他任务)
- [Task 4] depends on [Task 6] (ResistanceCalculator 需要 DCConductionSolver 进行求解)
- [Task 5] depends on [Task 6] (JouleHeatCalculator 需要 DCConductionSolver 的结果)
- [Task 6] depends on [Task 2, Task 3] (DCConductionSolver 需要积分器和扩展后的 FieldData)
- [Task 7] depends on [Task 4, Task 5, Task 6] (集成测试需要全部核心模块)

# Parallelizable Work

- Task 1 (枚举扩展) 和 Task 3 (FieldData扩展) 可并行开发
- Task 2 (DCConductionIntegrator) 可与 Task 1/3 并行开发
- Task 4 (ResistanceCalculator) 和 Task 5 (JouleHeatCalculator) 可并行开发
