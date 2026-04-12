# 直流电流场求解器（DC Conduction Solver）Spec

## Why
项目已完成静电场求解器（ElectrostaticSolver）的完整开发，包括PhysicsField抽象基类、FieldData后处理、BoundaryConditionManager、ExcitationManager、SolverScheduler等全套框架。根据 `/Users/poofee/emagddm/docs/solver/dccond.md` 文档，直流电流场（DC Conduction）与静电场共享**90%以上的底层基础代码**（网格、DOF管理、矩阵装配、线性求解器、边界条件施加），核心差异仅在：
- 控制方程：∇·(σ∇V) = 0 （电流连续性方程，用σ替代ε）
- 材料参数：电导率σ（S/m）替代介电常数ε（F/m）
- 激励源：电流源激励（Neumann边界）是核心功能
- 后处理：电阻/电导/焦耳热/电流密度 替代 电容/静电能量/电位移矢量

对标ANSYS Maxwell DC Conduction Solver，需支持2D/3D/轴对称模型、电压/电流边界条件、电阻计算、焦耳热计算等工业核心功能。

## What Changes
- **新增** 直流电流场积分器 `DCConductionIntegrator`（或复用ElectrostaticIntegrator并调整材料系数映射策略）
- **新增** 直流电流场求解器类 `DCConductionSolver`，派生自 `PhysicsField`
- **新增** 电阻计算模块 `ResistanceCalculator`（任意两端子间直流电阻 R = V/I）
- **新增** 焦耳热计算模块 `JouleHeatingCalculator`（P = ∫ σ|E|² dV）
- **新增** 电流场专用后处理扩展到 `FieldData`（电流密度 J=σE、焦耳热功率密度、总电流）
- **修改** `tool::em_enums.hpp` 的 `SimulationType` 枚举添加 `DC_CONDUCTION`
- **修改** `SolverScheduler::createPhysicsField()` 支持 `DC_CONDITION` 类型创建
- **零修改** 已有底层模块（EMDOFManager、EMAssembly、EMLinearSolverBase、BoundaryConditionManager、ExcitationManager）

## Impact
- Affected specs: add-electrostatic-solver（复用其架构模式）
- Affected code:
  - `include/tool/em_enums.hpp` — 添加 DC_CONDUCTION 枚举值
  - `include/solver/dc_conduction_solver.hpp` — 直流电流场求解器头文件（新建）
  - `include/solver/dc_conduction_integrator.hpp` — DC电流场单元积分器（新建）
  - `include/solver/resistance_calculator.hpp` — 电阻计算器（新建）
  - `include/solver/joule_heating_calculator.hpp` — 焦耳热计算器（新建）
  - `src/solver/` — 所有实现文件（新建）
  - `src/solver/solver_scheduler.cpp` — 添加 DC_CONDUCTION 分支
  - `test/test_dc_conduction_solver.cpp` — 测试用例（新建）

## ADDED Requirements

### Requirement: SimulationType 枚举扩展

系统 SHALL 在 `tool::SimulationType` 枚举中添加 `DC_CONDUCTION` 值。

### Requirement: DCConductionIntegrator 直流电流场单元积分器

系统 SHALL 提供 `DCConductionIntegrator` 类，实现直流电流场的单元级有限元矩阵计算。

#### 控制方程
- **二维/三维直流电流场**：∇·(σ∇V) = 0 （拉普拉斯方程形式）
  - 未知量：标量电位 V（伏特）
  - 材料参数：电导率 σ（S/m），可各向异性
  - 无源项（稳态恒定电流，无体电荷注入）

#### 弱形式
∫Ω ∇N^T · σ · ∇V dΩ = ∫Γ N · J_n dΓ

其中 J_n = -σ·∂V/∂n 为边界法向电流密度（Neumann条件）。

#### 单元刚度矩阵公式
K_e = ∫_Ωe (∇N)^T · σ · (∇N) dΩ

与静电场的区别仅在于材料系数：**σ 替代 ε**。

#### 支持的维度与单元

| 维度 | 单元类型 | DOF类型 |
|------|---------|---------|
| 2D | TRI3, TRI6, QUAD4, QUAD8, QUAD9 | SCALAR_ONLY |
| 3D | TET4, TET10, HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13 | SCALAR_ONLY |
| 轴对称 | TRI3, TRI6, QUAD4, QUAD8, QUAD9 | SCALAR_ONLY |

#### 实现策略
两种可选方案（推荐方案A）：

**方案A（推荐）：新建独立积分器**
- 新建 `DCConductionIntegrator` 类，继承 `EMElementIntegratorBase`
- 刚度矩阵使用 σ 作为材料系数
- 与 ElectrostaticIntegrator 结构完全平行，便于后续独立演化

**方案B：复用 ElectrostaticIntegrator**
- 在调用时将 MaterialProperties 中的 epsilon 字段设为 σ 值
- 利用现有 K_e = ∫(∇N)^T·ε·(∇N)dΩ 公式，将ε解释为σ
- 优点：零新代码；缺点：语义不清晰，调试困难

#### Scenario: 导体材料刚度矩阵计算
- **GIVEN** 一个铜导体单元（σ = 5.8×10⁷ S/m）
- **WHEN** 调用 DCConductionIntegrator::computeStiffnessMatrix()
- **THEN** 返回的 K_e 中所有元素正比于 5.8×10⁷
- **AND** 矩阵对称正定

---

### Requirement: DCConductionSolver 直流电流场求解器

系统 SHALL 提供 `DCConductionSolver` 类，继承 `PhysicsField`，实现完整的直流电流场有限元求解流程，对标Maxwell DC Conduction Solver。

#### 核心求解流程

```
1. setup()阶段：
   a. 前处理检查（网格有效性、材料参数合理性——重点校验电导率σ>0）
   b. 构建MaterialProperties映射表（从tool::Material提取σ作为主参数）
   c. DOF分配（EMDOFManager，SCALAR_ONLY类型）
   d. 全局矩阵装配（EMAssembly，调用DCConductionIntegrator）
      - 关键：MaterialProperties.epsilon字段存储的是σ值（电导率）
   e. 边界条件施加（BoundaryConditionManager）
      - Dirichlet: 电压边界 V = V₀
      - Neumann: 电流边界 J_n = J₀ → 右端项面积分 ∫ N·J_n dΓ
   f. 激励施加（ExcitationManager）
      - VOLTAGE_SOURCE → Dirichlet边界
      - CURRENT_SOURCE → Neumann边界（右端项贡献）

2. solve()阶段：
   a. 选择线性求解器（SPD矩阵→Cholesky/CGSolver）
   b. 求解线性方程组 KV = F
   c. 解向量回扩

3. postProcess()阶段：
   a. 计算电场强度 E = -∇V
   b. 计算电流密度 J = σE
   c. 计算焦耳热功率 P = ∫ σ|E|² dV
   d. 存储到FieldData
```

#### Scenario: 简单电阻条求解
- **GIVEN** 一个均匀截面电阻条模型（长度L，截面积A，材料电导率σ）
- **AND** 一端施加Dirichlet边界 V=10V，另一端接地 V=0V
- **WHEN** 调用 setup() → solve() → postProcess()
- **THEN** 电位沿长度方向线性分布
- **AND** 电场强度 E ≈ 10/L（均匀）
- **AND** 电流密度 J = σE ≈ 10σ/L
- **AND** 总电流 I = J·A = 10σA/L
- **AND** 计算电阻 R = V/I = L/(σA)，与解析解误差 < 5%

#### Scenario: PCB走线电流分布分析
- **GIVEN** 一个PCB铜走线模型（TET4单元，铜σ=5.8×10⁷ S/m）
- **AND** 一端施加电压源 V=1V，另一端接地
- **WHEN** 执行完整求解流程
- **THEN** 求解成功，电流密度在拐角处集中（趋肤效应不考虑但几何效应存在）
- **AND** 总电流值合理（可手工估算验证）

---

### Requirement: ResistanceCalculator 电阻计算器

系统 SHALL 提供 `ResistanceCalculator` 类，计算任意两个端子之间的直流电阻。

#### 计算方法

采用**欧姆定律法**：
1. 对端子i施加单位电压 V_i = 1V，端子j接地 V_j = 0V，其余绝缘
2. 求解直流电流场问题
3. 计算端子i表面的总流出电流 I_i = ∫_Γi J·dS = ∫_Γi σE·ndS
4. 电阻 R_ij = V_i / I_i = 1 / I_i

```cpp
class ResistanceCalculator {
public:
    double computeResistance(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        int terminal1_id,       // 高压端端子ID
        int terminal2_id,       // 接地端端子ID
        tool::DimType dim_type = tool::DimType::D2
    );

    double computeTotalCurrent(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<Eigen::Vector3d>& elem_current_density,
        int terminal_region_id
    );
};
```

#### Scenario: 均匀电阻条电阻计算
- **GIVEN** 铜电阻条（L=0.1m, A=1mm², σ=5.8×10⁷ S/m）
- **WHEN** 调用 computeResistance()
- **THEN** R = L/(σA) ≈ 1.724 Ω，与解析解误差 < 5%

#### Scenario: PCB过孔电阻计算
- **GIVEN** PCB过孔模型（镀铜层，复杂几何）
- **WHEN** 调用 computeResistance()
- **THEN** 返回合理的过孔电阻值（毫欧级别）

---

### Requirement: JouleHeatingCalculator 焦耳热计算器

系统 SHALL 提供 `JouleHeatingCalculator` 类，计算导电区域内的焦耳热损耗。

#### 计算方法

焦耳热功率密度：p = J·E = σ|E|² （W/m³）
总焦耳热：P_total = ∫_Ω p dV = ∫_Ω σ|E|² dV = φ^T · K · φ

```cpp
class JouleHeatingCalculator {
public:
    double computeTotalPower(
        const numeric::CsrMatrix<double>& K,
        const Eigen::VectorXd& phi
    );

    std::vector<double> computeElementPowerDensity(
        const std::vector<Eigen::Vector3d>& elem_E_field,
        const std::map<int, numeric::MaterialProperties>& materials,
        const fe_em::EMMeshData& mesh_data
    );
};
```

#### Scenario: 均匀导体焦耳热计算
- **GIVEN** 均匀导体（σ=5.8×10⁷ S/m），施加电压差 ΔV=1V，电阻R=1Ω
- **WHEN** 调用 computeTotalPower()
- **THEN** P = ΔV²/R = 1W，误差 < 5%

#### Scenario: 功率密度分布输出
- **GIVEN** 求解完成的PCB走线模型
- **WHEN** 调用 computeElementPowerDensity()
- **THEN** 返回每个单元的焦耳热功率密度
- **AND** 可通过FieldData导出到VTK进行热分布可视化

---

### Requirement: FieldData 电流场后处理扩展

系统 SHALL 扩展 `FieldData` 类，增加电流场专用的后处理计算接口。

#### 新增接口

```cpp
// ========== 电流密度计算 ==========
void computeCurrentDensity(
    const fe_em::EMMeshData& mesh_data,
    const std::map<int, numeric::MaterialProperties>& materials
);
const std::vector<Eigen::Vector3d>& getElementCurrentDensity() const;

// ========== 焦耳热计算 ==========
double computeJouleHeating(
    const numeric::CsrMatrix<double>& K,
    const Eigen::VectorXd& phi
) const;
double computeTotalCurrent(
    const fe_em::EMMeshData& mesh_data,
    const std::vector<Eigen::Vector3d>& elem_J
) const;

// ========== 电压降计算 ==========
double computeVoltageDrop(
    int node_id_1,
    int node_id_2
) const;
```

#### 物理量公式
- **电流密度**：J = σE（A/m²）
- **焦耳热功率密度**：p = |J|²/σ = σ|E|²（W/m³）
- **总焦耳热**：P = ∫ p dV = φ^T K φ（W）
- **总电流**：I = ∫_S J·n dS（A）
- **电压降**：ΔV = V(node₁) - V(node₂)（V）

#### Scenario: 电流密度分布正确性
- **GIVEN** 求解完成的均匀电阻条模型
- **WHEN** 调用 computeCurrentDensity()
- **THEN** 每个单元的电流密度 J 为常数向量（一阶单元）
- **AND** J的方向从高电位指向低电位
- **AND** |J| = σ × |E|

---

### Requirement: SolverScheduler 扩展

系统 SHALL 修改 `SolverScheduler::createPhysicsField()` 方法，支持 `DC_CONDITION` 类型。

```cpp
case tool::SimulationType::DC_CONDUCTION:
    FEEM_INFO("创建直流电流场求解器 (DCConductionSolver)");
    return std::make_unique<DCConductionSolver>(dim_type);
```

---

### Requirement: 边界条件与激励扩展

系统 SHALL 扩展 BoundaryConditionManager 和 ExcitationManager 以支持电流场特有功能。

#### 电流边界条件（Neumann类型的新用法）
- **CURRENT 类型边界**：指定流入/流出的总电流值 I（A）
  - 转化为右端项：F_i += ∫_Γ N_i · (I/A_boundary) dΓ
  - 其中 A_boundary 为边界的总面积/长度
- **RESISTANCE 类型边界**（薄层电阻）：模拟表面接触电阻
  - 通过Robin边界条件实现：σ_s · V + J_n = 0
  - 其中 σ_s 为表面电导（S/m²）

#### 激励源扩展
- **CURRENT_SOURCE 激励**：转化为Neumann边界条件（电流流入/流出）
- **GROUND 激励**：特殊电压源 V=0

---

### Requirement: 测试用例

系统 SHALL 提供完整的直流电流场求解器测试套件。

#### 测试用例清单
1. **均匀电阻条**：验证电位线性分布、电流均匀、电阻精度（R=L/(σA)）
2. **并联/串联电阻网络**：验证等效电阻计算
3. **PCB走线电流分布**：验证3D求解流程和电流密度方向
4. **接地电极电阻**：验证半无限域近似
5. **焦耳热计算**：验证 P = V²/R = I²R
6. **边界条件施加**：验证Dirichlet电压/Neumann电流边界
7. **错误处理**：验证零电导率材料检测、缺少边界条件报警
8. **与解析解对比**：简单几何（圆柱/矩形）的精确对比

## MODIFIED Requirements

### Requirement: SimulationType 枚举

在 `tool::SimulationType` 中添加 `DC_CONDUCTION` 值，位于 `ELECTROSTATIC` 之后：

```cpp
enum class SimulationType {
    ELECTROSTATIC,     ///< 静电场
    DC_CONDUCTION,      ///< 直流电流场（新增）
    MAGNETOSTATIC,     ///< 静磁场
    EDDYCURRENT,       ///< 涡流场
    TRANSIENT          ///< 瞬态场
};
```

### Requirement: SolverScheduler 工厂方法

在 `createPhysicsField()` 的 switch-case 中添加 `DC_CONDUCTION` 分支。

## REMOVED Requirements

无

## 技术约束

1. **语言标准**：C++17及以上
2. **核心依赖**：STL + Eigen 3.4+
3. **代码复用原则**：最大化复用 ElectrostaticSolver 的架构模式（setup/solve/postProcess 三阶段），仅替换物理特化部分
4. **日志规范**：使用 FEEM_INFO/FEEM_DEBUG/FEEM_ERROR 宏
5. **头文件规范**：仅包含文件名，无路径
6. **命名规范**：PascalCase类名、camelCase函数、snake_case变量
7. **注释要求**：Doxygen风格，关键逻辑必须注释
8. **内存安全**：智能指针管理动态资源，禁止内存泄漏
9. **禁止硬编码**：所有参数支持JSON配置
10. **禁止全局变量**

## 架构设计决策

### 代码复用策略

```
                    PhysicsField (抽象基类)
                         |
         ┌───────────────┼───────────────┐
         ↓               ↓               ↓
  ElectrostaticSolver  DCConductionSolver  MagnetostaticSolver (待开发)
         |               |               |
  ElectrostaticInteg  DCConductionInteg   MagneticScalar2DInteg
  (K = ∫∇N^T·ε·∇NdΩ) (K = ∫∇N^T·σ·∇NdΩ)  (K = ∫∇N^T·(1/μ)·∇NdΩ)
         |               |               
  FieldData(E,D,W)   FieldData(E,J,P)    
  CapacitanceCalc    ResistanceCalc     
  ElectrostaticForce  JouleHeatCalc      
```

**关键洞察**：DCConductionSolver 与 ElectrostaticSolver 的区别仅为：
1. 积分器中的材料系数（σ vs ε）
2. 后处理的物理量（J/P vs D/W）
3. 特有的后处理计算器（电阻/焦耳热 vs 电容/静电力）

因此 DCConductionSolver 可以**几乎原样复制** ElectrostaticSolver 的骨架代码，仅做上述三处替换。
