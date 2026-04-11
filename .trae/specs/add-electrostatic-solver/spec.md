# 静电场求解器（Electrostatic Solver）Spec

## Why
当前项目已完成有限元底层模块（网格、DOF管理、形函数、单元积分器、全局装配、线性求解器），但缺少**物理场求解器编排层**——即从项目数据出发，驱动"网格→DOF→装配→边界条件施加→求解→后处理"完整闭环的求解器框架。静电场求解器是所有低频电磁求解器中最基础的一个，其架构设计将作为后续静磁场、涡流场、瞬态场求解器的模板。对标ANSYS Maxwell的Electrostatic Solver，需支持2D/3D/轴对称模型、Dirichlet/Neumann/混合边界条件、电压源激励、电容矩阵计算、静电力计算等核心功能。

## What Changes
- **新增** 物理场抽象基类 `PhysicsField`（solver命名空间），定义所有物理场的通用接口
- **新增** 静电场求解器类 `ElectrostaticSolver`，派生自 `PhysicsField`
- **新增** 求解器调度器 `SolverScheduler`，负责求解流程的统一调度
- **新增** 边界条件施加模块 `BoundaryConditionManager`，将tool::Boundary转化为矩阵/向量修改
- **新增** 激励施加模块 `ExcitationManager`，将tool::Excitation转化为右端项
- **新增** 场数据与后处理模块 `FieldData`，存储求解结果并计算派生物理量
- **新增** 轴对称静电场积分权重处理（r-z坐标系）
- **新增** 电容矩阵计算模块 `CapacitanceCalculator`
- **新增** 静电力计算模块 `ElectrostaticForceCalculator`（虚功法）
- **新增** 静电场求解器完整测试用例
- **修改** `SolverApp::run()` 以集成静电场求解流程
- **零修改** 已有底层模块（网格、DOF、形函数、积分器、装配、线性求解器）

## Impact
- Affected specs: 无破坏性修改，纯增量开发
- Affected code:
  - `include/solver/physics_field.hpp` — 物理场抽象基类（新建）
  - `include/solver/electrostatic_solver.hpp` — 静电场求解器（新建）
  - `include/solver/solver_scheduler.hpp` — 求解器调度器（新建）
  - `include/solver/boundary_condition_manager.hpp` — 边界条件管理器（新建）
  - `include/solver/excitation_manager.hpp` — 激励管理器（新建）
  - `include/solver/field_data.hpp` — 场数据与后处理（新建）
  - `include/solver/capacitance_calculator.hpp` — 电容计算器（新建）
  - `include/solver/electrostatic_force_calculator.hpp` — 静电力计算器（新建）
  - `src/solver/` — 所有实现文件（新建目录）
  - `test/test_electrostatic_solver.cpp` — 测试用例（新建）
  - `src/app/solver_app.cpp` — 修改run()方法集成静电场求解流程
  - `src/CMakeLists.txt` — 更新编译配置

## ADDED Requirements

### Requirement: 物理场抽象基类 PhysicsField

系统 SHALL 提供物理场抽象基类 `PhysicsField`，定义所有物理场求解器的通用接口契约。

```cpp
class PhysicsField {
public:
    virtual ~PhysicsField() = default;

    virtual bool setup(const fe_em::EMMeshData& mesh_data,
                       const std::map<int, numeric::MaterialProperties>& materials,
                       const std::vector<tool::Boundary>& boundaries,
                       const std::vector<tool::Excitation>& excitations) = 0;

    virtual bool solve() = 0;

    virtual bool postProcess() = 0;

    virtual const solver::FieldData& getFieldData() const = 0;

    virtual tool::SimulationType getSimulationType() const = 0;

    virtual std::string getSolverName() const = 0;

    virtual void clear() = 0;
};
```

#### Scenario: 派生类接口一致性
- **WHEN** 新的物理场求解器继承 PhysicsField
- **THEN** 必须实现全部纯虚方法，保证多态调用时接口统一

---

### Requirement: 静电场求解器 ElectrostaticSolver

系统 SHALL 提供 `ElectrostaticSolver` 类，实现完整的静电场有限元求解流程，对标Maxwell Electrostatic Solver。

#### 控制方程

- **二维/三维静电场**：∇·(ε∇φ) = -ρ
  - 未知量：电位 φ（标量节点自由度）
  - 材料参数：介电常数 ε（F/m），可各向异性
  - 源项：体电荷密度 ρ（C/m³）

- **轴对称静电场（r-z坐标系）**：1/r · ∂/∂r(r·ε·∂φ/∂r) + ∂/∂z(ε·∂φ/∂z) = -ρ
  - 积分权重需乘以 2πr（轴对称旋转体积分）
  - 边界条件需考虑轴对称几何的特殊性

#### 弱形式

∫Ω ∇N^T · ε · ∇φ dΩ = ∫Ω N · ρ dΩ + ∫Γ N · (ε·∂φ/∂n) dΓ

#### 支持的维度与单元

| 维度 | 单元类型 | DOF类型 |
|------|---------|---------|
| 2D | TRI3, TRI6, QUAD4, QUAD8, QUAD9 | SCALAR_ONLY |
| 3D | TET4, TET10, HEX8, HEX20, HEX27, PRISM6, PRISM15, PYRAMID5, PYRAMID13 | SCALAR_ONLY |
| 轴对称 | TRI3, TRI6, QUAD4, QUAD8, QUAD9 | SCALAR_ONLY（积分权重含2πr因子） |

#### 核心求解流程

```
1. setup()阶段：
   a. 前处理检查（网格有效性、材料参数合理性）
   b. 构建MaterialProperties映射表（从tool::Material提取ε）
   c. DOF分配（EMDOFManager，SCALAR_ONLY类型）
   d. 全局矩阵装配（EMAssembly，调用ElectrostaticIntegrator）
   e. 边界条件施加（BoundaryConditionManager）
   f. 激励施加（ExcitationManager）

2. solve()阶段：
   a. 选择线性求解器（对称正定→SymmetricDirectSolver或CGSolver）
   b. 求解线性方程组 Kφ = F
   c. 解向量回扩（将自由DOF解扩展到包含约束DOF的完整向量）

3. postProcess()阶段：
   a. 计算电场强度 E = -∇φ
   b. 计算电位移矢量 D = εE
   c. 计算静电能量 W = 0.5·φ^T·K·φ
   d. 计算电容矩阵（多导体系统）
   e. 计算静电力（虚功法）
   f. 存储到FieldData
```

#### Scenario: 2D平行板电容器求解
- **GIVEN** 一个2D平行板电容器模型（两块平行导体板，中间为电介质ε_r=4.0）
- **AND** 上板施加Dirichlet边界条件 φ=100V，下板施加 φ=0V
- **WHEN** 调用 setup() → solve() → postProcess()
- **THEN** 求解成功，电位在两板间呈线性分布
- **AND** 电场强度 E ≈ 100/d（d为板间距），方向从高电位指向低电位
- **AND** 电容值与解析解 C = ε₀·ε_r·A/d 的误差 < 5%

#### Scenario: 3D绝缘子电场分析
- **GIVEN** 一个3D绝缘子模型（包含TET4单元，多种电介质材料）
- **AND** 高压端Dirichlet边界 φ=10kV，接地端 φ=0V
- **WHEN** 调用完整求解流程
- **THEN** 求解成功，电场在绝缘子表面分布合理
- **AND** 最大电场强度位置与物理预期一致

#### Scenario: 轴对称同轴电缆电场求解
- **GIVEN** 一个轴对称同轴电缆模型（r-z坐标系，内导体r₁，外导体r₂）
- **AND** 内导体 φ=V₀，外导体 φ=0
- **WHEN** 调用完整求解流程
- **THEN** 电位分布与解析解 φ(r) = V₀·ln(r₂/r)/ln(r₂/r₁) 吻合，误差 < 2%
- **AND** 电场强度 E(r) = V₀/(r·ln(r₂/r₁)) 与解析解吻合

---

### Requirement: 求解器调度器 SolverScheduler

系统 SHALL 提供 `SolverScheduler` 类，负责求解流程的统一调度与管理。

```cpp
class SolverScheduler {
public:
    SolverScheduler();
    ~SolverScheduler();

    bool initialize(tool::ProjectManager& pm);

    bool run();

    void setPhysicsField(std::unique_ptr<PhysicsField> field);

    const solver::FieldData& getFieldData() const;

    void clear();
};
```

#### 调度流程

1. **初始化阶段**：从ProjectManager提取EMMeshData、Material、Boundary、Excitation
2. **物理场创建**：根据SimulationType创建对应的PhysicsField派生类实例
3. **求解执行**：调用PhysicsField::setup() → solve() → postProcess()
4. **结果输出**：从FieldData提取结果，写入输出文件

#### Scenario: 根据SimulationType自动选择求解器
- **WHEN** ProjectManager中的SolutionSetup指定SimulationType=ELECTROSTATIC
- **THEN** SolverScheduler自动创建ElectrostaticSolver实例
- **AND** 完整执行setup→solve→postProcess流程

---

### Requirement: 边界条件管理器 BoundaryConditionManager

系统 SHALL 提供 `BoundaryConditionManager` 类，将 `tool::Boundary` 转化为对系统矩阵和右端向量的修改。

#### 支持的边界条件类型（静电场场景）

| BndType | 物理含义 | 实现方式 |
|---------|---------|---------|
| DIRICHLET | 给定电位值 φ=V₀ | 修改DOF映射表（约束标记）+ 右端项修正 |
| NEUMANN | 给定法向电场分量 ε·∂φ/∂n=g | 自然边界条件（默认g=0）或右端项添加面积分 |
| ROBIN | 阻抗边界条件 aφ + b·∂φ/∂n = c | 矩阵对角线修改 + 右端项修改 |
| PERFECT_E | 理想电壁（等电位面） | 等价于Dirichlet φ=const |
| INSULATION | 绝缘边界（法向D=0） | 等价于Neumann g=0（自然满足） |
| BALLOON | 气球边界（远场近似） | 等价于Neumann g=0（开放边界近似） |
| EVEN_SYMMETRY | 偶对称边界 | 法向E=0 → Neumann g=0 |
| ODD_SYMMETRY | 奇对称边界 | 切向E=0 → Dirichlet φ=0 |
| PERIODIC / ANTIPERIODIC | 周期性边界 | 主从DOF约束（需额外处理） |
| MASTER_SLAVE | 主从边界 | 从DOF值等于主DOF值（约束方程） |

#### Dirichlet边界条件施加方法

采用**消去法（Elimination）**：
1. 在DOF管理阶段，将Dirichlet DOF标记为约束（Local2Global.indices[i] = -1）
2. 装配阶段自动跳过约束DOF的行和列
3. 对右端项进行修正：F_free = F_free - K_free·constrained · φ_constrained
4. 约束DOF的值存储在constrained_dof_values中，求解后回扩到完整解向量

#### Scenario: Dirichlet边界条件正确施加
- **GIVEN** 一个静电场模型，部分节点施加Dirichlet边界条件 φ=100V
- **WHEN** 调用BoundaryConditionManager::applyDirichlet()
- **THEN** 约束DOF在Local2Global映射表中标记为-1
- **AND** 右端项向量正确修正（减去约束DOF对自由DOF的贡献）
- **AND** 求解后约束DOF的值正确回扩为100V

#### Scenario: Neumann边界条件自然满足
- **GIVEN** 一个静电场模型，外边界施加Neumann边界条件 ε·∂φ/∂n=0
- **WHEN** 执行装配和求解
- **THEN** Neumann条件自动满足（有限元弱形式的自然边界条件）
- **AND** 无需额外修改矩阵或右端项

---

### Requirement: 激励管理器 ExcitationManager

系统 SHALL 提供 `ExcitationManager` 类，将 `tool::Excitation` 转化为右端项向量的贡献。

#### 静电场支持的激励类型

| ExcitationType | 物理含义 | 实现方式 |
|---------------|---------|---------|
| VOLTAGE_SOURCE | 电压源（给定电位差） | 转化为Dirichlet边界条件 |
| CURRENT_SOURCE | 电流源（给定电荷注入） | 转化为右端项体积分/面积分 |
| VOLUME_CURRENT | 体电荷密度源 | 右端项添加 ∫N·ρ dΩ |

#### Scenario: 电压源激励转化为Dirichlet边界条件
- **GIVEN** 一个Excitation，类型为VOLTAGE_SOURCE，值为100V，施加于导体表面
- **WHEN** 调用ExcitationManager::processExcitation()
- **THEN** 自动创建对应的Dirichlet边界条件（φ=100V）
- **AND** 导体表面所有节点的DOF被标记为约束

---

### Requirement: 场数据与后处理 FieldData

系统 SHALL 提供 `FieldData` 类，存储求解结果数据并计算派生物理量。

```cpp
class FieldData {
public:
    void setNodalPotential(const Eigen::VectorXd& phi, int num_total_dofs);

    const Eigen::VectorXd& getNodalPotential() const;

    void computeElectricField(const fe_em::EMMeshData& mesh_data,
                              const std::map<int, numeric::MaterialProperties>& materials);

    void computeElectricDisplacement(const fe_em::EMMeshData& mesh_data,
                                     const std::map<int, numeric::MaterialProperties>& materials);

    double computeElectrostaticEnergy(const numeric::CsrMatrix<double>& K,
                                      const Eigen::VectorXd& phi) const;

    const std::vector<Eigen::Vector3d>& getElementElectricField() const;
    const std::vector<Eigen::Vector3d>& getElementElectricDisplacement() const;

    bool exportVTK(const std::string& filename, const fe_em::EMMeshData& mesh_data) const;
    bool exportCSV(const std::string& filename, const fe_em::EMMeshData& mesh_data) const;

    void clear();
};
```

#### 后处理计算公式

- **电场强度**：E = -∇φ（在每个单元的积分点处，通过形函数梯度计算）
- **电位移矢量**：D = ε·E
- **静电能量**：W = 0.5 · φ^T · K · φ
- **体电荷密度**：ρ = ∇·D（散度计算）

#### Scenario: 电场强度计算正确性
- **GIVEN** 求解完成的平行板电容器模型，电位在板间线性分布
- **WHEN** 调用FieldData::computeElectricField()
- **THEN** 每个单元的电场强度E为常数（一阶单元）
- **AND** E的方向从高电位指向低电位
- **AND** |E| ≈ Δφ/d

#### Scenario: VTK格式导出
- **WHEN** 调用FieldData::exportVTK()导出结果
- **THEN** 生成合法的VTK文件，包含节点电位和单元电场强度数据
- **AND** 可被ParaView正确加载和可视化

---

### Requirement: 电容矩阵计算 CapacitanceCalculator

系统 SHALL 提供 `CapacitanceCalculator` 类，计算多导体系统的电容矩阵，对标Maxwell的电容矩阵输出功能。

#### 计算方法

采用**能量法**计算电容矩阵：
1. 对N个导体系统，依次对第i个导体施加单位电压V_i=1V，其余导体接地V_j=0
2. 求解N次静电场问题
3. 计算每次求解的静电能量 W_i = 0.5·φ_i^T·K·φ_i
4. 计算互电容：C_ij = 2·W_ij（当i=j时为自电容）
5. 或采用电荷法：C_ij = Q_i/V_j（通过高斯面积分计算导体电荷）

```cpp
class CapacitanceCalculator {
public:
    Eigen::MatrixXd computeCapacitanceMatrix(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        const std::vector<tool::Boundary>& boundaries,
        const std::vector<int>& conductor_ids,
        numeric::EMLinearSolverBase& solver);

    double computeCapacitance(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        int conductor1_id, int conductor2_id,
        numeric::EMLinearSolverBase& solver);
};
```

#### Scenario: 平行板电容器电容计算
- **GIVEN** 一个2D平行板电容器模型（板面积A，间距d，介电常数ε_r）
- **WHEN** 调用CapacitanceCalculator::computeCapacitance()
- **THEN** 计算的电容值与解析解 C = ε₀·ε_r·A/d 的误差 < 5%

#### Scenario: 三导体系统电容矩阵
- **GIVEN** 一个包含3个导体的系统
- **WHEN** 调用CapacitanceCalculator::computeCapacitanceMatrix()
- **THEN** 返回3×3对称电容矩阵
- **AND** 对角线元素为自电容（正值），非对角线元素为互电容（负值）
- **AND** 每行元素之和等于该导体对地的总电容

---

### Requirement: 静电力计算 ElectrostaticForceCalculator

系统 SHALL 提供 `ElectrostaticForceCalculator` 类，采用虚功法计算导体上的静电力，对标Maxwell的力计算功能。

#### 计算方法

采用**虚功法（Virtual Work Method）**：
- 对导体表面进行虚位移δs
- 计算能量对位移的偏导数：F = -∂W/∂s
- 数值实现：在导体表面两侧各做微小位移，计算能量差分

```cpp
class ElectrostaticForceCalculator {
public:
    Eigen::Vector3d computeForce(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        int conductor_id,
        numeric::EMLinearSolverBase& solver,
        double delta = 1e-6);

    double computeTorque(
        const fe_em::EMMeshData& mesh_data,
        const std::map<int, numeric::MaterialProperties>& materials,
        int conductor_id,
        const Eigen::Vector3d& axis_point,
        const Eigen::Vector3d& axis_direction,
        numeric::EMLinearSolverBase& solver,
        double delta = 1e-6);
};
```

#### Scenario: 平行板电容器静电力计算
- **GIVEN** 一个平行板电容器模型（板面积A，间距d，电压V）
- **WHEN** 调用ElectrostaticForceCalculator::computeForce()
- **THEN** 计算的吸引力与解析解 F = ε₀·ε_r·A·V²/(2d²) 的误差 < 10%

---

### Requirement: 轴对称模型特殊处理

系统 SHALL 支持轴对称静电场模型的求解，需在积分层面进行特殊处理。

#### 轴对称控制方程

在r-z坐标系下，静电场方程为：
1/r · ∂/∂r(r·ε·∂φ/∂r) + ∂/∂z(ε·∂φ/∂z) = -ρ

#### 弱形式修改

∫Ω ∇N^T · ε · ∇φ · r dΩ = ∫Ω N · ρ · r dΩ

关键修改点：
1. **积分权重**：每个积分点需乘以 2πr_i（r_i为积分点的径向坐标），将截面积分转化为旋转体积分
2. **刚度矩阵**：K_e = ∫∇N^T · ε · ∇N · 2πr · detJ · w dξ
3. **源项向量**：F_e = ∫N · ρ · 2πr · detJ · w dξ
4. **轴上处理**：r=0时需特殊处理（L'Hôpital规则），避免1/r奇异性
5. **边界条件**：对称轴上自然满足 φ/r 的正则性条件

#### 实现策略

在 `ElectrostaticIntegrator` 中增加轴对称模式标志：
- 当 `DimType::AXIS` 时，积分权重自动乘以 2πr
- 形函数和梯度计算不变（仍使用2D参考单元）
- 仅影响积分权重，不改变单元拓扑

#### Scenario: 轴对称模式积分权重正确性
- **GIVEN** 一个轴对称模型（DimType::AXIS），TRI3单元
- **WHEN** 计算单元刚度矩阵
- **THEN** 每个积分点的权重乘以 2πr_i（r_i为积分点径向坐标）
- **AND** r=0的积分点使用L'Hôpital规则处理奇异性

---

### Requirement: SolverApp集成

系统 SHALL 修改 `SolverApp::run()` 方法，集成静电场求解流程。

#### 集成方案

```cpp
bool SolverApp::run() {
    auto& pm = tool::ProjectManager::getInstance();

    solver::SolverScheduler scheduler;
    if (!scheduler.initialize(pm)) return false;
    if (!scheduler.run()) return false;

    const auto& field_data = scheduler.getFieldData();
    // 导出结果...
    return true;
}
```

#### Scenario: 从JSON配置文件运行静电场求解
- **GIVEN** 一个JSON配置文件，指定SimulationType=ELECTROSTATIC
- **WHEN** 调用SolverApp::initialize() → run()
- **THEN** 自动创建ElectrostaticSolver，执行完整求解流程
- **AND** 结果输出到指定目录

---

### Requirement: 错误处理与鲁棒性

系统 SHALL 实现分级错误处理机制，确保静电场求解器的工业级鲁棒性。

#### 前处理检查项

1. **网格有效性**：节点坐标非NaN、单元节点ID合法、单元体积/面积>0
2. **材料参数合理性**：介电常数ε > 0、无未分配材料的单元
3. **边界条件完整性**：至少存在一个Dirichlet边界条件（否则系统奇异）
4. **边界条件冲突检测**：同一DOF不能同时施加多个不同值的Dirichlet条件
5. **激励合法性**：电压源值非NaN、电流源值合理

#### 运行时错误处理

1. 矩阵奇异检测（求解器返回NUMERICAL_ERROR时提供诊断信息）
2. 非正定矩阵检测（静电场刚度矩阵必须对称正定）
3. 内存不足时的优雅退出

#### Scenario: 缺少Dirichlet边界条件时报警
- **GIVEN** 一个静电场模型，未施加任何Dirichlet边界条件
- **WHEN** 调用ElectrostaticSolver::setup()
- **THEN** 输出ERROR日志："No Dirichlet boundary condition found, system is singular"
- **AND** setup()返回false

#### Scenario: 边界条件冲突检测
- **GIVEN** 同一节点同时施加两个不同值的Dirichlet边界条件（φ=0V和φ=100V）
- **WHEN** 调用BoundaryConditionManager::applyDirichlet()
- **THEN** 输出WARNING日志："Conflicting Dirichlet BC on DOF X: value1=0, value2=100"
- **AND** 使用后施加的值（或优先级更高的值）

---

### Requirement: 测试用例

系统 SHALL 提供完整的静电场求解器测试用例，使用Google Test框架。

#### 测试用例清单

1. **2D平行板电容器**：验证电位线性分布、电场均匀、电容计算精度
2. **2D同轴电缆**：验证圆柱坐标解析解对比
3. **3D立方体电容器**：验证3D求解流程正确性
4. **轴对称同轴电缆**：验证轴对称积分权重处理
5. **多导体电容矩阵**：验证3导体系统电容矩阵对称性和正确性
6. **边界条件施加**：验证Dirichlet/Neumann/混合边界条件
7. **静电力计算**：验证虚功法力计算精度
8. **错误处理**：验证缺少边界条件、材料参数异常等场景

#### Scenario: 2D平行板电容器测试通过与Maxwell对比
- **WHEN** 运行test_electrostatic_solver中的平行板电容器测试
- **THEN** 电位分布误差 < 1%（与解析解对比）
- **AND** 电容计算误差 < 5%
- **AND** 电场强度误差 < 5%

## MODIFIED Requirements

无（纯新增模块，不修改已有代码逻辑，仅扩展SolverApp::run()的调用路径）

## REMOVED Requirements

无

## 技术约束

1. **语言标准**：C++17及以上，禁止使用C++20特性
2. **核心依赖**：STL + Eigen 3.4+（必需）
3. **可选依赖**：SuperLU/MUMPS（条件编译）
4. **内存管理**：禁止裸指针管理动态资源，优先使用智能指针
5. **线程安全**：无全局变量，所有状态封装在实例对象内
6. **日志规范**：使用FEEM_INFO/FEEM_DEBUG/FEEM_ERROR等宏输出日志
7. **头文件规范**：仅包含文件名，不包含路径（由CMakeLists.txt管理路径）
8. **命名规范**：类名PascalCase、函数camelCase、变量snake_case、常量UPPER_SNAKE_CASE
9. **注释要求**：关键逻辑必须注释，使用Doxygen风格文档注释
10. **输入校验**：所有公开接口必须进行参数合法性检查
11. **禁止硬编码**：所有可配置参数支持JSON配置
12. **禁止代码冗余**：重复逻辑封装为工具函数
