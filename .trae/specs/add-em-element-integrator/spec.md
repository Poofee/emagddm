# 电磁场单元矩阵计算模块（EM Element Integrator Module）Spec

## Why
形函数模块已完成，但缺乏单元级矩阵组装能力。当前项目无法进行刚度矩阵、质量矩阵、阻尼矩阵的数值积分计算，无法支撑静态/瞬态电磁场有限元求解。本模块将基于已有的 Lagrange 标量节点基（支持全部 **19 种**类型，含 PRISM6/PRISM15/PYRAMID5/PYRAMID13）和 Nedelec 一阶棱边基（支持全部 **6 种**类型，含 PRISM6_EDGE/PYRAMID5_EDGE），实现完整的单元矩阵计算功能，为后续全局系统组装和求解器提供核心数据。

## What Changes
- 在 `include/numeric/` 和 `src/numeric/` 下新增电磁场积分模块全套文件
- 新增抽象基类 `EMElementIntegratorBase`，定义统一接口契约
- 新增高斯积分点库 `GaussQuadrature`，支持 2D/3D 全部常用单元类型（**含 PRISM 和 PYRAMID**）
- 新增二维/三维静电场/瞬态电场积分器 `ElectrostaticIntegrator`（支持 **8 种 3D 单元**）
- 新增二维静磁场/瞬态磁场标量位积分器 `MagneticScalar2DIntegrator`
- 新增三维静磁场/瞬态磁场矢量位积分器 `MagneticVector3DIntegrator`（支持 **4 种 3D 棱单元**）
- 新增测试用例文件 `test_em_element_integrator.cpp`，验证对称性和 H(curl) 共形性
- 更新 CMakeLists.txt 编译配置

## Impact
- Affected specs: 无（全新模块，依赖已完成的形函数模块）
- Affected code:
  - `include/numeric/em_element_integrator_base.hpp` — 抽象基类
  - `include/numeric/gauss_quadrature.hpp` — 高斯积分点库
  - `include/numeric/electrostatic_integrator.hpp` — 静电场/瞬态电场积分器声明
  - `include/numeric/magnetic_scalar_2d_integrator.hpp` — 二维静磁场/瞬态磁场标量位积分器声明
  - `include/numeric/magnetic_vector_3d_integrator.hpp` — 三维静磁场/瞬态磁场矢量位积分器声明
  - `src/numeric/gauss_quadrature.cpp` — 高斯积分点库实现
  - `src/numeric/electrostatic_integrator.cpp` — 静电场/瞬态电场积分器实现
  - `src/numeric/magnetic_scalar_2d_integrator.cpp` — 二维磁标位积分器实现
  - `src/numeric/magnetic_vector_3d_integrator.cpp` — 三维磁矢量位积分器实现
  - `test/test_em_element_integrator.cpp` — 测试用例
  - `src/CMakeLists.txt` / `test/CMakeLists.txt` — 编译配置更新

## ADDED Requirements

### Requirement: 高斯积分点库（GaussQuadrature）

系统 SHALL 提供高斯积分点库 `GaussQuadrature`，支持工程常用一阶积分精度：

**2D 单元类型：**
| 单元类型 | 积分点数 | 积分点位置 | 权重 |
|---------|---------|-----------|------|
| TRI3 | 1 | 形心 (1/3, 1/3) | 1/2 |
| TRI3 | 3 | 三边中点 | 1/6 |
| QUAD4 | 4 (2×2) | (±1/√3, ±1/√3) | 1 |

**3D 单元类型：**
| 单元类型 | 积分点数 | 积分点位置 | 权重 |
|---------|---------|-----------|------|
| TET4 | 1 | 形心 (1/4, 1/4, 1/4) | 1/6 |
| TET4 | 4 | 四面体内部点 | 1/24 |
| HEX8 | 8 (2×2×2) | (±1/√3, ±1/√3, ±1/√3) | 1 |
| **PRISM6** | **6** | **棱柱内部点** | **1/6** |
| **PYRAMID5** | **5** | **金字塔内部点** | **变权重** |

#### Scenario: 获取 TRI3 一阶积分点
- **WHEN** 调用 `GaussQuadrature::getPoints(ElementType::TRI3, 1)`
- **THEN** 返回包含 1 个积分点的列表，坐标为 (1/3, 1/3)，权重为 0.5

#### Scenario: 获取 PRISM6 一阶积分点
- **WHEN** 调用 `GaussQuadrature::getPoints(ElementType::PRISM6, 6)`
- **THEN** 返回包含 6 个积分点的列表，覆盖三棱柱体积域

#### Scenario: 获取 PYRAMID5 一阶积分点
- **WHEN** 调用 `GaussQuadrature::getPoints(ElementType::PYRAMID5, 5)`
- **THEN** 返回包含 5 个积分点的列表，覆盖金字塔体积域

### Requirement: 电磁场积分器抽象基类（EMElementIntegratorBase）

系统 SHALL 提供纯虚抽象基类 `EMElementIntegratorBase`，统一所有电磁场景的接口契约：

| 接口方法 | 功能说明 |
|---------|---------|
| `computeStiffnessMatrix()` | 计算单元刚度矩阵 K（静态/瞬态共用） |
| `computeMassMatrix()` | 计算单元质量矩阵 M（瞬态专用） |
| `computeDampingMatrix()` | 计算单元阻尼矩阵 C（瞬态专用） |
| `computeSourceVector()` | 计算单元源项向量 F |
| `setMaterialProperties()` | 设置材料参数（ε, μ, σ） |
| `isTransient()` | 查询是否为瞬态模式 |

#### Scenario: 基类接口一致性
- **WHEN** 派生类继承 EMElementIntegratorBase
- **THEN** 必须实现全部纯虚方法，保证多态调用时接口统一

### Requirement: 二维/三维静电场/瞬态电场积分器（ElectrostaticIntegrator）

系统 SHALL 提供 `ElectrostaticIntegrator` 类，支持二维和三维静电场/瞬态电场的单元矩阵计算：

**控制方程：**
- 静态：∇·(ε∇φ_e) = -ρ_e
- 瞬态：∇·(ε∇φ_e) + σ∂φ_e/∂t = -ρ_e

**弱形式（瞬态）：**
∫∇N^T·ε·∇N φ_e dV + ∫N^T·σ·N ∂φ_e/∂t dV = ∫N·ρ_e dV

**支持的单元类型：**
- 2D：TRI3、QUAD4（使用 ScalarShapeFunctionBase）
- 3D：**TET4、HEX8、PRISM6、PRISM15、PYRAMID5、PYRAMID13**（使用 ScalarShapeFunctionBase）

**材料参数：**
- 线性各向同性介电常数 ε（F/m）
- 电导率 σ（S/m）
- 预留非线性虚函数接口

**计算公式：**
- 刚度矩阵：K_e = ∫∇N^T · ε · ∇N · detJ · w
- 阻尼矩阵（瞬态）：C_e = ∫N^T · σ · N · detJ · w
- 源项向量：F_e = ∫N · ρ_e · detJ · w（先实现零源项）

#### Scenario: TRI3 单元静电场刚度矩阵计算
- **WHEN** 给定 TRI3 单元的节点坐标和材料参数 ε=8.854e-12 F/m
- **THEN** 返回 3×3 对称正定的刚度矩阵 K_e，满足 K_e = K_e^T 且所有特征值 > 0

#### Scenario: PRISM6 单元静电场刚度矩阵计算
- **WHEN** 给定 PRISM6 单元的节点坐标和材料参数 ε=8.854e-12 F/m
- **THEN** 返回 6×6 对称正定的刚度矩阵 K_e，满足 K_e = K_e^T 且所有特征值 > 0

#### Scenario: PYRAMID5 单元静电场刚度矩阵计算
- **WHEN** 给定 PYRAMID5 单元的节点坐标和材料参数 ε=8.854e-12 F/m
- **THEN** 返回 5×5 对称正定的刚度矩阵 K_e，满足 K_e = K_e^T 且所有特征值 > 0

#### Scenario: QUAD4 单元瞬态电场阻尼矩阵计算
- **WHEN** 给定 QUAD4 单元的节点坐标、ε 和 σ=1e6 S/m
- **THEN** 返回 4×4 对称半正定的阻尼矩阵 C_e，满足 C_e = C_e^T

### Requirement: 二维静磁场/瞬态磁场标量位积分器（MagneticScalar2DIntegrator）

系统 SHALL 提供 `MagneticScalar2DIntegrator` 类，支持二维静磁场/瞬态磁场的标量磁位 φ_m 计算：

**控制方程：**
- 静态：∇·(μ∇φ_m) = -ρ_m
- 瞬态（磁扩散）：∇·(μ∇φ_m) + σ∂φ_m/∂t = -ρ_m

**弱形式（瞬态）：**
∫∇N^T·μ·∇N φ_m dV + ∫N^T·σ·N ∂φ_m/∂t dV = ∫N·ρ_m dV

**支持的单元类型：**
- 2D：TRI3、QUAD4（使用 ScalarShapeFunctionBase）

**材料参数：**
- 线性各向同性磁导率 μ（H/m）
- 电导率 σ（S/m）
- 预留非线性虚函数接口（用于非线性 B-H 曲线）

**计算公式：**
- 刚度矩阵：K_m = ∫∇N^T · μ · ∇N · detJ · w
- 质量矩阵（瞬态）：M_m = ∫N^T · σ · N · detJ · w（对应磁扩散项）
- 源项向量：F_m（先实现零源项）

#### Scenario: TRI3 单元静磁场刚度矩阵计算
- **WHEN** 给定 TRI3 单元的节点坐标和材料参数 μ=4π×10^-7 H/m
- **THEN** 返回 3×3 对称正定的刚度矩阵 K_m，满足 K_m = K_m^T

#### Scenario: TRI3 单元瞬态磁场质量矩阵计算
- **WHEN** 给定 TRI3 单元的节点坐标、μ 和 σ=5.8e7 S/m（铜的电导率）
- **THEN** 返回 3×3 对称半正定的质量矩阵 M_m，满足 M_m = M_m^T

### Requirement: 三维静磁场/瞬态磁场矢量位积分器（MagneticVector3DIntegrator）

系统 SHALL 提供 `MagneticVector3DIntegrator` 类，支持三维静磁场/瞬态磁场的矢量磁位 A 计算（H(curl) 弱形式）：

**控制方程（忽略位移电流）：**
- 静态：∇×(ν∇×A) = J_s，其中 ν=1/μ 为磁阻率
- 瞬态（涡流场 A-V 简化）：∇×(ν∇×A) + σ∂A/∂t = J_s

**弱形式（瞬态）：**
∫(∇×N_edge)^T·ν·(∇×A) dV + ∫N_edge^T·σ·N_edge ∂A/∂t dV = ∫N_edge·J_s dV

**支持的单元类型：**
- 3D：**TET4_EDGE**（6 个棱边自由度）、**HEX8_EDGE**（12 个棱边自由度）、**PRISM6_EDGE**（9 个棱边自由度）、**PYRAMID5_EDGE**（8 个棱边自由度）（使用 NedelecEdgeShapeFunctionBase）

**材料参数：**
- 线性各向同性磁阻率 ν=1/μ（m/H）
- 电导率 σ（S/m）
- 预留非线性虚函数接口（用于非线性 ν(B)）

**计算公式：**
- 旋度刚度矩阵：K_curl = ∫(∇×N_edge)^T · ν · (∇×N_edge) · detJ · w
- 质量矩阵（瞬态）：M_A = ∫N_edge^T · σ · N_edge · detJ · w（对应涡流项）
- 源项向量：F_J = ∫N_edge·J_s · detJ · w（先实现零源项）

#### Scenario: TET4_EDGE 单元静磁场旋度刚度矩阵计算
- **WHEN** 给定 TET4_EDGE 单元的节点坐标和材料参数 ν=1/(4π×10^-7) m/H
- **THEN** 返回 6×6 的旋度刚度矩阵 K_curl，满足 K_curl = K_curl^T（对称性验证）

#### Scenario: PRISM6_EDGE 单元静磁场旋度刚度矩阵计算
- **WHEN** 给定 PRISM6_EDGE 单元的节点坐标和材料参数 ν=1/(4π×10^-7) m/H
- **THEN** 返回 9×9 的旋度刚度矩阵 K_curl，满足 K_curl = K_curl^T（对称性验证）

#### Scenario: PYRAMID5_EDGE 单元静磁场旋度刚度矩阵计算
- **WHEN** 给定 PYRAMID5_EDGE 单元的节点坐标和材料参数 ν=1/(4π×10^-7) m/H
- **THEN** 返回 8×8 的旋度刚度矩阵 K_curl，满足 K_curl = K_curl^T（对称性验证）

#### Scenario: TET4_EDGE 单元瞬态磁场质量矩阵计算
- **WHEN** 给定 TET4_EDGE 单元的节点坐标、ν 和 σ=5.8e7 S/m
- **THEN** 返回 6×6 对称半正定的质量矩阵 M_A，满足 M_A = M_A^T

#### Scenario: PRISM6_EDGE 单元 H(curl) 共形性验证
- **WHEN** 在共享棱边上评估相邻两个 PRISM6_EDGE 单元的基函数
- **THEN** 切向分量连续（H(curl) 协调性验证）

#### Scenario: PYRAMID5_EDGE 单元 H(curl) 共形性验证
- **WHEN** 在共享棱边上评估相邻两个 PYRAMID5_EDGE 单元的基函数
- **THEN** 切向分量连续（H(curl) 协调性验证）

### Requirement: 测试用例

系统 SHALL 提供可运行的测试程序 `test_em_element_integrator.cpp`，使用 Google Test 框架，对以下内容进行单单元验证：

1. **高斯积分点测试：**
   - TRI3 一阶/三阶积分点坐标与权重验证
   - QUAD4 二阶积分点验证
   - TET4 一阶/四阶积分点验证
   - HEX8 二阶积分点验证
   - **PRISM6 一阶积分点验证**
   - **PYRAMID5 一阶积分点验证**

2. **静电场/瞬态电场测试：**
   - TRI3 单元 K_e 对称性与正定性验证
   - QUAD4 单元 K_e、C_e 对称性验证
   - TET4 单元 K_e 对称性与正定性验证
   - HEX8 单元 K_e、C_e 对称性验证
   - **PRISM6 单元 K_e 对称性与正定性验证**
   - **PRISM15 单元 K_e 对称性验证**
   - **PYRAMID5 单元 K_e 对称性与正定性验证**
   - **PYRAMID13 单元 K_e 对称性验证**

3. **二维静磁场/瞬态磁场测试：**
   - TRI3 单元 K_m、M_m 对称性验证
   - QUAD4 单元 K_m、M_m 对称性验证

4. **三维静磁场/瞬态磁场测试：**
   - TET4_EDGE 单元 K_curl、M_A 对称性验证
   - HEX8_EDGE 单元 K_curl、M_A 对称性验证
   - **PRISM6_EDGE 单元 K_curl、M_A 对称性验证（9×9 矩阵）**
   - **PYRAMID5_EDGE 单元 K_curl、M_A 对称性验证（8×8 矩阵）**
   - **H(curl) 共形性验证（TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE 切向连续性检验）**

5. **材料参数设置测试：**
   - 各向同性材料参数正确应用验证
   - 非线性接口预留验证（默认返回线性值）

## MODIFIED Requirements
无（全新模块，不修改现有代码）。

## REMOVED Requirements
无。
