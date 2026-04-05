# Tasks

- [x] Task 1: 创建高斯积分点库 `GaussQuadrature`
  - [x] Task 1.1: 定义积分点数据结构（GaussPoint：坐标 + 权重）
  - [x] Task 1.2: 实现 2D 单元积分点（TRI3 一阶/三阶、QUAD4 二阶）
  - [x] Task 1.3: 实现 3D 单元积分点（TET4 一阶/四阶、HEX8 二阶、**PRISM6 一阶**、**PYRAMID5 一阶**）
  - [x] Task 1.4: 实现工厂方法 `getPoints(element_type, order)` 返回积分点列表

- [x] Task 2: 创建电磁场积分器抽象基类 `EMElementIntegratorBase`
  - [x] Task 2.1: 定义材料参数结构体（MaterialProperties：ε, μ, σ, 非线性标志）
  - [x] Task 2.2: 定义纯虚接口方法（computeStiffnessMatrix/computeMassMatrix/computeDampingMatrix/computeSourceVector/setMaterialProperties/isTransient）
  - [x] Task 2.3: 定义单元矩阵结果结构体（ElementMatrices：K, M, C, F）

- [x] Task 3: 实现二维/三维静电场/瞬态电场积分器 `ElectrostaticIntegrator`
  - [x] Task 3.1: 继承 EMElementIntegratorBase，支持 TRI3/QUAD4/**TET4/HEX8/PRISM6/PRISM15/PYRAMID5/PYRAMID13** 全部单元类型
  - [x] Task 3.2: 实现刚度矩阵计算 K_e = ∫∇N^T · ε · ∇N · detJ · w
  - [x] Task 3.3: 实现阻尼矩阵计算 C_e = ∫N^T · σ · N · detJ · w（瞬态模式）
  - [x] Task 3.4: 实现源项向量计算 F_e = ∫N · ρ_e · detJ · w（零源项实现）
  - [x] Task 3.5: 预留非线性材料接口（虚函数，默认返回线性 ε）

- [x] Task 4: 实现二维静磁场/瞬态磁场标量位积分器 `MagneticScalar2DIntegrator`
  - [x] Task 4.1: 继承 EMElementIntegratorBase，支持 TRI3/QUAD4 单元类型（仅 2D）
  - [x] Task 4.2: 实现刚度矩阵计算 K_m = ∫∇N^T · μ · ∇N · detJ · w
  - [x] Task 4.3: 实现质量矩阵计算 M_m = ∫N^T · σ · N · detJ · w（瞬态磁扩散项）
  - [x] Task 4.4: 实现源项向量计算 F_m（零源项实现）
  - [x] Task 4.5: 预留非线性材料接口（虚函数，默认返回线性 μ）

- [x] Task 5: 实现三维静磁场/瞬态磁场矢量位积分器 `MagneticVector3DIntegrator`
  - [x] Task 5.1: 继承 EMElementIntegratorBase，支持 **TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE** 全部 3D Nedelec 棱单元类型
  - [x] Task 5.2: 实现旋度刚度矩阵计算 K_curl = ∫(∇×N_edge)^T · ν · (∇×N_edge) · detJ · w
  - [x] Task 5.3: 实现质量矩阵计算 M_A = ∫N_edge^T · σ · N_edge · detJ · w（涡流项）
  - [x] Task 5.4: 实现源项向量计算 F_J = ∫N_edge·J_s · detJ · w（零源项实现）
  - [x] Task 5.5: 预留非线性材料接口（虚函数，默认返回线性 ν=1/μ）

- [x] Task 6: 编写测试用例并配置编译
  - [x] Task 6.1: 编写高斯积分点测试（TRI3/QUAD4/TET4/HEX8/**PRISM6**/**PYRAMID5** 积分点坐标与权重验证）
  - [x] Task 6.2: 编写静电场/瞬态电场测试（TRI3/QUAD4/TET4/HEX8/**PRISM6/PRISM15/PYRAMID5/PYRAMID13** 的 K_e、C_e 对称性验证）
  - [x] Task 6.3: 编写二维静磁场/瞬态磁场测试（TRI3/QUAD4 的 K_m、M_m 对称性验证）
  - [x] Task 6.4: 编写三维静磁场/瞬态磁场测试（**TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE** 的 K_curl、M_A 对称性与 H(curl) 共形性验证）
  - [x] Task 6.5: 更新 src/CMakeLists.txt 和 test/CMakeLists.txt 编译配置

# Task Dependencies
- [Task 2] depends on [Task 1] — 积分器基类依赖高斯积分点库
- [Task 3] depends on [Task 1, Task 2] — 静电场积分器依赖积分点和基类接口
- [Task 4] depends on [Task 1, Task 2] — 二维磁标位积分器依赖积分点和基类接口
- [Task 5] depends on [Task 1, Task 2] — 三维磁矢量位积分器依赖积分点和基类接口
- [Task 6] depends on [Task 3, Task 4, Task 5] — 测试需要全部积分器实现完成
