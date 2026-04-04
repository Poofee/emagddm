# 工业级有限元形函数计算模块（Shape Function Module）Spec

## Why
自研电磁有限元求解器（FETI-DP）缺乏统一的形函数计算基础设施。当前项目已有网格管理（EMMesh）和稀疏矩阵模块，但单元级形函数计算、雅可比变换、物理空间梯度映射等核心数值能力尚未实现，无法支撑刚度矩阵组装、载荷向量构造等下游求解流程。

## What Changes
- 在 `include/numeric/` 和 `src/numeric/` 下新增形函数计算模块全套文件
- 新增抽象基类定义纯虚接口（ShapeFunctionBase）
- 新增 Lagrange 标量节点单元完整实现（1D/2D/3D 共 **19 种**单元类型）
- 新增 Nedelec 1阶棱边矢量单元实现（2D/3D 共 **8 种**单元类型）
- 新增单元工厂类（ShapeFunctionFactory），支持字符串自动创建
- 新增测试用例文件（test_shape_function.cpp），单积分点验证
- 更新 CMakeLists.txt 编译配置

## Impact
- Affected specs: 无（全新模块）
- Affected code:
  - `include/numeric/shape_function_base.hpp` — 抽象基类
  - `include/numeric/lagrange_element.hpp` — Lagrange 单元声明
  - `include/numeric/nedelec_element.hpp` — Nedelec 棱单元声明
  - `include/numeric/shape_function_factory.hpp` — 工厂类声明
  - `src/numeric/lagrange_element.cpp` — Lagrange 单元实现
  - `src/numeric/nedelec_element.cpp` — Nedelec 棱单元实现
  - `src/numeric/shape_function_factory.cpp` — 工厂类实现
  - `test/test_shape_function.cpp` — 测试用例
  - `src/CMakeLists.txt` / `test/CMakeLists.txt` — 编译配置更新

## ADDED Requirements

### Requirement: 形函数抽象基类（ShapeFunctionBase）

系统 SHALL 提供纯虚抽象基类 `ShapeFunctionBase`，统一所有单元类型的接口契约：

| 接口方法 | 功能说明 |
|---------|---------|
| `getNodeType()` | 返回单元类型枚举（ElementType） |
| `getNodeCount()` | 返回节点/自由度数量 |
| `getDim()` | 返回单元维度（1/2/3） |
| `evalN(xi)` | 计算局部坐标下的形函数值向量 N |
| `evalGradN(xi)` | 计算形函数对局部坐标的偏导数矩阵 dN/dξ |
| `calcJacobian(xi, node_coords)` | 计算雅可比矩阵 J、detJ、invJ |
| `calcPhysicalGradN(xi, node_coords)` | 计算物理坐标系下的梯度 ∇N = invJ · dN/dξ |

#### Scenario: 基类接口一致性
- **WHEN** 派生类继承 ShapeFunctionBase
- **THEN** 必须实现全部纯虚方法，保证多态调用时接口统一

### Requirement: Lagrange 标量节点单元

系统 SHALL 提供完整的 Lagrange 标量节点单元实现类 `LagrangeElement<T>`，模板参数 T 为单元维度（支持 1D/2D/3D），覆盖以下 **19 种单元类型**：

**1D 单元（2种）：**
- `LINE2`：2节点线性单元，N₁=(1-ξ)/2, N₂=(1+ξ)/2
- `LINE3`：3节点二次单元，N₁=ξ(ξ-1)/2, N₂=1-ξ², N₃=ξ(ξ+1)/2

**2D 单元（5种）：**
- `TRI3`：3节点线性三角形，面积坐标插值
- `TRI6`：6节点二次三角形，含边中节点
- `QUAD4`：4节点双线性四边形，双线性拉格朗日插值
- `QUAD8`：8节点Serendipity四边形（缺中心节点）
- `QUAD9`：9节点双二次四边形（含中心节点）

**3D 单元（12种）：**
- `TET4`：4节点线性四面体，体积坐标插值
- `TET10`：10节点二次四面体
- `HEX8`：8节点三线性六面体
- `HEX20`：20节点Serendipity六面体
- `HEX27`：27节点三二次六面体
- **`PRISM6`**：6节点线性三棱柱（楔形体），三角形底面 × ζ方向线性插值
- **`PRISM15`**：15节点二次三棱柱，含边中节点与面心节点
- **`PYRAMID5`**：5节点线性金字塔，四边形底面 + 锥顶，退化等参形函数
- **`PYRAMID13`**：13节点二次金字塔，含边中节点与面心节点

每个 Lagrange 单元 SHALL 实现：
1. 局部坐标形函数 N(ξ) 的精确解析公式
2. 形函数对局部坐标的偏导数 ∂N/∂ξᵢ 解析公式
3. 通过 `ShapeFunctionBase` 接口提供雅可比变换与物理梯度计算

> **三棱柱（Prism）形函数要点**：
> - 局部坐标：(ξ, η, ζ)，其中 ξ≥0, η≥0, ξ+η≤1, -1≤ζ≤1
> - 线性棱柱 PRISM6：底面用 TRI3 面积坐标，ζ 方向用线性插值，Nᵢ = Lᵢ(ξ,η)·(1-ζ)/2（底面节点），Nᵢ = Lᵢ(ξ,η)·(1+ζ)/2（顶面节点）
>
> **金字塔（Pyramid）形函数要点**：
> - 局部坐标：(ξ, η, ζ)，底面为 ξ∈[-1,1], η∈[-1,1]，锥顶在 ζ=1
> - 线性金字塔 PYRAMID5：采用退化等参映射，N₁~N₄ 为底面角点，N₅ 为锥顶
> - 形函数需处理 ζ→1 时的奇异性（0/0 型），使用有理函数或内部变量替换

#### Scenario: TRI3 单元形函数验证
- **WHEN** 给定 TRI3 单元的三个物理节点坐标和局部坐标 ξ=(1/3, 1/3)（形心）
- **THEN** N = [1/3, 1/3, 1/3]，满足单位分解性 ΣNᵢ = 1

#### Scenario: QUAD4 单元形函数验证
- **WHEN** 给定 QUAD4 单元在 ξ=(0,0) 处求值
- **THEN** N = [0.25, 0.25, 0.25, 0.25]，四个角点贡献相等

#### Scenario: PRISM6 单元形函数验证
- **WHEN** 给定 PRISM6 单元在几何中心 ξ=(1/3, 1/3, 0) 处求值
- **THEN** 底面三个节点各贡献 1/6，顶面三个节点各贡献 1/6，ΣNᵢ = 1

#### Scenario: PYRAMID5 单元形函数验证
- **WHEN** 给定 PYRAMID5 单元在 ξ=(0, 0, 0.25)（靠近底面中心）处求值
- **THEN** 锥顶节点贡献较小，四个底面节点贡献较大，ΣNᵢ = 1

### Requirement: Nedelec 1阶棱边矢量单元

系统 SHALL 提供 Nedelec 第一类（H(curl) 协调）1阶棱边矢量单元实现类 `NedelecElement<T>`，覆盖以下 **8 种单元类型**：

**2D 棱单元（2种）：**
- `TRI3_EDGE`：对应 TRI3 三角形的 3 条棱边，每条棱边 1 个自由度
- `QUAD4_EDGE`：对应 QUAD4 四边形的 4 条棱边，每条棱边 1 个自由度

**3D 棱单元（6种）：**
- `TET4_EDGE`：对应 TET4 四面体的 6 条棱边，每条棱边 1 个自由度
- `HEX8_EDGE`：对应 HEX8 六面体的 12 条棱边，每条棱边 1 个自由度
- **`PRISM6_EDGE`**：对应 PRISM6 三棱柱的 **9 条棱边**（3条底棱 + 3条顶棱 + 3条侧棱），每条棱边 1 个自由度
- **`PYRAMID5_EDGE`**：对应 PYRAMID5 金字塔的 **8 条棱边**（4条底棱 + 4条侧棱），每条棱边 1 个自由度

每个 Nedelec 单元 SHALL 额外实现：
1. **棱边基函数** N_edge(ξ)：矢量值函数，沿棱边切向具有单位环流
2. **旋度 curl(N_edge)**：局部坐标系下的旋度（标量/矢量）
3. **物理坐标系映射**：
   - 矢量基函数物理映射：N_phys = invJ · N_edge / detJ （Piola 变换）
   - 旋度物理映射：curl_phys = curl(N_edge) / detJ

> **三棱柱棱单元（PRISM6_EDGE）基函数要点**：
> - 棱边编号：底面 3 条（1-2, 2-3, 3-1）、顶面 3 条（4-5, 5-6, 6-4）、侧棱 3 条（1-4, 2-5, 3-6）
> - 基函数由对应的楔形体矢量形状函数构造，保证切向连续性
>
> **金字塔棱单元（PYRAMID5_EDGE）基函数要点**：
> - 棱边编号：底面 4 条（1-2, 2-3, 3-4, 4-1）、侧棱 4 条（1-5, 2-5, 3-5, 4-5）
> - 金字塔为退化单元，基函数需特殊处理锥顶处的奇异性

#### Scenario: TRI3_EDGE 基函数切向连续性
- **WHEN** 在共享棱边上评估相邻两个 TRI3_EDGE 单元的基函数
- **THEN** 切向分量连续（H(curl) 协调性）

#### Scenario: TET4_EDGE 自由度数验证
- **WHEN** 查询 TET4_EDGE 单元的自由度数量
- **THEN** 返回 6（四面体有 6 条棱边）

#### Scenario: PRISM6_EDGE 自由度数验证
- **WHEN** 查询 PRISM6_EDGE 单元的自由度数量
- **THEN** 返回 9（三棱柱有 9 条棱边）

#### Scenario: PYRAMID5_EDGE 自由度数验证
- **WHEN** 查询 PYRAMID5_EDGE 单元的自由度数量
- **THEN** 返回 8（金字塔有 8 条棱边）

### Requirement: 单元工厂类（ShapeFunctionFactory）

系统 SHALL 提供工厂类 `ShapeFunctionFactory`，根据字符串标识符自动创建对应的单元对象：

| 字符串标识符 | 创建的单元类型 |
|------------|--------------|
| `"LINE2"` | LagrangeElement<1>（2节点线单元） |
| `"LINE3"` | LagrangeElement<1>（3节点二次线单元） |
| `"TRI3"` | LagrangeElement<2>（3节点三角形） |
| `"TRI6"` | LagrangeElement<2>（6节点二次三角形） |
| `"QUAD4"` | LagrangeElement<2>（4节点四边形） |
| `"QUAD8"` | LagrangeElement<2>（8节点Serendipity四边形） |
| `"QUAD9"` | LagrangeElement<2>（9节点双二次四边形） |
| `"TET4"` | LagrangeElement<3>（4节点四面体） |
| `"TET10"` | LagrangeElement<3>（10节点二次四面体） |
| `"HEX8"` | LagrangeElement<3>（8节点六面体） |
| `"HEX20"` | LagrangeElement<3>（20节点Serendipity六面体） |
| `"HEX27"` | LagrangeElement<3>（27节点三二次六面体） |
| **`"PRISM6"`** | **LagrangeElement<3>（6节点三棱柱）** |
| **`"PRISM15"`** | **LagrangeElement<3>（15节点二次三棱柱）** |
| **`"PYRAMID5"`** | **LagrangeElement<3>（5节点金字塔）** |
| **`"PYRAMID13"`** | **LagrangeElement<3>（13节点二次金字塔）** |
| `"TRI3_EDGE"` | NedelecElement<2>（三角形棱单元） |
| `"QUAD4_EDGE"` | NedelecElement<2>（四边形棱单元） |
| `"TET4_EDGE"` | NedelecElement<3>（四面体棱单元） |
| `"HEX8_EDGE"` | NedelecElement<3>（六面体棱单元） |
| **`"PRISM6_EDGE"`** | **NedelecElement<3>（三棱柱棱单元）** |
| **`"PYRAMID5_EDGE"`** | **NedelecElement<3>（金字塔棱单元）** |

工厂方法签名：
```cpp
static std::unique_ptr<ShapeFunctionBase> create(const std::string& element_type);
```

返回 `std::unique_ptr<ShapeFunctionBase>`，调用方通过基类指针操作。

#### Scenario: 工厂创建 TRI3 单元
- **WHEN** 调用 `ShapeFunctionFactory::create("TRI3")`
- **THEN** 返回指向 LagrangeElement<2>（TRI3 类型）的智能指针，非空且 getNodeCount() == 3

#### Scenario: 工厂创建 PRISM6 单元
- **WHEN** 调用 `ShapeFunctionFactory::create("PRISM6")`
- **THEN** 返回指向 LagrangeElement<3>（PRISM6 类型）的智能指针，非空且 getNodeCount() == 6

#### Scenario: 工厂创建未知类型
- **WHEN** 调用 `ShapeFunctionFactory::create("UNKNOWN_TYPE")`
- **THEN** 输出 ERROR 日志并返回 nullptr

### Requirement: 测试用例

系统 SHALL 提供可运行的测试程序 `test_shape_function.cpp`，使用 Google Test 框架，对以下内容进行单积分点验证：

1. **Lagrange 单元测试：**
   - TRI3 形心处形函数值与梯度验证
   - QUAD4 中心处形函数值验证
   - TET4 形心处形函数值验证
   - HEX8 中心处形函数值验证
   - **PRISM6 几何中心处形函数值验证**
   - **PYRAMID5 近底面中心处形函数值验证**
   - 单位分解性检验（ΣNᵢ = 1）
   - 雅可比行列式正定性检验（detJ > 0）

2. **Nedelec 棱单元测试：**
   - TRI3_EDGE 基函数切向分量验证
   - TET4_EDGE 自由度数量验证
   - **PRISM6_EDGE 自由度数量验证（应为 9）**
   - **PYRAMID5_EDGE 自由度数量验证（应为 8）**
   - 旋度计算正确性验证

3. **工厂类测试：**
   - 所有 **23 种**单元类型的创建成功验证
   - 未知类型返回 nullptr 验证

## MODIFIED Requirements
无（全新模块，不修改现有代码）。

## REMOVED Requirements
无。
