# Tasks

- [x] Task 1: 创建形函数抽象基类 `ShapeFunctionBase`
  - [x] Task 1.1: 定义单元类型枚举 `ElementType`（覆盖全部 **22 种**单元类型：16 Lagrange + 6 Nedelec）
  - [x] Task 1.2: 定义纯虚接口方法（getNodeType/getNodeCount/getDim/evalN/evalGradN/calcJacobian/calcPhysicalGradN）
  - [x] Task 1.3: 为 Nedelec 单元定义扩展虚接口（evalEdgeFunction/evalCurlEdge/calcPhysicalEdgeFunction/calcPhysicalCurl）
  - [x] Task 1.4: 定义数据结构（局部坐标点、节点坐标矩阵、雅可比结果结构体）

- [x] Task 2: 实现 Lagrange 标量节点单元（**16 种**类型）
  - [x] Task 2.1: 实现 1D 线单元（LINE2、LINE3）的形函数与偏导数
  - [x] Task 2.2: 实现 2D 三角形单元（TRI3、TRI6）的形函数与偏导数
  - [x] Task 2.3: 实现 2D 四边形单元（QUAD4、QUAD8、QUAD9）的形函数与偏导数
  - [x] Task 2.4: 实现 3D 四面体单元（TET4、TET10）的形函数与偏导数
  - [x] Task 2.5: 实现 3D 六面体单元（HEX8、HEX20、HEX27）的形函数与偏导数
  - [x] Task 2.6: 实现 3D 三棱柱单元（PRISM6、PRISM15）的形函数与偏导数
  - [x] Task 2.7: 实现 3D 金字塔单元（PYRAMID5、PYRAMID13）的形函数与偏导数

- [x] Task 3: 实现 Nedelec 1阶棱边矢量单元（**6 种**类型）
  - [x] Task 3.1: 实现 TRI3_EDGE 棱边基函数与旋度公式
  - [x] Task 3.2: 实现 QUAD4_EDGE 棱边基函数与旋度公式
  - [x] Task 3.3: 实现 TET4_EDGE 棱边基函数与旋度公式
  - [x] Task 3.4: 实现 HEX8_EDGE 棱边基函数与旋度公式
  - [x] Task 3.5: 实现 PRISM6_EDGE 棱边基函数（9 条棱边）与旋度公式
  - [x] Task 3.6: 实现 PYRAMID5_EDGE 棱边基函数（8 条棱边）与旋度公式
  - [x] Task 3.7: 实现物理坐标系 Piola 变换映射（矢量基函数 + 旋度）

- [x] Task 4: 实现单元工厂类 `ShapeFunctionFactory`
  - [x] Task 4.1: 实现字符串到单元类型的注册表映射（22 种）
  - [x] Task 4.2: 实现 create() 工厂方法，返回 unique_ptr<ShapeFunctionBase>
  - [x] Task 4.3: 实现错误处理（未知类型日志 + 返回 nullptr）

- [x] Task 5: 编写测试用例并配置编译
  - [x] Task 5.1: 编写 Lagrange 单元测试（TRI3/QUAD4/TET4/HEX8/**PRISM6**/**PYRAMID5** 形心验证 + 单位分解性 + 雅可比正定性）
  - [x] Task 5.2: 编写 Nedelec 棱单元测试（TRI3_EDGE/TET4_EDGE/**PRISM6_EDGE**/**PYRAMID5_EDGE** 基函数切向分量 + 自由度数量 + 旋度正确性）
  - [x] Task 5.3: 编写工厂类测试（**22 种**类型创建成功 + 未知类型返回 nullptr）
  - [x] Task 5.4: 更新 src/CMakeLists.txt 和 test/CMakeLists.txt 编译配置

# Task Dependencies
- [Task 2] depends on [Task 1] — Lagrange 单元依赖基类接口定义
- [Task 3] depends on [Task 1] — Nedelec 单元依赖基类接口定义
- [Task 4] depends on [Task 2, Task 3] — 工厂需要所有具体单元实现完成
- [Task 5] depends on [Task 2, Task 3, Task 4] — 测试需要全部模块就绪
