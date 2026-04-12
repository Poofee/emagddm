# Tree Gauge模块开发任务列表

## 任务概览
基于`docs/solver/tree.md`技术文档，实现高性能的Tree Gauge（树-余树规范）模块，用于3D静磁场求解器的规范条件施加。

# Tasks

- [x] Task 1: 实现图结构提取器（MeshGraphExtractor）
  - [x] 1.1 创建 `include/fe_em/mesh_graph_extractor.hpp` 头文件
    - 定义MeshGraphExtractor类/命名空间
    - 实现从EMMeshData提取邻接表的接口
    - 实现棱边到全局ID的映射结构（PairHash支持）
  - [x] 1.2 创建 `src/fe_em/mesh_graph_extractor.cpp` 实现文件
    - 遍历VECTOR_EDGE_ONLY和MIXED_AV单元
    - 构建无向邻接表（去重）
    - 构建edge_to_global_id映射
  - [x] 1.3 创建 `test/test_mesh_graph_extractor.cpp` 测试文件
    - 测试PRISM6单元的图提取正确性
    - 测试HEX8单元的多单元邻接表合并
    - 测试空网格和异常输入的处理

- [x] Task 2: 实现BFS/DFS生成树构建算法（SpanningTreeBuilder）
  - [x] 2.1 创建 `include/fe_em/spanning_tree.hpp` 头文件
    - 定义SpanningTreeBuilder类
    - 定义Result结构体（tree_edges, cotree_edges, node_parent, root_node）
    - 声明buildBFSTree()和buildDFSTree()静态方法
    - 声明buildBFSTreeWithBoundaryPriority()方法
    - 定义RootNodeStrategy枚举
  - [x] 2.2 创建 `src/fe_em/spanning_tree.cpp` 实现文件
    - 实现BFS算法（使用std::queue）
    - 实现DFS算法（使用递归或显式栈）
    - 实现边界优先策略（优先遍历Dirichlet边）
    - 实现智能根节点选择（MAX_DEGREE策略）
  - [x] 2.3 创建 `test/test_spanning_tree.cpp` 测试文件
    - 测试简单图的BFS树正确性（V=4, E=5）
    - 测试HEX8单元网格的生成树（8节点12棱边→7树边5余树边）
    - 测试边界优先策略的正确性
    - 测试孤立节点的处理

- [x] Task 3: 实现多连通域检测器（MultiConnectedDomainDetector）
  - [x] 3.1 创建 `include/fe_em/multi_connected_domain.hpp` 头文件
    - 定义MultiConnectedDomainDetector类
    - 定义Result结构体（num_holes, hole_loop_edges, circulation_dofs）
    - 声明detect()静态方法
  - [x] 3.2 创建 `src/fe_em/multi_connected_domain.cpp` 实现文件
    - 实现基于并查集的连通分量检测
    - 实现孔洞边界识别算法
    - 实现环量自由度选择逻辑
  - [x] 3.3 创建 `test/test_multi_connected_domain.cpp` 测试文件
    - 测试单连通域（num_holes=0）
    - 测试环形线圈模型（num_holes=1, circulation_dofs=1）
    - 测试多孔洞场景（num_holes>1）

- [x] Task 4: 实现Tree Gauge核心应用类（TreeGauge）
  - [x] 4.1 创建 `include/fe_em/tree_gauge.hpp` 头文件
    - 定义TreeGauge类（完整接口见spec.md Requirement 5）
    - 包含所有公共方法和私有成员
    - 包含与MeshGraphExtractor、SpanningTreeBuilder、MultiConnectedDomainDetector的组合关系
  - [x] 4.2 创建 `src/fe_em/tree_gauge.cpp` 实现文件
    - 实现extractGraphStructure()：调用MeshGraphExtractor
    - 实现buildSpanningTree()：调用SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority()
    - 实现detectMultiConnectedDomain()：调用MultiConnectedDomainDetector::detect()
    - 实现classifyEdges()：合并树边、约束边、余树边、环量自由度
    - 实现buildReducedMapping()：为余树边分配约化DOF编号
    - 实现mapReducedSolutionToFull()：解向量恢复
    - 完整实现build()主流程
    - 实现所有getter方法
  - [x] 4.3 创建 `test/test_tree_gauge.cpp` 测试文件
    - 测试单连通域完整流程（10节点模型）
    - 测试带Dirichlet边界模型的正确性
    - 测试多连通域模型的环量自由度保留
    - 测试mapReducedSolutionToFull()解恢复正确性
    - 测试性能基准（中型模型<50ms）

- [x] Task 5: 集成Tree Gauge到DOF管理器
  - [x] 5.1 扩展 `include/fe_em/em_dof_manager.hpp`
    - 新增applyTreeGauge()公共方法声明
    - 新增getTreeGauge()公共方法声明
    - 新增私有成员：tree_gauge_（unique_ptr<TreeGauge>）、original_num_free_dofs_、free_to_reduced_
  - [x] 5.2 扩展 `src/fe_em/em_dof_manager.cpp`
    - 实现applyTreeGauge()：创建TreeGauge对象、重新编号矢量DOF、更新映射表
    - 实现getTreeGauge()：返回tree_gauge_指针
  - [x] 5.3 创建 `test/test_tree_gauge_integration.cpp` 测试文件
    - 测试applyTreeGauge()后DOF数量减少
    - 测试Local2Global映射表中树边标记为-1
    - 测试getTreeGauge()返回有效指针
    - 测试不调用applyTreeGauge时getTreeGauge()返回nullptr

- [x] Task 6: 更新构建系统与最终验证
  - [x] 6.1 更新CMakeLists.txt
    - 添加新源文件到fe_em库
    - 添加新测试目标
  - [x] 6.2 运行完整测试套件
    - 执行所有新增测试用例
    - 确保通过率100%（75/75测试全部通过）
  - [x] 6.3 性能基准测试
    - 使用中型规模模型（400节点网格）测试性能
    - 验证构建时间<3ms（远优于50ms目标）
  - [x] 6.4 日志输出验证
    - 确认INFO级别日志包含完整的统计摘要
    - 确认DEBUG级别日志包含详细的算法步骤

# Task Dependencies
- [Task 1] 无依赖（可最先开始）
- [Task 2] 依赖 [Task 1]（需要图结构数据结构）
- [Task 3] 可与 [Task 2] 并行（独立算法）
- [Task 4] 依赖 [Task 1], [Task 2], [Task 3]（整合所有组件）
- [Task 5] 依赖 [Task 4]（需要TreeGauge完整实现）
- [Task 6] 依赖 [Task 4], [Task 5]（需要所有代码完成）

## 并行执行建议
**Phase 1（并行）**：
- Sub-Agent A: Task 1（图结构提取器）
- Sub-Agent B: Task 2 + Task 3（生成树 + 多连通域检测，Task 3可在Task 2完成后开始或稍后启动）

**Phase 2（串行）**：
- Task 4（TreeGauge核心类，依赖Phase 1的所有组件）

**Phase 3（串行）**：
- Task 5（DOF管理器集成）

**Phase 4（串行）**：
- Task 6（构建系统与验证）
