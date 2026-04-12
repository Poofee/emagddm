# Tree Gauge模块验证清单

## 验证策略
基于spec.md中的需求定义，系统化验证每个功能点的正确性、性能和鲁棒性。

- [ ] 图结构提取器（MeshGraphExtractor）正确性验证
  - [ ] PRISM6单元（6节点9棱边）的邻接表正确反映单元拓扑
  - [ ] 多个HEX8单元共享节点时，邻接表正确合并
  - [ ] edge_to_global_id映射覆盖所有棱边，无重复无遗漏
  - [ ] 空网格输入时返回空结构而不崩溃
  - [ ] 仅含标量单元的网格返回空图结构

- [ ] BFS生成树算法正确性验证
  - [ ] 简单图（4节点5边）：树边数=3（V-1），余树边数=2（E-V+1）
  - [ ] HEX8单元（8节点12边）：树边数=7，余树边数=5
  - [ ] 所有节点都被访问到（无孤立节点遗漏）
  - [ ] 树中无回路（通过node_parent可重建唯一路径）
  - [ ] DFS算法与BFS算法产生相同数量的树边/余树边（数量一致，具体边可能不同）

- [ ] Dirichlet边界优先策略验证
  - [ ] 提供的Dirichlet边界边被包含在tree_edges中
  - [ ] 根节点选在Dirichlet边界上（使用DIRICHLET_BOUNDARY策略）
  - [ ] 无Dirichlet边界时回退到默认策略（FIRST_NODE或MAX_DEGREE）
  - [ ] 混合边界条件（部分Dirichlet + 其他类型）处理正确

- [ ] 多连通域检测器验证
  - [ ] 单连通域模型：num_holes=0, circulation_dofs为空集
  - [ ] 环形线圈（toroid）模型：num_holes=1, circulation_dofs包含1条边
  - [ ] 双孔洞模型：num_holes=2, circulation_dofs包含2条独立边
  - [ ] hole_loop_edges中的回路确实围绕对应的孔洞
  - [ ] circulation_dofs是cotree_edges的子集（环量自由度来自余树边）

- [ ] TreeGauge核心类完整流程验证
  - [ ] build()方法按正确顺序执行5个步骤（提取图→构建树→检测多连通域→分类边→构建映射）
  - [ ] getReducedNumDofs() = cotree_edges.size() - constrained_cotree_edges.size()
  - [ ] getTreeEdges() ∪ getCotreeEdges() ∪ getCirculationDOFs() = 所有棱边的完备划分
  - [ ] getTreeEdges() ∩ getCotreeEdges() = 空集（互斥）
  - [ ] isTreeEdge()和isCotreeEdge()对所有棱边返回正确的布尔值
  - [ ] cotree_edge_to_reduced_dof_映射是双射（一一对应，连续编号从0开始）

- [ ] 解向量恢复验证（mapReducedSolutionToFull）
  - [ ] full_solution长度等于总棱边数
  - [ ] tree_edges对应的位置值为0.0
  - [ ] cotree_edges对应的位置值与reduced_solution一致
  - [ ] circulation_dofs对应的位置值与reduced_solution中对应位置一致
  - [ ] constrained_edges对应的位置保持约束值（不被覆盖）

- [ ] DOF管理器集成验证
  - [ ] applyTreeGauge()调用后getNumFreeDOFs()减少（减少量≈树边数）
  - [ ] applyTreeGauge()后Local2Global映射表中树边DOF标记为-1
  - [ ] applyTreeGauge()后Local2Global映射表中余树边DOF重新编号（连续从0开始）
  - [ ] getTreeGauge()在applyTreeGauge()后返回非空指针
  - [ ] getTreeGauge()在未调用applyTreeGauge()时返回nullptr
  - [ ] 标量DOF不受影响（仅矢量棱边DOF被压缩）
  - [ ] 原始约束DOF仍然标记为-1（不被Tree Gauge覆盖）

- [ ] 性能基准验证
  - [ ] 小型模型（1K节点5K棱边）：构建时间<1ms
  - [ ] 中型模型（100K节点500K棱边）：构建时间<50ms
  - [ ] 内存占用随规模线性增长（无二次或指数爆炸）
  - [ ] 邻接表使用vector<int>而非set<int>（缓存友好验证）

- [ ] 日志输出验证
  - [ ] INFO级别输出完整的统计摘要（节点数、棱边数、树边数、余树边数、压缩率）
  - [ ] DEBUG级别输出根节点选择信息
  - [ ] WARN级别对孤立节点和退化情况发出警告
  - [ ] 使用FEEM_INFO/FEEM_DEBUG/FEEM_WARN宏（符合项目规范）

- [ ] 边界情况与鲁棒性验证
  - [ ] 所有棱边都被约束的极端情况（reduced_num_dofs_=0）
  - [ ] 无约束边的理想情况（最大压缩率）
  - [ ] 单个单元的退化网格
  - [ ] 断开连接的多组件网格（多个 disconnected components）
  - [ ] 大度数节点（>100个邻居）的性能不退化

- [ ] 编码规范遵循验证
  - [ ] 所有头文件包含完整的Doxygen注释（@file @brief @class @param @return @note）
  - [ ] 使用项目规定的日志宏（FEEM_INFO/FEEM_DEBUG/FEEM_WARN/FEEM_ERROR）
  - [ ] 头文件引用无路径前缀（仅文件名，路径在CMakeLists.txt中配置）
  - [ ] 无全局变量（所有状态封装在类成员中）
  - [ ] 无内存泄漏（使用智能指针管理动态资源）
  - [ ] 异常处理完善（std::invalid_argument用于无效输入）
  - [ ] const正确性（只读方法标记为const）
