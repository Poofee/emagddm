# Tree Gauge（树-余树规范）模块 Spec

## Why
3D静磁场求解器使用矢量磁矢位A求解时存在**非唯一性问题**：对于任意标量函数φ，A' = A + ∇φ 都能得到相同的B场。这种非唯一性导致：
- 全局刚度矩阵**奇异**（存在零特征值）
- 直接求解器报错或结果错误
- 迭代求解器收敛极慢甚至不收敛

Tree Gauge是工业级软件（ANSYS Maxwell、JMAG等）的**标准解决方案**，通过图论方法消除冗余自由度，使矩阵可正常求解。当前项目已有完整的DOF管理器和网格数据结构，但缺少Tree Gauge实现，3D静磁场求解器无法正常运行。

## What Changes
- 新增核心类 `tree_gauge.hpp/.cpp`（TreeGauge类）：构建生成树并应用树-余树规范
- 新增图论工具类 `spanning_tree.hpp/.cpp`（SpanningTreeBuilder类）：BFS/DFS生成树构建算法
- 新增多连通域检测器 `multi_connected_domain.hpp/.cpp`（MultiConnectedDomainDetector类）：自动识别孔洞和环量自由度
- 扩展 `em_dof_manager.hpp`：新增Tree Gauge集成接口
- 更新CMakeLists.txt：添加新模块到构建系统
- 新增单元测试验证单连通域和多连通域场景的正确性

## Impact
- Affected specs: DOF管理模块（add-em-dof-manager）、网格模块（add-em-mesh-module）、线性求解器模块（add-em-linear-solver）
- Affected code:
  - `include/fe_em/` — 新增4个头文件
  - `src/fe_em/` — 新增4个源文件
  - `test/` — 新增4个测试文件
  - `CMakeLists.txt` — 更新构建配置

## ADDED Requirements

### Requirement: 图结构提取与邻接表构建
系统 SHALL 提供`MeshGraphExtractor`工具类或静态方法，从EMMeshData中提取图结构用于Tree Gauge计算。

#### 输入：
- `const EMMeshData& mesh_data`：网格拓扑数据
- `const std::vector<std::vector<int>>& elem_local_to_global_edge`：全局棱边ID映射表

#### 输出：
- `std::vector<std::vector<int>> adjacency_list`：节点邻接表（adjacency_list[node_id] = 相邻节点ID列表）
- `std::unordered_map<std::pair<int,int>, int, PairHash> edge_to_global_id`：无序节点对 → 全局棱边ID的映射
- `int num_nodes`：图中总节点数
- `int num_edges`：图中总棱边数

#### 实现逻辑：
1. 遍历所有VECTOR_EDGE_ONLY和MIXED_AV类型单元
2. 对每个单元的每条局部棱边，获取其两个节点的全局ID
3. 构建无向图的邻接表（避免重复添加）
4. 记录每条棱边对应的全局棱边ID（用于后续树边/余树边分类）

#### Scenario: PRISM6单元图提取
- **WHEN** 提取1个PRISM6单元（6节点9棱边）的图结构
- **THEN** 邻接表包含6个节点，每个节点的邻居列表正确反映单元拓扑，edge_to_global_id包含9条映射

### Requirement: BFS/DFS生成树构建算法
系统 SHALL 提供`SpanningTreeBuilder`类，使用广度优先搜索（BFS）或深度优先搜索（DFS）构建生成树。

#### 核心接口：
```cpp
class SpanningTreeBuilder {
public:
    struct Result {
        std::set<int> tree_edges;           ///< 树边的全局棱边ID集合
        std::set<int> cotree_edges;         ///< 余树边的全局棱边ID集合
        std::unordered_map<int, int> node_parent;  ///< 节点→父节点映射（根节点父节点为-1）
        int root_node;                      ///< 根节点ID
    };

    // 使用BFS构建生成树（推荐，缓存友好）
    static Result buildBFSTree(const std::vector<std::vector<int>>& adjacency_list,
                               const std::unordered_map<std::pair<int,int>, int, PairHash>& edge_to_global_id,
                               int root_node = 0);

    // 使用DFS构建生成树（备选，递归深度可控）
    static Result buildDFSTree(const std::vector<std::vector<int>>& adjacency_list,
                               const std::unordered_map<std::pair<int,int>, int, PairHash>& edge_to_global_id,
                               int root_node = 0);
};
```

#### 算法步骤（以BFS为例）：
1. 初始化：选择root_node作为根节点，标记为已访问
2. 入队：将root_node加入队列
3. 循环：当队列非空时
   - 出队当前节点u
   - 遍历u的所有相邻节点v
   - 若v未被访问：
     - 标记v为已访问
     - 将边(u,v)加入tree_edges（通过edge_to_global_id查找全局棱边ID）
     - 设置node_parent[v] = u
     - 将v入队
4. 结束后，所有不在tree_edges中的边即为cotree_edges

#### 性能要求：
- 时间复杂度：O(V + E)，其中V为节点数，E为棱边数
- 空间复杂度：O(V + E)
- 必须支持百万级节点的网格（内存效率优化）

#### Scenario: 简单立方体网格生成树
- **WHEN** 8节点12棱边的HEX8单元网格构建BFS生成树
- **THEN** tree_edges包含7条边（V-1），cotree_edges包含5条边（E-V+1），所有节点都被访问到

### Requirement: Dirichlet边界优先策略
系统 SHALL 支持**智能根节点选择**和**边界边优先策略**，确保Dirichlet边界条件与Tree Gauge规范条件一致。

#### 实现逻辑：
1. 从约束棱边集合（constraint_data_.constrained_edges）中选取一条边的一个节点作为根节点
2. 在BFS/DFS遍历时，**优先遍历Dirichlet边界边**（将它们放入树中）
3. 这样可以保证Dirichlet边界边属于tree_edges，其A值被强制置0，与Dirichlet边界条件一致

#### 接口扩展：
```cpp
// 带边界优先的BFS生成树
static Result buildBFSTreeWithBoundaryPriority(
    const std::vector<std::vector<int>>& adjacency_list,
    const std::unordered_map<std::pair<int,int>, int, PairHash>& edge_to_global_id,
    const std::set<int>& dirichlet_edges,  ///< Dirichlet边界棱边集合
    int root_node = 0);
```

#### Scenario: 带Dirichlet边界的模型
- **WHEN** 模型有2条Dirichlet边界棱边
- **THEN** 这2条边被优先包含在tree_edges中，确保A=0与边界条件一致

### Requirement: 多连通域检测与环量自由度保留
系统 SHALL 提供`MultiConnectedDomainDetector`类，自动检测多连通域（孔洞）并保留对应的环量自由度。

#### 数学背景：
对于有g个孔洞的多连通域：
- 生成树的边数 = V - 1 + g（而非简单的V-1）
- 每个孔洞对应一个**环量自由度**（不可置为0）
- 这些环量自由度对应于围绕孔洞的独立回路

#### 核心接口：
```cpp
class MultiConnectedDomainDetector {
public:
    struct Result {
        int num_holes;                     ///< 孔洞数量（亏格genus）
        std::vector<std::vector<int>> hole_loop_edges;  ///< 每个孔洞对应的回路棱边集合
        std::set<int> circulation_dofs;    ///< 需要保留的环量自由度（余树边子集）
    };

    // 检测多连通域并识别环量自由度
    static Result detect(const EMMeshData& mesh_data,
                         const SpanningTreeBuilder::Result& tree_result);
};
```

#### 检测算法（简化版）：
1. 基于网格拓扑构建对偶图
2. 使用并查集（Union-Find）检测连通分量
3. 识别内部边界（孔洞边界）
4. 为每个孔洞选择一条独立的余树边作为环量自由度

#### Scenario: 环形线圈（toroid）模型
- **WHEN** 环形线圈模型有1个中心孔洞
- **THEN** num_holes=1，circulation_dofs包含1条余树边（代表穿过孔洞的磁通）

### Requirement: Tree Gauge核心应用类
系统 SHALL 提供`TreeGauge`类，整合上述组件，执行完整的Tree Gauge流程。

#### 核心接口：
```cpp
class TreeGauge {
public:
    /**
     * @brief 构造函数
     * @param mesh_data 网格数据
     * @param elem_local_to_global_edge 全局棱边ID映射表
     * @param constrained_edges 已约束的全局棱边ID集合（来自DOF管理器）
     */
    explicit TreeGauge(const EMMeshData& mesh_data,
                       const std::vector<std::vector<int>>& elem_local_to_global_edge,
                       const std::set<int>& constrained_edges = {});

    /**
     * @brief 执行完整的Tree Gauge流程
     * 1. 提取图结构
     * 2. 构建生成树（带边界优先）
     * 3. 检测多连通域
     * 4. 分类树边/余树边/环量自由度
     * 5. 生成约化映射
     */
    void build();

    /**
     * @brief 获取约化后的自由DOF数量
     * @return int 余树边数量（减去被约束的余树边）
     */
    int getReducedNumDofs() const;

    /**
     * @brief 获取树边集合（这些边的A将被置0）
     * @return const std::set<int>& 树边的全局棱边ID集合
     */
    const std::set<int>& getTreeEdges() const;

    /**
     * @brief 获取余树边集合（需要求解的自由度）
     * @return const std::set<int>& 余树边的全局棱边ID集合（不包含环量自由度）
     */
    const std::set<int>& getCotreeEdges() const;

    /**
     * @brief 获取环量自由度集合（多连通域的特殊自由度）
     * @return const std::set<int>& 环量自由度的全局棱边ID集合
     */
    const std::set<int>& getCirculationDOFs() const;

    /**
     * @brief 获取余树边到约化DOF编号的映射
     * @return const std::unordered_map<int, int>& 全局棱边ID → 约化DOF编号
     */
    const std::unordered_map<int, int>& getCotreeEdgeToReducedDOF() const;

    /**
     * @brief 将约化解向量映射回完整的棱边解向量
     * @param reduced_solution 约化系统的解向量（仅含余树边自由度）
     * @param full_solution 输出：完整的棱边解向量（树边位置为0）
     */
    void mapReducedSolutionToFull(const std::vector<double>& reduced_solution,
                                  std::vector<double>& full_solution) const;

    /**
     * @brief 判断某条棱边是否为树边
     * @param global_edge_id 全局棱边ID
     * @return true 如果是树边
     */
    bool isTreeEdge(int global_edge_id) const;

    /**
     * @brief 判断某条棱边是否为余树边
     * @param global_edge_id 全局棱边ID
     * @return true 如果是余树边（包括环量自由度）
     */
    bool isCotreeEdge(int global_edge_id) const;

private:
    const EMMeshData& mesh_data_;
    const std::vector<std::vector<int>> elem_local_to_global_edge_;
    const std::set<int> constrained_edges_;

    // 图结构
    std::vector<std::vector<int>> adjacency_list_;
    std::unordered_map<std::pair<int,int>, int, PairHash> edge_to_global_id_;
    int num_nodes_ = 0;
    int num_edges_ = 0;

    // 生成树结果
    SpanningTreeBuilder::Result tree_result_;

    // 多连通域检测结果
    MultiConnectedDomainDetector::Result domain_result_;

    // 最终分类
    std::set<int> final_tree_edges_;       ///< 最终树边（包含被强制置0的边）
    std::set<int> final_cotree_edges_;     ///< 最终余树边（待求解的自由度）
    std::unordered_map<int, int> cotree_edge_to_reduced_dof_;  ///< 余树边→约化DOF编号
    int reduced_num_dofs_ = 0;             ///< 约化后的DOF数量

    // 内部方法
    void extractGraphStructure();
    void buildSpanningTree();
    void detectMultiConnectedDomain();
    void classifyEdges();
    void buildReducedMapping();
};
```

#### build()流程详细说明：
1. **extractGraphStructure()**：从网格数据提取邻接表和棱边映射
2. **buildSpanningTree()**：调用SpanningTreeBuilder::buildBFSTreeWithBoundaryPriority()
3. **detectMultiConnectedDomain()**：调用MultiConnectedDomainDetector::detect()
4. **classifyEdges()**：
   - tree_edges = 生成树边 ∪ 被约束的余树边（这些边A=0）
   - cotree_edges = 未被约束的余树边 - 环量自由度
   - circulation_dofs = 多连通域检测到的环量自由度
5. **buildReducedMapping()**：为每个最终余树边分配连续的约化DOF编号

#### Scenario: 单连通域完整流程
- **WHEN** 单连通域模型（无孔洞），10条棱边，2条Dirichlet边界棱边
- **THEN** tree_edges包含9条边（V-1=9，假设10节点），cotree_edges包含1条边，reduced_num_dofs_=1

### Requirement: 与DOF管理器的集成
系统 SHALL 扩展`EMDOFManager`类，提供Tree Gauge集成接口。

#### 新增公共方法：
```cpp
class EMDOFManager {
public:
    // ... 现有方法 ...

    /**
     * @brief 应用Tree Gauge规范（可选调用）
     * @details 在build()之后调用，进一步压缩矢量棱边DOF：
     *          - 树边DOF从系统中剔除（A=0）
     *          - 仅保留余树边DOF参与求解
     *          - 显著减少系统规模（通常减少30-50%）
     *
     * @note 调用后，getNumFreeDOFs()返回的是约化后的DOF数量
     * @note getElemLocalToGlobal()中的树边DOF标记为-1
     * @warning 此方法仅对VECTOR_EDGE_ONLY和MIXED_AV类型的矢量DOF有效
     */
    void applyTreeGauge();

    /**
     * @brief 获取Tree Gauge对象（用于查询树边/余树边信息）
     * @return const TreeGauge* Tree Gauge对象的常量指针（未调用applyTreeGauge时返回nullptr）
     */
    const TreeGauge* getTreeGauge() const;

private:
    std::unique_ptr<TreeGauge> tree_gauge_;  ///< Tree Gauge对象（延迟初始化）
    int original_num_free_dofs_ = 0;         ///< applyTreeGauge前的原始自由DOF数
    std::unordered_map<int, int> free_to_reduced_;  ///< 原始自由编号→约化编号映射
};
```

#### applyTreeGauge()实现逻辑：
1. 创建TreeGauge对象并调用build()
2. 保存original_num_free_dofs_
3. 遍历elem_local_to_global_中的所有矢量DOF：
   - 若该棱边是树边或被约束：indices[i] = -1
   - 若该棱边是余树边：indices[i] = new_reduced_number
4. 更新num_free_dofs_为reduced_num_dofs_
5. 构建free_to_reduced_映射供解恢复使用

#### Scenario: DOF管理器集成示例
- **WHEN** 原始自由DOF=100（70标量+30矢量），applyTreeGauge后树边=20，余树边=10
- **THEN** 新的自由DOF=80（70标量+10余树边矢量），20个树边DOF被置为-1

### Requirement: 高性能优化特性
系统 SHALL 实现以下性能优化以确保大规模模型的实时响应：

#### 1. 内存布局优化
- 邻接表使用`std::vector<int>`而非`std::set<int>`（缓存友好）
- 棱边映射使用开放寻址哈希表（比std::unordered_map更快）
- 预分配容器容量（reserve）避免动态扩容

#### 2. 并行化准备（Phase 1）
- 图结构提取阶段支持OpenMP并行（遍历单元时可并行）
- 生成树构建本身为串行算法（BFS/DFS难以并行），但可并行处理多个 disconnected components
- 多连通域检测可并行处理多个候选孔洞

#### 3. 智能根节点选择
- 优先选择**度数最大**的节点作为根（减少BFS树的高度）
- 或选择**Dirichlet边界上**的节点（保证边界兼容性）
- 可配置的策略枚举：
```cpp
enum class RootNodeStrategy {
    FIRST_NODE,           ///< 第一个节点（默认，简单快速）
    MAX_DEGREE,           ///< 度数最大的节点（推荐，树更平衡）
    DIRICHLET_BOUNDARY,   ///< Dirichlet边界上的节点（保证边界兼容性）
    USER_SPECIFIED        ///< 用户指定的节点ID
};
```

#### 性能基准目标：
| 模型规模 | 节点数 | 棱边数 | Tree Gauge构建时间 | 内存占用 |
|---------|--------|--------|-------------------|----------|
| 小型     | 1K     | 5K     | < 1 ms            | < 1 MB   |
| 中型     | 100K   | 500K   | < 50 ms           | < 50 MB  |
| 大型     | 1M     | 5M     | < 500 ms          | < 500 MB |
| 超大型   | 10M    | 50M    | < 5 s             | < 5 GB   |

### Requirement: 日志与诊断输出
系统 SHALL 提供详细的日志输出，便于调试和性能分析：

#### 日志级别：
- **INFO**：Tree Gauge构建完成摘要（树边数、余树边数、环量自由度数、压缩率）
- **DEBUG**：详细的算法步骤信息（根节点选择、BFS队列状态、每步分类结果）
- **WARN**：异常情况警告（孤立节点、未连接的子图、退化情况）

#### 统计信息输出示例：
```
[INFO] ========== Tree Gauge构建完成 ==========
[INFO] 图结构: 节点数=10000, 棱边数=55000
[INFO] 生成树: 根节点=42(策略=MAX_DEGREE), 树边数=9999, 余树边数=45001
[INFO] 多连通域: 孔洞数=0, 环量自由度数=0
[INFO] 最终分类: 树边=9999(含约束边12), 余树边=44989, 约化DOF数=44989
[INFO] 压缩率: 81.8% (原55000→44989)
[INFO] ========================================
```

## MODIFIED Requirements
无。本模块为全新功能，不修改现有代码逻辑（仅扩展EMDOFManager的可选接口）。

## REMOVED Requirements
无。
