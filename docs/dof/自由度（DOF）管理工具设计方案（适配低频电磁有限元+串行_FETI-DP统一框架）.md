# 自由度（DOF）管理工具设计方案（适配低频电磁有限元+串行/FETI-DP统一框架）

该工具是稀疏矩阵模块与网格/物理场模块的**核心桥梁**，承接网格拓扑信息、物理场DOF类型要求，输出稀疏矩阵组装所需的**DOF编号映射、批量DOF操作、界面DOF拓扑关系**，同时为FETI-DP区域分解提供**子域局部DOF/全局DOF/界面DOF**的三层索引体系，是实现矩阵自动化组装、分布式耦合的关键底层模块。

设计遵循**「电磁场景深度适配、串行/分布式统一、拓扑解耦、操作高效、工程化易用」**核心原则，聚焦**DOF编号管理、DOF集操作、界面DOF提取**三大核心功能，同时明确与网格、稀疏矩阵、内存管理工具的协同逻辑，以及全阶段开发落地要求，可直接交给AI按此实现。

## 一、工具整体定位与核心设计原则

### 1. 整体定位

- 作为**通用底层支撑模块**，贯穿稀疏矩阵开发全6个阶段，阶段1实现基础DOF编号，阶段3适配电磁标量/矢量元DOF类型，阶段4-5适配分布式/FETI-DP的界面DOF提取，阶段6完成工程化集成；

- 承上：从网格模块获取**节点/边/面拓扑信息**，从物理场模块获取**DOF类型（节点/边/面）、边界条件、场类型（静电/静磁/涡流）**要求；

- 启下：为稀疏矩阵模块提供**DOF索引映射、批量DOF索引、界面DOF拓扑关系**，直接支撑COO/CSR/块稀疏矩阵的组装、分布式子域矩阵构建、FETI-DP界面约束矩阵Bᵢ生成。

### 2. 核心设计原则

|原则|具体要求|
|---|---|
|电磁场景定制|原生支持**标量元（静电）-节点DOF**、**矢量元（静磁/涡流）-边DOF**，预留面DOF扩展；支持实/复DOF标记、边界DOF（Dirichlet/Neumann）标记|
|串行/分布式统一|采用**「全局DOF为基准，局部DOF为映射」**架构，1进程时局部DOF等价于全局DOF，无代码冗余|
|拓扑解耦|不侵入网格模块核心逻辑，通过**拓扑接口**获取网格信息，网格模块替换时仅需适配接口|
|操作高效|大规模DOF（百万/千万级）下，编号映射、集操作、界面提取的时间复杂度控制在O(N)或O(logN)|
|工程化易用|封装高层接口，隐藏底层拓扑计算、编号映射细节，AI开发矩阵模块时直接调用，无需关注网格拓扑|
|模块协同|与内存管理工具协同（DOF数组预分配/缓存对齐）、与稀疏矩阵模块协同（DOF索引直接赋值）|
|鲁棒性|内置DOF编号越界、重复、拓扑冲突检测，调试阶段抛出明确异常，生产阶段快速失败|
## 二、核心功能模块设计（三大核心）

### 模块1：DOF编号管理（DofNumberer）

#### 核心目标

实现**「单元局部DOF → 子域局部DOF → 全局DOF」**三层编号的自动映射与管理，支持DOF重排序（减少矩阵带宽）、DOF类型/属性标记、边界DOF处理，为矩阵组装提供**唯一、无冲突的DOF索引**。

#### 核心设计方案

##### 1. 三层DOF编号体系（适配串行+FETI-DP）

是整个工具的核心架构，从下到上支撑有限元单元组装→子域矩阵构建→全局矩阵/分布式耦合，**所有编号均为非负整数，从0开始连续编号**（适配C++数组索引）。

- **单元局部DOF**：单个单元内的DOF编号（如四面体标量元4个节点→03，四面体矢量元6条边→05），由**单元类型+物理场**决定，与全局/子域无关；

- **子域局部DOF**：某一子域内的全局DOF子集的本地编号（仅分布式场景有效），范围`[0, local_dof_size-1]`，通过**局部→全局映射表**关联全局DOF；

- **全局DOF**：整个计算域的唯一DOF编号，范围`[0, global_dof_size-1]`，串行场景下直接使用，分布式场景下作为子域间通信的统一基准。

##### 2. 核心数据结构（封装为模板类`DofNumberer<DoFType>`，DoFType为NODE/EDGE/FACE）

```C++

// DOF类型枚举（电磁核心适配NODE/EDGE）
enum class DofType { NODE, EDGE, FACE };
// DOF属性枚举（电磁场景必备）
enum class DofAttr { NONE, DIRICHLET, NEUMANN, INTERFACE, PERIODIC };

template<DofType T>
class DofNumberer {
public:
    // 全局/局部DOF规模
    size_t global_dof_size() const;
    size_t local_dof_size() const;
    // 编号映射：单元局部→全局，子域局部→全局，全局→子域局部（分布式）
    size_t elem2global(size_t elem_id, size_t elem_local_dof) const;
    size_t local2global(size_t local_dof) const;
    size_t global2local(size_t global_dof) const;
    // DOF属性操作
    void set_dof_attr(size_t global_dof, DofAttr attr);
    DofAttr get_dof_attr(size_t global_dof) const;
    // DOF重排序
    void reorder(ReorderType type = ReorderType::CMK); // CMK=Cuthill-McKee
private:
    // 核心映射表（由内存管理工具的DynamicArrayPreallocator管理，预分配+缓存对齐）
    std::vector<size_t> elem2global_map_;  // 单元局部→全局：[elem_id*elem_dof_num + elem_local_dof] → global_dof
    std::vector<size_t> local2global_map_; // 子域局部→全局
    std::vector<size_t> global2local_map_; // 全局→子域局部（分布式，全局DOF不在当前子域则为INVALID_DOF）
    // DOF属性与标记（电磁场景）
    std::vector<DofAttr> dof_attr_;        // 全局DOF属性
    std::vector<bool> is_complex_;         // 是否为复DOF（涡流场）
    // 网格拓扑依赖（通过接口解耦，仅存索引，不存网格对象）
    const MeshTopo* mesh_topo_;            // 网格拓扑接口指针
    DofType dof_type_;                     // 当前DOF类型（NODE/EDGE）
    size_t elem_dof_num_;                  // 单个单元的DOF数（由单元类型+DOF类型决定）
};
```

##### 3. 核心功能实现

1. **自动编号生成**：输入网格拓扑+DOF类型+子域拆分信息，自动遍历网格**节点/边**，生成无冲突的全局DOF，再按子域拆分生成局部DOF映射表；<br>

例：静磁场（边DOF）+ 四面体单元 + 2个子域拆分，工具自动遍历所有网格边生成全局边DOF，再将边按子域归属拆分为两个局部DOF集，生成映射表。

1. **单元局部→全局映射**：预计算`elem2global_map_`，将**单元ID+单元局部DOF**拼接为一维索引，直接查表得到全局DOF，**O(1)时间复杂度**（矩阵组装核心需求）；

2. **DOF重排序**：内置**Cuthill-McKee（CMK）**和**AMD（近似最小度）**重排序算法，重排序后更新所有映射表，**减少稀疏矩阵带宽**（提升缓存命中率和求解效率），支持开启/关闭；

3. **边界DOF处理**：为Dirichlet边界DOF标记属性，后续稀疏矩阵组装时可直接根据标记修改矩阵（行/列置零+对角置1），无需手动筛选；

4. **复DOF标记**：涡流场（复矩阵）时标记所有DOF为复DOF，为后续块稀疏矩阵的复数据存储提供依据。

##### 4. 电磁场景适配细节

- **静电场（标量元-节点DOF）**：遍历网格所有节点，每个节点对应1个全局DOF，四面体单元4个节点→4个DOF，六面体8个节点→8个DOF；

- **静磁场/涡流场（矢量元-边DOF）**：遍历网格所有边，每条边对应1个全局DOF，四面体单元6条边→6个DOF，六面体12条边→12个DOF；

- **周期性边界**：为周期性界面的DOF标记`DofAttr::PERIODIC`，并建立**周期DOF对映射表**，后续矩阵组装时实现周期约束的自动赋值。

### 模块2：DOF集操作（DofSet）

#### 核心目标

封装DOF编号的**集合化操作**，支持大规模DOF下的**并、交、差、补、筛选、批量查询**，为稀疏矩阵的**子矩阵提取、DOF批量筛选、界面DOF聚合**提供高效操作能力。

#### 核心设计方案

##### 1. 双容器底层实现（兼顾效率与易用性）

根据DOF集的**规模和操作类型**自动选择底层容器，无需用户手动指定，平衡**插入/删除/查询**效率和**内存占用**：

- **小规模DOF集（<1e4）**：使用`std::unordered_set<size_t>`，支持O(1)的插入、删除、存在性查询；

- **大规模DOF集（≥1e4）**：使用**位图（BitMask）**（封装为`DynamicBitMask`，基于内存管理工具的预分配数组），支持O(1)的按位操作（并/交/差），内存占用仅为有序容器的1/64（1个bit表示1个DOF）。

##### 2. 核心类设计（`DofSet`，与`DofNumberer`强关联）

```C++

class DofSet {
public:
    // 构造：关联DOF编号器，指定DOF类型（NODE/EDGE）
    DofSet(const DofNumbererBase* dof_num, DofType dof_type);
    // 基本操作：插入、删除、存在性查询
    void insert(size_t global_dof);
    void erase(size_t global_dof);
    bool contains(size_t global_dof) const;
    // 集合操作：并、交、差、补（返回新的DofSet）
    DofSet union_with(const DofSet& other) const;
    DofSet intersect_with(const DofSet& other) const;
    DofSet subtract(const DofSet& other) const;
    DofSet complement() const;
    // 筛选操作（电磁场景定制）：按属性/类型筛选
    DofSet filter_by_attr(DofAttr attr) const; // 如筛选所有界面DOF
    DofSet filter_by_complex(bool is_complex) const;
    // 批量操作：转为有序数组、批量获取局部DOF（分布式）
    std::vector<size_t> to_sorted_vector() const; // 按DOF编号升序，适配矩阵组装
    std::vector<size_t> to_local_dof() const;     // 转换为当前子域的局部DOF
    // 规模查询
    size_t size() const;
    bool empty() const;
private:
    const DofNumbererBase* dof_num_;    // 关联的DOF编号器，获取全局/局部映射
    DofType dof_type_;                  // DOF类型，避免跨类型操作
    std::unique_ptr<DofSetImpl> impl_;  // 底层实现指针（unordered_set/BitMask）
};
```

##### 3. 核心优化点

1. **批量插入/筛选**：提供`insert_range(const std::vector<size_t>& dofs)`接口，避免循环单插入的开销，适配**单元DOF批量加入、子域DOF批量筛选**；

2. **有序输出**：`to_sorted_vector()`默认按DOF编号升序输出，直接适配稀疏矩阵的**行/列索引有序性**（避免COO矩阵转换CSR时的排序开销）；

3. **分布式适配**：`to_local_dof()`自动过滤当前子域无权限的DOF，返回有效局部DOF数组，无无效值；

4. **内存协同**：底层位图/数组由**内存管理工具**管理，实现预分配和缓存行对齐，避免大规模DOF集的内存碎片化。

##### 4. 典型应用场景

- 稀疏矩阵**子矩阵提取**：从全局DOF集中筛选目标区域的DOF集，转为有序数组后直接提取子矩阵；

- **界面DOF聚合**：将多个子域的边界DOF集做并操作，得到全局界面DOF集（FETI-DP核心）；

- **边界DOF筛选**：筛选所有Dirichlet边界DOF集，后续矩阵组装时批量修改矩阵；

- **子域DOF管理**：每个子域对应一个DofSet，管理该子域的所有有效DOF，避免跨子域操作。

### 模块3：界面DOF提取（InterfaceDofExtractor）

#### 核心目标

**基于网格拓扑和子域拆分信息，自动提取子域间的界面DOF**（标量元-节点DOF、矢量元-边DOF），为FETI-DP生成**界面约束矩阵Bᵢ**提供**界面DOF集、界面DOF的主/从子域标记、界面拓扑邻接关系**，是区域分解的核心功能之一。

#### 核心设计方案

##### 1. 界面定义（适配低频电磁网格）

- **物理界面**：两个不同子域的网格**面/边**相交形成的界面，界面上的**节点/边**即为界面DOF（标量元/矢量元）；

- **周期性界面**：物理上重合的周期性边界，其对应的**节点/边**为周期性界面DOF，需单独标记；

- **界面单元**：仅包含界面DOF的虚拟单元，用于描述界面DOF的拓扑邻接关系（适配FETI-DP的界面约束）。

##### 2. 核心类设计（`InterfaceDofExtractor`，依赖`DofNumberer`和`DofSet`）

```C++

// 界面DOF信息（FETI-DP核心，包含拓扑和耦合信息）
struct InterfaceDofInfo {
    size_t global_dof;        // 全局DOF编号
    DofType dof_type;         // DOF类型（NODE/EDGE）
    std::pair<int, int> subdomains; // 归属的两个子域（主/从，-1表示边界）
    size_t interface_id;      // 界面ID，同一界面的DOF标记为同一ID
    DofAttr attr;             // 界面类型（物理/周期性）
};

class InterfaceDofExtractor {
public:
    // 构造：输入网格拓扑、DOF编号器、子域拆分信息
    InterfaceDofExtractor(const MeshTopo* mesh_topo, const DofNumbererBase* dof_num, const SubdomainSplit* sub_split);
    // 核心提取：提取所有子域的界面DOF
    void extract_all_interface_dofs();
    // 批量获取：按子域/界面ID/类型获取界面DOF信息
    std::vector<InterfaceDofInfo> get_interface_dofs_by_subdomain(int sub_id) const;
    std::vector<InterfaceDofInfo> get_interface_dofs_by_type(DofAttr attr = DofAttr::INTERFACE) const;
    // 获取界面DOF集（直接返回DofSet，与稀疏矩阵衔接）
    DofSet get_interface_dof_set(int sub_id) const;
    DofSet get_global_interface_dof_set() const;
    // 获取界面拓扑邻接关系：界面DOF→相邻的界面DOF（适配块稀疏矩阵）
    std::unordered_map<size_t, std::vector<size_t>> get_interface_adjacency() const;
private:
    // 核心提取逻辑（分DOF类型）
    void extract_node_interface_dofs();  // 标量元-节点DOF提取
    void extract_edge_interface_dofs();  // 矢量元-边DOF提取（电磁核心）
    void extract_periodic_interface_dofs(); // 周期性界面DOF提取
    // 主/从子域标记：按子域ID大小或网格拓扑标记，保证分布式一致性
    void mark_master_slave_subdomain(InterfaceDofInfo& info);
    // 存储：界面DOF信息和集合
    std::vector<InterfaceDofInfo> all_interface_dofs_;
    std::unordered_map<int, DofSet> subdomain_interface_sets_; // 子域→界面DOF集
    DofSet global_interface_set_;                              // 全局界面DOF集
};
```

##### 3. 核心提取算法（适配电磁矢量元边DOF）

低频电磁中**矢量元边DOF的界面提取**是重点（静磁/涡流场核心），算法基于**网格面的邻接关系**（界面由子域间的公共面构成），步骤如下：

1. 遍历网格所有**面**，筛选出**跨子域的公共面**（物理界面）或**周期性面对**（周期性界面）；

2. 对每个公共面，提取其**所有关联的边**（矢量元DOF为边），这些边即为**界面边DOF**；

3. 对每个界面边DOF，标记其归属的两个子域，按规则标记**主/从子域**（保证分布式多进程一致）；

4. 对界面边DOF去重（避免公共面的边被重复提取），生成界面DOF信息；

5. 构建界面边DOF的**邻接关系**（基于公共面的拓扑），为块稀疏矩阵的界面块组装提供依据。

##### 4. 关键优化与适配

1. **拓扑遍历优化**：基于网格的**面-边-节点**邻接表遍历，避免暴力遍历所有DOF，时间复杂度O(F)（F为网格面数），远低于O(N)（N为DOF数）；

2. **分布式一致性**：主/从子域标记、界面DOF编号、邻接关系在**所有MPI进程中完全一致**，避免跨进程耦合时的拓扑冲突；

3. **块稀疏适配**：提取界面DOF的**邻接关系**，直接生成界面块稀疏矩阵的**块索引**，无需稀疏矩阵模块再做拓扑计算；

4. **边界DOF过滤**：自动过滤仅属于一个子域的**物理边界DOF**（如模型外边界），不标记为界面DOF，避免FETI-DP耦合错误。

## 三、与现有模块的协同设计

DOF管理工具并非独立模块，需与之前设计的**内存管理工具、稀疏矩阵模块、分布式MPI模块**深度协同，确保整个框架的一致性和高效性，核心协同逻辑如下：

### 1. 与内存管理工具的协同

- **DOF映射表/数组的预分配**：`DofNumberer`中的`elem2global_map_`、`local2global_map_`等核心数组，由`DynamicArrayPreallocator`管理，根据网格规模/子域大小**提前预分配容量**，避免频繁扩容；

- **缓存行对齐**：DOF数组（如映射表、DOF集的位图）由`CacheLineAligner`做**64字节缓存行对齐**，提升DOF查表/集合操作的缓存命中率；

- **内存复用**：分布式场景下，子域DOF映射表在多轮FETI-DP迭代中通过`reset()`接口**复用内存**，减少内存分配/释放次数。

### 2. 与稀疏矩阵模块的协同

- **直接索引赋值**：`DofNumberer`的`elem2global()`接口直接返回全局DOF编号，稀疏矩阵组装时可**直接赋值给COO矩阵的row_idx/col_idx**，无需中间转换；

- **DOF集直接生成子矩阵**：`DofSet`的`to_sorted_vector()`接口返回有序DOF数组，稀疏矩阵模块可根据该数组**直接提取子矩阵**（如界面子矩阵、子域局部矩阵）；

- **界面DOF与块稀疏衔接**：`InterfaceDofExtractor`提取的界面DOF邻接关系，直接作为**块稀疏矩阵的块索引**，支撑界面约束矩阵Bᵢ的块组装。

### 3. 与分布式MPI模块的协同

- **DOF编号同步**：分布式场景下，通过MPI的`Allgather`接口**同步全局DOF规模和界面DOF信息**，确保所有进程的DOF编号一致；

- **界面DOF通信**：`InterfaceDofExtractor`生成的界面DOF集，直接作为MPI通信的**数据索引**，仅传输界面DOF对应的矩阵数据，减少通信量；

- **单进程兼容**：1进程时，MPI通信接口自动屏蔽，DOF管理工具的行为与串行场景完全一致，无代码分支。

## 四、关键设计考虑因素（AI开发需重点关注）

### 1. 电磁DOF类型的强适配，避免通用化过度

- 优先实现**节点DOF（标量元）和边DOF（矢量元）**，无需过早实现面DOF，避免增加复杂度；

- 矢量元边DOF的提取**必须基于网格面的邻接关系**，而非直接遍历边，确保界面DOF的拓扑正确性（电磁矢量元的切向连续约束依赖面拓扑）。

### 2. 大规模DOF下的性能与内存平衡

- 千万级DOF时，**位图（BitMask）**是DOF集的最优选择，内存占用远低于有序容器（如1e8个DOF仅需12.5MB内存）；

- DOF映射表采用**一维连续数组**，避免二维数组的内存碎片化，提升查表效率。

### 3. 拓扑解耦，降低与网格模块的耦合度

- 仅通过**网格拓扑接口**（如`mesh_topo->get_elem_edges(elem_id)`、`mesh_topo->get_face_adjacent_subdomains(face_id)`）获取拓扑信息，不直接访问网格模块的私有成员；

- 封装**网格拓扑适配层**，若后续替换网格模块（如从自研网格换为CGAL/OpenMesh），仅需修改适配层，无需修改DOF管理工具核心逻辑。

### 4. 分布式场景下的一致性保障

- **全局DOF编号的唯一性**：分布式场景下，全局DOF编号由**根进程统一生成**，再通过MPI广播到所有进程，避免各进程独立生成导致的编号冲突；

- **界面DOF信息的同步**：界面DOF的主/从子域标记、邻接关系在根进程生成后广播，确保所有进程的界面信息完全一致；

- **DOF越界检测**：分布式场景下，新增`is_global_dof_in_subdomain()`接口，检测全局DOF是否属于当前子域，避免访问无效DOF。

### 5. 鲁棒性与调试性

- 内置**DOF编号合法性检测**：检测编号是否为非负整数、是否越界、是否重复，检测失败时抛出`std::out_of_range`/`std::invalid_argument`异常，附带明确的错误信息（如“单元123的局部DOF4超出范围，该单元最大局部DOF为3”）；

- 提供**DOF信息打印接口**：打印全局/局部DOF规模、界面DOF数量、边界DOF数量，方便调试；

- 支持**DOF映射表的序列化/反序列化**：将DOF映射表写入文件，方便离线调试和结果复现。

### 6. 扩展性，预留未来功能接口

- 预留**面DOF**接口，适配后续高频电磁或弹性力学的面DOF需求；

- 预留**多物理场DOF耦合**接口，支持电磁-热、电磁-结构的多物理场DOF编号管理；

- 预留**自定义重排序算法**接口，支持后续添加更高效的重排序算法（如METIS）。

## 五、落地要求（融入稀疏矩阵全阶段开发）

### 1. 阶段开发适配（与原6阶段一一对应）

DOF管理工具与稀疏矩阵模块同步开发，分阶段实现功能，前一阶段为后一阶段打基础，无超前开发：

|稀疏矩阵开发阶段|DOF管理工具实现内容|
|---|---|
|阶段1：基础COO/CSR|实现**基础DofNumberer**（节点/边DOF）、单元局部→全局映射、核心数据结构，适配串行标量元DOF|
|阶段2：串行矩阵操作|完善DOF属性标记、边界DOF处理，实现**基础DofSet**（插入/删除/查询），支持子矩阵提取|
|阶段3：电磁场景适配|深度适配**矢量元边DOF**，实现复DOF标记、周期性DOF映射，完善DofSet的电磁定制筛选操作|
|阶段4：分布式矩阵基础|实现**分布式DOF编号**（局部→全局映射），DofSet适配分布式，实现批量局部DOF转换|
|阶段5：FETI-DP核心|实现**InterfaceDofExtractor**全功能，提取界面DOF、生成界面拓扑邻接关系，支撑Bᵢ/F矩阵组装|
|阶段6：集成与工程化|完成工具的工程化封装、Python绑定、文档完善，与有限元框架集成，实现DOF管理的自动化|
### 2. 交付物（独立模块，头文件+源文件+测试用例）

AI开发完成后，交付独立的DOF管理工具模块，与稀疏矩阵模块解耦，可单独编译、测试、维护，核心交付物如下：

```Plain Text

dof_manager/
├── include/
│   ├── dof_type.h        // DOF类型/属性枚举
│   ├── dof_numberer.h    // DOF编号管理类（基础+模板特化）
│   ├── dof_set.h         // DOF集操作类
│   ├── interface_dof_extractor.h // 界面DOF提取类
│   ├── dof_manager.h     // 汇总头文件，一键引入
│   └── mesh_topo_interface.h // 网格拓扑接口（解耦用）
├── src/
│   ├── dof_numberer.cpp
│   ├── dof_set.cpp
│   ├── interface_dof_extractor.cpp
│   └── mesh_topo_interface.cpp
└── test/
    ├── test_dof_numberer.cpp  // 测试DOF编号映射、重排序
    ├── test_dof_set.cpp       // 测试DOF集操作
    ├── test_interface_dof.cpp // 测试界面DOF提取（串行+分布式）
    └── test_dof_em_adapter.cpp// 测试电磁场景DOF适配
```

### 3. 测试验证标准

每个阶段完成后，需通过对应的测试用例，确保功能正确性、性能和鲁棒性，核心测试标准如下：

1. **功能正确性**：DOF映射表查表结果与理论值一致，界面DOF提取与网格拓扑一致，分布式场景下各进程DOF信息同步；

2. **性能**：1e6级DOF下，elem2global查表耗时<1ms，DofSet并操作耗时<5ms，界面DOF提取耗时<10ms；

3. **鲁棒性**：输入非法DOF编号、无效网格拓扑时，能正确抛出异常，无崩溃；

4. **电磁适配**：矢量元边DOF的界面提取结果与Maxwell的界面DOF一致，周期性DOF映射正确；

5. **分布式一致性**：4/8进程下，所有进程的全局DOF规模、界面DOF信息完全一致，无冲突。

## 六、核心API示例（AI开发稀疏矩阵时的典型调用）

为降低AI开发难度，工具提供**高层简洁API**，稀疏矩阵组装时无需关注底层拓扑细节，典型调用示例如下（以**静磁场矢量元边DOF**为例）：

```C++

// 1. 初始化网格拓扑、DOF编号器、界面DOF提取器
MeshTopo mesh_topo("mesh.vtk"); // 读取网格
DofNumberer<DofType::EDGE> dof_num(&mesh_topo); // 边DOF编号器
dof_num.generate_dofs(); // 生成全局DOF
dof_num.reorder(ReorderType::CMK); // CMK重排序
InterfaceDofExtractor if_extractor(&mesh_topo, &dof_num, &sub_split); // 子域拆分信息
if_extractor.extract_all_interface_dofs(); // 提取界面DOF

// 2. 稀疏矩阵组装：单元局部DOF→全局DOF，直接赋值
CooMatrix mat;
for (size_t elem_id = 0; elem_id < mesh_topo.elem_num(); ++elem_id) {
    Matrix Ke = compute_elem_matrix(elem_id); // 计算单元矩阵
    for (size_t i = 0; i < 6; ++i) { // 四面体矢量元6个局部DOF
        size_t global_i = dof_num.elem2global(elem_id, i); // 局部→全局
        for (size_t j = 0; j < 6; ++j) {
            size_t global_j = dof_num.elem2global(elem_id, j);
            mat.add_value(global_i, global_j, Ke(i, j)); // 直接组装COO矩阵
        }
    }
}

// 3. 提取子域界面DOF集，生成界面约束矩阵B_i
DofSet interface_dofs = if_extractor.get_interface_dof_set(sub_id);
std::vector<size_t> if_dof_vec = interface_dofs.to_sorted_vector();
CsrMatrix B_i = generate_interface_constraint(if_dof_vec); // 生成B_i矩阵
```
> （注：文档部分内容可能由 AI 生成）