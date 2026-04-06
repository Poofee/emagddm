# 全局装配模块（Assembly Module）Spec

## Why
有限元求解器的核心流程是：网格→形函数→单元矩阵→**全局装配**→线性求解。当前已完成前三个模块的开发，但缺少将单元级矩阵组装为全局系统矩阵的关键桥梁。没有装配模块，无法构建完整的离散代数方程组 Kx = F，也无法进入求解环节。本模块是实现"端到端闭环"的必要步骤。

## What Changes
- **新增** 全局装配核心类 `EMAssembly`（位于 numeric 命名空间）
  - 支持 COO 格式组装 + CSR 格式输出（工业级标准）
  - 支持静态场（K + F）和瞬态场（K + M + C + F）
  - 自动适配所有已支持的单元类型（标量/棱边/混合AV、一阶/二阶）
  - 完全复用已有的 EMMeshData、Local2Global、EMElementIntegratorBase 接口
- **新增** 装配性能统计接口（耗时、非零元素数量）
- **预留** OpenMP 多线程并行接口（单元级并行）
- **零修改** 已有代码（网格、形函数、DOF管理、单元矩阵模块）

## Impact
- Affected specs: 无（独立新模块）
- Affected code:
  - `include/numeric/em_assembly.hpp` (新建)
  - `src/numeric/em_assembly.cpp` (新建)
  - `test/test_assembly.cpp` (新建测试)
  - 复用: `fe_em::EMMeshData`, `fe_em::Local2Global`, `numeric::EMElementIntegratorBase`, `numeric::CooMatrix`, `numeric::CsrMatrix`

## ADDED Requirements

### Requirement: 全局装配核心类 EMAssembly

系统 SHALL 提供 `EMAssembly` 类，实现从单元级矩阵到全局系统矩阵的组装功能。

#### 核心输入接口

```cpp
class EMAssembly {
public:
    /**
     * @brief 执行全局装配主流程
     * @param mesh_data 网格拓扑数据（节点、单元）
     * @param elem_l2g 每个单元的局部-全局DOF映射表
     * @param materials 材料参数列表（按material_id索引）
     * @param num_free_dofs 全局自由DOF总数
     * @param is_transient 是否为瞬态场（默认false=静态场）
     */
    void assemble(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<fe_em::Local2Global>& elem_l2g,
        const std::map<int, numeric::MaterialProperties>& materials,
        int num_free_dofs,
        bool is_transient = false
    );
    
    /**
     * @brief 获取全局刚度矩阵（CSR格式）
     * @return const CsrMatrix<double>& 刚度矩阵引用
     */
    const CsrMatrix<double>& getStiffnessMatrix() const;
    
    /**
     * @brief 获取全局质量矩阵（CSR格式，仅瞬态场有效）
     * @return const CsrMatrix<double>& 质量矩阵引用
     */
    const CsrMatrix<double>& getMassMatrix() const;
    
    /**
     * @brief 获取全局阻尼矩阵（CSR格式，仅瞬态场有效）
     * @return const CsrMatrix<double>& 阻尼矩阵引用
     */
    const CsrMatrix<double>& getDampingMatrix() const;
    
    /**
     * @brief 获取全局右端项向量
     * @return const Eigen::VectorXd& 右端项向量引用
     */
    const Eigen::VectorXd& getSourceVector() const;
    
    /**
     * @brief 获取装配性能统计信息
     * @return AssemblyStats 性能统计结构体
     */
    AssemblyStats getStats() const;
    
    /**
     * @brief 清空所有装配结果
     */
    void clear();
};
```

#### Scenario: 静态场单单元无约束装配成功案例

- **GIVEN** 一个 PRISM6 单元（MIXED_AV类型），6个节点+9条棱边，共15个局部DOF
- **AND** Local2Global 映射表所有 indices[i] >= 0（无约束DOF）
- **AND** 材料参数已设置（μ_r=1.0, ε_r=1.0, σ=0）
- **WHEN** 调用 assemble() 方法（is_transient=false）
- **THEN** 输出刚度矩阵 K 的尺寸为 15×15
- **AND** K 的非零元素数量与单元 Ke 矩阵一致
- **AND** K 矩阵对称（|K(i,j) - K(j,i)| < 1e-10）
- **AND** 右端项向量 F 的长度为 15
- **AND** 统计信息显示装配完成且无错误

#### Scenario: 瞬态场多单元混合装配成功案例

- **GIVEN** 多个不同类型的单元（PRISM6_MIXED_AV + PYRAMID5_EDGE + TET4_SCALAR）
- **AND** 每个单元有对应的 Local2Global 映射表（部分DOF被约束，indices=-1）
- **AND** is_transient=true
- **WHEN** 调用 assemble() 方法
- **THEN** 输出 K、M、C 矩阵的尺寸均为 num_free_dofs × num_free_dofs
- **AND** 所有约束DOF对应的行/列在全局矩阵中不存在（已被跳过）
- **AND** K、M、C 矩阵均保持对称性
- **AND** F 向量的长度为 num_free_dofs

### Requirement: COO格式组装与CSR转换

系统 SHALL 采用 "COO坐标格式累加 → 合并重复项 → CSR压缩输出" 的标准工业流程。

#### 核心实现逻辑

```
1. 初始化阶段：
   - 创建 CooMatrix<double> K_coo(num_free_dofs, num_free_dofs)
   - 若瞬态：创建 M_coo, C_coo（同尺寸）
   - 初始化 Eigen::VectorXd F = Eigen::VectorXd::Zero(num_free_dofs)
   - 预分配COO内存：capacity = 预估总非零元数（基于单元数×平均DOF²）

2. 单元遍历循环（for each element in mesh_data.elements）：
   a. 获取 Local2Global l2g = elem_l2g[elem.id]
   b. 提取单元节点坐标：node_coords = buildNodeCoords(elem, mesh_data.nodes)
   c. 获取材料参数：mat = materials[elem.material_id]
   d. 调用单元积分器：
      - 静态场：Ke = integrator->computeStiffnessMatrix(node_coords)
                 Fe = integrator->computeSourceVector(node_coords)
      - 瞬态场：mats = integrator->computeAllMatrices(node_coords)
                 （包含Ke, Me, Ce, Fe）
   
   e. 矩阵 Scatter 双重循环（i=0..n_dof-1, j=0..n_dof-1）：
      i.   global_row = l2g.indices[i]
      ii.  global_col = l2g.indices[j]
      iii. 约束跳过规则：若 global_row == -1 或 global_col == -1 → continue
      iv.  K_coo.add_value(global_row, global_col, Ke(i,j))
      v.   瞬态场额外：M_coo.add_value(global_row, global_col, Me(i,j))
                       C_coo.add_value(global_row, global_col, Ce(i,j))
   
   f. 向量 Scatter 单循环（i=0..n_dof-1）：
      i.   global_row = l2g.indices[i]
      ii.  约束跳过：若 global_row == -1 → continue
      iii. F(global_row) += Fe(i)

3. 后处理阶段：
   - K_csr.build_from_coo(K_coo)  // 合并重复项并转换为CSR
   - 瞬态场：M_csr.build_from_coo(M_coo), C_csr.build_from_coo(C_coo)

4. 输出结果
```

#### Scenario: COO到CSR转换正确性验证

- **GIVEN** 两个相邻单元共享某些全局DOF
- **WHEN** 分别将两个单元的矩阵scatter到同一个COO矩阵
- **THEN** COO矩阵中相同位置的值自动累加
- **AND** 转换后的CSR矩阵在共享DOF位置包含两个单元的贡献之和
- **AND** CSR矩阵的 row_ptr、col_indices、values 数组格式正确

### Requirement: 约束DOF处理机制

系统 SHALL 严格遵循 Local2Global 映射表的约定，正确处理受约束的自由度。

#### 约束处理规则

- **规则1**: `l2g.indices[i] == -1` 表示第 i 个局部DOF被 Dirichlet 边界条件约束
- **规则2**: 约束DOF不参与任何矩阵元素的累加（行和列均跳过）
- **规则3**: 约束DOF不参与右端项向量的累加
- **规则4**: 最终输出的全局矩阵尺寸严格等于 num_free_dofs（不含约束DOF）

#### Scenario: 约束DOF完全跳过验证

- **GIVEN** 一个 PYRAMID5_EDGE 单元（8条棱边DOF）
- **AND** Local2Global 中标记第3条棱边为约束（indices[3] = -1）
- **WHEN** 执行装配
- **THEN** 全局矩阵的尺寸为 7×7（8-1=7个自由DOF）
- **AND** 全局矩阵中不存在任何与第3条棱边相关的行或列
- **AND** 对称性保持完好

### Requirement: 内存预分配优化

系统 SHALL 支持基于预估的非零元素数量预分配COO矩阵内存，避免运行时动态扩容。

#### 预估算法

```cpp
int estimateNNZ(const std::vector<Element>& elements, int avg_dofs_per_elem) {
    // 保守估计：每个单元贡献 avg_dofs_per_elem² 个非零元
    // 考虑到对称性和约束跳过，实际可能更少
    return static_cast<int>(elements.size()) * avg_dofs_per_elem * avg_dofs_per_elem;
}
```

#### Scenario: 大规模模型内存预分配有效性

- **GIVEN** 包含10000个单元的中等规模模型
- **WHEN** 使用预分配策略执行装配
- **THEN** COO矩阵在装配过程中无内存重新分配（reserve生效）
- **AND** 装配性能相比未预分配提升20%以上（减少realloc开销）

### Requirement: 装配性能统计接口

系统 SHALL 提供详细的装配过程统计信息，用于性能分析和调试。

#### 统计结构体定义

```cpp
struct AssemblyStats {
    double assembly_time_ms;          // 总装配耗时（毫秒）
    int total_elements;               // 处理的单元总数
    int free_dofs;                    // 自由DOF数量
    int constrained_dofs;             // 约束DOF数量
    int nnz_stiffness;                // 刚度矩阵非零元数
    int nnz_mass;                     // 质量矩阵非零元数（瞬态场）
    int nnz_damping;                  // 阻尼矩阵非零元数（瞬态场）
    bool is_symmetric_k;              // 刚度矩阵是否对称
    bool is_symmetric_m;              // 质量矩阵是否对称（瞬态场）
    bool is_symmetric_c;              // 阻尼矩阵是否对称（瞬态场）
    std::string error_message;        // 错误信息（空字符串表示成功）
};
```

#### Scenario: 统计信息完整性验证

- **WHEN** 装配完成后调用 getStats()
- **THEN** 返回完整的统计信息
- **AND** assembly_time_ms > 0（记录了实际耗时）
- **AND** total_elements 与输入单元数一致
- **AND** nnz_stiffness > 0（至少有对角线元素）
- **AND** is_symmetric_k == true（刚度矩阵必须对称）

### Requirement: 单元积分器适配层

系统 SHALL 通过工厂模式或策略模式适配不同类型的单元积分器。

#### 适配接口设计

```cpp
/**
 * @brief 根据单元类型和DOF类型创建对应的积分器实例
 * @param elem_type 单元类型枚举
 * @param dof_type DOF类型枚举
 * @return unique_ptr<EMElementIntegratorBase> 积分器智能指针
 */
static std::unique_ptr<EMElementIntegratorBase> createIntegrator(
    ElemType elem_type,
    DOFType dof_type
);
```

#### 支持的单元-积分器映射表

| 单元类型 | DOFType | 积分器类型 |
|---------|---------|-----------|
| TRI3/TRI6/QUAD4/.../HEX27/PRISM6/PRISM15/PYRAMID5/PYRAMID13 | SCALAR_ONLY | LagrangeEMIntegrator |
| TET4_EDGE/HEX8_EDGE/PRISM6_EDGE/PYRAMID5_EDGE | VECTOR_EDGE_ONLY | NedelecEMIntegrator |
| 任意标量单元 | MIXED_AV | LagrangeEMIntegrator(标量部分) + NedelecEMIntegrator(矢量部分) |

#### Scenario: 自动选择正确的积分器

- **GIVEN** 一个 HEX8 单元，dof_type = SCALAR_ONLY
- **WHEN** 调用 createIntegrator(HEX8, SCALAR_ONLY)
- **THEN** 返回 LagrangeEMIntegrator 实例
- **AND** 该实例可正确计算 8 个节点的单元矩阵

### Requirement: OpenMP并行预留接口

系统 SHALL 预留OpenMP多线程并行装配的接口框架（可选实现）。

#### 并行设计要点

- **粒度**: 单元级并行（每个线程处理一批单元）
- **冲突避免**: 
  - 方案A（推荐）: 每个线程维护私有COO矩阵，最后合并（线程安全）
  - 方案B（备选）: 单元颜色重排（Graph coloring，确保同色单元无共享DOF）
- **接口预留**:

```cpp
/**
 * @brief 设置并行线程数（默认1=串行）
 * @param num_threads 线程数量，<=1时退化为串行
 */
void setNumThreads(int num_threads);

/**
 * @brief 获取当前并行线程数
 * @return int 当前线程数
 */
int getNumThreads() const;
```

#### Scenario: 并行装配向后兼容

- **GIVEN** 默认配置（num_threads=1）
- **WHEN** 执行装配
- **THEN** 行为与串行版本完全一致（结果数值误差 < 1e-12）
- **AND** 无数据竞争或死锁风险

## MODIFIED Requirements

无（纯新增模块，不修改已有代码）

## REMOVED Requirements

无
