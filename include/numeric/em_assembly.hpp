/**
 * @file em_assembly.hpp
 * @brief 数值计算层 - 电磁场有限元全局组装器头文件
 * @details 实现电磁场有限元分析的全局系统矩阵/向量组装功能，
 *          负责将单元级矩阵（K_e, M_e, C_e）和源项向量（F_e）通过scatter操作
 *          组装为全局刚度矩阵（K）、质量矩阵（M）、阻尼矩阵（C）和右端项向量（F）。
 *
 *          核心组装流程：
 *          1. 预估非零元素数量并预分配COO矩阵内存
 *          2. 遍历所有单元，调用EMElementIntegratorBase计算单元矩阵
 *          3. 通过Local2Global映射将单元矩阵scatter到全局COO矩阵
 *          4. 将COO格式转换为CSR格式以支持高效求解
 *          5. 统计装配信息并进行对称性检查
 *
 *          适用场景：
 *          - 静态电磁场：组装K和F矩阵/向量
 *          - 瞬态涡流场：组装K、M、C、F全量矩阵/向量
 *          - 频域波动场：组装K、M、C、F全量矩阵/向量（复数形式）
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 *
 * @note 与其他模块的协作关系：
 *       - 上游依赖：EMMeshData（网格拓扑）、EMElementIntegratorBase（单元积分器）、Local2Global（DOF映射）
 *       - 下游服务：线性求解器（提供CSR格式矩阵）、后处理模块（提供解向量）
 *       - 并行策略：预留OpenMP多线程接口，支持单元级并行组装
 */

#pragma once

#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "em_element_integrator_base.hpp"
#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace numeric {

// ==================== 组装统计信息结构体 ====================

/**
 * @struct AssemblyStats
 * @brief 有限元全局组装统计信息结构体
 * @details 记录一次完整组装过程的性能指标和结果统计，
 *          用于性能分析、调试诊断和质量验证。
 *
 *          统计维度说明：
 *          - 时间维度：总装配耗时（毫秒级精度）
 *          - 规模维度：单元总数、自由/约束DOF数量
 *          - 稀疏性维度：各矩阵的非零元素数量
 *          - 对称性维度：各矩阵的对称性标志（用于选择最优求解器）
 *          - 错误维度：错误信息字符串（空表示成功）
 */
struct AssemblyStats {
    double assembly_time_ms;          ///< 总装配耗时（毫秒），从assemble()开始到结束的挂钟时间
    int total_elements;               ///< 处理的单元总数（等于网格单元数，跳过无效单元时可能小于此值）
    int free_dofs;                    ///< 自由DOF数量（参与求解的未知量数，决定系统矩阵维度）
    int constrained_dofs;             ///< 约束DOF数量（受Dirichlet等边界条件约束的DOF数）
    int nnz_stiffness;                ///< 刚度矩阵K的非零元素数量（COO去重后的实际非零数）
    int nnz_mass;                     ///< 质量矩阵M的非零元素数量（瞬态场专用，静态场为0）
    int nnz_damping;                  ///< 阻尼矩阵C的非零元素数量（瞬态场专用，静态场为0）
    bool is_symmetric_k;              ///< 刚度矩阵K是否对称（true时可使用Cholesky分解加速求解）
    bool is_symmetric_m;              ///< 质量矩阵M是否对称（瞬态场专用）
    bool is_symmetric_c;              ///< 阻尼矩阵C是否对称（瞬态场专用）
    std::string error_message;        ///< 错误信息（空字符串表示组装成功，否则包含具体错误描述）

    /**
     * @brief 默认构造函数，初始化为零/默认状态
     * @details 所有数值字段初始化为0，布尔字段初始化为false，
     *          error_message初始化为空字符串（表示成功状态）
     */
    AssemblyStats()
        : assembly_time_ms(0.0)
        , total_elements(0)
        , free_dofs(0)
        , constrained_dofs(0)
        , nnz_stiffness(0)
        , nnz_mass(0)
        , nnz_damping(0)
        , is_symmetric_k(false)
        , is_symmetric_m(false)
        , is_symmetric_c(false)
        , error_message()
    {
    }
};

// ==================== 全局组装器类定义 ====================

/**
 * @class EMAssembly
 * @brief 电磁场有限元全局组装器类
 * @details 实现从单元级矩阵到全局系统矩阵/向量的完整组装流程，
 *          是数值计算层连接单元积分器和线性求解器的核心桥梁。
 *
 *          设计原则与架构定位：
 *          - 单一职责：专注于矩阵/向量组装，不涉及求解和后处理
 *          - 策略模式：通过EMElementIntegratorBase虚函数分派到具体单元类型的积分算法
 *          - 内存优化：采用COO→CSR两阶段组装，避免直接操作CSR的插入开销
 *          - 可扩展性：预留OpenMP并行接口，支持大规模问题的多线程加速
 *
 *          组装算法流程（assemble方法）：
 *          ```
 *          输入: 网格数据(mesh)、DOF映射(dof_mapping)、材料参数(materials)
 *          输出: 全局矩阵(K,M,C) + 全局向量(F)
 *
 *          Step 1: 预处理
 *            - 估算非零元素数量 estimateNNZ()
 *            - 预分配COO矩阵内存（避免动态扩容）
 *
 *          Step 2: 单元循环（可并行化）
 *            for each element e in mesh:
 *              - 创建积分器 createIntegrator(elem_type, dof_type)
 *              - 构建节点坐标 buildNodeCoords(e, nodes)
 *              - 计算单元矩阵 integrator->computeAllMatrices(coords)
 *              - Scatter到全局COO scatterElementMatrix(coo, K_e, l2g)
 *              - Scatter源项向量 scatterElementVector(F, F_e, l2g)
 *
 *          Step 3: 后处理
 *            - COO → CSR转换 csr.build_from_coo(coo)
 *            - 对称性检查 checkSymmetry()
 *            - 统计信息收集
 *          ```
 *
 *          时间复杂度分析：
 *          - 设n_elem为单元数，n_dof_per_elem为平均每单元DOF数
 *          - 单元矩阵计算: O(n_elem × n_dof_per_elem³) （高斯积分主导）
 *          - Scatter操作: O(n_elem × n_dof_per_elem²) （稀疏累加）
 *          - COO→CSR排序: O(nnz × log(nnz)) （快速排序）
 *          - 总体复杂度: O(n_elem × n_dof_per_elem³) （单元积分主导）
 *
 *          内存占用估算：
 *          - COO临时存储: O(nnz) （3个数组：row, col, value）
 *          - CSR最终存储: O(nnz + n_rows + 1) （压缩行存储）
 *          - 峰值内存: 约2×nnz（COO+CSR共存阶段）
 *
 * @note 此类位于numeric命名空间下，属于数值计算层的核心模块
 * @see EMElementIntegratorBase 单元积分器抽象基类
 * @see AssemblyStats 组装统计信息结构体
 * @see CooMatrix COO格式稀疏矩阵
 * @see CsrMatrix CSR格式稀疏矩阵
 * @see fe_em::Local2Global 局部-全局DOF映射
 */
class EMAssembly {
public:
    // ==================== 构造与析构 ====================

    /**
     * @brief 默认构造函数
     * @details 初始化所有内部状态为空/零值，
     *          设置默认线程数为1（单线程模式）
     */
    EMAssembly();

    /**
     * @brief 析构函数
     * @details 释放所有内部资源（矩阵、向量等自动析构）
     */
    ~EMAssembly() = default;

    // ==================== 核心组装接口 ====================

    /**
     * @brief 执行完整的全局系统矩阵/向量组装
     * @param mesh_data 网格拓扑数据（节点、单元、边界标记），常量引用避免拷贝
     * @param elem_local_to_global 所有单元的局部-全局DOF映射表，由EMDOFManager构建
     * @param materials 所有材料的属性参数列表（索引对应Element::material_id）
     * @param is_transient 是否启用瞬态模式（true时额外组装M和C矩阵）
     * @return bool 组装成功返回true，失败返回false（详细错误记录在stats_.error_message）
     *
     * @details 这是组装器的核心入口方法，完成以下完整流程：
     *
     *          **第一阶段：预处理**
     *          1. 清空上次的组装结果（K_csr_, M_csr_, C_csr_, F_重置为空）
     *          2. 重置统计信息stats_为初始状态
     *          3. 估算全局矩阵的非零元素数量（基于单元连通性启发式算法）
     *          4. 预分配COO矩阵内存（避免组装过程中的动态扩容开销）
     *
     *          **第二阶段：单元级循环组装**
     *          对每个单元e执行：
     *          a. 查询单元的Local2Global映射l2g[e.id]
     *          b. 根据elem_type和dof_type创建对应的EMElementIntegratorBase派生类实例
     *          c. 从mesh_data.nodes提取该单元的节点坐标，构建dim×n_nodes坐标矩阵
     *          d. 调用integrator->setMaterialProperties(materials[e.material_id])设置材料参数
     *          e. 若is_transient=true，调用integrator->setTransientMode(true)启用瞬态模式
     *          f. 调用integrator->computeAllMatrices(node_coords)获取单元矩阵(K_e, M_e, C_e, F_e)
     *          g. 将K_e scatter到全局COO刚度矩阵（跳过indices[i]==-1的约束DOF）
     *          h. 若瞬态模式，将M_e和C_e分别scatter到对应的COO矩阵
     *          i. 将F_e scatter到全局右端项向量F_
     *
     *          **第三阶段：后处理与格式转换**
     *          1. 将COO格式的K转换为CSR格式（排序+压缩）
     *          2. 若瞬态模式，同样转换M和C
     *          3. 对每个CSR矩阵进行对称性检查（用于后续求解器选择）
     *          4. 填充stats_统计信息（时间、规模、稀疏性、对称性）
     *          5. 记录日志输出（INFO级别汇总，DEBUG级别详细信息）
     *
     *          **约束DOF处理规则**（关键逻辑）：
     *          - Local2Global::indices[i] == -1 表示第i个局部DOF为约束DOF
     *          - Scatter时跳过所有涉及约束DOF的矩阵元素（不写入COO）
     *          - 这实现了Dirichlet边界条件的"消去法"处理（而非罚函数法）
     *          - 最终CSR矩阵的维度 = 自由DOF数（不含约束DOF）
     *
     * @note 典型调用示例：
     *       @code
     *       EMAssembly assembler;
     *       assembler.setNumThreads(4);  // 启用4线程并行
     *       bool success = assembler.assemble(
     *           mesh_data,           // 网格拓扑
     *           dof_mappings,        // DOF映射表
     *           material_list,       // 材料参数列表
     *           true                // 瞬态模式
     *       );
     *       if (!success) {
     *           auto stats = assembler.getStats();
     *           std::cerr << "组装失败: " << stats.error_message << std::endl;
     *       }
     *       const auto& K = assembler.getStiffnessMatrix();  // 获取CSR刚度矩阵
     *       @endcode
     *
     * @warning 调用前需确保：
     *          - mesh_data不为空且已正确填充节点和单元数据
     *          - elem_local_to_global的大小等于mesh_data.getElementCount()
     *          - materials列表覆盖了所有单元引用的material_id
     *          - 若启用瞬态模式，材料参数中的sigma（电导率）应合理设置
     *
     * @exception 若输入数据非法（如空网格、material_id越界），返回false并设置error_message
     * @exception 若内存分配失败（预估NNZ过大），返回false并设置error_message
     */
    bool assemble(
        const fe_em::EMMeshData& mesh_data,
        const std::vector<fe_em::Local2Global>& elem_local_to_global,
        const std::vector<MaterialProperties>& materials,
        bool is_transient = false
    );

    // ==================== 结果访问接口 ====================

    /**
     * @brief 获取全局刚度矩阵（只读访问）
     * @return const CsrMatrix<double>& CSR格式的全局刚度矩阵常量引用
     *
     * @details 返回组装完成后的全局刚度矩阵K，
     *          矩阵维度为 free_dofs × free_dofs（已消去约束DOF）。
     *
     *          物理含义（以A-φ公式体系为例）：
     *          - K_ij = ∫_Ω (1/μ)(∇×N_i)·(∇×N_j) dΩ  （磁阻率curl-curl算子）
     *          - 或 K_ij = ∫_Ω ε∇N_i·∇N_j dΩ          （介电常数梯度算子）
     *
     * @note 必须在assemble()成功调用后方可使用，否则返回空矩阵
     * @warning 返回值为const引用，外部代码不得修改矩阵内容
     */
    const CsrMatrix<double>& getStiffnessMatrix() const;

    /**
     * @brief 获取全局质量矩阵（只读访问）
     * @return const CsrMatrix<double>& CSR格式的全局质量矩阵常量引用
     *
     * @details 返回组装完成后的全局质量矩阵M，
     *          仅在瞬态模式下（assemble()传入is_transient=true）有有效数据。
     *
     *          物理含义：
     *          - M_ij = ∫_Ω σN_i·N_j dΩ  （电导率质量项，关联∂A/∂t或∂²φ/∂t²）
     *          - 或 M_ij = ∫_Ω εN_i·N_j dΩ （介电常数质量项，波动方程惯性项）
     *
     * @note 仅瞬态模式下此矩阵有效，静态模式下返回空矩阵
     * @warning 返回值为const引用，外部代码不得修改矩阵内容
     */
    const CsrMatrix<double>& getMassMatrix() const;

    /**
     * @brief 获取全局阻尼矩阵（只读访问）
     * @return const CsrMatrix<double>& CSR格式的全局阻尼矩阵常量引用
     *
     * @details 返回组装完成后的全局阻尼矩阵C，
     *          仅在瞬态模式下有有效数据。
     *
     *          物理含义：
     *          - C_ij = ∫_Ω σN_i·N_j dΩ  （涡流耗散阻尼项）
     *          - 在某些公式体系中，C可能与M使用相同核函数但具有不同系数
     *
     * @note 仅瞬态模式下此矩阵有效，静态模式下返回空矩阵
     * @warning 返回值为const引用，外部代码不得修改矩阵内容
     */
    const CsrMatrix<double>& getDampingMatrix() const;

    /**
     * @brief 获取全局右端项源向量（只读访问）
     * @return const Eigen::VectorXd& 全局源项向量常量引用
     *
     * @details 返回组装完成后的全局右端项向量F，
     *          向量长度为free_dofs（已消去约束DOF的贡献）。
     *
     *          物理含义：
     *          - F_i = ∫_Ω N_i·J_src dΩ  （源电流密度投影）
     *          - 或 F_i = ∫_Γ N_i·q_src dΓ （边界源通量投影）
     *
     * @note 必须在assemble()成功调用后方可使用，否则返回空向量
     * @warning 返回值为const引用，外部代码不得修改向量内容
     */
    const Eigen::VectorXd& getSourceVector() const;

    // ==================== 统计与状态查询 ====================

    /**
     * @brief 获取本次组装的统计信息
     * @return AssemblyStats 组装统计信息结构体副本
     *
     * @details 返回最近一次assemble()调用产生的完整统计数据，
     *          包括时间性能、系统规模、稀疏特性、对称性和错误信息。
     *
     *          典型应用场景：
     *          - 性能分析：比较不同网格规模下的assembly_time_ms
     *          - 质量验证：检查is_symmetric_k确认刚度矩阵是否符合理论预期
     *          - 调试诊断：通过error_message定位组装失败原因
     *          - 日志记录：将stats序列化输出到文件供后续分析
     *
     * @note 若尚未调用过assemble()，返回默认构造的全零统计信息
     */
    AssemblyStats getStats() const;

    // ==================== 资源管理 ====================

    /**
     * @brief 清空所有内部数据，释放矩阵和向量占用的内存
     * @details 重置组装器到初始状态，包括：
     *          - 清空K_csr_、M_csr_、C_csr_三个CSR矩阵（释放内部vector内存）
     *          - 清空F_源项向量（释放Eigen向量内存）
     *          - 重置stats_为默认构造状态
     *          - 保持num_threads_不变（线程数设置独立于数据生命周期）
     *
     * @note 典型使用场景：
     *          - 重新组装前释放上次结果的内存（避免峰值内存过高）
     *          - 组装器对象复用（同一对象组装多个不同问题）
     *          - 异常恢复（组装失败后清理中间状态）
     *
     * @warning clear()后需重新调用assemble()才能获取有效的矩阵/向量
     */
    void clear();

    // ==================== 并行控制接口（OpenMP预留）====================

    /**
     * @brief 设置并行组装的线程数量
     * @param num_threads 线程数（必须>=1，1表示单线程模式）
     *
     * @details 配置OpenMP并行区域的最大线程数，
     *          影响assemble()中单元循环的并行度。
     *
     *          当前实现说明（预留接口）：
     *          - 内部存储num_threads_成员变量
     *          - 实际并行化需配合OpenMP pragma实现
     *          - 单元级并行需确保COO矩阵的线程安全写入（原子操作或归约）
     *
     *          性能预期（理想情况）：
     *          - 加速比 ≈ min(num_threads, n_elem / threshold)
     *          - 其中threshold为线程启动开销平衡点（通常100~1000个单元）
     *          - 内存带宽受限时加速比会下降（NUMA效应、缓存竞争）
     *
     * @note 推荐设置值：
     *          - 小规模问题（<10k单元）：1（单线程，避免并行开销）
     *          - 中等规模（10k~1M单元）：物理核心数
         *          - 大规模（>1M单元）：物理核心数（考虑NUMA绑定）
     *
     * @warning 若num_threads < 1，输出WARN日志并强制设置为1
     */
    void setNumThreads(int num_threads);

    /**
     * @brief 获取当前配置的线程数量
     * @return int 当前线程数（>=1）
     *
     * @details 返回通过setNumThreads()设置的线程数，
     *          或默认值1（若未显式设置）。
     */
    int getNumThreads() const;

private:
    // ==================== 私有成员变量 ====================

    CsrMatrix<double> K_csr_;      ///< 全局刚度矩阵（CSR格式，组装最终结果）
    CsrMatrix<double> M_csr_;      ///< 全局质量矩阵（CSR格式，瞬态场专用）
    CsrMatrix<double> C_csr_;      ///< 全局阻尼矩阵（CSR格式，瞬态场专用）
    Eigen::VectorXd F_;            ///< 全局右端项源向量（长度=free_dofs）
    AssemblyStats stats_;          ///< 本次组装的统计信息
    int num_threads_;              ///< 并行线程数（默认1，OpenMP预留）

    // ==================== 私有辅助方法 ====================

    /**
     * @brief 工厂方法：根据单元类型和DOF类型创建对应的单元积分器实例
     * @param elem_type 单元类型枚举（fe_em::ElemType，即numeric::ElementType）
     * @param dof_type DOF类型枚举（fe_em::DOFType，决定标量/矢量/混合）
     * @return std::unique_ptr<EMElementIntegratorBase> 指向新创建的积分器实例的智能指针
     *
     * @details 根据单元几何类型和物理场DOF类型的具体组合，
     *          实例化对应的EMElementIntegratorBase派生类：
     *
     *          | 单元类型 | DOFType | 创建的积分器类 |
     *          |---------|---------|---------------|
     *          | TRIANGLE / QUAD | SCALAR_ONLY | LagrangeEMIntegrator (2D标量) |
     *          | TETRA / HEX | SCALAR_ONLY | LagrangeEMIntegrator (3D标量) |
     *          | TRIANGLE / QUAD | VECTOR_EDGE_ONLY | NedelecEMIntegrator (2D矢量) |
     *          | TETRA / HEX | VECTOR_EDGE_ONLY | NedelecEMIntegrator (3D矢量) |
     *          | TETRA / HEX | MIXED_AV | MixedAVIntegrator (3D混合) |
     *
     *          设计模式：采用简单工厂模式（Simple Factory），
     *          通过if-else分支根据类型组合创建具体对象。
     *          未来可重构为工厂注册表模式以支持用户自定义积分器。
     *
     * @note 返回的智能指针自动管理积分器生命周期（RAII）
     * @warning 若传入不支持的类型组合，返回nullptr并输出ERROR日志
     */
    static std::unique_ptr<EMElementIntegratorBase> createIntegrator(
        fe_em::ElemType elem_type,
        fe_em::DOFType dof_type
    );

    /**
     * @brief 从全局节点列表构建单个单元的节点坐标矩阵
     * @param element 单元数据的常量引用（包含node_ids列表）
     * @param nodes 全局节点列表的常量引用
     * @return Eigen::MatrixXd 节点坐标矩阵，维度为 dim × node_count
     *
     * @details 根据单元的node_ids从全局nodes数组中提取对应节点的坐标，
     *          组装成列优先存储的坐标矩阵供积分器使用。
     *
     *          输出矩阵格式约定：
     *          - 行数 = 问题维数（2D问题为2，3D问题为3）
     *          - 列数 = 单元节点数（如TETRA为4，HEX为8）
     *          - 第j列 = 第j个节点的(x, y, z)^T坐标向量
     *          - 示例（3D四节点四面体）：
     *            ```
     *            coords = [x1 x2 x3 x4]
     *                      [y1 y2 y3 y4]
     *                      [z1 z2 z3 z4]
     *            ```
     *
     * @note 坐标单位为米（m），与材料参数的SI单位制一致
     * @warning 若element.node_ids包含无效索引（越界或负数），抛出std::out_of_range异常
     */
    static Eigen::MatrixXd buildNodeCoords(
        const fe_em::Element& element,
        const std::vector<fe_em::Node>& nodes
    );

    /**
     * @brief 估算全局矩阵的非零元素数量（用于COO矩阵预分配）
     * @param elements 全局单元列表的常量引用
     * @param avg_dof_per_elem 平均每个单元的自由度数（估计值）
     * @return int 估算的非零元素数量（上界，实际NNZ <= 此值）
     *
     * @details 基于单元连通性的启发式算法估算全局稀疏矩阵的NNZ，
     *          用于预先分配COO矩阵的内部vector容量，避免动态扩容的性能损失。
     *
     *          估算公式（保守上界）：
     *          NNZ_estimate = n_elements × (avg_dof_per_elem)² × overlap_factor
     *
     *          其中：
     *          - n_elements: 单元总数
     *          - avg_dof_per_elem: 平均每单元DOF数（如TETRA_SCALAR=4, HEX_VECTOR=12）
     *          - overlap_factor: 重叠因子（0.1~0.3，考虑共享节点的重复计数）
     *
     *          更精确的估算可基于图模型：
     *          NNZ_exact = Σ_{i=1}^{n_free_dofs} (1 + degree(i))
     *          其中degree(i)为第i个自由度的图邻接度（需构建连通图）
     *
     * @note 返回值为上界估计，实际COO→CSR转换后会因重复坐标合并而减少
     * @note 预分配过大不会导致错误（仅浪费少量内存），过小会导致多次扩容
     */
    static int estimateNNZ(
        const std::vector<fe_em::Element>& elements,
        int avg_dof_per_elem
    );

    /**
     * @brief 将单元刚度/质量/阻尼矩阵scatter到全局COO矩阵
     * @param coo_matrix 目标COO矩阵的引用（累积写入）
     * @param elem_matrix 单元级方阵（维度: local_dofs × local_dofs）
     * @param l2g 该单元的局部-全局DOF映射（用于确定scatter目标位置）
     *
     * @details 实现有限元组装的核心scatter操作，
     *          将单元矩阵的每个有效元素累加到全局矩阵的对应位置。
     *
     *          Scatter算法（双重循环遍历单元矩阵）：
     *          ```
     *          for i in [0, local_dofs):
     *              global_i = l2g.indices[i]
     *              if global_i == -1: continue  // 跳过约束DOF（行方向）
     *              for j in [0, local_dofs):
     *                  global_j = l2g.indices[j]
     *                  if global_j == -1: continue  // 跳过约束DOF（列方向）
     *                  coo_matrix.add_value(global_i, global_j, elem_matrix(i,j))
     *          ```
     *
     *          关键特性：
     *          - 自动跳过约束DOF（indices==-1），实现消去法边界条件处理
     *          - COO格式允许重复坐标（同一位置多次add_value将在CSR转换时合并）
     *          - 支持非对称矩阵组装（完整存储上三角和下三角，不利用对称性）
     *
     * @note 对于对称矩阵，可优化为仅存储上三角部分（节省50%内存和计算量）
     * @warning 调用者需确保coo_matrix已预分配足够容量（通过estimateNNZ）
     */
    void scatterElementMatrix(
        CooMatrix<double>& coo_matrix,
        const Eigen::MatrixXd& elem_matrix,
        const fe_em::Local2Global& l2g
    );

    /**
     * @brief 将单元源项向量scatter到全局右端项向量
     * @param global_vector 目标全局向量的引用（累积写入）
     * @param elem_vector 单元级源项向量（维度: local_dofs × 1）
     * @param l2g 该单元的局部-全局DOF映射（用于确定scatter目标位置）
     *
     * @details 实现有限元载荷向量的scatter操作，
     *          将单元源项向量的每个有效分量累加到全局向量的对应位置。
     *
     *          Scatter算法（单循环遍历单元向量）：
     *          ```
     *          for i in [0, local_dofs):
     *              global_i = l2g.indices[i]
     *              if global_i == -1: continue  // 跳过约束DOF
     *              global_vector(global_i) += elem_vector(i)
     *          ```
     *
     *          与scatterElementMatrix的区别：
     *          - 向量scatter仅需单层循环（O(n_dof)复杂度）
     *          - 无需考虑对称性（向量无对称概念）
     *          - 直接写入Eigen::VectorXd（无需经过COO中间格式）
     *
     * @note 全局向量F_应在assemble()开始时初始化为零向量
     * @warning 调用者需确保global_vector已正确初始化且尺寸匹配free_dofs
     */
    void scatterElementVector(
        Eigen::VectorXd& global_vector,
        const Eigen::VectorXd& elem_vector,
        const fe_em::Local2Global& l2g
    );

    /**
     * @brief 检查CSR矩阵是否满足数值对称性
     * @param matrix 待检查的CSR格式矩阵常量引用
     * @return true 如果矩阵在给定容差内对称，false否则
     *
     * @details 验证矩阵A是否满足 A_ij ≈ A_ji （对所有非零元素），
     *          用于判断是否可选择专用的对称矩阵求解器（如Cholesky分解）。
     *
     *          检查算法：
     *          ```
     *          for each nonzero (i, j, A_ij) in CSR matrix:
     *              if i == j: continue  // 对角线元素自动对称
     *              find A_ji at position (j, i) using binary search on col_indices
     *              if A_ji not found or |A_ij - A_ji| > tolerance:
     *                  return false
     *          return true
     *          ```
     *
     *          容差设置：
     *          - 绝对容差：tol_abs = 1e-10 （适用于矩阵元素值在1附近的场景）
     *          - 相对容差：tol_rel = 1e-8  （适用于元素值跨度大的场景）
     *          - 实际使用：|A_ij - A_ji| < tol_abs + tol_rel * max(|A_ij|, |A_ji|)
     *
     *          物理意义验证：
     *          - 刚度矩阵K：自伴算子离散结果应理论上对称（若形函数和积分规则合适）
     *          - 质量矩阵M：质量算子恒对称（正定对称）
     *          - 阻尼矩阵C：若材料各向同性则对称，否则可能非对称
     *
     * @note 时间复杂度: O(nnz × log(avg_row_nnz)) （二分查找转置元素）
     * @note 此方法标记为const，不修改矩阵内容
     * @warning 对于接近奇异的病态矩阵，浮点误差可能导致误判为非对称
     */
    bool checkSymmetry(const CsrMatrix<double>& matrix) const;
};

} // namespace numeric
