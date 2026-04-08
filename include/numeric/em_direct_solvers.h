/**
 * @file em_direct_solvers.h
 * @brief 核心数值层 - 多后端直接求解器实现（工业级鲁棒性）
 * @details 提供3个直接求解器实现，支持Eigen/SuperLU/MUMPS多后端策略模式，
 *          用于电磁场有限元线性方程组 Ax = b 的高性能求解。
 *
 * @par 支持的求解器类型：
 * - SymmetricDirectSolver：对称正定(SPD)矩阵Cholesky分解求解器
 * - SymmetricIndefiniteDirectSolver：对称不定(SID)矩阵LDL^T分解求解器
 * - GeneralDirectSolver：通用非对称矩阵LU分解求解器
 *
 * @par 多后端架构：
 * - Eigen后端（默认）：无需额外依赖，适合中小规模问题（<10万DOF）
 * - SuperLU后端：高性能LU/Cholesky分解，支持列选主元（需编译期启用）
 * - MUMPS后端：MPI并行分布式求解，适合大规模计算（需编译期启用）
 *
 * @par 核心特性：
 * - 瞬态复用优化：set_matrix()一次分解缓存，solve()多次复用（关键性能优化）
 * - 优雅降级机制：请求的后端不可用时自动回退到Eigen并输出WARNING
 * - 工业级错误处理：完整的输入校验、状态报告、错误信息反馈
 * - 对称性检测与正定性验证：自动检测矩阵属性并给出警告
 * - 条件数估计与病态警告：对病态矩阵输出FEEM_WARN级别日志
 *
 * @par 典型使用流程：
 * @code
 * // 创建对称正定求解器（默认使用Eigen后端）
 * numeric::SymmetricDirectSolver spd_solver;
 *
 * // 设置系数矩阵并触发Cholesky分解（结果缓存）
 * spd_solver.set_matrix(stiffness_matrix);
 *
 * // 瞬态分析中多次求解（复用缓存的分解因子，显著加速）
 * for (int t = 0; t < time_steps; t++) {
 *     auto result = spd_solver.solve(rhs_vector[t]);
 *     if (result.status == numeric::SolverStatus::SUCCESS) {
 *         FEEM_INFO("时间步{}求解完成, 残差={:.6e}", t, result.residual_norm);
 *     }
 * }
 *
 * // 动态切换后端（若SuperLU可用则使用，否则优雅降级到Eigen）
 * spd_solver.set_backend(numeric::DirectBackendType::SUPERLU);
 *
 * // 查询可用后端列表
 * auto backends = numeric::DirectBackendManager::getAvailableBackends();
 * @endcode
 *
 * @see EMLinearSolverBase 抽象基类定义
 * @see SparseConverter CSR-Eigen矩阵转换工具
 * @see DirectBackendManager 后端可用性管理器
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

#include "em_linear_solver.h"
#include "em_sparse_converter.h"
#include "em_solver_backends.hpp"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <memory>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>

#if EM_SOLVER_HAS_SUPERLU
// SuperLU头文件已在em_solver_backends.hpp中包含
#endif

#if EM_SOLVER_HAS_MUMPS
// MUMPS头文件已在em_solver_backends.hpp中包含
#endif

namespace numeric {

#if EM_SOLVER_HAS_SUPERLU

/**
 * @class SuperluContext
 * @brief SuperLU_MT 数据结构 RAII 封装，管理完整的 LU 分解生命周期
 * @details 封装 SuperLU_MT 所有的 C 语言数据结构和资源，提供类型安全的 C++ 接口。
 *          严格遵循官方示例（pdgssv.c / pdlinsolx.c / pdspmd.c）的正确调用顺序：
 *
 * @par 生命周期（严格遵循官方流程）：
 * @code
 * // 1. 初始化阶段（setup）
 * ctx.initialize(n, nnz, values, rowind, colptr);
 *
 * // 2. 分解阶段（factorize）—— 仅需一次
 * ctx.factorize(nprocs);
 *
 * // 3. 求解阶段（solve） —— 可重复调用
 * for (each rhs) {
 *     Eigen::VectorXd x = ctx.solve(b);
 * }
 *
 * // 4. 自动析构（RAII），或手动调用 reset()
 * @endcode
 *
 * @par 官方调用顺序对照：
 * - pdgssv.c: get_perm_c → pdgstrf_init → pdgstrf → dgstrs → pdgstrf_finalize → Destroy_*
 * - pdspmd.c: get_perm_c → pdgstrf_init → (threaded) pdgstrf → dgstrs → pxgstrf_finalize → Destroy_*
 * - 本封装严格按照此顺序实现，修复了原始代码中缺少 init/finalize 的致命问题
 *
 * @par 内存安全保证：
 * - 所有 SuperLU_MT 分配的内存通过对应的 Destroy/SUPERLU_FREE 正确释放
 * - L 因子使用 Destroy_SuperNode_SCP（SCP格式专用）
 * - U 因子使用 Destroy_CompCol_NCP（NCP格式专用）
 * - A 矩阵使用 Destroy_CompCol_Matrix（NC格式专用）
 * - 析构函数保证无资源泄漏
 *
 * @warning 此类非线程安全，外部需同步保护 factorize/solve 调用
 * @note 依赖 slu_mt_ddefs.h 中声明的所有 SuperLU_MT API
 */
class SuperluContext {
public:
    /**
     * @brief 默认构造函数，将所有内部状态初始化为未就绪
     */
    SuperluContext();

    /**
     * @brief 析构函数，自动释放所有 SuperLU_MT 资源
     * @details 按照正确的逆序释放所有分配的内存和数据结构，
     *          确保无内存泄漏。等价于调用 reset()。
     */
    ~SuperluContext();

    /** 禁止拷贝（管理唯一资源） */
    SuperluContext(const SuperluContext&) = delete;
    SuperluContext& operator=(const SuperluContext&) = delete;

    /** 允许移动语义（支持 std::unique_ptr 等场景） */
    SuperluContext(SuperluContext&& other) noexcept;
    SuperluContext& operator=(SuperluContext&& other) noexcept;

    /**
     * @brief 初始化 CSC 格式系数矩阵并预分配所有辅助数据结构
     * @param n   矩阵维度（方阵，行数=列数）
     * @param nnz 非零元素个数
     * @param nzval 非零值数组（长度 nnz，所有权转移至本对象）
     * @param rowind 行索引数组（长度 nnz，0-based，所有权转移至本对象）
     * @param colptr 列指针数组（长度 n+1，0-based，所有权转移至本对象）
     * @return true 初始化成功
     * @return false 内存分配失败或参数非法
     *
     * @details 执行以下操作（对应官方示例的初始化部分）：
     * 1. 调用 dCreate_CompCol_Matrix 创建 SuperMatrix A（CSC/NC 格式）
     * 2. 分配 perm_c[] 和 perm_r[] 置换向量
     * 3. 调用 get_perm_c() 计算列置换（最小度排序）
     * 4. 分配 etree[], colcnt_h[], part_super_h[] 辅助数组
     * 5. 配置 superlumt_options_t 默认参数
     *
     * @pre nzval/rowind/colptr 必须由 doubleMalloc/intMalloc 分配（SuperLU_MT 内存管理一致性）
     * @post 内部状态为 INITIALIZED，可调用 factorize()
     * @note 调用后传入的数组指针所有权转移至本对象，不应再由调用方释放
     */
    bool initialize(int_t n, int_t nnz, double* nzval, int_t* rowind, int_t* colptr);

    /**
     * @brief 执行并行 LU 分解（PA = LU）
     * @param nprocs 并行线程数（推荐值：1~硬件线程数）
     * @return 0 成功；>0 矩阵奇异（info=列号）；<0 参数错误或内存不足
     *
     * @details 严格按官方三步流程执行：
     * 1. **pdgstrf_init**: 初始化选项、计算消去树和列计数、创建置换后矩阵 AC
     * 2. **pdgstrf**: 并行 LU 分解，生成 L（SCP超节点格式）和 U（NCP列置换格式）
     * 3. **pdgstrf_finalize**: 清理 AC 矩阵和分解过程中的临时数据结构
     *
     * @pre 必须已成功调用 initialize()
     * @post 内部状态为 FACTORIZED，可重复调用 solve()
     * @warning 重复调用会先自动清理旧的 L/U 因子
     */
    int_t factorize(int_t nprocs = 1);

    /**
     * @brief 使用缓存的 L/U 因子执行三角回代求解
     * @param b 右端项向量（维度必须与 initialize 时的 n 一致）
     * @return SolverResult 包含解向量和状态信息
     *
     * @details 执行流程：
     * 1. 将 Eigen 向量转换为 SuperLU_DenseMatrix 格式
     * 2. 调用 dgstrs() 执行三角系统求解：Ly=Pb, Ux=y
     * 3. 从结果 DenseMatrix 提取解向量到 Eigen::VectorXd
     * 4. 释放临时 DenseMatrix 及其数据数组
     *
     * @pre 必须已成功调用 factorize()
     * @note 可多次调用以复用 L/U 因子（瞬态分析核心优化路径）
     */
    SolverResult solve(const Eigen::VectorXd& b);

    /**
     * @brief 重置所有内部状态并释放全部资源
     * @details 按照与初始化相反的顺序释放所有 SuperLU_MT 数据结构，
     *          使对象回到刚构造后的初始状态。可安全重复调用。
     *
     * @par 释放顺序（严格遵循官方示例的清理模式）：
     * 1. 若已分解：Destroy_SuperNode_SCP(L), Destroy_CompCol_NCP(U)
     * 2. 若有 AC 矩阵：pxgstrf_finalize(options, &AC)
     * 3. 若有 A 矩阵：Destroy_CompCol_Matrix(A)
     * 4. SUPERLU_FREE(perm_r, perm_c)
     * 5. SUPERLU_FREE(etree, colcnt_h, part_super_h)
     * 6. StatFree(Gstat) 若已分配
     */
    void reset();

    /**
     * @brief 检查是否已完成 LU 分解且可进行求解
     * @return true 已分解，solve() 可用
     */
    bool is_factored() const { return state_ == State::FACTORIZED; }

    /**
     * @brief 检查是否已初始化（至少调用了 initialize）
     * @return true 已初始化
     */
    bool is_initialized() const { return state_ == State::INITIALIZED || state_ == State::FACTORIZED; }

    /**
     * @brief 获取矩阵维度
     * @return 矩阵阶数 n（initialize 后有效）
     */
    int_t matrix_size() const { return n_; }

private:
    /**
     * @brief 内部状态枚举，严格控制操作顺序
     */
    enum class State {
        UNINITIALIZED,   ///< 未初始化（初始状态 / reset 后）
        INITIALIZED,     ///< 已初始化系数矩阵和辅助结构（initialize 后）
        FACTORIZED       ///< 已完成 LU 分解（factorize 后）
    };

    State state_ = State::UNINITIALIZED;  ///< 当前生命周期状态

    int_t n_ = 0;       ///< 矩阵维度
    int_t nnz_ = 0;     ///< 非零元数

    // ========== SuperLU_MT 核心数据结构 ==========
    // 对应官方示例中的局部变量，统一封装为成员变量管理生命周期

    std::unique_ptr<SuperMatrix> A_;       ///< 系数矩阵（CSC/NC 格式，dCreate_CompCol_Matrix 创建）
    std::unique_ptr<SuperMatrix> L_;       ///< L 因子（SCP 超节点格式，pdgstrf 输出）
    std::unique_ptr<SuperMatrix> U_;       ///< U 因子（NCP 列置换格式，pdgstrf 输出）
    std::unique_ptr<SuperMatrix> AC_;      ///< 列置换后矩阵（pdgstrf_init 创建，pdgstrf_finalize 销毁）

    std::unique_ptr<int_t[]> perm_c_;      ///< 列置换向量（get_perm_c 计算 + pdgstrf 修改）
    std::unique_ptr<int_t[]> perm_r_;      ///< 行置换向量（pdgstrf 输出）

    std::unique_ptr<int_t[]> etree_;       ///< 消去树（pdgstrf_init 计算，维度 n）
    std::unique_ptr<int_t[]> colcnt_h_;    ///< 列计数（pdgstrf_init 计算，维度 n）
    std::unique_ptr<int_t[]> part_super_h_;///< 超节点划分（pdgstrf_init 计算，维度 n）

    std::unique_ptr<superlumt_options_t> options_;  ///< 求解选项配置（pdgstrf_init 填充）
    std::unique_ptr<Gstat_t> Gstat_;                ///< 统计信息结构（StatAlloc/StatInit/StatFree）

    /**
     * @brief 释放 L 和 U 因子的内存（使用正确的格式专属函数）
     * @details 根据官方示例的清理模式：
     * - L 使用 Destroy_SuperNode_SCP（SCP 超节点列存储格式）
     * - U 使用 Destroy_CompCol_NCP（NCP 列置换压缩格式）
     * 注意：不能用通用的 Destroy_SuperMatrix_Store，会导致 NCP 格式内存泄漏
     */
    void release_factors();

    /**
     * @brief 释放 AC 置换矩阵（若 pdgstrf_init 已创建）
     * @details 调用 pxgstrf_finalize 或 pdgstrf_finalize 清理 AC 相关内部结构
     */
    void release_ac_matrix();
};

#endif  // EM_SOLVER_HAS_SUPERLU

/**
 * @class DirectBackendManager
 * @brief 直接求解器后端可用性管理器
 * @details 提供编译期和运行时的后端可用性查询功能，基于CMake宏检测第三方库是否已集成。
 *          所有方法均为静态方法，无需实例化即可调用。
 *
 * @note 设计目标：
 *       - 为上层提供统一的后端查询接口
 *       - 支持运行时动态检查后端可用性
 *       - 提供友好的后端名称字符串用于日志输出
 *
 * @par 编译期配置依赖：
 * - HAVE_SUPERLU：CMake编译选项，启用SuperLU后端支持
 * - HAVE_MUMPS：CMake编译选项，启用MUMPS后端支持
 * - Eigen后端始终可用（作为默认回退方案）
 */
class DirectBackendManager {
public:
    /**
     * @brief 检查指定后端类型在当前编译配置下是否可用
     * @param type 待检查的后端类型枚举值
     * @return true 后端可用（已编译进二进制文件）
     * @return false 后端不可用（未编译或未安装对应库）
     *
     * @details 基于编译期宏EM_SOLVER_HAS_SUPERLU/EM_SOLVER_HAS_MUMPS进行判断，
     *          Eigen后端始终返回true。
     */
    static bool isBackendAvailable(DirectBackendType type);

    /**
     * @brief 获取后端类型的友好名称字符串
     * @param type 后端类型枚举值
     * @return std::string 后端名称（如"Eigen"、"SuperLU"、"MUMPS"）
     *
     * @note 用于日志输出、调试信息和用户界面显示
     */
    static std::string getBackendName(DirectBackendType type);

    /**
     * @brief 获取当前编译配置下所有可用的后端类型列表
     * @return std::vector<DirectBackendType> 可用后端向量（至少包含EIGEN）
     *
     * @details 返回顺序：EIGEN始终在首位，随后是SUPERLU和MUMPS（若可用）
     */
    static std::vector<DirectBackendType> getAvailableBackends();

private:
    DirectBackendManager() = delete;  // 禁止实例化（纯静态工具类）
};

/**
 * @class SymmetricDirectSolver
 * @brief 对称正定(SPD)矩阵直接求解器 - Cholesky分解 (A = LL^T)
 * @details 继承EMLinearSolverBase抽象基类，专为电磁场有限元分析中的对称正定系统设计，
 *          如静磁场矢量磁位A-φ formulation的刚度矩阵求解、静电场标量电位Poisson方程等。
 *
 * @par 算法原理：
 * - 对于SPD矩阵A，执行Cholesky分解 A = L * L^T（L为下三角矩阵）
 * - 求解过程分两步：前代 Ly = b，回代 L^T x = y
 * - 时间复杂度：分解O(n³)，求解O(n²)（n为矩阵维度）
 * - 数值稳定性：对于良态SPD矩阵，Cholesky分解无条件稳定
 *
 * @par 适用场景：
 * - 静磁场分析（μ恒定介质，无涡流效应）
 * - 静电场分析（ε恒定介质的Laplace/Poisson方程）
 * - 结构力学刚度矩阵求解（对称正定性保证）
 * - 热传导稳态分析
 *
 * @par 性能优化特性：
 * - 瞬态复用：set_matrix()执行一次符号+数值分解，后续solve()仅做前代回代
 * - 对称性利用：仅需存储一半矩阵元素，内存占用减半
 * - 填入优化：SimplicialLLT使用近似最小度(AMD)排序减少填入
 *
 * @par 多后端支持：
 * - Eigen::SimplicialLLT（默认）：单线程，适合中小规模（推荐<5万DOF）
 * - SuperLU（可选）：支持列选主元的Cholesky变体，数值更稳健
 * - MUMPS（可选）：MPI并行，支持分布式内存计算（推荐>10万DOF）
 *
 * @code
 * // 示例：静磁场SPD系统求解
 * numeric::SymmetricDirectSolver solver;
 * solver.set_symmetry_tolerance(1e-8);  // 设置对称性校验容差
 *
 * // 设置刚度矩阵K（触发Cholesky分解并缓存L因子）
 * solver.set_matrix(K_spd);
 *
 * // 瞬态循环中多次求解（复用L因子，每次仅需O(n²)前代回代）
 * for (int step = 0; step < n_steps; step++) {
 *     auto result = solver.solve(F[step]);  // F为载荷向量
 *     if (result.status != numeric::SolverStatus::SUCCESS) {
 *         FEEM_ERROR("求解失败: {}", result.error_msg);
 *         break;
 *     }
 *     // 使用result.x更新场量...
 * }
 *
 * // 切换到MUMPS并行后端（若可用）
 * if (numeric::DirectBackendManager::isBackendAvailable(numeric::DirectBackendType::MUMPS)) {
 *     solver.set_backend(numeric::DirectBackendType::MUMPS);
 * }
 * @endcode
 *
 * @warning 输入矩阵必须为对称正定(SP)D，否则分解可能失败或产生无意义结果
 * @note 线程安全：set_matrix()/solve()需外部同步保护
 */
class SymmetricDirectSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数，初始化求解器并指定初始后端类型
     * @param backend 初始后端类型（默认EIGEN），若不可用则自动降级到EIGEN
     *
     * @details 构造时执行后端可用性检查，若请求的后端未编译则自动降级到Eigen
     *          并输出FEEM_WARN级别日志提示用户。
     */
    explicit SymmetricDirectSolver(DirectBackendType backend = DirectBackendType::EIGEN);

    /**
     * @brief 析构函数，自动清理所有内部资源
     * @details 调用clear()释放矩阵分解结果和临时工作空间
     */
    ~SymmetricDirectSolver() override = default;

    /**
     * @brief 动态切换求解后端类型
     * @param backend 目标后端类型（EIGEN/SUPERLU/MUMPS）
     *
     * @details 执行后端切换操作：
     * - 检查目标后端是否可用（编译期宏检测）
     * - 若不可用，自动降级到Eigen并输出WARNING日志
     * - 若当前已有缓存的矩阵分解结果，清除旧缓存
     * - 切换后需重新调用set_matrix()才能进行求解
     *
     * @note 优雅降级策略确保程序不会因缺少第三方库而崩溃
     * @warning 切换后端会丢失当前的矩阵分解缓存
     */
    void set_backend(DirectBackendType backend);

    /**
     * @brief 设置对称性校验容差阈值
     * @param tol 容差值（默认1e-10），用于判断矩阵是否"足够"对称
     *
     * @details 在set_matrix()中会检查 ||A - A^T||_F / ||A||_F < tol，
     *          若超过容差则输出WARNING但不阻止分解（某些浮点误差可接受）。
     *
     * @note 推荐范围：1e-8 ~ 1e-12（过松可能掩盖错误，过严可能导致误报）
     */
    void set_symmetry_tolerance(double tol);

    /**
     * @brief 设置系数矩阵并触发Cholesky分解（实现基类纯虚函数）
     * @param A CSR格式的稀疏系数矩阵（必须为方阵且对称正定）
     *
     * @details 执行完整的前置处理和分解流程：
     * 1. 输入合法性校验（维度、方阵检查）
     * 2. CSR → Eigen格式转换（通过SparseConverter工具类）
     * 3. 可选的对称性校验（基于symmetry_tol_阈值）
     * 4. 调用当前后端的分解接口（符号分析 + 数值分解）
     * 5. 缓存分解结果用于后续多次solve()复用
     * 6. 记录分解耗时和统计信息到日志
     *
     * @pre A必须是已构建完成的CSR矩阵（is_built() == true）
     * @post 内部缓存L因子，matrix_set_标志置为true
     * @note 此方法是瞬态优化的核心：一次分解，多次求解复用
     * @exception 若矩阵非法或分解失败，设置内部错误状态并记录ERROR日志
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 求解线性系统 Ax = b（实现基类纯虚函数）
     * @param b 右端项向量（载荷向量），维度必须与矩阵行数一致
     * @return SolverResult 求解结果结构体（含解向量和状态信息）
     *
     * @details 利用缓存的Cholesky分解因子执行高效求解：
     * 1. 输入校验（矩阵已设置、维度匹配）
     * 2. 调用当前后端的solve接口（前代 Ly = b + 回代 L^T x = y）
     * 3. 计算残差范数 ||Ax - b||_2 用于质量评估
     * 4. 记录求解耗时和性能指标
     *
     * @pre 必须先成功调用set_matrix()
     * @post 返回的result.x仅在status == SUCCESS时有效
     * @note 时间复杂度O(n²)（远低于重新分解的O(n³)）
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    /**
     * @brief 获取求解器名称标识（实现基类纯虚函数）
     * @return std::string 格式："SymmetricDirect_[后端名]"（如"SymmetricDirect_Eigen"）
     */
    std::string get_solver_name() const override;

    /**
     * @brief 清理所有内部资源并重置状态（实现基类纯虚函数）
     * @details 释放以下资源：
     * - Eigen求解器对象（L因子存储）
     * - SuperLU数据结构（若已分配）
     * - MUMPS数据结构（若已分配）
     * - Eigen格式矩阵缓存
     * - 重置matrix_set_标志为false
     *
     * @post 求解器恢复到刚构造后的初始状态
     * @note 可手动调用以提前释放内存，或在析构时自动调用
     */
    void clear() override;

private:
    DirectBackendType backend_type_;   ///< 当前激活的后端类型
    double symmetry_tol_ = 1e-10;      ///< 对称性校验容差阈值
    bool matrix_set_ = false;          ///< 矩阵是否已设置并成功分解

    Eigen::SparseMatrix<double> eigen_matrix_;  ///< 缓存的Eigen格式系数矩阵

    // Eigen后端实例（SimplicialLLT: 稀疏Cholesky分解 LL^T）
    std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> eigen_solver_;

#if EM_SOLVER_HAS_SUPERLU
    std::unique_ptr<SuperluContext> superlu_ctx_;  ///< SuperLU_MT RAII 封装（管理完整生命周期）
#endif

#if EM_SOLVER_HAS_MUMPS
    // MUMPS后端数据结构
    struct DMUMPS_STRUC_C mumps_data_;                   ///< MUMPS工作数据结构
    bool mumps_initialized_ = false;                      ///< MUMPS是否已初始化
#endif

    /**
     * @brief 将请求的后端降级到实际可用的后端（优雅降级策略）
     * @param requested 用户请求的后端类型
     * @return DirectBackendType 实际使用的后端类型（若requested不可用则返回EIGEN）
     *
     * @details 若requested后端未编译进二进制文件，输出WARNING日志并返回EIGEN
     */
    DirectBackendType fallback_to_eigen_if_unavailable(DirectBackendType requested);

    /**
     * @brief 使用Eigen后端执行Cholesky分解
     * @return SolverResult 分解结果（SUCCESS或NUMERICAL_ERROR）
     */
    SolverResult decompose_with_eigen();

    /**
     * @brief 使用Eigen后端执行求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve_with_eigen(const Eigen::VectorXd& b);

#if EM_SOLVER_HAS_SUPERLU
    /**
     * @brief 使用SuperLU后端执行分解（对称优化模式）
     * @return SolverResult 分解结果
     */
    SolverResult decompose_with_superlu();

    /**
     * @brief 使用SuperLU后端执行求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve_with_superlu(const Eigen::VectorXd& b);
#endif

#if EM_SOLVER_HAS_MUMPS
    /**
     * @brief 使用MUMPS后端执行分解（对称正定并行模式）
     * @return SolverResult 分解结果
     */
    SolverResult decompose_with_mumps();

    /**
     * @brief 使用MUMPS后端执行求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve_with_mumps(const Eigen::VectorXd& b);
#endif

    /**
     * @brief 校验矩阵对称性（辅助方法）
     * @param mat 待检验的Eigen稀疏矩阵
     * @param tol 容差阈值
     * @return true 矩阵在容差范围内对称
     * @return false 矩阵不对称性超过阈值
     *
     * @details 计算 ||mat - mat^T||_F / ||mat||_F 并与tol比较
     */
    bool check_symmetry(const Eigen::SparseMatrix<double>& mat, double tol) const;
};

/**
 * @class SymmetricIndefiniteDirectSolver
 * @brief 对称不定(SID)矩阵直接求解器 - LDL^T分解（带正则化）
 * @details 继承EMLinearSolverBase，专为对称但未必正定的线性系统设计，
 *          如Saddle Point问题、约束优化KKT系统、混合有限元离散化等。
 *
 * @par 算法原理：
 * - 执行LDL^T分解：A = L * D * L^T（L为单位下三角，D为块对角1x1/2x2）
 * - 相比Cholesky分解，允许D的对角元出现负值或零（处理不定矩阵）
 * - 通过正则化技术处理奇异矩阵（添加微小扰动 εI 保证可逆）
 *
 * @par 适用场景：
 * - 混合有限元公式（速度-压力耦合，如Stokes方程）
 * - Saddle Point系统（带约束的优化问题）
 * - 半正定矩阵（存在零特征值，如悬空节点）
 * - 对称但不保证正定的一般情况
 *
 * @par 正则化机制：
 * - 当检测到接近奇异的矩阵时，自动添加 regularization_epsilon_ * I 扰动
 * - 正则化参数可通过 set_regularization_epsilon() 配置（默认1e-12）
 * - 正则化前后输出INFO级别日志说明扰动幅度
 *
 * @par 零空间检测：
 * - 求解完成后分析D因子的对角元，统计接近零的主元个数
 * - 报告零空间维度（nullity）用于诊断矩阵秩亏
 * - 当零空间维度>0时输出WARNING日志提示可能的物理建模问题
 *
 * @code
 * // 示例：混合有限元Saddle Point系统求解
 * numeric::SymmetricIndefiniteDirectSolver sid_solver;
 * sid_solver.set_regularization_epsilon(1e-14);  // 设置微小正则化参数
 *
 * // 设置Saddle Point矩阵 [K  B^T; B  0]（对称但不定）
 * sid_solver.set_matrix(saddle_point_matrix);
 *
 * auto result = sid_solver.solve(rhs_vector);
 * if (result.status == numeric::SolverStatus::SUCCESS) {
 *     FEEM_INFO("SID求解成功, 解向量范数={:.6e}", result.x.norm());
 * } else if (result.status == numeric::SolverStatus::NUMERICAL_ERROR) {
 *     FEEM_ERROR("矩阵可能奇异: {}", result.error_msg);
 * }
 * @endcode
 *
 * @note Eigen::SimplicialLDLT自动处理2x2主元旋转，数值稳定性优于朴素LDL^T
 * @warning 对于强不定矩阵，建议优先测试GeneralDirectSolver（LU分解更稳定）
 */
class SymmetricIndefiniteDirectSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数，初始化SID求解器
     * @param backend 初始后端类型（目前仅支持EIGEN，其他后端预留接口）
     */
    explicit SymmetricIndefiniteDirectSolver(DirectBackendType backend = DirectBackendType::EIGEN);

    ~SymmetricIndefiniteDirectSolver() override = default;

    /**
     * @brief 设置正则化参数ε（用于处理奇异/近奇异矩阵）
     * @param epsilon 正则化强度（默认1e-12），添加到矩阵对角线：A_reg = A + εI
     *
     * @details 当矩阵接近奇异时（最小特征值接近0），添加微小扰动可保证可逆性。
     *          参数选择原则：
     *          - 过小（<1e-15）：无法有效消除奇异性
     *          - 过大（>1e-6）：显著改变原问题的物理意义
     *          - 推荐范围：1e-14 ~ 1e-10（根据问题规模调整）
     *
     * @note 正则化会影响解的精度，应在精度和稳定性间权衡
     */
    void set_regularization_epsilon(double epsilon);

    /**
     * @brief 设置系数矩阵并触发LDL^T分解（实现基类纯虚函数）
     * @param A CSR格式稀疏系数矩阵（必须为方阵且对称）
     *
     * @details 执行流程：
     * 1. 维度和对称性校验
     * 2. CSR→Eigen格式转换
     * 3. 尝试标准LDL^T分解
     * 4. 若分解失败（检测到数值零主元），自动应用正则化并重试
     * 5. 缓存分解结果
     *
     * @post matrix_set_标志置为true（即使使用了正则化）
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 求解线性系统 Ax = b（实现基类纯虚函数）
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    std::string get_solver_name() const override;
    void clear() override;

private:
    DirectBackendType backend_type_;
    double regularization_epsilon_ = 1e-12;  ///< 正则化参数（默认1e-12）
    double symmetry_tol_ = 1e-10;            ///< 对称性校验容差
    bool matrix_set_ = false;
    bool regularization_applied_ = false;    ///< 是否已施加正则化

    Eigen::SparseMatrix<double> eigen_matrix_;
    Eigen::SparseMatrix<double> regularized_matrix_;  ///< 正则化后的矩阵备份

    // Eigen LDL^T求解器（支持2x2 pivot，处理不定矩阵）
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> eigen_solver_;

    DirectBackendType fallback_to_eigen_if_unavailable(DirectBackendType requested);
    SolverResult decompose_with_eigen();
    SolverResult solve_with_eigen(const Eigen::VectorXd& b);
    bool check_symmetry(const Eigen::SparseMatrix<double>& mat, double tol) const;
};

/**
 * @class GeneralDirectSolver
 * @brief 通用非对称矩阵直接求解器 - LU分解（带部分选主元）
 * @details 继承EMLinearSolverBase，适用于一般方形矩阵（不要求对称或正定），
 *          是最通用的直接求解器实现，覆盖所有可直接求解的线性系统。
 *
 * @par 算法原理：
 * - 执行LU分解（带部分选主元PA = LU）：P为置换矩阵，L为下三角，U为上三角
 * - 求解过程：前代 Ly = Pb，回代 Ux = y
 * - 选主元策略：列选主元（partial pivoting）保证数值稳定性
 * - 时间复杂度：分解O(n³/3)，求解O(n²)
 *
 * @par 适用场景：
 * - 非对称矩阵系统（如含对流项的流体力学）
 * - 复数矩阵系统（经过实部/虚部拆分后）
 * - 无法保证对称性的通用稀疏系统
 * - 作为其他专用求解器的最终回退方案
 *
 * @par 病态矩阵处理：
 * - 自动估算条件数（通过U因子的对角元比值）
 * - 当条件数 > 1e10 时输出FEEM_WARN警告（解可能不准确）
 * - 当条件数 > 1e15 时返回NUMERICAL_ERROR（矩阵几乎奇异）
 * - 可通过 set_pivot_threshold() 调整选主元策略平衡稳定性和填充
 *
 * @par 性能考虑：
 * - LU分解的fill-in通常高于Cholesky/LDL^T（不对称性导致更多非零元）
 * - 对于大型非对称系统，建议考虑迭代法（如GMRES + ILU预条件子）
 * - Eigen::SparseLU使用AMD + COLAMD双重排序优化fill-in
 *
 * @code
 * // 示例：通用非对称系统求解
 * numeric::GeneralDirectSolver general_solver;
 * general_solver.set_pivot_threshold(1.0);  // 完全选主元（最稳定但较慢）
 *
 * general_solver.set_matrix(general_matrix);  // 触发LU分解
 *
 * auto result = general_solver.solve(rhs);
 * if (result.status == numeric::SolverStatus::SUCCESS) {
 *     if (result.residual_norm > 1e-6) {
 *         FEEM_WARN("残差较大({:.2e})，矩阵可能病态", result.residual_norm);
 *     }
 *     // 使用result.x...
 * }
 * @endcode
 *
 * @note 对于已知对称的正定矩阵，优先使用SymmetricDirectSolver（更快更稳定）
 * @warning 大型非对称系统的LU分解可能消耗大量内存（fill-in效应）
 */
class GeneralDirectSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数，初始化通用LU求解器
     * @param backend 初始后端类型（目前仅支持EIGEN）
     */
    explicit GeneralDirectSolver(DirectBackendType backend = DirectBackendType::EIGEN);

    ~GeneralDirectSolver() override = default;

    /**
     * @brief 设置LU分解的选主元阈值
     * @param threshold 选主元阈值（范围[0, 1]，默认1.0）
     *
     * @details 控制列选主元的激进程度：
     * - threshold=1.0：完全选主元（最稳定，最大fill-in）
     * - threshold=0.5：中等选主元（平衡稳定性和效率）
     * - threshold=0.0：不选主元（最快，但数值不稳定）
     *
     * @note 推荐保持默认值1.0以确保数值稳定性
     * @warning 降低阈值可能导致病态矩阵求解失败
     */
    void set_pivot_threshold(double threshold);

    /**
     * @brief 设置系数矩阵并触发LU分解（实现基类纯虚函数）
     * @param A CSR格式稀疏系数矩阵（必须为方阵）
     *
     * @details 执行流程：
     * 1. 方阵维度校验（不要求对称）
     * 2. CSR→Eigen格式转换
     * 3. 执行带选主元的LU分解 PA = LU
     * 4. 估算条件数并输出病态警告（若适用）
     * 5. 缓存P、L、U因子用于后续求解
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 求解线性系统 Ax = b（实现基类纯虚函数）
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    std::string get_solver_name() const override;
    void clear() override;

private:
    DirectBackendType backend_type_;
    double pivot_threshold_ = 1.0;      ///< 选主元阈值（默认1.0完全选主元）
    bool matrix_set_ = false;

    Eigen::SparseMatrix<double> eigen_matrix_;

    // Eigen LU求解器（带部分选主元）
    std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> eigen_solver_;

    DirectBackendType fallback_to_eigen_if_unavailable(DirectBackendType requested);
    SolverResult decompose_with_eigen();
    SolverResult solve_with_eigen(const Eigen::VectorXd& b);

    /**
     * @brief 估算矩阵条件数（基于LU分解的U因子）
     * @return double 条件数估计值 κ(A)
     *
     * @details 使用U因子对角元的最大/最小绝对值比作为条件数的粗略估计：
     *          κ(A) ≈ max(|U_ii|) / min(|U_ii|)
     *          这是一个下界估计，真实条件数可能更大。
     */
    double estimate_condition_number() const;
};

} // namespace numeric
