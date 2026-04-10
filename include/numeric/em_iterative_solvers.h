/**
 * @file em_iterative_solvers.h
 * @brief 核心数值层 - 高性能迭代求解器与AMG预条件子框架
 * @details 提供面向大规模电磁场问题（百万/千万自由度）的迭代求解系统，
 *          包括共轭梯度法(CG)、稳定双共轭梯度法(BiCGSTAB)和代数多重网格(AMG)预条件子。
 *          设计目标：
 *          - 支持对称正定/半正定矩阵（标量场、棱边单元静磁场）
 *          - 支持非对称矩阵（A-V混合涡流场、谐波场）
 *          - 提供多种预条件子选项（Jacobi、ILU0、AMG）
 *          - 完整的收敛监控和发散检测机制
 *
 * @par 系统架构：
 * @code
 * EMLinearSolverBase (抽象基类)
 *     ├── CGSolver (共轭梯度求解器)
 *     │   └── 预条件子: NONE / JACOBI / ILU0
 *     ├── BiCGSTABSolver (稳定双共轭梯度求解器)
 *     │   └── 预条件子: NONE / JACOBI / ILU0
 *     └── ScalarAMG (代数多重网格预条件子)
 *         └── 多层V/W循环 + Jacobi/GS光滑器
 *
 * HiptmairAMG (H(curl)专用AMG接口预留)
 * @endcode
 *
 * @par 性能特征：
 * - CGSolver: O(n) 每次迭代，适合SPD矩阵
 * - BiCGSTABSolver: O(n) 每次迭代，适合一般稀疏矩阵
 * - ScalarAMG: O(n) 预计算，O(n) 每次应用，接近最优收敛速度
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

#include "em_linear_solver.h"
#include "em_sparse_converter.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <memory>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace numeric {

// ==================== 公共配置和工具类 ====================

/**
 * @struct IterativeSolverConfig
 * @brief 迭代求解器统一配置参数结构体
 * @details 定义CG/BiCGSTAB等Krylov子空间方法的通用配置参数，
 *          包括收敛准则、最大迭代次数、预条件子类型选择等。
 */
struct IterativeSolverConfig {
    double tolerance = 1e-8;       ///< 相对残差收敛阈值：||r||/||b|| < tolerance
    int max_iterations = 1000;      ///< 最大允许迭代次数
    int residual_monitor_interval = 10; ///< 残差记录间隔（每N步记录一次）

    /**
     * @enum PreconditionerType
     * @brief 预条件子类型枚举
     */
    enum class PreconditionerType {
        NONE,       ///< 无预条件子（纯迭代）
        JACOBI,     ///< Jacobi预条件子（对角缩放）：M^{-1} = D^{-1}
        ILU0        ///< ILU(0)预条件子（不完全LU分解，零填充策略）
    };

    PreconditionerType preconditioner = PreconditionerType::NONE; ///< 预条件子类型选择
};

/**
 * @class ConvergenceHistory
 * @brief 迭代求解器收敛历史记录器
 * @details 记录每次迭代的残差范数，用于收敛曲线绘制、性能分析和调试诊断。
 *          支持查询总迭代次数、最终残差、是否已收敛等信息。
 *
 * @note 内存开销：O(iterations)，对于max_iterations=1000的情况可忽略不计
 * @note 线程安全：此类的实例不应在多线程间共享使用
 */
class ConvergenceHistory {
public:
    /**
     * @brief 构造函数
     */
    ConvergenceHistory() = default;

    /**
     * @brief 记录单次迭代信息
     * @param iter 当前迭代次数（从1开始计数）
     * @param residual_norm 当前相对残差范数 ||r||/||b||
     */
    void record_iteration(int iter, double residual_norm);

    /**
     * @brief 获取已记录的迭代次数序列
     * @return std::vector<int> 迭代次数列表 [1, 2, 3, ..., n]
     */
    std::vector<int> get_iterations() const;

    /**
     * @brief 获取对应的残差范数序列
     * @return std::vector<double> 残差范数列表 [r_1, r_2, r_3, ..., r_n]
     */
    std::vector<double> get_residuals() const;

    /**
     * @brief 判断是否已达到收敛标准
     * @param tolerance 收敛阈值
     * @return bool 最后一次记录的残差<tolerance返回true
     */
    bool has_converged(double tolerance) const;

    /**
     * @brief 获取总迭代次数
     * @return int 已记录的迭代轮数
     */
    int total_iterations() const;

    /**
     * @brief 获取最终残差范数
     * @return double 最后一次记录的残差值，无记录时返回-1.0
     */
    double final_residual() const;

    /**
     * @brief 清除所有历史记录
     */
    void clear();

private:
    std::vector<int> iterations_;      ///< 迭代次数序列
    std::vector<double> residuals_;    ///< 对应残差范数序列
};

// ==================== CGSolver 共轭梯度求解器 ====================

/**
 * @struct CGConfig
 * @brief CG求解器专用配置参数
 * @details 继承IterativeSolverConfig基础配置，添加CG特定参数。
 *          适用场景：对称正定(SPD)/半正定矩阵，如：
 *          - 标量电位有限元方程
 *          - 棱边单元静磁场（磁矢势A formulation）
 *          - 热传导、扩散方程离散系统
 */
struct CGConfig : public IterativeSolverConfig {
    // CG求解器无需额外特殊参数，继承基础配置即可
    using IterativeSolverConfig::PreconditionerType;
};

/**
 * @class CGSolver
 * @brief 共轭梯度法（Conjugate Gradient）迭代求解器
 * @details 基于Eigen::ConjugateGradient实现的SPD矩阵高效求解器。
 *
 * @par 算法原理：
 * CG方法通过在Krylov子空间 K_k(A, r_0) = span{r_0, Ar_0, ..., A^{k-1}r_0}
 * 中搜索最优解来求解线性系统 Ax = b。关键特性：
 * - 对于n×n SPD矩阵，理论上最多n步精确收敛（忽略舍入误差）
 * - 每次迭代仅需一次矩阵-向量乘法和少量向量操作
 * - 残差单调递减，保证稳定性
 *
 * @par 时间复杂度：
 * - 单次迭代: O(nnz)，nnz为矩阵非零元素数
 * - 总复杂度: O(k * nnz)，k为实际迭代次数
 * - 典型收敛: k << n（良好条件下对数级）
 *
 * @par 空间复杂度：
 * O(n) 用于存储工作向量（p, q, r, x等）
 *
 * @par 适用场景：
 * - 对称正定矩阵（刚度矩阵、质量矩阵）
 * - 条件数不太大的良态系统（κ(A) < 10^6）
 * - 大规模稀疏系统（百万DOF级别）
 *
 * @par 不适用场景：
 * - 非对称或不定矩阵 → 使用BiCGSTABSolver
 * - 极病态系统（κ > 10^8）→ 需要强预条件子（如AMG）
 *
 * @code
 * // 使用示例：求解静电场Poisson方程
 * numeric::CGConfig config;
 * config.tolerance = 1e-10;
 * config.max_iterations = 2000;
 * config.preconditioner = CGConfig::PreconditionerType::JACOBI;
 *
 * numeric::CGSolver solver(config);
 * solver.set_matrix(stiffness_matrix);
 *
 * auto result = solver.solve(rhs_vector);
 * if (result.status == numeric::SolverStatus::SUCCESS) {
 *     FEEM_INFO("CG求解成功, 迭代{}次, 残差={:.2e}", 
 *               result.iterations, result.residual_norm);
 * }
 * @endcode
 *
 * @see BiCGSTABSolver ScalarAMG EMLinearSolverBase
 */
class CGSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数
     * @param config CG求解器配置参数（可选，使用默认值）
     */
    explicit CGSolver(const CGConfig& config = CGConfig{});

    /**
     * @brief 虚析构函数
     */
    ~CGSolver() override = default;

    /**
     * @brief 更新求解器配置
     * @param config 新的配置参数
     * @note 调用后需重新调用set_matrix()以重新设置预条件子
     */
    void set_config(const CGConfig& config);

    /**
     * @brief 获取当前配置
     * @return CGConfig 当前使用的配置参数副本
     */
    CGConfig get_config() const;

    /**
     * @brief 获取最近一次求解的收敛历史
     * @return ConvergenceHistory 收敛历史对象（包含迭代次数和残差序列）
     */
    ConvergenceHistory get_convergence_history() const;

    /**
     * @brief 设置系数矩阵并预计算预条件子
     * @param A CSR格式系数矩阵（必须是对称正定的方阵）
     *
     * @details 内部执行：
     * 1. 将CsrMatrix转换为Eigen::SparseMatrix
     * 2. 若选择了Jacobi预条件子，提取对角线元素 D^{-1}
     * 3. 若选择了ILU0预条件子，执行不完全LU分解
     *
     * @note 此方法会触发预条件子的预计算，可能耗时较长（特别是ILU0）
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 执行CG迭代求解 Ax = b
     * @param b 右端项向量
     * @return SolverResult 包含解向量、状态、迭代次数、残差、耗时等完整信息
     *
     * @details 求解流程：
     * 1. 输入校验（维度匹配、矩阵已设置）
     * 2. 根据配置选择预条件子分支
     * 3. 执行CG迭代直到收敛或达到最大迭代次数
     * 4. 记录收敛历史、统计性能指标
     * 5. 返回完整结果
     *
     * @par 收敛判断：
     * relative_residual = ||Ax - b||_2 / ||b||_2 < tolerance
     *
     * @par 发散检测：
     * 连续多次残差增长超过1000倍时提前终止并标记DIVERGED
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    /**
     * @brief 设置复数系数矩阵并配置预条件子（复数版本）
     * @param A CSR格式复数稀疏系数矩阵（必须是对称正定的方阵）
     *
     * @details 内部执行：
     * 1. 将CsrMatrix<std::complex<double>>转换为Eigen::SparseMatrix<std::complex<double>>
     * 2. 创建或重置复数Eigen ConjugateGradient求解器实例
     * 3. 配置预条件子（使用IncompleteCholesky预条件子，与实数版本一致）
     * 4. 设置容差和最大迭代次数
     * 5. 调用eigen_solver_complex_->compute()分析矩阵
     *
     * @note 此方法会触发预条件子的预计算，可能耗时较长
     * @note 复数版本与实数版本维护独立的内部状态，可独立使用
     */
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;

    /**
     * @brief 执行CG迭代求解复数线性系统 Ax = b（复数版本）
     * @param b 复数右端项向量
     * @return SolverResult 求解结果（复数解存储在 x_complex 字段）
     *
     * @details 求解流程：
     * 1. 输入校验（检查matrix_set_complex_标志和维度匹配）
     * 2. 调用eigen_solver_complex_->solve(b)执行迭代求解
     * 3. 将结果存储到result.x_complex
     * 4. 记录迭代次数和误差信息
     * 5. 判断收敛状态并返回完整结果
     *
     * @par 收敛判断：
     * relative_residual = ||Ax - b||_2 / ||b||_2 < tolerance
     *
     * @note 必须先调用set_matrix(const CsrMatrix<std::complex<double>>&)设置复数矩阵
     */
    SolverResult solve(const Eigen::VectorXcd& b) override;

    /**
     * @brief 获取求解器名称标识
     * @return std::string "CGSolver_[预条件子类型]"
     */
    std::string get_solver_name() const override;

    /**
     * @brief 清理所有资源并重置状态
     * @details 释放矩阵存储、预条件子数据、清空收敛历史
     */
    void clear() override;

private:
    CGConfig config_;                                    ///< 配置参数
    Eigen::SparseMatrix<double> eigen_matrix_;           ///< Eigen格式的系数矩阵
    bool matrix_set_ = false;                            ///< 矩阵是否已设置标志
    ConvergenceHistory history_;                         ///< 收敛历史记录器

    Eigen::VectorXd jacobi_diag_;                        ///< Jacobi预条件子对角线 D^{-1}
    Eigen::IncompleteLUT<double> ilu0_precond_;          ///< ILU(0)预条件子

    // 复数版本内部状态
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;  ///< 缓存的复数Eigen格式系数矩阵
    bool matrix_set_complex_ = false;                    ///< 复数矩阵是否已设置

    /// 复数Eigen ConjugateGradient求解器实例（使用IncompleteCholesky预条件子）
    std::unique_ptr<Eigen::ConjugateGradient<
        Eigen::SparseMatrix<std::complex<double>>,
        Eigen::Lower|Eigen::Upper,
        Eigen::IncompleteCholesky<std::complex<double>>>> eigen_solver_complex_;

    /**
     * @brief 设置预条件子（根据config_.preconditioner类型）
     */
    void setup_preconditioner();

    /**
     * @brief 无预条件子的纯CG求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     */
    SolverResult solve_no_preconditioner(const Eigen::VectorXd& b);

    /**
     * @brief Jacobi预条件子的CG求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     * @note 预条件操作: z = M^{-1} * r = D^{-1} * r （逐元素除法）
     */
    SolverResult solve_with_jacobi(const Eigen::VectorXd& b);

    /**
     * @brief ILU0预条件子的CG求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     * @note 使用Eigen::IncompleteLUT实现ILU分解
     */
    SolverResult solve_with_ilu0(const Eigen::VectorXd& b);
};

// ==================== BiCGSTABSolver 双共轭梯度求解器 ====================

/**
 * @struct BiCGSTABConfig
 * @brief BiCGSTAB求解器专用配置参数
 * @details 继承IterativeSolverConfig基础配置。适用场景：
 * - 非对称矩阵（对流扩散、涡流场）
 * - 复对称矩阵（谐波电磁场）
 * - 一般稀疏线性系统
 */
struct BiCGSTABConfig : public IterativeSolverConfig {
    using IterativeSolverConfig::PreconditionerType;
};

/**
 * @class BiCGSTABSolver
 * @brief 稳定双共轭梯度法（Bi-Conjugate Gradient Stabilized）迭代求解器
 * @details 基于Eigen::BiCGSTAB实现的一般稀疏矩阵求解器，专门优化了标准BiCG的不稳定性问题。
 *
 * @par 算法原理：
 * BiCGSTAB是BiCG方法的改进版本，通过引入额外的稳定化步骤避免标准BiCG的振荡行为：
 * - 同时在Krylov子空间 K_k(A, r_0) 和 K_k(A^T, r*_0) 中搜索
 * - 每两步进行一次GMRES(1)风格的极小化，平滑残差曲线
 * - 对非对称矩阵具有更稳定的收敛特性
 *
 * @par 时间复杂度：
 * - 单次迭代: O(2*nnz)（两次矩阵-向量乘法）
 * - 总复杂度: O(k * nnz)，k通常比CG稍大
 *
 * @par 空间复杂度：
 * O(n) 存储额外的工作向量（r*, p, v, t, s, x等7-8个向量）
 *
 * @par 适用场景：
 * - 非对称/非正定矩阵（Navier-Stokes、涡流场A-V formulation）
 * - 复对称矩阵（时谐电磁场复数系统）
 * - 无法保证对称性的通用稀疏系统
 * - 中等规模到大规模问题（万级到百万级DOF）
 *
 * @par 与CGSolver对比：
 * | 特性 | CGSolver | BiCGSTABSolver |
 * |------|----------|----------------|
 * | 矩阵要求 | SPD/SPSD | 一般稀疏 |
 * | 每次迭代成本 | 1次SpMV | 2次SpMV |
 * | 收敛保证 | 单调递减 | 可能振荡 |
 * | 内存占用 | 4个工作向量 | 7-8个工作向量 |
 * | 稳定性 | 高 | 中等（比BiCG好）|
 *
 * @code
 * // 使用示例：求解涡流场非对称系统
 * numeric::BiCGSTABConfig config;
 * config.tolerance = 1e-6;
 * config.max_iterations = 3000;
 * config.preconditioner = BiCGSTABConfig::PreconditionerType::ILU0;
 *
 * numeric::BiCGSTABSolver solver(config);
 * solver.set_matrix(eddy_current_matrix);
 *
 * auto result = solver.solve(j_source_vector);
 * if (result.status == numeric::SolverStatus::SUCCESS) {
 *     FEEM_INFO("BiCGSTAB求解成功");
 * }
 * @endcode
 *
 * @see CGSolver ScalarAMG EMLinearSolverBase
 */
class BiCGSTABSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数
     * @param config BiCGSTAB求解器配置参数（可选，使用默认值）
     */
    explicit BiCGSTABSolver(const BiCGSTABConfig& config = BiCGSTABConfig{});

    /**
     * @brief 虚析构函数
     */
    ~BiCGSTABSolver() override = default;

    /**
     * @brief 更新求解器配置
     * @param config 新的配置参数
     */
    void set_config(const BiCGSTABConfig& config);

    /**
     * @brief 获取当前配置
     * @return BiCGSTABConfig 当前配置参数副本
     */
    BiCGSTABConfig get_config() const;

    /**
     * @brief 获取收敛历史
     * @return ConvergenceHistory 收敛历史对象
     */
    ConvergenceHistory get_convergence_history() const;

    /**
     * @brief 设置系数矩阵
     * @param A CSR格式系数矩阵（可以是任意方阵，不要求对称）
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 执行BiCGSTAB迭代求解
     * @param b 右端项向量
     * @return SolverResult 求解结果
     *
     * @details 额外的发散检测机制：
     * - 监控连续多次迭代的残差变化趋势
     * - 如果残差持续快速增长（>1000倍），判定为发散并提前终止
     * - 记录详细的警告日志以便调试
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    /**
     * @brief 设置复数系数矩阵（复数版本）
     * @param A CSR格式复数系数矩阵（可以是任意方阵，不要求对称）
     *
     * @details 内部执行：
     * 1. 将CsrMatrix<std::complex<double>>转换为Eigen::SparseMatrix<std::complex<double>>
     * 2. 创建或重置复数Eigen BiCGSTAB求解器实例
     * 3. 配置预条件子（使用IncompleteLUT预条件子，与实数版本一致）
     * 4. 设置容差和最大迭代次数
     * 5. 调用eigen_solver_complex_->compute()分析矩阵
     *
     * @note 此方法会触发预条件子的预计算，可能耗时较长（特别是ILU0）
     * @note 复数版本与实数版本维护独立的内部状态，可独立使用
     */
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;

    /**
     * @brief 执行BiCGSTAB迭代求解复数线性系统（复数版本）
     * @param b 复数右端项向量
     * @return SolverResult 求解结果（复数解存储在 x_complex 字段）
     *
     * @details 求解流程：
     * 1. 输入校验（检查matrix_set_complex_标志和维度匹配）
     * 2. 调用eigen_solver_complex_->solve(b)执行迭代求解
     * 3. 将结果存储到result.x_complex
     * 4. 记录迭代次数和误差信息
     * 5. 判断收敛状态并返回完整结果
     *
     * @par 额外的发散检测机制：
     * - 监控连续多次迭代的残差变化趋势
     * - 如果残差持续快速增长（>1000倍），判定为发散并提前终止
     *
     * @note 必须先调用set_matrix(const CsrMatrix<std::complex<double>>&)设置复数矩阵
     */
    SolverResult solve(const Eigen::VectorXcd& b) override;

    /**
     * @brief 获取求解器名称
     * @return std::string "BiCGSTABSolver_[预条件子类型]"
     */
    std::string get_solver_name() const override;

    /**
     * @brief 清理资源
     */
    void clear() override;

private:
    BiCGSTABConfig config_;                               ///< 配置参数
    Eigen::SparseMatrix<double> eigen_matrix_;             ///< Eigen格式系数矩阵
    bool matrix_set_ = false;                              ///< 矩阵设置标志
    ConvergenceHistory history_;                           ///< 收敛历史

    Eigen::VectorXd jacobi_diag_;                          ///< Jacobi对角线
    Eigen::IncompleteLUT<double> ilu0_precond_;           ///< ILU0预条件子

    // 复数版本内部状态
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;  ///< 缓存的复数Eigen格式系数矩阵
    bool matrix_set_complex_ = false;                      ///< 复数矩阵是否已设置

    /// 复数Eigen BiCGSTAB求解器实例（使用IncompleteLUT预条件子）
    std::unique_ptr<Eigen::BiCGSTAB<
        Eigen::SparseMatrix<std::complex<double>>,
        Eigen::IncompleteLUT<std::complex<double>>>> eigen_solver_complex_;

    void setup_preconditioner();
    SolverResult solve_no_preconditioner(const Eigen::VectorXd& b);
    SolverResult solve_with_jacobi(const Eigen::VectorXd& b);
    SolverResult solve_with_ilu0(const Eigen::VectorXd& b);
};

// ==================== ScalarAMG 代数多重网格预条件子 ====================

/**
 * @class ScalarAMG
 * @brief 标量代数多重网格（Algebraic Multigrid）预条件子
 * @details 为Eigen迭代求解器提供高性能AMG预条件子，通过多层网格层次结构加速Krylov方法的收敛。
 *
 * @par 算法原理：
 * AMG利用矩阵本身的代数结构自动构建粗细网格层次，无需几何网格信息：
 * 1. **粗化（Coarsening）**：基于强连接关系将变量分为C点（粗网点）和F点（细网点）
 * 2. **插值（Interpolation）**：构造从粗网格到细网格的插值算子 P
 * 3. **限制（Restriction）**：R = P^T（Galerkin方法）
 * 4. **光滑（Smoothing）**：在各层使用Jacobi/Gauss-Seidel消除高频误差
 * 5. **多层循环**：V循环或W循环递归地在各层之间传递误差
 *
 * @par 时间复杂度：
 * - 预计算（compute）: O(nnz) 构建层次结构
 * - 应用（apply）: O(nnz) 每次预条件子操作（常数层因子）
 * - 总体: 接近最优 O(n) 求解复杂度
 *
 * @par 空间复杂度：
 * O(nnz * levels) 存储各层矩阵和插值/限制算子
 * 通常levels = O(log n)，总内存约为原始矩阵的1.5-2倍
 *
 * @par 适用场景：
 * - 大规模椭圆型PDE离散系统（泊松方程、扩散方程）
 * - 二阶椭圆型算子的有限元/有限差分离散
 * - 病态系统（条件数 κ >> 10^6）作为预条件子
 * - 百万到千万自由度的超大规模问题
 *
 * @par 性能优势：
 * - 相比Jacobi/ILU0，可减少迭代次数50%-90%
 * - 对网格尺寸不敏感（最优性质）
 * - 特别适合多尺度问题（不同尺度误差分量在不同层级消除）
 *
 * @par 当前实现特点（简化版）：
 * - 采用简单的交替粗化策略（而非Ruge-Stüben经典算法）
 * - 支持直接插值算子构造
 * - 提供Jacobi光滑器和Gauss-Seidel光滑器选项
 * - 完整的V循环和W循环支持
 * - 最粗层使用直接求解器（稠密QR分解）
 *
 * @code
 * // 使用示例：作为CG求解器的预条件子
 * Eigen::SparseMatrix<double> A = ...;  // 大规模稀疏矩阵
 *
 * numeric::ScalarAMG amg(A);
 * amg.set_strong_threshold(0.25);
 * amg.set_smoother_iterations(2);
 * amg.set_cycle_type(numeric::ScalarAMG::CycleType::V_CYCLE);
 * amg.compute(A);  // 预计算层次结构
 *
 * // 在自定义迭代中使用AMG预条件子
 * Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
 * for (int iter = 0; iter < max_iter; ++iter) {
 *     Eigen::VectorXd r = b - A * x;
 *     Eigen::VectorXd z = amg.apply(r);  // AMG预条件子应用
 *     // ... 更新x ...
 * }
 * @endcode
 *
 * @note 未来改进方向：
 * - 实现Ruge-Stüben粗化算法提升粗化质量
 * - 添加Kaczmarz光滑器选项
 * - 支持并行AMG（并行光滑、并行粗化）
 * - 集成到Eigen预条件子框架作为即插即用组件
 *
 * @see CGSolver BiCGSTABSolver EMLinearSolverBase
 */
class ScalarAMG {
public:
    /**
     * @enum CycleType
     * @brief AMG循环类型枚举
     */
    enum class CycleType {
        V_CYCLE,    ///< V循环：细→粗→最粗层精确求解→粗→细（标准选择）
        W_CYCLE     ///< W循环：更激进的误差衰减，每层访问两次粗网格
    };

    /**
     * @enum SmootherType
     * @brief 光滑器类型枚举
     */
    enum class SmootherType {
        JACOBI,         ///< Jacobi（加权雅可比）光滑器
        GAUSS_SEIDEL    ///< Gauss-Seidel（高斯-赛德尔）光滑器
    };

    /**
     * @brief 构造函数
     * @param max_levels 最大网格层数（默认20，通常足够）
     * @param min_coarse_size 最粗层最小节点数（默认50，低于此值停止粗化）
     *
     * @details 创建AMG预条件子实例，但不立即构建层次结构。
     *          必须后续调用 compute(A) 来初始化具体的矩阵相关数据。
     */
    explicit ScalarAMG(int max_levels = 20, int min_coarse_size = 50);

    /**
     * @brief 析构函数
     */
    ~ScalarAMG() = default;

    /**
     * @brief 构建AMG多层次结构（预计算阶段）
     * @param A 输入的细网格系数矩阵（Eigen压缩稀疏矩阵格式）
     *
     * @details 执行以下步骤：
     * 1. 清除旧的层次数据
     * 2. 从最细层（level 0）开始逐层粗化
     * 3. 每层执行：C/F点分离 → 构建插值算子P → Galerkin粗网格 A_c = R*A_f*P
     * 4. 当粗网格尺寸 <= min_coarse_size 或层数 >= max_levels 时停止
     * 5. 记录层次统计信息（各层尺寸、非零元数）
     *
     * @note 此方法可能耗时较长（取决于矩阵规模和粗化质量）
     * @note 调用后可通过 apply(v) 方法使用预条件子
     *
     * @par 粗化策略说明：
     * 当前实现采用简化的交替粗化：
     * - 偶数索引节点选为C点（粗网点）
     * - 奇数索引节点选为F点（细网点）
     * - 这种简单策略在实际应用中效果尚可，但不如Ruge-Stüben算法鲁棒
     */
    void compute(const Eigen::SparseMatrix<double>& A);

    /**
     * @brief 应用AMG预条件子（求解 M^{-1} * v）
     * @param v 输入向量（通常是残差向量 r = b - Ax）
     * @return Eigen::VectorXd 预条件后的向量 z ≈ M^{-1}v
     *
     * @details 执行V循环或W循环：
     * 1. 从最细层（level 0）开始递归
     * 2. 各层执行前光滑 → 残差限制 → 递归粗层求解 → 插值校正 → 后光滑
     * 3. 最粗层使用直接求解器（colPivHouseholderQr）精确求解
     *
     * @note 这是AMG的核心操作，在Krylov方法中每步迭代都会调用
     * @note 时间复杂度: O(nnz)（常数因子约5-10倍单次SpMV）
     */
    Eigen::VectorXd apply(const Eigen::VectorXd& v) const;

    // ==================== 配置方法 ====================

    /**
     * @brief 设置强连接阈值θ
     * @param theta 强连接阈值（典型值0.25，范围[0.1, 0.5]）
     *
     * @details 强连接定义：|a_ij| >= θ * max_{k≠i} |a_ik|
     * - 较小的θ → 更多连接被视为"强连接" → 更多C点 → 粗化较慢
     * - 较大的θ → 只有最强的连接保留 → 粗化更快但可能丢失重要耦合
     *
     * @note 必须在compute()之前调用才能生效
     */
    void set_strong_threshold(double theta);

    /**
     * @brief 设置每层光滑器迭代次数
     * @param iter 前光滑+后光滑的总迭代次数（每方向iter/2次）
     *
     * @details 光滑器的作用是消除当前层的高频误差分量：
     * - iter=1: 快速但光滑效果较弱
     * - iter=2-3: 平衡性能和效果（推荐默认值）
     * - iter>5: 收益递减，增加计算成本
     */
    void set_smoother_iterations(int iter);

    /**
     * @brief 设置循环类型
     * @param cycle V_CYCLE 或 W_CYCLE
     *
     * @details 循环类型影响收敛速度和计算成本：
     * - V_CYCLE: 标准选择，每层只访问一次，成本低
     * - W_CYCLE: 在每个中间层递归两次，收敛更快但成本约2倍
     *
     * @par 选择建议：
     * - 一般问题：V_CYCLE（性价比最高）
     * - 高度各向异性或困难问题：W_CYCLE
     */
    void set_cycle_type(CycleType cycle);

    /**
     * @brief 设置光滑器类型
     * @param smoother JACOBI 或 GAUSS_SEIDEL
     */
    void set_smoother_type(SmootherType smoother);

    /**
     * @brief 设置Jacobi光滑器阻尼因子ω
     * @param omega 阻尼因子（典型值0.67-1.0，默认0.8）
     *
     * @details Jacobi更新公式：x_new = x_old + ω * D^{-1}(b - A*x_old)
     * - ω=1.0: 标准Jacobi（可能不稳定）
     * - ω=0.67: 保证对M矩阵收敛（理论值2/3）
     * - ω=0.8: 经验最佳值（平衡收敛速度和稳定性）
     */
    void set_jacobi_damping(double omega);

    // ==================== 查询方法 ====================

    /**
     * @brief 获取实际构建的网格层数
     * @return int 层数（>=1，至少包含最细层）
     */
    int get_num_levels() const;

    /**
     * @brief 获取各层网格尺寸
     * @return std::vector<int> 第i个元素表示第i层的节点数
     */
    std::vector<int> get_level_sizes() const;

    /**
     * @brief 获取AMG是否已初始化
     * @return bool 已调用compute()且成功构建层次结构返回true
     */
    bool is_initialized() const;

private:
    // ==================== 配置参数 ====================
    int max_levels_;              ///< 最大允许层数
    int min_coarse_size_;         ///< 最粗层最小尺寸阈值
    double strong_threshold_;     ///< 强连接阈值θ
    int smoother_iterations_;     ///< 每层光滑器迭代次数
    CycleType cycle_type_;        ///< 循环类型（V/W）
    SmootherType smoother_type_;  ///< 光滑器类型
    double jacobi_damping_;       ///< Jacobi阻尼因子ω

    // ==================== 多层次数据结构 ====================
    std::vector<Eigen::SparseMatrix<double>> P_;  ///< 插值算子序列 P[level]: fine -> coarse
    std::vector<Eigen::SparseMatrix<double>> R_;  ///< 限制算子序列 R[level] = P[level]^T
    std::vector<Eigen::SparseMatrix<double>> A_;  ///< 各层系数矩阵 A_[0] = 原始矩阵

    bool initialized_ = false;  ///< 是否已完成初始化

    // ==================== 私有辅助方法 ====================

    /**
     * @brief 构建插值算子和C/F点分离
     * @param A 当前层的系数矩阵
     * @return std::pair<P, coarse_indices> 插值算子和粗网点索引列表
     *
     * @details 简化的交替粗化策略：
     * 1. 遍历所有节点，偶数索引选为C点，奇数为F点
     * 2. F点的插值权重基于邻居C点的连接强度
     * 3. C点使用单位插值（保持原值）
     */
    std::pair<Eigen::SparseMatrix<double>, std::vector<int>>
    build_interpolation(const Eigen::SparseMatrix<double>& A);

    /**
     * @brief 执行V循环（核心递归函数）
     * @param level 当前层级（0=最细层）
     * @param x 输入输出：当前近似解（初始为0或上一步结果）
     * @param b 右端项（残差方程的右端项）
     * @param depth 递归深度（用于W循环控制）
     */
    void apply_v_cycle(int level, Eigen::VectorXd& x,
                       const Eigen::VectorXd& b, int depth) const;

    /**
     * @brief Jacobi光滑器一步迭代
     * @param A 当前层矩阵
     * @param x 当前近似解
     * @param b 右端项
     * @return Eigen::VectorXd 光滑后的新近似解
     *
     * @details x_new = x + ω * D^{-1} * (b - A*x)
     */
    Eigen::VectorXd jacobi_smooth(const Eigen::SparseMatrix<double>& A,
                                  const Eigen::VectorXd& x,
                                  const Eigen::VectorXd& b) const;

    /**
     * @brief Gauss-Seidel光滑器一步迭代（前向）
     * @param A 当前层矩阵
     * @param x 当前近似解（就地更新）
     * @param b 右端项
     *
     * @details 逐行更新：x[i] = (b[i] - Σ_{j<i} a[i][j]*x[j] - Σ_{j>i} a[i][j]*x[j]) / a[i][i]
     */
    void gauss_seidel_smooth_forward(const Eigen::SparseMatrix<double>& A,
                                     Eigen::VectorXd& x,
                                     const Eigen::VectorXd& b) const;

    /**
     * @brief Gauss-Seidel光滑器一步迭代（后向）
     * @param A 当前层矩阵
     * @param x 当前近似解（就地更新）
     * @param b 右端项
     *
     * @details 从最后一行向前逐行更新
     */
    void gauss_seidel_smooth_backward(const Eigen::SparseMatrix<double>& A,
                                      Eigen::VectorXd& x,
                                      const Eigen::VectorXd& b) const;
};

// ==================== HiptmairAMG 接口预留 ====================

/**
 * @class HiptmairAMG
 * @brief H(curl)空间专用代数多重网格预条件子（接口预留）
 * @details 适配三维棱边单元H(curl)大规模问题的专用AMG预条件子。
 *          利用∇算子的特殊结构加速Maxwell方程组求解。
 *
 * @par 设计思路（参考论文）：
 * - Hiptmair, "Multigrid method for Maxwell's equations", SIAM J. Numer. Anal., 2000
 * - 核心思想：将H(curl)问题分解为标量部分和梯度修复部分
 * - 标量部分：使用标准ScalarAMG处理旋度算子的主体
 * - 梯度修复：针对null space（梯度场）的特殊光滑器
 * - 通过辅助的H^1（标量）空间的AMG提供梯度场的有效预处理
 *
 * @par 关键技术组件（待实现）：
 * 1. **辅助标量AMG**: 对H^1协调空间（节点元）构建标准AMG
 * 2. **梯度光滑器**: ∇*M_1^{-1}*∇^T 形式的特殊预条件子
 * 3. **嵌套迭代**: 主AMG + 梯度修复的复合结构
 * 4. **边缘聚合粗化**: 基于棱边连接强度的特殊粗化策略
 *
 * @par 适用场景：
 * - 三维时谐Maxwell方程（curl-curl算子）
 * - 棱边单元（Nédélec元）离散系统
 * - 低频/中频电磁场问题（避免高频渐近锁定）
 * - 大规模三维电磁仿真（百万级棱边未知量）
 *
 * @note 当前版本为空实现框架，仅声明纯虚函数接口。
 *       完整实现需要：
 *       - 棱边元-节点元的离散梯度算子
 *       - 辅助标量AMG实例
 *       - 特殊的光滑器和粗化策略
 *       - 大量的测试验证
 *
 * @todo 完整实现以下接口：
 *       - compute(): 构建主AMG + 辅助标量AMG + 梯度算子
 *       - apply(): 执行Hiptmair复合V循环
 *       - 支持配置参数（辅助AMG层数、梯度光滑器迭代次数等）
 *
 * @see ScalarAMG CGSolver EMLinearSolverBase
 */
class HiptmairAMG {
public:
    /**
     * @brief 虚析构函数
     */
    virtual ~HiptmairAMG() = default;

    /**
     * @brief 构建Hiptmair AMG层次结构（TODO: 完整实现）
     * @param A H(curl)空间的系数矩阵（curl-curl算子离散）
     *
     * @details 未来实现应包括：
     * 1. 分析矩阵的null space结构（梯度场模式）
     * 2. 构建主AMG（对标量部分，类似ScalarAMG但适配棱边元）
     * 3. 构建辅助标量AMG（对H^1空间，用于梯度修复）
     * 4. 预计算离散梯度算子及其转置
     * 5. 设置复合光滑器参数
     */
    virtual void compute(const Eigen::SparseMatrix<double>& A) = 0;

    /**
     * @brief 应用Hiptmair AMG预条件子（TODO: 完整实现）
     * @param v 输入向量（残差）
     * @return Eigen::VectorXd 预条件后的向量
     *
     * @details 未来实现应执行：
     * 1. 主AMG V循环（处理旋度算子的高频误差）
     * 2. 梯度修复步骤（消除null space方向的慢收敛分量）
     * 3. 可选的第二轮光滑
     */
    virtual Eigen::VectorXd apply(const Eigen::VectorXd& v) const = 0;

protected:
    /**
     * @brief 保护构造函数（禁止直接实例化，只能通过派生类创建）
     */
    HiptmairAMG() = default;
};

// ==================== 工厂函数 ====================

/**
 * @brief 创建预条件子对象的工厂函数
 * @param type 预条件子类型枚举
 * @param A 系数矩阵（用于预计算ILU0/AMG等需要矩阵信息的预条件子）
 * @return std::unique_ptr<预条件子基类指针> 动态分配的预条件子对象
 *
 * @details 根据枚举类型创建对应的预条件子实例：
 * - NONE: 返回nullptr（表示无预条件子）
 * - JACOBI: 返回对角线缩放预条件子（实际在求解器内部实现）
 * - ILU0: 返回Eigen IncompleteLUT包装
 *
 * @note 当前主要用于内部使用，未来可扩展支持更多预条件子类型
 *
 * @code
 * auto precond = create_preconditioner(
 *     IterativeSolverConfig::PreconditionerType::ILU0, matrix);
 * if (precond) {
 *     // 使用预条件子...
 * }
 * @endcode
 */
std::unique_ptr<int> create_preconditioner(
    IterativeSolverConfig::PreconditionerType type,
    const Eigen::SparseMatrix<double>& A);

} // namespace numeric
