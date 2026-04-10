/**
 * @file solver_backend.h
 * @brief 求解器后端策略接口 - 统一所有求解器后端的契约
 * @details 定义 SolverBackend 抽象基类，作为策略模式的核心接口。
 *          所有求解器后端（Eigen/MUMPS/SuperLU/Pardiso等）必须实现此接口，
 *          实现后端无关的统一调用。
 *
 * @par 设计模式：
 * - **策略模式**：将不同求解算法封装为可互换的策略对象
 * - **开闭原则**：新增求解器后端只需添加新的 SolverBackend 子类，无需修改现有代码
 *
 * @par 架构层次：
 * ```
 * 上层调用代码（不变）
 *     ↓
 * UnifiedDirectSolver（统一调度层，模板方法）
 *     ↓
 * SolverBackend（本文件，策略接口）
 *     ↓
 * EigenLLTBackend / EigenLDLTBackend / EigenLUBackend
 * MUMPSBackend / SuperLUBackend / ...（具体后端实现）
 * ```
 *
 * @par 核心优势：
 * 1. **零代码重复**：公共逻辑在 UnifiedDirectSolver 中只写一次
 * 2. **极高扩展性**：新增后端仅需实现 SolverBackend 接口
 * 3. **运行时灵活**：可动态切换求解器后端而无需重新编译
 * 4. **完全兼容**：上层调用接口完全不变
 *
 * @see UnifiedDirectSolver 统一调度类
 * @see EMSolverFactory 工厂类
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版)
 */

#pragma once

#include "csr_matrix.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace numeric {

/**
 * @class SolverBackend
 * @brief 求解器后端策略抽象基类
 * @details 定义所有求解器后端必须实现的统一接口契约。
 *          采用策略模式封装不同的数值求解算法（Eigen/SuperLU/MUMPS等），
 *          使 UnifiedDirectSolver 能够以统一方式调用任意后端。
 *
 * @par 接口设计原则：
 * - **最小接口**：仅包含核心操作（set_matrix、solve、clear），避免过度设计
 * - **实数/复数统一**：同时支持实数和复数线性系统求解
 * - **状态透明**：通过属性查询方法暴露后端能力（是否支持复数、对称性要求等）
 * - **资源安全**：clear() 方法保证无内存泄漏
 *
 * @par 典型使用流程：
 * @code
 * // 1. 创建后端实例（通过工厂或手动创建）
 * auto backend = std::make_unique<EigenLLTBackend>();
 *
 * // 2. 设置系数矩阵（触发内部分解并缓存）
 * backend->set_matrix(csr_matrix);
 *
 * // 3. 多次求解（复用缓存的分解结果）
 * for (auto& rhs : rhs_list) {
 *     Eigen::VectorXd x = backend->solve_real(rhs);
 *     // 使用解向量...
 * }
 *
 * // 4. 清理资源（或依赖RAII自动清理）
 * backend->clear();
 * @endcode
 *
 * @note 线程安全：此接口的所有方法都不是线程安全的，多线程环境需外部同步
 * @note 内存管理：派生类负责管理内部缓存资源，clear() 必须释放所有动态分配的内存
 */
class SolverBackend {
public:
    /**
     * @brief 虚析构函数
     * @details 确保通过基类指针删除派生类对象时正确调用派生类析构函数，
     *          防止内存泄漏。这是 C++ 多态类的基本要求。
     */
    virtual ~SolverBackend() = default;

    // ========================================================================
    // 核心操作接口（纯虚函数，必须由派生类实现）
    // ========================================================================

    /**
     * @brief 设置实数稀疏系数矩阵并执行预计算（分解/分析）
     * @param A CSR格式的实数稀疏系数矩阵（必须是方阵且已构建完成）
     *
     * @details 此方法是瞬态优化的核心：
     * - 对于直接求解器：触发矩阵分解（LU/Cholesky/LDL^T等），分解结果被缓存
     * - 对于迭代求解器：可能预计算预条件子或提取对角线元素
     * - 后续多次 solve_real() 调用可复用缓存的分解结果，显著提升性能
     *
     * @pre A 必须是已构建完成的 CSR 矩阵（is_built() == true）
     * @post 内部状态更新，后续 solve_real() 调用可用
     *
     * @exception 若矩阵非法或分解失败，应抛出异常或设置内部错误状态
     *
     * @par 典型实现（Eigen LLT 后端）：
     * @code
     * void EigenLLTBackend::set_matrix(const CsrMatrix<double>& A) override {
     *     eigen_matrix_ = SparseConverter::to_eigen(A);
     *     solver_.compute(eigen_matrix_);  // Cholesky 分解并缓存
     * }
     * @endcode
     */
    virtual void set_matrix(const CsrMatrix<double>& A) = 0;

    /**
     * @brief 求解实数线性系统 Ax = b
     * @param b 右端项向量（载荷向量），维度必须与矩阵行数一致
     * @return Eigen::VectorXd 解向量 x
     *
     * @details 利用缓存的分解结果执行高效求解：
     * - 直接求解器：前代 + 回代，时间复杂度 O(n²)
     * - 迭代求解器：执行迭代算法直到收敛
     *
     * @pre 必须先成功调用 set_matrix()
     * @post 返回的解向量满足 Ax ≈ b（在数值精度范围内）
     *
     * @exception 若未设置矩阵或维度不匹配，应抛出 std::invalid_argument
     * @exception 若求解失败（矩阵奇异等），应抛出 std::runtime_error
     */
    virtual Eigen::VectorXd solve_real(const Eigen::VectorXd& b) = 0;

    /**
     * @brief 设置复数稀疏系数矩阵并执行预计算（分解/分析）
     * @param A CSR格式的复数稀疏系数矩阵（必须是方阵且已构建完成）
     *
     * @details 复数版本的 set_matrix，用于时谐场电磁分析等场景：
     * - 支持复数稀疏矩阵的输入和内部表示
     * - 触发复数矩阵分解（LU/QR等），分解结果被缓存
     * - 后续多次 solve_complex() 调用可复用缓存的分解结果
     *
     * @pre A 必须是已构建完成的 CSR 复数矩阵
     * @post 内部复数状态更新，后续 solve_complex() 调用可用
     *
     * @note 复数矩阵常用于频域电磁分析、涡流计算等时谐场问题
     */
    virtual void set_matrix(const CsrMatrix<std::complex<double>>& A) = 0;

    /**
     * @brief 求解复数线性系统 Ax = b
     * @param b 复数右端项向量（载荷向量），维度必须与矩阵行数一致
     * @return Eigen::VectorXcd 复数解向量 x
     *
     * @details 复数版本的 solve，利用缓存的复数分解结果执行求解。
     *          返回的解向量的实部和虚部分别对应原问题的实部和虚部解。
     *
     * @pre 必须先成功调用 set_matrix(const CsrMatrix<std::complex<double>>&)
     * @post 返回的复数解向量满足 Ax ≈ b（在数值精度范围内）
     */
    virtual Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) = 0;

    /**
     * @brief 清理求解器资源并重置状态
     * @details 释放所有内部缓存的资源，将求解器恢复到初始未初始化状态：
     * - 释放矩阵分解结果（L/U因子、符号分析信息等）
     * - 释放预条件子数据结构
     * - 重置临时工作空间
     * - 清除所有内部状态标志
     *
     * @note 调用后必须重新执行 set_matrix() 才能进行后续求解
     * @note 在析构函数中自动调用，也可手动调用以提前释放内存
     */
    virtual void clear() = 0;

    // ========================================================================
    // 属性查询接口（提供后端元信息）
    // ========================================================================

    /**
     * @brief 获取后端的唯一名称标识
     * @return std::string 后端名称字符串（如 "EigenLLT"、"MUMPS"、"SuperLU_MT" 等）
     *
     * @details 用于日志输出、性能统计和调试诊断。
     *          名称应能清晰标识求解器类型和关键参数配置。
     *
     * @note 建议命名格式："库名_算法"（如 "Eigen_SimplicialLLT"、"MUMPS_LU"）
     */
    virtual std::string get_backend_name() const = 0;

    /**
     * @brief 查询后端是否支持复数运算
     * @return true 支持复数矩阵和复数求解
     * @return false 仅支持实数运算
     *
     * @details 某些专用求解器可能仅支持实数运算（如某些定制的实数Cholesky实现）。
     *          此方法允许 UnifiedDirectSolver 在运行时检查能力并给出友好错误提示。
     *
     * @note 大多数现代求解器后端都应返回 true（Eigen/SuperLU/MUMPS 均支持复数）
     */
    virtual bool supports_complex() const = 0;

    /**
     * @brief 查询后端是否仅支持对称矩阵
     * @return true 仅接受对称（或Hermitian）矩阵作为输入
     * @return false 支持一般非对称矩阵
     *
     * @details 此方法用于输入校验和错误提示：
     * - 对于 Cholesky (LLT) 和 LDL^T 后端，返回 true（要求对称输入）
     * - 对于 LU 分解后端，返回 false（支持一般方阵）
     *
     * @note UnifiedDirectSolver 可据此在上层进行对称性预检
     */
    virtual bool is_symmetric_only() const = 0;

    // ========================================================================
    // 矩阵访问接口（用于残差计算等高级功能）
    // ========================================================================

    /**
     * @brief 获取缓存的实数系数矩阵（Eigen CSC 格式）
     * @return const Eigen::SparseMatrix<double>* 指向内部缓存矩阵的常量指针
     * @return nullptr 如果实数矩阵未设置
     *
     * @details 此方法允许外部代码（如 UnifiedDirectSolver）访问缓存的矩阵，
     *          用于计算精确的残差范数 ||Ax - b||_2 等高级功能。
     *
     * @warning 返回的是内部存储的指针，调用者不应修改或释放此内存
     * @note 默认实现返回 nullptr（派生类可根据需要重写）
     */
    virtual const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const { return nullptr; }

    /**
     * @brief 获取缓存的复数系数矩阵（Eigen CSC 格式）
     * @return const Eigen::SparseMatrix<std::complex<double>>* 指向内部缓存矩阵的常量指针
     * @return nullptr 如果复数矩阵未设置
     *
     * @details 复数版本的矩阵访问接口。
     *
     * @warning 返回的是内部存储的指针，调用者不应修改或释放此内存
     */
    virtual const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const { return nullptr; }
};

} // namespace numeric
