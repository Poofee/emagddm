/**
 * @file unified_direct_solver.h
 * @brief 统一直接求解器调度类 - 模板方法模式的核心实现
 * @details 实现 UnifiedDirectSolver 类，作为策略模式 + 模板方法模式的统一调度层。
 *          **核心价值**：将原3个求解器类中 70%+ 的重复代码（公共逻辑）合并到此类中，
 *          所有直接求解器都通过此类统一调度，零代码重复！
 *
 * @par 设计模式：
 * - **模板方法模式**：本类定义统一的求解算法骨架（set_matrix → solve → 残差计算）
 * - **策略模式**：通过 SolverBackend 接口调用具体后端实现
 * - **组合模式**：持有 SolverBackend 智针，运行时可动态切换后端
 *
 * @par 架构层次：
 * ```
 * 上层调用代码（完全不变）
 *     ↓
 * UnifiedDirectSolver（本文件，统一调度层）
 *   ├── 输入校验（统一实现）
 *   ├── 耗时统计（统一实现）
 *   ├── 残差计算（统一实现）
 *   ├── 错误处理（统一实现）
 *   └── 委托给 backend_->solve_real() / solve_complex()
 *     ↓
 * SolverBackend（策略接口）
 *     ↓
 * EigenLLTBackend / MUMPSBackend / SuperLUBackend / ...（具体后端）
 * ```
 *
 * @par 消除的重复代码（对比重构前）：
 * ✅ SymmetricDirectSolver::set_matrix() 的公共逻辑 → 统一到本类
 * ✅ SymmetricIndefiniteDirectSolver::set_matrix() 的公共逻辑 → 同上
 * ✅ GeneralDirectSolver::set_matrix() 的公共逻辑 → 同上
 * ✅ 3个类的 solve() 公共逻辑（输入校验、残差、计时）→ 统一到本类
 * ✅ 3个类的 clear() 公共逻辑 → 统一到本类
 * ✅ 复数版本的上述所有方法 → 统一到本类
 *
 * @par 兼容性保证：
 * - 继承 EMLinearSolverBase 基类（接口完全不变）
 * - 返回值类型 SolverResult 完全兼容
 * - 上层调用代码无需任何修改
 *
 * @par 典型使用方式（通过工厂创建）：
 * @code
 * // 方式1：通过工厂创建（推荐）
 * auto solver = EMSolverFactory::create_solver(EMSolverFactory::SolverType::SYMMETRIC_DIRECT);
 * solver->set_matrix(A);
 * auto result = solver->solve(b);
 *
 * // 方式2：手动创建（高级用法）
 * auto backend = std::make_unique<MUMPSBackend>(MUMPSBackend::MatrixType::SYMMETRIC_POSITIVE_DEFINITE);
 * auto solver = std::make_unique<UnifiedDirectSolver>(std::move(backend));
 * @endcode
 *
 * @see EMLinearSolverBase 抽象基类
 * @see SolverBackend 策略接口
 * @see EMSolverFactory 工厂类
 *
 * @author Poofee
 * @date 2026-04-10
 * @version 2.0 (重构版)
 */

#pragma once

#include "em_linear_solver.h"
#include "solver_backend.h"
#include <chrono>
#include <memory>
#include <string>

namespace numeric {

/**
 * @class UnifiedDirectSolver
 * @brief 统一直接求解器调度类（模板方法模式）
 * @details 继承 EMLinearSolverBase，封装所有直接求解器的公共逻辑。
 *          通过组合模式持有 SolverBackend 对象，将具体求解算法委托给后端处理。
 *
 * @par 核心职责：
 * 1. **输入校验**：统一验证矩阵和向量的合法性（维度、方阵检查等）
 * 2. **后端委托**：将 set_matrix/solve 调用转发给具体的 SolverBackend 实现
 * 3. **性能监控**：统一记录求解耗时、计算残差范数、输出详细日志
 * 4. **错误处理**：统一捕获异常并转换为 SolverResult 状态码
 * 5. **资源管理**：统一管理后端对象的生命周期
 *
 * @par 瞬态复用优化：
 * 本类完整支持瞬态分析场景的性能优化：
 * @code
 * UnifiedDirectSolver solver(std::make_unique<EigenLLTBackend>());
 *
 * solver.set_matrix(K);      // 触发分解（耗时操作，仅执行一次）
 *
 * for (int t = 0; t < time_steps; t++) {
 *     auto result = solver.solve(F[t]);  // 快速求解（复用缓存的分解因子）
 *     // 使用 result.x ...
 * }
 * @endcode
 */
class UnifiedDirectSolver : public EMLinearSolverBase {
public:
    /**
     * @brief 构造函数：传入求解器后端策略对象
     * @param backend 求解器后端的智能指针（所有权转移至本对象）
     *
     * @details 通过组合模式关联具体的求解器后端。
     *          本构造函数接受任意 SolverBackend 派生类的智能指针，
     *          实现运行时灵活选择求解算法而无需修改本类代码。
     *
     * @pre backend 不能为 nullptr（否则抛出 std::invalid_argument）
     * @post 本对象拥有 backend 的唯一所有权
     *
     * @par 示例：
     * @code
     * // 创建 Eigen LLT 后端的统一求解器
     * auto solver1 = UnifiedDirectSolver(std::make_unique<EigenLLTBackend>());
     *
     * // 创建 MUMPS 后端的统一求解器
     * auto solver2 = UnifiedDirectSolver(
     *     std::make_unique<MUMPSBackend>(MUMPSBackend::MatrixType::SYMMETRIC_POSITIVE_DEFINITE)
     * );
     * @endcode
     */
    explicit UnifiedDirectSolver(std::unique_ptr<SolverBackend> backend);

    /**
     * @brief 虚析构函数（默认实现）
     * @details 自动释放 backend_ 智能指针管理的资源
     */
    ~UnifiedDirectSolver() override = default;

    // ===== 实现 EMLinearSolverBase 纯虚函数接口 =====

    /**
     * @brief 设置实数系数矩阵并触发分解（统一实现）
     * @param A CSR格式的实数稀疏系数矩阵
     *
     * @details 统一的矩阵设置流程：
     * 1. 输入合法性校验（维度、方阵检查）
     * 2. 委托给 backend_->set_matrix(A) 执行实际的分解
     * 3. 记录耗时和日志信息
     *
     * @note 此方法是瞬态优化的核心：一次分解，多次求解复用
     */
    void set_matrix(const CsrMatrix<double>& A) override;

    /**
     * @brief 求解实数线性系统 Ax = b（统一实现）
     * @param b 右端项向量
     * @return SolverResult 包含解向量、状态码、残差范数、耗时等信息
     *
     * @details 统一的求解流程：
     * 1. 输入校验（矩阵已设置、维度匹配）
     * 2. 委托给 backend_->solve_real(b) 执行实际求解
     * 3. 计算残差范数 ||Ax - b|| 用于质量评估
     * 4. 记录求解耗时和性能指标
     * 5. 统一错误处理和日志输出
     */
    SolverResult solve(const Eigen::VectorXd& b) override;

    /**
     * @brief 设置复数系数矩阵并触发分解（统一实现）
     * @param A CSR格式的复数稀疏系数矩阵
     *
     * @details 复数版本的 set_matrix，用于时谐场电磁分析等场景。
     *          流程与实数版本完全对称。
     */
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override;

    /**
     * @brief 求解复数线性系统 Ax = b（统一实现）
     * @param b 复数右端项向量
     * @return SolverResult 复数解存储在 x_complex 字段
     *
     * @details 复数版本的 solve，流程与实数版本完全对称。
     */
    SolverResult solve(const Eigen::VectorXcd& b) override;

    /**
     * @brief 获取求解器名称标识
     * @return std::string 格式："UnifiedDirect_[后端名]"（如 "UnifiedDirect_EigenLLT"）
     */
    std::string get_solver_name() const override;

    /**
     * @brief 清理所有资源并重置状态（统一实现）
     * @details 委托给 backend_->clear() 并重置内部状态标志
     */
    void clear() override;

private:
    std::unique_ptr<SolverBackend> backend_;  ///< 求解器后端策略对象（组合模式）

    bool matrix_set_ = false;       ///< 实数矩阵是否已设置
    bool matrix_set_complex_ = false;  ///< 复数矩阵是否已设置

    /**
     * @brief 计算实数解的残差范数 ||Ax - b||_2
     * @param x 解向量
     * @param b 右端项向量
     * @return double 残差的2-范数
     *
     * @note 此方法需要访问缓存的 Eigen 矩阵，因此由本类而非后端实现
     */
    double compute_residual_norm_real(const Eigen::VectorXd& x, const Eigen::VectorXd& b);

    /**
     * @brief 计算复数解的残差范数 ||Ax - b||_2
     * @param x 复数解向量
     * @param b 复数右端项向量
     * @return double 残差的2-范数
     */
    double compute_residual_norm_complex(const Eigen::VectorXcd& x, const Eigen::VectorXcd& b);
};

} // namespace numeric
