/**
 * @file eigen_solver_backends.h
 * @brief Eigen 求解器后端策略实现 - LLT/LDLT/LU 三种分解算法（CRTP重构版）
 * @details 使用 CRTP (Curiously Recurring Template Pattern) 模板基类消除三个后端的重复代码。
 *          核心设计：EigenSolverBackendBase<Derived, RealSolver, ComplexSolver> 封装所有公共逻辑，
 *          派生类仅需指定模板参数和元信息，实现零重复的目标。
 *
 * @par 重构前问题：
 * - 原始实现中三个类（LLT/LDLT/LU）的 set_matrix/solve/clear 方法存在 ~95% 的代码重复
 * - 仅差异点为：求解器类型（模板参数）、错误消息中的类名
 * - clear() 方法完全相同，违反 DRY 原则
 *
 * @par 重构方案（CRTP + 策略模式）：
 * 1. **模板基类 EigenSolverBackendBase**：
 *    - 封装所有成员变量（eigen_matrix_, solver_, matrix_set_ 等）
 *    - 实现所有公共方法逻辑（set_matrix, solve_real/complex, clear）
 *    - 通过 CRTP 获取派生类的静态信息（backend_name, is_symmetric_only）
 *    - 模板参数 RealSolver/ComplexSolver 具体化求解器类型
 *
 * 2. **派生类薄包装**（每个仅 ~15 行）：
 *    - 继承基类，指定具体的 Eigen 求解器类型作为模板参数
 *    - 提供静态元信息（名称、对称性要求）
 *    - 自动获得完整的 SolverBackend 接口实现
 *
 * @par 设计优势：
 * - **零代码重复**：公共逻辑在基类中只写一次
 * - **编译期安全**：模板实例化时检查类型正确性，无运行时开销
 * - **极高扩展性**：新增 Eigen 求解器后端只需添加 ~10 行代码
 * - **完美兼容**：对外接口（SolverBackend）保持不变，上层代码无需修改
 *
 * @par 性能特征：
 * | 后端类型     | 时间复杂度    | 数值稳定性  | 内存占用   | 对称性要求 |
 * |------------|--------------|-----------|----------|----------|
 * | LLT        | O(n³/3)      | 最优       | ~1.5n²    | 必须SPD   |
 * | LDLT       | O(n³/3)      | 良好       | ~1.5n²    | 对称即可  |
 * | LU         | O(n³/3)      | 良好(选主元)| ~2n²     | 无要求    |
 *
 * @see SolverBackend 策略接口基类
 * @see UnifiedDirectSolver 统一调度类
 *
 * @author Poofee
 * @date 2026-04-11
 * @version 3.0 (CRTP重构版)
 */

#pragma once

#include "solver_backend.h"
#include "sparse_converter.h"
#include "logger_factory.hpp"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <memory>
#include <string>

namespace numeric {

// ============================================================================
// CRTP 模板基类：封装所有 Eigen 求解器的公共逻辑
// ============================================================================

/**
 * @class EigenSolverBackendBase
 * @brief Eigen 求解器后端 CRTP 基类（消除重复代码的核心）
 * @details 使用 CRTP 模式将三个 Eigen 后端（LLT/LDLT/LU）的公共逻辑抽象到基类，
 *          派生类通过模板参数指定具体的求解器类型，并通过 CRTP 提供特定的元信息。
 *
 * @tparam Derived 派生类类型（用于 CRTP 静态多态）
 * @tparam RealSolver 实数 Eigen 求解器类型（如 SimplicialLLT<SparseMatrix<double>>）
 * @tparam ComplexSolver 复数 Eigen 求解器类型（如 SimplicialLLT<SparseMatrix<complex>>）
 *
 * @par CRTP 机制说明：
 * - Derived 类必须提供两个静态方法：
 *   - `static std::string backend_name()` 返回后端名称字符串
 *   - `static bool symmetric_only()` 返回是否仅支持对称矩阵
 * - 基类通过 `static_cast<Derived*>(this)` 调用这些方法（编译期解析，零开销）
 *
 * @par 成员变量（统一管理）：
 * - eigen_matrix_ / eigen_matrix_complex_：缓存的实数/复数系数矩阵
 * - solver_ / solver_complex_：实数/复数分解器实例
 * - matrix_set_ / matrix_set_complex_：状态标志位
 *
 * @par 公共方法实现：
 * - set_matrix()：CSR → Eigen 格式转换 + 触发分解 + 错误处理
 * - solve_real/complex()：利用缓存的分解结果执行高效求解
 * - clear()：释放所有资源并重置状态
 * - get_eigen_matrix_real/complex()：提供矩阵访问接口
 */
template<typename Derived, typename RealSolver, typename ComplexSolver>
class EigenSolverBackendBase : public SolverBackend {
public:
    /**
     * @brief 默认构造函数，初始化内部状态
     */
    EigenSolverBackendBase() = default;

    /**
     * @brief 析构函数，自动调用 clear() 清理资源（RAII）
     */
    ~EigenSolverBackendBase() override { clear(); }

    // ========================================================================
    // 核心操作接口实现（从 SolverBackend 继承的虚函数）
    // ========================================================================

    /**
     * @brief 设置实数稀疏系数矩阵并执行预计算（分解）
     * @param A CSR格式的实数稀疏系数矩阵
     *
     * @details 统一的实数矩阵设置流程（适用于所有 Eigen 直接求解器）：
     * 1. 日志记录输入维度信息
     * 2. CSR → Eigen 稀疏矩阵格式转换（使用 SparseConverter 工具类）
     * 3. 压缩存储优化（makeCompressed）
     * 4. 创建具体类型的求解器实例（由模板参数 RealSolver 决定）
     * 5. 执行矩阵分解并缓存分解结果
     * 6. 检查分解成功与否（Eigen::Success 状态码）
     * 7. 异常处理和资源清理
     *
     * @note 此方法利用了 C++ 模板的静态多态特性：
     *       - RealSolver 模板参数在编译期确定具体求解器类型
     *       - Derived::backend_name() 在编译期内联展开（零运行时开销）
     */
    void set_matrix(const CsrMatrix<double>& A) override {
        FEEM_DEBUG("{}::set_matrix - 开始设置实数矩阵, 维度: {}x{}",
                   Derived::backend_name(), A.rows(), A.cols());

        try {
            eigen_matrix_ = SparseConverter::to_eigen(A);
            eigen_matrix_.makeCompressed();

            solver_ = std::make_unique<RealSolver>();
            solver_->compute(eigen_matrix_);

            if (solver_->info() != Eigen::Success) {
                FEEM_ERROR("{}::set_matrix - 分解失败（矩阵可能非正定或奇异）",
                           Derived::backend_name());
                solver_.reset();
                matrix_set_ = false;
                return;
            }

            matrix_set_ = true;
            FEEM_DEBUG("{}::set_matrix - 分解成功", Derived::backend_name());
        } catch (const std::exception& e) {
            FEEM_ERROR("{}::set_matrix - 异常: {}", Derived::backend_name(), e.what());
            solver_.reset();
            matrix_set_ = false;
        }
    }

    /**
     * @brief 求解实数线性系统 Ax = b
     * @param b 右端项向量
     * @return Eigen::VectorXd 解向量 x
     *
     * @details 利用缓存的实数分解结果执行高效前代+回代求解。
     *          包含前置条件检查（矩阵已设置、维度匹配）。
     */
    Eigen::VectorXd solve_real(const Eigen::VectorXd& b) override {
        if (!matrix_set_ || !solver_) {
            throw std::runtime_error(
                Derived::backend_name() + ": 矩阵未设置或分解失败，请先调用 set_matrix()");
        }
        if (b.size() != eigen_matrix_.rows()) {
            throw std::invalid_argument(
                Derived::backend_name() + ": 右端项向量维度不匹配");
        }

        return solver_->solve(b);
    }

    /**
     * @brief 设置复数稀疏系数矩阵并执行预计算（分解）
     * @param A CSR格式的复数稀疏系数矩阵
     *
     * @details 复数版本的 set_matrix，用于时谐场电磁分析等场景。
     *          流程与实数版本完全对称，仅数据类型不同。
     */
    void set_matrix(const CsrMatrix<std::complex<double>>& A) override {
        FEEM_DEBUG("{}::set_matrix(复数) - 开始设置复数矩阵, 维度: {}x{}",
                   Derived::backend_name(), A.rows(), A.cols());

        try {
            eigen_matrix_complex_ = SparseConverter::to_eigen_complex(A);
            eigen_matrix_complex_.makeCompressed();

            solver_complex_ = std::make_unique<ComplexSolver>();
            solver_complex_->compute(eigen_matrix_complex_);

            if (solver_complex_->info() != Eigen::Success) {
                FEEM_ERROR("{}::set_matrix(复数) - 复数分解失败",
                           Derived::backend_name());
                solver_complex_.reset();
                matrix_set_complex_ = false;
                return;
            }

            matrix_set_complex_ = true;
            FEEM_DEBUG("{}::set_matrix(复数) - 复数分解成功", Derived::backend_name());
        } catch (const std::exception& e) {
            FEEM_ERROR("{}::set_matrix(复数) - 异常: {}",
                       Derived::backend_name(), e.what());
            solver_complex_.reset();
            matrix_set_complex_ = false;
        }
    }

    /**
     * @brief 求解复数线性系统 Ax = b
     * @param b 复数右端项向量
     * @return Eigen::VectorXcd 复数解向量 x
     *
     * @details 复数版本的 solve，利用缓存的复数分解结果执行求解。
     */
    Eigen::VectorXcd solve_complex(const Eigen::VectorXcd& b) override {
        if (!matrix_set_complex_ || !solver_complex_) {
            throw std::runtime_error(
                Derived::backend_name() + ": 复数矩阵未设置或分解失败，请先调用 set_matrix(复数)");
        }
        if (b.size() != eigen_matrix_complex_.rows()) {
            throw std::invalid_argument(
                Derived::backend_name() + ": 复数右端项向量维度不匹配");
        }

        return solver_complex_->solve(b);
    }

    /**
     * @brief 清理求解器资源并重置状态
     * @details 释放所有内部缓存（矩阵、分解器、标志位），恢复到初始未初始化状态。
     *          此方法在析构函数中自动调用（RAII），也可手动调用以提前释放内存。
     */
    void clear() override {
        FEEM_DEBUG("{}::clear - 释放所有资源", Derived::backend_name());

        solver_.reset();
        solver_complex_.reset();
        eigen_matrix_ = Eigen::SparseMatrix<double>();
        eigen_matrix_complex_ = Eigen::SparseMatrix<std::complex<double>>();
        matrix_set_ = false;
        matrix_set_complex_ = false;
    }

    // ========================================================================
    // 属性查询接口实现
    // ========================================================================

    /**
     * @brief 获取后端唯一名称标识
     * @return std::string 后端名称（如 "EigenLLT"、"EigenLDLT"、"EigenLU"）
     */
    std::string get_backend_name() const override { return Derived::backend_name(); }

    /**
     * @brief 查询是否支持复数运算
     * @return true 所有 Eigen 后端均支持复数
     */
    bool supports_complex() const override { return true; }

    /**
     * @brief 查询是否仅支持对称矩阵
     * @return true/false 由派生类的 static symmetric_only() 决定
     */
    bool is_symmetric_only() const override { return Derived::symmetric_only(); }

    // ========================================================================
    // 矩阵访问接口实现
    // ========================================================================

    /**
     * @brief 获取缓存的实数系数矩阵指针
     * @return const Eigen::SparseMatrix<double>* 内部矩阵常量指针
     */
    const Eigen::SparseMatrix<double>* get_eigen_matrix_real() const override {
        return &eigen_matrix_;
    }

    /**
     * @brief 获取缓存的复数系数矩阵指针
     * @return const Eigen::SparseMatrix<std::complex<double>>* 内部复数矩阵常量指针
     */
    const Eigen::SparseMatrix<std::complex<double>>* get_eigen_matrix_complex() const override {
        return &eigen_matrix_complex_;
    }

protected:
    /**
     * @brief 禁止拷贝构造（防止资源双重释放）
     */
    EigenSolverBackendBase(const EigenSolverBackendBase&) = delete;

    /**
     * @brief 禁止拷贝赋值（防止资源双重释放）
     */
    EigenSolverBackendBase& operator=(const EigenSolverBackendBase&) = delete;

    // ========================================================================
    // 成员变量（统一管理，避免在每个派生类中重复声明）
    // ========================================================================

    Eigen::SparseMatrix<double> eigen_matrix_;  ///< 缓存的实数系数矩阵（Eigen CSC格式）
    Eigen::SparseMatrix<std::complex<double>> eigen_matrix_complex_;  ///< 缓存的复数系数矩阵

    std::unique_ptr<RealSolver> solver_;  ///< 实数分解器实例（由模板参数决定具体类型）
    std::unique_ptr<ComplexSolver> solver_complex_;  ///< 复数分解器实例

    bool matrix_set_ = false;       ///< 实数矩阵是否已设置并分解成功
    bool matrix_set_complex_ = false;  ///< 复数矩阵是否已设置并分解成功
};


// ============================================================================
// 派生类定义（薄包装，仅指定模板参数和元信息）
// ============================================================================

/**
 * @class EigenLLTBackend
 * @brief Eigen Cholesky 分解后端（对称正定矩阵）- CRTP 版本
 * @details 使用 SimplicialLLT 执行稀疏 Cholesky 分解 A = LL^T，
 *          专为对称正定(SPD)矩阵设计，具有最优的数值稳定性和效率。
 *
 * @par 从基类继承的功能（无需重复实现）：
 * - set_matrix() / solve_real() / solve_complex() / clear()
 * - get_eigen_matrix_real() / get_eigen_matrix_complex()
 * - supports_complex() → 返回 true
 *
 * @par 本类仅需提供的信息：
 * - backend_name() = "EigenLLT"
 * - symmetric_only() = true（仅接受 SPD 矩阵）
 * - 模板参数：RealSolver=SimplicialLLT<double>, ComplexSolver=SimplicialLLT<complex>
 *
 * @code
 * auto backend = std::make_unique<EigenLLTBackend>();
 * backend->set_matrix(spd_stiffness_matrix);  // Cholesky 分解
 * Eigen::VectorXd x = backend->solve_real(rhs_vector);
 * @endcode
 */
class EigenLLTBackend : public EigenSolverBackendBase<
    EigenLLTBackend,
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>,
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>>> {
public:
    using Base = EigenSolverBackendBase<
        EigenLLTBackend,
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>,
        Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>>>;

    EigenLLTBackend() = default;
    ~EigenLLTBackend() override = default;

    /**
     * @brief CRTP 静态方法：返回后端名称
     * @return "EigenLLT"
     */
    static std::string backend_name() { return "EigenLLT"; }

    /**
     * @brief CRTP 静态方法：查询是否仅支持对称矩阵
     * @return true LLT 要求输入必须是对称正定(SPD)矩阵
     */
    static bool symmetric_only() { return true; }
};

/**
 * @class EigenLDLTBackend
 * @brief Eigen LDL^T 分解后端（对称不定矩阵）- CRTP 版本
 * @details 使用 SimplicialLDLT 执行稀疏 LDL^T 分解 A = LDL^T，
 *          专为对称但未必正定的矩阵设计，可处理接近奇异的情况。
 *
 * @par 与 EigenLLTBackend 的区别：
 * - LLT 要求输入必须 SPD（否则分解失败）
 * - LDLT 允许输入对称但不定（D 对角元可为负或零）
 * - LDLT 数值稳定性略低于 LLT（对于 SPD 矩阵）
 *
 * @par 本类仅需提供的信息：
 * - backend_name() = "EigenLDLT"
 * - symmetric_only() = true（要求对称，但不要求正定）
 * - 模板参数：RealSolver=SimplicialLDLT<double>, ComplexSolver=SimplicialLDLT<complex>
 */
class EigenLDLTBackend : public EigenSolverBackendBase<
    EigenLDLTBackend,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>> {
public:
    using Base = EigenSolverBackendBase<
        EigenLDLTBackend,
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>,
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>>;

    EigenLDLTBackend() = default;
    ~EigenLDLTBackend() override = default;

    static std::string backend_name() { return "EigenLDLT"; }
    static bool symmetric_only() { return true; }
};

/**
 * @class EigenLUBackend
 * @brief Eigen LU 分解后端（通用非对称矩阵）- CRTP 版本
 * @details 使用 SparseLU 执行带部分选主元的 LU 分解 PA = LU，
 *          是最通用的直接求解器后端，适用于所有类型的方阵（包括非对称矩阵）。
 *
 * @par 适用场景：
 * - 非对称系统（一般方阵）
 * - 复数系统（时谐场分析）
 * - 无法保证对称性的通用稀疏系统
 *
 * @par 内存开销提示：
 * - LU 分解的 fill-in 通常高于 Cholesky/LDL^T（不对称性导致更多非零元）
 * - 大规模非对称问题建议考虑迭代法（GMRES/BiCGSTAB + 预条件子）
 *
 * @par 本类仅需提供的信息：
 * - backend_name() = "EigenLU"
 * - symmetric_only() = false（支持任意方阵）
 * - 模板参数：RealSolver=SparseLU<double>, ComplexSolver=SparseLU<complex>
 */
class EigenLUBackend : public EigenSolverBackendBase<
    EigenLUBackend,
    Eigen::SparseLU<Eigen::SparseMatrix<double>>,
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>> {
public:
    using Base = EigenSolverBackendBase<
        EigenLUBackend,
        Eigen::SparseLU<Eigen::SparseMatrix<double>>,
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>;

    EigenLUBackend() = default;
    ~EigenLUBackend() override = default;

    static std::string backend_name() { return "EigenLU"; }
    static bool symmetric_only() { return false; }
};

} // namespace numeric
