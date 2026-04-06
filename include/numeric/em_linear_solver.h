/**
 * @file em_linear_solver.h
 * @brief 电磁场线性求解器抽象基类与数据结构定义
 * @details 定义统一的线性求解器接口规范，支持直接法和迭代法多后端，
 *          为FETI-DP子域求解提供标准化的求解器抽象层。
 *          包含求解状态枚举、后端类型枚举、结果结构体和抽象基类定义。
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

#include "csr_matrix.hpp"
#include <Eigen/Dense>
#include <string>

namespace numeric {

/**
 * @enum SolverStatus
 * @brief 线性求解器状态枚举
 * @details 定义求解过程中可能出现的各种状态，用于标识求解结果的类型和错误原因
 */
enum class SolverStatus {
    SUCCESS,                ///< 求解成功，得到有效解向量
    NUMERICAL_ERROR,        ///< 数值错误：矩阵奇异、LU/Cholesky分解失败、数值溢出等
    DIVERGED,               ///< 迭代求解发散：残差范数持续增大或超过阈值
    MAX_ITER_REACHED,       ///< 达到最大迭代次数但未收敛，返回当前近似解
    INVALID_INPUT           ///< 输入非法：矩阵维度不匹配、向量长度不一致、矩阵未初始化等
};

/**
 * @enum DirectBackendType
 * @brief 直接求解器后端类型枚举
 * @details 支持多种高性能直接求解后端，可根据问题规模和硬件配置灵活选择
 */
enum class DirectBackendType {
    EIGEN,       ///< Eigen内置稀疏直接求解器（默认）：无需额外依赖，适合中小规模问题（<10万DOF），
                 ///< 支持LU、Cholesky、QR分解，单线程性能优秀
    SUPERLU,     ///< SuperLU开源直接求解器：支持并行LU分解，适合中大规模稀疏系统，
                 ///< 提供SuperLU_MT多线程版本和SuperLU_DIST分布式版本
    MUMPS        ///< MUMPS并行分布式直接求解器：支持MPI并行，适合大规模分布式计算场景，
                 ///< 具有优秀的内存管理和负载均衡能力，特别适用于FETI-DP子域求解
};

/**
 * @struct SolverResult
 * @brief 线性求解器结果结构体
 * @details 封装求解器的完整输出信息，包括解向量、状态码、性能指标和错误信息
 */
struct SolverResult {
    Eigen::VectorXd x;      ///< 求解结果：未知量向量，维度与输入矩阵行数一致
    SolverStatus status;    ///< 求解状态：指示求解是否成功及失败原因
    int iterations = 0;     ///< 迭代次数：直接求解器为0，迭代求解器记录实际迭代轮数
    double residual_norm = 0.0; ///< 最终残差范数：||Ax - b||_2，用于评估解的质量
    double solve_time_ms = 0.0; ///< 求解耗时：从开始求解到返回结果的 wall-clock 时间（毫秒）
    std::string error_msg;  ///< 错误信息：当status非SUCCESS时，包含详细的错误描述
};

/**
 * @class EMLinearSolverBase
 * @brief 电磁场线性求解器抽象基类
 * @details 定义统一的线性系统 Ax = b 求解接口，支持瞬态复用机制（矩阵不变时避免重复分解）。
 *          所有具体求解器实现（直接法/迭代法）必须继承此类并实现纯虚函数接口。
 *          设计目标：
 *          - 为FETI-DP子域求解提供标准化接口
 *          - 支持多种数值后端的无缝切换
 *          - 提供完整的求解状态反馈和性能监控
 *
 * @note 线程安全性：set_matrix()和solve()方法在多线程环境下需外部同步保护
 * @note 内存管理：派生类负责管理内部缓存资源，clear()方法必须释放所有动态分配的内存
 */
class EMLinearSolverBase {
public:
    /**
     * @brief 虚析构函数
     * @details 确保通过基类指针删除派生类对象时正确调用派生类析构函数
     */
    virtual ~EMLinearSolverBase() = default;

    /**
     * @brief 设置系数矩阵A
     * @param A CSR格式稀疏系数矩阵，必须是方阵且已构建完成
     *
     * @details 此方法实现瞬态复用机制：
     * - 当矩阵结构或数值发生变化时调用此方法更新内部表示
     * - 对于直接求解器，触发矩阵分解（LU/Cholesky等），分解结果被缓存
     * - 对于迭代求解器，可能预计算预条件子或提取对角线元素
     * - 后续多次solve()调用可复用缓存的分解结果，显著提升瞬态求解效率
     *
     * @note 调用前确保矩阵A已正确构建（is_built() == true）
     * @note 对于对称正定矩阵，建议使用Cholesky分解以获得更好性能
     * @exception 若矩阵为空、非方阵或未构建，应设置status为INVALID_INPUT并记录错误日志
     *
     * @par 典型使用流程（瞬态分析）：
     * @code
     * // 时间步循环中，仅当矩阵变化时重新设置
     * if (matrix_changed) {
     *     solver->set_matrix(A_new);  // 触发分解并缓存
     * }
     * // 多次求解复用同一分解结果
     * for (auto& rhs : rhs_list) {
     *     auto result = solver->solve(rhs);
     * }
     * @endcode
     */
    virtual void set_matrix(const CsrMatrix<double>& A) = 0;

    /**
     * @brief 求解线性系统 Ax = b
     * @param b 右端项向量（载荷向量），维度必须与矩阵行数一致
     * @return SolverResult 求解结果结构体，包含解向量和状态信息
     *
     * @details 执行实际的线性求解操作：
     * - 直接求解器：利用缓存的分解结果进行前代/回代，时间复杂度O(n²)
     * - 迭代求解器：执行选代算法直到收敛或达到最大迭代次数
     * - 自动记录迭代次数、残差范数、求解耗时等性能指标
     *
     * @note 必须先调用set_matrix()设置系数矩阵，否则行为未定义
     * @note 返回的result.x仅在status == SUCCESS时有效
     * @note 线程安全：此方法不是线程安全的，多线程并发调用需要外部加锁
     *
     * @par 错误处理策略：
     * - 矩阵未设置：返回INVALID_INPUT
     * - 维度不匹配：返回INVALID_INPUT
     * - 分解失败：返回NUMERICAL_ERROR
     * - 迭代发散：返回DIVERGED
     * - 超过最大迭代：返回MAX_ITER_REACHED（返回当前最佳近似解）
     */
    virtual SolverResult solve(const Eigen::VectorXd& b) = 0;

    /**
     * @brief 获取求解器名称标识
     * @return std::string 求解器的唯一名称字符串（如"Eigen_LU"、"GMRES_ILUT"等）
     *
     * @details 用于日志输出、性能统计和调试诊断。
     * 名称应能清晰标识求解器类型和关键参数配置。
     *
     * @note 建议命名格式："后端类型_算法_预处理子"
     * @note 示例返回值："Eigen_SimplicialLDLT"、"GMRES_ILU(0)"、"MUMPS_LU"
     */
    virtual std::string get_solver_name() const = 0;

    /**
     * @brief 清理求解器资源并重置状态
     * @details 释放所有内部缓存的资源，将求解器恢复到初始未初始化状态：
     * - 释放矩阵分解结果（L/U因子、符号分析信息等）
     * - 释放预条件子数据结构
     * - 重置临时工作空间
     * - 清除所有内部状态标志
     *
     * @note 调用后必须重新执行set_matrix()才能进行后续求解
     * @note 在析构函数中自动调用，也可手动调用以提前释放内存
     * @note 线程安全：此方法不是线程安全的
     *
     * @par 典型应用场景：
     * @code
     * // 场景1：显式释放内存
     * solver->clear();  // 释放当前矩阵相关资源
     * // ... 其他操作 ...
     * solver->set_matrix(A2);  // 可安全地设置新矩阵
     *
     * // 场景2：对象生命周期结束时的RAII清理
     * {
     *     auto solver = std::make_unique<ConcreteSolver>();
     *     solver->set_matrix(A);
     *     auto result = solver->solve(b);
     * }  // 析构时自动调用clear()
     * @endcode
     */
    virtual void clear() = 0;
};

} // namespace numeric
