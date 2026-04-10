/**
 * @file em_solver_factory.h
 * @brief 求解器工厂类 - 根据物理场景自动创建最优求解器实例
 * @details 提供工厂模式工具类，封装求解器的创建逻辑，简化上层调用代码。
 *          支持两种使用方式：
 *          1. 直接指定求解器类型枚举（create_solver方法）
 *          2. 根据物理场景自动选择最优求解器（create_solver_for_* 方法）
 *
 * @par 设计目标：
 * - 解耦上层业务与具体求解器实现的耦合
 * - 提供统一的求解器创建入口点
 * - 根据矩阵属性智能选择最优求解算法
 * - 支持大规模问题的迭代/直接法自适应切换
 *
 * @par 可创建的求解器类型映射表：
 * | SolverType枚举值        | 具体求解器类              | 适用场景                          |
 * |--------------------------|---------------------------|-----------------------------------|
 * | SYMMETRIC_DIRECT         | SymmetricDirectSolver     | 对称正定矩阵（Cholesky分解）       |
 * | SYMMETRIC_INDEFINITE     | SymmetricIndefiniteDirectSolver | 对称不定矩阵（LDL^T分解）    |
 * | GENERAL_DIRECT           | GeneralDirectSolver       | 非对称通用矩阵（LU分解）          |
 * | CG                       | CGSolver                  | 大规模SPD矩阵（共轭梯度法）      |
 * | BICGSTAB                 | BiCGSTABSolver            | 大规模非对称矩阵（双共轭梯度法） |
 *
 * @par 场景适配决策树：
 * @code
 * 物理场景 → 矩阵属性分析 → 最优求解器选择
 * ├─ 静电场（标量电位）→ SPD对称矩阵 → SymmetricDirectSolver / CGSolver
 * ├─ 静磁场（标量磁位）→ SPSD对称矩阵 → SymmetricIndefiniteDirectSolver
 * ├─ 静磁场（棱边元）→ 半正定奇异矩阵 → SymmetricIndefiniteDirectSolver
 * └─ 涡流场（A-V格式）→ 不定/非对称矩阵 → GeneralDirectSolver / BiCGSTABSolver
 * @endcode
 *
 * @par 使用示例：
 * @code
 * // 方式1：直接指定求解器类型
 * auto solver = numeric::EMSolverFactory::create_solver(
 *     numeric::EMSolverFactory::SolverType::SYMMETRIC_DIRECT);
 * solver->set_matrix(stiffness_matrix);
 * auto result = solver->solve(rhs_vector);
 *
 * // 方式2：基于物理场景自动选择
 * auto solver = numeric::EMSolverFactory::create_solver_for_electrostatic();
 * solver->set_matrix(electrostatic_matrix);
 * auto result = solver->solve(charge_density);
 *
 * // 方式3：基于矩阵属性智能选择
 * numeric::MatrixAttribute attr(
 *     numeric::MatrixSymmetry::SYMMETRIC,
 *     numeric::MatrixDefiniteness::POSITIVE_DEFINITE,
 *     numeric::MatrixDataType::REAL,
 *     numeric::MatrixElementType::SCALAR,
 *     numeric::PhysicalFieldType::ELECTROSTATIC
 * );
 * auto solver = numeric::EMSolverFactory::create_solver_for_attribute(attr);
 * @endcode
 *
 * @see EMLinearSolverBase 求解器抽象基类
 * @see SymmetricDirectSolver 对称正定直接求解器
 * @see SymmetricIndefiniteDirectSolver 对称不定直接求解器
 * @see GeneralDirectSolver 通用直接求解器
 * @see CGSolver 共轭梯度迭代求解器
 * @see BiCGSTABSolver 双共轭梯度迭代求解器
 * @see MatrixAttribute 矩阵属性标记结构体
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#pragma once

#include "em_linear_solver.h"
#include "em_direct_solvers.h"
#include "em_iterative_solvers.h"
#include "matrix_attribute.hpp"
#include "logger_factory.hpp"
#include <memory>
#include <string>
#include <stdexcept>

namespace numeric {

/**
 * @class EMSolverFactory
 * @brief 电磁场线性系统求解器工厂类
 * @details 采用静态工厂模式，提供统一的求解器创建接口。
 *          所有方法均为静态方法，禁止实例化。
 *          通过智能指针（unique_ptr）管理求解器生命周期，确保无内存泄漏。
 *
 * 本工厂类负责创建和管理各种类型的线性系统求解器实例。
 * 支持的求解器类型包括：
 * - 直接求解器：对称正定(LLT)、对称不定(LDLT)、通用非对称(LU)
 * - 迭代求解器：共轭梯度(CG)、稳定双共轭梯度(BiCGSTAB)、代数多重网格(AMG)
 * - 外部求解器：MUMPS、SuperLU
 *
 * **实数/复数统一架构支持：**
 * 工厂创建的所有求解器实例均同时支持实数矩阵和复数矩阵求解。
 * - 实数场景：使用 set_matrix(CsrMatrix<double>&) 和 solve(Eigen::VectorXd&)
 * - 复数场景：使用 set_matrix(CsrMatrix<std::complex<double>>&) 和 solve(Eigen::VectorXcd&)
 * - 两种模式可独立使用，互不干扰
 *
 * 使用示例：
 * @code
 * // 创建求解器实例
 * auto solver = EMSolverFactory::create(SolverType::SYMMETRIC_DIRECT);
 *
 * // 实数场景使用
 * CsrMatrix<double> K_real = ...;
 * Eigen::VectorXd F_real = ...;
 * solver->set_matrix(K_real);
 * auto result_real = solver->solve(F_real);
 *
 * // 复数场景使用（同一求解器实例）
 * CsrMatrix<std::complex<double>> K_complex = ...;
 * Eigen::VectorXcd F_complex = ...;
 * solver->set_matrix(K_complex);  // 自动切换到复数模式
 * auto result_complex = solver->solve(F_complex);
 * @endcode
 *
 * @note 线程安全：此类的所有方法都是线程安全的（无共享状态）
 * @note 内存安全：返回的unique_ptr自动管理对象生命周期，无需手动delete
 */
class EMSolverFactory {
public:
    /**
     * @enum SolverType
     * @brief 求解器类型枚举
     * @details 定义所有可由工厂创建的求解器类型，每种类型对应一个具体的求解器实现类
     */
    enum class SolverType {
        SYMMETRIC_DIRECT,       ///< 对称正定直接求解器（Eigen::SimplicialLLT Cholesky分解）
                                ///< 适用场景：静电场刚度矩阵、静磁场SPD系统、热传导方程
                                ///< 数值特性：无条件稳定，时间复杂度O(n³)分解+O(n²)求解
                                ///< 内存占用：存储L因子，约为原始矩阵的1.5-2倍

        SYMMETRIC_INDEFINITE,   ///< 对称不定直接求解器（Eigen::SimplicialLDLT LDL^T分解）
                                ///< 适用场景：Saddle Point问题、混合有限元、半正定奇异矩阵
                                ///< 数值特性：支持2x2主元旋转，可处理接近奇异的矩阵
                                ///< 特殊功能：内置正则化机制（添加εI扰动）

        GENERAL_DIRECT,         ///< 通用非对称直接求解器（Eigen::SparseLU LU分解）
                                ///< 适用场景：非对称矩阵系统、复数矩阵、无法保证对称性的通用系统
                                ///< 数值特性：带部分选主元PA=LU，数值稳定性好
                                ///< 性能考虑：fill-in通常高于Cholesky/LDL^T

        CG,                     ///< 共轭梯度迭代求解器（CGSolver）
                                ///< 适用场景：大规模SPD矩阵（百万DOF级别）、条件数适中的良态系统
                                ///< 收敛特性：残差单调递减，理论上n步精确收敛
                                ///< 时间复杂度：O(k*nnz)，k为迭代次数（通常k<<n）
                                ///< 预条件子支持：Jacobi、ILU(0)

        BICGSTAB                ///< 稳定双共轭梯度迭代求解器（BiCGSTABSolver）
                                ///< 适用场景：大规模非对称/非正定矩阵、涡流场A-V formulation
                                ///< 收敛特性：比标准BiCG更稳定，残差可能振荡
                                ///< 时间复杂度：O(k*2*nnz)，每次迭代两次SpMV
                                ///< 预条件子支持：Jacobi、ILU(0)
    };

    /**
     * @brief 根据类型枚举创建求解器实例
     * @param type 求解器类型枚举值（SolverType）
     * @return std::unique_ptr<EMLinearSolverBase> 指向具体求解器实例的智能指针
     *
     * @details 工厂核心方法，通过switch-case分发到对应的求解器构造函数：
     * - 直接求解器（前3种）：使用默认构造函数，可通过返回指针调用set_backend()切换后端
     * - 迭代求解器（后2种）：使用默认配置构造，可通过get_config()/set_config()调整参数
     *
     * @note 返回的指针可直接调用基类统一接口：
     *       - set_matrix(A): 设置系数矩阵并触发分解/预计算
     *       - solve(b): 执行线性求解
     *       - get_solver_name(): 获取求解器标识字符串
     *       - clear(): 释放资源并重置状态
     *
     * @note 使用std::unique_ptr管理生命周期，无需手动delete，避免内存泄漏
     *
     * @exception std::invalid_argument 当type为未知枚举值时抛出
     *
     * @par 典型用法：
     * @code
     * try {
     *     auto solver = EMSolverFactory::create_solver(SolverType::CG);
     *     solver->set_matrix(spd_matrix);
     *     auto result = solver->solve(rhs);
     *     if (result.status == SolverStatus::SUCCESS) {
     *         FEEM_INFO("求解成功, 迭代{}次", result.iterations);
     *     }
     * } catch (const std::invalid_argument& e) {
     *     FEEM_ERROR("求解器创建失败: {}", e.what());
     * }
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver(SolverType type) {
        FEEM_DEBUG("EMSolverFactory::create_solver 创建求解器类型: {}",
                   static_cast<int>(type));

        switch (type) {
            case SolverType::SYMMETRIC_DIRECT:
                FEEM_DEBUG("  -> 创建 SymmetricDirectSolver (Cholesky分解)");
                return std::make_unique<SymmetricDirectSolver>();

            case SolverType::SYMMETRIC_INDEFINITE:
                FEEM_DEBUG("  -> 创建 SymmetricIndefiniteDirectSolver (LDL^T分解)");
                return std::make_unique<SymmetricIndefiniteDirectSolver>();

            case SolverType::GENERAL_DIRECT:
                FEEM_DEBUG("  -> 创建 GeneralDirectSolver (LU分解)");
                return std::make_unique<GeneralDirectSolver>();

            case SolverType::CG: {
                FEEM_DEBUG("  -> 创建 CGSolver (共轭梯度法)");
                CGConfig config;
                return std::make_unique<CGSolver>(config);
            }

            case SolverType::BICGSTAB: {
                FEEM_DEBUG("  -> 创建 BiCGSTABSolver (双共轭梯度法)");
                BiCGSTABConfig config;
                return std::make_unique<BiCGSTABSolver>(config);
            }

            default:
                FEEM_ERROR("未知的求解器类型: {}", static_cast<int>(type));
                throw std::invalid_argument("Unknown solver type in EMSolverFactory");
        }
    }

    /**
     * @brief 获取求解器类型的友好名称字符串
     * @param type 求解器类型枚举值
     * @return std::string 求解器的可读名称（如"SymmetricDirectSolver"、"CGSolver"等）
     *
     * @details 用于日志输出、用户界面显示和调试信息。
     *          返回的名称与实际C++类名保持一致，便于代码追踪。
     *
     * @note 此方法不会抛出异常，对于未知类型返回"UnknownSolver"
     */
    static std::string get_solver_name(SolverType type) {
        switch (type) {
            case SolverType::SYMMETRIC_DIRECT:
                return "SymmetricDirectSolver";
            case SolverType::SYMMETRIC_INDEFINITE:
                return "SymmetricIndefiniteDirectSolver";
            case SolverType::GENERAL_DIRECT:
                return "GeneralDirectSolver";
            case SolverType::CG:
                return "CGSolver";
            case SolverType::BICGSTAB:
                return "BiCGSTABSolver";
            default:
                return "UnknownSolver";
        }
    }

    // ==================== 场景自动适配方法 ====================

    /**
     * @brief 为静电场场景创建最优求解器
     * @return std::unique_ptr<EMLinearSolverBase> 配置好的求解器实例
     *
     * @details 自动选择 SymmetricDirectSolver（Cholesky分解），
     *          因为静电场有限元离散后的刚度矩阵通常是对称正定(SP)D矩阵。
     *
     * @par 物理背景：
     * 静电场控制方程为 Poisson/Laplace 方程：∇·(ε∇φ) = -ρ
     * - 使用标量拉格朗日元离散化
     * - 刚度矩阵 K_ij = ∫ ε ∇N_i · ∇N_j dΩ （对称正定）
     * - 适用于导体电容计算、电场分布仿真等
     *
     * @par 选择理由：
     * - Cholesky分解对SPD矩阵具有最优数值稳定性
     * - 分解结果可复用（瞬态分析中多次右端项求解）
     * - 对于中小规模问题（<10万DOF），直接法效率高于迭代法
     *
     * @note 若问题规模特别大（>100万DOF），建议手动选择CGSolver
     *
     * @code
     * auto solver = EMSolverFactory::create_solver_for_electrostatic();
     * solver->set_matrix(electrostatic_stiffness);
     * auto potential = solver->solve(charge_density).x;
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver_for_electrostatic() {
        FEEM_INFO("EMSolverFactory: 为静电场场景创建求解器");
        FEEM_DEBUG("  场景特征: SPD对称矩阵, 标量电位, Laplace/Poisson方程");
        FEEM_DEBUG("  选择: SymmetricDirectSolver (Cholesky分解)");

        return create_solver(SolverType::SYMMETRIC_DIRECT);
    }

    /**
     * @brief 为静磁场标量位场景创建最优求解器
     * @param prefer_iterative 是否优先使用迭代求解器（默认false，优先直接法）
     * @return std::unique_ptr<EMLinearSolverBase> 配置好的求解器实例
     *
     * @details 根据问题规模和用户偏好自动选择：
     * - prefer_iterative=false（默认）：选择 SymmetricDirectSolver 或 SymmetricIndefiniteDirectSolver
     * - prefer_iterative=true：选择 CGSolver（适合大规模问题）
     *
     * @par 物理背景：
     * 静磁场标量磁位 ψ 的控制方程：∇·(μ⁻¹∇ψ) = 0
     * - 使用标量拉格朗日元离散化
     * - 刚度矩阵通常为对称半正定(SPSD)，可能存在零空间（常数模式）
     * - 适用于无源静磁场、磁屏蔽效应等
     *
     * @par 决策逻辑：
     * @code
     * if (prefer_iterative) {
     *     // 大规模问题或内存受限环境
     *     return CGSolver;  // O(n)内存，适合百万DOF
     * } else {
     *     // 中小规模问题，追求精度和稳定性
     *     return SymmetricIndefiniteDirectSolver;  // 处理可能的半正定性
     * }
     * @endcode
     *
     * @param prefer_iterative true表示优先迭代法（适合大规模问题），false优先直接法
     *
     * @code
     * // 中小规模静磁场问题（推荐）
     * auto solver = EMSolverFactory::create_solver_for_magnetostatic_scalar();
     *
     * // 大规模静磁场问题（>50万DOF）
     * auto solver = EMSolverFactory::create_solver_for_magnetostatic_scalar(true);
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver_for_magnetostatic_scalar(bool prefer_iterative = false) {
        FEEM_INFO("EMSolverFactory: 为静磁场标量位场景创建求解器");
        FEEM_DEBUG("  场景特征: SPSD对称矩阵, 标量磁位, 可能存在零空间");

        if (prefer_iterative) {
            FEEM_DEBUG("  用户偏好: 迭代法 → 选择 CGSolver (适合大规模问题)");
            return create_solver(SolverType::CG);
        } else {
            FEEM_DEBUG("  默认策略: 直接法 → 选择 SymmetricIndefiniteDirectSolver (处理半正定性)");
            return create_solver(SolverType::SYMMETRIC_INDEFINITE);
        }
    }

    /**
     * @brief 为静磁场矢量位（棱边单元）场景创建最优求解器
     * @return std::unique_ptr<EMLinearSolverBase> 配置好的求解器实例
     *
     * @details 自动选择 SymmetricIndefiniteDirectSolver（LDL^T分解），
     *          因为此场景下的刚度矩阵通常是半正定甚至奇异的。
     *
     * @par 物理背景：
     * 静磁场矢量磁位 A 的控制方程：∇×(μ⁻¹∇×A) = J
     * - 使用Nédélec棱边单元（H(curl)协调空间）离散化
     * - 刚度矩阵为对称半正定，存在非平凡null space（梯度场模式）
     * - 必须施加规范条件（gauge condition）消除奇异性的影响
     *
     * @par 矩阵特性：
     * - 对称性：对称 ✓
     * - 正定性：半正定（存在零特征值对应梯度场）
     * - 奇异性：可能奇异（取决于边界条件和规范处理）
     * - Null space维度：至少等于连通区域数
     *
     * @par 选择理由：
     * LDL^T分解相比Cholesky更能处理半正定矩阵：
     * - 允许D因子的对角元为零或负值
     * - 内置正则化机制可处理近奇异情况
     * - 2x2主元旋转提升数值稳定性
     *
     * @warning 使用前确保已正确施加树规范（tree gauge）或其他规范条件
     *
     * @code
     * auto solver = EMSolverFactory::create_solver_for_magnetostatic_edge();
     * solver->set_matrix(curl_curl_stiffness);  // ∇×(μ⁻¹∇×A) 矩阵
     * auto A = solver->solve(J_source).x;       // 磁矢势分布
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver_for_magnetostatic_edge() {
        FEEM_INFO("EMSolverFactory: 为静磁场矢量位（棱边元）场景创建求解器");
        FEEM_DEBUG("  场景特征: SPSD对称矩阵, H(curl)空间, 存在梯度场null space");
        FEEM_DEBUG("  可能奇异 → 需要 SymmetricIndefiniteDirectSolver (带正则化的LDL^T)");

        return create_solver(SolverType::SYMMETRIC_INDEFINITE);
    }

    /**
     * @brief 为涡流场A-V混合格式创建最优求解器
     * @param prefer_iterative 是否优先使用迭代求解器（默认false）
     * @return std::unique_ptr<EMLinearSolverBase> 配置好的求解器实例
     *
     * @details 涡流场的A-V格式产生非对称或不定的块矩阵系统，
     *          根据规模自动选择 GeneralDirectSolver 或 BiCGSTABSolver。
     *
     * @par 物理背景：
     * 时谐涡流场控制方程（A-V formulation）：
     * - ∇×(ν∇×A) + jωσA + σ∇V = J_s  （旋度旋度方程）
     * - ∇·(jωσA + σ∇V) = 0           （电荷守恒）
     *
     * 离散后形成Saddle Point结构的块矩阵系统：
     * @code
     * [K_AA   K_AV] [A]   [J_s]
     * [K_VA   0   ] [V] = [ 0 ]
     * @endcode
     *
     * @par 矩阵特性：
     * - 对称性：复对称（Hermitian）或非对称（取决于公式变体）
     * - 正定性：不定（Saddle Point结构导致）
     * - 块结构：左下角为零块的典型KKT系统
     * - 条件数：通常较大，预条件子至关重要
     *
     * @par 决策逻辑：
     * @code
     * if (prefer_iterative) {
     *     // 大规模涡流问题（3D复杂几何，数十万棱边未知量）
     *     return BiCGSTABSolver;  // 支持非对称矩阵，配合ILU预条件子
     * } else {
     *     // 中小规模或高精度要求
     *     return GeneralDirectSolver;  // LU分解保证稳定性
     * }
     * @endcode
     *
     * @param prefer_iterative true优先BiCGSTAB迭代法，false优先LU直接法
     *
     * @code
     * // 中小规模涡流问题
     * auto solver = EMSolverFactory::create_solver_for_eddy_current();
     *
     * // 大规模3D涡流仿真（推荐迭代法+ILU预条件子）
     * auto solver = EMSolverFactory::create_solver_for_eddy_current(true);
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver_for_eddy_current(bool prefer_iterative = false) {
        FEEM_INFO("EMSolverFactory: 为涡流场A-V格式场景创建求解器");
        FEEM_DEBUG("  场景特征: 不定/非对称块矩阵, Saddle Point结构, 复数系统");

        if (prefer_iterative) {
            FEEM_DEBUG("  用户偏好: 迭代法 → 选择 BiCGSTABSolver (支持非对称矩阵)");
            return create_solver(SolverType::BICGSTAB);
        } else {
            FEEM_DEBUG("  默认策略: 直接法 → 选择 GeneralDirectSolver (LU分解，稳定性最优)");
            return create_solver(SolverType::GENERAL_DIRECT);
        }
    }

    /**
     * @brief 基于MatrixAttribute属性自动选择最优求解器
     * @param attr 矩阵属性标记（包含对称性、正定性、元素类型、物理场类型等完整信息）
     * @return std::unique_ptr<EMLinearSolverBase> 根据属性智能选择的求解器实例
     *
     * @details 核心智能决策方法，通过分析矩阵的多维属性自动选择最合适的求解器。
     *          决策树综合考虑以下因素：
     *
     * @par 决策逻辑（优先级从高到低）：
     *
     * **1. 正定性判断（最高优先级）**
     * @code
     * if (attr.is_spd) {
     *     // 对称正定矩阵：Cholesky或CG均可
     *     // 小规模 → Cholesky（精度高，稳定性好）
     *     // 大规模 → CG（内存省，可扩展性强）
     * }
     * @endcode
     *
     * **2. 奇异性/半正定性检测**
     * @code
     * if (attr.is_singular || definiteness == POSITIVE_SEMIDEFINITE) {
     *     // 接近奇异或半正定：必须使用LDL^T（支持零主元）
     *     return SymmetricIndefiniteDirectSolver;
     * }
     * @endcode
     *
     * **3. 对称性检查**
     * @code
     * if (symmetry == UNSYMMETRIC) {
     *     // 非对称矩阵：只能用LU或BiCGSTAB
     *     // 小规模 → LU（稳定但慢）
     *     // 大规模 → BiCGSTAB（快但可能不收敛）
     * }
     * @endcode
     *
     * **4. 物理场类型辅助决策**
     * @code
     * switch (field_type) {
     *     case ELECTROSTATIC:  // 通常SPD → SymmetricDirectSolver
     *     case MAGNETOSTATIC:  // 通常SPSD → SymmetricIndefiniteDirectSolver
     *     case EDDY_CURRENT:   // 通常不定/非对称 → GeneralDirect/BiCGSTAB
     * }
     * @endcode
     *
     * **5. 规模自适应（未来扩展）**
     * 当前版本暂不支持基于矩阵行数的规模判断，
     * 未来可通过新增参数 rows > threshold 来触发迭代法优选。
     *
     * @par 属性组合示例：
     * | symmetry | definiteness | is_singular | field_type   | 推荐求解器                    |
     * |----------|--------------|-------------|--------------|-------------------------------|
     * | SYMMETRIC| POS_DEF      | false       | ELECTROSTATIC| SymmetricDirectSolver         |
     * | SYMMETRIC| POS_SEMIDEF  | true        | MAGNETOSTATIC| SymmetricIndefiniteDirectSolver|
     * | HERMITIAN| POS_DEF      | false       | EDDY_CURRENT | GeneralDirectSolver           |
     * | UNSYMMETRIC| INDEFINITE | false       | EDDY_CURRENT | GeneralDirectSolver           |
     * | SYMMETRIC| POS_DEF      | false       | MAGNETOSTATIC| SymmetricDirectSolver (or CG) |
     *
     * @note 此方法是场景适配方法的底层实现，其他create_solver_for_*方法最终都调用此方法或其简化版
     * @note 日志输出详细的决策过程，便于调试和性能优化分析
     *
     * @code
     * // 示例1：静电场矩阵属性
     * numeric::MatrixAttribute electrostatic_attr =
     *     numeric::MatrixAttribute::create_electrostatic();
     * auto solver = EMSolverFactory::create_solver_for_attribute(electrostatic_attr);
     * // 输出日志: "SPD对称矩阵 → SymmetricDirectSolver (Cholesky)"
     *
     * // 示例2：自定义属性（静磁场棱边元）
     * numeric::MatrixAttribute edge_attr(
     *     numeric::MatrixSymmetry::SYMMETRIC,
     *     numeric::MatrixDefiniteness::POSITIVE_SEMIDEFINITE,
     *     numeric::MatrixDataType::REAL,
     *     numeric::MatrixElementType::VECTOR_3D,
     *     numeric::PhysicalFieldType::MAGNETOSTATIC,
     *     true  // is_singular=true
     * );
     * auto solver = EMSolverFactory::create_solver_for_attribute(edge_attr);
     * // 输出日志: "半正定奇异矩阵 → SymmetricIndefiniteDirectSolver (LDL^T+正则化)"
     * @endcode
     */
    static std::unique_ptr<EMLinearSolverBase> create_solver_for_attribute(const MatrixAttribute& attr) {
        FEEM_INFO("EMSolverFactory: 基于矩阵属性自动选择求解器");
        FEEM_DEBUG("  完整属性: {}", attr.to_string());

        // ========== 第一优先级：对称正定矩阵 ==========
        if (attr.is_spd) {
            FEEM_DEBUG("  检测到: SPD对称正定矩阵");

            // 判断是否适合使用CG迭代法
            bool suitable_for_cg = attr.suitable_for_cg();

            if (suitable_for_cg) {
                // TODO: 未来可根据矩阵规模（rows > threshold）自动选择迭代法
                // 当前默认使用直接法以保证精度和稳定性
                FEEM_DEBUG("  适合CG但默认选择直接法（更高精度）");
                FEEM_DEBUG("  最终选择: SymmetricDirectSolver (Cholesky分解)");
                return create_solver(SolverType::SYMMETRIC_DIRECT);
            } else {
                // 复数SPD矩阵不适合标准CG
                FEEM_DEBUG("  复数SPD矩阵 → SymmetricDirectSolver (Cholesky支持复数)");
                return create_solver(SolverType::SYMMETRIC_DIRECT);
            }
        }

        // ========== 第二优先级：奇异/半正定矩阵 ==========
        if (attr.is_singular ||
            attr.definiteness == MatrixDefiniteness::POSITIVE_SEMIDEFINITE) {
            FEEM_DEBUG("  检测到: 半正定奇异矩阵（is_singular={}, definiteness={})",
                      attr.is_singular ? "是" : "否",
                      attr.definiteness_string());
            FEEM_DEBUG("  最终选择: SymmetricIndefiniteDirectSolver (LDL^T分解+正则化)");
            FEEM_DEBUG("  说明: LDL^T允许零主元，内置正则化机制处理近奇异情况");
            return create_solver(SolverType::SYMMETRIC_INDEFINITE);
        }

        // ========== 第三优先级：非对称矩阵 ==========
        if (attr.symmetry == MatrixSymmetry::UNSYMMETRIC) {
            FEEM_DEBUG("  检测到: 非对称矩阵");

            // 非对称矩阵只能用通用求解器
            // TODO: 未来可根据矩阵规模选择BiCGSTAB（大规模）或LU（小规模）
            FEEM_DEBUG("  最终选择: GeneralDirectSolver (LU分解，部分选主元)");
            FEEM_DEBUG("  替代方案: 大规模问题可考虑 BiCGSTABSolver");
            return create_solver(SolverType::GENERAL_DIRECT);
        }

        // ========== 第四优先级：埃尔米特（复对称）矩阵 ==========
        if (attr.symmetry == MatrixSymmetry::HERMITIAN) {
            FEEM_DEBUG("  检测到: 埃尔米特矩阵（复对称，常见于时谐场）");

            // 复对称矩阵的处理取决于正定性
            if (attr.definiteness == MatrixDefiniteness::POSITIVE_DEFINITE) {
                FEEM_DEBUG("  复正定 → GeneralDirectSolver (LU分解支持复数)");
            } else {
                FEEM_DEBUG("  复不定 → GeneralDirectSolver (LU分解最通用)");
            }
            return create_solver(SolverType::GENERAL_DIRECT);
        }

        // ========== 第五优先级：对称但不定矩阵 ==========
        if (attr.symmetry == MatrixSymmetry::SYMMETRIC &&
            attr.definiteness == MatrixDefiniteness::INDEFINITE) {
            FEEM_DEBUG("  检测到: 对称不定矩阵（Saddle Point结构？）");
            FEEM_DEBUG("  最终选择: SymmetricIndefiniteDirectSolver (LDL^T分解)");
            FEEM_DEBUG("  说明: 相比Cholesky，LDL^T能处理负特征值和零主元");
            return create_solver(SolverType::SYMMETRIC_INDEFINITE);
        }

        // ========== 默认回退方案 ==========
        FEEM_WARN("  无法精确匹配属性，使用默认通用求解器");
        FEEM_DEBUG("  默认选择: GeneralDirectSolver (最通用的LU分解)");
        FEEM_DEBUG("  建议: 检查MatrixAttribute设置是否正确，或手动指定求解器类型");
        return create_solver(SolverType::GENERAL_DIRECT);
    }

    /**
     * @brief 根据矩阵数据类型判断是否需要使用复数求解模式
     * @param data_type 矩阵数据类型（REAL 或 COMPLEX）
     * @return true 如果需要使用复数模式（data_type == COMPLEX），false 如果使用实数模式
     *
     * @details 此方法可用于上层代码根据问题类型自动选择正确的接口调用方式。
     *          在电磁场仿真中，不同物理场景可能需要不同的数值表示：
     *
     * @par 典型应用场景：
     * - **实数模式（MatrixDataType::REAL）**：
     *   静电场、静磁场等静态/准静态问题
     *   系统矩阵为实对称或实非对称
     *
     * - **复数模式（MatrixDataType::COMPLEX）**：
     *   时谐涡流场、波导问题、高频电磁仿真
     *   系统矩阵包含虚部（由 jωσ 或 k₀²ε 项引入）
     *
     * @par 使用示例：
     * @code
     * // 根据问题属性自动选择接口
     * numeric::MatrixAttribute attr = ...;  // 从配置或分析获得
     * bool use_complex = EMSolverFactory::is_complex_mode(attr.data_type);
     *
     * auto solver = EMSolverFactory::create_solver(SolverType::GENERAL_DIRECT);
     *
     * if (use_complex) {
     *     // 复数求解路径
     *     CsrMatrix<std::complex<double>> K_complex = build_complex_matrix();
     *     Eigen::VectorXcd F_complex = build_complex_rhs();
     *     solver->set_matrix(K_complex);
     *     auto result = solver->solve(F_complex);
     * } else {
     *     // 实数求解路径
     *     CsrMatrix<double> K_real = build_real_matrix();
     *     Eigen::VectorXd F_real = build_real_rhs();
     *     solver->set_matrix(K_real);
     *     auto result = solver->solve(F_real);
     * }
     * @endcode
     *
     * @note 此方法是纯逻辑判断工具，不涉及任何计算或状态修改
     * @note 可与 create_solver_for_attribute() 配合使用，实现完全自动化的求解流程
     */
    static bool is_complex_mode(MatrixDataType data_type) {
        return data_type == MatrixDataType::COMPLEX;
    }

private:
    EMSolverFactory() = delete;  // 禁止实例化（纯静态工厂类）
};

} // namespace numeric
