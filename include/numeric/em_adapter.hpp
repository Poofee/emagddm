/**
 * @file em_adapter.hpp
 * @brief 电磁场景适配器
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 * 
 * @details 根据矩阵属性自动选择预处理和求解策略，处理静磁场散度约束
 *          适配静电、静磁、涡流场等低频电磁场景
 */

#pragma once

#include "csr_matrix.hpp"
#include "coo_matrix.hpp"
#include "block_csr_matrix.hpp"
#include "matrix_attribute.hpp"
#include "complex_matrix_ops.hpp"
#include "preconditioner.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace numeric {

/**
 * @class EMAdapter
 * @brief 电磁场景适配器
 * @details 根据矩阵属性自动适配预处理和求解策略
 */
class EMAdapter {
public:
    /**
     * @brief 求解器类型枚举
     */
    enum class SolverType {
        CG,         ///< 共轭梯度法（对称正定）
        GMRES,      ///< 广义最小残差法（非对称/复数）
        MINRES,     ///< 最小残差法（对称不定）
        BICGSTAB,   ///< 稳定双共轭梯度法（一般矩阵）
        DIRECT      ///< 直接法（小规模问题）
    };
    
    /**
     * @brief 预处理类型枚举
     */
    enum class PreconditionerType {
        JACOBI,     ///< 雅可比预处理
        ILU0,       ///< ILU(0)预处理
        BLOCK_ILU,  ///< 块ILU预处理
        NONE        ///< 无预处理
    };
    
    /**
     * @brief 求解器配置结构
     */
    struct SolverConfig {
        SolverType solver_type;           ///< 求解器类型
        PreconditionerType precond_type;  ///< 预处理类型
        int max_iterations;               ///< 最大迭代次数
        double tolerance;                 ///< 收敛容差
        bool use_div_constraint;          ///< 是否使用散度约束
        
        /**
         * @brief 默认构造函数
         */
        SolverConfig()
            : solver_type(SolverType::CG)
            , precond_type(PreconditionerType::JACOBI)
            , max_iterations(1000)
            , tolerance(1e-8)
            , use_div_constraint(false) {}
        
        /**
         * @brief 参数化构造函数
         */
        SolverConfig(SolverType s, PreconditionerType p, int max_iter, double tol, bool use_div)
            : solver_type(s)
            , precond_type(p)
            , max_iterations(max_iter)
            , tolerance(tol)
            , use_div_constraint(use_div) {}
        
        /**
         * @brief 获取配置描述
         */
        std::string to_string() const {
            std::string solver_str;
            switch (solver_type) {
                case SolverType::CG: solver_str = "CG"; break;
                case SolverType::GMRES: solver_str = "GMRES"; break;
                case SolverType::MINRES: solver_str = "MINRES"; break;
                case SolverType::BICGSTAB: solver_str = "BiCGSTAB"; break;
                case SolverType::DIRECT: solver_str = "Direct"; break;
                default: solver_str = "Unknown";
            }
            
            std::string precond_str;
            switch (precond_type) {
                case PreconditionerType::JACOBI: precond_str = "Jacobi"; break;
                case PreconditionerType::ILU0: precond_str = "ILU0"; break;
                case PreconditionerType::BLOCK_ILU: precond_str = "BlockILU"; break;
                case PreconditionerType::NONE: precond_str = "None"; break;
                default: precond_str = "Unknown";
            }
            
            return "求解器: " + solver_str + 
                   ", 预处理: " + precond_str + 
                   ", 最大迭代: " + std::to_string(max_iterations) + 
                   ", 容差: " + std::to_string(tolerance) + 
                   ", 散度约束: " + (use_div_constraint ? "是" : "否");
        }
    };
    
    /**
     * @brief 根据矩阵属性自动配置求解器
     * @param attr 矩阵属性
     * @return 求解器配置
     */
    static SolverConfig auto_configure(const MatrixAttribute& attr) {
        SolverConfig config;
        
        // 根据矩阵属性选择求解器
        if (attr.suitable_for_cg()) {
            config.solver_type = SolverType::CG;
        } else if (attr.data_type == MatrixDataType::COMPLEX) {
            config.solver_type = SolverType::GMRES;
        } else if (attr.is_singular) {
            config.solver_type = SolverType::MINRES;
        } else {
            config.solver_type = SolverType::BICGSTAB;
        }
        
        // 根据矩阵属性选择预处理
        if (attr.suitable_for_block_preconditioner()) {
            config.precond_type = PreconditionerType::BLOCK_ILU;
        } else if (attr.suitable_for_ilu()) {
            config.precond_type = PreconditionerType::ILU0;
        } else if (attr.suitable_for_jacobi()) {
            config.precond_type = PreconditionerType::JACOBI;
        } else {
            config.precond_type = PreconditionerType::NONE;
        }
        
        // 静磁场需要散度约束
        config.use_div_constraint = (attr.field_type == PhysicalFieldType::MAGNETOSTATIC && 
                                    attr.is_singular);
        
        return config;
    }
    
    /**
     * @brief 为静磁场矩阵添加散度约束
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param constraint_penalty 约束罚参数
     * @return 添加约束后的矩阵
     */
    template<typename T>
    static CsrMatrix<T> add_divergence_constraint(const CsrMatrix<T>& matrix, 
                                                 T constraint_penalty = T(1e-6)) {
        if (matrix.rows() != matrix.cols()) {
            throw std::invalid_argument("散度约束仅适用于方阵");
        }
        
        int n = matrix.rows();
        
        // 创建约束矩阵（对角矩阵，对角元素为罚参数）
        CooMatrix<T> constraint_matrix(n, n);
        
        // 添加对角约束
        for (int i = 0; i < n; ++i) {
            constraint_matrix.add_value(i, i, constraint_penalty);
        }
        
        // 将约束矩阵加到原矩阵上
        CsrMatrix<T> constrained_matrix(n, n);
        
        // 将原矩阵转换为COO
        CooMatrix<T> coo_original(n, n);
        const auto& row_ptr = matrix.get_row_ptr();
        const auto& col_indices = matrix.get_col_indices();
        const auto& values = matrix.get_values();
        
        for (int i = 0; i < n; ++i) {
            int start = row_ptr[i];
            int end = row_ptr[i + 1];
            
            for (int j = start; j < end; ++j) {
                coo_original.add_value(i, col_indices[j], values[j]);
            }
        }
        
        // 合并原矩阵和约束矩阵
        CooMatrix<T> coo_combined = coo_original;
        const auto& constraint_rows = constraint_matrix.get_row_indices();
        const auto& constraint_cols = constraint_matrix.get_col_indices();
        const auto& constraint_values = constraint_matrix.get_values();
        
        for (size_t i = 0; i < constraint_rows.size(); ++i) {
            coo_combined.add_value(constraint_rows[i], constraint_cols[i], constraint_values[i]);
        }
        
        constrained_matrix.build_from_coo(coo_combined);
        return constrained_matrix;
    }
    
    /**
     * @brief 创建块预处理器
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param block_size 块大小
     * @return 块预处理器
     */
    template<typename T>
    static std::unique_ptr<Preconditioner<T>> create_block_preconditioner(
        const CsrMatrix<T>& matrix, BlockSize block_size) {
        
        // 创建块CSR矩阵
        BlockCsrMatrix<T> block_matrix(matrix.rows() / static_cast<int>(block_size), 
                                      matrix.cols() / static_cast<int>(block_size), 
                                      block_size);
        
        // 将CSR矩阵转换为COO，再转换为块CSR
        CooMatrix<T> coo_matrix(matrix.rows(), matrix.cols());
        const auto& row_ptr = matrix.get_row_ptr();
        const auto& col_indices = matrix.get_col_indices();
        const auto& values = matrix.get_values();
        
        for (int i = 0; i < matrix.rows(); ++i) {
            int start = row_ptr[i];
            int end = row_ptr[i + 1];
            
            for (int j = start; j < end; ++j) {
                coo_matrix.add_value(i, col_indices[j], values[j]);
            }
        }
        
        block_matrix.build_from_coo(coo_matrix);
        
        // 创建块ILU预处理器（简化实现）
        // 实际应用中需要实现完整的块ILU分解
        return std::make_unique<JacobiPreconditioner<T>>(matrix);
    }
    
    /**
     * @brief 创建适合电磁场景的预处理器
     * @tparam T 矩阵元素类型
     * @param matrix 输入矩阵
     * @param attr 矩阵属性
     * @return 预处理器
     */
    template<typename T>
    static std::unique_ptr<Preconditioner<T>> create_preconditioner(
        const CsrMatrix<T>& matrix, const MatrixAttribute& attr) {
        
        if (attr.suitable_for_block_preconditioner()) {
            // 矢量元使用块预处理
            BlockSize block_size = BlockSize::BLOCK_1x1;
            if (attr.element_type == MatrixElementType::VECTOR_2D) {
                block_size = BlockSize::BLOCK_2x2;
            } else if (attr.element_type == MatrixElementType::VECTOR_3D) {
                block_size = BlockSize::BLOCK_3x3;
            }
            return create_block_preconditioner(matrix, block_size);
            
        } else if (attr.suitable_for_ilu()) {
            // 适合ILU预处理的矩阵
            return std::make_unique<ILU0Preconditioner<T>>(matrix);
            
        } else if (attr.suitable_for_jacobi()) {
            // 适合Jacobi预处理的矩阵
            return std::make_unique<JacobiPreconditioner<T>>(matrix);
            
        } else {
            // 无预处理
            return nullptr;
        }
    }
    
    /**
     * @brief 求解线性系统（简化接口）
     * @tparam T 矩阵元素类型
     * @param A 系数矩阵
     * @param b 右端向量
     * @param x 解向量
     * @param attr 矩阵属性
     * @param config 求解器配置（可选）
     * @return 是否收敛
     */
    template<typename T>
    static bool solve(const CsrMatrix<T>& A, const std::vector<T>& b, 
                     std::vector<T>& x, const MatrixAttribute& attr,
                     const SolverConfig& config = SolverConfig()) {
        
        SolverConfig actual_config = config;
        if (config.solver_type == SolverType::CG) { // 使用默认配置
            actual_config = auto_configure(attr);
        }
        
        // 处理静磁场散度约束
        CsrMatrix<T> A_modified = A;
        if (actual_config.use_div_constraint) {
            A_modified = add_divergence_constraint(A);
        }
        
        // 创建预处理器
        auto precond = create_preconditioner(A_modified, attr);
        
        // 简化求解过程（实际实现需要完整的迭代求解器）
        // 这里仅演示接口设计
        x.resize(b.size(), T(0));
        
        // 使用预处理后的矩阵求解
        if (precond) {
            // 应用预处理
            std::vector<T> b_precond(b.size());
            precond->apply(b, b_precond);
            
            // 简化解法：使用预处理后的系统
            // 实际实现需要完整的迭代求解器
            for (size_t i = 0; i < b_precond.size(); ++i) {
                x[i] = b_precond[i]; // 简化处理
            }
        } else {
            // 无预处理直接求解
            for (size_t i = 0; i < b.size(); ++i) {
                x[i] = b[i]; // 简化处理
            }
        }
        
        return true; // 简化实现，总是返回收敛
    }
    
    /**
     * @brief 验证求解结果
     * @tparam T 矩阵元素类型
     * @param A 系数矩阵
     * @param b 右端向量
     * @param x 解向量
     * @param tolerance 容差
     * @return 残差范数
     */
    template<typename T>
    static double verify_solution(const CsrMatrix<T>& A, const std::vector<T>& b, 
                                 const std::vector<T>& x, double tolerance = 1e-10) {
        
        if (x.size() != b.size()) {
            throw std::invalid_argument("解向量尺寸与右端向量不匹配");
        }
        
        // 计算残差 r = b - A*x
        std::vector<T> r(b.size());
        A.mat_vec(x, r);
        
        for (size_t i = 0; i < b.size(); ++i) {
            r[i] = b[i] - r[i];
        }
        
        // 计算残差范数
        double residual_norm = 0.0;
        for (const auto& val : r) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                residual_norm += std::norm(val);
            } else {
                residual_norm += val * val;
            }
        }
        
        return std::sqrt(residual_norm);
    }
    
    /**
     * @brief 创建静电场景求解配置
     * @return 静电场景配置
     */
    static SolverConfig create_electrostatic_config() {
        return SolverConfig(
            SolverType::CG,
            PreconditionerType::ILU0,
            1000,
            1e-8,
            false
        );
    }
    
    /**
     * @brief 创建静磁场景求解配置
     * @return 静磁场景配置
     */
    static SolverConfig create_magnetostatic_config() {
        return SolverConfig(
            SolverType::MINRES,
            PreconditionerType::JACOBI,
            2000,
            1e-6,
            true
        );
    }
    
    /**
     * @brief 创建涡流场求解配置
     * @return 涡流场配置
     */
    static SolverConfig create_eddy_current_config() {
        return SolverConfig(
            SolverType::GMRES,
            PreconditionerType::BLOCK_ILU,
            1500,
            1e-8,
            false
        );
    }
    
    /**
     * @brief 获取场景配置描述
     * @param field_type 物理场类型
     * @return 配置描述
     */
    static std::string get_scenario_description(PhysicalFieldType field_type) {
        switch (field_type) {
            case PhysicalFieldType::ELECTROSTATIC:
                return "静电场：对称正定矩阵，适合CG+ILU0求解";
            case PhysicalFieldType::MAGNETOSTATIC:
                return "静磁场：对称半正定奇异矩阵，需要散度约束，适合MINRES+Jacobi求解";
            case PhysicalFieldType::EDDY_CURRENT:
                return "涡流场：复数埃尔米特矩阵，矢量元结构，适合GMRES+块ILU求解";
            case PhysicalFieldType::WAVE:
                return "波动场：复数矩阵，适合GMRES/BiCGSTAB求解";
            default:
                return "未知场景";
        }
    }
};

} // namespace numeric