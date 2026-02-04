/**
 * @file matrix_attribute.hpp
 * @brief 矩阵属性标记结构定义
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 * 
 * @details 定义矩阵的各种属性标记，用于自动适配预处理和操作逻辑
 *          支持对称性、正定性、数据类型、元素类型、物理场类型等属性
 */

#pragma once

#include <string>
#include <cstdint>

namespace numeric {

/**
 * @brief 矩阵对称性枚举
 */
enum class MatrixSymmetry {
    UNSYMMETRIC = 0,    ///< 非对称矩阵
    SYMMETRIC = 1,      ///< 对称矩阵
    HERMITIAN = 2,      ///< 埃尔米特矩阵（复对称）
    SKEW_SYMMETRIC = 3  ///< 反对称矩阵
};

/**
 * @brief 矩阵正定性枚举
 */
enum class MatrixDefiniteness {
    INDEFINITE = 0,     ///< 不定矩阵
    POSITIVE_DEFINITE = 1,  ///< 正定矩阵
    POSITIVE_SEMIDEFINITE = 2, ///< 半正定矩阵
    NEGATIVE_DEFINITE = 3,   ///< 负定矩阵
    NEGATIVE_SEMIDEFINITE = 4 ///< 半负定矩阵
};

/**
 * @brief 矩阵数据类型枚举
 */
enum class MatrixDataType {
    REAL = 0,           ///< 实数矩阵
    COMPLEX = 1         ///< 复数矩阵
};

/**
 * @brief 矩阵元素类型枚举（有限元类型）
 */
enum class MatrixElementType {
    SCALAR = 0,         ///< 标量元（拉格朗日元）
    VECTOR_2D = 1,      ///< 二维矢量元（Nedelec边元）
    VECTOR_3D = 2       ///< 三维矢量元（Nedelec边元）
};

/**
 * @brief 物理场类型枚举
 */
enum class PhysicalFieldType {
    ELECTROSTATIC = 0,  ///< 静电场
    MAGNETOSTATIC = 1,  ///< 静磁场
    EDDY_CURRENT = 2,   ///< 涡流场
    WAVE = 3            ///< 波动场
};

/**
 * @struct MatrixAttribute
 * @brief 矩阵属性标记结构体
 * @details 用于标记矩阵的各种属性，指导预处理和求解策略选择
 */
struct MatrixAttribute {
    MatrixSymmetry symmetry;          ///< 矩阵对称性
    MatrixDefiniteness definiteness;  ///< 矩阵正定性
    MatrixDataType data_type;         ///< 矩阵数据类型（实/复）
    MatrixElementType element_type;   ///< 矩阵元素类型（标量/矢量元）
    PhysicalFieldType field_type;     ///< 物理场类型
    bool is_singular;                 ///< 是否奇异矩阵
    bool is_spd;                      ///< 是否对称正定（快捷标记）
    bool is_hermitian;                ///< 是否埃尔米特矩阵（快捷标记）
    
    /**
     * @brief 默认构造函数
     */
    MatrixAttribute()
        : symmetry(MatrixSymmetry::UNSYMMETRIC)
        , definiteness(MatrixDefiniteness::INDEFINITE)
        , data_type(MatrixDataType::REAL)
        , element_type(MatrixElementType::SCALAR)
        , field_type(PhysicalFieldType::ELECTROSTATIC)
        , is_singular(false)
        , is_spd(false)
        , is_hermitian(false) {}
    
    /**
     * @brief 参数化构造函数
     * @param sym 对称性
     * @param def 正定性
     * @param dtype 数据类型
     * @param etype 元素类型
     * @param ftype 物理场类型
     * @param singular 是否奇异
     */
    MatrixAttribute(MatrixSymmetry sym, MatrixDefiniteness def, 
                   MatrixDataType dtype, MatrixElementType etype,
                   PhysicalFieldType ftype, bool singular = false)
        : symmetry(sym)
        , definiteness(def)
        , data_type(dtype)
        , element_type(etype)
        , field_type(ftype)
        , is_singular(singular)
        , is_spd(sym == MatrixSymmetry::SYMMETRIC && def == MatrixDefiniteness::POSITIVE_DEFINITE)
        , is_hermitian(sym == MatrixSymmetry::HERMITIAN) {}
    
    /**
     * @brief 更新快捷标记
     */
    void update_flags() {
        is_spd = (symmetry == MatrixSymmetry::SYMMETRIC && 
                  definiteness == MatrixDefiniteness::POSITIVE_DEFINITE);
        is_hermitian = (symmetry == MatrixSymmetry::HERMITIAN);
    }
    
    /**
     * @brief 获取对称性字符串描述
     * @return 对称性描述字符串
     */
    std::string symmetry_string() const {
        switch (symmetry) {
            case MatrixSymmetry::UNSYMMETRIC: return "非对称";
            case MatrixSymmetry::SYMMETRIC: return "对称";
            case MatrixSymmetry::HERMITIAN: return "埃尔米特";
            case MatrixSymmetry::SKEW_SYMMETRIC: return "反对称";
            default: return "未知";
        }
    }
    
    /**
     * @brief 获取正定性字符串描述
     * @return 正定性描述字符串
     */
    std::string definiteness_string() const {
        switch (definiteness) {
            case MatrixDefiniteness::INDEFINITE: return "不定";
            case MatrixDefiniteness::POSITIVE_DEFINITE: return "正定";
            case MatrixDefiniteness::POSITIVE_SEMIDEFINITE: return "半正定";
            case MatrixDefiniteness::NEGATIVE_DEFINITE: return "负定";
            case MatrixDefiniteness::NEGATIVE_SEMIDEFINITE: return "半负定";
            default: return "未知";
        }
    }
    
    /**
     * @brief 获取数据类型字符串描述
     * @return 数据类型描述字符串
     */
    std::string data_type_string() const {
        switch (data_type) {
            case MatrixDataType::REAL: return "实数";
            case MatrixDataType::COMPLEX: return "复数";
            default: return "未知";
        }
    }
    
    /**
     * @brief 获取元素类型字符串描述
     * @return 元素类型描述字符串
     */
    std::string element_type_string() const {
        switch (element_type) {
            case MatrixElementType::SCALAR: return "标量元";
            case MatrixElementType::VECTOR_2D: return "二维矢量元";
            case MatrixElementType::VECTOR_3D: return "三维矢量元";
            default: return "未知";
        }
    }
    
    /**
     * @brief 获取物理场类型字符串描述
     * @return 物理场类型描述字符串
     */
    std::string field_type_string() const {
        switch (field_type) {
            case PhysicalFieldType::ELECTROSTATIC: return "静电场";
            case PhysicalFieldType::MAGNETOSTATIC: return "静磁场";
            case PhysicalFieldType::EDDY_CURRENT: return "涡流场";
            case PhysicalFieldType::WAVE: return "波动场";
            default: return "未知";
        }
    }
    
    /**
     * @brief 获取完整属性描述
     * @return 完整属性描述字符串
     */
    std::string to_string() const {
        return "对称性: " + symmetry_string() + 
               ", 正定性: " + definiteness_string() + 
               ", 数据类型: " + data_type_string() + 
               ", 元素类型: " + element_type_string() + 
               ", 物理场: " + field_type_string() + 
               ", 奇异: " + (is_singular ? "是" : "否") + 
               ", SPD: " + (is_spd ? "是" : "否") + 
               ", 埃尔米特: " + (is_hermitian ? "是" : "否");
    }
    
    /**
     * @brief 检查是否适合使用CG求解器
     * @return 是否适合CG
     */
    bool suitable_for_cg() const {
        return is_spd && data_type == MatrixDataType::REAL;
    }
    
    /**
     * @brief 检查是否适合使用块预处理
     * @return 是否适合块预处理
     */
    bool suitable_for_block_preconditioner() const {
        return element_type != MatrixElementType::SCALAR;
    }
    
    /**
     * @brief 检查是否适合使用ILU预处理
     * @return 是否适合ILU预处理
     */
    bool suitable_for_ilu() const {
        return !is_singular && definiteness != MatrixDefiniteness::INDEFINITE;
    }
    
    /**
     * @brief 检查是否适合使用Jacobi预处理
     * @return 是否适合Jacobi预处理
     */
    bool suitable_for_jacobi() const {
        return definiteness != MatrixDefiniteness::INDEFINITE;
    }
    
    /**
     * @brief 创建静电场景属性
     * @return 静电场景矩阵属性
     */
    static MatrixAttribute create_electrostatic() {
        return MatrixAttribute(
            MatrixSymmetry::SYMMETRIC,
            MatrixDefiniteness::POSITIVE_DEFINITE,
            MatrixDataType::REAL,
            MatrixElementType::SCALAR,
            PhysicalFieldType::ELECTROSTATIC
        );
    }
    
    /**
     * @brief 创建静磁场景属性
     * @return 静磁场景矩阵属性
     */
    static MatrixAttribute create_magnetostatic() {
        return MatrixAttribute(
            MatrixSymmetry::SYMMETRIC,
            MatrixDefiniteness::POSITIVE_SEMIDEFINITE, // 静磁场通常半正定
            MatrixDataType::REAL,
            MatrixElementType::SCALAR,
            PhysicalFieldType::MAGNETOSTATIC,
            true  // 静磁场通常奇异
        );
    }
    
    /**
     * @brief 创建涡流场场景属性
     * @return 涡流场场景矩阵属性
     */
    static MatrixAttribute create_eddy_current() {
        return MatrixAttribute(
            MatrixSymmetry::HERMITIAN,
            MatrixDefiniteness::POSITIVE_DEFINITE,
            MatrixDataType::COMPLEX,
            MatrixElementType::VECTOR_3D, // 涡流场通常使用矢量元
            PhysicalFieldType::EDDY_CURRENT
        );
    }
    
    /**
     * @brief 创建矢量元属性
     * @param dim 维度（2或3）
     * @return 矢量元矩阵属性
     */
    static MatrixAttribute create_vector_element(int dim) {
        MatrixElementType etype = (dim == 2) ? MatrixElementType::VECTOR_2D : MatrixElementType::VECTOR_3D;
        return MatrixAttribute(
            MatrixSymmetry::SYMMETRIC,
            MatrixDefiniteness::POSITIVE_DEFINITE,
            MatrixDataType::REAL,
            etype,
            PhysicalFieldType::MAGNETOSTATIC
        );
    }
};

} // namespace numeric