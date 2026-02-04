/**
 * @file test_sparse_matrix_gtest_stage3_simple.cpp
 * @brief 稀疏矩阵模块阶段3简化测试用例（Google Test版本）
 * @author Poofee
 * @date 2026-02-04
 * @version 1.0
 * 
 * @details 测试阶段3基础功能：矩阵属性标记、简化测试
 */

#include <gtest/gtest.h>
#include "numeric/matrix_attribute.hpp"
#include <complex>

using namespace numeric;

/**
 * @brief 矩阵属性标记测试类
 */
class MatrixAttributeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 设置测试数据
    }
    
    void TearDown() override {
        // 清理测试数据
    }
};

// ============================================================================
// 矩阵属性标记测试
// ============================================================================

TEST_F(MatrixAttributeTest, DefaultConstructor) {
    MatrixAttribute attr;
    
    EXPECT_EQ(attr.symmetry, MatrixSymmetry::UNSYMMETRIC);
    EXPECT_EQ(attr.definiteness, MatrixDefiniteness::INDEFINITE);
    EXPECT_EQ(attr.data_type, MatrixDataType::REAL);
    EXPECT_EQ(attr.element_type, MatrixElementType::SCALAR);
    EXPECT_EQ(attr.field_type, PhysicalFieldType::ELECTROSTATIC);
    EXPECT_FALSE(attr.is_singular);
    EXPECT_FALSE(attr.is_spd);
    EXPECT_FALSE(attr.is_hermitian);
}

TEST_F(MatrixAttributeTest, ParameterizedConstructor) {
    MatrixAttribute attr(
        MatrixSymmetry::SYMMETRIC,
        MatrixDefiniteness::POSITIVE_DEFINITE,
        MatrixDataType::COMPLEX,
        MatrixElementType::VECTOR_3D,
        PhysicalFieldType::EDDY_CURRENT,
        false
    );
    
    EXPECT_EQ(attr.symmetry, MatrixSymmetry::SYMMETRIC);
    EXPECT_EQ(attr.definiteness, MatrixDefiniteness::POSITIVE_DEFINITE);
    EXPECT_EQ(attr.data_type, MatrixDataType::COMPLEX);
    EXPECT_EQ(attr.element_type, MatrixElementType::VECTOR_3D);
    EXPECT_EQ(attr.field_type, PhysicalFieldType::EDDY_CURRENT);
    EXPECT_FALSE(attr.is_singular);
    EXPECT_TRUE(attr.is_spd);
    EXPECT_FALSE(attr.is_hermitian);
}

TEST_F(MatrixAttributeTest, StringConversion) {
    MatrixAttribute attr(
        MatrixSymmetry::HERMITIAN,
        MatrixDefiniteness::POSITIVE_SEMIDEFINITE,
        MatrixDataType::COMPLEX,
        MatrixElementType::VECTOR_2D,
        PhysicalFieldType::MAGNETOSTATIC,
        true
    );
    
    EXPECT_EQ(attr.symmetry_string(), "埃尔米特");
    EXPECT_EQ(attr.definiteness_string(), "半正定");
    EXPECT_EQ(attr.data_type_string(), "复数");
    EXPECT_EQ(attr.element_type_string(), "二维矢量元");
    EXPECT_EQ(attr.field_type_string(), "静磁场");
    
    std::string desc = attr.to_string();
    EXPECT_NE(desc.find("埃尔米特"), std::string::npos);
    EXPECT_NE(desc.find("半正定"), std::string::npos);
    EXPECT_NE(desc.find("复数"), std::string::npos);
}

TEST_F(MatrixAttributeTest, SolverSuitability) {
    // 测试静电场景（适合CG）
    MatrixAttribute electrostatic = MatrixAttribute::create_electrostatic();
    EXPECT_TRUE(electrostatic.suitable_for_cg());
    EXPECT_TRUE(electrostatic.suitable_for_ilu());
    EXPECT_TRUE(electrostatic.suitable_for_jacobi());
    EXPECT_FALSE(electrostatic.suitable_for_block_preconditioner());
    
    // 测试静磁场景（适合MINRES）
    MatrixAttribute magnetostatic = MatrixAttribute::create_magnetostatic();
    EXPECT_FALSE(magnetostatic.suitable_for_cg());
    EXPECT_FALSE(magnetostatic.suitable_for_ilu()); // 奇异矩阵不适合ILU
    EXPECT_TRUE(magnetostatic.suitable_for_jacobi());
    EXPECT_FALSE(magnetostatic.suitable_for_block_preconditioner());
    
    // 测试涡流场场景（适合块预处理）
    MatrixAttribute eddy_current = MatrixAttribute::create_eddy_current();
    EXPECT_FALSE(eddy_current.suitable_for_cg()); // 复数矩阵不适合CG
    EXPECT_TRUE(eddy_current.suitable_for_ilu());
    EXPECT_TRUE(eddy_current.suitable_for_jacobi());
    EXPECT_TRUE(eddy_current.suitable_for_block_preconditioner());
}

TEST_F(MatrixAttributeTest, StaticFactoryMethods) {
    // 测试静电场景创建
    MatrixAttribute electrostatic = MatrixAttribute::create_electrostatic();
    EXPECT_EQ(electrostatic.symmetry, MatrixSymmetry::SYMMETRIC);
    EXPECT_EQ(electrostatic.definiteness, MatrixDefiniteness::POSITIVE_DEFINITE);
    EXPECT_EQ(electrostatic.data_type, MatrixDataType::REAL);
    EXPECT_EQ(electrostatic.element_type, MatrixElementType::SCALAR);
    EXPECT_EQ(electrostatic.field_type, PhysicalFieldType::ELECTROSTATIC);
    EXPECT_FALSE(electrostatic.is_singular);
    
    // 测试静磁场景创建
    MatrixAttribute magnetostatic = MatrixAttribute::create_magnetostatic();
    EXPECT_EQ(magnetostatic.symmetry, MatrixSymmetry::SYMMETRIC);
    EXPECT_EQ(magnetostatic.definiteness, MatrixDefiniteness::POSITIVE_SEMIDEFINITE);
    EXPECT_TRUE(magnetostatic.is_singular);
    
    // 测试涡流场场景创建
    MatrixAttribute eddy_current = MatrixAttribute::create_eddy_current();
    EXPECT_EQ(eddy_current.symmetry, MatrixSymmetry::HERMITIAN);
    EXPECT_EQ(eddy_current.data_type, MatrixDataType::COMPLEX);
    EXPECT_EQ(eddy_current.element_type, MatrixElementType::VECTOR_3D);
    
    // 测试矢量元创建
    MatrixAttribute vector_2d = MatrixAttribute::create_vector_element(2);
    EXPECT_EQ(vector_2d.element_type, MatrixElementType::VECTOR_2D);
    
    MatrixAttribute vector_3d = MatrixAttribute::create_vector_element(3);
    EXPECT_EQ(vector_3d.element_type, MatrixElementType::VECTOR_3D);
}

// ============================================================================
// 主函数
// ============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}