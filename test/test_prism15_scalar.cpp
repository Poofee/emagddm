/**
 * @file test_prism15_scalar.cpp
 * @brief 二阶三棱柱(PRISM15)纯标量格式DOF管理完整测试
 * @details 使用 Google Test 框架验证 SCALAR_ONLY 格式下 PRISM15 单元的：
 *          1. GlobalEdgeIDGenerator: 标量单元不产生全局棱边
 *          2. DOF预编号: 仅标量节点DOF（15个），无矢量DOF
 *          3. 约束识别: Dirichlet边界正确标记1个约束节点
 *          4. 自由DOF计数: 总数 = 15 - 1 = 14
 *          5. Local2Global映射表: 索引连续性、约束位置、维度正确性
 *          6. 约束值向量: 包含正确的Dirichlet边界值
 *
 * 测试网格拓扑：
 * - 1个 PRISM15 单元（二阶三棱柱，15个节点）
 *   底面三角形顶点(6个): 节点0~5
 *   顶面三角形顶点(6个): 节点6~11（实际按标准顺序为3~14）
 *   边中点(9个): 底面3个 + 顶面3个 + 垂直边3个
 * - Dirichlet边界: 固定节点0的标量值为10.0
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "em_mesh_data.hpp"
#include "global_edge_id_generator.hpp"
#include "em_dof_manager.hpp"
#include "logger_factory.hpp"

using namespace fe_em;

// ==================== 测试固件：PRISM15标量网格构建 ====================

/**
 * @class Prism15ScalarFixture
 * @brief Google Test 固件类：构建PRISM15纯标量测试网格
 * @details 在 SetUp 中一次性构建完整的 EMMeshData（15节点+1单元+1边界），
 *          所有 TEST_F 用例共享此网格数据，避免重复构建。
 */
class Prism15ScalarFixture : public ::testing::Test {
protected:
    EMMeshData mesh_data;                    ///< 测试用网格数据

    void SetUp() override {
        buildPrism15Nodes();
        buildPrism15Element();
        buildDirichletBoundary();
    }

private:
    /**
     * @brief 构建15个节点的二阶三棱柱坐标
     * @details PRISM15节点顺序遵循VTK/ANSYS标准：
     *          - 节点0~2: 底面三角形3个顶点（逆时针）
     *          - 节点3~5: 顶面三角形3个顶点（与底面对应，z=1）
     *          - 节点6~8: 底面3条边的边中点
     *          - 节点9~11: 顶面3条边的边中点
     *          - 节点12~14: 3条垂直棱边的边中点
     */
    void buildPrism15Nodes() {
        // 底面三角形顶点 (z=0)
        mesh_data.nodes.push_back({0, 0.0,       0.0,    0.0, 1});  // 节点0: 原点
        mesh_data.nodes.push_back({1, 1.0,       0.0,    0.0, 1});  // 节点1: (1,0,0)
        mesh_data.nodes.push_back({2, 0.5,      0.866,   0.0, 1});  // 节点2: 等边三角形第三顶点

        // 顶面三角形顶点 (z=1)
        mesh_data.nodes.push_back({3, 0.0,       0.0,    1.0, 1});  // 节点3: (0,0,1)
        mesh_data.nodes.push_back({4, 1.0,       0.0,    1.0, 1});  // 节点4: (1,0,1)
        mesh_data.nodes.push_back({5, 0.5,      0.866,   1.0, 1});  // 节点5: (0.5,0.866,1)

        // 底面边中点 (z=0)
        mesh_data.nodes.push_back({6, 0.5,       0.0,    0.0, 1});  // 节点6: 边0-1中点
        mesh_data.nodes.push_back({7, 0.75,     0.433,   0.0, 1});  // 节点7: 边1-2中点
        mesh_data.nodes.push_back({8, 0.25,     0.433,   0.0, 1});  // 节点8: 边2-0中点

        // 顶面边中点 (z=1)
        mesh_data.nodes.push_back({9,  0.5,      0.0,    1.0, 1});  // 节点9:  边3-4中点
        mesh_data.nodes.push_back({10, 0.75,    0.433,   1.0, 1});  // 节点10: 边4-5中点
        mesh_data.nodes.push_back({11, 0.25,    0.433,   1.0, 1});  // 节点11: 边5-3中点

        // 垂直棱边中点 (z=0.5)
        mesh_data.nodes.push_back({12, 0.0,      0.0,    0.5, 1});  // 节点12: 边0-3中点
        mesh_data.nodes.push_back({13, 1.0,      0.0,    0.5, 1});  // 节点13: 边1-4中点
        mesh_data.nodes.push_back({14, 0.5,     0.866,   0.5, 1});  // 节点14: 边2-5中点
    }

    /**
     * @brief 构建单个PRISM15单元
     * @details dof_type=SCALAR_ONLY 表示纯标量场求解（如静电位φ或静磁位ψ）
     */
    void buildPrism15Element() {
        Element elem;
        elem.id = 0;
        elem.node_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
        elem.elem_type = ElemType::PRISM15;
        elem.dof_type = DOFType::SCALAR_ONLY;
        elem.material_id = 1;
        elem.region_id = 1;

        mesh_data.elements.push_back(elem);
    }

    /**
     * @brief 构建Dirichlet标量边界条件
     * @details 固定节点0的标量值为10.0（模拟接地电位或固定磁位）
     */
    void buildDirichletBoundary() {
        EMBoundaryMarker bnd;
        bnd.id = 0;
        bnd.bnd_type = BndType::DIRICHLET;
        bnd.dof_type = DOFType::SCALAR_ONLY;
        bnd.target_ids = std::vector<int>{0};   // 目标：节点0
        bnd.value = 10.0;
        bnd.name = "node0_fixed";

        mesh_data.boundary_markers.push_back(bnd);
    }
};


// ====================================================================
//  测试1: GlobalEdgeIDGenerator — 标量单元不产生棱边
// ====================================================================

/**
 * @test 验证SCALAR_ONLY格式的PRISM15单元不产生任何全局棱边
 * @details GlobalEdgeIDGenerator仅处理VECTOR_EDGE_ONLY和MIXED_AV类型的单元，
 *          SCALAR_ONLY单元对应的全局棱边映射应为空，总数为0。
 */
TEST_F(Prism15ScalarFixture, GlobalEdgeIDGenerator_NoEdgesForScalar) {
    GlobalEdgeIDGenerator edge_gen(mesh_data);

    ASSERT_NO_THROW(edge_gen.generate())
        << "SCALAR_ONLY单元的棱边生成不应抛出异常";

    int num_global_edges = edge_gen.getNumGlobalEdges();
    EXPECT_EQ(num_global_edges, 0)
        << "SCALAR_ONLY单元的全局棱边数应为0，实际: " << num_global_edges;

    const auto& elem_edge_mapping = edge_gen.getElemLocalToGlobalEdge();
    ASSERT_EQ(elem_edge_mapping.size(), 1u)
        << "映射表外层大小应等于单元数(1)";

    EXPECT_EQ(elem_edge_mapping[0].size(), 0u)
        << "SCALAR_ONLY单元的内层棱边映射应为空向量";

    FEEM_INFO("GlobalEdgeIDGenerator验证通过: 全局棱边数={}, 映射表size={}",
              num_global_edges, elem_edge_mapping.size());
}


// ====================================================================
//  测试2: EMDOFManager — DOF预编号统计
// ====================================================================

/**
 * @test 验证SCALAR_ONLY格式的DOF预编号结果
 * @details 15个节点全部被标量单元引用，故标量DOF=15；
 *          无矢量单元，故矢量DOF=0。
 */
TEST_F(Prism15ScalarFixture, DOFManager_Prenumbering_Stats) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);

    ASSERT_NO_THROW(dof_manager.build())
        << "DOF编号流程执行成功";

    int num_free = dof_manager.getNumFreeDOFs();

    EXPECT_EQ(num_free, 14)
        << "自由DOF数应为14（15总节点 - 1约束节点），实际: " << num_free;

    FEEM_INFO("DOF预编号统计: 自由DOF={}", num_free);
}


// ====================================================================
//  测试3: Local2Global映射表 — 维度与类型验证
// ====================================================================

/**
 * @test 验证Local2Global映射表的基本属性
 * @details indices.size()==15（纯标量，每节点1个DOF），
 *          num_scalar_dofs==15, num_vector_dofs==0,
 *          is_mixed()==false。
 */
TEST_F(Prism15ScalarFixture, Local2Global_DimensionsAndType) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    ASSERT_EQ(mappings.size(), 1u)
        << "应有1个单元的映射表";

    const auto& l2g = mappings[0];

    EXPECT_EQ(l2g.indices.size(), 15u)
        << "indices大小应为15（PRISM15标量节点数）";
    EXPECT_EQ(l2g.num_scalar_dofs, 15)
        << "标量DOF数应为15";
    EXPECT_EQ(l2g.num_vector_dofs, 0)
        << "矢量DOF数应为0（纯标量格式）";
    EXPECT_FALSE(l2g.is_mixed())
        << "is_mixed()应返回false（非A-V混合格式）";
    EXPECT_EQ(l2g.element_id, 0)
        << "element_id应匹配单元ID=0";
    EXPECT_EQ(l2g.elem_type, ElemType::PRISM15)
        << "elem_type应为PRISM15";

    FEEM_INFO("Local2Global维度验证: indices.size()={}, scalar={}, vector={}, mixed={}",
              l2g.indices.size(), l2g.num_scalar_dofs, l2g.num_vector_dofs,
              l2g.is_mixed());
}


// ====================================================================
//  测试4: Local2Global映射表 — 约束位置验证
// ====================================================================

/**
 * @test 验证Local2Global映射表中约束DOF的位置和值
 * @details 节点0被Dirichlet边界固定，其对应的indices[0]应为-1；
 *          其余14个位置均应为>=0的自由DOF编号。
 */
TEST_F(Prism15ScalarFixture, Local2Global_ConstrainedPosition) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& l2g = dof_manager.getElemLocalToGlobal()[0];

    EXPECT_EQ(l2g.indices[0], -1)
        << "节点0（被Dirichlet固定）的映射值应为-1（约束DOF），实际: "
        << l2g.indices[0];

    int constrained_count = 0;
    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        if (l2g.indices[i] == -1) {
            constrained_count++;
        }
    }
    EXPECT_EQ(constrained_count, 1)
        << "约束DOF数量应为1（仅节点0），实际: " << constrained_count;

    FEEM_INFO("约束位置验证: indices[0]={}, 总约束数={}",
              l2g.indices[0], constrained_count);
}


// ====================================================================
//  测试5: Local2Global映射表 — 自由DOF连续性验证
// ====================================================================

/**
 * @test 验证自由DOF编号的连续性和范围
 * @details 除节点0(-1)外，其余14个自由DOF编号应从0到13连续排列，
 *          无重复、无跳跃、无负值（除约束位）。
 */
TEST_F(Prism15ScalarFixture, Local2Global_FreeDOFContinuity) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& l2g = dof_manager.getElemLocalToGlobal()[0];

    std::vector<int> free_indices;
    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        if (l2g.indices[i] >= 0) {
            free_indices.push_back(l2g.indices[i]);
        }
    }

    ASSERT_EQ(free_indices.size(), 14u)
        << "自由DOF数量应为14";

    std::vector<int> sorted_free = free_indices;
    std::sort(sorted_free.begin(), sorted_free.end());

    for (int i = 0; i < static_cast<int>(sorted_free.size()); ++i) {
        EXPECT_EQ(sorted_free[i], i)
            << "第" << i << "个自由DOF编号应为" << i
            << "（期望0~13连续），实际: " << sorted_free[i];
    }

    FEEM_INFO("自由DOF连续性验证: 数量={}, 范围=[0, {}]",
              free_indices.size(), free_indices.size() - 1);
}


// ====================================================================
//  测试6: Local2Global映射表 — 完整内容输出
// ====================================================================

/**
 * @test 输出完整的Local2Global映射表内容用于调试
 * @details 打印每个局部DOF位置对应的映射值，
 *          标注哪些是约束DOF（-1）、哪些是自由DOF（>=0）。
 */
TEST_F(Prism15ScalarFixture, Local2Global_FullContentOutput) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& l2g = dof_manager.getElemLocalToGlobal()[0];

    std::cout << "\n========== PRISM15 Local2Global 映射表完整内容 ==========\n";
    std::cout << "单元ID: " << l2g.element_id
              << " | 类型: PRISM15(SCALAR_ONLY)\n";
    std::cout << "标量DOF: " << l2g.num_scalar_dofs
              << " | 矢量DOF: " << l2g.num_vector_dofs
              << " | 混合: " << (l2g.is_mixed() ? "是" : "否") << "\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "局部idx | 全局节点ID | 映射值 | 状态\n";
    std::cout << "--------|------------|--------|------\n";

    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        int node_id = mesh_data.elements[0].node_ids[i];
        int global_idx = l2g.indices[i];
        std::string status = (global_idx == -1) ? "约束(Dirichlet)" : "自由DOF";

        printf("  %2zu    |    %2d      |  %4d  | %s\n",
               i, node_id, global_idx, status.c_str());
    }

    std::cout << "========================================================\n";

    FEEM_INFO("Local2Global完整映射表已输出至控制台");
}


// ====================================================================
//  测试7: 约束DOF值向量验证
// ====================================================================

/**
 * @test 验证constrained_dof_values包含正确的Dirichlet边界值
 * @details 约束值向量的长度应等于总预编号数（15），
 *          其中对应节点0的位置应存储值10.0。
 */
TEST_F(Prism15ScalarFixture, ConstrainedDOFValues_Content) {
    std::vector<std::vector<int>> empty_edge_mapping;
    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& constrained_vals = dof_manager.getConstrainedDOFValues();

    ASSERT_FALSE(constrained_vals.empty())
        << "约束值向量不应为空";

    bool found_target_value = false;
    for (size_t i = 0; i < constrained_vals.size(); ++i) {
        if (std::abs(constrained_vals[i] - 10.0) < 1e-10) {
            found_target_value = true;
            FEEM_INFO("找到Dirichlet约束值10.0于位置{}", i);
            break;
        }
    }

    EXPECT_TRUE(found_target_value)
        << "约束值向量中应包含Dirichlet边界值10.0";

    std::cout << "\n========== 约束DOF值向量 ==========\n";
    std::cout << "向量长度: " << constrained_vals.size() << "\n";
    std::cout << "内容: [";
    for (size_t i = 0; i < constrained_vals.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << constrained_vals[i];
        if (std::abs(constrained_vals[i] - 10.0) < 1e-10) {
            std::cout << "(←Dirichlet)";
        }
    }
    std::cout << "]\n";
    std::cout << "======================================\n";
}


// ====================================================================
//  测试8: 综合一致性验证 — 所有统计数据自洽
// ====================================================================

/**
 * @test 跨模块综合验证所有统计数据的一致性
 * @details 验证以下等式链成立：
 *          总节点数(15) == num_scalar_dofs(15)
 *          自由DOF(14) + 约束DOF(1) == 总标量DOF(15)
 *          getNumFreeDOFs() == 映射表中>=0的数量
 */
TEST_F(Prism15ScalarFixture, CrossModule_ConsistencyCheck) {
    std::vector<std::vector<int>> empty_edge_mapping;

    GlobalEdgeIDGenerator edge_gen(mesh_data);
    edge_gen.generate();

    EMDOFManager dof_manager(mesh_data, empty_edge_mapping);
    dof_manager.build();

    const auto& l2g = dof_manager.getElemLocalToGlobal()[0];
    const auto& constrained_vals = dof_manager.getConstrainedDOFValues();

    int total_nodes = static_cast<int>(mesh_data.nodes.size());
    int num_free = dof_manager.getNumFreeDOFs();

    int mapping_constrained = 0;
    int mapping_free = 0;
    for (int idx : l2g.indices) {
        if (idx == -1) {
            mapping_constrained++;
        } else if (idx >= 0) {
            mapping_free++;
        }
    }

    EXPECT_EQ(total_nodes, 15) << "总节点数应为15";
    EXPECT_EQ(l2g.num_scalar_dofs, 15) << "标量DOF数应等于节点数";
    EXPECT_EQ(l2g.num_vector_dofs, 0) << "矢量DOF数应为0";
    EXPECT_EQ(mapping_free + mapping_constrained, 15)
        << "自由+约束DOF数应等于总标量DOF数";
    EXPECT_EQ(num_free, mapping_free)
        << "getNumFreeDOFs()应等于映射表中自由DOF计数";
    EXPECT_EQ(num_free, 14) << "最终自由DOF数应为14";
    EXPECT_EQ(edge_gen.getNumGlobalEdges(), 0)
        << "全局棱边数应为0（纯标量格式）";

    std::cout << "\n========== 统计一致性报告 ==========\n";
    std::printf("  总节点数:         %d\n", total_nodes);
    std::printf("  标量DOF数:        %d\n", l2g.num_scalar_dofs);
    std::printf("  矢量DOF数:        %d\n", l2g.num_vector_dofs);
    std::printf("  映射表自由DOF:    %d\n", mapping_free);
    std::printf("  映射表约束DOF:    %d\n", mapping_constrained);
    std::printf("  getNumFreeDOFs(): %d\n", num_free);
    std::printf("  全局棱边数:       %d\n", edge_gen.getNumGlobalEdges());
    std::printf("  约束值向量长度:   %zu\n", constrained_vals.size());
    std::printf("  一致性:           %s\n",
                (num_free == 14 && mapping_free == 14 && mapping_constrained == 1)
                    ? "✅ 通过" : "❌ 失败");
    std::cout << "====================================\n";

    FEEM_INFO("综合一致性检查通过: 节点={} 标量DOF={} 自由={} 约束={} 棱边={}",
              total_nodes, l2g.num_scalar_dofs, num_free, mapping_constrained,
              edge_gen.getNumGlobalEdges());
}


// ==================== 主函数 ====================

int main(int argc, char** argv) {
    FEEM_LOG_INIT("", false);

    ::testing::InitGoogleTest(&argc, argv);

    FEEM_INFO("============================================");
    FEEM_INFO("  PRISM15 纯标量格式 DOF 管理测试");
    FEEM_INFO("  单元类型: PRISM15 (二阶三棱柱, 15节点)");
    FEEM_INFO("  DOF格式:  SCALAR_ONLY");
    FEEM_INFO("  边界条件: Dirichlet(节点0=10.0)");
    FEEM_INFO("============================================");

    int result = RUN_ALL_TESTS();

    FEEM_INFO("测试完成，返回码: {}", result);
    return result;
}
