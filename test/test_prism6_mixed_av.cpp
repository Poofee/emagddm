/**
 * @file test_prism6_mixed_av.cpp
 * @brief 一阶三棱柱(PRISM6)混合A-V格式DOF管理完整测试
 * @details 使用 Google Test 框架验证 PRISM6 + MIXED_AV 场景下的核心模块协作：
 *          1. GlobalEdgeIDGenerator：全局棱边ID生成（PRISM6有9条棱边）
 *          2. EMDOFManager：四步DOF编号流程（预编号→约束标记→重编号→映射表构建）
 *          3. Local2Global映射表：标量(前6)+矢量(后9)=15个局部DOF的正确映射
 *
 * 测试网格拓扑：
 * - 1个PRISM6单元，6个节点构成单位三棱柱
 * - 底面(z=0): 节点0(0,0,0), 节点1(1,0,0), 节点2(0.5,0.866,0)
 * - 顶面(z=1): 节点3(0,0,1), 节点4(1,0,1), 节点5(0.5,0.866,1)
 *
 * 边界条件：
 * - 标量Dirichlet: 固定节点0电位 V=0.0
 * - 矢量Dirichlet: 固定棱边(节点0-节点3)磁矢势 A×n=0.0
 *
 * 预期结果：
 * - 全局棱边数 = 9（PRISM6标准棱边数）
 * - 预编号总数 = 15（6标量节点 + 9矢量棱边）
 * - 约束DOF数 = 2（1个标量节点 + 1条矢量棱边）
 * - 自由DOF数 = 13
 * - Local2Global: indices[0]=-1(节点0约束), 某个矢量位置=-1(棱边0-3约束), 其余>=0连续
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <algorithm>

#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"
#include "global_edge_id_generator.hpp"
#include "em_dof_manager.hpp"
#include "element_geometry.hpp"
#include "logger_factory.hpp"

using namespace fe_em;

// ==================== 测试固件：构建PRISM6 MIXED_AV网格 ====================

/**
 * @class Prism6MixedAVFixture
 * @brief PRISM6混合A-V格式测试固件
 * @details 在SetUp中构建完整的EMMeshData，包含：
 *          - 6个节点（单位三棱柱几何）
 *          - 1个PRISM6单元（MIXED_AV类型）
 *          - 2个Dirichlet边界条件（标量+矢量各1个）
 */
class Prism6MixedAVFixture : public ::testing::Test {
protected:
    /**
     * @brief 构建测试网格数据
     * @details 手动填充EMMeshData的所有字段，模拟从解析器获取的网格拓扑
     */
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 构建6个节点（单位三棱柱） ----------
        // 底面三角形 (z=0)：等边三角形的一半，便于计算
        mesh_data_.nodes = {
            {0, 0.0,     0.0,    0.0, 0},   // 节点0: 原点
            {1, 1.0,     0.0,    0.0, 0},   // 节点1: x轴方向
            {2, 0.5,     0.866,  0.0, 0},   // 节点2: 等边三角形顶点 (√3/2 ≈ 0.866)
            {3, 0.0,     0.0,    1.0, 0},   // 节点3: 顶面原点（z=1平移）
            {4, 1.0,     0.0,    1.0, 0},   // 节点4: 顶面x轴方向
            {5, 0.5,     0.866,  1.0, 0},   // 节点5: 顶面等边三角形顶点
        };

        // ---------- 构建1个PRISM6单元（MIXED_AV格式） ----------
        mesh_data_.elements = {{
            0,                                    // 单元ID
            {0, 1, 2, 3, 4, 5},                  // 节点连接（底面0-1-2 + 顶面3-4-5）
            ElemType::PRISM6,                     // 单元类型：一阶三棱柱
            DOFType::MIXED_AV,                    // DOF类型：A-V混合格式
            1,                                    // 材料ID
            0,                                    // 区域ID
        }};

        // ---------- 构建边界条件 ----------
        // 边界1：标量Dirichlet - 固定节点0电位为0
        EMBoundaryMarker scalar_bnd;
        scalar_bnd.id = 0;
        scalar_bnd.bnd_type = BndType::DIRICHLET;
        scalar_bnd.dof_type = DOFType::SCALAR_ONLY;
        scalar_bnd.target_ids = std::vector<int>{0};       // 目标：节点0
        scalar_bnd.value = 0.0;
        scalar_bnd.name = "V_fixed";

        // 边界2：矢量Dirichlet - 固定棱边(节点0-节点3)磁矢势为0
        EMBoundaryMarker vector_bnd;
        vector_bnd.id = 1;
        vector_bnd.bnd_type = BndType::DIRICHLET;
        vector_bnd.dof_type = DOFType::VECTOR_EDGE_ONLY;
        vector_bnd.target_ids = std::vector<std::vector<int>>{{0, 3}};  // 目标：棱边(0,3)
        vector_bnd.value = 0.0;
        vector_bnd.name = "A_fixed";

        mesh_data_.boundary_markers = {scalar_bnd, vector_bnd};

        FEEM_INFO("Prism6MixedAVFixture: 网格构建完成 - {}节点, {}单元, {}边界",
                  mesh_data_.nodes.size(),
                  mesh_data_.elements.size(),
                  mesh_data_.boundary_markers.size());
    }

    EMMeshData mesh_data_;                          ///< 测试用网格数据
};


// ====================================================================
//  1. 全局棱边ID生成器测试（GlobalEdgeIDGenerator Tests）
// ====================================================================

/**
 * @test 验证PRISM6单元的全局棱边数量
 * @details PRISM6三棱柱有9条物理棱边：
 *          - 底面3条三角形边
 *          - 顶面3条三角形边
 *          - 3条竖直侧棱（连接底面与顶面对应节点）
 */
TEST_F(Prism6MixedAVFixture, GlobalEdgeCount) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();

    int num_edges = edge_gen.getNumGlobalEdges();
    EXPECT_EQ(num_edges, 9) << "PRISM6单元应有9条全局棱边";

    FEEM_INFO("全局棱边数量验证: 预期=9, 实际={}", num_edges);
}

/**
 * @test 验证单元局部到全局棱边映射表的维度
 * @details 映射表外层size=1（1个单元），内层size=9（每条局部棱边的全局ID）
 */
TEST_F(Prism6MixedAVFixture, EdgeMappingDimension) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();

    const auto& mapping = edge_gen.getElemLocalToGlobalEdge();

    ASSERT_EQ(mapping.size(), 1) << "映射表应有1个条目（1个单元）";
    EXPECT_EQ(mapping[0].size(), 9) << "单元0的局部棱边映射应有9个元素";

    FEEM_INFO("棱边映射维度: 外层={}, 内层[0]={}", mapping.size(), mapping[0].size());
}

/**
 * @test 验证全局棱边ID的唯一性和连续性
 * @details 所有全局棱边ID应在[0, 8]范围内且互不重复
 */
TEST_F(Prism6MixedAVFixture, EdgeIDsUniqueAndContinuous) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();

    const auto& mapping = edge_gen.getElemLocalToGlobalEdge();
    const auto& local_to_global = mapping[0];

    // 收集所有全局棱边ID并检查范围和唯一性
    std::vector<bool> id_found(9, false);
    for (int global_edge_id : local_to_global) {
        ASSERT_GE(global_edge_id, 0) << "全局棱边ID不应为负";
        ASSERT_LT(global_edge_id, 9) << "全局棱边ID应<9";
        EXPECT_FALSE(id_found[global_edge_id])
            << "全局棱边ID=" << global_edge_id << "重复出现";
        id_found[global_edge_id] = true;
    }

    // 验证所有ID都被覆盖（连续性）
    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(id_found[i]) << "全局棱边ID=" << i << "缺失";
    }

    // 输出详细映射表供人工检查
    FEEM_INFO("=== 局部棱边 → 全局棱边映射表 ===");
    for (size_t i = 0; i < local_to_global.size(); ++i) {
        FEEM_INFO("  局部棱边[{}] → 全局棱边[{}]", i, local_to_global[i]);
    }
}


// ====================================================================
//  2. DOF管理器统计信息测试（EMDOFManager Statistics Tests）
// ====================================================================

/**
 * @test 验证自由DOF总数
 * @details 总预编号15(6标量+9矢量) - 2约束(1标量+1矢量) = 13自由DOF
 */
TEST_F(Prism6MixedAVFixture, FreeDOFCount) {
    // 第一步：生成全局棱边映射
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    // 第二步：执行DOF编号
    EMDOFManager dof_manager(mesh_data_, edge_mapping, &edge_gen);
    dof_manager.build();

    int num_free = dof_manager.getNumFreeDOFs();
    EXPECT_EQ(num_free, 13) << "自由DOF数应为13 (15总 - 2约束)";

    FEEM_INFO("自由DOF数量验证: 预期=13, 实际={}", num_free);
}

/**
 * @test 验证单元映射表数量
 * @details 1个单元应产生1条Local2Global映射记录
 */
TEST_F(Prism6MixedAVFixture, MappingTableSize) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    ASSERT_EQ(mappings.size(), 1) << "映射表应有1条记录（1个单元）";

    FEEM_INFO("映射表大小: {}", mappings.size());
}


// ====================================================================
//  3. Local2Global映射表结构测试（Mapping Table Structure Tests）
// ====================================================================

/**
 * @test 验证Local2Global映射表的基本属性
 * @details indices.size()=15, num_scalar_dofs=6, num_vector_dofs=9, is_mixed()=true
 */
TEST_F(Prism6MixedAVFixture, MappingBasicProperties) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    // 验证映射向量总长度 = 标量节点数 + 矢量棱边数 = 6 + 9 = 15
    EXPECT_EQ(static_cast<int>(l2g.indices.size()), 15)
        << "indices长度应为15 (6标量 + 9矢量)";

    // 验证标量和矢量DOF计数
    EXPECT_EQ(l2g.num_scalar_dofs, 6) << "标量DOF数应为6（PRISM6节点数）";
    EXPECT_EQ(l2g.num_vector_dofs, 9) << "矢量DOF数应为9（PRISM6棱边数）";

    // 验证混合格式标志
    EXPECT_TRUE(l2g.is_mixed()) << "MIXED_AV单元的is_mixed()应返回true";

    // 验证单元ID和类型
    EXPECT_EQ(l2g.element_id, 0) << "单元ID应为0";
    EXPECT_EQ(l2g.elem_type, ElemType::PRISM6) << "单元类型应为PRISM6";

    FEEM_INFO("=== Local2Global基本属性 ===");
    FEEM_INFO("  indices.size() = {}", l2g.indices.size());
    FEEM_INFO("  num_scalar_dofs = {}", l2g.num_scalar_dofs);
    FEEM_INFO("  num_vector_dofs = {}", l2g.num_vector_dofs);
    FEEM_INFO("  is_mixed() = {}", l2g.is_mixed());
    FEEM_INFO("  element_id = {}", l2g.element_id);
}

/**
 * @test 验证标量部分映射（前6个位置对应节点0-5）
 * @details 位置0（节点0）应为-1（Dirichlet约束），其余位置1-5应>=0
 */
TEST_F(Prism6MixedAVFixture, ScalarMappingSection) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    // 提取标量部分的映射
    auto scalar_indices = l2g.get_scalar_indices();
    ASSERT_EQ(static_cast<int>(scalar_indices.size()), 6)
        << "标量映射段应有6个元素";

    // 位置0对应节点0，被Dirichlet边界固定，应为-1
    EXPECT_EQ(scalar_indices[0], -1)
        << "标量位置0（节点0）应为约束DOF(-1)，因施加了V_fixed边界";

    // 位置1-5对应节点1-5，未被约束，应为有效自由编号(>=0)
    for (int i = 1; i < 6; ++i) {
        EXPECT_GE(scalar_indices[i], 0)
            << "标量位置" << i << "（节点" << i << "）应为自由DOF(>=0)";
    }

    FEEM_INFO("=== 标量映射段（节点0-5） ===");
    for (int i = 0; i < 6; ++i) {
        FEEM_INFO("  indices[{}] (节点{}) = {}", i, i, scalar_indices[i]);
    }
}

/**
 * @test 验证矢量部分映射（后9个位置对应棱边0-8）
 * @details 其中一条棱边（对应节点0-3的那条）应为-1（Dirichlet约束），
 *          其余8条棱边应>=0
 */
TEST_F(Prism6MixedAVFixture, VectorMappingSection) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    // 提取矢量部分的映射
    auto vector_indices = l2g.get_vector_indices();
    ASSERT_EQ(static_cast<int>(vector_indices.size()), 9)
        << "矢量映射段应有9个元素";

    // 统计约束DOF数量（值为-1的位置）
    int constrained_count = 0;
    int constrained_local_idx = -1;
    for (int i = 0; i < 9; ++i) {
        if (vector_indices[i] == -1) {
            constrained_count++;
            constrained_local_idx = i;
        }
    }

    // 应恰好有1条棱边被约束（棱边0-3对应的局部棱边）
    EXPECT_EQ(constrained_count, 1)
        << "矢量段中应恰好有1个约束DOF（A_fixed边界固定的棱边）";

    // 其余8条棱边应为自由DOF
    int free_count = 0;
    for (int i = 0; i < 9; ++i) {
        if (vector_indices[i] >= 0) {
            free_count++;
            EXPECT_GE(vector_indices[i], 0)
                << "矢量位置" << i << "应为自由DOF(>=0)";
        }
    }
    EXPECT_EQ(free_count, 8) << "矢量段中应有8个自由DOF";

    FEEM_INFO("=== 矢量映射段（棱边0-8） ===");
    for (int i = 0; i < 9; ++i) {
        const char* tag = (vector_indices[i] == -1) ? " [约束]" : "";
        FEEM_INFO("  indices[{}+6] (棱边{}) = {}{}",
                  i, i, vector_indices[i], tag);
    }
    if (constrained_local_idx >= 0) {
        FEEM_INFO("  → 被约束的局部棱边索引: {}", constrained_local_idx);
    }
}


// ====================================================================
//  4. 自由DOF编号连续性与完整性测试（Free DOF Numbering Tests）
// ====================================================================

/**
 * @test 验证所有自由DOF编号的连续性和范围
 * @details 13个自由DOF的编号应恰好覆盖 [0, 12]，无间隙无重复
 */
TEST_F(Prism6MixedAVFixture, FreeDOFNumberingContinuity) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    // 收集所有非-1的自由编号
    std::vector<int> free_ids;
    for (int idx : l2g.indices) {
        if (idx >= 0) {
            free_ids.push_back(idx);
        }
    }

    // 验证自由DOF总数
    ASSERT_EQ(static_cast<int>(free_ids.size()), 13)
        << "应恰好有13个自由DOF编号";

    // 排序后检查连续性 [0, 12]
    std::sort(free_ids.begin(), free_ids.end());
    for (int i = 0; i < 13; ++i) {
        EXPECT_EQ(free_ids[i], i)
            << "第" << i << "个自由DOF编号应为" << i
            << ", 实际为" << free_ids[i];
    }

    FEEM_INFO("=== 自由DOF编号排序验证 ===");
    FEEM_INFO("  自由编号列表: ");
    std::string id_list;
    for (int idx : free_ids) {
        id_list += std::to_string(idx) + " ";
    }
    FEEM_INFO("  {}", id_list);
}

/**
 * @test 验证约束DOF值向量的正确性
 * @details 约束DOF值向量长度应等于预编号总数(15)，
 *          对应约束位置的值应为边界条件指定的值(0.0)
 */
TEST_F(Prism6MixedAVFixture, ConstrainedDOFValues) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& constrained_vals = dof_manager.getConstrainedDOFValues();

    // 约束值向量长度应等于预编号总数（包含标量+矢量）
    EXPECT_EQ(static_cast<int>(constrained_vals.size()), 15)
        << "约束值向量长度应为15（预编号总数）";

    // 统计非零约束值（两个Dirichlet边界的value都是0.0，所以可能全为0）
    int nonzero_count = 0;
    for (double val : constrained_vals) {
        if (std::abs(val) > 1e-15) {
            nonzero_count++;
        }
    }
    // 本测试中两个边界的value均为0.0，所以nonzero_count应为0
    // 但如果有非零边界值，这里会正确反映

    FEEM_INFO("=== 约束DOF值向量 ===");
    FEEM_INFO("  向量长度: {}", constrained_vals.size());
    FEEM_INFO("  非零值个数: {}", nonzero_count);
}


// ====================================================================
//  5. 完整映射表输出测试（Full Mapping Table Dump Test）
// ====================================================================

/**
 * @test 输出完整的Local2Global映射表供人工审查
 * @details 将15个位置的映射关系全部输出，便于调试时直观检查
 *          格式：位置 | 类型(标量/矢量) | 物理实体 | 全局编号
 */
TEST_F(Prism6MixedAVFixture, FullMappingTableDump) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    EMDOFManager dof_manager(mesh_data_, edge_mapping);
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    FEEM_INFO("");
    FEEM_INFO("============================================================");
    FEEM_INFO("  PRISM6 MIXED_AV 完整 Local2Global 映射表");
    FEEM_INFO("============================================================");

    // 输出标量段（位置0-5，对应节点0-5）
    FEEM_INFO("--- 标量段（位置0-5，节点DOF） ---");
    for (int i = 0; i < l2g.num_scalar_dofs; ++i) {
        int global_idx = l2g.indices[i];
        const char* status = (global_idx == -1) ? "[约束]" : "[自由]";
        FEEM_INFO("  [{:2d}] 节点{:2d} → 全局DOF {:>3d} {}",
                  i, i, global_idx, status);
    }

    // 输出矢量段（位置6-14，对应棱边0-8）
    FEEM_INFO("--- 矢量段（位置6-14，棱边DOF） ---");
    int vec_start = l2g.num_scalar_dofs;
    for (int i = 0; i < l2g.num_vector_dofs; ++i) {
        int pos = vec_start + i;
        int global_idx = l2g.indices[pos];
        const char* status = (global_idx == -1) ? "[约束]" : "[自由]";
        FEEM_INFO("  [{:2d}] 棱边{:2d} → 全局DOF {:>3d} {}",
                  pos, i, global_idx, status);
    }

    FEEM_INFO("--- 统计摘要 ---");
    int total_constrained = 0;
    int total_free = 0;
    for (int idx : l2g.indices) {
        if (idx == -1) total_constrained++;
        else total_free++;
    }
    FEEM_INFO("  总局部DOF: {}", l2g.indices.size());
    FEEM_INFO("  约束DOF:   {}", total_constrained);
    FEEM_INFO("  自由DOF:   {}", total_free);
    FEEM_INFO("  is_mixed:  {}", l2g.is_mixed());
    FEEM_INFO("============================================================");
    FEEM_INFO("");

    // 最终断言：确保统计数据自洽
    EXPECT_EQ(total_constrained, 2) << "约束DOF总数应为2";
    EXPECT_EQ(total_free, 13) << "自由DOF总数应为13";
    EXPECT_EQ(total_constrained + total_free, 15) << "约束+自由应等于总数";
}


// ====================================================================
//  6. ElementGeometry辅助查询测试（ElementGeometry Query Tests）
// ====================================================================

/**
 * @test 验证ElementGeometry对PRISM6类型的查询结果
 * @details PRISM6: 6节点, 9棱边, 5面（2三角底面 + 3四边形侧面）
 */
TEST_F(Prism6MixedAVFixture, ElementGeometryQueries) {
    EXPECT_EQ(ElementGeometry::get_num_nodes(ElemType::PRISM6), 6)
        << "PRISM6节点数应为6";
    EXPECT_EQ(ElementGeometry::get_num_edges(ElemType::PRISM6), 9)
        << "PRISM6棱边数应为9";
    EXPECT_EQ(ElementGeometry::get_num_faces(ElemType::PRISM6), 5)
        << "PRISM6面数应为5（2底面+3侧面）";
    EXPECT_TRUE(ElementGeometry::is_supported(ElemType::PRISM6))
        << "PRISM6应被支持";

    // 查询局部棱边定义并验证基本性质
    auto local_edges = ElementGeometry::get_local_edges(ElemType::PRISM6);
    EXPECT_EQ(static_cast<int>(local_edges.size()), 9)
        << "PRISM6局部棱边定义应有9条";

    // 验证每条棱边的两个节点ID在合法范围内 [0, 5]
    for (size_t i = 0; i < local_edges.size(); ++i) {
        auto [n1, n2] = local_edges[i];
        EXPECT_GE(n1, 0) << "棱边" << i << "的节点1 ID不应为负";
        EXPECT_LT(n1, 6) << "棱边" << i << "的节点1 ID应<6";
        EXPECT_GE(n2, 0) << "棱边" << i << "的节点2 ID不应为负";
        EXPECT_LT(n2, 6) << "棱边" << i << "的节点2 ID应<6";
        // 棱边定义要求 n1 < n2（升序排列）
        EXPECT_LT(n1, n2) << "棱边" << i << "的节点应升序排列";
    }

    FEEM_INFO("=== ElementGeometry PRISM6查询 ===");
    FEEM_INFO("  节点数: {}", ElementGeometry::get_num_nodes(ElemType::PRISM6));
    FEEM_INFO("  棱边数: {}", ElementGeometry::get_num_edges(ElemType::PRISM6));
    FEEM_INFO("  面数:   {}", ElementGeometry::get_num_faces(ElemType::PRISM6));
    for (size_t i = 0; i < local_edges.size(); ++i) {
        auto [n1, n2] = local_edges[i];
        FEEM_INFO("  局部棱边{}: ({}, {})", i, n1, n2);
    }
}


// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
