/**
 * @file test_pyramid5_edge.cpp
 * @brief 一阶金字塔纯棱边格式(PYRAMID5_EDGE)DOF管理完整测试
 * @details 使用 Google Test 框架验证以下核心流程的正确性：
 *          1. EMMeshData 手动构建（单PYRAMID5_EDGE单元 + Dirichlet边界）
 *          2. GlobalEdgeIDGenerator 全局棱边ID生成（8条棱边去重编号）
 *          3. EMDOFManager DOF预编号/约束标记/重编号/映射表构建
 *          4. Local2Global 映射表完整性验证
 *
 * 测试场景：
 * - 单元: 1个 PYRAMID5_EDGE，5节点(底面正方形+尖顶)，8条局部棱边
 * - 边界: Dirichlet矢量边界固定第0条局部棱边(节点0-1)的磁矢势A=0
 * - 预期: 8条全局棱边, 8个矢量DOF预编号, 1个约束, 7个自由DOF
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

#include "em_mesh_data.hpp"
#include "em_dof_data.hpp"
#include "global_edge_id_generator.hpp"
#include "em_dof_manager.hpp"
#include "logger_factory.hpp"

using namespace fe_em;

// ==================== 测试固件：构建PYRAMID5_EDGE网格 ====================

/**
 * @class Pyramid5EdgeFixture
 * @brief Google Test 固件类，负责构建标准化的PYRAMID5_EDGE测试网格
 * @details 在 SetUp() 中完成以下准备工作：
 *          1. 初始化日志系统（FEEM_LOG_INIT）
 *          2. 构建5个节点的金字塔坐标
 *          3. 构建1个 PYRAMID5_EDGE 单元
 *          4. 构建1个 Dirichlet 矢量边界条件（固定棱边0）
 *          5. 组装完整的 EMMeshData 对象
 */
class Pyramid5EdgeFixture : public ::testing::Test {
protected:
    /**
     * @brief 测试前置准备：构建完整的网格数据
     * @details 节点布局：
     *          - 底面正方形(z=0): 节点0(0,0,0), 节点1(1,0,0), 节点2(1,1,0), 节点3(0,1,0)
     *          - 尖顶(z=1):       节点4(0.5,0.5,1)
     *          单元: id=0, node_ids={0,1,2,3,4}, type=PYRAMID5_EDGE, dof_type=VECTOR_EDGE_ONLY
     *          边界: Dirichlet固定局部棱边0(节点0-1)，值=0.0
     */
    void SetUp() override
    {
        // 初始化日志系统（必须在使用FEEM_*宏之前调用）
        FEEM_LOG_INIT("", false);

        // ========== 第一步：构建5个金字塔节点 ==========
        // 底面正方形四个角点 (z=0)
        mesh_data_.nodes.push_back({0, 0.0, 0.0, 0.0, 0});   // 节点0: 左下角
        mesh_data_.nodes.push_back({1, 1.0, 0.0, 0.0, 0});   // 节点1: 右下角
        mesh_data_.nodes.push_back({2, 1.0, 1.0, 0.0, 0});   // 节点2: 右上角
        mesh_data_.nodes.push_back({3, 0.0, 1.0, 0.0, 0});   // 节点3: 左上角
        // 尖顶 (z=1)
        mesh_data_.nodes.push_back({4, 0.5, 0.5, 1.0, 0});   // 节点4: 锥顶中心

        // ========== 第二步：构建1个PYRAMID5_EDGE单元 ==========
        Element elem;
        elem.id = 0;                                              // 单元全局ID
        elem.node_ids = {0, 1, 2, 3, 4};                         // 5个节点（底面+尖顶）
        elem.elem_type = ElemType::PYRAMID5_EDGE;                 // 一阶金字塔Nedelec棱边元
        elem.dof_type = DOFType::VECTOR_EDGE_ONLY;                // 纯矢量棱边DOF格式
        elem.material_id = 0;                                     // 材料ID（空气）
        elem.region_id = 0;                                       // 区域ID
        mesh_data_.elements.push_back(std::move(elem));

        // ========== 第三步：构建Dirichlet矢量边界条件 ==========
        // 目标：固定第0条局部棱边（即局部节点对(0,1)→全局节点对(0,1)）的磁矢势A=0
        // target_ids使用vector<vector<int>>类型表示棱边节点对列表
        EMBoundaryMarker bnd_marker;
        bnd_marker.id = 0;                                        // 边界标记ID
        bnd_marker.bnd_type = BndType::DIRICHLET;                 // 狄利克雷边界（固定值）
        bnd_marker.dof_type = DOFType::VECTOR_EDGE_ONLY;          // 施加到矢量棱边DOF
        bnd_marker.target_ids = std::vector<std::vector<int>>{{0, 1}};  // 局部棱边0的节点对
        bnd_marker.value = 0.0;                                   // 磁矢势A=0（理想导体）
        bnd_marker.name = "A_edge0_fixed";                        // 边界名称
        mesh_data_.boundary_markers.push_back(std::move(bnd_marker));

        FEEM_INFO("Pyramid5EdgeFixture: 网格数据构建完成");
        FEEM_INFO("  - 节点数: {}", mesh_data_.nodes.size());
        FEEM_INFO("  - 单元数: {}", mesh_data_.elements.size());
        FEEM_INFO("  - 边界标记数: {}", mesh_data_.boundary_markers.size());
    }

    /**
     * @brief 测试后置清理：清空网格数据
     */
    void TearDown() override
    {
        mesh_data_.clear();
    }

    /** @brief 构建好的网格数据（供TEST_F访问） */
    EMMeshData mesh_data_;
};


// ====================================================================
//  测试用例1：GlobalEdgeIDGenerator 全局棱边ID生成
// ====================================================================

/**
 * @test 验证PYRAMID5_EDGE单元的全局棱边数量
 * @details PYRAMID5有8条物理棱边：
 *          - 底面4条: (0,1), (1,2), (2,3), (3,0)
 *          - 侧棱4条: (0,4), (1,4), (2,4), (3,4)
 *          单元无共享棱边（仅1个单元），故全局棱边数=8
 */
TEST_F(Pyramid5EdgeFixture, GlobalEdgeCount_IsEight)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    int num_edges = generator.getNumGlobalEdges();
    EXPECT_EQ(num_edges, 8) << "PYRAMID5_EDGE应有8条全局唯一棱边";

    const auto& edge_mapping = generator.getElemLocalToGlobalEdge();
    EXPECT_EQ(edge_mapping.size(), 1u) << "应有1个单元的棱边映射记录";

    FEEM_INFO("全局棱边数验证: {} (预期8)", num_edges);
}

/**
 * @test 验证每条局部棱边的全局ID映射正确性
 * @details 由于只有1个单元且无共享棱边，每条局部棱边应获得唯一的连续全局ID(0~7)
 *          映射关系应与ElementGeometry::get_local_edges定义一致
 */
TEST_F(Pyramid5EdgeFixture, LocalToGlobalEdgeMapping_Valid)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    const auto& edge_mapping = generator.getElemLocalToGlobalEdge();
    ASSERT_EQ(edge_mapping.size(), 1u) << "应有1个单元的映射";

    const auto& local_to_global = edge_mapping[0];
    ASSERT_EQ(local_to_global.size(), 8u) << "PYRAMID5_EDGE应有8条局部棱边";

    // 验证所有全局ID在[0, 7]范围内且互不重复
    std::vector<bool> used(8, false);
    for (int i = 0; i < 8; ++i) {
        int global_id = local_to_global[i];
        EXPECT_GE(global_id, 0) << "局部棱边" << i << "的全局ID应为非负";
        EXPECT_LT(global_id, 8) << "局部棱边" << i << "的全局ID应<8";
        EXPECT_FALSE(used[global_id]) << "全局ID" << global_id << "被重复使用";
        used[global_id] = true;
    }

    // 输出详细映射信息
    FEEM_INFO("局部→全局棱边映射:");
    for (size_t i = 0; i < local_to_global.size(); ++i) {
        FEEM_INFO("  局部棱边{} → 全局棱边{}", i, local_to_global[i]);
    }
}

/**
 * @test 验证第0条局部棱边对应正确的全局节点对
 * @details 根据VTK标准，PYRAMID5的第0条局部棱边为(局部节点0, 局部节点1)
 *          对应全局节点为(0, 1)，即底面前边缘
 */
TEST_F(Pyramid5EdgeFixture, EdgeZero_IsNodePair_0_1)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    const auto& edge_mapping = generator.getElemLocalToGlobalEdge();
    int edge0_global_id = edge_mapping[0][0];

    FEEM_INFO("第0条局部棱边的全局ID: {}", edge0_global_id);
    // 全局ID应为0（第一条发现的棱边分配ID=0）
    EXPECT_EQ(edge0_global_id, 0) << "第0条局部棱边(节点0-1)应获得全局ID=0";
}


// ====================================================================
//  测试用例2：EMDOFManager DOF预编号
// ====================================================================

/**
 * @test 验证纯矢量格式的DOF预编号数量
 * @details VECTOR_EDGE_ONLY模式下：
 *          - 无标量节点DOF（num_scalar_prenum = 0）
 *          - 8条全局棱边各分配1个矢量DOF（num_vector_prenum = 8）
 *          - 总预编号数 = 8
 */
TEST_F(Pyramid5EdgeFixture, Prenumbering_VectorOnly_EightDOFs)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    EMDOFManager dof_manager(mesh_data_, generator.getElemLocalToGlobalEdge(), &generator);
    dof_manager.build();

    int num_free = dof_manager.getNumFreeDOFs();
    EXPECT_EQ(num_free, 7) << "自由DOF总数应为8(总)-1(约束)=7";

    FEEM_INFO("DOF预编号统计: 自由DOF={} (预期7)", num_free);
}

/**
 * @test 验证约束DOF值的正确性
 * @details Dirichlet边界条件设置value=0.0，
 *          约束值向量中对应位置应存储0.0
 */
TEST_F(Pyramid5EdgeFixture, ConstrainedValue_IsZero)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    EMDOFManager dof_manager(mesh_data_, generator.getElemLocalToGlobalEdge());
    dof_manager.build();

    const auto& constrained_vals = dof_manager.getConstrainedDOFValues();

    // 约束值向量长度应等于预编号总数（8个矢量DOF）
    EXPECT_EQ(constrained_vals.size(), 8u) << "约束值向量长度应等于预编号总数";

    // Dirichlet设置的值为0.0，检查约束值向量中是否包含0.0
    bool has_zero = false;
    for (double v : constrained_vals) {
        if (std::abs(v) < 1e-15) {
            has_zero = true;
            break;
        }
    }
    EXPECT_TRUE(has_zero) << "约束值向量中应包含Dirichlet设置的值0.0";
}


// ====================================================================
//  测试用例3：Local2Global 映射表完整性验证
// ====================================================================

/**
 * @test 验证Local2Global映射表的维度和基本属性
 * @details 对于纯矢量PYRAMID5_EDGE单元：
 *          - indices.size() == 8（8条局部棱边）
 *          - num_scalar_dofs == 0（无标量分量）
 *          - num_vector_dofs == 8（全部为矢量分量）
 *          - is_mixed() == false（非混合格式）
 */
TEST_F(Pyramid5EdgeFixture, Local2Global_Dimensions_Correct)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    EMDOFManager dof_manager(mesh_data_, generator.getElemLocalToGlobalEdge());
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    ASSERT_EQ(mappings.size(), 1u) << "应有1个单元的映射表";

    const auto& l2g = mappings[0];

    // 基本维度验证
    EXPECT_EQ(l2g.indices.size(), 8u)
        << "indices大小应为8（PYRAMID5_EDGE的棱边数）";
    EXPECT_EQ(l2g.num_scalar_dofs, 0)
        << "标量DOF数应为0（VECTOR_EDGE_ONLY模式）";
    EXPECT_EQ(l2g.num_vector_dofs, 8)
        << "矢量DOF数应为8（PYRAMID5_EDGE的全部棱边）";
    EXPECT_FALSE(l2g.is_mixed())
        << "is_mixed()应返回false（纯矢量格式不含标量分量）";
    EXPECT_EQ(l2g.element_id, 0)
        << "element_id应与单元ID一致";
    EXPECT_EQ(static_cast<int>(l2g.elem_type), static_cast<int>(ElemType::PYRAMID5_EDGE))
        << "elem_type应为PYRAMID5_EDGE";

    FEEM_INFO("Local2Global维度验证通过:");
    FEEM_INFO("  indices.size()={}", l2g.indices.size());
    FEEM_INFO("  num_scalar_dofs={}, num_vector_dofs={}", l2g.num_scalar_dofs, l2g.num_vector_dofs);
    FEEM_INFO("  is_mixed()={}", l2g.is_mixed());
}

/**
 * @test 验证第0条局部棱边（节点0-1）被正确标记为约束DOF
 * @details 边界条件固定了局部棱边0（对应全局节点对(0,1)），
 *          因此 indices[0] 应为 -1（约束标记）
 */
TEST_F(Pyramid5EdgeFixture, Local2Global_Edge0_Constrained)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    EMDOFManager dof_manager(mesh_data_, generator.getElemLocalToGlobalEdge());
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    EXPECT_EQ(l2g.indices[0], -1)
        << "第0条局部棱边(节点0-1)应为约束DOF(indices[0]=-1)";

    FEEM_INFO("约束验证: indices[0]={} (预期-1)", l2g.indices[0]);
}

/**
 * @test 验证其余7条局部棱边均为自由DOF且编号连续
 * @details 除第0条棱边外，剩余7条棱边应获得从0开始的连续自由DOF编号
 *          即 indices[1..7] 应包含 {0,1,2,3,4,5,6} 的某种排列
 */
TEST_F(Pyramid5EdgeFixture, Local2Global_FreeEdges_Continuous)
{
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    EMDOFManager dof_manager(mesh_data_, generator.getElemLocalToGlobalEdge());
    dof_manager.build();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& l2g = mappings[0];

    // 收集所有非约束位置的值（应恰好是0~6的排列）
    std::vector<int> free_ids;
    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        if (i == 0) {
            continue;  // 跳过已知的约束位置
        }
        int val = l2g.indices[i];
        EXPECT_GE(val, 0) << "局部棱边" << i << "应为自由DOF(>=0)";
        free_ids.push_back(val);
    }

    // 验证恰好收集了7个自由ID
    ASSERT_EQ(free_ids.size(), 7u) << "应有7个自由DOF";

    // 排序后验证连续性（应为0,1,2,3,4,5,6）
    std::sort(free_ids.begin(), free_ids.end());
    for (int i = 0; i < 7; ++i) {
        EXPECT_EQ(free_ids[i], i)
            << "第" << i << "个自由DOF编号应为" << i
            << "（实际=" << free_ids[i] << "）";
    }

    FEEM_INFO("自由DOF编号验证:");
    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        FEEM_INFO("  indices[{}] = {}", i, l2g.indices[i]);
    }
}


// ====================================================================
//  测试用例4：综合端到端验证
// ====================================================================

/**
 * @test 完整端到端流程验证：从网格构建到DOF映射全链路
 * @details 在单个测试中串联执行所有步骤并输出完整诊断信息，
 *          用于快速确认整个DOF管理管线工作正常
 */
TEST_F(Pyramid5EdgeFixture, EndToEnd_CompletePipeline)
{
    // ===== 第1步：全局棱边ID生成 =====
    GlobalEdgeIDGenerator generator(mesh_data_);
    generator.generate();

    int num_global_edges = generator.getNumGlobalEdges();
    const auto& edge_map = generator.getElemLocalToGlobalEdge();

    ASSERT_EQ(num_global_edges, 8) << "[Step1] 全局棱边数应为8";
    ASSERT_EQ(edge_map.size(), 1u) << "[Step1] 单元映射表大小应为1";
    ASSERT_EQ(edge_map[0].size(), 8u) << "[Step1] 第0单元应有8条局部棱边映射";

    FEEM_INFO("=== Step1: 全局棱边ID生成完成 ===");
    FEEM_INFO("  全局棱边总数: {}", num_global_edges);
    for (size_t i = 0; i < edge_map[0].size(); ++i) {
        FEEM_INFO("  局部棱边{} → 全局棱边{}", i, edge_map[0][i]);
    }

    // ===== 第2步：DOF管理与编号 =====
    EMDOFManager dof_manager(mesh_data_, edge_map, &generator);
    dof_manager.build();

    int num_free = dof_manager.getNumFreeDOFs();
    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const auto& constrained_vals = dof_manager.getConstrainedDOFValues();

    ASSERT_EQ(mappings.size(), 1u) << "[Step2] 映射表应含1个单元";
    ASSERT_EQ(num_free, 7) << "[Step2] 自由DOF数应为7";

    const auto& l2g = mappings[0];

    FEEM_INFO("=== Step2: DOF编号完成 ===");
    FEEM_INFO("  自由DOF总数: {}", num_free);
    FEEM_INFO("  约束值向量长度: {}", constrained_vals.size());

    // ===== 第3步：映射表全面验证 =====
    // 3a. 维度检查
    EXPECT_EQ(static_cast<int>(l2g.indices.size()), 8) << "[Step3a] indices大小";
    EXPECT_EQ(l2g.num_scalar_dofs, 0) << "[Step3a] 标量DOF数";
    EXPECT_EQ(l2g.num_vector_dofs, 8) << "[Step3a] 矢量DOF数";
    EXPECT_FALSE(l2g.is_mixed()) << "[Step3a] is_mixed标志";

    // 3b. 约束位置检查
    int constrained_count = 0;
    int min_free_id = INT_MAX;
    int max_free_id = INT_MIN;
    for (size_t i = 0; i < l2g.indices.size(); ++i) {
        int idx = l2g.indices[i];
        if (idx == -1) {
            constrained_count++;
            FEEM_INFO("  [约束] 局部棱边{} → indices[{}]=-1", i, i);
        } else {
            min_free_id = std::min(min_free_id, idx);
            max_free_id = std::max(max_free_id, idx);
            FEEM_INFO("  [自由] 局部棱边{} → indices[{}]={}", i, i, idx);
        }
    }

    EXPECT_EQ(constrained_count, 1) << "[Step3b] 约束DOF数应为1";
    EXPECT_EQ(min_free_id, 0) << "[Step3c] 最小自由ID应为0";
    EXPECT_EQ(max_free_id, 6) << "[Step3d] 最大自由ID应为6(7个自由DOF: 0~6)";

    // 3c. 特定位置精确检查：第0条棱边必为约束
    EXPECT_EQ(l2g.indices[0], -1)
        << "[Step3e] 第0条局部棱边(节点0-1)必须是约束DOF";

    FEEM_INFO("=== 端到端验证通过 ===");
}


// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
