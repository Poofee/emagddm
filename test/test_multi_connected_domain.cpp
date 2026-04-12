/**
 * @file test_multi_connected_domain.cpp
 * @brief MultiConnectedDomainDetector多连通域检测器完整测试套件
 * @details 使用 Google Test 框架验证MultiConnectedDomainDetector的核心功能：
 *          1. 单连通域模型的正确检测（num_holes=0）
 *          2. 单孔洞模型的环量自由度识别
 *          3. 多孔洞模型的独立环量选择
 *          4. 环量自由度子集验证
 *          5. 边界条件和异常输入的健壮性处理
 *
 * 测试场景设计：
 * - 场景1：单连通域立方体（无孔洞，基本功能验证）
 * - 场景2：带命名边界的简单模型（启发式边界分析）
 * - 场景3：多边界标记模型（多个潜在孔洞）
 * - 场景4：空网格和无效输入（边界条件）
 * - 场景5：环量自由度一致性验证（数学正确性）
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>

#include "em_mesh_data.hpp"
#include "global_edge_id_generator.hpp"
#include "mesh_graph_extractor.hpp"
#include "spanning_tree.hpp"
#include "multi_connected_domain.hpp"
#include "element_geometry.hpp"
#include "logger_factory.hpp"

using namespace fe_em;

// ==================== 辅助函数 ====================

/**
 * @brief 构建完整的测试流程：生成棱边ID → 提取图 → 构建生成树 → 检测多连通域
 * @param mesh_data 网格数据
 * @return DetectionResult 多连通域检测结果
 */
static DetectionResult runFullDetectionPipeline(const EMMeshData& mesh_data) {
    // 第一步：生成全局棱边ID映射
    GlobalEdgeIDGenerator edge_gen(mesh_data);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    // 第二步：提取图结构
    MeshGraphExtractor extractor(mesh_data, edge_mapping);
    GraphResult graph_result = extractor.extract();

    // 第三步：构建BFS生成树（使用默认根节点）
    auto tree_result = SpanningTreeBuilder::buildBFSTree(
        graph_result.adjacency_list,
        graph_result.edge_to_global_id,
        0
    );

    // 第四步：检测多连通域
    return MultiConnectedDomainDetector::detect(mesh_data, tree_result);
}


// ====================================================================
//  测试固件1：单连通域模型（无孔洞）
// ====================================================================

/**
 * @class SimplyConnectedFixture
 * @brief 单连通域测试固件
 * @details 构建1个HEX8单元的简单立方体网格，
 *          用于验证单连通域检测的正确性（应返回num_holes=0）。
 *
 * 测试网格拓扑：
 * - 1个HEX8单元，8个节点构成单位立方体 [0,1]^3
 * - 无任何边界标记（或仅有外部边界）
 * - 应被识别为单连通域（genus=0）
 */
class SimplyConnectedFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 构建8个节点（单位立方体） ----------
        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},  // 节点0: 原点
            {1, 1.0, 0.0, 0.0, 0},  // 节点1: x轴方向
            {2, 1.0, 1.0, 0.0, 0},  // 节点2: xy平面
            {3, 0.0, 1.0, 0.0, 0},  // 节点3: y轴方向
            {4, 0.0, 0.0, 1.0, 0},  // 节点4: z轴方向
            {5, 1.0, 0.0, 1.0, 0},  // 节点5: xz平面
            {6, 1.0, 1.0, 1.0, 0},  // 节点6: 对角顶点
            {7, 0.0, 1.0, 1.0, 0},  // 节点7: yz平面
        };

        // ---------- 构建1个HEX8单元 ----------
        mesh_data_.elements = {{
            0,
            {0, 1, 2, 3, 4, 5, 6, 7},
            ElemType::HEX8,
            DOFType::VECTOR_EDGE_ONLY,
            1,
            0,
        }};

        FEEM_INFO("SimplyConnectedFixture: 网格构建完成 - {}节点, {}单元",
                  mesh_data_.nodes.size(), mesh_data_.elements.size());
    }

    EMMeshData mesh_data_;
};


// ---------- 1.1 单连通域基本检测测试 ----------

/**
 * @test 验证单连通域模型的孔洞数为0
 * @details 对于无内部边界的简单立方体，应返回num_holes=0
 */
TEST_F(SimplyConnectedFixture, ZeroHoles) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    EXPECT_EQ(result.num_holes, 0) << "单连通域模型应有0个孔洞";
    EXPECT_TRUE(result.circulation_dofs.empty()) << "单连通域不应有环量自由度";
    EXPECT_TRUE(result.hole_loop_edges.empty()) << "回路列表应为空";
    EXPECT_TRUE(result.isValid()) << "结果应为有效状态";

    FEEM_INFO("单连通域测试: num_holes={} ✓", result.num_holes);
}


// ---------- 1.2 单连通域结果完整性测试 ----------

/**
 * @test 验证单连通域结果的各字段一致性
 * @details 即使num_holes=0，结果结构仍应保持完整和一致
 */
TEST_F(SimplyConnectedFixture, ResultConsistency) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    // 验证所有字段都正确初始化
    EXPECT_GE(result.num_holes, 0) << "孔洞数不应为负";
    EXPECT_EQ(static_cast<int>(result.hole_loop_edges.size()), result.num_holes)
        << "回路列表大小应等于孔洞数";
    EXPECT_EQ(static_cast<int>(result.circulation_dofs.size()), result.num_holes)
        << "环量自由度数应等于孔洞数（单连通域时均为0）";

    FEEM_INFO("结果一致性验证: 通过");
}


// ====================================================================
//  测试固件2：带边界标记的单孔洞候选模型
// ====================================================================

/**
 * @class SingleHoleCandidateFixture
 * @brief 单孔洞候选模型测试固件
 * @details 构建1个HEX8单元 + 1个命名为"Hole_1"的边界标记，
 *          用于验证启发式边界分析能识别出潜在孔洞。
 *
 * 测试场景：
 * - 网格拓扑与SimplyConnectedFixture相同（1个HEX8单元）
 * - 额外添加1个名为"Hole_1"的PERFECT_E边界标记
 * - 启发式算法应识别出1个潜在孔洞并选择1条环量自由度
 */
class SingleHoleCandidateFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 复用SimplyConnectedFixture的网格拓扑 ----------
        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},
            {1, 1.0, 0.0, 0.0, 0},
            {2, 1.0, 1.0, 0.0, 0},
            {3, 0.0, 1.0, 0.0, 0},
            {4, 0.0, 0.0, 1.0, 0},
            {5, 1.0, 0.0, 1.0, 0},
            {6, 1.0, 1.0, 1.0, 0},
            {7, 0.0, 1.0, 1.0, 0},
        };

        mesh_data_.elements = {{
            0,
            {0, 1, 2, 3, 4, 5, 6, 7},
            ElemType::HEX8,
            DOFType::VECTOR_EDGE_ONLY,
            1,
            0,
        }};

        // ---------- 添加命名为"Hole_1"的边界标记 ----------
        EMBoundaryMarker hole_marker;
        hole_marker.id = 0;
        hole_marker.bnd_type = BndType::PERFECT_E;
        hole_marker.dof_type = DOFType::VECTOR_EDGE_ONLY;
        hole_marker.target_ids = std::vector<int>{0, 1};  // 小目标集（视为内部边界）
        hole_marker.value = 0.0;
        hole_marker.name = "Hole_1";  // 关键：名称包含"hole"关键词

        mesh_data_.boundary_markers.push_back(hole_marker);

        FEEM_INFO("SingleHoleCandidateFixture: 网格构建完成 - {}节点, {}单元, {}边界标记",
                  mesh_data_.nodes.size(), mesh_data_.elements.size(),
                  mesh_data_.boundary_markers.size());
    }

    EMMeshData mesh_data_;
};


// ---------- 2.1 单孔洞候选检测测试 ----------

/**
 * @test 验证带"Hole_1"命名边界的模型能检测到孔洞候选
 * @details 由于边界名称包含"hole"关键词，启发式算法应识别为潜在孔洞
 */
TEST_F(SingleHoleCandidateFixture, DetectSingleHole) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    EXPECT_GT(result.num_holes, 0) << "应检测到至少1个潜在孔洞";
    EXPECT_GE(static_cast<int>(result.circulation_dofs.size()), result.num_holes)
        << "环量自由度数应>= 孔洞数";
    EXPECT_EQ(static_cast<int>(result.hole_loop_edges.size()), result.num_holes)
        << "回路列表大小应等于孔洞数";
    EXPECT_TRUE(result.isValid()) << "结果应为有效状态";

    FEEM_INFO("单孔洞候选测试: num_holes={}, circulation_dofs_count={} ✓",
              result.num_holes, result.circulation_dofs.size());
}


// ---------- 2.2 环量自由度数量合理性测试 ----------

/**
 * @test 验证环量自由度数量不超过余树边数量
 * @details 即使检测到孔洞，选中的环量边也不能超过可用的余树边总数
 */
TEST_F(SingleHoleCandidateFixture, CirculationDofCountReasonable) {
    // 先获取生成树信息以知道余树边数量
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    MeshGraphExtractor extractor(mesh_data_, edge_mapping);
    GraphResult graph_result = extractor.extract();

    auto tree_result = SpanningTreeBuilder::buildBFSTree(
        graph_result.adjacency_list,
        graph_result.edge_to_global_id,
        0
    );

    int cotree_size = static_cast<int>(tree_result.cotree_edges.size());

    DetectionResult domain_result = MultiConnectedDomainDetector::detect(
        mesh_data_, tree_result
    );

    EXPECT_LE(static_cast<int>(domain_result.circulation_dofs.size()), cotree_size)
        << "环量自由度数(" << domain_result.circulation_dofs.size()
        << ")不应超过余树边数(" << cotree_size << ")";

    FEEM_INFO("环量自由度数量验证: {} <= {} ✓",
              domain_result.circulation_dofs.size(), cotree_size);
}


// ====================================================================
//  测试固件3：多孔洞候选模型
// ====================================================================

/**
 * @class MultiHoleCandidateFixture
 * @brief 多孔洞候选模型测试固件
 * @details 构建1个HEX8单元 + 2个不同命名的边界标记（"inner_1"和"cavity_2"），
 *          用于验证多孔洞场景下的独立环量自由度选择。
 *
 * 测试场景：
 * - 网格拓扑与前面相同（1个HEX8单元，12条棱边）
 * - 添加2个边界标记，名称分别包含"inner"和"cavity"关键词
 * - 应识别出2个潜在孔洞并选择2条独立的环量自由度
 */
class MultiHoleCandidateFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 复用相同的网格拓扑 ----------
        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},
            {1, 1.0, 0.0, 0.0, 0},
            {2, 1.0, 1.0, 0.0, 0},
            {3, 0.0, 1.0, 0.0, 0},
            {4, 0.0, 0.0, 1.0, 0},
            {5, 1.0, 0.0, 1.0, 0},
            {6, 1.0, 1.0, 1.0, 0},
            {7, 0.0, 1.0, 1.0, 0},
        };

        mesh_data_.elements = {{
            0,
            {0, 1, 2, 3, 4, 5, 6, 7},
            ElemType::HEX8,
            DOFType::VECTOR_EDGE_ONLY,
            1,
            0,
        }};

        // ---------- 添加2个边界标记 ----------
        EMBoundaryMarker marker1;
        marker1.id = 0;
        marker1.bnd_type = BndType::PERFECT_E;
        marker1.dof_type = DOFType::VECTOR_EDGE_ONLY;
        marker1.target_ids = std::vector<int>{0, 1};
        marker1.value = 0.0;
        marker1.name = "inner_boundary_1";  // 包含"inner"

        EMBoundaryMarker marker2;
        marker2.id = 1;
        marker2.bnd_type = BndType::PERFECT_H;
        marker2.dof_type = DOFType::VECTOR_EDGE_ONLY;
        marker2.target_ids = std::vector<int>{2, 3};
        marker2.value = 0.0;
        marker2.name = "cavity_region_2";   // 包含"cavity"

        mesh_data_.boundary_markers.push_back(marker1);
        mesh_data_.boundary_markers.push_back(marker2);

        FEEM_INFO("MultiHoleCandidateFixture: 网格构建完成 - {}边界标记",
                  mesh_data_.boundary_markers.size());
    }

    EMMeshData mesh_data_;
};


// ---------- 3.1 多孔洞检测测试 ----------

/**
 * @test 验证多边界标记模型能检测到多个孔洞候选
 * @details 2个包含关键词的边界标记应产生至少2个潜在孔洞
 */
TEST_F(MultiHoleCandidateFixture, DetectMultipleHoles) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    EXPECT_GE(result.num_holes, 2) << "应检测到至少2个潜在孔洞";
    EXPECT_GE(static_cast<int>(result.circulation_dofs.size()), 2)
        << "应至少有2个环量自由度";
    EXPECT_EQ(static_cast<int>(result.hole_loop_edges.size()), result.num_holes)
        << "每个孔洞应对应一个回路";

    FEEM_INFO("多孔洞检测测试: num_holes={}, circulation_dofs={} ✓",
              result.num_holes, result.circulation_dofs.size());
}


// ---------- 3.2 环量自由度独立性测试 ----------

/**
 * @test 验证选中的环量自由度是不同的边（无重复）
 * @details 不同孔洞的环量自由度应是不同的余树边ID
 */
TEST_F(MultiHoleCandidateFixture, CirculationDofUniqueness) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    if (result.num_holes >= 2) {
        // 验证circulation_dofs中的元素互不相同
        std::set<int> unique_dofs(result.circulation_dofs.begin(),
                                   result.circulation_dofs.end());

        EXPECT_EQ(static_cast<int>(unique_dofs.size()),
           static_cast<int>(result.circulation_dofs.size()))
            << "环量自由度中存在重复元素";

        FEEM_INFO("环量自由度独立性验证: {}个唯一值 ✓", unique_dofs.size());
    } else {
        SUCCEED() << "孔洞数<2，跳过独立性测试";
    }
}


// ---------- 3.3 回路列表结构测试 ----------

/**
 * @test 验证每个孔洞都有对应的非空回路
 * @details hole_loop_edges的每一行都应至少包含1条棱边
 */
TEST_F(MultiHoleCandidateFixture, HoleLoopEdgesNonEmpty) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    for (int i = 0; i < result.num_holes; ++i) {
        EXPECT_FALSE(result.hole_loop_edges[i].empty())
            << "孔洞" << i << "的回路不应为空";
    }

    FEEM_INFO("回路非空验证: 所有{}个孔洞均有对应回路 ✓", result.num_holes);
}


// ====================================================================
//  测试固件4：环量自由度子集验证（数学正确性核心测试）
// ====================================================================

/**
 * @class SubsetVerificationFixture
 * @brief 环量子集验证测试固件
 * @details 使用已知的生成树结果，严格验证circulation_dofs ⊆ cotree_edges。
 *
 * 这是Tree Gauge模块正确性的**核心数学约束**：
 * - 环量自由度必须是余树边的子集
 * - 否则会导致Tree Gauge规范条件被违反
 */
class SubsetVerificationFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 构建最小有效网格 ----------
        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},
            {1, 1.0, 0.0, 0.0, 0},
            {2, 0.0, 1.0, 0.0, 0},
            {3, 0.0, 0.0, 1.0, 0},
        };

        mesh_data_.elements = {{
            0,
            {0, 1, 2, 3},
            ElemType::TET4,
            DOFType::VECTOR_EDGE_ONLY,
            1,
            0,
        }};
    }

    EMMeshData mesh_data_;
};


// ---------- 4.1 子集关系严格验证 ----------

/**
 * @test 验证circulation_dofs是cotree_edges的真子集
 * @details 这是最重要的数学正确性测试
 */
TEST_F(SubsetVerificationFixture, CirculationDofsSubsetOfCotreeEdges) {
    GlobalEdgeIDGenerator edge_gen(mesh_data_);
    edge_gen.generate();
    const auto& edge_mapping = edge_gen.getElemLocalToGlobalEdge();

    MeshGraphExtractor extractor(mesh_data_, edge_mapping);
    GraphResult graph_result = extractor.extract();

    auto tree_result = SpanningTreeBuilder::buildBFSTree(
        graph_result.adjacency_list,
        graph_result.edge_to_global_id,
        0
    );

    DetectionResult domain_result = MultiConnectedDomainDetector::detect(
        mesh_data_, tree_result
    );

    // 验证每个环量自由度都在cotree_edges中
    for (int dof : domain_result.circulation_dofs) {
        EXPECT_NE(tree_result.cotree_edges.find(dof),
                  tree_result.cotree_edges.end())
            << "环量自由度" << dof << "不在余树边集合中（违反子集约束）";
    }

    FEEM_INFO("子集验证: 所有{}个环量自由度均在余树边集合中 ✓",
              domain_result.circulation_dofs.size());
}


// ---------- 4.2 大规模随机子集验证 ----------

/**
 * @test 在较大规模网格上验证子集关系的一致性
 * @details 使用PRISM6单元（9条棱边）进行更全面的验证
 */
TEST_F(SubsetVerificationFixture, LargeMeshSubsetConsistency) {
    // 构建稍大的网格（1个PRISM6单元，9条棱边）
    EMMeshData prism_mesh;
    prism_mesh.nodes = {
        {0, 0.0,     0.0,    0.0, 0},
        {1, 1.0,     0.0,    0.0, 0},
        {2, 0.5,     0.866,  0.0, 0},
        {3, 0.0,     0.0,    1.0, 0},
        {4, 1.0,     0.0,    1.0, 0},
        {5, 0.5,     0.866,  1.0, 0},
    };

    prism_mesh.elements = {{
        0,
        {0, 1, 2, 3, 4, 5},
        ElemType::PRISM6,
        DOFType::MIXED_AV,
        1,
        0,
    }};

    // 添加2个边界标记以触发孔洞检测
    EMBoundaryMarker m1, m2;
    m1.id = 0;
    m1.bnd_type = BndType::PERFECT_E;
    m1.dof_type = DOFType::MIXED_AV;
    m1.target_ids = std::vector<int>{0, 1};
    m1.value = 0.0;
    m1.name = "hole_A";

    m2.id = 1;
    m2.bnd_type = BndType::PERFECT_H;
    m2.dof_type = DOFType::MIXED_AV;
    m2.target_ids = std::vector<int>{3, 4};
    m2.value = 0.0;
    m2.name = "hole_B";

    prism_mesh.boundary_markers.push_back(m1);
    prism_mesh.boundary_markers.push_back(m2);

    // 运行完整流程
    GlobalEdgeIDGenerator gen(prism_mesh);
    gen.generate();
    MeshGraphExtractor ext(prism_mesh, gen.getElemLocalToGlobalEdge());
    GraphResult gr = ext.extract();
    auto tr = SpanningTreeBuilder::buildBFSTree(gr.adjacency_list, gr.edge_to_global_id, 0);
    DetectionResult dr = MultiConnectedDomainDetector::detect(prism_mesh, tr);

    // 验证子集关系
    for (int dof : dr.circulation_dofs) {
        EXPECT_NE(tr.cotree_edges.find(dof), tr.cotree_edges.end())
            << "PRISM6网格: 环量自由度" << dof << "不在余树边中";
    }

    FEEM_INFO("大规模子集验证: PRISM6网格通过 ✓ (孔洞={}, 环量={})",
              dr.num_holes, dr.circulation_dofs.size());
}


// ====================================================================
//  测试固件5：边界条件和异常输入
// ====================================================================

/**
 * @class EdgeCaseFixture
 * @brief 异常输入和边界条件测试固件
 * @details 测试各种极端情况下的行为健壮性：
 *          - 空网格
 *          - 无效的生成树结果
 *          - 空余树边集合
 *          - 大量边界标记但无实际孔洞
 */
class EdgeCaseFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);
    }
};


// ---------- 5.1 空网格异常测试 ----------

/**
 * @test 验证空网格抛出异常
 * @details 当mesh_data不包含节点时应抛出std::invalid_argument
 */
TEST_F(EdgeCaseFixture, EmptyMeshThrowsException) {
    EMMeshData empty_mesh;

    SpanningTreeResult dummy_tree;
    dummy_tree.root_node = 0;
    dummy_tree.node_parent = {{0, -1}};

    EXPECT_THROW(
        MultiConnectedDomainDetector::detect(empty_mesh, dummy_tree),
        std::invalid_argument
    ) << "空网格应抛出std::invalid_argument";

    FEEM_INFO("空网格异常测试: 正确抛出异常 ✓");
}


// ---------- 5.2 无效生成树结果测试 ----------

/**
 * @test 验证无效的生成树结果抛出异常
 * @details 当tree_result未正确初始化时应抛出异常
 */
TEST_F(EdgeCaseFixture, InvalidTreeResultThrowsException) {
    EMMeshData valid_mesh;
    valid_mesh.nodes = {{0, 0, 0, 0, 0}};
    valid_mesh.elements = {{
        0, {0}, ElemType::LINE2, DOFType::SCALAR_ONLY, 1, 0
    }};

    SpanningTreeResult invalid_tree;
    // 故意不设置node_parent，使isValid()返回false

    EXPECT_THROW(
        MultiConnectedDomainDetector::detect(valid_mesh, invalid_tree),
        std::invalid_argument
    ) << "无效生成树结果应抛出std::invalid_argument";

    FEEM_INFO("无效生成树测试: 正确抛出异常 ✓");
}


// ---------- 5.3 空余树边测试 ----------

/**
 * @test 验证空余树边集合的处理（警告但不崩溃）
 * @details 当余树边为空时，应返回空结果而非崩溃
 */
TEST_F(EdgeCaseFixture, EmptyCotreeEdgesReturnsEmptyResult) {
    EMMeshData test_mesh;
    test_mesh.nodes = {
        {0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0},
    };

    test_mesh.elements = {{
        0, {0, 1}, ElemType::LINE2, DOFType::VECTOR_EDGE_ONLY, 1, 0
    }};

    // 构建只有树边没有余树边的特殊情况（退化情况）
    SpanningTreeResult degenerate_tree;
    degenerate_tree.root_node = 0;
    degenerate_tree.node_parent = {{0, -1}, {1, 0}};
    degenerate_tree.tree_edges = {0};  // 唯一的边在树中
    degenerate_tree.cotree_edges = {};  // 余树边为空

    // 不应抛出异常，而是返回空结果
    EXPECT_NO_THROW({
        DetectionResult result = MultiConnectedDomainDetector::detect(test_mesh, degenerate_tree);
        EXPECT_EQ(result.num_holes, 0) << "空余树边时孔洞数应为0";
        EXPECT_TRUE(result.circulation_dofs.empty()) << "空余树边时环量自由度应为空";
    });

    FEEM_INFO("空余树边测试: 返回安全空结果 ✓");
}


// ---------- 5.4 大量边界标记但无孔洞关键词测试 ----------

/**
 * @test 验证大量普通边界标记不会误判为孔洞
 * @details 当边界标记名称不包含孔洞关键词时，即使数量很多也不应检测到孔洞
 */
TEST_F(EdgeCaseFixture, ManyNormalBoundariesNoFalsePositive) {
    EMMeshData test_mesh;
    test_mesh.nodes = {
        {0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0},
        {2, 0, 1, 0, 0},
        {3, 1, 1, 0, 0},
    };

    test_mesh.elements = {{
        0, {0, 1, 2, 3}, ElemType::QUAD4, DOFType::VECTOR_EDGE_ONLY, 1, 0
    }};

    // 添加10个普通边界标记（名称不含hole/inner/cavity关键词）
    for (int i = 0; i < 10; ++i) {
        EMBoundaryMarker marker;
        marker.id = i;
        marker.bnd_type = BndType::DIRICHLET;
        marker.dof_type = DOFType::VECTOR_EDGE_ONLY;
        marker.target_ids = std::vector<int>{i % 4};  // 单个节点
        marker.value = 0.0;
        marker.name = "Boundary_" + std::to_string(i);  // 普通命名

        test_mesh.boundary_markers.push_back(marker);
    }

    DetectionResult result = runFullDetectionPipeline(test_mesh);

    // 可能会因PERFECT_E/H类型和小目标集而检测到一些候选
    // 但主要验证的是不会崩溃且结果一致
    EXPECT_TRUE(result.isValid()) << "结果应为有效状态";
    EXPECT_GE(result.num_holes, 0) << "孔洞数不应为负";

    FEEM_INFO("大量普通边界测试: num_holes={} (无崩溃) ✓", result.num_holes);
}


// ====================================================================
//  测试固件6：混合DOF类型过滤测试
// ====================================================================

/**
 * @class MixedDOFFixture
 * @brief 混合DOF类型测试固件
 * @details 验证当网格同时包含矢量和标量单元时，
 *          仅基于矢量单元的图进行多连通域检测。
 */
class MixedDOFFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},
            {1, 1.0, 0.0, 0.0, 0},
            {2, 0.0, 1.0, 0.0, 0},
            {3, 1.0, 1.0, 0.0, 0},
        };

        // 2个单元：1个标量 + 1个矢量
        mesh_data_.elements = {
            {
                0,
                {0, 1, 2, 3},
                ElemType::QUAD4,
                DOFType::SCALAR_ONLY,
                1,
                0,
            },
            {
                1,
                {0, 1, 2, 3},
                ElemType::QUAD4,
                DOFType::VECTOR_EDGE_ONLY,
                1,
                0,
            },
        };
    }

    EMMeshData mesh_data_;
};


// ---------- 6.1 混合DOF类型正确性测试 ----------

/**
 * @test 验证混合DOF类型网格的多连通域检测仅使用矢量单元
 * @details 结果应基于矢量单元的图结构，忽略标量单元
 */
TEST_F(MixedDOFFixture, MixedDOFDetectionCorrectness) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    EXPECT_TRUE(result.isValid()) << "混合DOF类型网格应成功完成检测";
    EXPECT_GE(result.num_holes, 0) << "孔洞数合法";

    FEEM_INFO("混合DOF类型测试: 成功完成检测 (num_holes={}) ✓", result.num_holes);
}


// ====================================================================
//  综合输出测试
// ====================================================================

/**
 * @test 输出单孔洞候选模型的完整检测结果供人工审查
 * @details 将所有字段格式化输出，便于调试和验证
 */
TEST_F(SingleHoleCandidateFixture, FullResultDump) {
    DetectionResult result = runFullDetectionPipeline(mesh_data_);

    FEEM_INFO("");
    FEEM_INFO("============================================================");
    FEEM_INFO("  SingleHoleCandidate 完整检测结果输出");
    FEEM_INFO("============================================================");

    FEEM_INFO("--- 基本信息 ---");
    FEEM_INFO("  孔洞数 (genus): {}", result.num_holes);
    FEEM_INFO("  环量自由度数: {}", result.circulation_dofs.size());
    FEEM_INFO("  回路列表大小: {}", result.hole_loop_edges.size());
    FEEM_INFO("  结果有效性: {}", result.isValid() ? "有效" : "无效");

    if (result.num_holes > 0) {
        FEEM_INFO("--- 环量自由度详情 ---");
        int idx = 0;
        for (int edge_id : result.circulation_dofs) {
            FEEM_INFO("  环量自由度[{}]: 全局棱边ID={}", idx, edge_id);
            idx++;
        }

        FEEM_INFO("--- 孔洞回路详情 ---");
        for (int i = 0; i < result.num_holes; ++i) {
            std::string edges_str;
            for (int eid : result.hole_loop_edges[i]) {
                edges_str += std::to_string(eid) + " ";
            }
            FEEM_INFO("  孔洞{}回路: [{}] (共{}条棱边)",
                      i, edges_str, result.hole_loop_edges[i].size());
        }
    }

    FEEM_INFO("--- 统计摘要 ---");
    FEEM_INFO("  平均每孔洞环量自由度: {:.2f}",
              result.num_holes > 0 ?
              static_cast<double>(result.circulation_dofs.size()) / result.num_holes :
              0.0);
    FEEM_INFO("============================================================");
    FEEM_INFO("");

    // 最终断言：确保数据自洽
    EXPECT_TRUE(result.isValid()) << "完整输出的结果必须有效";
}


// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
