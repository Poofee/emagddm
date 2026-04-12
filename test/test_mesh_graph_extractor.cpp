/**
 * @file test_mesh_graph_extractor.cpp
 * @brief MeshGraphExtractor图结构提取器完整测试套件
 * @details 使用 Google Test 框架验证MeshGraphExtractor的核心功能：
 *          1. PRISM6单元（6节点9棱边）的图提取正确性
 *          2. HEX8单元的多单元邻接表合并与棱边去重
 *          3. 空网格和异常输入的健壮性处理
 *          4. 节点索引映射的正确性（非连续ID处理）
 *
 * 测试场景设计：
 * - 场景1：单个PRISM6 MIXED_AV单元（基本功能验证）
 * - 场景2：两个相邻HEX8 VECTOR_EDGE_ONLY单元（棱边共享与去重）
 * - 场景3：空网格和纯标量网格（边界条件）
 * - 场景4：输入尺寸不匹配（异常处理）
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
#include "element_geometry.hpp"
#include "logger_factory.hpp"

using namespace fe_em;

// ==================== 辅助函数 ====================

/**
 * @brief 构建GlobalEdgeIDGenerator并生成映射表的便捷函数
 * @param mesh_data 网格数据
 * @return std::vector<std::vector<int>> 全局棱边ID映射表
 */
static std::vector<std::vector<int>> generateEdgeMapping(const EMMeshData& mesh_data) {
    GlobalEdgeIDGenerator edge_gen(mesh_data);
    edge_gen.generate();
    return edge_gen.getElemLocalToGlobalEdge();
}


// ====================================================================
//  测试固件1：PRISM6混合A-V格式单元
// ====================================================================

/**
 * @class Prism6GraphFixture
 * @brief PRISM6单元图提取测试固件
 * @details 构建1个PRISM6 MIXED_AV单元的完整网格数据，
 *          用于验证基本图提取功能的正确性。
 *
 * 测试网格拓扑：
 * - 1个PRISM6单元，6个节点构成单位三棱柱
 * - 底面(z=0): 节点0(0,0,0), 节点1(1,0,0), 节点2(0.5,0.866,0)
 * - 顶面(z=1): 节点3(0,0,1), 节点4(1,0,1), 节点5(0.5,0.866,1)
 */
class Prism6GraphFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 构建6个节点（单位三棱柱） ----------
        mesh_data_.nodes = {
            {0, 0.0,     0.0,    0.0, 0},
            {1, 1.0,     0.0,    0.0, 0},
            {2, 0.5,     0.866,  0.0, 0},
            {3, 0.0,     0.0,    1.0, 0},
            {4, 1.0,     0.0,    1.0, 0},
            {5, 0.5,     0.866,  1.0, 0},
        };

        // ---------- 构建1个PRISM6单元（MIXED_AV格式） ----------
        mesh_data_.elements = {{
            0,
            {0, 1, 2, 3, 4, 5},
            ElemType::PRISM6,
            DOFType::MIXED_AV,
            1,
            0,
        }};

        // 生成全局棱边映射
        edge_mapping_ = generateEdgeMapping(mesh_data_);

        FEEM_INFO("Prism6GraphFixture: 网格构建完成 - {}节点, {}单元",
                  mesh_data_.nodes.size(), mesh_data_.elements.size());
    }

    EMMeshData mesh_data_;
    std::vector<std::vector<int>> edge_mapping_;
};


// ---------- 1.1 基本统计信息测试 ----------

/**
 * @test 验证PRISM6单元的图结构节点数
 * @details PRISM6有6个节点，所有节点都应出现在图中
 */
TEST_F(Prism6GraphFixture, NodeCount) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    EXPECT_EQ(result.num_nodes, 6) << "PRISM6单元应有6个节点";
    EXPECT_TRUE(result.isValid()) << "结果应为有效状态";

    FEEM_INFO("节点数验证: 预期=6, 实际={}", result.num_nodes);
}

/**
 * @test 验证PRISM6单元的图结构棱边数
 * @details PRISM6三棱柱有9条物理棱边（3底边 + 3顶边 + 3侧棱）
 */
TEST_F(Prism6GraphFixture, EdgeCount) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    EXPECT_EQ(result.num_edges, 9) << "PRISM6单元应有9条唯一棱边";
    EXPECT_EQ(static_cast<int>(result.edge_to_global_id.size()), 9)
        << "edge_to_global_id映射表大小应为9";

    FEEM_INFO("棱边数验证: 预期=9, 实际={}", result.num_edges);
}


// ---------- 1.2 邻接表结构测试 ----------

/**
 * @test 验证邻接表维度
 * @details adjacency_list应有6个条目（每个节点一个）
 */
TEST_F(Prism6GraphFixture, AdjacencyListDimension) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    ASSERT_EQ(static_cast<int>(result.adjacency_list.size()), 6)
        << "邻接表应有6个条目（对应6个节点）";

    FEEM_INFO("邻接表维度: {}", result.adjacency_list.size());
}

/**
 * @test 验证每个节点的邻居数量合理性
 * @details PRISM6单元中：
 * - 每个底面/顶面节点连接3条棱边（度数=3）
 * - 无孤立节点（所有节点的邻居数>0）
 */
TEST_F(Prism6GraphFixture, NodeDegrees) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    for (int i = 0; i < result.num_nodes; ++i) {
        int degree = static_cast<int>(result.adjacency_list[i].size());
        EXPECT_GT(degree, 0) << "节点" << i << "不应为孤立节点";

        // PRISM6中每个节点至少连接2条棱边，最多连接3条
        EXPECT_GE(degree, 2) << "节点" << i << "的度数不应<2";
        EXPECT_LE(degree, 3) << "节点" << i << "的度数不应>3";
    }

    FEEM_INFO("=== 节点度数分布 ===");
    for (int i = 0; i < result.num_nodes; ++i) {
        int global_id = result.index_to_global_node[i];
        FEEM_INFO("  节点{}(全局ID={}): 度数={}",
                  i, global_id, result.adjacency_list[i].size());
    }
}

/**
 * @test 验证邻接表的无向对称性
 * @details 若节点i在节点j的邻居列表中，则节点j也应在节点i的邻居列表中
 */
TEST_F(Prism6GraphFixture, UndirectedSymmetry) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    for (int i = 0; i < result.num_nodes; ++i) {
        for (int neighbor : result.adjacency_list[i]) {
            bool found = false;
            for (int neighbor_of_neighbor : result.adjacency_list[neighbor]) {
                if (neighbor_of_neighbor == i) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found)
                << "无向对称性破坏: 节点" << i << "指向节点" << neighbor
                << ", 但反向链接不存在";
        }
    }

    FEEM_INFO("无向对称性验证: 通过");
}


// ---------- 1.3 棱边映射表测试 ----------

/**
 * @test 验证edge_to_global_id的键值范围
 * @details 所有全局棱边ID应在[0, 8]范围内且互不重复
 */
TEST_F(Prism6GraphFixture, EdgeMappingValueRange) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    std::vector<bool> id_found(9, false);
    for (const auto& [edge_key, global_id] : result.edge_to_global_id) {
        EXPECT_GE(global_id, 0) << "全局棱边ID不应为负";
        EXPECT_LT(global_id, 9) << "全局棱边ID应<9";
        EXPECT_FALSE(id_found[global_id])
            << "全局棱边ID=" << global_id << "重复出现";
        id_found[global_id] = true;
    }

    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(id_found[i]) << "全局棱边ID=" << i << "缺失";
    }

    FEEM_INFO("棱边ID连续性验证: 通过（覆盖[0, 8]全部9个ID）");
}

/**
 * @test 验证edge_to_global_id的键为排序后的节点对
 * @details 每个键的first <= second（升序排列）
 */
TEST_F(Prism6GraphFixture, EdgeKeyOrdering) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    for (const auto& [edge_key, global_id] : result.edge_to_global_id) {
        EXPECT_LE(edge_key.first, edge_key.second)
            << "棱边键(" << edge_key.first << ", " << edge_key.second
            << ")未按升序排列";
    }

    FEEM_INFO("棱边键排序验证: 通过（所有键均为升序）");
}


// ---------- 1.4 节点索引映射测试 ----------

/**
 * @test 验证全局节点ID到紧凑索引的双向映射一致性
 * @details global_node_to_index和index_to_global_node应为互逆映射
 */
TEST_F(Prism6GraphFixture, NodeIndexMappingConsistency) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    ASSERT_EQ(static_cast<int>(result.global_node_to_index.size()), 6)
        << "global_node_to_index应有6个条目";
    ASSERT_EQ(static_cast<int>(result.index_to_global_node.size()), 6)
        << "index_to_global_node应有6个元素";

    for (int i = 0; i < result.num_nodes; ++i) {
        int global_id = result.index_to_global_node[i];
        auto it = result.global_node_to_index.find(global_id);
        ASSERT_NE(it, result.global_node_to_index.end())
            << "全局节点ID=" << global_id << "未在global_node_to_index中找到";
        EXPECT_EQ(it->second, i)
            << "双向映射不一致: index=" << i
            << " → global=" << global_id
            << " → index=" << it->second;
    }

    FEEM_INFO("节点索引映射一致性验证: 通过");
}


// ====================================================================
//  测试固件2：两个相邻HEX8单元（棱边共享测试）
// ====================================================================

/**
 * @class Hex8MultiElementFixture
 * @brief HEX8多单元图提取测试固件
 * @details 构建2个共享一个面的HEX8 VECTOR_EDGE_ONLY单元，
 *          用于验证棱边去重和邻接表合并的正确性。
 *
 * 测试网格拓扑：
 * - 单元0: 节点0-7（标准HEX8）
 * - 单元1: 节点4-11（与单元0共享节点4,5,6,7构成的面）
 * - 总共12个节点，共享4个节点
 * - 共享面的4条棱边应只记录一次
 */
class Hex8MultiElementFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);

        // ---------- 构建12个节点 ----------
        // 单元0的8个节点（单位立方体 [0,1]^3）
        mesh_data_.nodes = {
            {0, 0.0, 0.0, 0.0, 0},  // 节点0
            {1, 1.0, 0.0, 0.0, 0},  // 节点1
            {2, 1.0, 1.0, 0.0, 0},  // 节点2
            {3, 0.0, 1.0, 0.0, 0},  // 节点3
            {4, 0.0, 0.0, 1.0, 0},  // 节点4（共享面左下前）
            {5, 1.0, 0.0, 1.0, 0},  // 节点5（共享面右下前）
            {6, 1.0, 1.0, 1.0, 0},  // 节点6（共享面右上后）
            {7, 0.0, 1.0, 1.0, 0},  // 节点7（共享面左上后）
            // 单元1的新增4个节点（在z=2平面）
            {8, 0.0, 0.0, 2.0, 0},  // 节点8
            {9, 1.0, 0.0, 2.0, 0},  // 节点9
            {10, 1.0, 1.0, 2.0, 0}, // 节点10
            {11, 0.0, 1.0, 2.0, 0}, // 节点11
        };

        // ---------- 构建2个HEX8单元 ----------
        mesh_data_.elements = {
            {
                0,
                {0, 1, 2, 3, 4, 5, 6, 7},  // 单元0: 底部立方体
                ElemType::HEX8,
                DOFType::VECTOR_EDGE_ONLY,
                1,
                0,
            },
            {
                1,
                {4, 5, 6, 7, 8, 9, 10, 11},  // 单元1: 顶部立方体（共享z=1面）
                ElemType::HEX8,
                DOFType::VECTOR_EDGE_ONLY,
                1,
                0,
            },
        };

        edge_mapping_ = generateEdgeMapping(mesh_data_);

        FEEM_INFO("Hex8MultiElementFixture: 网格构建完成 - {}节点, {}单元",
                  mesh_data_.nodes.size(), mesh_data_.elements.size());
    }

    EMMeshData mesh_data_;
    std::vector<std::vector<int>> edge_mapping_;
};


// ---------- 2.1 多单元棱边去重测试 ----------

/**
 * @test 验证两个HEX8单元的棱边总数
 * @details 2个独立HEX8应有2×12=24条局部棱边，
 *          但共享4条棱边后，唯一棱边数=24-4=20
 */
TEST_F(Hex8MultiElementFixture, EdgeDeduplication) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    EXPECT_EQ(result.num_nodes, 12) << "12个节点应全部参与图";
    EXPECT_EQ(result.num_edges, 20)
        << "2个HEX8共享4条棱边后，唯一棱边数应为20（而非24）";

    FEEM_INFO("棱边去重验证: 2×12={}条局部棱边 → {}条唯一棱边（去重{}条）",
              24, result.num_edges, 24 - result.num_edges);
}

/**
 * @test 验证共享面上节点的度数增加
 * @details 共享面上的节点（4,5,6,7）因连接两个单元，
 *          其度数应大于非共享节点
 */
TEST_F(Hex8MultiElementFixture, SharedNodeDegreeIncrease) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    // 获取共享节点的索引
    std::set<int> shared_node_globals = {4, 5, 6, 7};
    std::vector<int> shared_indices;
    std::vector<int> non_shared_indices;

    for (int i = 0; i < result.num_nodes; ++i) {
        if (shared_node_globals.count(result.index_to_global_node[i])) {
            shared_indices.push_back(i);
        } else {
            non_shared_indices.push_back(i);
        }
    }

    ASSERT_EQ(static_cast<int>(shared_indices.size()), 4)
        << "应有4个共享节点";
    ASSERT_EQ(static_cast<int>(non_shared_indices.size()), 8)
        << "应有8个非共享节点";

    // 计算平均度数
    double avg_shared_degree = 0.0;
    for (int idx : shared_indices) {
        avg_shared_degree += result.adjacency_list[idx].size();
    }
    avg_shared_degree /= shared_indices.size();

    double avg_non_shared_degree = 0.0;
    for (int idx : non_shared_indices) {
        avg_non_shared_degree += result.adjacency_list[idx].size();
    }
    avg_non_shared_degree /= non_shared_indices.size();

    EXPECT_GT(avg_shared_degree, avg_non_shared_degree)
        << "共享节点平均度数(" << avg_shared_degree
        << ")应 > 非共享节点平均度数(" << avg_non_shared_degree << ")";

    FEEM_INFO("度数对比: 共享节点平均={:.2f}, 非共享节点平均={:.2f}",
              avg_shared_degree, avg_non_shared_degree);
}


// ====================================================================
//  测试固件3：边界条件和异常输入
// ====================================================================

/**
 * @class EdgeCaseFixture
 * @brief 异常输入测试固件
 * @details 测试各种边界条件下的行为：
 *          - 空网格
 *          - 纯标量单元网格
 *          - 输入尺寸不匹配
 */
class EdgeCaseFixture : public ::testing::Test {
protected:
    void SetUp() override {
        FEEM_LOG_INIT("", false);
    }
};


// ---------- 3.1 空网格测试 ----------

/**
 * @test 验证空网格返回无效结果
 * @details 当mesh_data不包含任何单元时，应返回空的GraphResult
 */
TEST_F(EdgeCaseFixture, EmptyMesh) {
    EMMeshData empty_mesh;
    std::vector<std::vector<int>> empty_mapping;

    MeshGraphExtractor extractor(empty_mesh, empty_mapping);
    GraphResult result = extractor.extract();

    EXPECT_FALSE(result.isValid()) << "空网格应返回无效结果";
    EXPECT_EQ(result.num_nodes, 0) << "节点数应为0";
    EXPECT_EQ(result.num_edges, 0) << "棱边数应为0";
    EXPECT_TRUE(result.adjacency_list.empty()) << "邻接表应为空";

    FEEM_INFO("空网格测试: 返回无效结果 ✓");
}


// ---------- 3.2 纯标量单元测试 ----------

/**
 * @test 验证纯标量单元网格返回空结果
 * @details 当所有单元都是SCALAR_ONLY类型时，应返回空图（无矢量棱边）
 */
TEST_F(EdgeCaseFixture, ScalarOnlyMesh) {
    EMMeshData scalar_mesh;

    scalar_mesh.nodes = {
        {0, 0.0, 0.0, 0.0, 0},
        {1, 1.0, 0.0, 0.0, 0},
        {2, 0.0, 1.0, 0.0, 0},
    };

    scalar_mesh.elements = {{
        0,
        {0, 1, 2},
        ElemType::TRI3,
        DOFType::SCALAR_ONLY,
        1,
        0,
    }};

    auto mapping = generateEdgeMapping(scalar_mesh);

    MeshGraphExtractor extractor(scalar_mesh, mapping);
    GraphResult result = extractor.extract();

    EXPECT_FALSE(result.isValid()) << "纯标量网格应返回无效结果";
    EXPECT_EQ(result.num_nodes, 0) << "标量单元的节点不应包含在图中";

    FEEM_INFO("纯标量单元测试: 返回空结果 ✓");
}


// ---------- 3.3 输入尺寸不匹配测试 ----------

/**
 * @test 验证elem_local_to_global_edge尺寸不匹配时抛出异常
 * @details 当映射表尺寸与单元数不一致时，应抛出std::invalid_argument
 */
TEST_F(EdgeCaseFixture, InputSizeMismatch) {
    EMMeshData test_mesh;

    test_mesh.nodes = {{0, 0, 0, 0, 0}, {1, 1, 0, 0, 0}};
    test_mesh.elements = {{
        0, {0, 1}, ElemType::LINE2, DOFType::VECTOR_EDGE_ONLY, 1, 0
    }};

    // 故意传入错误尺寸的映射表（空向量，而实际有1个单元）
    std::vector<std::vector<int>> wrong_size_mapping;

    MeshGraphExtractor extractor(test_mesh, wrong_size_mapping);
    EXPECT_THROW(extractor.extract(), std::invalid_argument)
        << "尺寸不匹配时应抛出std::invalid_argument";

    FEEM_INFO("输入尺寸不匹配测试: 正确抛出异常 ✓");
}


// ---------- 3.4 混合单元类型过滤测试 ----------

/**
 * @test 验证混合单元网格仅提取矢量单元的图结构
 * @details 网格同时包含SCALAR_ONLY和VECTOR_EDGE_ONLY单元时，
 *          应忽略标量单元，仅从矢量单元提取图
 */
TEST_F(EdgeCaseFixture, MixedElementTypeFiltering) {
    EMMeshData mixed_mesh;

    // 4个节点
    mixed_mesh.nodes = {
        {0, 0.0, 0.0, 0.0, 0},
        {1, 1.0, 0.0, 0.0, 0},
        {2, 0.0, 1.0, 0.0, 0},
        {3, 1.0, 1.0, 0.0, 0},
    };

    // 2个单元：1个标量 + 1个矢量
    mixed_mesh.elements = {
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

    auto mapping = generateEdgeMapping(mixed_mesh);

    MeshGraphExtractor extractor(mixed_mesh, mapping);
    GraphResult result = extractor.extract();

    EXPECT_TRUE(result.isValid()) << "应成功提取矢量单元的图";
    EXPECT_EQ(result.num_nodes, 4) << "4个节点应全部包含（来自矢量单元）";

    // QUAD4有4条棱边
    EXPECT_EQ(result.num_edges, 4) << "QUAD4矢量单元应有4条棱边";

    FEEM_INFO("混合单元过滤测试: 仅提取矢量单元 ✓ (节点={}, 棱边={})",
              result.num_nodes, result.num_edges);
}


// ====================================================================
//  综合输出测试
// ====================================================================

/**
 * @test 输出PRISM6单元的完整图结构信息供人工审查
 * @details 将邻接表、棱边映射、节点映射全部格式化输出
 */
TEST_F(Prism6GraphFixture, FullGraphDump) {
    MeshGraphExtractor extractor(mesh_data_, edge_mapping_);
    GraphResult result = extractor.extract();

    FEEM_INFO("");
    FEEM_INFO("============================================================");
    FEEM_INFO("  PRISM6 MIXED_AV 完整图结构输出");
    FEEM_INFO("============================================================");

    FEEM_INFO("--- 节点索引映射 ---");
    for (int i = 0; i < result.num_nodes; ++i) {
        FEEM_INFO("  索引{} ↔ 全局节点ID={}", i, result.index_to_global_node[i]);
    }

    FEEM_INFO("--- 邻接表 ---");
    for (int i = 0; i < result.num_nodes; ++i) {
        int global_id = result.index_to_global_node[i];
        std::string neighbors_str;
        for (int neighbor_idx : result.adjacency_list[i]) {
            int neighbor_global = result.index_to_global_node[neighbor_idx];
            neighbors_str += std::to_string(neighbor_global) + " ";
        }
        FEEM_INFO("  节点{}(全局{}): [{}]", i, global_id, neighbors_str);
    }

    FEEM_INFO("--- 棱边→全局ID映射 (共{}条) ---", result.num_edges);
    for (const auto& [edge_key, global_id] : result.edge_to_global_id) {
        int node1_global = result.index_to_global_node[edge_key.first];
        int node2_global = result.index_to_global_node[edge_key.second];
        FEEM_INFO("  棱边({},{}) [索引({},{})] → 全局棱边ID={}",
                  node1_global, node2_global,
                  edge_key.first, edge_key.second, global_id);
    }

    FEEM_INFO("--- 统计摘要 ---");
    int total_neighbor_count = 0;
    for (const auto& neighbors : result.adjacency_list) {
        total_neighbor_count += static_cast<int>(neighbors.size());
    }
    double avg_degree = static_cast<double>(total_neighbor_count) / result.num_nodes;

    FEEM_INFO("  总节点数: {}", result.num_nodes);
    FEEM_INFO("  总棱边数: {}", result.num_edges);
    FEEM_INFO("  总邻接关系数: {} (应为2×棱边数={})", total_neighbor_count, 2 * result.num_edges);
    FEEM_INFO("  平均节点度数: {:.2f}", avg_degree);
    FEEM_INFO("============================================================");
    FEEM_INFO("");

    // 最终断言：确保统计数据自洽
    EXPECT_EQ(total_neighbor_count, 2 * result.num_edges)
        << "总邻接关系数应为棱边数的2倍（无向图每条棱边贡献2个邻接关系）";
}


// ==================== 主函数 ====================

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
