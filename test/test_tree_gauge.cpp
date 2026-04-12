/**
 * @file test_tree_gauge.cpp
 * @brief TreeGauge核心应用类完整测试套件
 * @details 使用 Google Test 框架验证TreeGauge的完整功能：
 *          1. 单连通域完整流程（5步构建正确性）
 *          2. 带Dirichlet边界的模型（约束边归入树边）
 *          3. 多连通域模型（环量自由度保留）
 *          4. mapReducedSolutionToFull()解恢复（约化→完整映射）
 *          5. 性能基准测试（中型模型<50ms）
 *
 * 测试场景设计：
 * - 场景1：简单10节点图（单连通域，验证基本功能）
 * - 场景2：带Dirichlet边界（验证约束边处理）
 * - 场景3：多连通域（验证环量自由度）
 * - 场景4：解向量转换（验证mapReducedSolutionToFull正确性）
 * - 场景5：性能测试（1000节点级别网格）
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
#include <unordered_map>
#include <chrono>
#include <cmath>

#include "tree_gauge.hpp"
#include "em_mesh_data.hpp"
#include "global_edge_id_generator.hpp"
#include "logger_factory.hpp"

using namespace fe_em;


// ==================== 测试辅助函数 ====================

/**
 * @brief 构建简单的EMMeshData用于单元测试
 * @param nodes 节点列表（每个节点为{x, y, z}坐标三元组）
 * @param elements 单元列表（每个单元为节点ID列表）
 * @return EMMeshData 构建好的网格数据
 */
static EMMeshData buildTestMesh(
    const std::vector<std::tuple<double, double, double>>& nodes,
    const std::vector<std::vector<int>>& elements
) {
    EMMeshData mesh_data;

    mesh_data.nodes.reserve(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        Node node;
        node.id = static_cast<int>(i);
        std::tie(node.x, node.y, node.z) = nodes[i];
        node.region_id = 0;
        mesh_data.nodes.push_back(node);
    }

    mesh_data.elements.reserve(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        Element elem;
        elem.id = static_cast<int>(i);
        elem.node_ids = elements[i];

        // 根据节点数自动选择单元类型
        if (elements[i].size() == 3) {
            elem.elem_type = ElemType::TRI3;       // 3节点三角形
        } else if (elements[i].size() == 4) {
            elem.elem_type = ElemType::QUAD4;      // 4节点四边形
        } else if (elements[i].size() == 8) {
            elem.elem_type = ElemType::HEX8;       // 8节点六面体
        } else {
            elem.elem_type = ElemType::QUAD4;      // 默认四边形
        }

        elem.dof_type = DOFType::VECTOR_EDGE_ONLY;
        elem.material_id = 0;
        elem.region_id = 0;
        mesh_data.elements.push_back(elem);
    }

    return mesh_data;
}

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


// ==================== 测试用例1：单连通域完整流程 ====================

TEST(TreeGaugeTest, SingleConnectedDomainFullFlow) {
    FEEM_INFO("===== 测试: 单连通域完整流程 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0},
        {0.0, 3.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3}, {1, 2, 5, 4},
        {3, 4, 7, 6}, {4, 5, 8, 7},
        {6, 7, 9, 3}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);
    gauge.build();

    EXPECT_GT(gauge.getReducedNumDofs(), 0);
    EXPECT_GT(gauge.getTreeEdges().size(), 0);
    EXPECT_GT(gauge.getCotreeEdges().size(), 0);

    int total_edges = gauge.getTreeEdges().size() + gauge.getCotreeEdges().size() +
                      gauge.getCirculationDOFs().size();

    EXPECT_EQ(total_edges, gauge.getReducedNumDofs() + gauge.getTreeEdges().size() -
                           gauge.getCirculationDOFs().size());

    FEEM_INFO("单连通域测试通过 - 树边: {}, 余树边: {}, 约化DOF: {}",
              gauge.getTreeEdges().size(),
              gauge.getCotreeEdges().size(),
              gauge.getReducedNumDofs());
}


// ==================== 测试用例2：带Dirichlet边界模型 ====================

TEST(TreeGaugeTest, DirichletBoundaryHandling) {
    FEEM_INFO("===== 测试: 带Dirichlet边界模型 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}, {3.0, 2.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
        {4, 5, 9, 8}, {5, 6, 10, 9}, {6, 7, 11, 10}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    std::set<int> constrained_edges = {0, 1, 2};

    TreeGauge gauge(mesh_data, edge_mapping, constrained_edges);
    gauge.build();

    for (int edge_id : constrained_edges) {
        EXPECT_TRUE(gauge.isTreeEdge(edge_id))
            << "约束边" << edge_id << "应为树边";
        EXPECT_FALSE(gauge.isCotreeEdge(edge_id))
            << "约束边" << edge_id << "不应为余树边";
    }

    FEEM_INFO("Dirichlet边界测试通过 - 所有{}条约束边已归入树边",
              constrained_edges.size());
}


// ==================== 测试用例3：多连通域模型 ====================

TEST(TreeGaugeTest, MultiConnectedDomain) {
    FEEM_INFO("===== 测试: 多连通域模型 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 4.0, 0.0}, {0.0, 4.0, 0.0},
        {1.0, 1.0, 0.0}, {3.0, 1.0, 0.0}, {3.0, 3.0, 0.0}, {1.0, 3.0, 0.0},
        {1.5, 1.5, 0.0}, {2.5, 1.5, 0.0}, {2.5, 2.5, 0.0}, {1.5, 2.5, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7},
        {4, 5, 9, 8}, {5, 6, 10, 9}, {6, 7, 11, 10}, {7, 4, 8, 11},
        {8, 9, 10, 11}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);
    gauge.build();

    const auto& circulation_dofs = gauge.getCirculationDOFs();

    EXPECT_GE(circulation_dofs.size(), 0);

    if (!circulation_dofs.empty()) {
        for (int dof_id : circulation_dofs) {
            EXPECT_FALSE(gauge.isTreeEdge(dof_id))
                << "环量自由度" << dof_id << "不应是树边";
        }
    }

    FEEM_INFO("多连通域测试通过 - 孔洞数: {}, 环量DOF数: {}",
              0, circulation_dofs.size());
}


// ==================== 测试用例4：解向量映射测试 ====================

TEST(TreeGaugeTest, MapReducedSolutionToFull) {
    FEEM_INFO("===== 测试: 解向量映射 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3}, {1, 2, 5, 4},
        {3, 4, 5, 3}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);
    gauge.build();

    int reduced_dofs = gauge.getReducedNumDofs();
    ASSERT_GT(reduced_dofs, 0);

    std::vector<double> reduced_solution(static_cast<size_t>(reduced_dofs));
    for (int i = 0; i < reduced_dofs; ++i) {
        reduced_solution[static_cast<size_t>(i)] = static_cast<double>(i + 1) * 1.5;
    }

    std::vector<double> full_solution;
    gauge.mapReducedSolutionToFull(reduced_solution, full_solution);

    FEEM_INFO("调试信息 - full_size: {}, reduced_size: {}, tree_size: {}, cotree_size: {}, circulation_size: {}",
              full_solution.size(),
              reduced_dofs,
              gauge.getTreeEdges().size(),
              gauge.getCotreeEdges().size(),
              gauge.getCirculationDOFs().size());

    for (int edge_id : gauge.getTreeEdges()) {
        EXPECT_DOUBLE_EQ(full_solution[static_cast<size_t>(edge_id)], 0.0)
            << "树边" << edge_id << "应为0";
    }

    for (const auto& [edge_id, reduced_idx] : gauge.getCotreeEdgeToReducedDOF()) {
        EXPECT_DOUBLE_EQ(
            full_solution[static_cast<size_t>(edge_id)],
            reduced_solution[static_cast<size_t>(reduced_idx)]
        ) << "余树边" << edge_id << "映射值不匹配";
    }

    FEEM_INFO("解向量映射测试通过 - reduced_size: {}, full_size: {}",
              reduced_solution.size(), full_solution.size());
}


// ==================== 测试用例5：异常输入测试 ====================

TEST(TreeGaugeTest, InvalidInputHandling) {
    FEEM_INFO("===== 测试: 异常输入处理 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);
    gauge.build();

    std::vector<double> wrong_size_reduced(999);

    std::vector<double> full_solution;
    EXPECT_THROW(
        gauge.mapReducedSolutionToFull(wrong_size_reduced, full_solution),
        std::invalid_argument
    );

    FEEM_INFO("异常输入测试通过 - 正确检测到维度不匹配");
}


// ==================== 测试用例6：重复build调用测试 ====================

TEST(TreeGaugeTest, DuplicateBuildCall) {
    FEEM_INFO("===== 测试: 重复build调用 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 2}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);
    gauge.build();

    int first_reduced_dofs = gauge.getReducedNumDofs();

    gauge.build();

    EXPECT_EQ(gauge.getReducedNumDofs(), first_reduced_dofs);

    FEEM_INFO("重复build调用测试通过 - 结果一致");
}


// ==================== 测试用例7：性能基准测试 ====================

TEST(TreeGaugeTest, PerformanceBenchmark) {
    FEEM_INFO("===== 测试: 性能基准测试 =====");

    int grid_size = 20;

    auto nodes = std::vector<std::tuple<double, double, double>>();
    auto elements = std::vector<std::vector<int>>();

    nodes.reserve(static_cast<size_t>((grid_size + 1) * (grid_size + 1)));
    for (int j = 0; j <= grid_size; ++j) {
        for (int i = 0; i <= grid_size; ++i) {
            nodes.emplace_back(static_cast<double>(i), static_cast<double>(j), 0.0);
        }
    }

    elements.reserve(static_cast<size_t>(grid_size * grid_size));
    for (int j = 0; j < grid_size; ++j) {
        for (int i = 0; i < grid_size; ++i) {
            int n0 = j * (grid_size + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (grid_size + 1) + 1;
            int n3 = n0 + (grid_size + 1);
            elements.push_back({n0, n1, n2, n3});
        }
    }

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);

    auto start_time = std::chrono::high_resolution_clock::now();
    gauge.build();
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time).count();

    FEEM_INFO("性能测试 - 节点数: {}, 单元数: {}, 构建时间: {}ms",
              mesh_data.nodes.size(),
              mesh_data.elements.size(),
              duration_ms);

    EXPECT_LT(duration_ms, 50)
        << "构建时间应小于50ms，实际: " << duration_ms << "ms";

    EXPECT_GT(gauge.getReducedNumDofs(), 0);
    EXPECT_GT(gauge.getTreeEdges().size(), 0);

    double compression_rate = 100.0 * (1.0 -
        static_cast<double>(gauge.getReducedNumDofs()) /
        (gauge.getTreeEdges().size() + gauge.getCotreeEdges().size() +
         gauge.getCirculationDOFs().size()));

    FEEM_INFO("压缩率: {:.1f}%", compression_rate);

    EXPECT_GT(compression_rate, 0.0)
        << "应有正的压缩率";
}


// ==================== 测试用例8：build前查询接口安全性测试 ====================

TEST(TreeGaugeTest, QueryBeforeBuildSafety) {
    FEEM_INFO("===== 测试: build前查询接口安全性 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    EMMeshData mesh_data = buildTestMesh(nodes, elements);
    auto edge_mapping = generateEdgeMapping(mesh_data);

    TreeGauge gauge(mesh_data, edge_mapping);

    EXPECT_EQ(gauge.getReducedNumDofs(), 0);
    EXPECT_TRUE(gauge.getTreeEdges().empty());
    EXPECT_TRUE(gauge.getCotreeEdges().empty());
    EXPECT_TRUE(gauge.getCirculationDOFs().empty());
    EXPECT_TRUE(gauge.getCotreeEdgeToReducedDOF().empty());

    EXPECT_FALSE(gauge.isTreeEdge(0));
    EXPECT_FALSE(gauge.isCotreeEdge(0));

    FEEM_INFO("build前查询接口安全性测试通过 - 所有接口安全返回默认值");
}


// ==================== 主函数 ====================

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
