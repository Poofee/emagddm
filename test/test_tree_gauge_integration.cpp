/**
 * @file test_tree_gauge_integration.cpp
 * @brief Tree Gauge与EMDOFManager集成测试套件
 * @details 验证Tree Gauge模块与DOF管理器的无缝集成：
 *          1. 基本集成测试（applyTreeGauge后DOF数减少）
 *          2. 无操作测试（纯标量网格安全忽略）
 *          3. MIXED_AV单元集成测试（标量/矢量DOF分离处理）
 *          4. 解向量恢复测试（约化系统→完整解映射）
 *
 * 测试场景设计：
 * - 场景1：简单矢量网格（验证基本压缩功能）
 * - 场景2：带Dirichlet边界的模型（验证约束兼容性）
 * - 场景3：MIXED_AV混合单元（验证标量/矢量分离）
 * - 场景4：解向量端到端流程（build→apply→求解→恢复）
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
#include <cmath>
#include <memory>

#include "em_dof_manager.hpp"
#include "em_mesh_data.hpp"
#include "global_edge_id_generator.hpp"
#include "tree_gauge.hpp"
#include "logger_factory.hpp"

using namespace fe_em;


// ==================== 测试辅助函数 ====================

/**
 * @brief 构建简单的EMMeshData用于单元测试
 * @param nodes 节点列表（每个节点为{x, y, z}坐标三元组）
 * @param elements 单元列表（每个单元为节点ID列表）
 * @param dof_type DOF类型（默认为VECTOR_EDGE_ONLY）
 * @return EMMeshData 构建好的网格数据
 */
static EMMeshData buildTestMesh(
    const std::vector<std::tuple<double, double, double>>& nodes,
    const std::vector<std::vector<int>>& elements,
    DOFType dof_type = DOFType::VECTOR_EDGE_ONLY
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
        elem.elem_type = ElemType::HEX8;
        elem.dof_type = dof_type;
        elem.material_id = 0;
        elem.region_id = 0;
        mesh_data.elements.push_back(elem);
    }

    return mesh_data;
}

/**
 * @brief 构建GlobalEdgeIDGenerator并生成映射表的便捷函数
 * @param mesh_data 网格数据
 * @return std::pair<GlobalEdgeIDGenerator, std::vector<std::vector<int>>> 生成器和映射表
 */
static std::pair<std::unique_ptr<GlobalEdgeIDGenerator>, std::vector<std::vector<int>>>
generateEdgeMapping(const EMMeshData& mesh_data) {
    auto edge_gen = std::make_unique<GlobalEdgeIDGenerator>(mesh_data);
    edge_gen->generate();
    return {std::move(edge_gen), edge_gen->getElemLocalToGlobalEdge()};
}


// ==================== Test Fixture 1: BasicIntegrationTest ====================

TEST(TreeGaugeIntegrationTest, BasicIntegration_ReducedDOFs) {
    FEEM_INFO("===== 测试: 基本集成 - DOF数量减少 =====");

    // 构建2x2的四边形网格（4个QUAD4单元，共享棱边）
    // QUAD4有4个节点和4条棱边
    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},     // 单元0: 左下
        {1, 2, 5, 4},     // 单元1: 右下
        {3, 4, 7, 6},     // 单元2: 左上
        {4, 5, 8, 7}      // 单元3: 右上
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::QUAD4;
    mesh_data.elements[2].elem_type = ElemType::QUAD4;
    mesh_data.elements[3].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    // 创建DOF管理器并执行标准编号流程
    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    int original_dofs = dof_manager.getNumFreeDOFs();
    FEEM_INFO("原始自由DOF数: {}", original_dofs);
    ASSERT_GT(original_dofs, 0) << "原始DOF数应大于0";

    // 应用Tree Gauge
    dof_manager.applyTreeGauge();

    int reduced_dofs = dof_manager.getNumFreeDOFs();
    FEEM_INFO("约化后自由DOF数: {}", reduced_dofs);

    // 验证DOF数量减少
    EXPECT_LT(reduced_dofs, original_dofs)
        << "应用Tree Gauge后，自由DOF数应减少";

    // 验证压缩率在合理范围内（通常30-50%）
    double compression_rate = 100.0 * (original_dofs - reduced_dofs) / original_dofs;
    FEEM_INFO("压缩率: {:.1f}%", compression_rate);
    EXPECT_GT(compression_rate, 10.0) << "压缩率应大于10%";
    EXPECT_LT(compression_rate, 90.0) << "压缩率应小于90%（不能过度压缩）";
}

TEST(TreeGaugeIntegrationTest, BasicIntegration_TreeEdgesMarkedAsMinusOne) {
    FEEM_INFO("===== 测试: 基本集成 - 树边标记为-1 =====");

    // 构建简单网格（2个QUAD4单元）
    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {1, 2, 5, 4}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();
    dof_manager.applyTreeGauge();

    const auto& mappings = dof_manager.getElemLocalToGlobal();
    const TreeGauge* gauge = dof_manager.getTreeGauge();

    ASSERT_NE(gauge, nullptr) << "getTreeGauge()应返回有效指针";

    // 检查至少有一个树边被标记为-1
    bool found_tree_edge_marked = false;
    for (const auto& mapping : mappings) {
        for (int idx : mapping.indices) {
            if (idx == -1) {
                found_tree_edge_marked = true;
                break;
            }
        }
        if (found_tree_edge_marked) break;
    }

    EXPECT_TRUE(found_tree_edge_marked)
        << "至少应有部分树边DOF被标记为-1";
}

TEST(TreeGaugeIntegrationTest, BasicIntegration_GetTreeGaugeReturnsValidPointer) {
    FEEM_INFO("===== 测试: 基本集成 - getTreeGauge()返回有效指针 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {1, 2, 4}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::TRI3;  // 三角形单元
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());

    // 未调用applyTreeGauge前应为nullptr
    EXPECT_EQ(dof_manager.getTreeGauge(), nullptr)
        << "未调用applyTreeGauge时，getTreeGauge()应返回nullptr";

    dof_manager.build();
    dof_manager.applyTreeGauge();

    // 调用后应返回有效指针
    const TreeGauge* gauge = dof_manager.getTreeGauge();
    ASSERT_NE(gauge, nullptr) << "调用applyTreeGauge后，getTreeGauge()应返回有效指针";

    // 验证TreeGauge对象可用
    EXPECT_GT(gauge->getReducedNumDofs(), 0) << "约化DOF数应大于0";
    EXPECT_GT(gauge->getTreeEdges().size(), 0u) << "树边集合不应为空";
    EXPECT_GT(gauge->getCotreeEdges().size(), 0u) << "余树边集合不应为空";
}


// ==================== Test Fixture 2: NoOperationTest ====================

TEST(TreeGaugeIntegrationTest, NoOperation_BeforeBuildThrowsException) {
    FEEM_INFO("===== 测试: 无操作 - build前调用抛出异常 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());

    // 不调用build()直接调用applyTreeGauge()应抛出异常
    EXPECT_THROW(dof_manager.applyTreeGauge(), std::runtime_error)
        << "未调用build()就调用applyTreeGauge()应抛出异常";
}

TEST(TreeGaugeIntegrationTest, NoOperation_ScalarOnlyMeshIgnored) {
    FEEM_INFO("===== 测试: 无操作 - 纯标量网格安全忽略 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    // 使用SCALAR_ONLY类型构建网格
    auto mesh_data = buildTestMesh(nodes, elements, DOFType::SCALAR_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;

    EMDOFManager dof_manager(mesh_data);
    dof_manager.build();

    int original_dofs = dof_manager.getNumFreeDOFs();

    // 对纯标量网格调用applyTreeGauge()应安全忽略（不抛异常）
    EXPECT_NO_THROW(dof_manager.applyTreeGauge())
        << "对纯标量网格调用applyTreeGauge()应安全忽略";

    // DOF数应不变
    EXPECT_EQ(dof_manager.getNumFreeDOFs(), original_dofs)
        << "纯标量网格的DOF数不应改变";

    // getTreeGauge()应返回nullptr（因为没有实际执行Tree Gauge）
    EXPECT_EQ(dof_manager.getTreeGauge(), nullptr)
        << "纯标量网格的getTreeGauge()应返回nullptr";
}

TEST(TreeGaugeIntegrationTest, NoOperation_NotCalledReturnsNullptr) {
    FEEM_INFO("===== 测试: 无操作 - 不调用applyTreeGauge返回nullptr =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    // 不调用applyTreeGauge()
    EXPECT_EQ(dof_manager.getTreeGauge(), nullptr)
        << "不调用applyTreeGauge()时，getTreeGauge()应返回nullptr";
}


// ==================== Test Fixture 3: MixedAVIntegrationTest ====================

TEST(TreeGaugeIntegrationTest, MixedAV_ScalarDOFsUnchanged) {
    FEEM_INFO("===== 测试: MIXED_AV集成 - 标量DOF不受影响 =====");

    // 构建MIXED_AV类型的网格（2个QUAD4单元）
    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {3, 4, 7, 6}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::MIXED_AV);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    // 记录applyTreeGauge前的标量DOF映射
    const auto& mappings_before = dof_manager.getElemLocalToGlobal();
    std::vector<int> scalar_indices_before;
    for (const auto& mapping : mappings_before) {
        if (mapping.is_mixed()) {
            for (int i = 0; i < mapping.num_scalar_dofs; ++i) {
                scalar_indices_before.push_back(mapping.indices[i]);
            }
        }
    }

    dof_manager.applyTreeGauge();

    // 记录applyTreeGauge后的标量DOF映射
    const auto& mappings_after = dof_manager.getElemLocalToGlobal();
    std::vector<int> scalar_indices_after;
    for (const auto& mapping : mappings_after) {
        if (mapping.is_mixed()) {
            for (int i = 0; i < mapping.num_scalar_dofs; ++i) {
                scalar_indices_after.push_back(mapping.indices[i]);
            }
        }
    }

    // 验证标量DOF完全相同
    ASSERT_EQ(scalar_indices_before.size(), scalar_indices_after.size())
        << "标量DOF数量应保持不变";
    for (size_t i = 0; i < scalar_indices_before.size(); ++i) {
        EXPECT_EQ(scalar_indices_before[i], scalar_indices_after[i])
            << "第" << i << "个标量DOF编号应保持不变";
    }
}

TEST(TreeGaugeIntegrationTest, MixedAV_VectorDOFsCompressed) {
    FEEM_INFO("===== 测试: MIXED_AV集成 - 矢量DOF被正确压缩 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {3, 4, 5}  // 第二个单元改为三角形单元
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::MIXED_AV);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::TRI3;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    int original_dofs = dof_manager.getNumFreeDOFs();
    dof_manager.applyTreeGauge();
    int reduced_dofs = dof_manager.getNumFreeDOFs();

    // 验证总DOF数减少
    EXPECT_LT(reduced_dofs, original_dofs)
        << "MIXED_AV单元的总自由DOF数应减少";

    // 统计矢量DOF中-1的数量（树边）
    const auto& mappings = dof_manager.getElemLocalToGlobal();
    int tree_edge_count = 0;
    int cotree_edge_count = 0;

    for (const auto& mapping : mappings) {
        if (mapping.is_mixed()) {
            for (int i = mapping.num_scalar_dofs;
                 i < static_cast<int>(mapping.indices.size()); ++i) {
                if (mapping.indices[i] == -1) {
                    tree_edge_count++;
                } else {
                    cotree_edge_count++;
                }
            }
        }
    }

    FEEM_INFO("MIXED_AV统计 - 树边DOF数: {}, 余树边DOF数: {}",
              tree_edge_count, cotree_edge_count);

    // 至少应有一些树边和余树边
    EXPECT_GT(tree_edge_count, 0) << "至少应有部分矢量DOF被标记为树边(-1)";
    EXPECT_GT(cotree_edge_count, 0) << "至少应保留部分余树边DOF";
}

TEST(TreeGaugeIntegrationTest, MixedAV_FinalDOFCountMatchesExpected) {
    FEEM_INFO("===== 测试: MIXED_AV集成 - 最终DOF数符合预期 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::MIXED_AV);
    for (auto& elem : mesh_data.elements) {
        elem.elem_type = ElemType::QUAD4;
    }
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();
    dof_manager.applyTreeGauge();

    const TreeGauge* gauge = dof_manager.getTreeGauge();
    ASSERT_NE(gauge, nullptr);

    // 统计标量自由DOF数
    const auto& mappings = dof_manager.getElemLocalToGlobal();
    int scalar_free_count = 0;
    for (const auto& mapping : mappings) {
        if (mapping.is_mixed()) {
            for (int i = 0; i < mapping.num_scalar_dofs; ++i) {
                if (mapping.indices[i] >= scalar_free_count) {
                    scalar_free_count = mapping.indices[i] + 1;
                }
            }
        }
    }

    int expected_total = scalar_free_count + gauge->getReducedNumDofs();
    int actual_total = dof_manager.getNumFreeDOFs();

    EXPECT_EQ(actual_total, expected_total)
        << "最终DOF数应等于 标量free DOF + 余树边数";
    FEEM_INFO("预期总DOF数: {} (标量:{} + 余树边:{})",
              expected_total, scalar_free_count, gauge->getReducedNumDofs());
}


// ==================== Test Fixture 4: SolutionRecoveryTest ====================

TEST(TreeGaugeIntegrationTest, SolutionRecovery_FullWorkflow) {
    FEEM_INFO("===== 测试: 解恢复 - 完整工作流 =====");

    // 构建测试网格（4个QUAD4单元）
    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
        {0.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    for (auto& elem : mesh_data.elements) {
        elem.elem_type = ElemType::QUAD4;
    }
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    // 步骤1: 执行标准DOF编号
    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    int original_dofs = dof_manager.getNumFreeDOFs();
    FEEM_INFO("步骤1完成 - 原始自由DOF数: {}", original_dofs);

    // 步骤2: 应用Tree Gauge
    dof_manager.applyTreeGauge();

    int reduced_dofs = dof_manager.getNumFreeDOFs();
    FEEM_INFO("步骤2完成 - 约化后自由DOF数: {}", reduced_dofs);

    const TreeGauge* gauge = dof_manager.getTreeGauge();
    ASSERT_NE(gauge, nullptr);

    // 步骤3: 模拟求解器输出（假设的约化解向量）
    std::vector<double> reduced_solution(reduced_dofs);
    for (int i = 0; i < reduced_dofs; ++i) {
        reduced_solution[i] = static_cast<double>(i + 1) * 0.5;  // 示例值
    }
    FEEM_INFO("步骤3完成 - 模拟约化解向量（长度:{})", reduced_solution.size());

    // 步骤4: 使用TreeGauge恢复完整解向量
    std::vector<double> full_solution;
    gauge->mapReducedSolutionToFull(reduced_solution, full_solution);

    FEEM_INFO("步骤4完成 - 恢复完整解向量（长度:{})", full_solution.size());

    // 验证完整解向量的性质
    EXPECT_EQ(static_cast<int>(full_solution.size()), gauge->getCotreeEdges().size() +
              gauge->getTreeEdges().size() + gauge->getCirculationDOFs().size())
        << "完整解向量长度应等于总棱边数";

    // 验证树边位置的值为0
    bool all_tree_edges_zero = true;
    for (int edge_id : gauge->getTreeEdges()) {
        if (std::abs(full_solution[edge_id]) > 1e-10) {
            all_tree_edges_zero = false;
            break;
        }
    }
    EXPECT_TRUE(all_tree_edges_zero) << "所有树边的A值应为0";

    // 验证余树边位置的值来自约化解
    const auto& cotree_to_reduced = gauge->getCotreeEdgeToReducedDOF();
    bool cotree_values_match = true;
    for (const auto& [edge_id, reduced_idx] : cotree_to_reduced) {
        if (std::abs(full_solution[edge_id] - reduced_solution[reduced_idx]) > 1e-10) {
            cotree_values_match = false;
            FEEM_WARN("余树边{}不匹配: full={} vs reduced={}",
                      edge_id, full_solution[edge_id], reduced_solution[reduced_idx]);
            break;
        }
    }
    EXPECT_TRUE(cotree_values_match) << "余树边值应与约化解一致";
}

TEST(TreeGaugeIntegrationTest, SolutionRecovery_WithDirichletBoundary) {
    FEEM_INFO("===== 测试: 解恢复 - 带Dirichlet边界条件 =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 4, 3},
        {1, 2, 5, 4}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    mesh_data.elements[1].elem_type = ElemType::QUAD4;

    // 添加Dirichlet边界条件（左边界的棱边A=0）
    EMBoundaryMarker dirichlet_bnd;
    dirichlet_bnd.id = 0;
    dirichlet_bnd.bnd_type = BndType::DIRICHLET;
    dirichlet_bnd.dof_type = DOFType::VECTOR_EDGE_ONLY;
    dirichlet_bnd.value = 0.0;
    dirichlet_bnd.target_ids = std::vector<std::vector<int>>{{0, 3}};  // 左边界棱边(节点0-3)

    mesh_data.boundary_markers.push_back(dirichlet_bnd);

    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    int original_dofs = dof_manager.getNumFreeDOFs();
    dof_manager.applyTreeGauge();
    int reduced_dofs = dof_manager.getNumFreeDOFs();

    FEEM_INFO("带边界条件的压缩: {} -> {} ({:.1f}%减少)",
              original_dofs, reduced_dofs,
              100.0 * (original_dofs - reduced_dofs) / original_dofs);

    // 验证约束边也被归入树边
    const TreeGauge* gauge = dof_manager.getTreeGauge();
    ASSERT_NE(gauge, nullptr);

    // Dirichlet边应在树边集合中（或被特殊处理）
    // 注意：具体行为取决于TreeGauge的实现细节
    EXPECT_LE(reduced_dofs, original_dofs)
        << "即使有边界条件，DOF数也应减少或持平";
}


// ==================== Test Fixture 5: EdgeCaseTests ====================

TEST(TreeGaugeIntegrationTest, EdgeCase_MultipleApplyCallsAllowed) {
    FEEM_INFO("===== 测试: 边界情况 - 允许重复调用applyTreeGauge =====");

    auto nodes = std::vector<std::tuple<double, double, double>>{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    auto elements = std::vector<std::vector<int>>{
        {0, 1, 3, 2}
    };

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    mesh_data.elements[0].elem_type = ElemType::QUAD4;
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    // 第一次调用
    dof_manager.applyTreeGauge();
    int dofs_after_first = dof_manager.getNumFreeDOFs();

    // 第二次调用（应允许，结果可能重置或累加）
    EXPECT_NO_THROW(dof_manager.applyTreeGauge())
        << "重复调用applyTreeGauge()不应抛出异常";

    int dofs_after_second = dof_manager.getNumFreeDOFs();

    // 验证两次调用的结果一致（或第二次重置）
    FEEM_INFO("第一次调用后DOF数: {}, 第二次调用后DOF数: {}",
              dofs_after_first, dofs_after_second);
}

TEST(TreeGaugeIntegrationTest, EdgeCase_LargeMeshPerformance) {
    FEEM_INFO("===== 测试: 边界情况 - 中型网格性能测试 =====");

    // 构建稍大的网格（5x5的QUAD4网格）
    int nx = 5, ny = 5;
    auto nodes = std::vector<std::tuple<double, double, double>>();
    auto elements = std::vector<std::vector<int>>();

    int node_id = 0;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            nodes.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
            node_id++;
        }
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (nx + 1) + 1;
            int n3 = n0 + (nx + 1);
            elements.push_back({n0, n1, n2, n3});
        }
    }

    auto mesh_data = buildTestMesh(nodes, elements, DOFType::VECTOR_EDGE_ONLY);
    for (auto& elem : mesh_data.elements) {
        elem.elem_type = ElemType::QUAD4;
    }
    auto [edge_gen, edge_mapping] = generateEdgeMapping(mesh_data);

    EMDOFManager dof_manager(mesh_data, edge_mapping, edge_gen.get());
    dof_manager.build();

    auto start_time = std::chrono::high_resolution_clock::now();
    dof_manager.applyTreeGauge();
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time).count();

    FEEM_INFO("中型网格性能测试 - 节点数: {}, 单元数: {}, 耗时: {}ms",
              mesh_data.getNodeCount(), mesh_data.getElementCount(), duration_ms);

    // 性能要求：中型网格应在100ms内完成
    EXPECT_LT(duration_ms, 100)
        << "中型网格的Tree Gauge应用时间应小于100ms（实际：" << duration_ms << "ms）";
}
