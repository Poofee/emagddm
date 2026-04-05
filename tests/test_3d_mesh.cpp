/**
 * @file test_3d_mesh.cpp
 * @brief 三维涡流场A-V混合格式测试用例
 * @details 测试TET4单元的网格数据结构、A-V混合DOF类型、导体材料属性、
 *          理想导体边界条件、ElementGeometry几何定义、ProjectManager集成。
 * @author Poofee
 * @date 2026-04-05
 */

#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
#include <cmath>

#include "fe_em/em_mesh_data.hpp"
#include "fe_em/element_geometry.hpp"
#include "fe_em/mesh_query.hpp"
#include "tool/project_manager.hpp"
#include "tool/project_data.hpp"

using namespace fe_em;
// 注意：tool命名空间不全局using，避免歧义
// 使用 tool:: 前缀访问ProjectManager, Material等

/**
 * @brief 测试1: 创建TET4单元模型并验证A-V混合DOF类型
 */
void test_tet4_av_mixed_dof() {
    std::cout << "\n========== 测试1: TET4 A-V混合DOF类型 ==========\n";
    
    // 创建EMMeshData
    auto mesh_data = std::make_unique<EMMeshData>();
    
    // 添加4个节点（单位正四面体的一部分）
    mesh_data->nodes.push_back({1, 0.0, 0.0, 0.0, 1});     // 节点1: 原点
    mesh_data->nodes.push_back({2, 1.0, 0.0, 0.0, 1});     // 节点2: X轴上
    mesh_data->nodes.push_back({3, 0.5, 0.866, 0.0, 1});   // 节点3: XY平面
    mesh_data->nodes.push_back({4, 0.5, 0.289, 0.816, 1}); // 节点4: Z方向顶点
    
    // 添加1个TET4单元（A-V混合格式用于涡流场导体区）
    Element tet_elem;
    tet_elem.id = 1;
    tet_elem.node_ids = {1, 2, 3, 4};  // 使用全局节点ID
    tet_elem.elem_type = ElemType::TET4;
    tet_elem.dof_type = DOFType::MIXED_AV;  // ⭐ 关键：A-V混合格式
    tet_elem.material_id = 1;               // 导体材料ID
    tet_elem.region_id = 1;                 // 导体区域
    
    mesh_data->elements.push_back(tet_elem);
    
    // 验证基本数据
    assert(mesh_data->getNodeCount() == 4);
    assert(mesh_data->getElementCount() == 1);
    
    // 验证DOF类型
    const auto& elem = mesh_data->elements[0];
    assert(elem.dof_type == DOFType::MIXED_AV);
    assert(elem.elem_type == ElemType::TET4);
    assert(elem.node_ids.size() == 4);
    
    std::cout << "✓ TET4 A-V混合DOF类型验证通过\n";
    std::cout << "  - 单元类型: TET4 (枚举值=" << static_cast<int>(elem.elem_type) << ")\n";
    std::cout << "  - DOF类型: MIXED_AV (A-V混合格式)\n";
    std::cout << "  - 节点数: " << elem.node_ids.size() << "\n";
}

/**
 * @brief 测试2: 验证ElementGeometry的TET4几何定义
 */
void test_tet4_element_geometry() {
    std::cout << "\n========== 测试2: TET4几何定义 ==========\n";
    
    // 检查支持性
    assert(ElementGeometry::is_supported(ElemType::TET4));
    
    // 获取拓扑信息
    int num_nodes = ElementGeometry::get_num_nodes(ElemType::TET4);
    int num_edges = ElementGeometry::get_num_edges(ElemType::TET4);
    int num_faces = ElementGeometry::get_num_faces(ElemType::TET4);
    
    assert(num_nodes == 4);    // TET4有4个节点
    assert(num_edges == 6);    // TET4有6条棱边 C(4,2)=6
    assert(num_faces == 4);    // TET4有4个三角形面
    
    // 获取并验证棱边定义
    auto edges = ElementGeometry::get_local_edges(ElemType::TET4);
    assert(edges.size() == 6);
    
    // 打印所有棱边
    std::cout << "  TET4棱边定义 (" << edges.size() << "条):\n";
    for (size_t i = 0; i < edges.size(); ++i) {
        auto [n1, n2] = edges[i];
        std::cout << "    Edge[" << i << "]: (" << n1 << ", " << n2 << ")\n";
        
        // 验证升序排列
        assert(n1 < n2);
    }
    
    // 验证标准6条棱边
    assert(std::get<0>(edges[0]) == 0 && std::get<1>(edges[0]) == 1);
    assert(std::get<0>(edges[1]) == 1 && std::get<1>(edges[1]) == 2);
    assert(std::get<0>(edges[2]) == 0 && std::get<1>(edges[2]) == 2);
    assert(std::get<0>(edges[3]) == 0 && std::get<1>(edges[3]) == 3);
    assert(std::get<0>(edges[4]) == 1 && std::get<1>(edges[4]) == 3);
    assert(std::get<0>(edges[5]) == 2 && std::get<1>(edges[5]) == 3);
    
    // 获取并验证面定义
    auto faces = ElementGeometry::get_local_faces(ElemType::TET4);
    assert(faces.size() == 4);
    
    std::cout << "  TET4面定义 (" << faces.size() << "个):\n";
    for (size_t i = 0; i < faces.size(); ++i) {
        std::cout << "    Face[" << i << "]: [";
        for (size_t j = 0; j < faces[i].size(); ++j) {
            std::cout << faces[i][j];
            if (j < faces[i].size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // 验证每个面都是三角形
    for (const auto& face : faces) {
        assert(face.size() == 3);  // TET4的所有面都是三角形
    }
    
    std::cout << "✓ TET4几何定义验证通过\n";
}

/**
 * @brief 测试3: 导体材料和理想导体边界条件
 */
void test_conductor_material_and_boundary() {
    std::cout << "\n========== 测试3: 导体材料与理想导体边界 ==========\n";
    
    auto& pm = tool::ProjectManager::getInstance();
    pm.createProject("Test3DProject");
    
    // 创建铜导体材料（高电导率）
    auto conductor_mat = std::make_shared<tool::Material>("Copper");
    conductor_mat->setRelativePermeability(1.0);         // μ_r = 1.0 (非磁性)
    conductor_mat->setPermittivity(1.0);                // ε_r ≈ 1
    conductor_mat->setConductivity(1e6);                // σ = 1×10^6 S/m (良导体)
    
    pm.addMaterial(conductor_mat);
    
    // 创建TET4单元网格
    auto em_mesh = std::make_unique<EMMeshData>();
    em_mesh->nodes = {
        {1, 0.0, 0.0, 0.0, 1},
        {2, 1.0, 0.0, 0.0, 1},
        {3, 0.5, 0.866, 0.0, 1},
        {4, 0.5, 0.289, 0.816, 1}
    };
    
    Element tet_elem;
    tet_elem.id = 1;
    tet_elem.node_ids = {1, 2, 3, 4};  // 使用全局节点ID
    tet_elem.elem_type = ElemType::TET4;
    tet_elem.dof_type = DOFType::MIXED_AV;  // 涡流场导体区必须用A-V格式
    tet_elem.material_id = 1;               // 引用Copper材料
    tet_elem.region_id = 1;
    
    em_mesh->elements.push_back(tet_elem);
    
    // 添加理想导体边界条件（Dirichlet棱边边界）- 组合设计
    // 物理意义: A×n = 0 （理想导体表面切向磁场为零）
    EMBoundaryMarker pec_boundary;
    pec_boundary.id = 1;
    pec_boundary.bnd_type = BndType::PERFECT_E;          // 边界类型：理想电壁（复用tool::BndType）
    pec_boundary.dof_type = DOFType::VECTOR_EDGE_ONLY;  // DOF类型：矢量棱边
    
    // 标记棱边(1,2)为理想导体边界（即节点ID=1和ID=2之间的边）
    pec_boundary.target_ids = std::vector<std::vector<int>>{{1, 2}};
    pec_boundary.value = 0.0;  // A_tangential = 0
    pec_boundary.name = "Ideal_Conductor_Edge";
    
    em_mesh->boundary_markers.push_back(pec_boundary);
    
    // 设置到ProjectManager
    pm.setEMMeshData(std::move(em_mesh));
    
    // 验证材料参数访问
    auto* mesh_data = pm.getEMMeshData();
    assert(mesh_data != nullptr);
    
    // 通过material_id索引获取材料参数（零拷贝！）
    auto copper_opt = pm.getMaterial("Copper");
    assert(copper_opt.has_value());
    
    auto* copper = copper_opt.value().get();
    double mu_r = copper->getRelativePermeability();
    double sigma = copper->getConductivity();
    
    assert(std::abs(mu_r - 1.0) < 1e-10);          // 铜是非磁性材料
    assert(std::abs(sigma - 1e6) < 1.0);           // 高电导率导体
    
    // 验证边界条件
    assert(mesh_data->getBoundaryMarkerCount() == 1);
    const auto& bnd = mesh_data->boundary_markers[0];
    assert(bnd.bnd_type == BndType::PERFECT_E);           // 验证BndType
    assert(bnd.dof_type == DOFType::VECTOR_EDGE_ONLY);     // 验证DOFType
    assert(bnd.name == "Ideal_Conductor_Edge");
    assert(std::abs(bnd.value - 0.0) < 1e-15);
    
    // 验证目标棱边
    if (std::holds_alternative<std::vector<std::vector<int>>>(bnd.target_ids)) {
        auto edge_list = std::get<std::vector<std::vector<int>>>(bnd.target_ids);
        assert(edge_list.size() == 1);
        assert(edge_list[0].size() == 2);
        assert(edge_list[0][0] == 1 && edge_list[0][1] == 2);  // 全局节点ID (1,2)
    }
    
    std::cout << "✓ 导体材料与理想导体边界验证通过\n";
    std::cout << "  - 材料: Copper (μ_r=" << mu_r << ", σ=" << sigma << " S/m)\n";
    std::cout << "  - DOF类型: MIXED_AV (涡流场导体区)\n";
    std::cout << "  - 边界组合: PERFECT_E + VECTOR_EDGE_ONLY (理想导体 A×n=0)\n";
    std::cout << "  - 约束棱边: (节点1, 节点2)\n";
    
    pm.closeProject();
}

/**
 * @brief 测试4: 三维MeshQuery查询
 */
void test_3d_mesh_query() {
    std::cout << "\n========== 测试4: 三维MeshQuery查询 ==========\n";
    
    EMMeshData test_mesh;
    test_mesh.nodes = {
        {1, 0.0, 0.0, 0.0, 1},     // 区域1
        {2, 1.0, 0.0, 0.0, 1},     // 区域1
        {3, 0.5, 0.866, 0.0, 2},   // 区域2
        {4, 0.5, 0.289, 0.816, 2}  // 区域2
    };
    
    test_mesh.elements = {
        {1, {1, 2, 3, 4}, ElemType::TET4, DOFType::MIXED_AV, 1, 1}  // 使用全局节点ID
    };
    
    // 获取单元节点坐标
    auto nodes = MeshQuery::get_element_nodes(test_mesh.elements[0], test_mesh.nodes);
    assert(nodes.size() == 4);
    
    // 验证三维坐标
    assert(nodes[3].z > 0);  // 节点4有Z坐标
    
    // 统计各区域的节点数
    auto region1_nodes = MeshQuery::get_nodes_in_region(1, test_mesh.nodes);
    auto region2_nodes = MeshQuery::get_nodes_in_region(2, test_mesh.nodes);
    
    assert(region1_nodes.size() == 2);
    assert(region2_nodes.size() == 2);
    
    std::cout << "✓ 三维MeshQuery查询验证通过\n";
    std::cout << "  - 区域1节点数: " << region1_nodes.size() << "\n";
    std::cout << "  - 区域2节点数: " << region2_nodes.size() << "\n";
}

/**
 * @brief 测试5: BndType + DOFType 组合设计验证
 * @details 验证组合设计替代独立EMBoundaryType枚举的正确性
 */
void test_boundary_combination() {
    std::cout << "\n========== 测试5: BndType + DOFType 组合设计 ==========\n";
    
    // 验证BndType枚举可用（复用tool层）
    assert(static_cast<int>(BndType::DIRICHLET) == 0);
    assert(static_cast<int>(BndType::NEUMANN) == 1);
    assert(static_cast<int>(BndType::PERFECT_E) == 12);
    assert(static_cast<int>(BndType::PERFECT_H) == 13);
    
    // 验证DOFType枚举可用（fe_em层）
    assert(DOFType::SCALAR_ONLY != DOFType::VECTOR_EDGE_ONLY);
    assert(DOFType::MIXED_AV != DOFType::SCALAR_ONLY);
    
    // 验证典型组合语义正确性
    EMBoundaryMarker test_marker;
    test_marker.bnd_type = BndType::DIRICHLET;
    test_marker.dof_type = DOFType::SCALAR_ONLY;
    test_marker.target_ids = std::vector<int>{1, 2};  // 节点标记
    
    assert(test_marker.bnd_type == BndType::DIRICHLET);      // 狄利克雷边界
    assert(test_marker.dof_type == DOFType::SCALAR_ONLY);     // 标量DOF
    assert(std::holds_alternative<std::vector<int>>(test_marker.target_ids));  // 节点列表
    
    std::cout << "✓ BndType + DOFType 组合设计验证通过\n";
    std::cout << "  - 复用tool::BndType (23种通用类型)\n";
    std::cout << "  - 复用fe_em::DOFType (3种DOF类型)\n";
    std::cout << "  - 组合示例: DIRICHLET + SCALAR_ONLY + 节点标记\n";
}

// ==================== 主函数 ====================

int main() {
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║  三维涡流场A-V混合网格模块测试 (TET4)        ║\n";
    std::cout << "║  方案C架构 + 导体材料 + 理想导体边界        ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n";
    
    try {
        test_tet4_av_mixed_dof();                 // 测试1: A-V混合DOF
        test_tet4_element_geometry();              // 测试2: TET4几何定义
        test_conductor_material_and_boundary();     // 测试3: 导体+PEC边界（核心）
        test_3d_mesh_query();                       // 测试4: 3D查询
        test_boundary_combination();                  // 测试5: BndType+DOFType组合设计
        
        std::cout << "\n╔══════════════════════════════════════════════╗\n";
        std::cout << "║  ✅ 所有三维测试用例通过！                   ║\n";
        std::cout << "║  A-V混合格式 + 导体材料 + PEC边界验证成功   ║\n";
        std::cout << "╚══════════════════════════════════════════════╝\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}
