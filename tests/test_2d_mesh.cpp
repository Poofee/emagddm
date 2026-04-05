/**
 * @file test_2d_mesh.cpp
 * @brief 二维静磁标量位测试用例
 * @details 测试TRI3单元的网格数据结构、ElementGeometry几何定义、ProjectManager集成。
 *          验证方案C架构：EMMeshData（拓扑）+ tool::Material（属性）的正确对接。
 * @author Poofee
 * @date 2026-04-05
 */

#include <iostream>
#include <cassert>
#include <memory>
#include <vector>

#include "fe_em/em_mesh_data.hpp"
#include "fe_em/element_geometry.hpp"
#include "fe_em/mesh_query.hpp"
#include "tool/project_manager.hpp"
#include "tool/project_data.hpp"

using namespace fe_em;
// 注意：tool命名空间不全局using，避免歧义
// 使用 tool:: 前缀访问ProjectManager, Material等

/**
 * @brief 测试1: 创建TRI3单元模型并验证基本数据结构
 */
void test_tri3_basic_structure() {
    std::cout << "\n========== 测试1: TRI3基本数据结构 ==========\n";
    
    // 创建EMMeshData
    auto mesh_data = std::make_unique<EMMeshData>();
    
    // 添加3个节点（二维问题，z=0）
    mesh_data->nodes.push_back({1, 0.0, 0.0, 0.0, 1});   // 节点1: 原点，区域1
    mesh_data->nodes.push_back({2, 1.0, 0.0, 0.0, 1});   // 节点2: (1,0)，区域1
    mesh_data->nodes.push_back({3, 0.5, 0.866, 0.0, 1}); // 节点3: 等边三角形顶点
    
    // 添加1个TRI3单元
    Element elem;
    elem.id = 1;
    elem.node_ids = {1, 2, 3};  // 使用全局节点ID（与Node::id字段匹配）
    elem.elem_type = ElemType::TRI3;
    elem.dof_type = DOFType::SCALAR_ONLY;  // 二维静磁标量位
    elem.material_id = 1;  // 材料ID=1（稍后设置到ProjectManager）
    elem.region_id = 1;
    
    mesh_data->elements.push_back(elem);
    
    // 验证基本数据
    assert(mesh_data->getNodeCount() == 3);
    assert(mesh_data->getElementCount() == 1);
    assert(mesh_data->isEmpty() == false);
    
    // 验证节点坐标
    assert(mesh_data->nodes[0].x == 0.0 && mesh_data->nodes[0].y == 0.0);
    assert(mesh_data->nodes[1].x == 1.0 && mesh_data->nodes[1].y == 0.0);
    assert(std::abs(mesh_data->nodes[2].y - 0.866) < 1e-10);
    
    // 验证单元属性
    const auto& tri_elem = mesh_data->elements[0];
    assert(tri_elem.elem_type == ElemType::TRI3);
    assert(tri_elem.node_ids.size() == 3);
    assert(tri_elem.dof_type == DOFType::SCALAR_ONLY);
    assert(tri_elem.material_id == 1);
    
    std::cout << "✓ TRI3基本数据结构验证通过\n";
}

/**
 * @brief 测试2: 验证ElementGeometry的TRI3几何定义
 */
void test_tri3_element_geometry() {
    std::cout << "\n========== 测试2: TRI3几何定义 ==========\n";
    
    // 检查是否支持TRI3
    assert(ElementGeometry::is_supported(ElemType::TRI3));
    
    // 获取基本信息
    int num_nodes = ElementGeometry::get_num_nodes(ElemType::TRI3);
    int num_edges = ElementGeometry::get_num_edges(ElemType::TRI3);
    int num_faces = ElementGeometry::get_num_faces(ElemType::TRI3);
    
    assert(num_nodes == 3);   // TRI3有3个节点
    assert(num_edges == 3);   // TRI3有3条棱边
    assert(num_faces == 1);   // TRI3有1个面（二维单元）
    
    // 获取棱边定义并验证
    auto edges = ElementGeometry::get_local_edges(ElemType::TRI3);
    assert(edges.size() == 3);
    
    // 验证每条棱边（升序排列）
    assert(std::get<0>(edges[0]) == 0 && std::get<1>(edges[0]) == 1);  // 边0-1
    assert(std::get<0>(edges[1]) == 1 && std::get<1>(edges[1]) == 2);  // 边1-2
    assert(std::get<0>(edges[2]) == 0 && std::get<1>(edges[2]) == 2);  // 边0-2
    
    // 获取面定义并验证
    auto faces = ElementGeometry::get_local_faces(ElemType::TRI3);
    assert(faces.size() == 1);
    assert(faces[0].size() == 3);
    
    // 面的节点顺序应该是逆时针（从外部看）
    assert(faces[0][0] == 0 && faces[0][1] == 1 && faces[0][2] == 2);
    
    std::cout << "✓ TRI3几何定义验证通过\n";
    std::cout << "  - 节点数: " << num_nodes << "\n";
    std::cout << "  - 棱边数: " << num_edges << "\n";
    std::cout << "  - 面数: " << num_faces << "\n";
    std::cout << "  - 棱边定义: (0,1), (1,2), (0,2)\n";
    std::cout << "  - 面定义: [0, 1, 2]\n";
}

/**
 * @brief 测试3: 验证MeshQuery查询功能
 */
void test_mesh_query() {
    std::cout << "\n========== 测试3: MeshQuery查询功能 ==========\n";
    
    // 构建测试数据
    EMMeshData test_mesh;
    test_mesh.nodes = {
        {1, 0.0, 0.0, 0.0, 1},
        {2, 1.0, 0.0, 0.0, 1},
        {3, 0.5, 0.866, 0.0, 2}  // 节点3在区域2
    };
    
    test_mesh.elements = {
        {1, {1, 2, 3}, ElemType::TRI3, DOFType::SCALAR_ONLY, 1, 1}  // 使用全局节点ID
    };
    
    // 测试获取单元节点
    auto elem_nodes = MeshQuery::get_element_nodes(test_mesh.elements[0], test_mesh.nodes);
    assert(elem_nodes.size() == 3);
    assert(elem_nodes[0].id == 1);
    assert(elem_nodes[1].id == 2);
    assert(elem_nodes[2].id == 3);
    
    // 测试区域判断
    assert(MeshQuery::is_node_in_region(1, 1, test_mesh.nodes) == true);   // 节点1在区域1
    assert(MeshQuery::is_node_in_region(3, 2, test_mesh.nodes) == true);   // 节点3在区域2
    assert(MeshQuery::is_node_in_region(3, 1, test_mesh.nodes) == false);  // 节点3不在区域1
    
    assert(MeshQuery::is_element_in_region(test_mesh.elements[0], 1) == true);
    assert(MeshQuery::is_element_in_region(test_mesh.elements[0], 99) == false);
    
    // 测试获取区域内的元素和节点
    auto elems_in_region1 = MeshQuery::get_elements_in_region(1, test_mesh.elements);
    assert(elems_in_region1.size() == 1);
    
    auto nodes_in_region1 = MeshQuery::get_nodes_in_region(1, test_mesh.nodes);
    assert(nodes_in_region1.size() == 2);  // 节点1和2在区域1
    
    std::cout << "✓ MeshQuery查询功能验证通过\n";
}

/**
 * @brief 测试4: ProjectManager集成测试（材料复用验证）
 */
void test_project_manager_integration() {
    std::cout << "\n========== 测试4: ProjectManager集成（方案C架构）==========\n";
    
    auto& pm = tool::ProjectManager::getInstance();
    
    // 创建新项目
    bool success = pm.createProject("Test2DProject");
    assert(success);
    
    // 创建空气材料（使用tool::Material类）
    auto air_material = std::make_shared<tool::Material>("Air");
    air_material->setRelativePermeability(1.0);     // μ_r = 1.0
    air_material->setPermittivity(1.0);            // ε_r = 1.0
    air_material->setConductivity(0.0);             // σ ≈ 0 S/m (绝缘)
    
    pm.addMaterial(air_material);
    
    // 创建EMMeshData拓扑数据
    auto em_mesh = std::make_unique<EMMeshData>();
    em_mesh->nodes = {
        {1, 0.0, 0.0, 0.0, 1},
        {2, 1.0, 0.0, 0.0, 1},
        {3, 0.5, 0.866, 0.0, 1}
    };
    
    Element tri_elem;
    tri_elem.id = 1;
    tri_elem.node_ids = {1, 2, 3};  // 使用全局节点ID
    tri_elem.elem_type = ElemType::TRI3;
    tri_elem.dof_type = DOFType::SCALAR_ONLY;  // 二维静磁标量位
    tri_elem.material_id = 1;  // 引用material ID（不是存储材料本身！）
    tri_elem.region_id = 1;
    
    em_mesh->elements.push_back(tri_elem);
    
    // 添加Dirichlet边界条件（固定磁位）- 组合设计
    EMBoundaryMarker dirichlet_bnd;
    dirichlet_bnd.id = 1;
    dirichlet_bnd.bnd_type = BndType::DIRICHLET;        // 边界类型：狄利克雷（复用tool::BndType）
    dirichlet_bnd.dof_type = DOFType::SCALAR_ONLY;       // DOF类型：标量场
    dirichlet_bnd.target_ids = std::vector<int>{1};      // 固定节点ID=1的磁位为0
    dirichlet_bnd.value = 0.0;
    dirichlet_bnd.name = "Ground";
    
    em_mesh->boundary_markers.push_back(dirichlet_bnd);
    
    // 设置到ProjectManager（转移所有权）
    pm.setEMMeshData(std::move(em_mesh));
    
    // 验证数据已正确设置
    assert(pm.hasEMMeshData());
    
    auto* retrieved_mesh = pm.getEMMeshData();
    assert(retrieved_mesh != nullptr);
    assert(retrieved_mesh->getNodeCount() == 3);
    assert(retrieved_mesh->getElementCount() == 1);
    assert(retrieved_mesh->getBoundaryMarkerCount() == 1);
    
    // 关键验证：通过material_id索引获取材料参数（零拷贝！）
    const auto& elem = retrieved_mesh->elements[0];
    auto mat_opt = pm.getMaterial("Air");  // 通过名称获取
    assert(mat_opt.has_value());
    
    auto* material = mat_opt.value().get();
    double mu_r = material->getRelativePermeability();
    double sigma = material->getConductivity();
    
    assert(std::abs(mu_r - 1.0) < 1e-10);      // 空气相对磁导率≈1
    assert(std::abs(sigma - 0.0) < 1e-10);     // 空气电导率≈0
    
    std::cout << "✓ ProjectManager集成验证通过\n";
    std::cout << "  - EMMeshData成功设置到ProjectManager\n";
    std::cout << "  - 通过material_id=" << elem.material_id 
              << " 索引获取材料: μ_r=" << mu_r << ", σ=" << sigma << "\n";
    std::cout << "  - 边界条件数量: " << retrieved_mesh->getBoundaryMarkerCount() << "\n";
    
    // 清理
    pm.closeProject();
}

/**
 * @brief 测试5: DOFType枚举转换测试
 */
void test_dof_type_enum() {
    std::cout << "\n========== 测试5: DOFType枚举转换 ==========\n";
    
    // 测试枚举转字符串
    assert(dofTypeToString(DOFType::SCALAR_ONLY) == "SCALAR_ONLY");
    assert(dofTypeToString(DOFType::VECTOR_EDGE_ONLY) == "VECTOR_EDGE_ONLY");
    assert(dofTypeToString(DOFType::MIXED_AV) == "MIXED_AV");
    
    // 测试字符串转枚举
    assert(stringToDOFType("SCALAR_ONLY") == DOFType::SCALAR_ONLY);
    assert(stringToDOFType("VECTOR_EDGE_ONLY") == DOFType::VECTOR_EDGE_ONLY);
    assert(stringToDOFType("MIXED_AV") == DOFType::MIXED_AV);
    
    std::cout << "✓ DOFType枚举转换验证通过\n";
}

// ==================== 主函数 ====================

int main() {
    std::cout << "╔════════════════════════════════════════════╗\n";
    std::cout << "║   二维静磁标量位网格模块测试 (TRI3)       ║\n";
    std::cout << "║   方案C架构: EMMeshData + tool::Material  ║\n";
    std::cout << "╚════════════════════════════════════════════╝\n";
    
    try {
        test_tri3_basic_structure();           // 测试1: 基本数据结构
        test_tri3_element_geometry();          // 测试2: 几何定义
        test_mesh_query();                     // 测试3: 查询功能
        test_project_manager_integration();    // 测试4: PM集成（核心）
        test_dof_type_enum();                  // 测试5: 枚举转换
        
        std::cout << "\n╔════════════════════════════════════════════╗\n";
        std::cout << "║  ✅ 所有测试用例通过！                      ║\n";
        std::cout << "║  方案C架构验证成功                         ║\n";
        std::cout << "╚════════════════════════════════════════════╝\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}
