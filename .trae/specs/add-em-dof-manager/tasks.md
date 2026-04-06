# Tasks

- [x] Task 1: 创建核心公共数据结构 `em_dof_data.hpp`
  - [x] 在 `include/fe_em/` 下创建头文件，定义 Local2Global 结构体
  - [x] 实现 indices 向量、element_id、num_scalar_dofs、num_vector_dofs、elem_type 字段
  - [x] 实现 is_mixed()、get_scalar_indices()、get_vector_indices() 辅助方法
  - [x] 添加完整 Doxygen 注释，遵循项目编码规范

- [x] Task 2: 实现全局棱边ID生成器 `GlobalEdgeIDGenerator`
  - [x] 在 `include/fe_em/` 下创建 `global_edge_id_generator.hpp` 头文件
  - [x] 在 `src/fe_em/` 下创建 `global_edge_id_generator.cpp` 实现文件
  - [x] 实现 generate() 方法：遍历单元棱边，用 unordered_map 去重分配连续ID
  - [x] 仅处理 VECTOR_EDGE_ONLY 和 MIXED_AV 类型单元
  - [x] 棱边唯一表示：std::tuple<int,int>（两个节点ID升序排列）
  - [x] 输出 num_global_edges 和 elem_local_to_global_edge 映射表

- [x] Task 3: 实现DOF管理核心类 `EMDOFManager`
  - [x] 在 `include/fe_em/` 下创建 `em_dof_manager.hpp` 头文件
  - [x] 在 `src/fe_em/` 下创建 `em_dof_manager.cpp` 实现文件
  - [x] 实现四步编号流程：预编号 → 标记约束 → 分块重编号 → 生成映射表
  - [x] 预编号：标量节点DOF（SCALAR_ONLY/MIXED_AV）+ 矢量棱边DOF（VECTOR_EDGE_ONLY/MIXED_AV）
  - [x] 约束标记：处理 DIRICHLET_SCALAR_NODE 和 DIRICHLET_EDGE_EDGE 边界条件
  - [x] 分块重编号：消去法，自由DOF连续编号，约束DOF置-1
  - [x] 映射表生成：按 DOFType 分别处理 SCALAR_ONLY / VECTOR_EDGE_ONLY / MIXED_AV 三类

- [x] Task 4: 编写测试用例1 — PRISM6 MIXED_AV混合格式测试
  - [x] 在 `test/` 下创建 `test_prism6_mixed_av.cpp`
  - [x] 手动构建1个 PRISM6 单元（6节点，MIXED_AV类型）
  - [x] 标记1个节点为Dirichlet边界（标量V固定），1条棱边为Dirichlet边界（矢量A固定）
  - [x] 验证：全局棱边ID生成正确性（9条棱边）
  - [x] 验证：DOF编号正确性（15个局部DOF → 自由+约束分布正确）
  - [x] 验证：Local2Global 映射表拼接正确性（前6标量+后9矢量）

- [x] Task 5: 编写测试用例2 — PYRAMID5_EDGE纯棱边测试
  - [x] 在 `test/` 下创建 `test_pyramid5_edge.cpp`
  - [x] 手动构建1个 PYRAMID5_EDGE 单元（5节点8棱边，VECTOR_EDGE_ONLY类型）
  - [x] 标记1条棱边为Dirichlet边界（A固定）
  - [x] 验证：全局棱边ID生成正确性（8条棱边）
  - [x] 验证：DOF编号正确性（8个局部DOF → 7自由+1约束）
  - [x] 验证：Local2Global 映射表正确性

- [x] Task 6: 编写测试用例3 — PRISM15纯标量测试
  - [x] 在 `test/` 下创建 `test_prism15_scalar.cpp`
  - [x] 手动构建1个 PRISM15 单元（15节点，SCALAR_ONLY类型）
  - [x] 标记1个节点为Dirichlet边界（固定值）
  - [x] 验证：DOF编号正确性（15个局部DOF → 14自由+1约束）
  - [x] 验证：Local2Global 映射表正确性（全标量，无矢量部分）

- [x] Task 7: 更新CMake构建配置
  - [x] 更新 `test/CMakeLists.txt`，添加3个新测试目标的编译配置
  - [x] 配置 include_directories 和 target_link_libraries（链接 fe_em_lib, tool_lib, spdlog_lib）
  - [x] 添加 enable_testing() 和 add_test() 注册

# Task Dependencies
- [Task 2] 依赖 [Task 1]（GlobalEdgeIDGenerator 使用 EMMeshData 数据结构）
- [Task 3] 依赖 [Task 1] 和 [Task 2]（EMDOFManager 依赖 Local2Global 和 GlobalEdgeIDGenerator）
- [Task 4] 依赖 [Task 3]（测试依赖完整的DOF管理功能）
- [Task 5] 依赖 [Task 3]
- [Task 6] 依赖 [Task 3]
- [Task 7] 依赖 [Task 4], [Task 5], [Task 6]
