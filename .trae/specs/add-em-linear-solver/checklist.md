# Checklist - 工业级电磁场有限元线性求解器模块

## 基础架构验证

- [ ] **C1**: em_linear_solver.h 头文件存在且包含完整的 SolverStatus 枚举定义（5个状态值）
- [ ] **C2**: SolverResult 结构体字段完整（x, status, iterations, residual_norm, solve_time_ms, error_msg）
- [ ] **C3**: EMLinearSolverBase 抽象基类声明4个纯虚函数（set_matrix/solve/get_solver_name/clear）
- [ ] **C4**: 所有头文件使用 #pragma once 保护
- [ ] **C5**: 头文件仅包含文件名，无路径前缀

## 矩阵转换工具验证

- [ ] **C6**: SparseConverter::to_eigen() 方法实现正确，支持 CsrMatrix<double> → Eigen::SparseMatrix<double>
- [ ] **C7**: SparseConverter::from_eigen() 方法实现正确，支持反向转换
- [ ] **C8**: 往返转换精度测试通过：||A_original - A_roundtrip||_F < 1e-15
- [ ] **C9**: 维度校验功能正常工作（不匹配时返回错误或抛异常）
- [ ] **C10**: 大规模矩阵转换性能可接受（10000×10000 矩阵 < 500ms）

## 直接求解器验证

### 多后端基础设施
- [ ] **C10**: DirectBackendType 枚举定义完整（EIGEN/SUPERLU/MUMPS 三个值）
- [ ] **C10a**: DirectBackendManager 类实现完整（isBackendAvailable/getBackendName/getAvailableBackends 三个静态方法）
- [ ] **C10b**: em_solver_backends.hpp 条件编译配置文件存在，正确处理 #ifdef HAVE_SUPERLU / #ifdef HAVE_MUMPS
- [ ] **C10c**: DirectBackendManager::isBackendAvailable(DirectBackendType::EIGEN) 始终返回 true
- [ ] **C10d**: 未编译 SuperLU/MUMPS 时，对应的 isBackendAvailable 返回 false
- [ ] **C10e**: getBackendName() 对每个后端类型返回友好的中英文字符串名称

### SymmetricDirectSolver（多后端支持）
- [ ] **C11**: 类正确继承 EMLinearSolverBase 并实现所有纯虚函数
- [ ] **C12**: EIGEN 后端使用 Eigen::SimplicialLLT 作为底层求解器（或有 MKL Pardiso 条件编译分支）
- [ ] **C12a**: 构造函数支持 DirectBackendType 参数指定初始后端
- [ ] **C12b**: set_backend() 方法可动态切换后端（运行时切换）
- [ ] **C12c**: 切换到不可用后端时自动回退到 EIGEN 并输出 WARNING 日志
- [ ] **C13**: 15×15 SPD 测试矩阵求解成功，误差 < 1e-10（EIGEN 后端验证）
- [ ] **C13a**: SuperLU 后端启用时（#ifdef HAVE_SUPERLU），SuperLU 求解结果与 Eigen 一致（误差 < 1e-10）
- [ ] **C13b**: MUMPS 后端启用时（#ifdef HAVE_MUMPS），MUMPS 求解结果与 Eigen 一致（误差 < 1e-10）
- [ ] **C14**: set_matrix() 执行分解并缓存，后续 solve() 复用结果（所有后端一致）
- [ ] **C15**: 正定性检测失败时返回 NUMERICAL_ERROR 状态（所有后端一致）
- [ ] **C16**: 对称性校验选项可用且正常工作（set_symmetry_tolerance 方法存在）
- [ ] **C17**: 求解时间统计使用 std::chrono 高精度计时（所有后端均统计）

### SymmetricIndefiniteDirectSolver（多后端支持）
- [ ] **C18**: 类正确继承 EMLinearSolverBase 并实现所有纯虚函数
- [ ] **C19**: EIGEN 后端使用 Eigen::SimplicialLDLT 作为底层求解器
- [ ] **C19a**: 支持多后端切换（同 C12b/C12c 验证标准）
- [ ] **C20**: 8×8 半正定奇异矩阵处理正常（不崩溃，返回合理状态）（EIGEN 后端）
- [ ] **C20a**: SuperLU/MUMPS 后端对奇异矩阵的鲁棒性测试通过（如已编译）
- [ ] **C21**: 零空间维度检测/报告机制工作正常（get_null_space_dimension 方法）
- [ ] **C22**: 正则化参数配置接口可用（set_regularization_epsilon 方法）

### GeneralDirectSolver（多后端支持）
- [ ] **C23**: 类正确继承 EMLinearSolverBase 并实现所有纯虚函数
- [ ] **C24**: EIGEN 后端使用 Eigen::SparseLU 作为底层求解器
- [ ] **C24a**: 支持多后端切换（同 C12b/C12c 验证标准）
- [ ] **C25**: 非对称测试矩阵求解成功（EIGEN 后端）
- [ ] **C25a**: SuperLU/MUMPS 后端非对称矩阵求解成功且结果一致（如已编译）
- [ ] **C26**: 病态矩阵警告机制工作（条件数 > 1e10 时输出 WARNING 日志）

## 迭代求解器验证

### CGSolver
- [ ] **C27**: 类正确继承 EMLinearSolverBase 并实现所有纯虚函数
- [ ] **C28**: 使用 Eigen::ConjugateGradient 作为底层算法
- [ ] **C29**: tolerance 和 max_iterations 参数可配置
- [ ] **C30**: Jacobi 预条件子实现有效（加速收敛 >30%）
- [ ] **C31**: ILU0 预条件子实现有效
- [ ] **C32**: SPD 测试矩阵收敛至 tolerance=1e-8
- [ ] **C33**: 迭代次数和最终残差正确记录在 SolverResult 中
- [ ] **C34**: 达到最大迭代次数时返回 MAX_ITER_REACHED 状态

### BiCGSTABSolver
- [ ] **C35**: 类正确继承 EMLinearSolverBase 并实现所有纯虚函数
- [ ] **C36**: 使用 Eigen::BiCGSTAB 作为底层算法
- [ ] **C37**: 不定矩阵测试用例稳定收敛
- [ ] **C38**: 发散检测机制工作正常（残差异常增长时提前终止）

### AMG 预条件子
- [ ] **C39**: ScalarAMG 标量AMG预条件子框架实现
- [ ] **C40**: AMG 核心组件完整（粗化、插值、限制、光滑、V循环/W循环）
- [ ] **C41**: AMG 接口兼容 Eigen 迭代求解器的预条件子规范
- [ ] **C42**: HiptmairAMG 接口预留（类声明 + TODO 注释 + 设计文档说明）

## 工厂模式验证

- [ ] **C43**: EMSolverFactory::create_solver() 正确创建5种求解器类型
- [ ] **C44**: 返回 std::unique_ptr<EMLinearSolverBase> 智能指针
- [ ] **C45**: 创建的实例可直接调用基类接口（多态性验证）
- [ ] **C46**: 场景自动适配方法（可选）基于 MatrixAttribute 正确选择求解器

## 测试用例验证

### Test Case 1: 标量场泊松方程（PRISM15单元）
- [ ] **C47**: test_linear_solver_scalar.cpp 文件存在且可编译
- [ ] **C48**: PRISM15 单元刚度矩阵构建正确（15×15 对称正定）
- [ ] **C49**: SymmetricDirectSolver 求解状态为 SUCCESS
- [ ] **C50**: 解向量误差 ||x - x_exact|| < 1e-6
- [ ] **C51**: 残差范数 < 1e-10
- [ ] **C52**: 多次 solve() 调用结果一致（分解复用验证）

### Test Case 2: 棱边单元静磁场（PYRAMID5_EDGE单元）
- [ ] **C53**: test_linear_solver_edge.cpp 文件存在且可编译
- [ ] **C54**: PYRAMID5_EDGE 单元刚度矩阵构建正确（8×8 半正定）
- [ ] **C55**: SymmetricIndefiniteDirectSolver 求解完成（状态 SUCCESS 或明确错误信息）
- [ ] **C56**: 解向量满足残差要求
- [ ] **C57**: 奇异性处理日志输出正常

### Test Case 3: A-V混合涡流场（PRISM6混合单元）
- [ ] **C58**: test_linear_solver_mixed_av.cpp 文件存在且可编译
- [ ] **C59**: PRISM6_MIXED_AV 单元刚度矩阵构建正确（18×18 不定）
- [ ] **C60**: GeneralDirectSolver 求解状态为 SUCCESS
- [ ] **C61**: BiCGSTABSolver 求解状态为 SUCCESS
- [ ] **C62**: 两个求解器结果一致性 ||x_direct - x_iterative|| < 1e-6
- [ ] **C63**: BiCGSTAB 迭代次数 < 500 且残差 < 1e-8

## 编译与集成验证

- [ ] **C64**: CMakeLists.txt 已添加3个新测试目标配置
- [ ] **C65**: cmake 配置无错误或警告
- [ ] **C66**: 全量编译成功（包括 numeric_lib 更新后的源文件）
- [ ] **C67**: 3个测试目标均可独立编译链接
- [ ] **C68**: ctest 运行全部 PASS（0 failures）
- [ ] **C69**: 无编译器警告（-Wall -Wextra -Werror 通过）

## 代码质量验证

- [ ] **C70**: 所有公开函数有 Doxygen 风格注释（@brief @param @return @details）
- [ ] **C71**: 关键算法步骤有行内注释说明设计意图
- [ ] **C72**: 使用 FEEM_INFO/FEEM_DEBUG/FEEM_ERROR 宏输出日志（禁止 printf/cout 用于日志）
- [ ] **C73**: 无全局变量（所有状态封装在类实例中）
- [ ] **C74**: 内存管理安全（智能指针管理动态资源，无裸指针泄漏）
- [ ] **C75**: 输入参数合法性校验覆盖所有公开接口
- [ ] **C76**: 错误处理完善（使用 SolverStatus 枚举，禁止静默失败）
- [ ] **C77**: 命名规范符合项目标准：
  - 类名 PascalCase (如 SymmetricDirectSolver)
  - 函数 camelCase (如 set_matrix, solve)
  - 变量 snake_case (如 tolerance_, max_iterations_)
  - 常量 UPPER_SNAKE_CASE (如 MAX_DEFAULT_ITERATIONS)
- [ ] **C78**: 单个文件代码长度合理（头文件 < 500 行，源文件按需拆分）
- [ ] **C79**: 单个函数长度 < 100 行（过长函数已拆分）
- [ ] **C80**: 无平台专属语法或编译器特定特性（保证跨平台兼容性）

## 性能基准验证（可选但推荐）

- [ ] **C81**: SymmetricDirectSolver 求解 1000×1000 SPD 矩阵耗时 < 100ms
- [ ] **C82**: CGSolver 配合 Jacobi 预条件子在 10000×1000 矩阵上 500 次迭代内收敛
- [ ] **C83**: 瞬态多步复用场景下，第2-N次 solve 耗时比首次降低 >80%（分解缓存生效）

---

## 验证执行顺序建议

### Phase 1: 基础架构检查 (C1-C5)
快速确认头文件结构和基类定义是否完整。

### Phase 2: 核心功能验证 (C6-C26)
依次验证矩阵转换工具 → 3个直接求解器的基础功能。

### Phase 3: 迭代求解器验证 (C27-C42)
验证 CG/BiCGSTAB 求解器和 AMG 预条件子框架。

### Phase 4: 工厂与集成验证 (C43-C46, C64-C69)
验证工厂模式和 CMake 编译集成。

### Phase 5: 测试用例全量验证 (C47-C63)
运行3个测试用例，确认全电磁场景覆盖。

### Phase 6: 代码质量审查 (C70-C80)
静态代码分析，确保符合项目编码规范。

### Phase 7: 性能验证 (C81-C83) [可选]
性能基准测试，确认满足工业级性能要求。
