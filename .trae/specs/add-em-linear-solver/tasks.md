# Tasks - 工业级电磁场有限元线性求解器模块

## 任务概览
本模块实现完整的电磁场线性求解器系统，包含抽象基类、直接/迭代求解器实现、矩阵转换工具、工厂模式和测试用例。

---

## Task 1: 创建求解器抽象基类与数据结构定义 ✅
**文件**: `include/numeric/em_linear_solver.h`
**优先级**: 高（所有后续任务依赖此基类）
**预估工作量**: 小
**状态**: 已完成 (2026-04-06)

### 子任务
- [x] 1.1 定义 SolverStatus 枚举（SUCCESS/NUMERICAL_ERROR/DIVERGED/MAX_ITER_REACHED/INVALID_INPUT）
- [x] 1.2 定义 SolverResult 结构体（解向量x、状态码、迭代次数、残差范数、耗时、错误信息）
- [x] 1.3 定义 EMLinearSolverBase 抽象基类（纯虚函数：set_matrix/solve/get_solver_name/clear）
- [x] 1.4 添加 Doxygen 风格完整文档注释
- [x] 1.5 确保头文件仅包含文件名（无路径），符合项目规范

### 验证标准
- [x] 头文件可被其他模块正常 include
- [x] 编译无警告无错误
- [x] 抽象基类接口清晰，符合用户需求规格

---

## Task 2: 实现CSR矩阵与Eigen稀疏矩阵转换工具 ✅
**文件**:
- `include/numeric/em_sparse_converter.h`（声明）
- `src/numeric/em_sparse_converter.cpp`（实现）
**优先级**: 高（直接/迭代求解器均依赖转换功能）
**预估工作量**: 中
**状态**: 已完成 (2026-04-06)

### 子任务
- [x] 2.1 创建 SparseConverter 工具类（静态方法设计）
- [x] 2.2 实现 to_eigen() 方法：CsrMatrix<double> → Eigen::SparseMatrix<double>
  - 使用 Eigen::Triplets 构建稀疏矩阵（零拷贝优化）
  - 添加维度校验和异常处理
  - 支持对称矩阵优化存储（仅存上三角或下三角）
- [x] 2.3 实现 from_eigen() 方法：Eigen::SparseMatrix<double> → CsrMatrix<double>
  - 提取 Eigen 内部 CSR 数据结构（rowPtr/colIdx/values）
  - 保证数据一致性，支持压缩存储模式
- [x] 2.4 添加辅助方法：
  - validate_matrix_dimensions() - 维度一致性校验
  - get_conversion_stats() - 转换统计信息（耗时、非零元数量等）
- [x] 2.5 性能优化：对于大型矩阵（>10000 DOF）使用 reserve() 预分配空间
- [x] 2.6 单元测试验证往返转换精度（误差 < 1e-15）

### 验证标准
- [x] to_eigen/from_eigen 往返转换数据完全一致
- [x] 支持 15×15 小规模到 10000×10000 大规模矩阵
- [x] 转换耗时在可接受范围内（<100ms for 10K×10K matrix）
- [x] 异常输入返回明确错误信息

---

## Task 3: 实现直接求解器（工业级鲁棒性 + 多后端策略模式）✅
**文件**:
- `include/numeric/em_direct_solvers.h`（DirectBackendType枚举、DirectBackendManager、3个派生类声明）
- `src/numeric/em_direct_solvers.cpp`（实现，含条件编译分支）
**优先级**: 高（核心求解能力）
**预估工作量**: 大（需实现多后端架构）
**状态**: 已完成 (2026-04-06)

### 架构设计：策略模式（Strategy Pattern）
每个直接求解器内部通过 DirectBackendType 枚举选择底层算法后端：
- **EIGEN**: 默认后端，使用 Eigen 内置稀疏求解器（无需额外依赖）
- **SUPERLU**: 可选后端，高性能开源直接求解器（#ifdef HAVE_SUPERLU）
- **MUMPS**: 可选后端，并行分布式求解器（#ifdef HAVE_MUMPS）

### 子任务

#### 3.0 基础设施：后端类型枚举和管理器
- [ ] 3.0.1 定义 DirectBackendType 枚举（EIGEN/SUPERLU/MUMPS）
  - 放置在 em_linear_solver.h 或 em_direct_solvers.h 中
  - 添加 Doxygen 注释说明各后端特点
- [ ] 3.0.2 实现 DirectBackendManager 辅助类
  - isBackendAvailable(DirectBackendType) - 运行时检测后端可用性
    - 基于 #ifdef HAVE_SUPERLU / #ifdef HAVE_MUMPS 编译期宏
    - 未定义时返回 false
  - getBackendName(DirectBackendType) - 返回友好名称字符串
  - getAvailableBackends() - 返回所有已启用后端的 vector
  - 日志输出当前可用后端列表（FEEM_INFO 级别）

#### 3.1 实现 SymmetricDirectSolver 类（支持多后端）
- [ ] 3.1.1 类结构设计
  - 继承 EMLinearSolverBase
  - 成员变量：
    ```cpp
    DirectBackendType backend_type_ = DirectBackendType::EIGEN;
    std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>> eigen_solver_;
    #ifdef HAVE_SUPERLU
    // SuperLU 求解器实例指针
    #endif
    #ifdef HAVE_MUMPS
    // MUMPS 求解器实例指针
    #endif
    bool matrix_set_ = false;
    double symmetry_tol_ = 1e-10;
    ```
- [ ] 3.1.2 构造与配置方法
  - 构造函数支持指定 backend_type 参数
  - set_backend(DirectBackendType) 方法：动态切换后端
    - 切换时清空已有分解结果（clear()）
    - 如果目标后端不可用，输出 WARNING 并回退到 EIGEN
  - set_symmetry_tolerance(double tol) 配置对称性校验容差
- [ ] 3.1.3 EIGEN 后端实现（默认）
  - 使用 Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>
  - set_matrix():
    - 转换 CsrMatrix → Eigen::SparseMatrix（使用 SparseConverter::to_eigen）
    - 执行符号分解 + 数值分解并缓存
    - 捕获 Eigen 分解失败异常，返回 NUMERICAL_ERROR
  - solve():
    - 复用已有分解结果
    - 计算解向量 x = A^{-1} b
    - 统计求解时间（std::chrono）
  - 支持 Intel MKL Pardiso（条件编译 #ifdef EIGEN_USE_MKL_ALL）
- [ ] 3.1.4 SUPERLU 后端实现（可选，#ifdef HAVE_SUPERLU）
  - 使用 SuperLU 的对称矩阵优化接口
  - set_matrix():
    - 将 CSR 数据转换为 SuperLU 内部格式（SuperMatrix）
    - 执行列分析 + 符号分解 + 数值分解
    - 对称矩阵仅存储上三角（节省内存）
  - solve():
    - 调用 SuperLU 求解接口
    - 处理 SuperLU 返回的错误码，映射到 SolverStatus
- [ ] 3.1.5 MUMPS 后端实现（可选，#ifdef HAVE_MUMPS）
  - 使用 MUMPS 分布式并行求解器
  - set_matrix():
    - 初始化 MUMPS 数据结构
    - 设置对称正定矩阵参数 (ICNTL(28)=1, ICNTL(29)=2)
    - 执行分析 + 分解阶段
  - solve():
    - 调用 MUMPS 求解阶段
    - 处理 MUMPS 信息码 (INFO(1), INFO(2))
  - 需 MPI 初始化/终结管理（如果 USE_MPI 启用）
- [ ] 3.1.6 公共功能（与后端无关）
  - 自动检测矩阵正定性（尝试 Cholesky 分解，失败则返回错误）
  - 对称性校验（遍历非零元检查 |A(i,j) - A(j,i)| < tolerance）
  - 性能日志输出（分解时间、求解时间、后端类型标识）

#### 3.2 实现 SymmetricIndefiniteDirectSolver 类（支持多后端）
- [ ] 3.2.1 类结构设计（同 3.1.1，但使用 LDL^T/LU 分解）
- [ ] 3.2.2 EIGEN 后端实现
  - 使用 Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>
  - 处理零主元问题（添加微小扰动 regularization_epsilon_）
  - 奇异矩阵检测（检查 LDL^T 分解中的零对角元数量）
- [ ] 3.2.3 SUPERLU 后端实现（可选）
  - 使用 SuperLU 完全选主元 LU 分解
  - 支持不定对称矩阵的稳定求解
- [ ] 3.2.4 MUMPS 后端实现（可选）
  - 使用 MUMPS 非对称通用模式 + 对称预处理
  - 强大的奇异矩阵处理能力（秩检测、零空间计算）
- [ ] 3.2.5 核心功能
  - 正则化参数配置（set_regularization_epsilon(double)）
  - 零空间维度报告（get_null_space_dimension() 方法）
  - 奇异矩阵警告日志输出

#### 3.3 实现 GeneralDirectSolver 类（支持多后端）
- [ ] 3.3.1 类结构设计（同上，面向非对称矩阵）
- [ ] 3.3.2 EIGEN 后端实现
  - 使用 Eigen::SparseLU<Eigen::SparseMatrix<double>>
  - 列主序 LU 分解，阈值部分选主元（pivot threshold）
  - 病态矩阵检测（估算条件数或检查 pivot 大小）
- [ ] 3.3.3 SUPERLU 后端实现（可选）
  - 使用 SuperLU 完全选主元模式（保证数值稳定性）
  - 适合高病态矩阵（条件数可达 1e15）
- [ ] 3.3.4 MUMPS 后端实现（可选）
  - 使用 MUMPS 非对称并行模式
  - 支持大规模分布式非对称求解
  - 适合复杂耦合场（如涡流场复数矩阵）
- [ ] 3.3.5 核心功能
  - 病态矩阵警告机制（条件数 > 1e10 时 FEEM_WARN 输出）
  - 选主元策略配置（set_pivot_threshold(double)）
  - 数值稳定性增强选项

#### 3.4 公共功能封装（避免代码冗余）
- [ ] 3.4.1 extract_solver_result() 辅助函数
  - 统一构建 SolverResult 结构体
  - 封装耗时统计、残差计算、状态码设置逻辑
- [ ] 3.4.2 validate_input() 输入校验函数
  - 检查矩阵维度合法性（行数 == 列数 > 0）
  - 检查向量长度匹配（b.size() == A.rows()）
  - 检查矩阵是否已设置（matrix_set_ 标志）
  - 不合法时返回 SolverResult{status: INVALID_INPUT, error_msg: "详细原因"}
- [ ] 3.4.3 measure_solve_time() 高精度计时工具
  - 使用 std::chrono::high_resolution_clock
  - 返回毫秒级精度耗时（double solve_time_ms）
  - 封装为 RAII 计时器对象（构造开始计时，析构自动记录）
- [ ] 3.4.4 fallback_to_eigen() 优雅降级函数
  - 当请求的后端不可用时自动切换到 EIGEN
  - 输出 WARNING 日志说明降级原因
  - 返回实际使用的后端类型

#### 3.5 条件编译配置文件（新增 em_solver_backends.hpp）
- [ ] 3.5.1 创建编译期宏定义头文件
  ```cpp
  // em_solver_backends.hpp
  #pragma once
  
  // SuperLU 后端支持（在CMakeLists.txt中通过 -DHAVE_SUPERLU 启用）
  #ifdef HAVE_SUPERLU
    #define EM_SOLVER_HAS_SUPERLU 1
    // 包含SuperLU头文件
    #include "slu_ddefs.h"  // SuperLU Double Precision
  #else
    #define EM_SOLVER_HAS_SUPERLU 0
  #endif
  
  // MUMPS 后端支持（在CMakeLists.txt中通过 -DHAVE_MUMPS 启用）
  #ifdef HAVE_MUMPS
    #define EM_SOLVER_HAS_MUMPS 1
    // 包含MUMPS头文件
    #include "dmumps_c.h"  // MUMPS Double Precision
  #else
    #define EM_SOLVER_HAS_MUMPS 0
  #endif
  ```
- [ ] 3.5.2 文档说明如何在 CMakeLists.txt 中启用各后端
  - SuperLU: option(USE_SUPERLU "Enable SuperLU backend" OFF)
  - MUMPS: option(USE_MUMPS "Enable MUMPS backend" OFF) + 需要 USE_MPI=ON

### 验证标准
- [ ] **EIGEN 后端验证**:
  - SymmetricDirectSolver 正确求解 15×15 SPD 矩阵（误差 < 1e-10）
  - SymmetricIndefiniteDirectSolver 正确处理 8×8 半正定奇异矩阵
  - GeneralDirectSolver 正确求解非对称测试矩阵
- [ ] **多后端切换验证**:
  - set_backend(SUPERLU) 在未编译 SuperLU 时优雅回退到 EIGEN
  - set_backend(MUMPS) 在未编译 MUMPS 时优雅回退到 EIGEN
  - 回退时输出正确的 WARNING 日志
  - 回退后求解结果仍然正确
- [ ] **后端管理器验证**:
  - DirectBackendManager::isBackendAvailable(EIGEN) 始终返回 true
  - DirectBackendManager::getAvailableBackends() 返回正确列表
  - DirectBackendManager::getBackendName() 返回友好字符串
- [ ] **瞬态复用验证**:
  - 所有求解器支持 set_matrix一次，solve多次
  - 第2-N次 solve 耗时显著降低（分解缓存生效）
- [ ] **错误处理验证**:
  - 错误场景正确返回对应状态码（NUMERICAL_ERROR/INVALID_INPUT）
  - 日志输出规范（FEEM_INFO/FEEM_DEBUG/FEEM_ERROR）
- [ ] GeneralDirectSolver 正确求解非对称测试矩阵
- [ ] 所有求解器支持瞬态多步复用（set_matrix一次，solve多次）
- [ ] 错误场景正确返回对应状态码（NUMERICAL_ERROR/INVALID_INPUT）
- [ ] 日志输出规范（FEEM_INFO/FEEM_DEBUG/FEEM_ERROR）

---

## Task 4: 实现迭代求解器（高性能大规模问题）
**文件**:
- `include/numeric/em_iterative_solvers.h`（CGSolver/BiCGSTABSolver/AMG预条件子声明）
- `src/numeric/em_iterative_solvers.cpp`（实现）
**优先级**: 高（大规模问题必备）
**预估工作量**: 大

### 子任务
- [ ] 4.1 实现 CGSolver 类
  - 继承 EMLinearSolverBase
  - 内部使用 Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>
  - 可配置参数：
    - tolerance_ (double, 默认 1e-8)
    - max_iterations_ (int, 默认 1000)
    - preconditioner_type_ (枚举: NONE/JACOBI/ILU0)
  - Jacobi预条件子实现（对角缩放）
  - ILU0预条件子实现（不完全LU分解，零填充）
  - 迭代过程残差监控（每N步记录一次残差历史）
  - 收敛判断逻辑：||r|| / ||b|| < tolerance
  
- [ ] 4.2 实现 BiCGSTABSolver 类
  - 继承 EMLinearSolverBase
  - 内部使用 Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>
  - 可配置参数同 CGSolver
  - ILU0预条件子支持
  - 发散检测机制（残差增长超过阈值时提前终止）
  
- [ ] 4.3 实现标量AMG预条件子 ScalarAMG
  - 继承 Eigen 预条件子基类接口（或提供兼容包装器）
  - 核心算法组件：
    - 粗化策略（Classic Ruge-Stüben 或 Aggregation-based AMG）
    - 插值算子构造（直接插值或多步插值）
    - 限制算子（插值算子的转置）
    - 光滑器选择（Jacobi/Gauss-Seidel）
    - V循环/W循环多层求解
  - 适配 Eigen 迭代求解器的预条件子 apply() 接口
  - 最大层级限制（默认 20 层，防止过度粗化）
  - 最小粗网格尺寸阈值（默认 50 DOF）
  
- [ ] 4.4 预留 HiptmairAMG 接口（H(curl)专用AMG）
  - 声明 HiptmairAMG 类框架（纯虚函数或空实现）
  - 文档说明设计意图和使用场景
  - 预留标量AMG + 梯度修复的关键接口签名
  - 标注 TODO 注释等待后续实现
  
- [ ] 4.5 迭代求解器公共功能
  - IterativeSolverConfig 结构体（统一配置参数）
  - ConvergenceHistory 类（记录收敛曲线数据）
  - 预条件子工厂函数 create_preconditioner()

### 验证标准
- [ ] CGSolver 在 SPD 测试矩阵上正确收敛（误差 < tolerance）
- [ ] BiCGSTABSolver 在不定测试矩阵上稳定收敛
- [ ] Jacobi/ILU0 预条件子有效加速收敛（减少迭代次数 >30%）
- [ ] AMG 预条件子在大型测试矩阵上表现良好（可选基准测试）
- [ ] 所有迭代求解器正确返回迭代次数和最终残差
- [ ] 发散场景正确识别并返回 DIVERGED 状态

---

## Task 5: 实现求解器工厂类
**文件**: `include/numeric/em_solver_factory.h`
**优先级**: 中（便利上层调用）
**预估工作量**: 小

### 子任务
- [ ] 5.1 创建 EMSolverFactory 工厂类
  - 定义 SolverType 枚举（SYMMETRIC_DIRECT/SYMMETRIC_INDEFINITE/GENERAL_DIRECT/CG/BICGSTAB）
  - 实现 static create_solver(SolverType type) 工厂方法
  - 返回 std::unique_ptr<EMLinearSolverBase>
  
- [ ] 5.2 添加便捷创建方法（可选增强）
  - create_solver_for_electrostatic() - 自动选择 SymmetricDirectSolver
  - create_solver_for_magnetostatic_edge() - 自动选择 SymmetricIndefiniteDirectSolver
  - create_solver_for_eddy_current() - 自动选择 BiCGSTABSolver 或 GeneralDirectSolver
  
- [ ] 5.3 场景自动适配逻辑（基于 MatrixAttribute）
  - 根据 MatrixAttribute::is_spd 选择 CG 或 Direct
  - 根据 MatrixAttribute::is_singular 选择 IndefiniteDirectSolver
  - 根据 MatrixAttribute::symmetry 选择对称/非对称求解器

### 验证标准
- [ ] 工厂方法正确创建对应类型的求解器实例
- [ ] 返回的智能指针可直接使用基类接口
- [ ] 无内存泄漏（unique_ptr 自动管理生命周期）

---

## Task 6: 实现测试用例1 - 标量场泊松方程求解（PRISM15单元）
**文件**: `test/test_linear_solver_scalar.cpp`
**优先级**: 高（验证全链路正确性）
**预估工作量**: 中

### 测试内容
- **场景**: 手动构建1个 PRISM15 二阶三棱柱单元（15个节点，DOFType=SCALAR_ONLY）
- **矩阵**: 装配生成 15×15 对称正定刚度矩阵 K
- **右端项**: 构造已知解析解 x_exact = [1, 2, 3, ..., 15]^T，计算 b = K * x_exact
- **求解器**: SymmetricDirectSolver
- **验证点**:
  - [ ] 求解状态为 SUCCESS
  - [ ] ||x_computed - x_exact|| < 1e-6
  - [ ] 残差范数 ||b - A*x_computed|| < 1e-10
  - [ ] 求解时间 > 0 ms（计时正常工作）
  - [ ] 多次 solve() 调用结果一致（分解复用验证）

### 辅助函数
- build_prism15_stiffness_matrix() - 构建 PRISM15 刚度矩阵
- verify_symmetric_positive_definite() - 验证 SPD 属性

---

## Task 7: 实现测试用例2 - 棱边单元静磁场求解（PYRAMID5_EDGE单元）
**文件**: `test/test_linear_solver_edge.cpp`
**优先级**: 高（验证半正定/奇异矩阵处理）
**预估工作量**: 中

### 测试内容
- **场景**: 手动构建1个 PYRAMID5_EDGE 一阶金字塔棱边单元（8条棱边，DOFType=VECTOR_EDGE_ONLY）
- **矩阵**: 装配生成 8×8 对称半正定刚度矩阵 K（含零空间，秩可能 < 8）
- **右端项**: 构造相容右端项 F（确保在值域空间内）
- **求解器**: SymmetricIndefiniteDirectSolver
- **验证点**:
  - [ ] 求解状态为 SUCCESS（或 NUMERICAL_ERROR with clear message if truly singular）
  - [ ] 解向量满足 K*x ≈ F（残差合理）
  - [ ] 奇异性检测日志输出正常
  - [ ] 与 GeneralDirectSolver 结果对比一致性

### 辅助函数
- build_pyramid5_edge_stiffness_matrix() - 构建 PYRAMID5_EDGE 刚度矩阵
- check_null_space_dimension() - 检测零空间维度

---

## Task 8: 实现测试用例3 - A-V混合涡流场求解（PRISM6混合单元）
**文件**: `test/test_linear_solver_mixed_av.cpp`
**优先级**: 高（验证不定矩阵和多求解器对比）
**预估工作量**: 中

### 测试内容
- **场景**: 手动构建1个 PRISM6 一阶三棱柱混合单元（6节点+9棱边=18DOF，DOFType=MIXED_AV）
- **矩阵**: 装配生成 18×18 对称不定刚度矩阵 K（A-V耦合块结构）
- **右端项**: 构造已知解的右端项 F = K * x_exact
- **求解器1**: GeneralDirectSolver
- **求解器2**: BiCGSTABSolver（配合 ILU0 预条件子）
- **验证点**:
  - [ ] 两个求解器状态均为 SUCCESS
  - [ ] 两个求解器结果一致：||x_direct - x_iterative|| < 1e-6
  - [ ] 残差范数 < 1e-8
  - [ ] BiCGSTAB 迭代次数合理（< 500 次）
  - [ ] 不定矩阵处理稳定，无数值崩溃

### 辅助函数
- build_prism6_mixed_av_stiffness_matrix() - 构建 PRISM6_MIXED_AV 刚度矩阵
- compare_solver_results() - 对比两个求解器结果的一致性

---

## Task 9: 更新 CMakeLists.txt 构建配置
**文件**: `test/CMakeLists.txt`, `src/CMakeLists.txt`（如需要）
**优先级**: 高（编译集成必需）
**预估工作量**: 小

### 子任务
- [ ] 9.1 在 test/CMakeLists.txt 中添加3个新测试可执行目标
  - test_linear_solver_scalar（链接 numeric_lib, spdlog_lib, GTest）
  - test_linear_solver_edge（链接 numeric_lib, spdlog_lib, GTest）
  - test_linear_solver_mixed_av（链接 numeric_lib, spdlog_lib, GTest）
- [ ] 9.2 配置正确的 include_directories（Eigen路径、spdlog路径）
- [ ] 9.3 配置 link_libraries（numeric_lib 包含新的求解器源文件）
- [ ] 9.4 添加 add_test() 注册测试到 CTest 系统
- [ ] 9.5 验证 cmake && make 全流程编译成功

### 验证标准
- [ ] cmake 配置无错误
- [ ] 编译3个测试目标均成功
- [ ] ctest 运行全部通过

---

# Task Dependencies
- [Task 1] 必须最先完成（基类定义）
- [Task 2] 依赖于 [Task 1]（转换工具依赖基类的类型定义）
- [Task 3] 依赖于 [Task 1, Task 2]（直接求解器依赖基类和转换工具）
- [Task 4] 依赖于 [Task 1, Task 2]（迭代求解器依赖基类和转换工具）
- [Task 5] 依赖于 [Task 3, Task 4]（工厂类需要所有具体求解器实现完成）
- [Task 6, 7, 8] 可并行开发（测试用例相互独立，但都依赖 [Task 3, 4]）
- [Task 9] 依赖于 [Task 6, 7, 8]（CMake配置需要在测试文件就绪后添加）

## 并行执行建议
**Phase 1 (串行)**: Task 1 → Task 2
**Phase 2 (并行)**: Task 3 + Task 4 （直接求解器和迭代求解器可并行开发）
**Phase 3 (串行)**: Task 5 （工厂类整合）
**Phase 4 (并行)**: Task 6 + Task 7 + Task 8 （三个测试用例并行编写）
**Phase 5 (串行)**: Task 9 （CMake集成和最终验证）
