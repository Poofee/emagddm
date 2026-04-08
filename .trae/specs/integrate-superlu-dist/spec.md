# SuperLU 集成 Spec（修订版）

## Why

当前项目的 SuperLU_MT 库集成在 Windows/MSVC 平台存在 pthread 兼容性问题（SEH 异常 0xc0000005），导致所有实际求解功能崩溃。经过全面测试和评估，发现：
- SuperLU_DIST 需要 MPI 环境，不适合项目的共享内存使用模式
- SuperLU_MT 的 pthread 兼容层完善需要 8-16 小时，成功率仅 40%
- Eigen 后端已集成且性能优秀，可作为默认后端

本 Spec 修订为**采用渐进式 SuperLU 集成策略**：短期使用 Eigen 后端，中期评估 SuperLU serial 版本，长期保留 SuperLU_DIST/M 选项。

## What Changes

### 已完成的修改（需要回滚）
- ❌ `CMakeLists.txt` - USE_SUPERLU_DIST 配置（不适用，需回滚）
- ❌ `src/CMakeLists.txt` - 链接 SuperLU_DIST（不适用，需回滚）
- ❌ `include/numeric/em_solver_backends.hpp` - SuperLU_DIST 头文件（需回滚到 SuperLU_MT 或 Eigen-only）
- ❌ `cmake/FindSuperLUDIST.cmake` - CMake 查找模块（保留供未来使用）

### 实际实施的修改
- ✅ **回滚到 Eigen 后端作为默认后端**
- ✅ **保留 SuperLU 宏定义但默认不启用**
- ✅ **更新文档说明 SuperLU 集成状态**
- ✅ **保留测试套件供未来使用**

## Impact

- Affected specs: 线性求解器多后端支持能力（暂时使用单 Eigen 后端）
- Affected code: 
  - `src/numeric/em_direct_solvers.cpp` - 默认使用 Eigen
  - `include/numeric/em_solver_backends.hpp` - 简化为 Eigen-only 或保留 SuperLU 桩
  - `docs/superlu_integration_status.md` - 新增状态说明文档
- External dependencies: 
  - **移除**: SuperLU_DIST（暂不使用）
  - **保留**: Eigen（已集成）
  - **未来可选**: SuperLU serial / SuperLU_MT (vcpkg) / SuperLU_DIST (vcpkg)

## ADDED Requirements

### Requirement: Eigen 后端作为默认后端
系统 SHALL 默认使用 Eigen 稀疏求解器作为对称正定矩阵的求解后端：

```cpp
// em_direct_solvers.cpp - 默认实现
SolverResult SymmetricDirectSolver::decompose() {
    // 使用 Eigen 的 SimplicialLLT 或 SparseLU
    return decompose_with_eigen();
}
```

#### Scenario: Eigen 后端求解
- **WHEN** 用户调用 `SymmetricDirectSolver::set_matrix(A)`
- **THEN** 内部调用 Eigen 的 `SimplicialLLT<SparseMatrix<double>>` 执行 Cholesky 分解
- **THEN** 缓存 L 因子
- **WHEN** 用户调用 `solve(b)`
- **THEN** 执行三角回代 `x = L^{-T} L^{-1} b`
- **THEN** 返回解向量 x，满足 ||b - Ax|| / (||A||·||x|| + ||b||) < 1e-10

### Requirement: SuperLU 后端预留接口
系统 SHALL 保留 SuperLU 后端的宏定义和接口，供未来启用：

```cpp
// em_solver_backends.hpp
#ifdef HAVE_SUPERLU
    #define EM_SOLVER_HAS_SUPERLU 1
    // SuperLU 头文件包含（暂不启用）
#else
    #define EM_SOLVER_HAS_SUPERLU 0
#endif
```

#### Scenario: 未来启用 SuperLU
- **WHEN** 用户设置 CMake 选项 `USE_SUPERLU=ON`
- **THEN** 编译时定义 `HAVE_SUPERLU` 宏
- **THEN** `EM_SOLVER_HAS_SUPERLU` 为 1
- **THEN** SuperLU 后端代码被编译
- **THEN** 用户可通过 `DirectBackendType::SUPERLU` 选择 SuperLU 后端

### Requirement: 文档说明
系统 SHALL 提供详细的 SuperLU 集成状态文档：

1. `docs/superlu_integration_status.md` - 当前状态和建议方案
2. `docs/superlu_mt_test_report.md` - SuperLU_MT 测试结果
3. README.md 中的依赖说明

## MODIFIED Requirements

### Requirement: CMake 配置选项
**修改前**:
```cmake
option(USE_SUPERLU_DIST "Enable SuperLU_DIST support" ON)
```

**修改后**:
```cmake
# SuperLU 集成状态：实验性（默认禁用，使用 Eigen 后端）
# 可选值：OFF(默认), SERIAL, MT, DIST
#   - SERIAL: SuperLU serial 版本（无 pthread 依赖）
#   - MT: SuperLU_MT 多线程版本（需要 pthread 兼容层）
#   - DIST: SuperLU_DIST 分布式版本（需要 MPI）
set(SUPERLU_BACKEND "OFF" CACHE STRING "SuperLU backend type")
set_property(CACHE SUPERLU_BACKEND PROPERTY STRINGS OFF SERIAL MT DIST)

if(SUPERLU_BACKEND STREQUAL "SERIAL")
    # SuperLU serial 配置
elseif(SUPERLU_BACKEND STREQUAL "MT")
    # SuperLU_MT 配置（实验性）
elseif(SUPERLU_BACKEND STREQUAL "DIST")
    # SuperLU_DIST 配置（需要 MPI）
endif()
```

### Requirement: 后端选择
**修改前**: 用户通过 `DirectBackendType::SUPERLU` 选择 SuperLU 后端

**修改后**: 
```cpp
// 默认使用 Eigen
SymmetricDirectSolver solver(DirectBackendType::EIGEN);

// SuperLU 后端（未来启用后使用）
// SymmetricDirectSolver solver(DirectBackendType::SUPERLU);
```

## REMOVED Requirements

### Requirement: SuperLU_DIST 预编译库集成
**Reason**: SuperLU_DIST 需要 MPI 环境，不适合项目的共享内存使用模式
**Migration**: 保留 `lib/superlu_dist` 目录供未来参考，但不主动集成

### Requirement: Windows pthread 兼容层完善
**Reason**: 工作量过大（8-16 小时），成功率低（40%），性价比不高
**Migration**: 
- 删除 `lib/superlu_mt/SRC/pthread.h`（如已创建）
- 删除 `lib/superlu_mt/SRC/superlu_mt_msvc_fix.h`（如已创建）
- 使用 Eigen 后端替代

## 技术约束

1. **默认后端**: Eigen SimplicialLLT 或 SparseLU
2. **精度要求**: 双精度 (double)
3. **编译选项**: 移除 `HAVE_SUPERLU` 定义（除非显式启用）
4. **性能目标**: 
   - 100x100 矩阵分解时间 < 50ms
   - 1000x1000 矩阵分解时间 < 2s
   - 残差范数 < 1e-10
5. **向后兼容**: 保持 DirectBackendType 枚举不变

## 实施步骤

1. **回滚 SuperLU_DIST 配置** (1 小时)
   - 恢复 `CMakeLists.txt` 到 SuperLU_MT 版本或删除 SuperLU 配置
   - 恢复 `em_solver_backends.hpp` 到原始版本

2. **简化为 Eigen-only 后端** (0.5 小时)
   - 移除 `em_direct_solvers.cpp` 中的 SuperLU 适配层（或保留为桩）
   - 确保 Eigen 后端完全功能正常

3. **更新文档** (0.5 小时)
   - 创建 `docs/superlu_integration_status.md`（已完成）
   - 更新 README.md

4. **测试验证** (1 小时)
   - 运行 `test_superlu_comprehensive`（跳过 SuperLU 相关测试）
   - 验证 Eigen 后端所有测试通过

**预计总工作量**: 3 小时

---

**版本**: 2.0 (Revised)  
**创建日期**: 2026-04-08  
**修订日期**: 2026-04-08  
**状态**: 待审核  
**建议方案**: 采用 Eigen 后端作为默认后端
