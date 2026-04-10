# SuperLU_DIST 集成任务列表

## 任务分解

- [ ] **Task 1: 准备 SuperLU_DIST 预编译库**
  - [ ] 从官网下载 SuperLU_DIST 6.x Windows 预编译包
  - [ ] 解压到 `lib/superlu_dist` 目录
  - [ ] 验证包含文件：`superlu_defs.h`、`superlu_enum_consts.h` 等
  - [ ] 验证库文件：`superludist.lib` (双精度版本)
  - [ ] 验证依赖：CBLAS 库已存在于 `lib/superlu_mt/CBLAS`

- [ ] **Task 2: 创建 CMake 查找模块**
  - [ ] 创建 `cmake/FindSuperLUDIST.cmake`
  - [ ] 实现 `find_path()` 定位头文件
  - [ ] 实现 `find_library()` 定位库文件
  - [ ] 设置 `SUPERLUDIST_FOUND`、`SUPERLUDIST_INCLUDE_DIRS`、`SUPERLUDIST_LIBRARIES`
  - [ ] 处理依赖查找 (BLAS/LAPACK)
  - [ ] 添加版本检查 (>= 6.0)

- [ ] **Task 3: 修改根 CMakeLists.txt**
  - [ ] 添加 `find_package(SuperLUDIST)` 调用
  - [ ] 设置 `HAVE_SUPERLU` 编译定义
  - [ ] 包含 SuperLU_DIST 头文件路径
  - [ ] 移除 SuperLU_MT 源码编译配置
  - [ ] 保留 CBLAS 库编译配置

- [ ] **Task 4: 修改 src/CMakeLists.txt**
  - [ ] 为 numeric_lib 链接 SuperLU_DIST 库
  - [ ] 添加 `target_link_libraries(numeric_lib PRIVATE ${SUPERLUDIST_LIBRARIES})`
  - [ ] 确保链接顺序正确 (SuperLU_DIST -> CBLAS)

- [ ] **Task 5: 修改 em_solver_backends.hpp**
  - [ ] 替换 `#include "slu_mt_ddefs.h"` 为 `#include "superlu_defs.h"`
  - [ ] 添加 `#include "superlu_enum_consts.h"`
  - [ ] 更新注释说明使用 SuperLU_DIST
  - [ ] 保持 `HAVE_SUPERLU` 宏定义逻辑不变

- [ ] **Task 6: 修改 em_direct_solvers.h**
  - [ ] 检查 SuperMatrix 等类型定义是否兼容
  - [ ] 更新成员变量注释 (从 SuperLU_MT 改为 SuperLU_DIST)
  - [ ] 确认 int_t 类型定义一致

- [ ] **Task 7: 重写 em_direct_solvers.cpp 的 SuperLU 适配层**
  - [ ] 修改 `decompose_with_superlu()`:
    - 使用 `set_default_options_dist()` 初始化选项
    - 配置 `SuperLU_dist_options_t` (ParSymbFact, ColPerm 等)
    - 调用 `pdgssv()` 执行 LU 分解 (替代 pdgstrf)
    - 缓存 L/U 因子和置换向量
  - [ ] 修改 `solve_with_superlu()`:
    - 调用 `pdgtrs()` 执行三角回代 (替代 dgstrs)
    - 从 DenseMatrix 提取解向量
  - [ ] 修改 `clear()`:
    - 调用 SuperLU_DIST 的清理函数 `Destroy_SuperMatrix_Store()`
    - 释放选项结构体资源
  - [ ] 更新日志输出 (从"SuperLU_MT"改为"SuperLU_DIST")

- [ ] **Task 8: 清理 SuperLU_MT 兼容层文件**
  - [ ] 删除 `lib/superlu_mt/SRC/pthread.h`
  - [ ] 删除 `lib/superlu_mt/SRC/superlu_mt_msvc_fix.h`
  - [ ] 删除 `lib/superlu_mt/CBLAS/dcomplex.h`
  - [ ] 删除 `lib/superlu_mt/CBLAS/scomplex.h`
  - [ ] 恢复 `lib/superlu_mt/SRC/dlamch.c` 原始内容 (撤销 min/max 重命名)
  - [ ] 恢复 `lib/superlu_mt/SRC/slu_mt_util.h` 原始内容 (撤销 FLOAT 重命名)

- [ ] **Task 9: 更新测试文件**
  - [ ] 修改 `test_superlu_comprehensive.cpp`:
    - 更新头文件包含
    - 验证 API 调用与 SuperLU_DIST 兼容
    - 添加 SuperLU_DIST 特定测试 (如 MPI 模式检测)
  - [ ] 修改 `test_superlu_diagnostic.cpp`:
    - 适配 SuperLU_DIST API
  - [ ] 确保所有测试用例编译通过

- [ ] **Task 10: 构建与调试**
  - [ ] 完全清理构建目录 (Remove-Item -Recurse -Force build)
  - [ ] 重新配置 CMake: `cmake -B build -G "Visual Studio 17 2022" -A x64 -DUSE_SUPERLU_DIST=ON`
  - [ ] 构建 SuperLU_DIST 测试目标
  - [ ] 构建 numeric_lib
  - [ ] 构建所有测试可执行文件
  - [ ] 修复编译错误 (如有)

- [ ] **Task 11: 执行全面测试**
  - [ ] 运行 `test_superlu_comprehensive` 测试套件
  - [ ] 验证所有测试通过 (预期：11/11 通过)
  - [ ] 记录测试结果 (分解时间、求解时间、残差范数)
  - [ ] 对比 Eigen 后端性能
  - [ ] 生成测试报告

- [ ] **Task 12: 文档更新**
  - [ ] 更新 `docs/superlu_mt_test_report.md` 添加解决方案章节
  - [ ] 创建 `docs/superlu_dist_integration_guide.md` 集成指南
  - [ ] 更新 `docs/third_party_libraries.md` 添加 SuperLU_DIST 说明
  - [ ] 在 `README.md` 添加 SuperLU_DIST 依赖说明

## 任务依赖关系

```
Task 1 (准备库) 
    ↓
Task 2 (CMake 模块) → Task 3 (根 CMake) → Task 4 (src CMake)
                                      ↓
Task 5 (头文件切换) → Task 6 (适配层头文件) → Task 7 (适配层实现)
                                              ↓
Task 8 (清理旧文件)                          ↓
                                              ↓
Task 9 (更新测试) ←──────────────────────────┘
    ↓
Task 10 (构建调试)
    ↓
Task 11 (执行测试)
    ↓
Task 12 (文档更新)
```

## 预计工作量

| 阶段 | 任务 | 预计时间 |
|------|------|----------|
| 准备 | Task 1-2 | 2 小时 |
| CMake 配置 | Task 3-4 | 1 小时 |
| 代码适配 | Task 5-7 | 3 小时 |
| 清理 | Task 8 | 0.5 小时 |
| 测试更新 | Task 9 | 1 小时 |
| 构建调试 | Task 10 | 2 小时 |
| 测试执行 | Task 11 | 1 小时 |
| 文档 | Task 12 | 1 小时 |
| **总计** | | **11.5 小时** |

## 验收标准

1. ✅ SuperLU_DIST 库正确链接，无编译错误
2. ✅ `HAVE_SUPERLU` 宏正确定义
3. ✅ `test_superlu_comprehensive` 所有测试通过 (11/11)
4. ✅ 无 SEH 异常 (0xc0000005)
5. ✅ 求解精度满足要求 (残差 < 1e-10)
6. ✅ 性能优于或等于 Eigen 后端
7. ✅ 文档完整，包含集成指南和 API 说明

---

**创建日期**: 2026-04-08  
**状态**: 待审核
