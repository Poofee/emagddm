---
name: "code-reviewer"
description: "Reviews industrial C++ electromagnetic solver code for quality, identifying good vs bad code patterns. Invoke when user requests code review, needs to evaluate code quality, or wants to improve existing C++ electromagnetic FEM code."
---

# 工业级C++电磁求解器代码审查Skill

## Skill概述

本Skill专门针对**低频电磁有限元求解器**这类高性能、长生命周期的工业软件进行代码审查，通过系统化的检查项快速区分"好代码"与"坏代码"，确保代码达到工业级标准。

## 触发条件

当以下情况发生时调用此Skill：
- 用户请求代码审查或代码质量评估
- 需要判断代码是否符合工业级电磁求解器的编码标准
- 需要改进现有C++电磁FEM代码的质量
- 代码提交前需要质量把关
- 发现潜在的性能问题或架构缺陷

## 审查维度（按优先级排序）

### 一、命名与可读性审查（第一印象）

#### 检查要点
1. **命名是否"见名知意"**
   - ✅ 禁止无意义变量名：`a`, `tmp`, `data`
   - ✅ 禁止过度缩写：`el` 应写为 `element`（`dof` 除外）
   - ✅ 变量名应清晰表达其用途和含义

2. **命名风格是否统一**
   - 类名：大驼峰（`DOFManager`, `MatrixAssembler`）
   - 函数名：小驼峰（`assembleMatrix`, `calculateCurrentDensity`）
   - 变量名：下划线风格（`num_dofs`）或小驼峰（`numDofs`），**二选一且全程一致**

3. **是否避免"魔法数"**
   - 所有常量必须有名字
   - 示例：`3.14159` → `constexpr double PI = 3.141592653589793;`

#### 对比示例

**❌ 坏代码**：
```cpp
void calc(std::vector<double>& a, const std::vector<int>& b, Mesh* m) {
    for (int i = 0; i < m->elem.size(); i++) {
        auto e = m->elem[i];
        for (int j = 0; j < 4; j++) { // 魔法数 4
            int id = e->node[j]->id;
            a[id] += b[i] * 0.5; // 魔法数 0.5
        }
    }
}
```

**✅ 好代码**：
```cpp
constexpr int TETRAHEDRAL_NODES = 4;
constexpr double GAUSS_WEIGHT = 0.5;

void assembleRightHandSide(
    std::vector<double>& rhs_vector,
    const std::vector<double>& element_source_terms,
    const Mesh& mesh
) {
    for (const auto& element : mesh.elements) {
        const int element_id = element->getId();
        const double source_term = element_source_terms[element_id];

        for (int node_idx = 0; node_idx < TETRAHEDRAL_NODES; ++node_idx) {
            const Node& node = element->getNode(node_idx);
            const int dof_id = node.getDofId();
            rhs_vector[dof_id] += source_term * GAUSS_WEIGHT;
        }
    }
}
```

### 二、函数与类设计审查（核心架构）

#### 检查要点
1. **函数是否"单一职责"**
   - 一个函数只做一件事
   - 长度原则上不超过 **50行**（性能热点代码除外，但需加注释说明）

2. **类是否"高内聚、低耦合"**
   - ❌ 禁止"上帝类"（一个类包含所有功能）
   - ✅ 使用依赖注入而非硬编码（不直接依赖 `MeshSingleton`，而是通过构造函数传 `std::shared_ptr<Mesh>`）

3. **参数数量控制**
   - 函数参数超过 **5个** 时，考虑用结构体封装

#### 对比示例（绕组建模场景）

**❌ 坏代码**：
```cpp
class WindingManager {
public:
    WindingManager() {
        mesh_ = MeshSingleton::getInstance(); // 硬编码依赖
    }

    void doEverything() { // 一个函数做了N件事
        initWinding();
        calcCurrent();
        assembleMatrix();
        calcInductance();
    }

private:
    Mesh* mesh_;
    std::vector<double> current_;
    std::vector<double> inductance_;
    // ... 几十上百个成员变量
};
```

**✅ 好代码**：
```cpp
struct WindingConfig {
    int num_turns = 100;
    double wire_diameter = 0.001;
    int parallel_branches = 1;
};

class WindingCurrentCalculator {
public:
    explicit WindingCurrentCalculator(std::shared_ptr<Mesh> mesh, WindingConfig config)
        : mesh_(std::move(mesh)), config_(std::move(config)) {}

    std::vector<double> calculateCurrentDensity() const;

private:
    std::shared_ptr<Mesh> mesh_;
    WindingConfig config_;
};

class WindingInductanceCalculator {
public:
    explicit WindingInductanceCalculator(std::shared_ptr<Mesh> mesh);
    double calculateSelfInductance(const std::vector<double>& current_density) const;
};
```

### 三、内存与资源管理审查（工业级底线）

#### 检查要点
1. **智能指针使用**
   - ❌ 禁止裸 `new`/`delete`（除非在极底层的性能热点，且封装在类内）
   - ✅ 独占所有权用 `std::unique_ptr`
   - ✅ 共享所有权用 `std::shared_ptr`

2. **避免不必要的拷贝**
   - 大对象（矩阵、网格）传递用 `const&` 或 `&&`（移动语义）
   - 循环内避免重复构造大对象

3. **内存泄漏风险检查**
   - 是否有未关闭的文件句柄
   - 是否有未释放的资源（智能指针能解决大部分问题）

#### 对比示例（矩阵传递场景）

**❌ 坏代码**：
```cpp
class MatrixAssembler {
public:
    void assemble(SparseMatrix mat) { // 大对象值传递，性能灾难！
        for (int i = 0; i < mat.rows(); i++) {
            // ...
        }
    }

private:
    SparseMatrix* mat_; // 裸指针，谁来 delete？
};

SparseMatrix m(10000, 10000);
assembler.assemble(m); // 拷贝了 10000x10000 的矩阵！
```

**✅ 好代码**：
```cpp
class MatrixAssembler {
public:
    void assemble(const SparseMatrix& mat) { // const引用，不拷贝
        for (int i = 0; i < mat.rows(); i++) {
            // ...
        }
    }

    void setMatrix(std::unique_ptr<SparseMatrix> mat) {
        mat_ = std::move(mat);
    }

private:
    std::unique_ptr<SparseMatrix> mat_; // 智能指针管理
};

// 使用时
SparseMatrix m(10000, 10000);
assembler.assemble(m); // 不拷贝，只读引用

auto mat_ptr = std::make_unique<SparseMatrix>(10000, 10000);
assembler.setMatrix(std::move(mat_ptr)); // 移动语义
```

### 四、错误处理与鲁棒性审查（数值稳定性）

#### 检查要点
1. **返回值检查**
   - 所有可能失败的函数必须检查返回值（文件读取、求解器收敛等）

2. **数值稳定性保护**
   - 除零保护：`if (std::abs(det) < 1e-12) throw std::runtime_error("Singular matrix");`
   - 迭代求解必须有收敛准则（不能硬编码迭代次数就结束）

3. **异常安全保证**
   - 至少保证"基本异常安全"（发生异常时不泄漏资源）

#### 对比示例（线性求解器场景）

**❌ 坏代码**：
```cpp
void solveLinearSystem(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b) {
    double det = A.determinant(); // 不检查 det 是否为 0
    MUMPS_solve(A, x, b); // 不检查求解是否收敛
    // 直接返回，用户不知道解对不对
}
```

**✅ 好代码**：
```cpp
struct SolverResult {
    bool converged = false;
    int iterations = 0;
    double residual_norm = 0.0;
};

SolverResult solveLinearSystem(
    const SparseMatrix& A,
    std::vector<double>& x,
    const std::vector<double>& b,
    const double tolerance = 1e-6,
    const int max_iterations = 1000
) {
    SolverResult result;

    const double det = A.determinant();
    if (std::abs(det) < 1e-12) {
        throw std::runtime_error("Singular matrix detected");
    }

    result = MUMPS_solve(A, x, b, tolerance, max_iterations);

    if (!result.converged) {
        throw std::runtime_error(
            "Solver did not converge after " + std::to_string(result.iterations) +
            " iterations. Residual norm: " + std::to_string(result.residual_norm)
        );
    }

    return result;
}
```

### 五、可扩展性与架构审查（面向未来）

#### 检查要点
1. **面向接口编程**
   - 用抽象基类定义接口（如 `PhysicsSolver`）
   - 具体实现（`ElectrostaticSolver`, `MagnetostaticSolver`）通过多态扩展

2. **开闭原则**
   - 加新功能时只需"加新类"，无需"改旧代码"

3. **消除硬编码 switch-case**
   - 大量 `switch-case` 是扩展性差的信号
   - 用多态或策略模式替代

#### 对比示例（物理场扩展场景）

**❌ 坏代码**：
```cpp
enum class PhysicsType { Electrostatic, Magnetostatic };

class Solver {
public:
    void solve(PhysicsType type) {
        switch (type) {
            case PhysicsType::Electrostatic:
                solveElectrostatic();
                break;
            case PhysicsType::Magnetostatic:
                solveMagnetostatic();
                break;
            // 加 EddyCurrent 要改这里，违反开闭原则
        }
    }
};
```

**✅ 好代码**：
```cpp
class PhysicsSolver {
public:
    virtual ~PhysicsSolver() = default;
    virtual void initialize() = 0;
    virtual void assemble() = 0;
    virtual SolverResult solve() = 0;
};

class ElectrostaticSolver : public PhysicsSolver {
public:
    void initialize() override;
    void assemble() override;
    SolverResult solve() override;
};

class MagnetostaticSolver : public PhysicsSolver {
public:
    void initialize() override;
    void assemble() override;
    SolverResult solve() override;
};

class SolverScheduler {
public:
    void registerSolver(std::unique_ptr<PhysicsSolver> solver) {
        solvers_.push_back(std::move(solver));
    }

    void runAll() {
        for (auto& solver : solvers_) {
            solver->initialize();
            solver->assemble();
            const auto result = solver->solve();
        }
    }

private:
    std::vector<std::unique_ptr<PhysicsSolver>> solvers_;
};
```

## 审查执行流程

### 步骤1：整体架构评估（优先级：⭐⭐⭐⭐⭐）
- [ ] 类图是否清晰
- [ ] 模块划分是否合理
- [ ] 是否存在"上帝类"
- [ ] 依赖关系是否合理

### 步骤2：核心模块深度审查（优先级：⭐⭐⭐⭐）
重点检查以下核心类/函数：
- `DOFManager` - 自由度管理
- `MatrixAssembler` - 矩阵装配
- `SolverScheduler` - 求解调度
- `Mesh` / `Element` / `Node` - 网格数据结构

### 步骤3：代码细节审查（优先级：⭐⭐⭐）
- [ ] 命名规范一致性
- [ ] 内存管理安全性
- [ ] 错误处理完备性
- [ ] 注释质量

### 步骤4：测试与文档检查（优先级：⭐⭐）
- [ ] 单元测试覆盖率
- [ ] 关键算法的文档说明
- [ ] API文档完整性

## 审查输出模板

### 评分标准

| 维度 | 权重 | 评分（1-10） | 说明 |
|------|------|-------------|------|
| 命名与可读性 | 15% | | |
| 函数与类设计 | 25% | | |
| 内存与资源管理 | 20% | | |
| 错误处理与鲁棒性 | 20% | | |
| 可扩展性与架构 | 20% | | |

### 输出格式示例

```markdown
# 代码审查报告

## 总体评价
[优秀/良好/需改进/不合格] - [一句话总结]

## 详细发现

### ✅ 优点
1. [具体优点描述]

### ⚠️ 需改进项
1. **[问题描述]**
   - 位置：[文件名:行号]
   - 问题类型：[命名/设计/内存/错误处理/扩展性]
   - 严重程度：🔴严重 / 🟡中等 / 🟢轻微
   - 当前代码：
     ```cpp
     [坏代码示例]
     ```
   - 建议修改为：
     ```cpp
     [好代码示例]
     ```

## 评分汇总
- 命名与可读性：X/10
- 函数与类设计：X/10
- 内存与资源管理：X/10
- 错误处理与鲁棒性：X/10
- 可扩展性与架构：X/10
- **加权总分：X.X/10**

## 改进建议优先级
1. 🔴 [立即修复] - ...
2. 🟡 [尽快修复] - ...
3. 🟢 [建议优化] - ...
```

## 特殊场景检查清单

### 绕组建模相关
- [ ] 电流密度计算是否正确处理匝数
- [ ] 电感矩阵计算是否有对称性验证
- [ ] 多支路绕组是否支持并行计算

### 矩阵装配相关
- [ ] 稀疏矩阵存储格式选择是否合理
- [ ] 边界条件处理是否正确
- [ ] 积分精度是否满足要求

### 求解器相关
- [ ] 收敛准则是否可配置
- [ ] 预条件子选择是否合适
- [ ] 并行效率是否经过测试

## 注意事项

1. **上下文理解**：审查前先了解代码的业务背景和设计意图
2. **平衡原则**：工业软件需要在性能、可维护性、扩展性之间取得平衡
3. **建设性反馈**：不仅指出问题，还要提供具体的改进方案
4. **尊重历史**：对于遗留代码，建议渐进式改进而非一次性重写
5. **关注安全**：电磁求解器涉及工程计算，数值稳定性至关重要
