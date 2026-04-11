# 工业级C++电磁求解器代码审查清单
这份清单专门针对**低频电磁有限元求解器**这类高性能、长生命周期的工业软件设计，通过可落地的检查项快速区分"好代码"与"坏代码"。

---

## 一、命名与可读性审查（第一印象）
### 检查要点
1. **命名是否"见名知意"**：禁止无意义变量名（`a`、`tmp`、`data`），禁止过度缩写（`el` 应写为 `element`，`dof` 除外）。
2. **命名风格是否统一**：类名用大驼峰（`DOFManager`），函数名用小驼峰（`assembleMatrix`），变量名用下划线（`num_dofs`）或小驼峰（`numDofs`），二选一且全程一致。
3. **是否避免"魔法数"**：所有常量必须有名字（`3.14159` → `const double PI = 3.141592653589793;`）。

### 对比示例
**坏代码**：
```cpp
// 坏代码：命名模糊，魔法数
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

**好代码**：
```cpp
// 好代码：命名清晰，常量有名字
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

---

## 二、函数与类设计审查（核心架构）
### 检查要点
1. **函数是否"单一职责"**：一个函数只做一件事，长度原则上不超过 50 行（性能热点代码除外，但需加注释）。
2. **类是否"高内聚、低耦合"**：
   - 禁止"上帝类"（一个类包含所有功能）。
   - 依赖注入而非硬编码（不直接依赖 `MeshSingleton`，而是通过构造函数传 `std::shared_ptr<Mesh>`）。
3. **参数是否过多**：函数参数超过 5 个时，考虑用结构体封装。

### 对比示例（绕组建模场景）
**坏代码**：
```cpp
// 坏代码：上帝类，职责混乱，硬编码依赖
class WindingManager {
public:
    WindingManager() {
        mesh_ = MeshSingleton::getInstance(); // 硬编码依赖
    }
    
    // 一个函数做了N件事：初始化、计算电流、装配矩阵、计算电感
    void doEverything() {
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

**好代码**：
```cpp
// 好代码：职责分离，依赖注入
struct WindingConfig {
    int num_turns = 100;
    double wire_diameter = 0.001;
    int parallel_branches = 1;
};

class WindingCurrentCalculator {
public:
    explicit WindingCurrentCalculator(std::shared_ptr<Mesh> mesh, WindingConfig config)
        : mesh_(std::move(mesh)), config_(std::move(config)) {}
    
    std::vector<double> calculateCurrentDensity() const; // 只做一件事
    
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

---

## 三、内存与资源管理审查（工业级底线）
### 检查要点
1. **是否用智能指针管理资源**：
   - 禁止裸 `new`/`delete`（除非在极底层的性能热点，且封装在类内）。
   - 独占所有权用 `std::unique_ptr`，共享所有权用 `std::shared_ptr`。
2. **是否避免不必要的拷贝**：
   - 大对象（矩阵、网格）传递用 `const&` 或 `&&`（移动语义）。
   - 循环内避免重复构造大对象。
3. **是否有内存泄漏风险**：检查是否有未关闭的文件、未释放的句柄（虽然智能指针能解决大部分问题）。

### 对比示例（矩阵传递场景）
**坏代码**：
```cpp
// 坏代码：裸指针，不必要的拷贝
class MatrixAssembler {
public:
    // 矩阵值传递，发生拷贝！
    void assemble(SparseMatrix mat) { // 大对象拷贝，性能灾难
        for (int i = 0; i < mat.rows(); i++) {
            // ...
        }
    }
    
private:
    SparseMatrix* mat_; // 裸指针，谁来 delete？
};

// 使用时
SparseMatrix m(10000, 10000);
assembler.assemble(m); // 拷贝了 10000x10000 的矩阵！
```

**好代码**：
```cpp
// 好代码：智能指针，移动语义，引用传递
class MatrixAssembler {
public:
    // 矩阵用 const& 传递，避免拷贝
    void assemble(const SparseMatrix& mat) {
        for (int i = 0; i < mat.rows(); i++) {
            // ...
        }
    }
    
    // 或者用移动语义，转移所有权
    void setMatrix(std::unique_ptr<SparseMatrix> mat) {
        mat_ = std::move(mat);
    }
    
private:
    std::unique_ptr<SparseMatrix> mat_; // 智能指针，自动管理内存
};

// 使用时
SparseMatrix m(10000, 10000);
assembler.assemble(m); // 不拷贝，只读引用

// 或者转移所有权
auto mat_ptr = std::make_unique<SparseMatrix>(10000, 10000);
assembler.setMatrix(std::move(mat_ptr)); // 移动，不拷贝
```

---

## 四、错误处理与鲁棒性审查（数值稳定性）
### 检查要点
1. **是否忽略返回值**：所有可能失败的函数（文件读取、求解器收敛）必须检查返回值。
2. **是否有数值稳定性保护**：
   - 检查除零（如 `if (std::abs(det) < 1e-12) throw std::runtime_error("Singular matrix");`）。
   - 迭代求解必须有收敛准则（不能硬编码迭代 100 次就结束）。
3. **是否有异常安全保证**：至少保证"基本异常安全"（发生异常时不泄漏资源）。

### 对比示例（线性求解器场景）
**坏代码**：
```cpp
// 坏代码：忽略返回值，无收敛检查，无除零保护
void solveLinearSystem(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b) {
    double det = A.determinant(); // 不检查 det 是否为 0
    MUMPS_solve(A, x, b); // 不检查求解是否收敛
    // 直接返回，用户不知道解对不对
}
```

**好代码**：
```cpp
// 好代码：检查返回值，有收敛准则，有数值保护
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
    
    // 检查矩阵奇异性
    const double det = A.determinant();
    if (std::abs(det) < 1e-12) {
        throw std::runtime_error("Singular matrix detected, determinant is nearly zero.");
    }
    
    // 调用求解器并检查收敛
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

---

## 五、可扩展性与架构审查（面向未来）
### 检查要点
1. **是否面向接口编程**：用抽象基类定义接口（如 `PhysicsSolver`），具体实现（`ElectrostaticSolver`、`MagnetostaticSolver`）通过多态扩展。
2. **是否遵循"开闭原则"**：加新功能时是否只需"加新类"，而无需"改旧代码"。
3. **是否有硬编码的 `switch-case`**：大量 `switch-case` 通常是扩展性差的信号，考虑用多态或策略模式替代。

### 对比示例（物理场扩展场景）
**坏代码**：
```cpp
// 坏代码：硬编码 switch-case，加新物理场要改旧代码
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
            // 加 EddyCurrent 要在这里加 case，违反开闭原则
        }
    }
};
```

**好代码**：
```cpp
// 好代码：面向接口，插件化，加新物理场只需加新类
class PhysicsSolver {
public:
    virtual ~PhysicsSolver() = default;
    virtual void initialize() = 0;
    virtual void assemble() = 0;
    virtual SolverResult solve() = 0;
};

// 静电场求解器（独立模块）
class ElectrostaticSolver : public PhysicsSolver {
public:
    void initialize() override;
    void assemble() override;
    SolverResult solve() override;
};

// 静磁场求解器（独立模块）
class MagnetostaticSolver : public PhysicsSolver {
public:
    void initialize() override;
    void assemble() override;
    SolverResult solve() override;
};

// 核心调度器（不需要修改，只需注册新 solver）
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
            // 处理结果
        }
    }

private:
    std::vector<std::unique_ptr<PhysicsSolver>> solvers_;
};
```

---

## 六、代码审查执行流程（按优先级）
1. **先看整体架构**：类图是否清晰，模块划分是否合理，有没有"上帝类"。
2. **再看核心类/函数**：`DOFManager`、`MatrixAssembler`、`SolverScheduler` 这些核心模块的设计是否符合上述原则。
3. **再看具体代码细节**：命名、内存管理、错误处理。
4. **最后看测试与文档**：有没有单元测试，关键函数有没有注释说明"做什么"以及"为什么这么做"（不用说明"怎么做"）。
