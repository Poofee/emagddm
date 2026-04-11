输入文件(.json/.aedt/.xml)
    │
    ▼
SolverApp::initialize()
    │
    ▼
ProjectManager::openProject()  ──→  IProjectLoader(AEDT/JSON/XML)
    │                                    │
    │                                    ▼
    │                               填充 Material/Boundary/Excitation/Mesh/EMMeshData等
    │
    ▼
SolverApp::run()
    │
    ├── 1. GlobalEdgeIDGenerator::generate()  ← EMMeshData
    │       └── 输出: elem_local_to_global_edge
    │
    ├── 2. EMDOFManager::build()  ← EMMeshData + edge_mapping
    │       └── 输出: Local2Global映射表, constrained_dof_values, num_free_dofs
    │
    ├── 3. EMAssembly::assemble()  ← EMMeshData + Local2Global + MaterialProperties
    │       └── 内部: createIntegrator() → computeAllMatrices() → scatter → COO→CSR
    │       └── 输出: K_csr, M_csr, C_csr, F_vector
    │
    ├── 4. EMLinearSolverBase::set_matrix(K) → solve(F)
    │       └── 输出: SolverResult(x, status, iterations, residual)
    │
    └── 5. FETI-DP迭代求解 (FETIDPContext)
            └── 输出: 最终解向量