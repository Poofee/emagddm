# 命名规范

所有命名需清晰、规范，见名知意，禁止使用拼音、简写（通用简写除外，如DOF、MPI、VTK），禁止大小写混用。

## 类命名
- **规则**：PascalCase（首字母大写，每个单词首字母大写），后缀可添加模块标识
- **示例**：EMMesh、APhiAssembler、FETIDPContext、SparseMatrix、DirectSolver
- **禁止**：em_mesh、aPhiAssembler、FETIDP_Context、sparseMatrix

## 函数/方法命名
- **规则**：camelCase（首字母小写，后续单词首字母大写），动词开头，明确函数功能
- **示例**：solveSubdomain、buildSchurComplement、extractInterfaceDOF、calcMagneticFlux
- **禁止**：subdomain_solve、schurBuild、getDOF、calcB

## 变量命名
- **规则**：snake_case（全小写，单词之间用下划线连接），见名知意
- **示例**：subdomain_id、convergence_threshold、schur_matrix、magnetic_flux_density
- **禁止**：subID、convThresh、S、B（注释说明除外）

## 常量/宏命名
- **规则**：全大写，单词之间用下划线连接，放在头文件顶部，集中管理
- **示例**：MAX_ITERATION、DEFAULT_CONV_THRESH、PI、CSR_FORMAT
- **禁止**：maxIter、defaultConv、Pi、csrFormat

## 结构体/枚举命名
- **结构体**：PascalCase，后缀可添加Data（如SubdomainData、CoarseGridData）
- **枚举**：PascalCase，枚举值全大写+下划线

```cpp
enum SolverType { 
    DIRECT_SOLVER, 
    ITERATIVE_SOLVER 
};
```

## 命名空间命名
- **规则**：全小写，按模块划分，避免命名冲突，禁止使用全局命名空间
- **示例**：namespace fe_em、namespace fetidp、namespace numeric、namespace tool