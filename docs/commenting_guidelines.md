# 注释规范

注释需清晰、简洁、准确，覆盖"是什么、做什么、为什么"，禁止无用注释、冗余注释，禁止注释与代码不一致，支持中文注释（编码UTF-8）。

## 头文件注释（文件顶部）

每个头文件顶部添加文件说明注释，包含：文件名称、模块归属、功能描述、作者、创建日期、修改记录。

```cpp
/**
 * @file fe_em_mesh.hpp
 * @brief 电磁物理层 - 电磁网格管理模块头文件
 * @details 基于Gmsh API实现2D/3D网格的读取、解析、子域拆分、界面DOF提取，适配低频电磁A-φ法的网格特性
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */
```

## 类注释

每个类定义前添加类注释，说明类的功能、核心职责、适用场景，若有继承/依赖，需说明继承关系、依赖模块。

```cpp
/**
 * @class APhiAssembler
 * @brief A-φ有限元组装类，负责低频电磁A-φ方程的单元级组装
 * @details 支持2D/3D、稳态/瞬态场景的刚度矩阵（K）、涡流阻尼矩阵（C）、载荷向量（F）组装，
 *          实现介质区纯A组装、导体区A-φ解耦组装，为FETI-DP子域求解提供子域矩阵/向量
 * @note 依赖 EMMesh、Material、LoadBoundary 类，属于电磁物理层核心模块
 */
class APhiAssembler { ... };
```

## 函数/方法注释

每个函数（尤其是公有接口）定义前添加函数注释，说明函数功能、参数含义、返回值、异常情况、使用注意事项。

```cpp
/**
 * @brief 求解单个子域的A-φ方程，得到子域解x_i=[A_i;φ_i]
 * @param subdomain_id 子域ID（需在[0, subdomain_num-1]范围内）
 * @param lambda 全局对偶变量（拉格朗日乘子），用于施加界面约束
 * @param x_i 输出参数，子域解向量（A和φ的自由度组合）
 * @return bool 求解成功返回true，失败返回false（失败原因通过日志输出）
 * @note 需提前初始化子域求解器并缓存分解结果，否则求解效率极低
 * @exception 若subdomain_id无效、lambda维度不匹配，输出ERROR日志并返回false
 */
bool solveSubdomain(int subdomain_id, const DenseVector& lambda, DenseVector& x_i);
```

## 代码内注释

### 复杂逻辑
每步关键操作添加注释，说明操作目的（如Schur补构造、粗矫正步骤）。

```cpp
// 构造约束矩阵B_i（0-1稀疏矩阵），描述子域界面DOF与全局对偶变量的关联
// B_i的维度：界面DOF数 × 全局对偶变量数，界面DOF对应位置为1，其余为0
SparseMatrix B_i = constructConstraintMatrix(interface_dof, global_dual_dof_num);

// 计算粗网格刚度矩阵Kc = P^T * S * P（P为粗投影算子，S为全局Schur补矩阵）
// Kc维度：粗节点数 × 粗节点数，仅需计算一次并缓存
SparseMatrix Kc = P.transpose() * S * P;
```

### 关键参数
定义关键参数时，添加注释说明参数含义、取值范围。

### 算法细节
FETI-DP核心算法（如约束矩阵B_i构造、粗投影算子P构造），需添加注释说明算法逻辑、公式对应关系。

### 禁止行为
- 注释与代码重复（如`// 定义变量i，int i = 0;`）
- 无用注释（如`// 开始循环`）