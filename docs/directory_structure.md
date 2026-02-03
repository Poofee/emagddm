# 目录与文件规范

## 目录结构（固定不变）

```plain text
FE-EM-FETIDP/  # 工程根目录（禁止修改名称）
├── CMakeLists.txt               # 根CMake（管理所有子目录、库依赖，禁止拆分）
├── config/                      # 配置文件目录（仅放JSON配置，文件名小写+下划线）
├── mesh/                        # 测试网格目录（仅放.msh格式，文件名小写+下划线）
├── include/                     # 头文件目录（按模块划分，与src目录一一对应）
│   ├── fe_em/                   # 电磁物理层头文件（.h/.hpp）
│   ├── fetidp/                  # FETI-DP核心层头文件（.h/.hpp）
│   ├── numeric/                 # 核心数值层头文件（.h/.hpp）
│   ├── tool/                    # 基础工具层头文件（.h/.hpp）
│   └── app/                     # 应用层头文件（.h/.hpp）
├── src/                         # 源文件目录（按模块划分，与include一一对应）
│   ├── fe_em/                   # 电磁物理层源文件（.cpp）
│   ├── fetidp/                  # FETI-DP核心层源文件（.cpp）
│   ├── numeric/                 # 核心数值层源文件（.cpp）
│   ├── tool/                    # 基础工具层源文件（.cpp）
│   ├── app/                     # 应用层源文件（.cpp）
│   └── main.cpp                 # 驱动层主函数（仅1个，禁止拆分）
├── test/                        # 测试程序目录（按模块命名测试文件）
├── lib/                         # 第三方库目录（按库名称分类，禁止随意修改结构）
└── output/                      # 结果输出目录（自动生成，禁止手动创建子目录）
    ├── field/                   # 场量VTK输出（自动生成）
    └── log/                     # 日志输出（自动生成）
```

## 文件命名规范

### 头文件/源文件
- 文件名小写+下划线，与类/模块名称对应
- 头文件后缀 .hpp（复杂类）/.h（简单工具）
- 源文件后缀 .cpp

**正确示例**：fe_em_mesh.hpp、fetidp_context.h、sparse_matrix.cpp  
**错误示例**：FEMesh.H、FETIDPContext.cpp、SparseMatrix.C

### 配置文件
- 文件名小写+下划线，后缀 .json
- 每个场景对应1个配置文件，包含详细参数注释

**正确示例**：motor2d_steady.json、transformer3d_trans.json

### 网格文件
- 文件名小写+下划线，后缀 .msh
- 仅存放测试用网格，按场景分类命名

**正确示例**：motor2d_simple.msh、airgap3d.msh

### 文档文件
- Markdown格式（.md）
- 文件名首字母大写+其余小写

**示例**：DevelopmentDoc.md、UserManual.md

## CMakeLists.txt 规范

### 根CMake
- 负责指定C++标准（C++17）、全局编译选项、第三方库依赖、子目录包含
- 禁止添加具体模块的编译逻辑

### 子目录CMake
- 每个子目录（include/下的模块、src/下的模块）单独编写 CMakeLists.txt
- 仅负责本模块的头文件路径、源文件编译、模块内依赖

### 库集成
- 第三方库的集成逻辑，统一放在根CMake或 numeric/ 模块的CMake中
- 禁止在高层模块（app/、fetidp/）中直接添加库依赖

### 编译输出
- 可执行文件输出到根目录的 bin/ 目录（自动生成）
- 测试程序输出到 test/bin/ 目录
- 禁止输出到其他目录