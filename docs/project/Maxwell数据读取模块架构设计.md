# Maxwell数据读取模块架构设计

## 1. 模块概述

Maxwell数据读取模块负责解析Maxwell Project文件（.aedt/.aedtz/.xml/.amat），将数据转换为自研数据模型，为电磁仿真提供统一的数据接口。

## 2. 设计原则

- **模块化设计**：各功能模块职责单一，便于维护和扩展
- **错误处理**：完善的异常处理机制，确保数据读取的稳定性
- **性能优化**：支持大文件流式读取，避免内存溢出
- **扩展性**：支持多种文件格式和版本兼容

## 3. 模块架构图

```
┌─────────────────────────────────────────────────────────────┐
│                    Maxwell数据读取模块                        │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │  文件解析器  │  │  数据转换器  │  │    数据验证器       │  │
│  │  Parser     │  │  Converter  │  │    Validator       │  │
│  └─────────────┘  └─────────────┘  └─────────────────────┘  │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │  格式适配器  │  │  版本管理器  │  │    缓存管理器       │  │
│  │  Adapter    │  │  VersionMgr │  │    CacheMgr        │  │
│  └─────────────┘  └─────────────┘  └─────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
```

## 4. 核心组件设计

### 4.1 文件解析器 (Parser)

**职责**：
- 识别文件格式和版本
- 解析XML/二进制文件结构
- 提取原始数据

**接口设计**：
```cpp
class IMaxwellParser {
public:
    virtual ~IMaxwellParser() = default;
    virtual bool canParse(const std::string& file_path) = 0;
    virtual MaxwellFileInfo parseFileInfo() = 0;
    virtual std::vector<MaterialData> parseMaterials() = 0;
    virtual std::vector<BoundaryData> parseBoundaries() = 0;
    virtual std::vector<ExcitationData> parseExcitations() = 0;
    virtual SolutionSetupData parseSolutionSetup() = 0;
    virtual GeometryData parseGeometry() = 0;
};
```

### 4.2 数据转换器 (Converter)

**职责**：
- 将Maxwell数据转换为自研数据模型
- 处理数据类型和单位转换
- 处理枚举值映射

**接口设计**：
```cpp
class IMaxwellConverter {
public:
    virtual ~IMaxwellConverter() = default;
    virtual std::shared_ptr<tool::Material> convertMaterial(const MaterialData& data) = 0;
    virtual std::shared_ptr<tool::Boundary> convertBoundary(const BoundaryData& data) = 0;
    virtual std::shared_ptr<tool::Excitation> convertExcitation(const ExcitationData& data) = 0;
    virtual std::shared_ptr<tool::SolutionSetup> convertSolutionSetup(const SolutionSetupData& data) = 0;
    virtual std::shared_ptr<fe_em::Geometry> convertGeometry(const GeometryData& data) = 0;
};
```

### 4.3 数据验证器 (Validator)

**职责**：
- 验证转换后的数据完整性
- 检查数据一致性
- 生成验证报告

**接口设计**：
```cpp
class IMaxwellValidator {
public:
    virtual ~IMaxwellValidator() = default;
    virtual ValidationResult validateMaterials(const std::vector<std::shared_ptr<tool::Material>>& materials) = 0;
    virtual ValidationResult validateBoundaries(const std::vector<std::shared_ptr<tool::Boundary>>& boundaries) = 0;
    virtual ValidationResult validateExcitations(const std::vector<std::shared_ptr<tool::Excitation>>& excitations) = 0;
    virtual ValidationResult validateSolutionSetup(const std::shared_ptr<tool::SolutionSetup>& setup) = 0;
    virtual ValidationResult validateGeometry(const std::shared_ptr<fe_em::Geometry>& geometry) = 0;
};
```

## 5. 文件格式支持

### 5.1 支持的文件格式

| 格式 | 扩展名 | 解析方式 | 特点 |
|------|--------|----------|------|
| Maxwell项目 | .aedt | XML解析 | 完整项目配置 |
| 压缩项目 | .aedtz | 解压+XML解析 | 压缩格式 |
| XML导出 | .xml | XML解析 | 标准化格式 |
| 材料库 | .amat | XML解析 | 材料数据 |

### 5.2 版本兼容性

| Maxwell版本 | 支持状态 | 备注 |
|-------------|----------|------|
| R15-R19 | 基础支持 | 核心功能支持 |
| R20-R24 | 完全支持 | 最新功能支持 |
| 未来版本 | 预留接口 | 扩展支持 |

## 6. 错误处理机制

### 6.1 异常类型

```cpp
enum class MaxwellErrorCode {
    FILE_NOT_FOUND,
    INVALID_FORMAT,
    VERSION_NOT_SUPPORTED,
    DATA_CORRUPTED,
    CONVERSION_FAILED,
    VALIDATION_FAILED
};

class MaxwellException : public std::exception {
public:
    MaxwellException(MaxwellErrorCode code, const std::string& message);
    MaxwellErrorCode getErrorCode() const;
    const char* what() const noexcept override;
};
```

### 6.2 错误恢复策略

- **文件级错误**：跳过损坏文件，继续处理其他文件
- **数据级错误**：使用默认值或跳过错误数据
- **转换错误**：记录错误日志，提供修复建议

## 7. 性能优化策略

### 7.1 内存管理

- 流式读取大文件，避免一次性加载
- 使用智能指针管理资源
- 实现对象池减少内存分配

### 7.2 缓存机制

- 文件解析结果缓存
- 数据转换结果缓存
- 支持缓存失效策略

## 8. 扩展性设计

### 8.1 插件架构

支持通过插件方式扩展新的文件格式和解析器：

```cpp
class IMaxwellPlugin {
public:
    virtual ~IMaxwellPlugin() = default;
    virtual std::string getPluginName() = 0;
    virtual std::vector<std::string> getSupportedFormats() = 0;
    virtual std::unique_ptr<IMaxwellParser> createParser() = 0;
};
```

### 8.2 配置驱动

支持通过配置文件调整解析参数和转换规则：

```json
{
    "maxwell_reader": {
        "parser": {
            "buffer_size": 8192,
            "encoding": "UTF-8"
        },
        "converter": {
            "unit_conversion": true,
            "strict_mode": false
        },
        "validator": {
            "enable_validation": true,
            "report_level": "WARNING"
        }
    }
}
```

## 9. 测试策略

### 9.1 单元测试

- 文件解析器测试
- 数据转换器测试
- 验证器测试

### 9.2 集成测试

- 完整数据流测试
- 性能测试
- 兼容性测试

### 9.3 回归测试

- 版本兼容性测试
- 边界条件测试
- 错误处理测试

## 10. 部署与集成

### 10.1 依赖管理

- 第三方库：pugixml（XML解析）、zlib（压缩支持）
- 内部依赖：项目数据模型、工具库

### 10.2 集成方式

- 静态库集成
- 动态库集成
- 头文件包含

---

**设计完成时间**：2026-02-06
**设计者**：AI Developer
**版本**：1.0