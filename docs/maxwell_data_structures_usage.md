# Maxwell AEDT 数据结构使用指南

## 概述

本文档介绍如何使用新设计的Maxwell AEDT文件数据结构进行数据读取、转换和存储。该数据结构提供了完整的C++17兼容的接口，用于处理Maxwell电磁仿真项目的配置信息。

## 核心数据结构

### 1. 项目数据容器 (MaxwellProjectData)

`MaxwellProjectData` 是顶层容器，包含完整的项目信息：

```cpp
#include "tool/maxwell_data_structures.hpp"
#include "tool/maxwell_data_converter.hpp"

using namespace fe_em::tool::maxwell_data;

// 创建项目数据容器
MaxwellProjectData project_data;

// 访问项目元数据
const auto& metadata = project_data.metadata;
std::cout << "项目创建时间: " << metadata.created_time << std::endl;
std::cout << "产品类型: " << metadata.product << std::endl;

// 访问材料库
for (const auto& [name, material] : project_data.materials) {
    std::cout << "材料: " << name << std::endl;
    if (material.properties.permittivity) {
        std::cout << "介电常数: " << material.properties.permittivity->value << std::endl;
    }
}

// 访问设计模型
for (const auto& design_model : project_data.design_models) {
    std::cout << "设计模型: " << design_model.name << std::endl;
    std::cout << "求解类型: " << static_cast<int>(design_model.solution_type) << std::endl;
}
```

### 2. 材料属性 (Material)

材料数据结构支持电磁、热学、结构等多种物理属性：

```cpp
// 创建新材料
Material copper;
copper.name = "copper";
copper.coord_system = CoordinateSystemType::CARTESIAN;
copper.is_bulk_material = true;
copper.physics_types = {PhysicsType::ELECTROMAGNETIC, PhysicsType::THERMAL};

// 设置材料属性
copper.properties.conductivity = PhysicalQuantity(5.8e7, UnitType::SIEMENS_PER_METER, "58000000");
copper.properties.thermal_conductivity = PhysicalQuantity(400, UnitType::WATT, "400");

// 设置外观属性
copper.properties.appearance = MaterialAppearance(242, 140, 102);
```

### 3. 几何模型 (GeometryPart)

几何模型支持多种操作类型和拓扑信息：

```cpp
// 创建几何部件
GeometryPart box_part;
box_part.attributes.name = "Box1";
box_part.attributes.material_name = "copper";
box_part.attributes.color = {143, 175, 143};

// 创建几何操作
GeometryOperation box_operation;
box_operation.type = GeometryOperationType::BOX;
box_operation.id = 5;

// 设置长方体参数
BoxParameters box_params;
box_params.position = {-0.6, -2.2, 0};  // mm
box_params.size = {1.4, 1.0, 1.8};      // mm
box_operation.parameters = box_params;

// 添加到部件
box_part.operations.push_back(box_operation);
```

### 4. 边界条件 (BoundaryCondition)

边界条件支持多种类型和激励源：

```cpp
// 创建电流激励边界
BoundaryCondition current_boundary;
current_boundary.id = 7;
current_boundary.type = BoundaryType::CURRENT;
current_boundary.current = PhysicalQuantity(0.01, UnitType::AMPERE, "10mA");
current_boundary.faces = {92};  // 关联的面ID

// 创建磁场边界
BoundaryCondition h_field_boundary;
h_field_boundary.id = 5;
h_field_boundary.type = BoundaryType::TANGENTIAL_H_FIELD;
h_field_boundary.faces = {10};
```

## 数据转换工具使用

### 1. 从文件直接转换

最简单的使用方式是从AEDT文件直接转换：

```cpp
#include "tool/maxwell_data_converter.hpp"

using namespace fe_em::tool::maxwell_converter;

// 从文件转换
std::string file_path = "docs/project/Project49.aedt";
auto project_data = MaxwellDataConverter::convert_from_file(file_path);

if (project_data) {
    std::cout << "成功转换项目数据" << std::endl;
    
    // 访问转换后的数据
    const auto& materials = project_data->materials;
    const auto& designs = project_data->design_models;
    
    std::cout << "材料数量: " << materials.size() << std::endl;
    std::cout << "设计模型数量: " << designs.size() << std::endl;
} else {
    std::cout << "转换失败" << std::endl;
}
```

### 2. 使用解析器手动转换

如果需要更精细的控制，可以使用解析器手动转换：

```cpp
#include "tool/maxwell_parser.hpp"
#include "tool/maxwell_data_converter.hpp"

// 创建解析器
auto parser = std::make_shared<maxwell_parser::MaxwellParser>();
parser->parse_file("docs/project/Project49.aedt");

// 创建转换器
MaxwellDataConverter converter(parser);

// 转换数据
auto project_data = converter.convert();

// 处理转换后的数据
// ...
```

### 3. 部分数据转换

如果只需要转换特定部分的数据，可以使用专门的转换器：

```cpp
// 只转换材料数据
auto materials_node = parser->get_root_node()->find_children("Materials")[0];
for (const auto& material_node : materials_node->children) {
    MaterialConverter material_converter(material_node);
    auto material = material_converter.convert();
    
    // 处理单个材料
    std::cout << "材料名称: " << material.name << std::endl;
}
```

## 与FETI-DP求解器的集成

### 1. 材料参数提取

将Maxwell材料属性转换为FETI-DP求解器所需的参数：

```cpp
// 提取材料电磁参数
std::unordered_map<std::string, Eigen::Matrix3d> material_tensors;

for (const auto& [name, material] : project_data.materials) {
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Identity();
    
    if (material.properties.permeability) {
        tensor(0,0) = material.properties.permeability->value;
        tensor(1,1) = material.properties.permeability->value;
        tensor(2,2) = material.properties.permeability->value;
    }
    
    material_tensors[name] = tensor;
}
```

### 2. 几何分区映射

根据几何模型生成子域划分：

```cpp
// 根据几何部件生成子域划分
std::vector<GeometryPart> subdomains;

for (const auto& design_model : project_data.design_models) {
    for (const auto& geometry_part : design_model.geometry_parts) {
        // 根据材料属性判断是否为磁性材料
        if (geometry_part.attributes.material_name.find("TDK") != std::string::npos) {
            subdomains.push_back(geometry_part);
        }
    }
}
```

### 3. 边界条件转换

将Maxwell边界条件转换为FETI-DP的界面条件：

```cpp
// 转换边界条件
std::vector<BoundaryCondition> interface_conditions;

for (const auto& design_model : project_data.design_models) {
    for (const auto& boundary : design_model.boundaries) {
        if (boundary.type == BoundaryType::SYMMETRY || 
            boundary.type == BoundaryType::INDEPENDENT) {
            interface_conditions.push_back(boundary);
        }
    }
}
```

## 数据验证和错误处理

### 1. 数据完整性检查

```cpp
bool validate_project_data(const MaxwellProjectData& data) {
    // 检查必需字段
    if (data.metadata.product.empty()) {
        LOG_CONVERTER(error, "项目产品类型为空");
        return false;
    }
    
    // 检查材料定义
    for (const auto& [name, material] : data.materials) {
        if (material.properties.conductivity && material.properties.conductivity->value <= 0) {
            LOG_CONVERTER(warn, "材料{}的电导率无效: {}", name, material.properties.conductivity->value);
        }
    }
    
    return true;
}
```

### 2. 异常处理

```cpp
try {
    auto project_data = MaxwellDataConverter::convert_from_file(file_path);
    
    if (!project_data) {
        throw std::runtime_error("转换失败");
    }
    
    if (!validate_project_data(*project_data)) {
        throw std::runtime_error("数据验证失败");
    }
    
    // 处理有效数据
    process_project_data(*project_data);
    
} catch (const std::exception& e) {
    LOG_CONVERTER(error, "处理Maxwell数据时出错: {}", e.what());
}
```

## 性能优化建议

### 1. 内存管理

- 使用 `std::optional` 避免不必要的内存分配
- 对于大型项目，考虑使用延迟加载策略
- 使用移动语义避免不必要的拷贝

### 2. 数据访问优化

```cpp
// 使用引用避免拷贝
const auto& materials = project_data.materials;

// 使用find避免重复查找
auto it = materials.find("copper");
if (it != materials.end()) {
    const auto& copper = it->second;
    // 使用材料数据
}
```

### 3. 缓存策略

对于频繁访问的数据，可以考虑缓存机制：

```cpp
class MaxwellDataCache {
private:
    std::unordered_map<std::string, MaxwellProjectData> cache_;
    
public:
    std::optional<MaxwellProjectData> get_cached_data(const std::string& file_path) {
        auto it = cache_.find(file_path);
        if (it != cache_.end()) {
            return it->second;
        }
        
        // 从文件加载并缓存
        auto data = MaxwellDataConverter::convert_from_file(file_path);
        if (data) {
            cache_[file_path] = *data;
        }
        
        return data;
    }
};
```

## 扩展和自定义

### 1. 添加新的几何操作类型

```cpp
// 扩展几何操作类型
enum class GeometryOperationType {
    // ... 现有类型
    TORUS,           // 圆环体
    CONE,            // 圆锥体
    PYRAMID,         // 棱锥体
    CUSTOM           // 自定义操作
};

// 添加对应的参数结构
struct TorusParameters {
    double major_radius;
    double minor_radius;
    Vector3D center;
    Vector3D axis;
};

// 扩展variant
using GeometryParameters = std::variant<
    BoxParameters,
    CylinderParameters,
    TorusParameters
>;
```

### 2. 自定义数据验证规则

```cpp
class CustomDataValidator {
public:
    bool validate_material(const Material& material) {
        // 自定义材料验证逻辑
        return true;
    }
    
    bool validate_boundary(const BoundaryCondition& boundary) {
        // 自定义边界条件验证逻辑
        return true;
    }
};
```

## 总结

新设计的Maxwell AEDT数据结构提供了：

1. **完整的类型安全**：使用现代C++特性确保类型安全
2. **灵活的扩展性**：支持新的几何操作和边界条件类型
3. **高效的性能**：优化的内存管理和数据访问
4. **易于集成**：与FETI-DP求解器无缝集成
5. **完善的错误处理**：详细的日志和异常处理机制

通过使用这些数据结构，可以轻松地将Maxwell项目配置转换为FETI-DP求解器所需的输入格式，大大提高开发效率和代码可维护性。