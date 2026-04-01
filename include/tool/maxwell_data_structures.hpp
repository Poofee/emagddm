#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <variant>
#include <optional>
#include <array>

namespace fe_em {
namespace tool {

/**
 * @brief Maxwell数据模型命名空间
 * 定义用于存储Maxwell AEDT文件数据的完整数据结构
 */
namespace maxwell_data {

// 基本数据类型定义
using Vector3D = std::array<double, 3>;
using RGBColor = std::array<int, 3>;

/**
 * @brief 物理量单位枚举
 */
enum class UnitType {
    METER,          // 米
    MILLIMETER,     // 毫米
    CENTIMETER,     // 厘米
    AMPERE,         // 安培
    MILLIAMPER,     // 毫安
    VOLT,           // 伏特
    MILLIVOLT,      // 毫伏
    TESLA,          // 特斯拉
    AMPERE_PER_METER, // 安培/米
    SIEMENS_PER_METER, // 西门子/米
    HENRY_PER_METER, // 亨利/米
    FARAD_PER_METER, // 法拉/米
    DEGREE,         // 度
    RADIAN,         // 弧度
    CELSIUS,        // 摄氏度
    KELVIN,         // 开尔文
    PASCAL,         // 帕斯卡
    NEWTON,         // 牛顿
    WATT,           // 瓦特
    HERTZ,          // 赫兹
    UNKNOWN         // 未知单位
};

/**
 * @brief 带单位的物理量
 */
struct PhysicalQuantity {
    double value;           // 数值
    UnitType unit;          // 单位
    std::string raw_string; // 原始字符串表示（如"10mA"）
    
    PhysicalQuantity(double v = 0.0, UnitType u = UnitType::UNKNOWN, 
                    const std::string& raw = "")
        : value(v), unit(u), raw_string(raw) {}
};

/**
 * @brief 坐标系类型
 */
enum class CoordinateSystemType {
    CARTESIAN,      // 笛卡尔坐标系
    CYLINDRICAL,    // 柱坐标系
    SPHERICAL,      // 球坐标系
    UNKNOWN         // 未知坐标系
};

/**
 * @brief 物理类型
 */
enum class PhysicsType {
    ELECTROMAGNETIC,    // 电磁
    THERMAL,            // 热学
    STRUCTURAL,         // 结构
    UNKNOWN             // 未知
};

/**
 * @brief 材料外观数据
 */
struct MaterialAppearance {
    RGBColor color;         // RGB颜色
    double transparency;    // 透明度 (0-1)
    
    MaterialAppearance() : color({255, 255, 255}), transparency(0.0) {}
    MaterialAppearance(int r, int g, int b, double t = 0.0)
        : color({r, g, b}), transparency(t) {}
};

/**
 * @brief 材料属性
 */
struct MaterialProperties {
    // 电磁属性
    std::optional<PhysicalQuantity> permittivity;          // 介电常数
    std::optional<PhysicalQuantity> permeability;          // 磁导率
    std::optional<PhysicalQuantity> conductivity;          // 电导率
    std::optional<Vector3D> magnetic_coercivity;           // 矫顽力
    
    // 热学属性
    std::optional<PhysicalQuantity> thermal_conductivity;  // 热导率
    std::optional<PhysicalQuantity> specific_heat;         // 比热容
    std::optional<PhysicalQuantity> thermal_expansion_coefficient; // 热膨胀系数
    
    // 结构属性
    std::optional<PhysicalQuantity> youngs_modulus;        // 杨氏模量
    std::optional<PhysicalQuantity> poissons_ratio;        // 泊松比
    std::optional<PhysicalQuantity> mass_density;          // 质量密度
    
    // 外观属性
    std::optional<MaterialAppearance> appearance;          // 外观设置
};

/**
 * @brief 材料定义
 */
struct Material {
    std::string name;                           // 材料名称
    CoordinateSystemType coord_system;          // 坐标系类型
    bool is_bulk_material;                      // 是否为体材料
    std::vector<PhysicsType> physics_types;     // 支持的物理类型
    MaterialProperties properties;              // 材料属性
    std::string library;                        // 所属库
    std::string lib_location;                   // 库位置
    bool modified_since_lib;                    // 是否相对于库有修改
    
    Material() : coord_system(CoordinateSystemType::CARTESIAN), 
                is_bulk_material(true), 
                modified_since_lib(false) {}
};

/**
 * @brief 几何操作参数 - 长方体
 */
struct BoxParameters {
    Vector3D position;      // 位置 (x, y, z)
    Vector3D size;          // 尺寸 (x_size, y_size, z_size)
    int kernel_version;     // 内核版本
    
    BoxParameters() : position({0, 0, 0}), size({0, 0, 0}), kernel_version(0) {}
};

/**
 * @brief 几何操作拓扑信息
 */
struct TopologyInfo {
    int num_lumps;          // 块数
    int num_shells;         // 壳数
    int num_faces;          // 面数
    int num_wires;          // 线数
    int num_loops;          // 环数
    int num_coedges;        // 共边数
    int num_edges;          // 边数
    int num_vertices;       // 顶点数
    
    TopologyInfo() : num_lumps(0), num_shells(0), num_faces(0), 
                    num_wires(0), num_loops(0), num_coedges(0), 
                    num_edges(0), num_vertices(0) {}
};

/**
 * @brief 几何操作身份标识
 */
struct OperationIdentity {
    TopologyInfo topology;          // 拓扑信息
    int body_id;                    // 体ID
    int start_face_id;              // 起始面ID
    int start_edge_id;              // 起始边ID
    int start_vertex_id;            // 起始顶点ID
    int num_new_faces;              // 新面数
    int num_new_edges;              // 新边数
    int num_new_vertices;           // 新顶点数
    
    OperationIdentity() : body_id(-1), start_face_id(-1), start_edge_id(-1),
                         start_vertex_id(-1), num_new_faces(0), 
                         num_new_edges(0), num_new_vertices(0) {}
};

/**
 * @brief 几何操作类型
 */
enum class GeometryOperationType {
    BOX,            // 长方体
    CYLINDER,       // 圆柱体
    SPHERE,         // 球体
    EXTRUDE,        // 拉伸
    REVOLVE,        // 旋转
    SWEEP,          // 扫掠
    UNKNOWN         // 未知操作
};

/**
 * @brief 几何操作
 */
struct GeometryOperation {
    GeometryOperationType type;         // 操作类型
    int id;                            // 操作ID
    int reference_coord_system_id;     // 参考坐标系ID
    bool is_suppressed;                // 是否被抑制
    
    // 操作参数（使用variant支持多种操作类型）
    std::variant<BoxParameters> parameters;
    
    OperationIdentity identity;        // 操作身份信息
    
    GeometryOperation() : type(GeometryOperationType::UNKNOWN), id(-1),
                         reference_coord_system_id(-1), is_suppressed(false) {}
};

/**
 * @brief 几何体属性
 */
struct GeometryAttributes {
    std::string name;                   // 几何体名称
    RGBColor color;                     // 颜色
    double transparency;                // 透明度
    int part_coordinate_system;         // 部件坐标系
    std::string material_name;          // 材料名称
    std::string surface_material_name;  // 表面材料名称
    bool solve_inside;                  // 是否求解内部
    bool is_shell_element;              // 是否为壳单元
    double shell_element_thickness;     // 壳单元厚度
    bool is_material_editable;          // 材料是否可编辑
    bool use_material_appearance;       // 是否使用材料外观
    bool is_lightweight;                // 是否为轻量级
    bool is_always_hidden;              // 是否总是隐藏
    
    GeometryAttributes() : color({143, 175, 143}), transparency(0.0),
                          part_coordinate_system(1), solve_inside(true),
                          is_shell_element(false), shell_element_thickness(0.0),
                          is_material_editable(true), use_material_appearance(false),
                          is_lightweight(false), is_always_hidden(false) {}
};

/**
 * @brief 几何部件
 */
struct GeometryPart {
    GeometryAttributes attributes;                      // 几何体属性
    std::vector<GeometryOperation> operations;          // 几何操作列表
    
    GeometryPart() = default;
};

/**
 * @brief 边界条件类型
 */
enum class BoundaryType {
    INSULATING,                     // 绝缘边界
    ZERO_TANGENTIAL_H_FIELD,        // 零切向磁场
    INTEGRATED_ZERO_TANGENTIAL_H_FIELD, // 积分零切向磁场
    SYMMETRY,                       // 对称边界
    INDEPENDENT,                    // 独立边界
    DEPENDENT,                      // 依赖边界
    TANGENTIAL_H_FIELD,             // 切向磁场
    CURRENT,                        // 电流激励
    CURRENT_DENSITY,                // 电流密度激励
    CURRENT_DENSITY_TERMINAL,       // 电流密度终端
    VOLTAGE,                        // 电压激励
    VOLTAGE_DROP,                   // 电压降
    UNKNOWN                         // 未知边界类型
};

/**
 * @brief 几何位置
 */
struct GeometryPosition {
    bool is_attached_to_entity;     // 是否附着到实体
    int entity_id;                  // 实体ID
    std::vector<int> parent_ids;    // 父ID列表
    int faceted_body_triangle_index; // 面体三角形索引
    int triangle_vertex_index;      // 三角形顶点索引
    bool has_xyz;                   // 是否有XYZ坐标
    std::string position_type;      // 位置类型
    double u_param;                 // U参数
    double v_param;                 // V参数
    Vector3D position;              // 位置坐标
    
    GeometryPosition() : is_attached_to_entity(false), entity_id(-1),
                        faceted_body_triangle_index(-1), triangle_vertex_index(-1),
                        has_xyz(false), u_param(0), v_param(0), position({0, 0, 0}) {}
};

/**
 * @brief 坐标系向量
 */
struct CoordinateSystemVector {
    GeometryPosition start_point;   // 起点
    GeometryPosition end_point;     // 终点
    bool reverse_v;                 // 是否反转V方向
    
    CoordinateSystemVector() : reverse_v(false) {}
};

/**
 * @brief 边界条件
 */
struct BoundaryCondition {
    int id;                         // 边界ID
    BoundaryType type;              // 边界类型
    bool is_component;              // 是否为组件
    std::vector<int> faces;         // 关联的面ID列表
    std::vector<int> objects;       // 关联的对象ID列表
    int parent_boundary_id;         // 父边界ID
    
    // 激励源相关属性
    std::optional<PhysicalQuantity> current;           // 电流值
    std::optional<PhysicalQuantity> voltage;           // 电压值
    std::optional<Vector3D> current_density;           // 电流密度
    std::optional<PhysicalQuantity> voltage_drop;      // 电压降
    
    // 坐标系相关
    std::optional<CoordinateSystemVector> coord_system; // 坐标系向量
    std::optional<int> dependent_boundary_id;          // 依赖边界ID
    bool relation_is_same;                             // 关系是否相同
    bool is_odd_symmetry;                              // 是否为奇对称
    
    BoundaryCondition() : id(-1), type(BoundaryType::UNKNOWN), 
                         is_component(false), parent_boundary_id(-1),
                         relation_is_same(false), is_odd_symmetry(false) {}
};

/**
 * @brief 求解器类型
 */
enum class SolverType {
    MAGNETOSTATIC,      // 静磁求解
    EDDY_CURRENT,       // 涡流求解
    TRANSIENT,          // 瞬态求解
    UNKNOWN             // 未知求解类型
};

/**
 * @brief 网格方法
 */
enum class MeshMethod {
    AUTO,               // 自动网格
    TAU,                // TAU网格
    CURVILINEAR,        // 曲线网格
    UNKNOWN             // 未知方法
};

/**
 * @brief 求解器设置
 */
struct SolverSetup {
    int id;                             // 设置ID
    SolverType type;                    // 求解器类型
    bool enabled;                       // 是否启用
    
    // 收敛控制
    int maximum_passes;                 // 最大迭代次数
    int minimum_passes;                 // 最小迭代次数
    int minimum_converged_passes;       // 最小收敛迭代次数
    double percent_refinement;          // 细化百分比
    double percent_error;               // 误差百分比
    
    // 求解器选项
    bool solve_field_only;              // 仅求解场
    bool solve_matrix_at_last;          // 最后求解矩阵
    bool use_iterative_solver;          // 使用迭代求解器
    double relative_residual;           // 相对残差
    double non_linear_residual;         // 非线性残差
    bool smooth_bh_curve;               // 平滑BH曲线
    bool mu_non_linear_bh;              // 非线性BH磁导率
    
    SolverSetup() : id(-1), type(SolverType::UNKNOWN), enabled(true),
                   maximum_passes(10), minimum_passes(2), 
                   minimum_converged_passes(1), percent_refinement(30.0),
                   percent_error(1.0), solve_field_only(false),
                   solve_matrix_at_last(true), use_iterative_solver(false),
                   relative_residual(1e-6), non_linear_residual(0.001),
                   smooth_bh_curve(false), mu_non_linear_bh(true) {}
};

/**
 * @brief 网格设置
 */
struct MeshSettings {
    MeshMethod method;                          // 网格方法
    int slider_mesh_settings;                   // 网格滑块设置
    bool use_auto_length;                       // 使用自动长度
    bool apply_curvilinear;                     // 应用曲线网格
    bool dynamic_surface_resolution;            // 动态表面分辨率
    bool use_alternative_mesh_methods;          // 使用替代网格方法
    
    MeshSettings() : method(MeshMethod::AUTO), slider_mesh_settings(5),
                    use_auto_length(true), apply_curvilinear(true),
                    dynamic_surface_resolution(false), 
                    use_alternative_mesh_methods(true) {}
};

/**
 * @brief 后处理扫描变量
 */
struct PostprocessSweep {
    std::string variable;               // 变量名
    int regular_sweep;                  // 常规扫描
    std::string units;                  // 单位
    double minimum;                     // 最小值
    double maximum;                     // 最大值
    double increment;                   // 增量
    bool create_indexed_subsweep;       // 创建索引子扫描
    
    PostprocessSweep() : regular_sweep(1), minimum(0), maximum(0), 
                        increment(0), create_indexed_subsweep(false) {}
};

/**
 * @brief 求解结果
 */
struct Solution {
    int id;                             // 解ID
    std::string name;                   // 解名称
    std::vector<PostprocessSweep> sweeps; // 后处理扫描
    
    Solution() : id(-1) {}
};

/**
 * @brief 温度设置
 */
struct TemperatureSettings {
    bool include_temperature_dependence;    // 包含温度依赖性
    bool enable_feedback;                   // 启用反馈
    std::vector<std::pair<int, std::string>> temperatures; // 温度列表
    
    TemperatureSettings() : include_temperature_dependence(false),
                           enable_feedback(false) {}
};

/**
 * @brief 设计模型
 */
struct DesignModel {
    std::string name;                           // 模型名称
    int design_id;                              // 设计ID
    SolverType solution_type;                   // 求解类型
    std::string geometry_mode;                  // 几何模式
    std::string background_material;            // 背景材料
    
    TemperatureSettings temperature_settings;   // 温度设置
    std::vector<GeometryPart> geometry_parts;   // 几何部件列表
    std::vector<BoundaryCondition> boundaries;  // 边界条件列表
    SolverSetup solver_setup;                   // 求解器设置
    MeshSettings mesh_settings;                 // 网格设置
    std::vector<Solution> solutions;            // 求解结果列表
    
    DesignModel() : design_id(-1), solution_type(SolverType::UNKNOWN) {}
};

/**
 * @brief 项目元数据
 */
struct ProjectMetadata {
    std::string created_time;           // 创建时间
    std::string product;                // 产品类型
    bool file_owned_by_workbench;       // 是否由Workbench拥有
    std::pair<int, int> version;        // 版本号 (主版本, 次版本)
    std::pair<int, int> infrastructure_version; // 基础设施版本
    bool uses_advanced_features;        // 是否使用高级功能
    int next_unique_id;                 // 下一个唯一ID
    bool move_backwards;                // 是否向后移动
    
    ProjectMetadata() : file_owned_by_workbench(false), 
                       version({0, 0}), infrastructure_version({0, 0}),
                       uses_advanced_features(false), next_unique_id(0),
                       move_backwards(false) {}
};

/**
 * @brief Maxwell项目数据容器
 * 包含完整的项目配置信息
 */
struct MaxwellProjectData {
    ProjectMetadata metadata;                   // 项目元数据
    std::unordered_map<std::string, Material> materials; // 材料库
    std::vector<DesignModel> design_models;     // 设计模型列表
    
    MaxwellProjectData() = default;
};

} // namespace maxwell_data
} // namespace tool
} // namespace fe_em