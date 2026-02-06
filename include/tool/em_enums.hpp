/**
 * @file em_enums.hpp
 * @brief 基础工具层 - 电磁模块全局枚举工具头文件（简化版）
 */

#pragma once

#include <string>

namespace tool {

enum class DimType { 
    D2,     ///< 二维问题
    D3,     ///< 三维问题
    AXIS    ///< 轴对称问题
};

enum class FieldType { 
    SCALAR, ///< 标量场
    VECTOR  ///< 矢量场
};

enum class SimulationType { 
    ELECTROSTATIC,  ///< 静电场
    MAGNETOSTATIC,  ///< 静磁场
    EDDYCURRENT,    ///< 涡流场
    TRANSIENT,      ///< 瞬态场
    HARMONIC        ///< 谐波场
};

enum class SolverType { 
    AUTO,           ///< 自动选择求解器
    ITERATIVE,      ///< 迭代求解器
    DIRECT,         ///< 直接求解器
    MOONEY_RIVLIN   ///< Mooney-Rivlin材料模型求解器
};

enum class ConvergenceType { 
    RESIDUAL,       ///< 残差收敛
    ENERGY,         ///< 能量收敛
    MIXED           ///< 混合收敛
};

enum class MatType {
    LINEAR_ISOTROPIC,      ///< 线性各向同性材料
    LINEAR_ANISOTROPIC,    ///< 线性各向异性材料
    NONLINEAR_ISOTROPIC,   ///< 非线性各向同性材料
    NONLINEAR_ANISOTROPIC, ///< 非线性各向异性材料
    PERMANENT_MAGNET,      ///< 永磁体材料
    CONDUCTOR,             ///< 导体材料
    DIELECTRIC,            ///< 电介质材料
    SUPERCONDUCTOR         ///< 超导材料
};

enum class CoreLossModelType {
    NONE,                   ///< 无铁损模型
    STEINMETZ,              ///< Steinmetz铁损模型
    GENERALIZED_STEINMETZ,  ///< 广义Steinmetz铁损模型
    IORIO_VARALDI,          ///< Iorio-Varaldi铁损模型
    PREISACH,               ///< Preisach磁滞模型
    Bertotti,               ///< Bertotti铁损模型
    CUSTOM                  ///< 自定义铁损模型
};

enum class BHCurveType { 
    NONE,           ///< 无BH曲线
    SINGLE_CURVE,   ///< 单条BH曲线
    TEMP_DEPENDENT, ///< 温度相关BH曲线
    FREQ_DEPENDENT, ///< 频率相关BH曲线
    CUSTOM_CURVE    ///< 自定义BH曲线
};

enum class BndType {
    DIRICHLET,          ///< 狄利克雷边界条件
    NEUMANN,            ///< 诺伊曼边界条件
    ROBIN,              ///< 罗宾边界条件
    PERIODIC,           ///< 周期性边界条件
    ANTIPERIODIC,       ///< 反周期性边界条件
    MASTER_SLAVE,       ///< 主从边界条件
    EVEN_SYMMETRY,      ///< 偶对称边界条件
    ODD_SYMMETRY,       ///< 奇对称边界条件
    BALLOON,            ///< 气球边界条件
    IMPEDANCE,          ///< 阻抗边界条件
    RADIATION,          ///< 辐射边界条件
    INFINITE_BOX,       ///< 无限远边界条件
    PERFECT_E,          ///< 理想电边界条件
    PERFECT_H,          ///< 理想磁边界条件
    LUMPED_PORT,        ///< 集总端口边界条件
    WAVE_PORT,          ///< 波导端口边界条件
    SURFACE_CURRENT,    ///< 表面电流边界条件
    VOLUME_CURRENT,     ///< 体电流边界条件
    SOURCE_SURFACE,     ///< 源表面边界条件
    FORCE_LOAD,         ///< 力载荷边界条件
    TORQUE_LOAD,        ///< 扭矩载荷边界条件
    INSULATION,         ///< 绝缘边界条件
    CONTACTS            ///< 接触边界条件
};

enum class BoundarySubType { 
    NONE,               ///< 无子类型
    SKIN_DEPTH,         ///< 趋肤深度子类型
    EDDY_CURRENT,       ///< 涡流子类型
    PROXIMITY_EFFECT,   ///< 邻近效应子类型
    EDGE_BASED,         ///< 基于边的子类型
    FACE_BASED          ///< 基于面的子类型
};

enum class PeriodicMappingType { 
    NONE,               ///< 无映射
    TRANSLATIONAL,      ///< 平移映射
    ROTATIONAL,         ///< 旋转映射
    CYCLIC              ///< 循环映射
};

enum class ExcitationType {
    CURRENT_DENSITY,     ///< 电流密度激励
    VOLTAGE_SOURCE,      ///< 电压源激励
    CURRENT_SOURCE,      ///< 电流源激励
    MAGNETIC_FLUX,       ///< 磁通激励
    COIL,                ///< 线圈激励
    WINDING,             ///< 绕组激励
    SOLID_CONDUCTOR,     ///< 实心导体激励
    EDDY_CURRENT,        ///< 涡流激励
    EXTERNAL_CIRCUIT,    ///< 外电路激励
    MOVING_PART,         ///< 运动部件激励
    FORCE,               ///< 力激励
    TORQUE,              ///< 扭矩激励
    THERMAL_HEAT,        ///< 热激励
    DISPLACEMENT_CURRENT ///< 位移电流激励
};

enum class ExcitationWaveformType {
    DC,             ///< 直流波形
    SINUSOIDAL,     ///< 正弦波形
    COSINE,         ///< 余弦波形
    SQUARE,         ///< 方波波形
    TRIANGULAR,     ///< 三角波形
    SAWTOOTH,       ///< 锯齿波形
    PULSE,          ///< 脉冲波形
    RAMP,           ///< 斜坡波形
    EXPONENTIAL,    ///< 指数波形
    CUSTOM,         ///< 自定义波形
    USER_DEFINED    ///< 用户定义波形
};

enum class CoilConnectionType { 
    SERIES,             ///< 串联连接
    PARALLEL,           ///< 并联连接
    SERIES_PARALLEL,    ///< 串并联连接
    DELTA,              ///< 三角形连接
    WYE                 ///< 星形连接
};

enum class WindingType { 
    SOLID,              ///< 实心绕组
    STRANDED,           ///< 绞线绕组
    LITZ,               ///< 利兹线绕组
    FOIL                ///< 箔绕组
};

enum class MotionType { 
    NONE,               ///< 无运动
    ROTATION,           ///< 旋转运动
    TRANSLATION,        ///< 平移运动
    CUSTOM_MOTION       ///< 自定义运动
};

enum class CoordinateSystemType { 
    GLOBAL,             ///< 全局坐标系
    LOCAL,              ///< 局部坐标系
    CYLINDRICAL,        ///< 柱坐标系
    SPHERICAL           ///< 球坐标系
};

enum class MeshElementType {
    NODE,           ///< 节点
    EDGE,           ///< 边
    FACE_TRI,       ///< 三角形面
    FACE_QUAD,      ///< 四边形面
    VOL_TET,        ///< 四面体体单元
    VOL_HEX,        ///< 六面体体单元
    VOL_PRISM,      ///< 棱柱体体单元
    VOL_PYRAMID     ///< 棱锥体体单元
};

enum class MeshGenerationType { 
    AUTOMATIC,          ///< 自动网格生成
    USER_DEFINED,       ///< 用户定义网格
    MAPPED,             ///< 映射网格
    COONS,              ///< Coons曲面网格
    ADVANCING_FRONT     ///< 前沿推进网格
};

enum class MeshRefinementType {
    NONE,               ///< 无网格加密
    BOUNDARY_LAYER,     ///< 边界层加密
    CURVATURE,          ///< 曲率加密
    ACTIVITY,          ///< 活动区域加密
    CORELOSS,          ///< 铁损区域加密
    SKIN_DEPTH,        ///< 趋肤深度加密
    USER_SPECIFIED     ///< 用户指定加密
};

enum class MeshQualityMetric { 
    SKEWNESS,          ///< 偏斜度
    ASPECT_RATIO,      ///< 长宽比
    JACOBIAN,          ///< 雅可比矩阵
    WARPAGE,           ///< 翘曲度
    TAPER,             ///< 锥度
    MIN_ANGLE,         ///< 最小角度
    MAX_ANGLE          ///< 最大角度
};

enum class HPCParallelMode { 
    SERIAL,         ///< 串行模式
    MPI,            ///< MPI并行
    OPENMP,         ///< OpenMP并行
    MPI_OPENMP,     ///< MPI+OpenMP混合并行
    CUDA,           ///< CUDA GPU并行
    HIP,            ///< HIP GPU并行
    ONEAPI          ///< oneAPI并行
};

enum class HPCSolverMode { 
    DISTRIBUTED,        ///< 分布式求解
    SHARED_MEMORY,      ///< 共享内存求解
    GPU_ACCELERATED,    ///< GPU加速求解
    CLOUD               ///< 云端求解
};

enum class DomainDecompositionType { 
    GEOMETRIC,      ///< 几何分解
    ALGEBRAIC,      ///< 代数分解
    RECURSIVE,      ///< 递归分解
    KWAY,           ///< K路分解
    SPECTRAL        ///< 谱分解
};

enum class ProjectFileType { 
    AEDT,           ///< ANSYS AEDT格式
    AEDTZ,          ///< ANSYS AEDT压缩格式
    XML,            ///< XML格式
    AMAT,           ///< 材料文件格式
    AEDT_AUTO,      ///< ANSYS AEDT自动检测格式
    EMF,            ///< 电磁场格式
    CUSTOM_XML,     ///< 自定义XML格式
    STEP,           ///< STEP几何格式
    IGES,           ///< IGES几何格式
    SAT,            ///< ACIS SAT几何格式
    STL             ///< STL几何格式
};

enum class MaxwellVersion { 
    UNKNOWN,        ///< 未知版本
    R15,            ///< Maxwell R15版本
    R16,            ///< Maxwell R16版本
    R17,            ///< Maxwell R17版本
    R18,            ///< Maxwell R18版本
    R19,            ///< Maxwell R19版本
    R20,            ///< Maxwell R20版本
    R21,            ///< Maxwell R21版本
    R22,            ///< Maxwell R22版本
    R23,            ///< Maxwell R23版本
    R24,            ///< Maxwell R24版本
    NEWER           ///< 更新版本
};

enum class ProjectState { 
    CREATED,        ///< 项目已创建
    LOADED,         ///< 项目已加载
    MODIFIED,       ///< 项目已修改
    SAVED,          ///< 项目已保存
    CLOSING         ///< 项目正在关闭
};

enum class IDCategory {
    MATERIAL,           ///< 材料ID
    GEOMETRY,           ///< 几何ID
    BOUNDARY,           ///< 边界条件ID
    EXCITATION,         ///< 激励ID
    MESH,               ///< 网格ID
    UNIT,               ///< 单位ID
    DOF,                ///< 自由度ID
    PROJECT,            ///< 项目ID
    RESULT,             ///< 结果ID
    UDP,                ///< 用户定义参数ID
    WIRE,               ///< 线缆ID
    COIL,               ///< 线圈ID
    WINDING,            ///< 绕组ID
    SCALING_GROUP,      ///< 缩放组ID
    PARAMETER,          ///< 参数ID
    VARIABLE            ///< 变量ID
};

std::string dimTypeToString(DimType dim_type);
DimType stringToDimType(const std::string& str);
std::string fieldTypeToString(FieldType field_type);
FieldType stringToFieldType(const std::string& str);
std::string matTypeToString(MatType mat_type);
MatType stringToMatType(const std::string& str);
std::string bndTypeToString(BndType bnd_type);
BndType stringToBndType(const std::string& str);
std::string excitationTypeToString(ExcitationType exc_type);
ExcitationType stringToExcitationType(const std::string& str);
std::string projectFileTypeToString(ProjectFileType file_type);
ProjectFileType stringToProjectFileType(const std::string& str);
std::string maxwellVersionToString(MaxwellVersion version);
MaxwellVersion stringToMaxwellVersion(const std::string& str);

} // namespace tool
