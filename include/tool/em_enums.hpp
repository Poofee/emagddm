/**
 * @file em_enums.hpp
 * @brief 基础工具层 - 电磁模块全局枚举工具头文件（简化版）
 */

#pragma once

#include <string>

namespace tool {

enum class DimType { D2, D3, AXIS };
enum class FieldType { SCALAR, VECTOR };
enum class SimulationType { ELECTROSTATIC, MAGNETOSTATIC, EDDYCURRENT, TRANSIENT, HARMONIC };
enum class SolverType { AUTO, ITERATIVE, DIRECT, MOONEY_RIVLIN };
enum class ConvergenceType { RESIDUAL, ENERGY, MIXED };

enum class MatType {
    LINEAR_ISOTROPIC, LINEAR_ANISOTROPIC, NONLINEAR_ISOTROPIC, NONLINEAR_ANISOTROPIC,
    PERMANENT_MAGNET, CONDUCTOR, DIELECTRIC, SUPERCONDUCTOR
};

enum class CoreLossModelType {
    NONE, STEINMETZ, GENERALIZED_STEINMETZ, IORIO_VARALDI, PREISACH, Bertotti, CUSTOM
};

enum class BHCurveType { NONE, SINGLE_CURVE, TEMP_DEPENDENT, FREQ_DEPENDENT, CUSTOM_CURVE };

enum class BndType {
    DIRICHLET, NEUMANN, ROBIN, PERIODIC, ANTIPERIODIC, MASTER_SLAVE,
    EVEN_SYMMETRY, ODD_SYMMETRY, BALLOON, IMPEDANCE, RADIATION, INFINITE_BOX,
    PERFECT_E, PERFECT_H, LUMPED_PORT, WAVE_PORT, SURFACE_CURRENT, VOLUME_CURRENT,
    SOURCE_SURFACE, FORCE_LOAD, TORQUE_LOAD, INSULATION, CONTACTS
};

enum class BoundarySubType { NONE, SKIN_DEPTH, EDDY_CURRENT, PROXIMITY_EFFECT, EDGE_BASED, FACE_BASED };
enum class PeriodicMappingType { NONE, TRANSLATIONAL, ROTATIONAL, CYCLIC };

enum class ExcitationType {
    CURRENT_DENSITY, VOLTAGE_SOURCE, CURRENT_SOURCE, MAGNETIC_FLUX, COIL, WINDING,
    SOLID_CONDUCTOR, EDDY_CURRENT, EXTERNAL_CIRCUIT, MOVING_PART, FORCE, TORQUE,
    THERMAL_HEAT, DISPLACEMENT_CURRENT
};

enum class ExcitationWaveformType {
    DC, SINUSOIDAL, COSINE, SQUARE, TRIANGULAR, SAWTOOTH, PULSE, RAMP, EXPONENTIAL, CUSTOM, USER_DEFINED
};

enum class CoilConnectionType { SERIES, PARALLEL, SERIES_PARALLEL, DELTA, WYE };
enum class WindingType { SOLID, STRANDED, LITZ, FOIL };
enum class MotionType { NONE, ROTATION, TRANSLATION, CUSTOM_MOTION };
enum class CoordinateSystemType { GLOBAL, LOCAL, CYLINDRICAL, SPHERICAL };

enum class MeshElementType {
    NODE, EDGE, FACE_TRI, FACE_QUAD, VOL_TET, VOL_HEX, VOL_PRISM, VOL_PYRAMID
};

enum class MeshGenerationType { AUTOMATIC, USER_DEFINED, MAPPED, COONS, ADVANCING_FRONT };
enum class MeshRefinementType {
    NONE, BOUNDARY_LAYER, CURVATURE, ACTIVITY, CORELOSS, SKIN_DEPTH, USER_SPECIFIED
};
enum class MeshQualityMetric { SKEWNESS, ASPECT_RATIO, JACOBIAN, WARPAGE, TAPER, MIN_ANGLE, MAX_ANGLE };

enum class HPCParallelMode { SERIAL, MPI, OPENMP, MPI_OPENMP, CUDA, HIP, ONEAPI };
enum class HPCSolverMode { DISTRIBUTED, SHARED_MEMORY, GPU_ACCELERATED, CLOUD };
enum class DomainDecompositionType { GEOMETRIC, ALGEBRAIC, RECURSIVE, KWAY, SPECTRAL };

enum class ProjectFileType { AEDT, AEDTZ, XML, AMAT, AEDT_AUTO, EMF, CUSTOM_XML, STEP, IGES, SAT, STL };
enum class MaxwellVersion { UNKNOWN, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, NEWER };
enum class ProjectState { CREATED, LOADED, MODIFIED, SAVED, CLOSING };

enum class IDCategory {
    MATERIAL, GEOMETRY, BOUNDARY, EXCITATION, MESH, UNIT, DOF,
    PROJECT, RESULT, UDP, WIRE, COIL, WINDING, SCALING_GROUP, PARAMETER, VARIABLE
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
