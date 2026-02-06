/**
 * @file em_enums.cpp
 * @brief 基础工具层 - 电磁模块全局枚举工具源文件
 */

#include "tool/em_enums.hpp"
#include <stdexcept>

namespace tool {

std::string dimTypeToString(DimType dim_type) {
    switch (dim_type) {
        case DimType::D2:     return "D2";
        case DimType::D3:     return "D3";
        case DimType::AXIS:   return "AXIS";
        default:                    return "UNKNOWN";
    }
}

DimType stringToDimType(const std::string& str) {
    if (str == "D2" || str == "2D") return DimType::D2;
    if (str == "D3" || str == "3D") return DimType::D3;
    if (str == "AXIS" || str == "Axisymmetric") return DimType::AXIS;
    throw std::invalid_argument("Invalid DimType: " + str);
}

std::string fieldTypeToString(FieldType field_type) {
    switch (field_type) {
        case FieldType::SCALAR:  return "SCALAR";
        case FieldType::VECTOR:   return "VECTOR";
        default:                        return "UNKNOWN";
    }
}

FieldType stringToFieldType(const std::string& str) {
    if (str == "SCALAR" || str == "Scalar") return FieldType::SCALAR;
    if (str == "VECTOR" || str == "Vector") return FieldType::VECTOR;
    throw std::invalid_argument("Invalid FieldType: " + str);
}

std::string matTypeToString(MatType mat_type) {
    switch (mat_type) {
        case MatType::LINEAR_ISOTROPIC:        return "LINEAR_ISOTROPIC";
        case MatType::LINEAR_ANISOTROPIC:      return "LINEAR_ANISOTROPIC";
        case MatType::NONLINEAR_ISOTROPIC:     return "NONLINEAR_ISOTROPIC";
        case MatType::NONLINEAR_ANISOTROPIC:   return "NONLINEAR_ANISOTROPIC";
        case MatType::PERMANENT_MAGNET:        return "PERMANENT_MAGNET";
        case MatType::CONDUCTOR:               return "CONDUCTOR";
        case MatType::DIELECTRIC:              return "DIELECTRIC";
        case MatType::SUPERCONDUCTOR:          return "SUPERCONDUCTOR";
        default:                                          return "UNKNOWN";
    }
}

MatType stringToMatType(const std::string& str) {
    if (str == "LINEAR_ISOTROPIC" || str == "Linear Isotropic") return MatType::LINEAR_ISOTROPIC;
    if (str == "LINEAR_ANISOTROPIC" || str == "Linear Anisotropic") return MatType::LINEAR_ANISOTROPIC;
    if (str == "NONLINEAR_ISOTROPIC" || str == "Nonlinear Isotropic") return MatType::NONLINEAR_ISOTROPIC;
    if (str == "PERMANENT_MAGNET" || str == "Permanent Magnet") return MatType::PERMANENT_MAGNET;
    if (str == "CONDUCTOR" || str == "Conductor") return MatType::CONDUCTOR;
    if (str == "DIELECTRIC" || str == "Dielectric") return MatType::DIELECTRIC;
    if (str == "SUPERCONDUCTOR" || str == "Superconductor") return MatType::SUPERCONDUCTOR;
    throw std::invalid_argument("Invalid MatType: " + str);
}

std::string bndTypeToString(BndType bnd_type) {
    switch (bnd_type) {
        case BndType::DIRICHLET:          return "DIRICHLET";
        case BndType::NEUMANN:            return "NEUMANN";
        case BndType::ROBIN:              return "ROBIN";
        case BndType::PERIODIC:           return "PERIODIC";
        case BndType::ANTIPERIODIC:       return "ANTIPERIODIC";
        case BndType::MASTER_SLAVE:       return "MASTER_SLAVE";
        case BndType::EVEN_SYMMETRY:      return "EVEN_SYMMETRY";
        case BndType::ODD_SYMMETRY:       return "ODD_SYMMETRY";
        case BndType::BALLOON:            return "BALLOON";
        case BndType::PERFECT_E:          return "PERFECT_E";
        case BndType::PERFECT_H:          return "PERFECT_H";
        default:                                        return "UNKNOWN";
    }
}

BndType stringToBndType(const std::string& str) {
    if (str == "DIRICHLET" || str == "Dirichlet" || str == "Fixed") return BndType::DIRICHLET;
    if (str == "NEUMANN" || str == "Neumann") return BndType::NEUMANN;
    if (str == "PERIODIC" || str == "Periodic") return BndType::PERIODIC;
    if (str == "PERFECT_E" || str == "Perfect E") return BndType::PERFECT_E;
    if (str == "PERFECT_H" || str == "Perfect H") return BndType::PERFECT_H;
    if (str == "BALLOON" || str == "Balloon") return BndType::BALLOON;
    throw std::invalid_argument("Invalid BndType: " + str);
}

std::string excitationTypeToString(ExcitationType exc_type) {
    switch (exc_type) {
        case ExcitationType::COIL:               return "COIL";
        case ExcitationType::WINDING:            return "WINDING";
        case ExcitationType::CURRENT_DENSITY:    return "CURRENT_DENSITY";
        case ExcitationType::VOLTAGE_SOURCE:      return "VOLTAGE_SOURCE";
        default:                                              return "UNKNOWN";
    }
}

ExcitationType stringToExcitationType(const std::string& str) {
    if (str == "COIL" || str == "Coil") return ExcitationType::COIL;
    if (str == "WINDING" || str == "Winding") return ExcitationType::WINDING;
    if (str == "CURRENT_DENSITY" || str == "Current Density") return ExcitationType::CURRENT_DENSITY;
    if (str == "VOLTAGE_SOURCE" || str == "Voltage Source") return ExcitationType::VOLTAGE_SOURCE;
    throw std::invalid_argument("Invalid ExcitationType: " + str);
}

std::string projectFileTypeToString(ProjectFileType file_type) {
    switch (file_type) {
        case ProjectFileType::AEDT:      return "AEDT";
        case ProjectFileType::AEDTZ:     return "AEDTZ";
        case ProjectFileType::AMAT:      return "AMAT";
        case ProjectFileType::EMF:       return "EMF";
        default:                                          return "UNKNOWN";
    }
}

ProjectFileType stringToProjectFileType(const std::string& str) {
    if (str == "AEDT" || str == ".aedt") return ProjectFileType::AEDT;
    if (str == "EMF" || str == ".emf") return ProjectFileType::EMF;
    if (str == "AMAT" || str == ".amat") return ProjectFileType::AMAT;
    throw std::invalid_argument("Invalid ProjectFileType: " + str);
}

std::string maxwellVersionToString(MaxwellVersion version) {
    switch (version) {
        case MaxwellVersion::R15:     return "R15";
        case MaxwellVersion::R16:     return "R16";
        case MaxwellVersion::R17:     return "R17";
        case MaxwellVersion::R18:     return "R18";
        case MaxwellVersion::R19:     return "R19";
        case MaxwellVersion::R20:     return "R20";
        case MaxwellVersion::R21:     return "R21";
        case MaxwellVersion::R22:     return "R22";
        case MaxwellVersion::R23:     return "R23";
        case MaxwellVersion::R24:     return "R24";
        case MaxwellVersion::NEWER:   return "NEWER";
        case MaxwellVersion::UNKNOWN: return "UNKNOWN";
        default: return "UNKNOWN";
    }
}

MaxwellVersion stringToMaxwellVersion(const std::string& str) {
    if (str == "R22" || str == "22" || str == "2022") return MaxwellVersion::R22;
    if (str == "R24" || str == "24" || str == "2024") return MaxwellVersion::R24;
    if (str == "NEWER" || str == "2025") return MaxwellVersion::NEWER;
    return MaxwellVersion::UNKNOWN;
}

} // namespace tool
