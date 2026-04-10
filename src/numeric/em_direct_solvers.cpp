/**
 * @file em_direct_solvers.cpp
 * @brief 求解器后端辅助类实现（SuperLU/MUMPS/DirectBackendManager）
 * @details 实现辅助数据结构封装类的全部方法：
 *          - DirectBackendManager：后端可用性管理器
 *
 * @note 原三个直接求解器类已迁移至新架构：
 *       - UnifiedDirectSolver + SolverBackend 策略模式
 *       - EigenLLTBackend / EigenLDLTBackend / EigenLUBackend
 *       - MUMPSBackend / SuperLUBackend
 *
 * @see UnifiedDirectSolver 统一调度层
 * @see SolverBackend 策略接口
 *
 * @author Poofee
 * @date 2026-04-06
 * @version 2.0 (重构版 - 移除旧求解器类)
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"
#include <algorithm>

namespace numeric {

// ============================================================================
// DirectBackendManager 实现
// ============================================================================

bool DirectBackendManager::isBackendAvailable(DirectBackendType type) {
    switch (type) {
        case DirectBackendType::EIGEN:
            return true;  // Eigen 始终可用

#ifdef HAVE_MUMPS
        case DirectBackendType::MUMPS:
            return true;
#endif

#ifdef HAVE_SUPERLU
        case DirectBackendType::SUPERLU_MT:
            return true;
#endif

        default:
            return false;
    }
}

std::string DirectBackendManager::getBackendName(DirectBackendType type) {
    switch (type) {
        case DirectBackendType::EIGEN:
            return "Eigen";
        case DirectBackendType::MUMPS:
            return "MUMPS";
        case DirectBackendType::SUPERLU_MT:
            return "SuperLU_MT";
        default:
            return "Unknown";
    }
}

std::vector<DirectBackendType> DirectBackendManager::getAvailableBackends() {
    std::vector<DirectBackendType> backends;
    backends.push_back(DirectBackendType::EIGEN);

#ifdef HAVE_MUMPS
    backends.push_back(DirectBackendType::MUMPS);
#endif

#ifdef HAVE_SUPERLU
    backends.push_back(DirectBackendType::SUPERLU_MT);
#endif

    return backends;
}

} // namespace numeric
