// OpenMP轻量封装头文件
// 提供低侵入的OpenMP宏封装，支持无感适配

#ifndef FE_EM_OMP_WRAPPER_H
#define FE_EM_OMP_WRAPPER_H

#include "config.h"
#include <stdexcept>
#include <string>

namespace fe_em {

/**
 * @brief OpenMP线程管理异常类
 */
class OMPException : public std::runtime_error {
public:
    explicit OMPException(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * @brief OpenMP轻量封装类
 * @note 仅头文件实现，无源文件依赖
 */
class OMPWrapper {
private:
    OMPWrapper() = delete;
    ~OMPWrapper() = delete;

public:
    /**
     * @brief 获取当前线程ID
     * @return 当前线程ID（0到线程数-1）
     */
    static int get_thread_id() {
#if USE_OPENMP
        return omp_get_thread_num();
#else
        return 0; // 单线程时始终返回0
#endif
    }

    /**
     * @brief 获取总线程数
     * @return 当前并行区域的总线程数
     */
    static int get_thread_num() {
#if USE_OPENMP
        return omp_get_num_threads();
#else
        return 1; // 单线程时始终返回1
#endif
    }

    /**
     * @brief 获取最大可用线程数
     * @return 系统支持的最大线程数
     */
    static int get_max_threads() {
#if USE_OPENMP
        return omp_get_max_threads();
#else
        return 1; // 单线程时返回1
#endif
    }

    /**
     * @brief 设置并行区域使用的线程数
     * @param num_threads 线程数（必须大于0且不超过最大线程数）
     * @throws OMPException 如果线程数非法
     */
    static void set_thread_num(int num_threads) {
        if (num_threads <= 0) {
            throw OMPException("线程数必须大于0");
        }
        
        int max_threads = get_max_threads();
        if (num_threads > max_threads) {
            throw OMPException("线程数超过系统最大支持数: " + std::to_string(num_threads) + 
                             ">" + std::to_string(max_threads));
        }

#if USE_OPENMP
        omp_set_num_threads(num_threads);
#endif
    }

    /**
     * @brief 判断当前是否在并行区域内
     * @return true如果在并行区域内，false否则
     */
    static bool in_parallel() {
#if USE_OPENMP
        return omp_in_parallel();
#else
        return false; // 单线程时始终不在并行区域
#endif
    }

    /**
     * @brief 获取线程数（带验证）
     * @param num_threads 请求的线程数
     * @return 实际使用的线程数（经过验证和调整）
     */
    static int get_valid_thread_num(int num_threads) {
        if (num_threads <= 0) {
            return 1; // 默认使用1个线程
        }
        
        int max_threads = get_max_threads();
        return (num_threads > max_threads) ? max_threads : num_threads;
    }
};

} // namespace fe_em

// ============================================================================
// OpenMP宏封装
// ============================================================================

#if USE_OPENMP

// OpenMP启用时的原生宏封装

/**
 * @brief 并行区域宏
 * @param ... 可选的线程数设置
 */
#define OMP_PARALLEL(...) _Pragma("omp parallel __VA_ARGS__")

/**
 * @brief 并行for循环宏
 * @param ... 循环参数和调度策略
 */
#define OMP_FOR(...) _Pragma("omp for __VA_ARGS__")

/**
 * @brief 并行区域+for循环组合宏
 * @param ... 线程数和循环参数
 */
#define OMP_PARALLEL_FOR(...) _Pragma("omp parallel for __VA_ARGS__")

/**
 * @brief 临界区宏
 */
#define OMP_CRITICAL _Pragma("omp critical")

/**
 * @brief 原子操作宏
 */
#define OMP_ATOMIC _Pragma("omp atomic")

/**
 * @brief 屏障同步宏
 */
#define OMP_BARRIER _Pragma("omp barrier")

/**
 * @brief 单线程执行宏（仅一个线程执行）
 */
#define OMP_SINGLE _Pragma("omp single")

/**
 * @brief 主线程执行宏（仅主线程执行）
 */
#define OMP_MASTER _Pragma("omp master")

/**
 * @brief 线程局部存储宏
 */
#define OMP_THREAD_LOCAL _Thread_local

#else

// OpenMP禁用时的空展开宏（无感适配）

/**
 * @brief 并行区域宏（空展开）
 */
#define OMP_PARALLEL(...)

/**
 * @brief 并行for循环宏（空展开）
 */
#define OMP_FOR(...)

/**
 * @brief 并行区域+for循环组合宏（空展开）
 */
#define OMP_PARALLEL_FOR(...)

/**
 * @brief 临界区宏（空展开）
 */
#define OMP_CRITICAL

/**
 * @brief 原子操作宏（空展开）
 */
#define OMP_ATOMIC

/**
 * @brief 屏障同步宏（空展开）
 */
#define OMP_BARRIER

/**
 * @brief 单线程执行宏（空展开）
 */
#define OMP_SINGLE

/**
 * @brief 主线程执行宏（空展开）
 */
#define OMP_MASTER

/**
 * @brief 线程局部存储宏（使用thread_local）
 */
#define OMP_THREAD_LOCAL thread_local

#endif

// ============================================================================
// 高级封装宏（简化使用）
// ============================================================================

/**
 * @brief 并行区域宏（带线程数设置）
 * @param num_threads 线程数（可选，默认使用系统最大线程数）
 */
#define OMP_PARALLEL_REGION(num_threads) \
    OMP_PARALLEL(if(fe_em::OMPWrapper::get_valid_thread_num(num_threads) > 1) \
                 num_threads(fe_em::OMPWrapper::get_valid_thread_num(num_threads)))

/**
 * @brief 并行for循环宏（带调度策略）
 * @param schedule_type 调度类型（static/dynamic/guided/auto）
 * @param chunk_size 块大小（可选）
 */
#define OMP_FOR_SCHEDULE(schedule_type, chunk_size) \
    OMP_FOR(schedule(schedule_type chunk_size))

/**
 * @brief 并行区域+for循环组合宏（完整参数）
 * @param num_threads 线程数
 * @param schedule_type 调度类型
 * @param chunk_size 块大小
 */
#define OMP_PARALLEL_FOR_FULL(num_threads, schedule_type, chunk_size) \
    OMP_PARALLEL_FOR(if(fe_em::OMPWrapper::get_valid_thread_num(num_threads) > 1) \
                     num_threads(fe_em::OMPWrapper::get_valid_thread_num(num_threads)) \
                     schedule(schedule_type chunk_size))

/**
 * @brief 简化并行for循环宏（自动线程数）
 */
#define OMP_PARALLEL_FOR_AUTO OMP_PARALLEL_FOR_FULL(0, "auto", )

/**
 * @brief 临界区保护宏（自动加锁解锁）
 * @param code_block 需要保护的代码块
 */
#define OMP_CRITICAL_BLOCK(code_block) \
    OMP_CRITICAL \
    code_block

/**
 * @brief 原子操作保护宏
 * @param operation 原子操作表达式
 */
#define OMP_ATOMIC_OPERATION(operation) \
    OMP_ATOMIC \
    operation

/**
 * @brief 线程局部变量定义宏
 * @param type 变量类型
 * @param name 变量名
 * @param initial_value 初始值
 */
#define OMP_THREAD_LOCAL_VAR(type, name, initial_value) \
    OMP_THREAD_LOCAL type name = initial_value

#endif // FE_EM_OMP_WRAPPER_H