// OpenMP基础测试用例
// 测试OpenMP多线程功能和无感适配

#include <iostream>
#include <vector>
#include <numeric>
#include <cassert>
#include "omp_wrapper.h"

/**
 * @brief 测试OpenMP线程信息功能
 */
void test_omp_info() {
    std::cout << "=== 测试OpenMP线程信息 ===" << std::endl;
    
    // 测试线程信息获取
    int max_threads = fe_em::OMPWrapper::get_max_threads();
    std::cout << "最大可用线程数: " << max_threads << std::endl;
    
    int thread_id = fe_em::OMPWrapper::get_thread_id();
    std::cout << "当前线程ID: " << thread_id << std::endl;
    
    bool in_parallel = fe_em::OMPWrapper::in_parallel();
    std::cout << "是否在并行区域内: " << (in_parallel ? "是" : "否") << std::endl;
    
    // 验证单线程模式下的默认值
    assert(thread_id == 0);
    assert(!in_parallel);
    
    std::cout << "线程信息测试通过" << std::endl;
}

/**
 * @brief 测试并行for循环
 */
void test_omp_parallel_for() {
    std::cout << "\n=== 测试并行for循环 ===" << std::endl;
    
    const int N = 1000;
    std::vector<double> data(N, 0.0);
    std::vector<double> expected(N);
    
    // 生成预期结果
    for (int i = 0; i < N; ++i) {
        expected[i] = i * 2.5;
    }
    
    // 使用OpenMP并行for循环
    OMP_PARALLEL_FOR_AUTO
    for (int i = 0; i < N; ++i) {
        data[i] = i * 2.5;
    }
    
    // 验证结果
    for (int i = 0; i < N; ++i) {
        assert(data[i] == expected[i]);
    }
    
    std::cout << "并行for循环测试通过 (N=" << N << ")" << std::endl;
}

/**
 * @brief 测试临界区保护
 */
void test_omp_critical() {
    std::cout << "\n=== 测试临界区保护 ===" << std::endl;
    
    int shared_counter = 0;
    const int NUM_ITERATIONS = 10000;
    
    // 使用临界区保护共享变量
    OMP_PARALLEL_FOR_AUTO
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        OMP_CRITICAL_BLOCK(
            shared_counter++;
        )
    }
    
    // 验证结果
    assert(shared_counter == NUM_ITERATIONS);
    
    std::cout << "临界区保护测试通过 (计数器=" << shared_counter << ")" << std::endl;
}

/**
 * @brief 测试原子操作
 */
void test_omp_atomic() {
    std::cout << "\n=== 测试原子操作 ===" << std::endl;
    
    double atomic_sum = 0.0;
    const int NUM_ITERATIONS = 1000;
    
    // 使用原子操作累加
    OMP_PARALLEL_FOR_AUTO
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        OMP_ATOMIC_OPERATION(
            atomic_sum += 1.5;
        )
    }
    
    // 验证结果
    double expected_sum = NUM_ITERATIONS * 1.5;
    assert(atomic_sum == expected_sum);
    
    std::cout << "原子操作测试通过 (总和=" << atomic_sum << ")" << std::endl;
}

/**
 * @brief 测试线程局部存储
 */
void test_omp_thread_local() {
    std::cout << "\n=== 测试线程局部存储 ===" << std::endl;
    
    // 定义线程局部变量
    OMP_THREAD_LOCAL_VAR(int, thread_local_counter, 0);
    
    const int NUM_ITERATIONS = 100;
    
    // 每个线程独立累加
    OMP_PARALLEL_FOR_AUTO
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        thread_local_counter++;
    }
    
    // 注意：线程局部变量在并行区域外不可访问
    // 这里只是验证编译通过和执行无异常
    
    std::cout << "线程局部存储测试通过" << std::endl;
}

/**
 * @brief 测试线程数设置
 */
void test_omp_thread_num() {
    std::cout << "\n=== 测试线程数设置 ===" << std::endl;
    
    // 测试线程数验证功能
    int valid_threads = fe_em::OMPWrapper::get_valid_thread_num(4);
    std::cout << "请求4个线程，实际使用: " << valid_threads << "个" << std::endl;
    
    // 测试非法线程数处理
    try {
        fe_em::OMPWrapper::set_thread_num(-1);
        assert(false); // 应该抛出异常
    } catch (const fe_em::OMPException& e) {
        std::cout << "非法线程数异常捕获成功: " << e.what() << std::endl;
    }
    
    std::cout << "线程数设置测试通过" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始OpenMP基础测试" << std::endl;
    std::cout << "编译模式: " << COMPILATION_MODE << std::endl;
    
    try {
        test_omp_info();
        test_omp_parallel_for();
        test_omp_critical();
        test_omp_atomic();
        test_omp_thread_local();
        test_omp_thread_num();
        
        std::cout << "\n✅ 所有OpenMP测试通过!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}