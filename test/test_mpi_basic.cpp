// MPI基础测试用例
// 测试MPI多进程功能和无感适配

#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>
#include "mpi_wrapper.h"

/**
 * @brief 测试MPI进程信息功能
 */
void test_mpi_info() {
    std::cout << "=== 测试MPI进程信息 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    
    // 测试进程信息获取
    int rank = mpi_comm.get_rank();
    int size = mpi_comm.get_size();
    bool is_root = mpi_comm.is_root();
    
    std::cout << "进程ID: " << rank << "/" << size 
              << ", 是否主进程: " << (is_root ? "是" : "否") << std::endl;
    
    // 验证进程信息
    assert(rank >= 0 && rank < size);
    assert(size > 0);
    assert(is_root == (rank == 0));
    
    std::cout << "进程信息测试通过" << std::endl;
}

/**
 * @brief 测试屏障同步
 */
void test_mpi_barrier() {
    std::cout << "=== 测试屏障同步 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    
    // 每个进程打印消息
    std::cout << "进程 " << mpi_comm.get_rank() << " 进入屏障前" << std::endl;
    
    // 执行屏障同步
    mpi_comm.barrier();
    
    std::cout << "进程 " << mpi_comm.get_rank() << " 离开屏障后" << std::endl;
    
    std::cout << "屏障同步测试通过" << std::endl;
}

/**
 * @brief 测试点对点通信
 */
void test_mpi_send_recv() {
    std::cout << "=== 测试点对点通信 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    int rank = mpi_comm.get_rank();
    int size = mpi_comm.get_size();
    
    if (size < 2) {
        std::cout << "进程数不足2，跳过点对点通信测试" << std::endl;
        return;
    }
    
    // 测试整数发送接收
    if (rank == 0) {
        int send_data = 42;
        mpi_comm.send(&send_data, 1, 1);
        std::cout << "进程0发送数据: " << send_data << std::endl;
    } else if (rank == 1) {
        int recv_data;
        int count = mpi_comm.recv(&recv_data, 1, 0);
        assert(count == 1);
        assert(recv_data == 42);
        std::cout << "进程1接收数据: " << recv_data << std::endl;
    }
    
    mpi_comm.barrier();
    std::cout << "点对点通信测试通过" << std::endl;
}

/**
 * @brief 测试广播通信
 */
void test_mpi_broadcast() {
    std::cout << "=== 测试广播通信 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    int rank = mpi_comm.get_rank();
    
    // 测试数组广播
    std::vector<double> data(5);
    
    if (rank == 0) {
        // 主进程初始化数据
        for (int i = 0; i < 5; ++i) {
            data[i] = i * 1.5;
        }
        std::cout << "进程0广播数据: ";
        for (double val : data) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    
    // 广播数据
    mpi_comm.broadcast(data.data(), 5, 0);
    
    // 所有进程验证数据
    for (int i = 0; i < 5; ++i) {
        assert(data[i] == i * 1.5);
    }
    
    std::cout << "进程" << rank << "接收广播数据: ";
    for (double val : data) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "广播通信测试通过" << std::endl;
}

/**
 * @brief 测试收集通信
 */
void test_mpi_gather() {
    std::cout << "=== 测试收集通信 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    int rank = mpi_comm.get_rank();
    int size = mpi_comm.get_size();
    
    // 每个进程准备数据
    std::vector<int> send_data(3);
    for (int i = 0; i < 3; ++i) {
        send_data[i] = rank * 10 + i;
    }
    
    std::cout << "进程" << rank << "发送数据: ";
    for (int val : send_data) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // 主进程准备接收缓冲区
    std::vector<int> recv_data;
    std::vector<int> recv_counts;
    
    if (rank == 0) {
        recv_data.resize(3 * size);
        recv_counts.resize(size, 3); // 每个进程发送3个数据
    }
    
    // 收集数据
    mpi_comm.gather(send_data.data(), 3, 
                    recv_data.data(), recv_counts.data(), 0);
    
    // 主进程验证收集结果
    if (rank == 0) {
        std::cout << "进程0收集数据: ";
        for (int val : recv_data) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        
        // 验证数据正确性
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < 3; ++j) {
                int expected = i * 10 + j;
                int actual = recv_data[i * 3 + j];
                assert(actual == expected);
            }
        }
    }
    
    std::cout << "收集通信测试通过" << std::endl;
}

/**
 * @brief 测试数据分域功能
 */
void test_mpi_split_data() {
    std::cout << "=== 测试数据分域功能 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    int rank = mpi_comm.get_rank();
    int size = mpi_comm.get_size();
    
    // 测试数据分域
    const int TOTAL_SIZE = 17;
    auto [start_index, chunk_size] = mpi_comm.split_data(TOTAL_SIZE);
    
    std::cout << "进程" << rank << "负责数据: 起始索引=" << start_index 
              << ", 数据个数=" << chunk_size << std::endl;
    
    // 验证数据分域的合理性
    assert(start_index >= 0);
    assert(chunk_size >= 0);
    assert(start_index + chunk_size <= TOTAL_SIZE);
    
    // 验证数据完整性（所有进程的数据覆盖整个范围）
    if (rank == 0) {
        int total_covered = 0;
        for (int i = 0; i < size; ++i) {
            auto [s, c] = fe_em::MPIComm::split_data(TOTAL_SIZE, i, size);
            total_covered += c;
        }
        assert(total_covered == TOTAL_SIZE);
        std::cout << "数据完整性验证通过: 总数据量=" << total_covered << std::endl;
    }
    
    std::cout << "数据分域测试通过" << std::endl;
}

/**
 * @brief 测试异常处理
 */
void test_mpi_exception() {
    std::cout << "=== 测试异常处理 ===" << std::endl;
    
    auto& mpi_comm = fe_em::MPIComm::get_instance();
    
    // 测试非法进程ID异常
    try {
        int data = 1;
        mpi_comm.send(&data, 1, -1); // 非法进程ID
        assert(false); // 应该抛出异常
    } catch (const fe_em::MPIException& e) {
        std::cout << "非法进程ID异常捕获成功: " << e.what() << std::endl;
    }
    
    // 测试非法数据个数异常
    try {
        int data = 1;
        mpi_comm.send(&data, -1, 0); // 非法数据个数
        assert(false); // 应该抛出异常
    } catch (const fe_em::MPIException& e) {
        std::cout << "非法数据个数异常捕获成功: " << e.what() << std::endl;
    }
    
    std::cout << "异常处理测试通过" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main(int argc, char** argv) {
    std::cout << "开始MPI基础测试" << std::endl;
    std::cout << "编译模式: " << COMPILATION_MODE << std::endl;
    
    try {
        // 初始化MPI
        auto& mpi_comm = fe_em::MPIComm::get_instance();
        mpi_comm.init(&argc, &argv);
        
        test_mpi_info();
        test_mpi_barrier();
        test_mpi_send_recv();
        test_mpi_broadcast();
        test_mpi_gather();
        test_mpi_split_data();
        test_mpi_exception();
        
        // 终止MPI
        mpi_comm.finalize();
        
        std::cout << "\n✅ 所有MPI测试通过!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ 测试失败: " << e.what() << std::endl;
        
        // 确保MPI正确终止
        try {
            auto& mpi_comm = fe_em::MPIComm::get_instance();
            if (mpi_comm.is_initialized() && !mpi_comm.is_finalized()) {
                mpi_comm.finalize();
            }
        } catch (...) {
            // 忽略终止异常
        }
        
        return 1;
    }
}