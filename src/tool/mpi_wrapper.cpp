// MPI面向对象封装实现文件

#include "mpi_wrapper.h"
#include <iostream>

namespace fe_em {

// 静态成员初始化
std::unique_ptr<MPIComm> MPIComm::instance_ = nullptr;

// 私有构造函数
MPIComm::MPIComm() 
    : initialized_(false), finalized_(false), rank_(0), size_(1) {
    // 单进程模式下的默认值
}

// 析构函数
MPIComm::~MPIComm() {
    if (initialized_ && !finalized_) {
        try {
            finalize();
        } catch (const MPIException& e) {
            // 析构函数中不能抛出异常，记录错误即可
            std::cerr << "MPI终止失败: " << e.what() << std::endl;
        }
    }
}

// 获取单例实例
MPIComm& MPIComm::get_instance() {
    if (!instance_) {
        instance_ = std::unique_ptr<MPIComm>(new MPIComm());
    }
    return *instance_;
}

// 初始化MPI
void MPIComm::init(int* argc, char*** argv) {
#if USE_MPI
    if (initialized_) {
        return; // 已经初始化过
    }
    
    int error_code;
    
    // 检查MPI是否已初始化
    int flag;
    error_code = MPI_Initialized(&flag);
    check_error(error_code, "MPI_Initialized");
    
    if (flag) {
        // MPI已经被其他代码初始化
        initialized_ = true;
    } else {
        // 初始化MPI
        error_code = MPI_Init(argc, argv);
        check_error(error_code, "MPI_Init");
        initialized_ = true;
    }
    
    // 获取进程信息
    error_code = MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    check_error(error_code, "MPI_Comm_rank");
    
    error_code = MPI_Comm_size(MPI_COMM_WORLD, &size_);
    check_error(error_code, "MPI_Comm_size");
    
    std::cout << "MPI初始化成功: 进程 " << rank_ << "/" << size_ << std::endl;
    
#else
    // 单进程模式：设置默认值
    initialized_ = true;
    rank_ = 0;
    size_ = 1;
    
    std::cout << "MPI单进程模式: 进程 0/1" << std::endl;
#endif
}

// 终止MPI
void MPIComm::finalize() {
#if USE_MPI
    if (finalized_ || !initialized_) {
        return; // 已经终止或未初始化
    }
    
    int error_code;
    
    // 检查MPI是否已终止
    int flag;
    error_code = MPI_Finalized(&flag);
    check_error(error_code, "MPI_Finalized");
    
    if (flag) {
        // MPI已经被其他代码终止
        finalized_ = true;
    } else {
        // 终止MPI
        error_code = MPI_Finalize();
        check_error(error_code, "MPI_Finalize");
        finalized_ = true;
    }
    
    std::cout << "MPI终止成功: 进程 " << rank_ << std::endl;
    
#else
    // 单进程模式：标记为已终止
    finalized_ = true;
    std::cout << "MPI单进程模式终止" << std::endl;
#endif
}

// 屏障同步
void MPIComm::barrier() const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法进行屏障同步");
    }
    
    int error_code = MPI_Barrier(MPI_COMM_WORLD);
    check_error(error_code, "MPI_Barrier");
#else
    // 单进程模式：无操作
#endif
}

// 检查MPI错误码
void MPIComm::check_error(int error_code, const std::string& operation_name) const {
#if USE_MPI
    if (error_code != MPI_SUCCESS) {
        char error_string[MPI_MAX_ERROR_STRING];
        int result_len;
        MPI_Error_string(error_code, error_string, &result_len);
        
        std::string msg = operation_name + "失败: " + std::string(error_string);
        throw MPIException(msg, error_code, rank_);
    }
#else
    // 单进程模式：忽略错误检查
    (void)error_code;
    (void)operation_name;
#endif
}

// 数据分域辅助函数
std::pair<int, int> MPIComm::split_data(int total_size, int rank, int size) {
    if (total_size <= 0) {
        return {0, 0};
    }
    
    if (rank < 0 || rank >= size) {
        throw MPIException("进程ID越界: " + std::to_string(rank));
    }
    
    // 计算每个进程的基本数据量
    int base_chunk = total_size / size;
    int remainder = total_size % size;
    
    // 前remainder个进程多分配一个数据
    int start_index, chunk_size;
    if (rank < remainder) {
        chunk_size = base_chunk + 1;
        start_index = rank * chunk_size;
    } else {
        chunk_size = base_chunk;
        start_index = remainder * (base_chunk + 1) + (rank - remainder) * base_chunk;
    }
    
    return {start_index, chunk_size};
}

} // namespace fe_em