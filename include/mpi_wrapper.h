// MPI面向对象封装头文件
// 提供C++风格的MPI接口封装，屏蔽原生C接口

#ifndef FE_EM_MPI_WRAPPER_H
#define FE_EM_MPI_WRAPPER_H

#include "config.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#if USE_MPI
#include <mpi.h>
#else
// 单进程模式下的MPI类型定义
typedef int MPI_Datatype;
#endif

namespace fe_em {

/**
 * @brief MPI异常类
 */
class MPIException : public std::runtime_error {
private:
    int error_code_;
    int rank_;

public:
    /**
     * @brief 构造函数
     * @param msg 错误信息
     * @param error_code MPI错误码
     * @param rank 发生错误的进程ID
     */
    explicit MPIException(const std::string& msg, int error_code = 0, int rank = -1)
        : std::runtime_error(msg), error_code_(error_code), rank_(rank) {}

    /**
     * @brief 获取MPI错误码
     * @return MPI错误码
     */
    int get_error_code() const { return error_code_; }

    /**
     * @brief 获取发生错误的进程ID
     * @return 进程ID
     */
    int get_rank() const { return rank_; }

    /**
     * @brief 获取完整的错误信息
     * @return 包含错误码和进程ID的完整错误信息
     */
    std::string get_full_message() const {
        std::string msg = what();
        if (error_code_ != 0) {
            msg += " (错误码: " + std::to_string(error_code_) + ")";
        }
        if (rank_ != -1) {
            msg += " (进程ID: " + std::to_string(rank_) + ")";
        }
        return msg;
    }
};

/**
 * @brief MPI通信器封装类（单例）
 */
class MPIComm {
private:
    static std::unique_ptr<MPIComm> instance_;
    bool initialized_;
    bool finalized_;
    int rank_;
    int size_;

    // 私有构造函数
    MPIComm();

public:
    /**
     * @brief 析构函数
     */
    ~MPIComm();

    // 禁止拷贝和移动
    MPIComm(const MPIComm&) = delete;
    MPIComm& operator=(const MPIComm&) = delete;
    MPIComm(MPIComm&&) = delete;
    MPIComm& operator=(MPIComm&&) = delete;

    /**
     * @brief 获取单例实例
     * @return MPIComm单例引用
     */
    static MPIComm& get_instance();

    /**
     * @brief 初始化MPI（如果尚未初始化）
     * @param argc 命令行参数个数指针
     * @param argv 命令行参数数组指针
     * @throws MPIException 如果初始化失败
     */
    void init(int* argc = nullptr, char*** argv = nullptr);

    /**
     * @brief 终止MPI（如果尚未终止）
     * @throws MPIException 如果终止失败
     */
    void finalize();

    /**
     * @brief 获取当前进程ID
     * @return 当前进程ID（0到进程数-1）
     */
    int get_rank() const { return rank_; }

    /**
     * @brief 获取总进程数
     * @return 总进程数
     */
    int get_size() const { return size_; }

    /**
     * @brief 判断当前进程是否为主进程
     * @return true如果是主进程（rank=0），false否则
     */
    bool is_root() const { return rank_ == 0; }

    /**
     * @brief 判断MPI是否已初始化
     * @return true如果已初始化，false否则
     */
    bool is_initialized() const { return initialized_; }

    /**
     * @brief 判断MPI是否已终止
     * @return true如果已终止，false否则
     */
    bool is_finalized() const { return finalized_; }

    /**
     * @brief 屏障同步
     * @throws MPIException 如果同步失败
     */
    void barrier() const;

    /**
     * @brief 点对点发送数据
     * @tparam T 数据类型
     * @param data 要发送的数据
     * @param count 数据个数
     * @param dest_rank 目标进程ID
     * @param tag 消息标签
     * @throws MPIException 如果发送失败
     */
    template<typename T>
    void send(const T* data, int count, int dest_rank, int tag = 0) const;

    /**
     * @brief 点对点接收数据
     * @tparam T 数据类型
     * @param data 接收数据的缓冲区
     * @param count 数据个数
     * @param source_rank 源进程ID
     * @param tag 消息标签
     * @return 实际接收到的数据个数
     * @throws MPIException 如果接收失败
     */
    template<typename T>
    int recv(T* data, int count, int source_rank, int tag = 0) const;

    /**
     * @brief 广播数据（从根进程到所有进程）
     * @tparam T 数据类型
     * @param data 数据缓冲区（根进程发送，所有进程接收）
     * @param count 数据个数
     * @param root_rank 根进程ID
     * @throws MPIException 如果广播失败
     */
    template<typename T>
    void broadcast(T* data, int count, int root_rank = 0) const;

    /**
     * @brief 收集数据（所有进程发送到根进程）
     * @tparam T 数据类型
     * @param send_data 发送数据缓冲区
     * @param send_count 每个进程发送的数据个数
     * @param recv_data 接收数据缓冲区（仅根进程有效）
     * @param recv_counts 每个进程接收的数据个数（仅根进程有效）
     * @param root_rank 根进程ID
     * @throws MPIException 如果收集失败
     */
    template<typename T>
    void gather(const T* send_data, int send_count, 
                T* recv_data, int* recv_counts, int root_rank = 0) const;

    /**
     * @brief 分发数据（根进程分发到所有进程）
     * @tparam T 数据类型
     * @param send_data 发送数据缓冲区（仅根进程有效）
     * @param send_counts 每个进程发送的数据个数（仅根进程有效）
     * @param recv_data 接收数据缓冲区
     * @param recv_count 每个进程接收的数据个数
     * @param root_rank 根进程ID
     * @throws MPIException 如果分发失败
     */
    template<typename T>
    void scatter(const T* send_data, const int* send_counts,
                 T* recv_data, int recv_count, int root_rank = 0) const;

    /**
     * @brief 数据分域辅助函数
     * @param total_size 总数据大小
     * @param rank 当前进程ID
     * @param size 总进程数
     * @return 当前进程负责的数据范围（起始索引，数据个数）
     */
    static std::pair<int, int> split_data(int total_size, int rank, int size);

    /**
     * @brief 数据分域辅助函数（自动使用当前进程信息）
     * @param total_size 总数据大小
     * @return 当前进程负责的数据范围（起始索引，数据个数）
     */
    std::pair<int, int> split_data(int total_size) const {
        return split_data(total_size, rank_, size_);
    }

private:
    /**
     * @brief 检查MPI错误码并抛出异常
     * @param error_code MPI错误码
     * @param operation_name 操作名称
     * @throws MPIException 如果错误码不是MPI_SUCCESS
     */
    void check_error(int error_code, const std::string& operation_name) const;

    /**
     * @brief 获取MPI数据类型
     * @tparam T C++数据类型
     * @return 对应的MPI数据类型
     */
    template<typename T>
    static MPI_Datatype get_mpi_datatype() {
        // 默认实现，需要在特化中定义具体类型
#if USE_MPI
        return MPI_DATATYPE_NULL;
#else
        // 单进程模式下返回一个占位符值
        return static_cast<MPI_Datatype>(0);
#endif
    }
};

// ============================================================================
// 模板函数实现
// ============================================================================

template<typename T>
void MPIComm::send(const T* data, int count, int dest_rank, int tag) const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法发送数据");
    }
    
    if (dest_rank < 0 || dest_rank >= size_) {
        throw MPIException("目标进程ID越界: " + std::to_string(dest_rank));
    }
    
    if (count <= 0) {
        throw MPIException("发送数据个数必须大于0");
    }
    
    MPI_Datatype mpi_type = get_mpi_datatype<T>();
    int error_code = MPI_Send(data, count, mpi_type, dest_rank, tag, MPI_COMM_WORLD);
    check_error(error_code, "MPI_Send");
#else
    // 单进程模式：忽略发送操作
    (void)data;
    (void)count;
    (void)dest_rank;
    (void)tag;
#endif
}

template<typename T>
int MPIComm::recv(T* data, int count, int source_rank, int tag) const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法接收数据");
    }
    
    if (source_rank < 0 || source_rank >= size_) {
        throw MPIException("源进程ID越界: " + std::to_string(source_rank));
    }
    
    if (count <= 0) {
        throw MPIException("接收数据个数必须大于0");
    }
    
    MPI_Datatype mpi_type = get_mpi_datatype<T>();
    MPI_Status status;
    int error_code = MPI_Recv(data, count, mpi_type, source_rank, tag, MPI_COMM_WORLD, &status);
    check_error(error_code, "MPI_Recv");
    
    int actual_count;
    MPI_Get_count(&status, mpi_type, &actual_count);
    return actual_count;
#else
    // 单进程模式：返回0表示未接收到数据
    (void)data;
    (void)count;
    (void)source_rank;
    (void)tag;
    return 0;
#endif
}

template<typename T>
void MPIComm::broadcast(T* data, int count, int root_rank) const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法广播数据");
    }
    
    if (root_rank < 0 || root_rank >= size_) {
        throw MPIException("根进程ID越界: " + std::to_string(root_rank));
    }
    
    if (count <= 0) {
        throw MPIException("广播数据个数必须大于0");
    }
    
    MPI_Datatype mpi_type = get_mpi_datatype<T>();
    int error_code = MPI_Bcast(data, count, mpi_type, root_rank, MPI_COMM_WORLD);
    check_error(error_code, "MPI_Bcast");
#else
    // 单进程模式：数据保持不变
    (void)data;
    (void)count;
    (void)root_rank;
#endif
}

template<typename T>
void MPIComm::gather(const T* send_data, int send_count, 
                     T* recv_data, int* recv_counts, int root_rank) const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法收集数据");
    }
    
    if (root_rank < 0 || root_rank >= size_) {
        throw MPIException("根进程ID越界: " + std::to_string(root_rank));
    }
    
    if (send_count <= 0) {
        throw MPIException("发送数据个数必须大于0");
    }
    
    MPI_Datatype mpi_type = get_mpi_datatype<T>();
    
    // 准备接收计数和位移数组
    std::vector<int> recv_displs(size_);
    if (rank_ == root_rank) {
        int total = 0;
        for (int i = 0; i < size_; ++i) {
            recv_displs[i] = total;
            total += recv_counts[i];
        }
    }
    
    int error_code = MPI_Gatherv(send_data, send_count, mpi_type,
                                 recv_data, recv_counts, recv_displs.data(), mpi_type,
                                 root_rank, MPI_COMM_WORLD);
    check_error(error_code, "MPI_Gatherv");
#else
    // 单进程模式：直接复制数据
    if (root_rank == 0) {
        std::copy(send_data, send_data + send_count, recv_data);
        recv_counts[0] = send_count;
    }
#endif
}

template<typename T>
void MPIComm::scatter(const T* send_data, const int* send_counts,
                      T* recv_data, int recv_count, int root_rank) const {
#if USE_MPI
    if (!initialized_) {
        throw MPIException("MPI未初始化，无法分发数据");
    }
    
    if (root_rank < 0 || root_rank >= size_) {
        throw MPIException("根进程ID越界: " + std::to_string(root_rank));
    }
    
    if (recv_count <= 0) {
        throw MPIException("接收数据个数必须大于0");
    }
    
    MPI_Datatype mpi_type = get_mpi_datatype<T>();
    
    // 准备发送位移数组
    std::vector<int> send_displs(size_);
    if (rank_ == root_rank) {
        int total = 0;
        for (int i = 0; i < size_; ++i) {
            send_displs[i] = total;
            total += send_counts[i];
        }
    }
    
    int error_code = MPI_Scatterv(send_data, send_counts, send_displs.data(), mpi_type,
                                 recv_data, recv_count, mpi_type,
                                 root_rank, MPI_COMM_WORLD);
    check_error(error_code, "MPI_Scatterv");
#else
    // 单进程模式：直接复制数据
    if (root_rank == 0) {
        std::copy(send_data, send_data + recv_count, recv_data);
    }
#endif
}

// MPI数据类型映射
#define DEFINE_MPI_TYPE(CppType, MpiType) \
    template<> \
    inline MPI_Datatype MPIComm::get_mpi_datatype<CppType>() { \
        return MpiType; \
    }

#if USE_MPI
DEFINE_MPI_TYPE(int, MPI_INT)
DEFINE_MPI_TYPE(unsigned int, MPI_UNSIGNED)
DEFINE_MPI_TYPE(long, MPI_LONG)
DEFINE_MPI_TYPE(unsigned long, MPI_UNSIGNED_LONG)
DEFINE_MPI_TYPE(long long, MPI_LONG_LONG)
DEFINE_MPI_TYPE(unsigned long long, MPI_UNSIGNED_LONG_LONG)
DEFINE_MPI_TYPE(float, MPI_FLOAT)
DEFINE_MPI_TYPE(double, MPI_DOUBLE)
DEFINE_MPI_TYPE(char, MPI_CHAR)
DEFINE_MPI_TYPE(unsigned char, MPI_UNSIGNED_CHAR)
DEFINE_MPI_TYPE(bool, MPI_C_BOOL)
#endif

#undef DEFINE_MPI_TYPE

} // namespace fe_em

#endif // FE_EM_MPI_WRAPPER_H