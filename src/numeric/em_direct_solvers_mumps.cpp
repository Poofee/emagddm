/**
 * @file em_direct_solvers_mumps.cpp
 * @brief MUMPS 后端辅助类实现（独立编译单元）
 * @details 将 MUMPS 相关代码从 em_direct_solvers.cpp 拆分出来，
 *          避免非 MUMPS 目标链接时引入 Intel Fortran 运行时依赖。
 *
 * 包含：
 * - MumpsContext RAII 封装类的完整实现
 *
 * @note 原三个直接求解器类已迁移至新架构：
 *       - UnifiedDirectSolver + SolverBackend 策略模式
 *       - MUMPSBackend 统一 MUMPS 后端
 *
 * @see MUMPSBackend 新架构 MUMPS 后端
 * @see MumpsContext RAII 资源管理
 *
 * @author Poofee
 * @date 2026-04-09
 * @version 2.0 (重构版 - 移除旧求解器类)
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"

#ifdef HAVE_MUMPS

namespace {

numeric::SolverResult create_success_result(Eigen::VectorXd x, int iterations = 0) {
    numeric::SolverResult result;
    result.x = std::move(x);
    result.status = numeric::SolverStatus::SUCCESS;
    result.iterations = iterations;
    return result;
}

numeric::SolverResult create_error_result(numeric::SolverStatus status, const std::string& error_msg) {
    numeric::SolverResult result;
    result.status = status;
    result.error_msg = error_msg;
    return result;
}

} // anonymous namespace

namespace numeric {

// ============================================================================
// MumpsContext 实现（MUMPS RAII 封装）
// ============================================================================

/**
 * @brief 构造函数，初始化为未初始化状态
 */
MumpsContext::MumpsContext() : state_(State::UNINITIALIZED), n_(0) {
    memset(&mumps_data_, 0, sizeof(mumps_data_));
}

/**
 * @brief 析构函数，自动释放 MUMPS 资源
 */
MumpsContext::~MumpsContext() {
    reset();
}

/**
 * @brief 移动构造函数
 */
MumpsContext::MumpsContext(MumpsContext&& other) noexcept
    : state_(other.state_), n_(other.n_),
      mumps_data_(other.mumps_data_),
      rhs_buffer_(std::move(other.rhs_buffer_)),
      irn_storage_(std::move(other.irn_storage_)),
      jcn_storage_(std::move(other.jcn_storage_)),
      a_storage_(std::move(other.a_storage_))
{
    other.state_ = State::UNINITIALIZED;
    other.n_ = 0;
    memset(&other.mumps_data_, 0, sizeof(other.mumps_data_));
}

/**
 * @brief 移动赋值运算符
 */
MumpsContext& MumpsContext::operator=(MumpsContext&& other) noexcept {
    if (this != &other) {
        reset();
        state_ = other.state_;
        n_ = other.n_;
        mumps_data_ = other.mumps_data_;
        rhs_buffer_ = std::move(other.rhs_buffer_);
        irn_storage_ = std::move(other.irn_storage_);
        jcn_storage_ = std::move(other.jcn_storage_);
        a_storage_ = std::move(other.a_storage_);
        other.state_ = State::UNINITIALIZED;
        other.n_ = 0;
        memset(&other.mumps_data_, 0, sizeof(other.mumps_data_));
    }
    return *this;
}

/**
 * @brief 初始化 MUMPS 实例
 * @param sym 对称性标志：0=非对称，1=正定对称，2=一般对称
 * @return true 初始化成功，false 初始化失败
 */
bool MumpsContext::initialize(int sym) {
    if (state_ != State::UNINITIALIZED) {
        reset();
    }

    memset(&mumps_data_, 0, sizeof(mumps_data_));

    mumps_data_.par = 1;
    mumps_data_.sym = sym;
    mumps_data_.job = -1;

    dmumps_c(&mumps_data_);

    if (mumps_data_.info[0] < 0) {
        FEEM_ERROR("MumpsContext::initialize - MUMPS 初始化失败: info[0]={}, info[1]={}",
                   mumps_data_.info[0], mumps_data_.info[1]);
        memset(&mumps_data_, 0, sizeof(mumps_data_));
        return false;
    }

    configure_default_icntl();

    state_ = State::INITIALIZED;
    FEEM_DEBUG("MumpsContext::initialize 成功, sym={}", sym);
    return true;
}

/**
 * @brief 执行矩阵分析+分解（JOB=1 分析 + JOB=2 分解）
 * @param n 矩阵维度
 * @param nz 非零元个数
 * @param irn 行索引数组（1-based，MUMPS 格式）
 * @param jcn 列索引数组（1-based，MUMPS 格式）
 * @param a 非零值数组
 * @return 0 成功，非零 MUMPS 错误码
 *
 * @note Windows+Intel Fortran 下 JOB=5 存在 NNZ 传递 bug，必须拆分为 JOB=1+JOB=2
 */
int MumpsContext::factorize(int n, int nz, int* irn, int* jcn, double* a) {
    if (state_ != State::INITIALIZED) {
        FEEM_ERROR("MumpsContext::factorize - 状态错误（需要先调用 initialize）");
        return -1;
    }

    n_ = n;

    mumps_data_.n = n;
    mumps_data_.nz = nz;
    mumps_data_.nnz = static_cast<MUMPS_INT8>(nz);
    mumps_data_.irn = irn;
    mumps_data_.jcn = jcn;
    mumps_data_.a = reinterpret_cast<DMUMPS_COMPLEX*>(a);

    mumps_data_.job = 1;
    dmumps_c(&mumps_data_);

    if (mumps_data_.info[0] < 0) {
        int err = mumps_data_.info[0];
        FEEM_ERROR("MumpsContext::factorize - MUMPS 分析阶段失败: info[0]={}, info[1]={}",
                   err, mumps_data_.info[1]);
        state_ = State::INITIALIZED;
        return err;
    }

    mumps_data_.job = 2;
    dmumps_c(&mumps_data_);

    if (mumps_data_.info[0] < 0) {
        int err = mumps_data_.info[0];
        FEEM_ERROR("MumpsContext::factorize - MUMPS 分析/分解失败: info[0]={}, info[1]={}, infog[1]={}, infog[2]={}",
                   err, mumps_data_.info[1], mumps_data_.infog[0], mumps_data_.infog[1]);
        state_ = State::INITIALIZED;
        return err;
    }

    rhs_buffer_.resize(n);

    state_ = State::FACTORIZED;
    FEEM_DEBUG("MumpsContext::factorize 成功, n={}, nz={}", n, nz);
    return 0;
}

/**
 * @brief 使用 MUMPS 执行求解（JOB=3）
 * @param b 右端项向量
 * @return SolverResult 求解结果
 */
SolverResult MumpsContext::solve(const Eigen::VectorXd& b) {
    if (state_ != State::FACTORIZED) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "MUMPS 未完成分解，请先调用 factorize()");
    }
    if (b.size() != n_) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量维度不匹配");
    }

    for (int i = 0; i < n_; ++i) {
        rhs_buffer_[i] = b(i);
    }

    mumps_data_.rhs = rhs_buffer_.data();
    mumps_data_.job = 3;
    dmumps_c(&mumps_data_);

    if (mumps_data_.info[0] < 0) {
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "MUMPS 求解失败 (info[0]=" +
                                   std::to_string(mumps_data_.info[0]) + ")");
    }

    Eigen::VectorXd x(n_);
    for (int i = 0; i < n_; ++i) {
        x(i) = rhs_buffer_[i];
    }

    return create_success_result(std::move(x));
}

/**
 * @brief 重置 MUMPS 实例，释放所有资源
 */
void MumpsContext::reset() {
    if (state_ != State::UNINITIALIZED) {
        mumps_data_.job = -2;
        dmumps_c(&mumps_data_);
    }

    memset(&mumps_data_, 0, sizeof(mumps_data_));
    rhs_buffer_.clear();
    irn_storage_.clear();
    jcn_storage_.clear();
    a_storage_.clear();
    state_ = State::UNINITIALIZED;
    n_ = 0;
}

/**
 * @brief 配置 MUMPS 默认 ICNTL 参数
 * @details 关闭所有输出，使用集中式输入、PORD 排序、自动符号分解
 */
void MumpsContext::configure_default_icntl() {
    mumps_data_.icntl[0] = 6;
    mumps_data_.icntl[1] = 0;
    mumps_data_.icntl[2] = 0;
    mumps_data_.icntl[3] = 0;

    mumps_data_.icntl[4] = 0;
    mumps_data_.icntl[6] = 7;
    mumps_data_.icntl[12] = 1;
    mumps_data_.icntl[17] = 0;
    mumps_data_.icntl[21] = 0;
    mumps_data_.icntl[27] = 0;
    mumps_data_.icntl[57] = 0;
}

/**
 * @brief 从 CSR 矩阵设置数据并执行分解（便捷方法）
 * @param csr CSR 格式稀疏矩阵（实数版本）
 * @param sym 对称性标志：0=不对称，1=对称正定，2=一般对称
 * @return 0 成功；<0 错误码
 *
 * @details 自动完成以下操作：
 * 1. 若未初始化则自动调用 initialize(sym)
 * 2. 将 CSR 格式转换为 COO 格式存储到内部数组（irn_storage_, jcn_storage_, a_storage_）
 * 3. 调用 factorize() 执行分析+分解
 */
int MumpsContext::factorize_from_csr(const CsrMatrix<double>& csr, int sym) {
    // 自动初始化（若尚未初始化）
    if (state_ == State::UNINITIALIZED) {
        if (!initialize(sym)) {
            FEEM_ERROR("MumpsContext::factorize_from_csr - 自动初始化失败");
            return -1;
        }
    }

    // 提取 CSR 矩阵数据
    const int n = csr.rows();
    const auto& row_ptr = csr.get_row_ptr();
    const auto& col_indices = csr.get_col_indices();
    const auto& values = csr.get_values();

    // CSR → COO 转换（MUMPS 使用 1-based 索引的 COO 格式）
    irn_storage_.clear();
    jcn_storage_.clear();
    a_storage_.clear();

    const int nz = row_ptr[n]; // 从 CSR 的 row_ptr 获取非零元总数
    irn_storage_.reserve(nz);
    jcn_storage_.reserve(nz);
    a_storage_.reserve(nz);

    for (int i = 0; i < n; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            irn_storage_.push_back(i + 1);           // 1-based 行索引
            jcn_storage_.push_back(col_indices[k] + 1); // 1-based 列索引
            a_storage_.push_back(values[k]);
        }
    }

    // 调用 factorize()（使用内部存储的数组指针）
    return factorize(n, nz, irn_storage_.data(), jcn_storage_.data(), a_storage_.data());
}

} // namespace numeric

/**
 * @brief Fortran 主程序入口点 stub
 * @details 当 C++ 程序链接 Intel Fortran 运行时库 (libifcoremd.lib) 时，
 *          for_main.obj 中的 main() 会引用 MAIN__ 符号。此 stub 满足链接器需求。
 *          由于 C++ 对象文件中的 main() 优先于库文件中的 main()，
 *          此函数永远不会被实际调用。
 */
extern "C" void MAIN__() {}

#endif // HAVE_MUMPS
