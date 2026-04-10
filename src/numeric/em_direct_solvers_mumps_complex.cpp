/**
 * @file em_direct_solvers_mumps_complex.cpp
 * @brief MUMPS 复数后端直接求解器实现（独立编译单元）
 * @details 封装 ZMUMPS_STRUC_C 和 zmumps_c C API，
 *          提供 ComplexMumpsContext RAII 类的完整实现。
 *
 * 复数类型映射：
 * - C++ 侧：std::complex<double>（16 字节，实部+虚部连续存储）
 * - MUMPS 侧：ZMUMPS_COMPLEX = mumps_double_complex {double r, double i}
 * - C++11 标准保证两者内存布局兼容，可通过 reinterpret_cast 互转
 *
 * @author Poofee
 * @date 2026-04-09
 */

#include "em_direct_solvers.h"
#include "logger_factory.hpp"

#if EM_SOLVER_HAS_MUMPS

#include <complex>
#include <cstring>

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
// ComplexMumpsContext 实现
// ============================================================================

ComplexMumpsContext::ComplexMumpsContext()
    : state_(State::UNINITIALIZED), n_(0) {
    memset(&zmumps_data_, 0, sizeof(zmumps_data_));
}

ComplexMumpsContext::~ComplexMumpsContext() {
    reset();
}

ComplexMumpsContext::ComplexMumpsContext(ComplexMumpsContext&& other) noexcept
    : state_(other.state_), n_(other.n_),
      zmumps_data_(other.zmumps_data_),
      rhs_buffer_(std::move(other.rhs_buffer_))
{
    other.state_ = State::UNINITIALIZED;
    other.n_ = 0;
    memset(&other.zmumps_data_, 0, sizeof(other.zmumps_data_));
}

ComplexMumpsContext& ComplexMumpsContext::operator=(ComplexMumpsContext&& other) noexcept {
    if (this != &other) {
        reset();
        state_ = other.state_;
        n_ = other.n_;
        zmumps_data_ = other.zmumps_data_;
        rhs_buffer_ = std::move(other.rhs_buffer_);
        other.state_ = State::UNINITIALIZED;
        other.n_ = 0;
        memset(&other.zmumps_data_, 0, sizeof(other.zmumps_data_));
    }
    return *this;
}

/**
 * @brief 初始化 MUMPS 复数实例（job=-1）
 * @param sym 对称性标志：0=非对称，1=Hermitian正定，2=一般Hermitian
 * @return true 初始化成功
 */
bool ComplexMumpsContext::initialize(int sym) {
    if (state_ != State::UNINITIALIZED) {
        reset();
    }

    memset(&zmumps_data_, 0, sizeof(zmumps_data_));

    zmumps_data_.par = 1;
    zmumps_data_.sym = sym;
    zmumps_data_.job = -1;

    zmumps_c(&zmumps_data_);

    if (zmumps_data_.info[0] < 0) {
        FEEM_ERROR("ComplexMumpsContext::initialize - MUMPS 复数初始化失败: info[0]={}, info[1]={}",
                   zmumps_data_.info[0], zmumps_data_.info[1]);
        memset(&zmumps_data_, 0, sizeof(zmumps_data_));
        return false;
    }

    configure_default_icntl();

    state_ = State::INITIALIZED;
    FEEM_DEBUG("ComplexMumpsContext::initialize 成功, sym={}", sym);
    return true;
}

/**
 * @brief 执行复数矩阵分析+分解（job=1 + job=2）
 * @param n 矩阵维度
 * @param nz 非零元个数
 * @param irn 行索引数组（1-based）
 * @param jcn 列索引数组（1-based）
 * @param a 复数非零值数组（std::complex<double> 格式）
 * @return 0 成功，非零 MUMPS 错误码
 * @note Windows+Intel Fortran 下 JOB=5 存在 NNZ 传递 bug，使用 JOB=1+JOB=2 替代
 */
int ComplexMumpsContext::factorize(int n, int nz, int* irn, int* jcn, std::complex<double>* a) {
    if (state_ != State::INITIALIZED) {
        FEEM_ERROR("ComplexMumpsContext::factorize - 状态错误（需要先调用 initialize）");
        return -1;
    }

    n_ = n;

    zmumps_data_.n = n;
    zmumps_data_.nz = nz;
    zmumps_data_.nnz = static_cast<MUMPS_INT8>(nz);
    zmumps_data_.irn = irn;
    zmumps_data_.jcn = jcn;
    zmumps_data_.a = reinterpret_cast<ZMUMPS_COMPLEX*>(a);

    zmumps_data_.job = 1;
    zmumps_c(&zmumps_data_);

    if (zmumps_data_.info[0] < 0) {
        int err = zmumps_data_.info[0];
        FEEM_ERROR("ComplexMumpsContext::factorize - MUMPS 复数分析阶段失败: info[0]={}, info[1]={}",
                   err, zmumps_data_.info[1]);
        state_ = State::INITIALIZED;
        return err;
    }

    zmumps_data_.job = 2;
    zmumps_c(&zmumps_data_);

    if (zmumps_data_.info[0] < 0) {
        int err = zmumps_data_.info[0];
        FEEM_ERROR("ComplexMumpsContext::factorize - MUMPS 复数分解失败: info[0]={}, info[1]={}",
                   err, zmumps_data_.info[1]);
        state_ = State::INITIALIZED;
        return err;
    }

    rhs_buffer_.resize(n);

    state_ = State::FACTORIZED;
    FEEM_DEBUG("ComplexMumpsContext::factorize 成功, n={}, nz={}", n, nz);
    return 0;
}

/**
 * @brief 使用 MUMPS 执行复数求解（job=3）
 * @param b 复数右端项向量
 * @return SolverResult 求解结果（result.x 为实数向量，长度 2n，
 *         按 [Re(x0), Im(x0), Re(x1), Im(x1), ...] 交错排列）
 */
SolverResult ComplexMumpsContext::solve(const Eigen::VectorXcd& b) {
    if (state_ != State::FACTORIZED) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "MUMPS 复数求解器未完成分解，请先调用 factorize()");
    }
    if (b.size() != n_) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "右端项向量维度不匹配");
    }

    for (int i = 0; i < n_; ++i) {
        rhs_buffer_[i].r = b(i).real();
        rhs_buffer_[i].i = b(i).imag();
    }

    zmumps_data_.rhs = rhs_buffer_.data();
    zmumps_data_.job = 3;
    zmumps_c(&zmumps_data_);

    if (zmumps_data_.info[0] < 0) {
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "MUMPS 复数求解失败 (info[0]=" +
                                   std::to_string(zmumps_data_.info[0]) + ")");
    }

    Eigen::VectorXd x(2 * n_);
    for (int i = 0; i < n_; ++i) {
        x(2 * i) = rhs_buffer_[i].r;
        x(2 * i + 1) = rhs_buffer_[i].i;
    }

    return create_success_result(std::move(x));
}

/**
 * @brief 重置 MUMPS 复数实例，释放所有资源
 */
void ComplexMumpsContext::reset() {
    if (state_ != State::UNINITIALIZED) {
        zmumps_data_.job = -2;
        zmumps_c(&zmumps_data_);
    }

    memset(&zmumps_data_, 0, sizeof(zmumps_data_));
    rhs_buffer_.clear();
    state_ = State::UNINITIALIZED;
    n_ = 0;
}

/**
 * @brief 配置 MUMPS 复数默认 ICNTL 参数
 * @details 关闭所有输出，使用集中式输入、PORD 排序、自动符号分解
 */
void ComplexMumpsContext::configure_default_icntl() {
    zmumps_data_.icntl[0] = 6;
    zmumps_data_.icntl[1] = 0;
    zmumps_data_.icntl[2] = 0;
    zmumps_data_.icntl[3] = 0;

    zmumps_data_.icntl[4] = 0;
    zmumps_data_.icntl[6] = 7;
    zmumps_data_.icntl[12] = 1;
    zmumps_data_.icntl[17] = 0;
    zmumps_data_.icntl[21] = 0;
    zmumps_data_.icntl[27] = 0;
    zmumps_data_.icntl[57] = 0;
}

} // namespace numeric

#endif // EM_SOLVER_HAS_MUMPS
