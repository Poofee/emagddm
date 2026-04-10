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

#ifdef HAVE_MUMPS

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
// MUMPS复数求解器通用辅助函数（消除三个求解器的重复代码）
// ============================================================================

/**
 * @brief MUMPS复数分解通用实现
 * @param eigen_matrix 复数Eigen稀疏矩阵引用
 * @param mumps_ctx 输出参数，用于存储创建的ComplexMumpsContext
 * @param sym 对称性标志：0=非对称, 1=Hermitian正定, 2=一般Hermitian
 * @param solver_name 求解器名称（用于日志）
 * @return SolverResult 分解结果
 *
 * @details 统一的MUMPS复数分解流程：
 *          1. 压缩Eigen矩阵
 *          2. 转换为COO格式
 *          3. 构建CSR矩阵
 *          4. 创建ComplexMumpsContext并执行分解
 */
SolverResult decompose_mumps_complex_helper(
    const Eigen::SparseMatrix<std::complex<double>>& eigen_matrix,
    std::unique_ptr<ComplexMumpsContext>& mumps_ctx,
    int sym,
    const std::string& solver_name) {

    FEEM_DEBUG("{}::decompose_mumps_complex_helper - 开始MUMPS复数分解 (sym={})", solver_name, sym);

    try {
        const int n = static_cast<int>(eigen_matrix.rows());
        const_cast<Eigen::SparseMatrix<std::complex<double>>&>(eigen_matrix).makeCompressed();

        // 将 Eigen 复数矩阵转换为 CsrMatrix
        CooMatrix<std::complex<double>> coo(n, n);
        for (int col = 0; col < n; ++col) {
            for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(eigen_matrix, col); it; ++it) {
                coo.add_value(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
            }
        }

        CsrMatrix<std::complex<double>> csr(n, n);
        csr.build_from_coo(coo);

        mumps_ctx = std::make_unique<ComplexMumpsContext>();

        // 使用便捷方法：自动处理初始化、CSR→COO转换、复数分解
        int info = mumps_ctx->factorize_from_csr(csr, sym);
        if (info != 0) {
            mumps_ctx.reset();
            return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                       "MUMPS复数分解失败 (info=" + std::to_string(info) + ")");
        }

        FEEM_DEBUG("{}::decompose_mumps_complex_helper - MUMPS复数分解成功, 维度: {}x{}", solver_name, n, n);
        return create_success_result(Eigen::VectorXd());  // 空向量占位（分解无解向量）
    } catch (const std::bad_alloc& e) {
        FEEM_ERROR("{}::decompose_mumps_complex_helper - 内存分配失败: {}", solver_name, e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   "内存不足，无法完成MUMPS复数分解");
    } catch (const std::exception& e) {
        FEEM_ERROR("{}::decompose_mumps_complex_helper - 异常: {}", solver_name, e.what());
        return create_error_result(SolverStatus::NUMERICAL_ERROR,
                                   std::string("MUMPS复数分解异常: ") + e.what());
    }
}

/**
 * @brief MUMPS复数求解通用实现
 * @param mumps_ctx ComplexMumpsContext智能指针引用
 * @param b 复数右端项向量
 * @param solver_name 求解器名称（用于日志）
 * @return SolverResult 求解结果（x_complex字段存储复数解）
 *
 * @details 统一的MUMPS复数求解流程：
 *          1. 验证分解状态
 *          2. 调用ComplexMumpsContext::solve()
 *          3. 将交错实数向量转换为复数向量（使用头文件中已定义的extract_complex）
 */
SolverResult solve_mumps_complex_helper(
    const std::unique_ptr<ComplexMumpsContext>& mumps_ctx,
    const Eigen::VectorXcd& b,
    const std::string& solver_name) {

    FEEM_DEBUG("{}::solve_mumps_complex_helper - 开始MUMPS复数求解", solver_name);

    if (!mumps_ctx || !mumps_ctx->is_factored()) {
        return create_error_result(SolverStatus::INVALID_INPUT,
                                   "MUMPS复数未完成分解，请先调用set_matrix()");
    }

    // ComplexMumpsContext::solve() 返回的 result.x 为实数向量（长度 2n），按 [Re,Im,...] 排列
    SolverResult raw_result = mumps_ctx->solve(b);

    if (raw_result.status != SolverStatus::SUCCESS) {
        return raw_result;
    }

    // 将交错排列的实数向量转换为 VectorXcd（使用头文件中定义的extract_complex）
    SolverResult result;
    result.x_complex = extract_complex(raw_result);
    result.status = SolverStatus::SUCCESS;
    result.iterations = raw_result.iterations;

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
      rhs_buffer_(std::move(other.rhs_buffer_)),
      irn_storage_(std::move(other.irn_storage_)),
      jcn_storage_(std::move(other.jcn_storage_)),
      a_storage_(std::move(other.a_storage_))
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
        irn_storage_ = std::move(other.irn_storage_);
        jcn_storage_ = std::move(other.jcn_storage_);
        a_storage_ = std::move(other.a_storage_);
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
    irn_storage_.clear();
    jcn_storage_.clear();
    a_storage_.clear();
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

/**
 * @brief 从 CSR 矩阵设置数据并执行复数分解（便捷方法）
 * @param csr CSR 格式稀疏矩阵（复数版本）
 * @param sym 对称性标志：0=不对称，1=Hermitian正定，2=一般Hermitian
 * @return 0 成功；<0 错误码
 *
 * @details 自动完成以下操作：
 * 1. 若未初始化则自动调用 initialize(sym)
 * 2. 将 CSR 格式转换为 COO 格式存储到内部数组
 * 3. 调用 factorize() 执行分析+分解
 */
int ComplexMumpsContext::factorize_from_csr(const CsrMatrix<std::complex<double>>& csr, int sym) {
    // 自动初始化（若尚未初始化）
    if (state_ == State::UNINITIALIZED) {
        if (!initialize(sym)) {
            FEEM_ERROR("ComplexMumpsContext::factorize_from_csr - 自动初始化失败");
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

// ============================================================================
// 三个求解器的 MUMPS 复数后端方法实现
// ============================================================================

namespace numeric {

// ============================================================================
// SymmetricDirectSolver MUMPS 复数后端
// ============================================================================

/**
 * @brief 使用MUMPS复数后端执行分解（Hermitian正定并行模式）
 * @return SolverResult 分解结果
 *
 * @details 委托给通用辅助函数decompose_mumps_complex_helper，sym=1: Hermitian正定
 */
SolverResult SymmetricDirectSolver::decompose_with_mumps_complex() {
    return decompose_mumps_complex_helper(eigen_matrix_complex_, mumps_ctx_complex_, 1,
                                         "SymmetricDirectSolver");
}

/**
 * @brief 使用MUMPS复数后端执行求解
 * @param b 复数右端项向量
 * @return SolverResult 求解结果（复数解存储在 x_complex 字段）
 *
 * @details 委托给通用辅助函数solve_mumps_complex_helper
 */
SolverResult SymmetricDirectSolver::solve_with_mumps_complex(const Eigen::VectorXcd& b) {
    return solve_mumps_complex_helper(mumps_ctx_complex_, b, "SymmetricDirectSolver");
}

// ============================================================================
// SymmetricIndefiniteDirectSolver MUMPS 复数后端
// ============================================================================

/**
 * @brief 使用MUMPS复数后端执行分解（一般Hermitian并行模式）
 * @return SolverResult 分解结果
 *
 * @details 委托给通用辅助函数decompose_mumps_complex_helper，sym=2: 一般Hermitian
 */
SolverResult SymmetricIndefiniteDirectSolver::decompose_with_mumps_complex() {
    return decompose_mumps_complex_helper(eigen_matrix_complex_, mumps_ctx_complex_, 2,
                                         "SymmetricIndefiniteDirectSolver");
}

/**
 * @brief 使用MUMPS复数后端执行求解
 * @param b 复数右端项向量
 * @return SolverResult 求解结果（复数解存储在 x_complex 字段）
 *
 * @details 委托给通用辅助函数solve_mumps_complex_helper
 */
SolverResult SymmetricIndefiniteDirectSolver::solve_with_mumps_complex(const Eigen::VectorXcd& b) {
    return solve_mumps_complex_helper(mumps_ctx_complex_, b, "SymmetricIndefiniteDirectSolver");
}

// ============================================================================
// GeneralDirectSolver MUMPS 复数后端
// ============================================================================

/**
 * @brief 使用MUMPS复数后端执行分解（非对称并行模式）
 * @return SolverResult 分解结果
 *
 * @details 委托给通用辅助函数decompose_mumps_complex_helper，sym=0: 非对称
 */
SolverResult GeneralDirectSolver::decompose_with_mumps_complex() {
    return decompose_mumps_complex_helper(eigen_matrix_complex_, mumps_ctx_complex_, 0,
                                         "GeneralDirectSolver");
}

/**
 * @brief 使用MUMPS复数后端执行求解
 * @param b 复数右端项向量
 * @return SolverResult 求解结果（复数解存储在 x_complex 字段）
 *
 * @details 委托给通用辅助函数solve_mumps_complex_helper
 */
SolverResult GeneralDirectSolver::solve_with_mumps_complex(const Eigen::VectorXcd& b) {
    return solve_mumps_complex_helper(mumps_ctx_complex_, b, "GeneralDirectSolver");
}

} // namespace numeric

#endif // HAVE_MUMPS
