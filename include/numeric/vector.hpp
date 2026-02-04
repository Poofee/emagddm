#pragma once

#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace emag {

/**
 * @brief 向量数据类型枚举
 */
enum class VectorDataType {
    REAL,      ///< 实数向量
    COMPLEX    ///< 复数向量
};

/**
 * @brief 向量类模板 - 支持实数和复数向量
 * @tparam T 数据类型 (double 或 std::complex<double>)
 */
template<typename T>
class Vector {
public:
    /**
     * @brief 默认构造函数
     */
    Vector() : size_(0) {}
    
    /**
     * @brief 构造函数，指定向量大小
     * @param size 向量大小
     */
    explicit Vector(int size) : size_(size), data_(size, T(0)) {}
    
    /**
     * @brief 构造函数，从初始值列表初始化
     * @param init_list 初始值列表
     */
    Vector(std::initializer_list<T> init_list) 
        : size_(static_cast<int>(init_list.size())), data_(init_list) {}
    
    /**
     * @brief 拷贝构造函数
     * @param other 另一个向量
     */
    Vector(const Vector& other) = default;
    
    /**
     * @brief 移动构造函数
     * @param other 另一个向量
     */
    Vector(Vector&& other) noexcept = default;
    
    /**
     * @brief 析构函数
     */
    ~Vector() = default;
    
    /**
     * @brief 拷贝赋值运算符
     * @param other 另一个向量
     * @return Vector& 当前向量引用
     */
    Vector& operator=(const Vector& other) = default;
    
    /**
     * @brief 移动赋值运算符
     * @param other 另一个向量
     * @return Vector& 当前向量引用
     */
    Vector& operator=(Vector&& other) noexcept = default;
    
    /**
     * @brief 获取向量大小
     * @return int 向量大小
     */
    int size() const { return size_; }
    
    /**
     * @brief 获取向量数据类型
     * @return VectorDataType 数据类型
     */
    VectorDataType get_data_type() const {
        if constexpr (std::is_same_v<T, double>) {
            return VectorDataType::REAL;
        } else {
            return VectorDataType::COMPLEX;
        }
    }
    
    /**
     * @brief 获取向量数据
     * @return const std::vector<T>& 向量数据引用
     */
    const std::vector<T>& get_data() const { return data_; }
    
    /**
     * @brief 获取向量数据（可修改）
     * @return std::vector<T>& 向量数据引用
     */
    std::vector<T>& get_data() { return data_; }
    
    /**
     * @brief 设置向量大小
     * @param size 新的向量大小
     */
    void resize(int size) {
        size_ = size;
        data_.resize(size, T(0));
    }
    
    /**
     * @brief 清空向量
     */
    void clear() {
        size_ = 0;
        data_.clear();
    }
    
    /**
     * @brief 向量元素访问（只读）
     * @param index 元素索引
     * @return const T& 元素值
     */
    const T& operator[](int index) const {
        if (index < 0 || index >= size_) {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[index];
    }
    
    /**
     * @brief 向量元素访问（可修改）
     * @param index 元素索引
     * @return T& 元素值引用
     */
    T& operator[](int index) {
        if (index < 0 || index >= size_) {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[index];
    }
    
    /**
     * @brief 向量加法
     * @param other 另一个向量
     * @return Vector 结果向量
     */
    Vector operator+(const Vector& other) const {
        if (size_ != other.size_) {
            throw std::invalid_argument("Vector sizes do not match for addition");
        }
        Vector result(size_);
        for (int i = 0; i < size_; ++i) {
            result[i] = data_[i] + other.data_[i];
        }
        return result;
    }
    
    /**
     * @brief 向量减法
     * @param other 另一个向量
     * @return Vector 结果向量
     */
    Vector operator-(const Vector& other) const {
        if (size_ != other.size_) {
            throw std::invalid_argument("Vector sizes do not match for subtraction");
        }
        Vector result(size_);
        for (int i = 0; i < size_; ++i) {
            result[i] = data_[i] - other.data_[i];
        }
        return result;
    }
    
    /**
     * @brief 向量数乘
     * @param scalar 标量
     * @return Vector 结果向量
     */
    Vector operator*(T scalar) const {
        Vector result(size_);
        for (int i = 0; i < size_; ++i) {
            result[i] = data_[i] * scalar;
        }
        return result;
    }
    
    /**
     * @brief 向量点积
     * @param other 另一个向量
     * @return T 点积结果
     */
    T dot(const Vector& other) const {
        if (size_ != other.size_) {
            throw std::invalid_argument("Vector sizes do not match for dot product");
        }
        T result = T(0);
        for (int i = 0; i < size_; ++i) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                result += std::conj(data_[i]) * other.data_[i];
            } else {
                result += data_[i] * other.data_[i];
            }
        }
        return result;
    }
    
    /**
     * @brief 向量2-范数
     * @return double 2-范数值
     */
    double norm() const {
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            return std::sqrt(std::real(dot(*this)));
        } else {
            return std::sqrt(dot(*this));
        }
    }
    
    /**
     * @brief 向量填充
     * @param value 填充值
     */
    void fill(T value) {
        std::fill(data_.begin(), data_.end(), value);
    }
    
    /**
     * @brief 向量置零
     */
    void set_zero() {
        fill(T(0));
    }
    
    /**
     * @brief 打印向量信息
     */
    void print_info() const {
        std::cout << "Vector size: " << size_ << std::endl;
        std::cout << "Data type: " << (get_data_type() == VectorDataType::REAL ? "REAL" : "COMPLEX") << std::endl;
        std::cout << "Norm: " << norm() << std::endl;
    }
    
    /**
     * @brief 打印向量内容
     * @param max_elements 最大显示元素数（默认显示全部）
     */
    void print(int max_elements = -1) const {
        int display_count = (max_elements > 0) ? std::min(max_elements, size_) : size_;
        std::cout << "Vector [" << size_ << "]: ";
        for (int i = 0; i < display_count; ++i) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                std::cout << "(" << std::real(data_[i]) << "+" << std::imag(data_[i]) << "i)";
            } else {
                std::cout << data_[i];
            }
            if (i < display_count - 1) std::cout << ", ";
        }
        if (display_count < size_) std::cout << "...";
        std::cout << std::endl;
    }

private:
    int size_;              ///< 向量大小
    std::vector<T> data_;   ///< 向量数据
};

// 常用类型别名
using VectorReal = Vector<double>;
using VectorComplex = Vector<std::complex<double>>;

} // namespace emag