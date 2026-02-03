/**
 * @file fetidp_context.hpp
 * @brief FETI-DP核心层 - 全局上下文管理模块头文件
 * @details 负责FETI-DP求解器的全局状态管理，包括子域信息、对偶变量、收敛控制等
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#pragma once

#include <vector>
#include <memory>

namespace fetidp {

/**
 * @struct SubdomainData
 * @brief 子域数据结构体
 */
struct SubdomainData {
    int id;
    int node_count;
    int element_count;
    int interface_dof_count;
};

/**
 * @class FETIDPContext
 * @brief FETI-DP全局上下文管理类
 */
class FETIDPContext {
public:
    /**
     * @brief 构造函数
     */
    FETIDPContext();

    /**
     * @brief 析构函数
     */
    ~FETIDPContext();

    /**
     * @brief 初始化FETI-DP求解器
     * @param subdomain_count 子域数量
     * @return bool 初始化成功返回true，失败返回false
     */
    bool initialize(int subdomain_count);

    /**
     * @brief 获取子域数量
     * @return int 子域数量
     */
    int getSubdomainCount() const;

    /**
     * @brief 获取子域数据
     * @param subdomain_id 子域ID
     * @return SubdomainData 子域数据
     */
    SubdomainData getSubdomainData(int subdomain_id) const;

private:
    int subdomain_count_;
    std::vector<SubdomainData> subdomains_;
};

} // namespace fetidp