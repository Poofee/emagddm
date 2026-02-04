/**
 * @file fetidp_context.cpp
 * @brief FETI-DP核心层 - 全局上下文管理模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "fetidp_context.hpp"

namespace fetidp {

FETIDPContext::FETIDPContext() : subdomain_count_(0) {
}

FETIDPContext::~FETIDPContext() {
}

bool FETIDPContext::initialize(int subdomain_count) {
    if (subdomain_count <= 0) {
        return false;
    }
    
    subdomain_count_ = subdomain_count;
    subdomains_.resize(subdomain_count);
    
    for (int i = 0; i < subdomain_count; ++i) {
        subdomains_[i].id = i;
        subdomains_[i].node_count = 0;
        subdomains_[i].element_count = 0;
        subdomains_[i].interface_dof_count = 0;
    }
    
    return true;
}

int FETIDPContext::getSubdomainCount() const {
    return subdomain_count_;
}

SubdomainData FETIDPContext::getSubdomainData(int subdomain_id) const {
    if (subdomain_id < 0 || subdomain_id >= subdomain_count_) {
        return SubdomainData{};
    }
    return subdomains_[subdomain_id];
}

} // namespace fetidp