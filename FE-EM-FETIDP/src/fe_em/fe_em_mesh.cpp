/**
 * @file fe_em_mesh.cpp
 * @brief 电磁物理层 - 电磁网格管理模块源文件
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include "fe_em_mesh.hpp"

namespace fe_em {

EMMesh::EMMesh() : node_count_(0), element_count_(0) {
}

EMMesh::~EMMesh() {
}

bool EMMesh::loadFromGmsh(const std::string& mesh_file) {
    // TODO: 实现Gmsh网格文件加载逻辑
    return false;
}

int EMMesh::getNodeCount() const {
    return node_count_;
}

int EMMesh::getElementCount() const {
    return element_count_;
}

} // namespace fe_em