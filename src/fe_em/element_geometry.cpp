/**
 * @file element_geometry.cpp
 * @brief 电磁物理层 - 有限元单元几何定义工具类实现
 * @details 硬编码标准有限元单元的局部棱边和面定义。
 *          遵循ANSYS Fluent/VTK标准顺序。
 * @author Poofee
 * @date 2026-04-05
 * @version 2.0 (重构为枚举类型索引)
 */

#include "element_geometry.hpp"
#include <unordered_map>
#include <algorithm>

namespace fe_em {

// ==================== 内部辅助：单元几何定义注册表 ====================

/**
 * @struct ElementGeometryDef
 * @brief 单元几何定义内部数据结构
 */
struct ElementGeometryDef {
    int num_nodes;                              ///< 节点数
    int num_edges;                              ///< 棱边数
    int num_faces;                              ///< 面数
    std::vector<std::tuple<int, int>> edges;     ///< 局部棱边定义（升序）
    std::vector<std::vector<int>> faces;        ///< 局部面定义（逆时针/右手定则）
};

/**
 * @brief 初始化所有支持的单元类型几何定义
 * @return std::unordered_map<ElemType, ElementGeometryDef> 几何定义映射表（枚举索引，O(1)查询）
 */
static std::unordered_map<ElemType, ElementGeometryDef> initialize_geometry_definitions() {
    std::unordered_map<ElemType, ElementGeometryDef> defs;
    
    // ==================== 二维单元 ====================
    
    // TRI3: 3节点线性三角形
    // 标准来源: VTK_TRIANGLE / ANSYS PLANE13
    defs[ElemType::TRI3] = {
        3,  // 节点数
        3,  // 棱边数
        1,  // 面数（二维单元只有1个面，即单元本身）
        {   // 棱边定义：(0,1), (1,2), (0,2)
            {0, 1},
            {1, 2},
            {0, 2}
        },
        {   // 面定义: [0, 1, 2] (逆时针)
            {0, 1, 2}
        }
    };
    
    // TRI6: 6节点二次三角形（顶点 + 边中点）
    // 标准来源: VTK_QUADRATIC_TRIANGLE
    defs[ElemType::TRI6] = {
        6,  // 节点数
        9,  // 棱边数（包括边中点形成的子棱边）
        3,  // 面数（3个子三角形）
        {   // 棱边定义
            {0, 1}, {1, 2}, {0, 2},           // 主棱边
            {0, 3}, {3, 1},                   // 边0-1的子棱边
            {1, 4}, {4, 2},                   // 边1-2的子棱边
            {0, 5}, {5, 2}                    // 边0-2的子棱边
        },
        {   // 面定义: 3个子三角形
            {0, 3, 5},
            {3, 1, 4},
            {5, 4, 2}
        }
    };
    
    // QUAD4: 4节点线性四边形
    // 标准来源: VTK_QUAD / ANSYS PLANE42
    defs[ElemType::QUAD4] = {
        4,  // 节点数
        4,  // 棱边数
        1,  // 面数
        {   // 棱边定义: (0,1), (1,2), (2,3), (0,3)
            {0, 1},
            {1, 2},
            {2, 3},
            {0, 3}
        },
        {   // 面定义: [0, 1, 2, 3] (逆时针)
            {0, 1, 2, 3}
        }
    };
    
    // QUAD8: 8节点二次四边形（顶点 + 边中点）
    // 标准来源: VTK_QUADRATIC_QUAD
    defs[ElemType::QUAD8] = {
        8,  // 节点数
        8,  // 棱边数（包括子棱边）
        1,  // 面数
        {   // 棱边定义
            {0, 1}, {1, 2}, {2, 3}, {0, 3},   // 主棱边
            {0, 4}, {4, 1},                    // 子棱边
            {1, 5}, {5, 2},
            {2, 6}, {6, 3},
            {0, 7}, {7, 3}
        },
        {   // 面定义
            {0, 1, 2, 3}
        }
    };

    // ==================== 三维单元 ====================
    
    // TET4: 4节点线性四面体
    // 标准来源: VTK_TETRA / ANSYS SOLID87
    // 节点顺序: 底面(0,1,2) + 顶点(3)
    defs[ElemType::TET4] = {
        4,  // 节点数
        6,  // 棱边数 (C(4,2) = 6)
        4,  // 面数 (4个三角形面)
        {   // 棱边定义: 所有6条可能的棱边
            {0, 1},  // 底面棱边
            {1, 2},  // 底面棱边
            {0, 2},  // 底面棱边
            {0, 3},  // 侧棱边
            {1, 3},  // 侧棱边
            {2, 3}   // 侧棱边
        },
        {   // 面定义: 4个三角形面（外法向朝外）
            {0, 1, 2},  // 底面（从外部看逆时针）
            {0, 3, 1},  // 侧面1
            {1, 3, 2},  // 侧面2
            {2, 3, 0}   // 侧面3
        }
    };
    
    // TET10: 10节点二次四面体（4顶点 + 6边中点）
    // 标准来源: VTK_QUADRATIC_TETRA
    defs[ElemType::TET10] = {
        10,  // 节点数
        18,  // 棱边数
        4,   // 面数
        {   // 棱边定义（主棱边 + 子棱边）
            {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3},  // 主棱边
            {0, 4}, {4, 1},                                    // 边0-1的子棱边
            {1, 5}, {5, 2},                                    // 边1-2的子棱边
            {0, 6}, {6, 2},                                    // 边0-2的子棱边
            {0, 7}, {7, 3},                                    // 边0-3的子棱边
            {1, 8}, {8, 3},                                    // 边1-3的子棱边
            {2, 9}, {9, 3}                                     // 边2-3的子棱边
        },
        {   // 面定义: 4个曲面三角形面
            {0, 1, 2, 4, 5, 6},
            {0, 3, 1, 7, 8, 4},
            {1, 3, 2, 8, 9, 5},
            {2, 3, 0, 9, 7, 6}
        }
    };
    
    // HEX8: 8节点线性六面体
    // 标准来源: VTK_HEXAHEDRON / ANSYS SOLID45/SOLID70
    // 节点顺序: 底面(0,1,2,3) + 顶面(4,5,6,7)，底面和顶面对应点垂直对齐
    defs[ElemType::HEX8] = {
        8,   // 节点数
        12,  // 棱边数
        6,   // 面数
        {    // 棱边定义: 12条棱边
             // 底面4条
            {0, 1}, {1, 2}, {2, 3}, {3, 0},
             // 顶面4条
            {4, 5}, {5, 6}, {6, 7}, {7, 4},
             // 垂直4条
            {0, 4}, {1, 5}, {2, 6}, {3, 7}
        },
        {    // 面定义: 6个四边形面（外法向朝外）
            {0, 1, 2, 3},  // 底面 (-Z方向)
            {4, 5, 6, 7},  // 顶面 (+Z方向)
            {0, 1, 5, 4},  // 前面 (+Y方向)
            {2, 3, 7, 6},  // 后面 (-Y方向)
            {0, 3, 7, 4},  // 左面 (-X方向)
            {1, 2, 6, 5}   // 右面 (+X方向)
        }
    };
    
    // HEX20: 20节点二次六面体（8顶点 + 12边中点）
    // 标准来源: VTK_QUADRATIC_HEXAHEDRON
    defs[ElemType::HEX20] = {
        20,  // 节点数
        32,  // 棵边数
        6,   // 面数
        {    // 棱边定义
             // 底面
            {0, 1}, {1, 2}, {2, 3}, {3, 0},
             // 顶面
            {4, 5}, {5, 6}, {6, 7}, {7, 4},
             // 垂直边
            {0, 4}, {1, 5}, {2, 6}, {3, 7},
             // 子棱边（边中点分割）
            {0, 8}, {8, 1}, {1, 9}, {9, 2},
            {2, 10}, {10, 3}, {3, 11}, {11, 0},
            {4, 12}, {12, 5}, {5, 13}, {13, 6},
            {6, 14}, {14, 7}, {7, 15}, {15, 4},
            {0, 16}, {16, 4}, {1, 17}, {17, 5},
            {2, 18}, {18, 6}, {3, 19}, {19, 7}
        },
        {    // 面定义: 6个曲面四边形面
            {0, 1, 2, 3, 8, 9, 10, 11},       // 底面
            {4, 5, 6, 7, 12, 13, 14, 15},      // 顶面
            {0, 1, 5, 4, 8, 17, 12, 16},       // 前面
            {2, 3, 7, 6, 10, 11, 15, 14},      // 后面
            {0, 3, 7, 4, 11, 19, 15, 16},      // 左面
            {1, 2, 6, 5, 9, 18, 13, 17}        // 右面
        }
    };
    
    // PRISM6: 6节点三棱柱（五面体）
    // 标准来源: VTK_WEDGE / ANSYS SOLID86
    // 节点顺序: 底面三角形(0,1,2) + 顶面三角形(3,4,5)
    defs[ElemType::PRISM6] = {
        6,   // 节点数
        9,   // 棱边数
        5,   // 面数（2个三角形底面 + 3个四边形侧面）
        {    // 棱边定义
             // 底面3条
            {0, 1}, {1, 2}, {0, 2},
             // 顶面3条
            {3, 4}, {4, 5}, {3, 5},
             // 垂直3条
            {0, 3}, {1, 4}, {2, 5}
        },
        {    // 面定义
            {0, 1, 2},        // 底面
            {3, 5, 4},        // 顶面（注意顺序保证外法向）
            {0, 1, 4, 3},     // 侧面1
            {1, 2, 5, 4},     // 侧面2
            {0, 2, 5, 3}      // 侧面3
        }
    };
    
    // PRISM15: 15节点二次三棱柱（6顶点 + 9边中点）
    // 标准来源: VTK_QUADRATIC_WEDGE
    defs[ElemType::PRISM15] = {
        15,  // 节点数
        24,  // 棱边数
        5,   // 面数
        {    // 棱边定义
             // 底面
            {0, 1}, {1, 2}, {0, 2},
             // 顶面
            {3, 4}, {4, 5}, {3, 5},
             // 垂直边
            {0, 3}, {1, 4}, {2, 5},
             // 子棱边
            {0, 6}, {6, 1}, {1, 7}, {7, 2},
            {0, 8}, {8, 2},
            {3, 9}, {9, 4}, {4, 10}, {10, 5},
            {3, 11}, {11, 5},
            {0, 12}, {12, 3}, {1, 13}, {13, 4},
            {2, 14}, {14, 5}
        },
        {    // 面定义
            {0, 1, 2, 6, 7, 8},              // 底面
            {3, 4, 5, 9, 10, 11},            // 顶面
            {0, 1, 4, 3, 6, 13, 9, 12},     // 侧面1
            {1, 2, 5, 4, 7, 14, 10, 13},    // 侧面2
            {0, 2, 5, 3, 8, 14, 11, 12}     // 侧面3
        }
    };
    
    // PYRAMID5: 5节点金字塔（四边形底面 + 尖顶）
    // 标准来源: VTK_PYRAMID
    // 节点顺序: 底面四边形(0,1,2,3) + 尖顶(4)
    defs[ElemType::PYRAMID5] = {
        5,   // 节点数
        8,   // 棱边数
        5,   // 面数（1个四边形底面 + 4个三角形侧面）
        {    // 棱边定义
             // 底面4条
            {0, 1}, {1, 2}, {2, 3}, {3, 0},
             // 侧棱4条（到底面各顶点）
            {0, 4}, {1, 4}, {2, 4}, {3, 4}
        },
        {    // 面定义
            {0, 1, 2, 3},     // 底面
            {0, 1, 4},         // 侧面1
            {1, 2, 4},         // 侧面2
            {2, 3, 4},         // 侧面3
            {3, 0, 4}          // 侧面4
        }
    };
    
    // PYRAMID13: 13节点二次金字塔（5顶点 + 8边中点）
    // 标准来源: VTK_QUADRATIC_PYRAMID
    defs[ElemType::PYRAMID13] = {
        13,  // 节点数
        20,  // 棱边数
        5,   // 面数
        {    // 棱边定义
             // 底面
            {0, 1}, {1, 2}, {2, 3}, {3, 0},
             // 侧棱
            {0, 4}, {1, 4}, {2, 4}, {3, 4},
             // 子棱边
            {0, 5}, {5, 1}, {1, 6}, {6, 2},
            {2, 7}, {7, 3}, {3, 8}, {8, 0},
            {0, 9}, {9, 4}, {1, 10}, {10, 4},
            {2, 11}, {11, 4}, {3, 12}, {12, 4}
        },
        {    // 面定义
            {0, 1, 2, 3, 5, 6, 7, 8},         // 底面
            {0, 1, 4, 5, 10, 9},               // 侧面1
            {1, 2, 4, 6, 11, 10},              // 侧面2
            {2, 3, 4, 7, 12, 11},              // 侧面3
            {3, 0, 4, 8, 9, 12}                // 侧面4
        }
    };
    
    // ==================== Nedelec棱边元几何定义（复用基础类型的拓扑）====================
    // Nedelec一阶棱边元与对应的Lagrange节点元共享完全相同的几何拓扑定义
    // 区别仅在于自由度定义在棱边上而非节点上
    
    // TET4_EDGE: 四面体一阶Nedelec棱边元（6条棱边自由度）
    defs[ElemType::TET4_EDGE] = defs[ElemType::TET4];
    
    // HEX8_EDGE: 六面体一阶Nedelec棱边元（12条棱边自由度）
    defs[ElemType::HEX8_EDGE] = defs[ElemType::HEX8];
    
    // PRISM6_EDGE: 三棱柱一阶Nedelec棱边元（9条棱边自由度）
    defs[ElemType::PRISM6_EDGE] = defs[ElemType::PRISM6];
    
    // PYRAMID5_EDGE: 金字塔一阶Nedelec棱边元（8条棱边自由度）
    defs[ElemType::PYRAMID5_EDGE] = defs[ElemType::PYRAMID5];
    
    return defs;
}

// ==================== 静态存储几何定义 ====================
static const std::unordered_map<ElemType, ElementGeometryDef> geometry_definitions = 
    initialize_geometry_definitions();

// ==================== 公共接口实现 ====================

std::vector<std::tuple<int, int>> ElementGeometry::get_local_edges(ElemType elem_type) {
    auto it = geometry_definitions.find(elem_type);
    if (it == geometry_definitions.end()) {
        throw std::invalid_argument("不支持的单元类型");
    }
    return it->second.edges;
}

std::vector<std::vector<int>> ElementGeometry::get_local_faces(ElemType elem_type) {
    auto it = geometry_definitions.find(elem_type);
    if (it == geometry_definitions.end()) {
        throw std::invalid_argument("不支持的单元类型");
    }
    return it->second.faces;
}

int ElementGeometry::get_num_nodes(ElemType elem_type) {
    auto it = geometry_definitions.find(elem_type);
    if (it == geometry_definitions.end()) {
        throw std::invalid_argument("不支持的单元类型");
    }
    return it->second.num_nodes;
}

int ElementGeometry::get_num_edges(ElemType elem_type) {
    auto it = geometry_definitions.find(elem_type);
    if (it == geometry_definitions.end()) {
        throw std::invalid_argument("不支持的单元类型");
    }
    return it->second.num_edges;
}

int ElementGeometry::get_num_faces(ElemType elem_type) {
    auto it = geometry_definitions.find(elem_type);
    if (it == geometry_definitions.end()) {
        throw std::invalid_argument("不支持的单元类型");
    }
    return it->second.num_faces;
}

bool ElementGeometry::is_supported(ElemType elem_type) {
    return geometry_definitions.find(elem_type) != geometry_definitions.end();
}

std::vector<ElemType> ElementGeometry::get_supported_types() {
    std::vector<ElemType> types;
    types.reserve(geometry_definitions.size());
    for (const auto& pair : geometry_definitions) {
        types.push_back(pair.first);
    }
    return types;
}

} // namespace fe_em
