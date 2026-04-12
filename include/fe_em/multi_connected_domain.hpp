/**
 * @file multi_connected_domain.hpp
 * @brief 电磁物理层 - 多连通域检测器头文件
 * @details 自动检测有限元网格中的多连通域（孔洞）并识别对应的环量自由度。
 *          用于Tree Gauge规范中保留不可置零的磁通回路自由度。
 *
 * @note 数学背景：
 *       对于有g个孔洞的多连通域：
 *       - 生成树的边数 = V - 1 + g（而非简单的V-1）
 *       - 每个孔洞对应一个环量自由度（不可置为0）
 *       - 这些环量自由度对应于围绕孔洞的独立回路
 *
 * @note 核心功能：
 *       - 基于网格拓扑分析识别内部边界（孔洞边界）
 *       - 使用并查集（Union-Find）检测连通分量
 *       - 为每个孔洞选择独立的余树边作为环量自由度
 *
 * @note 与后续模块对接：
 *       - TreeGauge：使用检测结果保留环量自由度
 *       - SpanningTreeBuilder：基于生成树的余树边选择环量边
 *
 * @author Poofee
 * @date 2026-04-12
 * @version 1.0
 */

#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include "em_mesh_data.hpp"
#include "spanning_tree.hpp"
#include "logger_factory.hpp"

namespace fe_em {

// ==================== 输出数据结构 ====================

/**
 * @struct DetectionResult
 * @brief 多连通域检测结果数据结构
 * @details 包含孔洞数量、每个孔洞的回路棱边、以及需要保留的环量自由度集合。
 *
 * @note 数据结构设计原则：
 *       - num_holes: 亏格genus，即独立的孔洞数量
 *       - hole_loop_edges: 二维向量，每行对应一个孔洞的闭合回路棱边ID列表
 *       - circulation_dofs: 余树边的子集，这些边代表穿过孔洞的独立磁通回路
 */
struct DetectionResult {
    int num_holes = 0;                                          ///< 孔洞数量（亏格genus）
    std::vector<std::vector<int>> hole_loop_edges;              ///< 每个孔洞对应的回路棱边集合（全局棱边ID）
    std::set<int> circulation_dofs;                             ///< 需要保留的环量自由度（余树边子集，全局棱边ID）

    /**
     * @brief 检查检测结果是否有效
     * @return bool 如果结果包含有效数据返回true
     */
    bool isValid() const {
        return num_holes >= 0 &&
               static_cast<int>(hole_loop_edges.size()) == num_holes &&
               static_cast<int>(circulation_dofs.size()) >= num_holes;
    }
};

// ==================== 核心工具类 ====================

/**
 * @class MultiConnectedDomainDetector
 * @brief 多连通域检测器
 * @details 基于网格拓扑结构和生成树信息，自动识别多连通域中的孔洞，
 *          并为每个孔洞选择合适的余树边作为环量自由度。
 *
 * @note 算法流程（Phase 1简化版）：
 *       1. **输入验证**：检查网格数据和生成树结果的完整性
 *       2. **边界分析**：从boundary_markers中提取边界信息
 *       3. **拓扑分析**：基于图论方法识别可能的孔洞候选
 *       4. **环量选择**：从余树边中选择独立的边作为环量自由度
 *       5. **结果验证**：确保circulation_dofs是cotree_edges的子集
 *
 * @note Phase 2优化方向（TODO）：
 *       - 基于对偶图的精确孔洞检测算法
 *       - 代数拓扑Betti数计算
 *       - 几何验证（使用节点坐标判断内外边界）
 *       - 支持复杂多孔洞拓扑（嵌套孔洞、非简单连通等）
 *
 * @note 设计原则：
 *       - 静态方法设计，无状态，线程安全
 *       - 输入仅依赖常量引用，不修改原始数据
 *       - 返回值语义清晰，便于单元测试
 *       - 算法复杂度可控：O(V + E)（V=节点数，E=棱边数）
 *
 * @note 典型使用流程：
 *       @code
 *       // 第一步：构建生成树（已在TreeGauge中完成）
 *       SpanningTreeBuilder::Result tree_result = SpanningTreeBuilder::buildBFSTree(...);
 *
 *       // 第二步：检测多连通域
 *       DetectionResult domain_result = MultiConnectedDomainDetector::detect(
 *           mesh_data, tree_result
 *       );
 *
 *       // 第三步：使用结果
 *       if (domain_result.num_holes > 0) {
 *           FEEM_INFO("检测到{}个孔洞，需保留{}个环量自由度",
 *                     domain_result.num_holes, domain_result.circulation_dofs.size());
 *           // 在TreeGauge中保留这些自由度不被置零
 *       }
 *       @endcode
 *
 * @note 性能特征：
 *       - 时间复杂度：O(V + E + B)，B为边界标记数
 *       - 空间复杂度：O(V + E)，主要用于并查集和临时容器
 *       - 内存优化：预分配容器容量，避免动态扩容
 */
class MultiConnectedDomainDetector {
public:
    /**
     * @brief 检测多连通域并识别环量自由度（静态工厂方法）
     * @param mesh_data 有限元网格数据的常量引用（用于边界信息）
     * @param tree_result SpanningTreeBuilder的生成树结果（用于余树边信息）
     * @return DetectionResult 包含孔洞数量、回路棱边、环量自由度的完整结果
     * @throws std::invalid_argument 当输入数据无效时抛出
     *
     * @details 核心算法步骤（Phase 1简化实现）：
     *          1. **输入验证**
     *             - 检查mesh_data是否包含有效的节点和单元数据
     *             - 检查tree_result是否包含有效的树边和余树边集合
     *
     *          2. **边界信息提取**（简化版）
     *             - 遍历boundary_markers收集所有边界类型
     *             - 统计不同类型的边界数量
     *             - 识别可能存在内部边界的区域
     *
     *          3. **拓扑分析**（基于图连通性）
     *             - 使用并查集分析图的连通性
     *             - 计算基本圈基（fundamental cycle basis）
     *             - 识别与内部边界相关的独立回路
     *
     *          4. **环量自由度选择**
     *             - 从cotree_edges中选择与孔洞相关的边
     *             - 确保选择的边线性无关（构成基本圈基）
     *             - 将选中的边加入circulation_dofs集合
     *
     *          5. **结果后处理**
     *             - 验证circulation_dofs ⊆ cotree_edges
     *             - 构建hole_loop_edges（每个孔洞的闭合回路）
     *             - 输出详细统计日志
     *
     * @warning 调用前确保：
     *          - mesh_data已正确填充节点、单元、边界标记
     *          - tree_result已由SpanningTreeBuilder正确构建
     *
     * @note 此方法可重复调用，每次调用会重新计算并返回新结果
     *
     * @note 日志输出示例：
     *       [INFO] ========== MultiConnectedDomainDetector: 多连通域检测完成 ==========
     *       [INFO] 输入统计: 节点数=10000, 单元数=50000, 边界标记数=10
     *       [INFO] 生成树信息: 树边数=9999, 余树边数=45001
     *       [INFO] 检测结果: 孔洞数=1, 环量自由度数=1
     *       [INFO] 环量边: [12345]
     *       [INFO] =======================================================
     */
    static DetectionResult detect(
        const EMMeshData& mesh_data,
        const SpanningTreeResult& tree_result
    );

private:
    /**
     * @brief 验证输入数据的有效性
     * @param mesh_data 网格引用
     * @param tree_result 生成树结果引用
     * @throws std::invalid_argument 当输入数据无效时抛出
     * @details 检查项包括：
     *          - 网格数据不为空（至少有节点和单元）
     *          - 生成树结果包含有效的树边和余树边集合
     *          - 余树边集合非空（否则无法选择环量自由度）
     */
    static void validateInput(
        const EMMeshData& mesh_data,
        const SpanningTreeResult& tree_result
    );

    /**
     * @brief 分析边界标记信息，提取潜在的孔洞候选
     * @param mesh_data 网格数据引用
     * @param out_boundary_types 输出：发现的边界类型集合
     * @param out_inner_boundary_count 输出：潜在内部边界数量
     * @details 基于boundary_markers进行启发式分析：
     *          - 识别PERFECT_E/PERFECT_H类型的边界（可能为内部导体边界）
     *          - 统计不同名称的边界（如"Hole_1", "Hole_2"等命名模式）
     *          - 基于目标实体数量估计可能的孔洞数
     *
     * @note Phase 1实现为简化版，主要依赖边界标记的元数据分析
     * @note Phase 2将结合几何计算进行精确判定
     */
    static void analyzeBoundaryMarkers(
        const EMMeshData& mesh_data,
        std::set<int>& out_boundary_types,
        int& out_inner_boundary_count
    );

    /**
     * @brief 基于图拓扑选择环量自由度（核心算法）
     * @param cotree_edges 余树边集合（全局棱边ID）
     * @param tree_edges 树边集合（全局棱边ID）
     * @param node_parent 节点→父节点映射（来自生成树）
     * @param estimated_hole_count 估计的孔洞数量（来自边界分析）
     * @return std::set<int> 选中的环量自由度集合（cotree_edges的子集）
     * @details 选择策略（Phase 1简化版）：
     *          1. 若estimated_hole_count == 0：返回空集（单连通域）
     *          2. 若estimated_hole_count > 0但无法精确识别：
     *             - 使用贪心策略从cotree_edges中选择前N条边
     *             - N = min(estimated_hole_count, cotree_edges.size())
     *             - 选择标准：优先选择连接不同子树的边（增加多样性）
     *
     * @note 选择原则：
     *       - 选中的边必须属于cotree_edges（保证是余树边）
     *       - 不同孔洞的环量边应尽量线性无关（避免冗余）
     *       - 优先选择"跨越性强"的边（连接距离较远的节点）
     *
     * @note Phase 2优化方向：
     *       - 基于基本圈基（fundamental cycles）的精确选择
     *       - 使用对偶图找到真正包围孔洞的回路
     *       - 结合几何信息验证选择的正确性
     */
    static std::set<int> selectCirculationDOFs(
        const std::set<int>& cotree_edges,
        const std::set<int>& tree_edges,
        const std::unordered_map<int, int>& node_parent,
        int estimated_hole_count
    );

    /**
     * @brief 构建每个孔洞的回路棱边列表
     * @param circulation_dofs 环量自由度集合
     * @param tree_edges 树边集合
     * @param node_parent 节点→父节点映射
     * @param adjacency_list 节点邻接表（可选，用于路径查找）
     * @return std::vector<std::vector<int>> 每个孔洞的回路棱边列表
     * @details 对每个环量自由度（余树边），构建其对应的基本圈（fundamental cycle）：
     *          1. 取余树边(u, v)
     *          2. 在生成树中查找u到v的唯一路径
     *          3. 该路径 + 边(u,v) 构成一个基本圈
     *          4. 这个基本圈即为该孔洞的回路
     *
     * @note Phase 1实现返回简化的回路（仅包含环量边本身）
     * @note Phase 2将返回完整的闭合回路棱边序列
     */
    static std::vector<std::vector<int>> buildHoleLoopEdges(
        const std::set<int>& circulation_dofs,
        const std::set<int>& tree_edges,
        const std::unordered_map<int, int>& node_parent,
        const std::vector<std::vector<int>>& adjacency_list = {}
    );

    /**
     * @brief 验证结果的一致性
     * @param result 待验证的结果
     * @param cotree_edges 原始余树边集合（用于子集验证）
     * @return bool 如果结果一致返回true
     * @details 验证规则：
     *          - num_holes >= 0
     *          - hole_loop_edges.size() == num_holes
     *          - circulation_dofs中的每个元素都属于cotree_edges
     *          - circulation_dofs.size() >= num_holes（每个孔洞至少一个环量自由度）
     */
    static bool validateResult(
        const DetectionResult& result,
        const std::set<int>& cotree_edges
    );
};

} // namespace fe_em
