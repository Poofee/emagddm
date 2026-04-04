/**
 * @file string_format.hpp
 * @brief 独立字符串格式化工具
 * @details 提供与日志库解耦的通用字符串格式化功能。
 *          格式化后端可通过此文件统一替换（如迁移至C++20 std::format），
 *          不影响上层调用方。当前基于项目已有的fmt库实现。
 * @author Poofee
 * @date 2026-04-04
 */

#pragma once

#include <string>
#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <spdlog/fmt/bundled/format.h>

namespace tool {

/**
 * @brief 格式化字符串（支持{}占位符）
 * @tparam Args 可变参数类型
 * @param fmt_str 格式字符串，使用{}作为占位符
 * @param args 要填充的参数列表
 * @return std::string 格式化后的结果字符串
 *
 * @code
 * auto s = formatString("值: {}, 名称: {}", 42, "test");
 * // s == "值: 42, 名称: test"
 * @endcode
 *
 * @note 格式化后端封装在此处，替换后端时仅需修改本文件实现，
 *       上层Logger等调用方无需任何改动
 */
template<typename... Args>
std::string formatString(const std::string& fmt_str, Args&&... args) {
    return fmt::format(fmt_str, std::forward<Args>(args)...);
}

} // namespace tool
