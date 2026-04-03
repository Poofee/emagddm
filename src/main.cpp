/**
 * @file main.cpp
 * @brief 驱动层主函数
 * @details FETI-DP求解器的主入口程序，负责命令行参数解析、应用初始化和运行
 *          支持两种运行模式：
 *          1. JSON配置文件模式：使用配置文件运行电磁场求解
 *          2. AEDT解析模式：从ANSYS Maxwell项目文件提取数据
 * @author Poofee
 * @date 2026-XX-XX
 * @version 1.0
 */

#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include "solver_app.hpp"
#include "maxwell_parser_impl.hpp"
#include "base64_utils.hpp"
#include "logger_factory.hpp"
#include "json.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

/**
 * @brief 显示程序使用帮助信息
 */
void showUsage() {
    std::cout << "FETI-DP电磁有限元求解器" << std::endl;
    std::cout << "用法:" << std::endl;
    std::cout << "  fetidp_solver <配置文件路径>          使用JSON配置文件运行求解器" << std::endl;
    std::cout << "  fetidp_solver --aedt <aedt文件路径>   解析AEDT项目文件并提取数据" << std::endl;
    std::cout << std::endl;
    std::cout << "示例:" << std::endl;
    std::cout << "  fetidp_solver config/motor2d_steady.json" << std::endl;
    std::cout << "  fetidp_solver --aedt docs/project/Project49.aedt" << std::endl;
}

/**
 * @brief 获取文件的扩展名（转换为小写）
 * @param file_path 文件路径
 * @return 小写扩展名
 */
std::string getFileExtension(const std::string& file_path) {
    size_t pos = file_path.find_last_of('.');
    if (pos == std::string::npos) {
        return "";
    }
    std::string ext = file_path.substr(pos);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

/**
 * @brief 从AEDT文件中提取Base64编码的图像数据
 * @param file_path AEDT文件路径
 * @return Base64编码的图像数据
 */
std::string extractImage64FromFile(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        return "";
    }
    
    std::string line;
    std::string base64_data;
    bool in_image_section = false;
    
    while (std::getline(file, line)) {
        if (line.find("Image64='") != std::string::npos) {
            in_image_section = true;
            size_t start_pos = line.find("Image64='") + 9;
            base64_data = line.substr(start_pos);
            if (!base64_data.empty() && base64_data.back() == '\'') {
                base64_data.pop_back();
                break;
            }
            continue;
        }
        
        if (in_image_section) {
            for (char c : line) {
                if (c == '\'') {
                    in_image_section = false;
                    break;
                }
                if (c != '\\' && c != '\n' && c != '\r') {
                    base64_data += c;
                }
            }
            if (!in_image_section) {
                break;
            }
        }
    }
    
    file.close();
    return base64_data;
}

/**
 * @brief 解析AEDT项目文件并保存提取的数据
 * @param aedt_path AEDT项目文件路径
 * @return 处理成功返回true
 */
bool processAedtFile(const std::string& aedt_path) {
    FEEM_INFO("========================================");
    FEEM_INFO("解析AEDT项目文件: {}", aedt_path);
    FEEM_INFO("========================================");
    
    // 验证输入文件是否存在
    if (!fs::exists(aedt_path)) {
        FEEM_ERROR("文件不存在: {}", aedt_path);
        return false;
    }
    
    // 创建Maxwell解析器实例
    tool::MaxwellParserImpl parser;
    
    // 验证文件格式是否受支持
    if (!parser.canParse(aedt_path)) {
        FEEM_ERROR("不支持的文件格式: {}", aedt_path);
        return false;
    }
    
    // 步骤1：解析文件基本信息
    FEEM_INFO("");
    FEEM_INFO("[1/6] 解析文件信息...");
    auto file_info = parser.parseFileInfo();
    FEEM_INFO("  项目名称: {}", file_info.project_name);
    FEEM_INFO("  Maxwell版本: {}", file_info.maxwell_version);
    FEEM_INFO("  仿真类型: {}", file_info.simulation_type);
    FEEM_INFO("  维度: {}", file_info.dimension);
    FEEM_INFO("  文件大小: {} 字节", file_info.file_size);
    
    // 步骤2：解析材料数据
    FEEM_INFO("");
    FEEM_INFO("[2/6] 解析材料数据...");
    auto materials = parser.parseMaterials();
    FEEM_INFO("  找到 {} 个材料:", materials.size());
    for (const auto& mat : materials) {
        FEEM_INFO("    - {}", mat.value("name", "Unknown"));
    }
    
    // 步骤3：解析边界条件
    FEEM_INFO("");
    FEEM_INFO("[3/6] 解析边界条件...");
    auto boundaries = parser.parseBoundaries();
    FEEM_INFO("  找到 {} 个边界条件:", boundaries.size());
    for (const auto& bnd : boundaries) {
        FEEM_INFO("    - {}", bnd.value("name", "Unknown"));
    }
    
    // 步骤4：解析激励源
    FEEM_INFO("");
    FEEM_INFO("[4/6] 解析激励源...");
    auto excitations = parser.parseExcitations();
    FEEM_INFO("  找到 {} 个激励源:", excitations.size());
    for (const auto& exc : excitations) {
        FEEM_INFO("    - {}", exc.value("name", "Unknown"));
    }
    
    // 步骤5：解析求解设置
    FEEM_INFO("");
    FEEM_INFO("[5/6] 解析求解设置...");
    auto solution_setup = parser.parseSolutionSetup();
    FEEM_INFO("  求解器名称: {}", solution_setup.value("name", "Unknown"));
    FEEM_INFO("  求解类型: {}", solution_setup.value("solution_type", "Unknown"));
    
    // 步骤6：提取并保存预览图像
    FEEM_INFO("");
    FEEM_INFO("[6/6] 提取预览图像...");
    std::string base64_image = extractImage64FromFile(aedt_path);
    std::string image_saved_path;
    
    if (!base64_image.empty()) {
        // 将Base64字符串解码为二进制图像数据
        std::vector<uint8_t> image_data = fe_em::tool::Base64Utils::decode(base64_image);
        
        // 根据二进制头自动检测图像格式
        std::string image_type = fe_em::tool::Base64Utils::detectImageType(image_data);
        
        // 创建输出目录
        fs::path output_dir = fs::path(aedt_path).parent_path() / "output";
        fs::create_directories(output_dir);
        
        // 构造输出文件名
        std::string base_name = fs::path(aedt_path).stem().string();
        image_saved_path = (output_dir / (base_name + "_preview" + image_type)).string();
        
        // 以二进制模式写入图像文件
        std::ofstream img_file(image_saved_path, std::ios::binary);
        if (img_file.is_open()) {
            img_file.write(reinterpret_cast<const char*>(image_data.data()), image_data.size());
            img_file.close();
            FEEM_INFO("  图像已保存: {}", image_saved_path);
            FEEM_INFO("  图像大小: {} 字节", image_data.size());
        } else {
            FEEM_ERROR("无法创建图像文件: {}", image_saved_path);
        }
    } else {
        FEEM_WARN("未找到预览图像");
    }
    
    // 构建输出JSON对象
    json output_json;
    output_json["file_info"] = file_info.toJson();
    output_json["materials"] = materials;
    output_json["boundaries"] = boundaries;
    output_json["excitations"] = excitations;
    output_json["solution_setup"] = solution_setup;
    
    if (!image_saved_path.empty()) {
        output_json["extracted_image"] = image_saved_path;
    }
    
    // 确保输出目录存在
    fs::path output_dir = fs::path(aedt_path).parent_path() / "output";
    fs::create_directories(output_dir);
    
    // 构造JSON输出文件名
    std::string base_name = fs::path(aedt_path).stem().string();
    std::string output_json_path = (output_dir / (base_name + "_data.json")).string();
    
    // 写入JSON文件
    std::ofstream out_file(output_json_path);
    if (out_file.is_open()) {
        out_file << output_json.dump(4);
        out_file.close();
        FEEM_INFO("");
        FEEM_INFO("========================================");
        FEEM_INFO("数据已保存到: {}", output_json_path);
        FEEM_INFO("========================================");
    } else {
        FEEM_ERROR("无法创建输出文件: {}", output_json_path);
    }
    
    return true;
}

/**
 * @brief 程序主入口函数
 * @param argc 命令行参数个数
 * @param argv 命令行参数数组
 * @return 程序退出码
 */
int main(int argc, char* argv[]) {
    // 初始化日志系统
    tool::LoggerFactory::initializeDefaultLogger("", true);
    
    int exit_code = 0;
    
    try {
        // 至少需要一个命令行参数
        if (argc < 2) {
            showUsage();
            return 1;
        }
        
        std::string arg1 = argv[1];
        
        // 处理--aedt选项
        if (arg1 == "--aedt" || arg1 == "-a") {
            if (argc < 3) {
                FEEM_ERROR("请指定AEDT文件路径");
                showUsage();
                return 1;
            }
            
            std::string aedt_path = argv[2];
            exit_code = processAedtFile(aedt_path) ? 0 : 1;
        } 
        // 处理--help选项
        else if (arg1 == "--help" || arg1 == "-h") {
            showUsage();
            exit_code = 0;
        } 
        // 默认处理
        else {
            std::string config_file = arg1;
            std::string ext = getFileExtension(config_file);
            
            // 如果是.aedt文件，自动切换到AEDT解析模式
            if (ext == ".aedt") {
                exit_code = processAedtFile(config_file) ? 0 : 1;
            } 
            // 否则视为JSON配置文件，启动求解器
            else {
                app::SolverApp solver;
                
                if (!solver.initialize(config_file)) {
                    FEEM_ERROR("求解器初始化失败: {}", config_file);
                    exit_code = 1;
                } else if (!solver.run()) {
                    FEEM_ERROR("求解器运行失败");
                    exit_code = 1;
                } else {
                    FEEM_INFO("求解器运行成功完成");
                    exit_code = 0;
                }
            }
        }
    } catch (const std::exception& e) {
        FEEM_ERROR("程序异常: {}", e.what());
        exit_code = 1;
    } catch (...) {
        FEEM_ERROR("未知异常");
        exit_code = 1;
    }
    
    // 刷新日志缓冲区
    tool::LoggerFactory::getDefaultLogger().flush();
    
    return exit_code;
}
