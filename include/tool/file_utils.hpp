/**
 * @file file_utils.hpp
 * @brief 基础工具层 - 文件操作工具头文件（简化版）
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include <chrono>
#include <fstream>
#include <filesystem>

namespace tool {

namespace file_utils {

std::string getExtension(const std::string& file_path);
std::string getFilename(const std::string& file_path);
std::string getFilenameWithoutExtension(const std::string& file_path);
std::string getDirectory(const std::string& file_path);
std::string combinePath(const std::string& directory, const std::string& filename);
std::string normalizePath(const std::string& path);
bool isAbsolute(const std::string& path);
bool isRelative(const std::string& path);
std::string changeExtension(const std::string& file_path, const std::string& new_extension);

bool exists(const std::string& path);
bool isFile(const std::string& path);
bool isDirectory(const std::string& path);
uint64_t fileSize(const std::string& file_path);
bool createDirectory(const std::string& dir_path);
bool createDirectories(const std::string& dir_path);
bool remove(const std::string& path);
bool copyFile(const std::string& source, const std::string& target);
bool moveFile(const std::string& source, const std::string& target);

std::vector<uint8_t> readBinary(const std::string& file_path);
bool writeBinary(const std::string& file_path, const std::uint8_t* data, std::size_t size);
bool readText(const std::string& file_path, std::string& content);
bool writeText(const std::string& file_path, const std::string& content);
bool appendText(const std::string& file_path, const std::string& content);

std::vector<std::string> listFiles(const std::string& directory);
std::vector<std::string> listFilesRecursive(const std::string& directory);
std::vector<std::string> findFilesByExtension(const std::string& directory, const std::string& extension);

class FileWatcher {
public:
    explicit FileWatcher(const std::string& file_path);
    bool hasChanged();
    std::chrono::system_clock::time_point lastWriteTime() const;

private:
    std::string file_path_;
    std::chrono::system_clock::time_point last_write_time_;
};

struct FileInfo {
    std::string path;
    std::string filename;
    std::string extension;
    uint64_t size;
    bool is_directory;
    std::chrono::system_clock::time_point creation_time;
    std::chrono::system_clock::time_point modification_time;
    std::chrono::system_clock::time_point last_access_time;
};

std::optional<FileInfo> getInfo(const std::string& path);

bool isPathValid(const std::string& path);
std::string getTempDirectory();
std::string getHomeDirectory();
std::string getCurrentDirectory();
std::string getExecutablePath();
bool setCurrentDirectory(const std::string& directory);

} // namespace file_utils

} // namespace tool
