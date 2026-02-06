/**
 * @file file_utils.cpp
 * @brief 基础工具层 - 文件操作工具源文件
 */

#include "tool/file_utils.hpp"
#include <algorithm>
#include <array>
#include <cstdio>
#include <sstream>

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

namespace tool::file_utils {

std::string getExtension(const std::string& file_path) {
    auto pos = file_path.find_last_of('.');
    if (pos != std::string::npos && pos < file_path.length() - 1) {
        return file_path.substr(pos);
    }
    return "";
}

std::string getFilename(const std::string& file_path) {
    auto pos = file_path.find_last_of("/\\");
    if (pos != std::string::npos) {
        return file_path.substr(pos + 1);
    }
    return file_path;
}

std::string getFilenameWithoutExtension(const std::string& file_path) {
    std::string filename = getFilename(file_path);
    auto pos = filename.find_last_of('.');
    if (pos != std::string::npos) {
        return filename.substr(0, pos);
    }
    return filename;
}

std::string getDirectory(const std::string& file_path) {
    auto pos = file_path.find_last_of("/\\");
    if (pos != std::string::npos) {
        return file_path.substr(0, pos);
    }
    return "";
}

std::string combinePath(const std::string& directory, const std::string& filename) {
    if (directory.empty()) {
        return filename;
    }
    char last_char = directory.back();
    if (last_char == '/' || last_char == '\\') {
        return directory + filename;
    }
    return directory + "\\" + filename;
}

std::string normalizePath(const std::string& path) {
    std::string result = path;
    std::replace(result.begin(), result.end(), '/', '\\');
    return result;
}

bool isAbsolute(const std::string& path) {
#if defined(_WIN32)
    return path.size() >= 2 && path[1] == ':';
#else
    return !path.empty() && path[0] == '/';
#endif
}

bool isRelative(const std::string& path) {
    return !isAbsolute(path);
}

std::string changeExtension(const std::string& file_path, const std::string& new_extension) {
    std::string ext = new_extension;
    if (!ext.empty() && ext[0] != '.') {
        ext = "." + ext;
    }
    auto pos = file_path.find_last_of('.');
    if (pos != std::string::npos) {
        return file_path.substr(0, pos) + ext;
    }
    return file_path + ext;
}

bool exists(const std::string& path) {
#if defined(_WIN32)
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES;
#else
    struct stat st;
    return stat(path.c_str(), &st) == 0;
#endif
}

bool isFile(const std::string& path) {
#if defined(_WIN32)
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && !(attr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return S_ISREG(st.st_mode);
#endif
}

bool isDirectory(const std::string& path) {
#if defined(_WIN32)
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return S_ISDIR(st.st_mode);
#endif
}

uint64_t fileSize(const std::string& file_path) {
#if defined(_WIN32)
    HANDLE hFile = CreateFileA(file_path.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) return 0;
    LARGE_INTEGER size;
    BOOL success = GetFileSizeEx(hFile, &size);
    CloseHandle(hFile);
    return success ? static_cast<uint64_t>(size.QuadPart) : 0;
#else
    struct stat st;
    if (stat(file_path.c_str(), &st) != 0) return 0;
    return static_cast<uint64_t>(st.st_size);
#endif
}

bool createDirectory(const std::string& dir_path) {
#if defined(_WIN32)
    return CreateDirectoryA(dir_path.c_str(), NULL) != 0 || GetLastError() == ERROR_ALREADY_EXISTS;
#else
    return mkdir(dir_path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

bool createDirectories(const std::string& dir_path) {
    if (dir_path.empty()) return false;
    if (exists(dir_path)) return isDirectory(dir_path);

    std::string parent = getDirectory(dir_path);
    if (!parent.empty() && !exists(parent)) {
        if (!createDirectories(parent)) return false;
    }

    return createDirectory(dir_path);
}

bool remove(const std::string& path) {
#if defined(_WIN32)
    return isDirectory(path) ? RemoveDirectoryA(path.c_str()) : DeleteFileA(path.c_str());
#else
    return ::remove(path.c_str()) == 0;
#endif
}

bool copyFile(const std::string& source, const std::string& target) {
#if defined(_WIN32)
    return CopyFileA(source.c_str(), target.c_str(), FALSE) != 0;
#else
    std::ifstream src(source, std::ios::binary);
    if (!src) return false;
    std::ofstream dst(target, std::ios::binary);
    if (!dst) return false;
    dst << src.rdbuf();
    return true;
#endif
}

bool moveFile(const std::string& source, const std::string& target) {
#if defined(_WIN32)
    return MoveFileA(source.c_str(), target.c_str()) != 0;
#else
    return rename(source.c_str(), target.c_str()) == 0;
#endif
}

std::vector<uint8_t> readBinary(const std::string& file_path) {
    std::vector<uint8_t> data;
    std::ifstream file(file_path, std::ios::binary);
    if (!file) return data;

    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    if (size <= 0) return data;

    file.seekg(0, std::ios::beg);
    data.resize(static_cast<size_t>(size));
    if (!file.read(reinterpret_cast<char*>(data.data()), size)) {
        data.clear();
    }
    return data;
}

bool writeBinary(const std::string& file_path, const std::uint8_t* data, std::size_t size) {
    std::ofstream file(file_path, std::ios::binary);
    if (!file) return false;
    file.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(size));
    return file.good();
}

bool readText(const std::string& file_path, std::string& content) {
    std::ifstream file(file_path);
    if (!file) return false;

    std::stringstream buffer;
    buffer << file.rdbuf();
    content = buffer.str();
    return true;
}

bool writeText(const std::string& file_path, const std::string& content) {
    std::ofstream file(file_path);
    if (!file) return false;
    file << content;
    return file.good();
}

bool appendText(const std::string& file_path, const std::string& content) {
    std::ofstream file(file_path, std::ios::app);
    if (!file) return false;
    file << content;
    return file.good();
}

std::vector<std::string> listFiles(const std::string& directory) {
    std::vector<std::string> files;
#if defined(_WIN32)
    std::string pattern = directory + "\\*";
    WIN32_FIND_DATAA find_data;
    HANDLE hFind = FindFirstFileA(pattern.c_str(), &find_data);
    if (hFind == INVALID_HANDLE_VALUE) return files;
    do {
        std::string name = find_data.cFileName;
        if (name != "." && name != "..") {
            files.push_back(combinePath(directory, name));
        }
    } while (FindNextFileA(hFind, &find_data));
    FindClose(hFind);
#else
    DIR* dir = opendir(directory.c_str());
    if (!dir) return files;
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name = entry->d_name;
        if (name != "." && name != "..") {
            files.push_back(combinePath(directory, name));
        }
    }
    closedir(dir);
#endif
    return files;
}

std::vector<std::string> listFilesRecursive(const std::string& directory) {
    std::vector<std::string> all_files;
    std::vector<std::string> dirs = {directory};

    while (!dirs.empty()) {
        std::string current = dirs.back();
        dirs.pop_back();

        auto files = listFiles(current);
        for (const auto& file : files) {
            if (isDirectory(file)) {
                dirs.push_back(file);
            } else {
                all_files.push_back(file);
            }
        }
    }
    return all_files;
}

std::vector<std::string> findFilesByExtension(const std::string& directory, const std::string& extension) {
    std::vector<std::string> result;
    std::string ext = extension;
    if (!ext.empty() && ext[0] != '.') ext = "." + ext;

    auto all_files = listFilesRecursive(directory);
    for (const auto& file : all_files) {
        if (getExtension(file) == ext) {
            result.push_back(file);
        }
    }
    return result;
}

FileWatcher::FileWatcher(const std::string& file_path) : file_path_(file_path) {
    if (exists(file_path)) {
        last_write_time_ = lastWriteTime();
    }
}

bool FileWatcher::hasChanged() {
    if (!exists(file_path_)) return true;
    return lastWriteTime() != last_write_time_;
}

std::chrono::file_time_type FileWatcher::lastWriteTime() const {
#if defined(_WIN32)
    HANDLE hFile = CreateFileA(file_path_.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) return {};
    FILETIME ft;
    BOOL success = GetFileTime(hFile, NULL, NULL, &ft);
    CloseHandle(hFile);
    if (!success) return {};
    return std::chrono::file_time_type(std::chrono::seconds(ft.dwLowDateTime));
#else
    struct stat st;
    if (stat(file_path_.c_str(), &st) != 0) return {};
    return std::chrono::file_time_type(std::chrono::seconds(st.st_mtime));
#endif
}

std::optional<FileInfo> getInfo(const std::string& path) {
    if (!exists(path)) return std::nullopt;

    FileInfo info;
    info.path = path;
    info.filename = getFilename(path);
    info.extension = getExtension(path);
    info.size = isFile(path) ? fileSize(path) : 0;
    info.is_directory = isDirectory(path);

    return info;
}

bool isPathValid(const std::string& path) {
    if (path.empty()) return false;
#if defined(_WIN32)
    std::string normalized = normalizePath(path);
    if (normalized.size() >= 2 && normalized[1] == ':') {
        return true;
    }
    return normalized.size() >= 3 && normalized[0] == '\\' && normalized[1] == '\\';
#else
    return !path.empty();
#endif
}

std::string getTempDirectory() {
#if defined(_WIN32)
    char temp_path[MAX_PATH];
    GetTempPathA(MAX_PATH, temp_path);
    return std::string(temp_path);
#else
    const char* temp = getenv("TMPDIR");
    return temp ? std::string(temp) : "/tmp";
#endif
}

std::string getHomeDirectory() {
#if defined(_WIN32)
    char home_path[MAX_PATH];
    GetEnvironmentVariableA("USERPROFILE", home_path, MAX_PATH);
    return std::string(home_path);
#else
    const char* home = getenv("HOME");
    return home ? std::string(home) : std::string();
#endif
}

std::string getCurrentDirectory() {
    char buffer[4096];
#if defined(_WIN32)
    GetCurrentDirectoryA(4096, buffer);
#else
    getcwd(buffer, 4096);
#endif
    return std::string(buffer);
}

std::string getExecutablePath() {
#if defined(_WIN32)
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    return std::string(buffer);
#else
    char buffer[4096];
    readlink("/proc/self/exe", buffer, 4096);
    return std::string(buffer);
#endif
}

bool setCurrentDirectory(const std::string& directory) {
#if defined(_WIN32)
    return SetCurrentDirectoryA(directory.c_str()) != 0;
#else
    return chdir(directory.c_str()) == 0;
#endif
}

} // namespace file_utils

} // namespace tool
