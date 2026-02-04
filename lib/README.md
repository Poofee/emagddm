# 第三方库目录说明

本目录用于存放项目所需的第三方库文件。

## 库目录结构

- `eigen/` - Eigen线性代数库（单头文件库）
- `googletest/` - Google Test单元测试框架
- `pugixml/` - pugixml XML解析库（单头文件库）
- `spdlog/` - spdlog日志库（单头文件库）

## 库集成说明

### 单头文件库（Eigen、pugixml、spdlog）
这些库是单头文件库，使用时只需：
1. 将库的头文件放入对应目录
2. 在CMakeLists.txt中添加包含路径
3. 直接通过`#include`引用

### 需要编译的库（googletest）
googletest需要编译，使用时：
1. 下载源码放入对应目录
2. 在CMakeLists.txt中通过`add_subdirectory`包含
3. 链接生成的库文件

## 库版本要求

- **Eigen**: 3.4+ 
- **googletest**: 1.12+
- **pugixml**: 1.12+
- **spdlog**: 1.11+

## 下载和安装

请从以下官方地址下载对应库的最新版本：

- Eigen: https://gitlab.com/libeigen/eigen
- googletest: https://github.com/google/googletest
- pugixml: https://github.com/zeux/pugixml
- spdlog: https://github.com/gabime/spdlog

下载后将库文件放入对应目录即可。