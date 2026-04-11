
日志输出使用emagddm\include\tool\logger_factory.hpp的宏。
头文件包含不要有任何路径，直接包含文件名，在cmake文件当中包含路径。

- 禁止使用平台专属语法、编译器专属特性
- 跨平台的语句需要使用单独的平台分支，防止在别的平台报错
- 禁止硬编码（所有可配置参数均需支持JSON配置）
- 禁止全局变量（除FETI-DP全局上下文结构体）
- 禁止代码冗余（重复逻辑封装为工具函数）
- 禁止不写注释
- 禁止内存泄漏

## 详细规范参考

详细的技术规范请查阅 `docs/` 目录下的相应文档：

- [编码风格规范](docs/coding_style.md)
- [内存管理规范](docs/memory_management.md)
- [目录与文件规范](docs/directory_structure.md)
- [命名规范](docs/naming_conventions.md)
- [注释规范](docs/commenting_guidelines.md)
- [第三方库集成规范](docs/third_party_libraries.md)
- [开发流程规范](docs/development_workflow.md)
- [交付物规范](docs/deliverables_requirements.md)
- [调试与日志规范](docs/debugging_logging.md)
- [头文件规范](docs/header_inclusion_guidelines.md)
