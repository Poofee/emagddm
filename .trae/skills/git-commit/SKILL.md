---
name: "git-commit"
description: "Automates git commit workflow with Conventional Commits format. Invoke when user asks to commit changes, create a commit, or submit code."
---

# Git Commit Skill

自动化Git提交流程，遵循Conventional Commits规范，确保代码质量和文档完整性。

## 功能特性

### 1. 自动生成提交信息
- 分析代码变更，自动生成符合Conventional Commits规范的提交信息
- 支持的类型：
  - `feat`: 新功能
  - `fix`: Bug修复
  - `docs`: 文档更新
  - `test`: 测试相关
  - `perf`: 性能优化
  - `refactor`: 代码重构
  - `style`: 代码格式
  - `chore`: 构建/工具变更

### 2. 代码质量检查
提交前自动执行：
- 编译检查（确保代码可编译）
- Lint检查（如果项目配置了lint命令）
- Typecheck检查（如果项目配置了typecheck命令）
- 测试运行（如果有相关测试）

### 3. 文档更新提醒
检查并提醒更新以下文档：
- `translation_log.md`: 记录翻译工作
- `todo_list.md`: 任务清单
- `daily_log.md`: 每日工作日志
- `design_decisions.md`: 设计决策

### 4. 执行Git命令
- 自动将变更文件添加到暂存区
- 执行git commit
- 可选推送到远程仓库

## 使用场景

**立即调用此skill当：**
- 用户说"提交代码"、"commit"、"创建提交"
- 用户说"帮我生成commit message"
- 用户说"提交这些更改"
- 用户完成了一个功能或修复，需要提交

## 工作流程

1. **检查工作区状态**
   ```
   git status
   git diff --cached
   git diff
   ```

2. **代码质量检查**
   - 检查项目配置文件（package.json、CMakeLists.txt等）
   - 执行lint和typecheck命令
   - 如有错误，停止并报告

3. **分析变更内容**
   - 查看修改的文件列表
   - 分析diff内容
   - 确定提交类型和范围

4. **生成提交信息**
   格式：
   ```
   <type>(<scope>): <subject>

   <body>

   <footer>
   ```

5. **文档检查**
   - 检查相关文档是否需要更新
   - 提醒用户更新必要文档

6. **执行提交**
   - `git add <files>`
   - `git commit -m "message"`
   - 可选：`git push`

## 提交信息模板

### Fortran到C++翻译提交
```
feat(translation): 翻译<module_name>模块

- 将Fortran子程序<name>翻译为C++17
- 添加单元测试验证数值精度
- 更新translation_log.md

性能: 单核运行时间XXms，为原Fortran的XX%
测试: 通过所有单元测试，相对误差<1e-12
```

### Bug修复提交
```
fix(<component>): 修复<bug_description>

- 问题原因：<root_cause>
- 解决方案：<solution>
- 影响范围：<scope>

Fixes #<issue_number>
```

### 功能添加提交
```
feat(<module>): 添加<feature_description>

- 实现<detail_1>
- 实现<detail_2>
- 添加相应测试

Refs: <reference_document>
```

## 注意事项

1. **禁止提交的内容**
   - 包含敏感信息的文件（密钥、密码等）
   - 未通过测试的代码
   - 编译失败的代码
   - 临时文件和IDE配置文件

2. **必须检查的项目**
   - 代码是否符合项目编码规范
   - 是否需要更新文档
   - 提交信息是否清晰准确

3. **提交粒度**
   - 一个提交解决一个问题
   - 避免大型提交（建议<500行）
   - 相关变更放在同一提交

4. **分支策略**
   - 新功能：feature/<feature-name>
   - Bug修复：fix/<bug-name>
   - 重构：refactor/<refactor-name>

## 示例对话

**用户**: 提交代码

**Assistant**:
1. 检查git状态，发现有3个文件修改
2. 运行lint和typecheck，全部通过
3. 分析变更：翻译了MaxwellSolver模块
4. 生成提交信息：
   ```
   feat(translation): 翻译MaxwellSolver模块

   - 将Fortran子程序maxwell_solver翻译为C++17
   - 添加单元测试验证数值精度
   - 更新translation_log.md

   性能: 单核运行时间45ms，为原Fortran的95%
   测试: 通过所有单元测试，相对误差<1e-12
   ```
5. 提醒：需要更新todo_list.md，标记任务完成
6. 执行git commit

## 错误处理

- **编译失败**: 停止提交，显示错误信息，建议修复方案
- **测试失败**: 停止提交，显示失败的测试用例
- **文档缺失**: 警告用户，但不阻止提交（除非用户明确要求）
- **合并冲突**: 提示用户先解决冲突

## 配置选项

可在项目根目录创建 `.trae/git-commit-config.json` 自定义配置：

```json
{
  "run_tests": true,
  "run_lint": true,
  "run_typecheck": true,
  "check_docs": true,
  "docs_to_check": [
    "translation_log.md",
    "todo_list.md",
    "daily_log.md"
  ],
  "auto_push": false,
  "max_diff_lines": 500
}
```
