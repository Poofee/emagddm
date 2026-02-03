# 编码风格规范

## 缩进与换行
- 使用 4 个空格缩进（禁止Tab）
- 每个语句单独一行
- 函数体、循环体、条件语句的大括号 `{}` 单独成行
- 左括号紧跟语句末尾，右括号与左括号对齐

```cpp
// 正确示例
if (isConverged) {
    mergeGlobalField();
    exportResult();
}

// 错误示例
if (isConverged){mergeGlobalField();exportResult();}
```

## 空行与空格
- 函数之间、类的成员函数之间保留1个空行
- 运算符（+、-、*、/、=、==等）前后各留1个空格
- 逗号、分号后留1个空格，前不留空格

```cpp
// 正确示例
for (int i = 0; i < subdomainNum; i++) {
    solveSubdomain(i, lambda, x_i);
}

// 错误示例
for(int i=0;i<subdomainNum;i++){
    solveSubdomain(i,lambda,x_i);
}
```

## 代码长度限制
- 单个函数代码长度不超过 100 行
- 单个类的头文件代码不超过 500 行
- 过长函数拆分为子函数，过大类拆分为多个相关类（遵循单一职责原则）

## 异常处理
- 禁止使用 try-catch 捕获所有异常
- 仅针对特定场景（如文件读取、库调用失败）捕获异常
- 捕获后通过 spdlog 输出 ERROR 级日志，明确异常原因，避免程序崩溃