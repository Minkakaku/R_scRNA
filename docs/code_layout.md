# 代码分层

`src/00_runtime.R` 到 `src/12_visualization.R` 是底层通用函数；`src/steps/` 把底层函数组合成可直接调用的分析点；`workflow/` 记录当前项目在哪一天、以什么顺序实际调用了哪些分析点。

新增代码时按下面判断：

- 换一个项目仍然成立：放 `src/`。
- 只与当前项目的 cluster、marker、对象或人工判断有关：放新的 workflow 步骤。
- 只是参数变化：放 `config/config.R`。

不要把生成的数据对象、图片或表格放进 `src/`。
