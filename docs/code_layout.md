# 代码分层

`src/R/00_runtime.R` 到 `12_visualization.R` 是底层通用函数；`src/R/steps/` 将它们组合为分析点；`projects/<项目>/workflow/` 只记录某个项目实际调用的顺序。

新增内容时按下面判断：

- 换一个项目仍成立：放 `src/R/`。
- 所有项目共用的参数或 marker：放根目录 `config/`。
- 只属于当前项目的分组、阈值、cluster、marker 或人工判断：放项目 `config/` 或新的 workflow 步骤。
- 中间对象放项目 `data/processed/`；最终图表放项目 `results/`。

不要把数据对象、图片、表格或内部绝对路径放进 `src/R/`。
