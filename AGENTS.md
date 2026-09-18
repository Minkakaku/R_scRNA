# 项目代码约定

本仓库是多步骤 R/Seurat 主流程，修改时遵守以下规则。

## 分层

- 可跨项目复用的函数放 `src/R/`，组合分析点放 `src/R/steps/`。
- 公共默认参数和 marker 参考放根目录 `config/`。
- 样本、分组、QC 阈值、比较、人工注释和一次性决定放 `projects/<项目>/config/`。
- 原始输入、中间对象、最终结果和日志分别放项目的 `data/raw`、`data/processed`、`results` 和 `logs`。
- `source_catalog/` 是来源档案，不作为正式 workflow 入口。

## Workflow

- workflow 文件按实际执行顺序编号，只加载项目并调用一个清楚的分析步骤。
- 已运行脚本不覆盖；重跑、改参数或补分析时新增编号脚本。
- 不用 notebook 作为正式入口，不建立跨越 QC 或人工注释检查点的 `run_all`。
- 不把项目路径、项目编号、样本名、固定 cluster 注释或一次性生物学判断写进 `src/R/`。

## 输出

- QS2 中间对象按步骤编号写入项目 `data/processed/`。
- 最终图表写入项目 `results/<编号_分析点>/`，图默认同名 PNG/PDF。
- 有解释价值的图保留支撑 TSV/CSV；正式项目补充 `README_results_CN.md`。
- 保持代码简洁，优先复用已有函数，不堆叠无必要的防御性封装。
