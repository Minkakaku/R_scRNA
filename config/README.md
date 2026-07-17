# config

`config.R` 保存项目级路径和参数。新项目开始时优先检查：

- `GROUP_LEVELS` 与 `PROPORTION_COMPARISONS`；
- `QC_READ` 与 `QC_DEFAULTS`；
- `NORMALIZE` 与 `CLUSTER`；
- `SUBCLUSTER$target_celltype`。

以下 CSV 由分析步骤生成，再由人工检查或填写：

- `qc_thresholds.csv`：逐样本 QC 阈值；
- `sample_groups.csv`：样本与分组；
- `cell_annotation.csv`：主聚类注释；
- `subcluster_annotation.csv`：子聚类注释。

这些 CSV 属于当前项目。建议在具体项目仓库中提交，用来记录 QC、分组和人工注释决策；如果样本名涉及敏感信息，再按项目需要单独忽略。
