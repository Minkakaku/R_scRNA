# 公共配置

- `default.yaml`：跨项目通用的 QC、归一化、聚类、marker 和子聚类默认值。
- `reference.yaml`：公共大类细胞 marker 参考。

每个项目通过 `projects/<项目>/config/config.yaml` 递归覆盖同名参数。样本分组、逐样本 QC 阈值和人工注释属于项目状态，分别保存在该项目的 `samples.tsv`、`qc_thresholds.tsv`、`cluster_annotation.tsv` 和 `subcluster_annotation.tsv`，不要写入根目录公共配置。
