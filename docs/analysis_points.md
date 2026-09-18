# R/Seurat 单细胞主流程分析点

| 顺序 | 正式代码 | 分析点 | 主要输出 |
|---|---|---|---|
| 01 | `src/R/01_sample_discovery.R`、`02_sample_metadata.R` | 样本发现和 metadata | 输入审计表 |
| 02 | `src/R/03_quality_control.R`、`04_qc_filter.R` | QC 指标、阈值和过滤 | 过滤前后图表与 QS2 |
| 03 | `src/R/05_normalization.R`、`06_batch_correction.R` | NormalizeData/SCTransform 与批次校正 | 标准化对象 |
| 04 | `src/R/07_clustering.R` | PCA、自动拐点、邻居图、聚类、UMAP | PCA 表、UMAP、聚类对象 |
| 05 | `src/R/08_markers.R` | cluster/cell-type marker | 完整 marker 表和 DotPlot |
| 06 | `src/R/09_annotation.R` | 人工生物学注释 | 注释表、UMAP、注释对象 |
| 07 | `src/R/10_cell_proportion.R` | 样本级细胞比例 | 样本比例表和组间图 |
| 08 | `src/R/11_subcluster.R` | 指定大类的亚群重分析 | 亚群聚类、marker、注释和比例 |
| 09 | `src/R/12_visualization.R` | 通用结果可视化 | 同名 PNG/PDF 与绘图表 |

复用实现位于 `src/R/`；每个项目通过自己的 `workflow/` 记录实际执行历史。
