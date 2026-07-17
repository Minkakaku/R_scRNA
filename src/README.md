# src

这里是可跨项目复用的单细胞分析代码库，可以单独上传到 GitHub。

## 通用函数

- `00_runtime.R`：日志、QS2 读写和图片保存。
- `01_sample_discovery.R`：发现 Cell Ranger 表达矩阵。
- `02_sample_metadata.R`：样本分组表和 Seurat 元数据。
- `03_quality_control.R`：QC 指标、双细胞、去污染和 QC 图。
- `04_qc_filter.R`：逐样本阈值表和细胞过滤。
- `05_normalization.R`：NormalizeData、SCTransform 和 SCT integration。
- `06_batch_correction.R`：Harmony、CCA 和批次校正选择。
- `07_clustering.R`：PCA 维度、邻接图、聚类和 UMAP。
- `08_markers.R`：marker 计算和 marker 汇总。
- `09_annotation.R`：人工注释模板、注释应用和 cluster 清理。
- `10_cell_proportion.R`：细胞比例统计和输出。
- `11_subcluster.R`：细胞类型提取和子聚类。
- `12_visualization.R`：比例图、Sankey 和热图。
- `load_all.R`：按依赖顺序加载全部通用代码。

## 可复用分析步骤

`steps/` 把上述函数组合成标准分析点，但不会自动执行：

- `01_qc_before_filter.R`
- `02_qc_filter.R`
- `03_global_clustering.R`
- `04_markers.R`
- `05_annotation.R`
- `06_proportion.R`
- `07_subcluster_clustering.R`
- `08_subcluster_markers.R`
- `09_subcluster_annotation.R`
- `10_subcluster_proportion.R`
- `11_stepwise.R`

项目运行脚本不要写进 `src/`。实际运行时在项目的 `workflow/` 中新建 `stepXX.名称.R`，文件只需要加载配置和 `src/load_all.R`，然后调用当前分析步骤。
