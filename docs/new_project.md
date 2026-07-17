# 新项目指南

## 标准分析顺序

| Workflow 文件 | 调用方法 | 主要输出 |
|---|---|---|
| `step01.QC_before_filter.R` | `step_qc_before_filter()` | QC 图、逐样本阈值表、合并对象 |
| `step02.QC_filter.R` | `step_qc_filter()` | 过滤后对象、保留率、样本分组表 |
| `step03.Global_clustering.R` | `step_global_clustering()` | PCA、批次校正、聚类、UMAP |
| `step04.Markers.R` | `step_find_markers()` | 全部 marker、注释模板 |
| `step05.Annotation.R` | `step_cell_annotation()` | 细胞类型注释对象与图表 |
| `step06.Proportion.R` | `step_cell_proportion()` | 细胞比例表与图 |
| `step07.Subcluster.R` | `step_subcluster_clustering()` | 指定细胞类型子聚类 |
| `step08.Subcluster_markers.R` | `step_subcluster_find_markers()` | 子聚类 marker 与注释模板 |
| `step09.Subcluster_annotation.R` | `step_subcluster_annotation()` | 子类型注释对象与图表 |
| `step10.Subcluster_proportion.R` | `step_subcluster_proportion()` | 子类型比例表与图 |

## 每个 workflow 的固定骨架

```r
rm(list = ls())
setwd("你的项目根目录")
source("config/config.R")
source("src/load_all.R")

# 本次只调用一个清楚的分析方法
step_global_clustering()
```

如需改变 PC 数、resolution 或分组，先修改 `config/config.R`，再新建步骤文件记录这次运行，例如 `step03b.Recluster_30PC.R`。不要覆盖已经运行过的 `step03.Global_clustering.R`。

## 数据布局

代码会自动发现以下两种 Cell Ranger 表达矩阵：

```text
data/<sample>/outs/filtered_feature_bc_matrix/
data/<sample>/outs/filtered_feature_bc_matrix.h5
```

`<sample>` 会成为 Seurat 对象中的 `orig.ident`。
