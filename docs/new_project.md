# 新项目指南

## 标准分析顺序

| Workflow 文件 | 调用方法 | 主要输出 |
|---|---|---|
| `01_QC_before_filter.R` | `step_qc_before_filter()` | QC 图、逐样本阈值表、合并对象 |
| `02_QC_filter.R` | `step_qc_filter()` | 过滤后对象、保留率、样本分组表 |
| `03_Global_clustering.R` | `step_global_clustering()` | PCA、批次校正、聚类、UMAP |
| `04_Markers.R` | `step_find_markers()` | 完整 marker、注释模板 |
| `05_Annotation.R` | `step_cell_annotation()` | 细胞类型注释对象与图表 |
| `06_Proportion.R` | `step_cell_proportion()` | 样本级细胞比例表与图 |
| `07_Subcluster.R` | `step_subcluster_clustering()` | 指定大类细胞子聚类 |
| `08_Subcluster_markers.R` | `step_subcluster_find_markers()` | 子聚类 marker 与注释模板 |
| `09_Subcluster_annotation.R` | `step_subcluster_annotation()` | 子类型注释对象与图表 |
| `10_Subcluster_proportion.R` | `step_subcluster_proportion()` | 子类型比例表与图 |

## 固定入口

从仓库根目录运行项目脚本。每个脚本只加载项目并调用一个方法：

```r
source("src/R/project_config.R")
load_project("projects/project1")
step_global_clustering()
```

如需改变 PC 数或 resolution，修改项目 `config/config.yaml` 后新建下一个编号脚本；不要覆盖已经执行的步骤。

## 输入布局

代码会在项目 `data/` 下递归发现以下 Cell Ranger 矩阵，因此推荐放到 `data/raw/<sample>/outs/`：

```text
filtered_feature_bc_matrix/
filtered_feature_bc_matrix.h5
```

`<sample>` 会成为 Seurat 对象中的 `orig.ident`。
