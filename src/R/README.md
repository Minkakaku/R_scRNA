# R 复用代码

## 通用函数

- `00_runtime.R`：日志、QS2 读写、NPG 配色、非栅格 UMAP、紧凑 DotPlot 和 600 dpi PNG/矢量 PDF 保存。
- `01_sample_discovery.R`：发现 Cell Ranger 表达矩阵。
- `02_sample_metadata.R`：样本分组表和 Seurat metadata。
- `03_quality_control.R`、`04_qc_filter.R`：QC 指标、逐样本阈值和过滤。
- `05_normalization.R`、`06_batch_correction.R`：NormalizeData、SCTransform 和批次校正。
- `07_clustering.R`：PCA 自动拐点、邻接图、聚类和 UMAP。
- `08_markers.R`、`09_annotation.R`：marker 与人工注释。
- `10_cell_proportion.R`：样本级组成分析。
- `11_subcluster.R`：细胞大类提取和子聚类。
- `12_visualization.R`：比例图、Sankey 和热图。
- `project_config.R`：合并公共/项目 YAML，并建立项目输入输出路径。
- `load_all.R`：按编号加载全部通用函数和 `steps/`。

`steps/` 将函数组合成可直接调用的标准分析点，但不会自动执行。具体项目的入口只放在 `projects/<项目>/workflow/`；不要在复用层保存样本名、绝对路径、固定 cluster 注释或一次性分析判断。

已有对象处理：`steps/12_existing_rds_overview.R` 检查 Seurat RDS 并绘制原始 UMAP；`steps/13_merge_existing_rds.R` 校验重叠细胞、保留一次并输出标准 QS，输入和输出路径均来自项目配置。

已有结果重绘：`steps/15_refresh_existing_figures.R` 按项目提供的对象、marker 和图名重绘 QC、PCA、UMAP、DotPlot、样本/组间组成、ROE、缓存 AUCell 与 PCD 小提琴。保留细胞、注释、坐标和原有评分/统计，不重存输入、不生成额外支撑表；组间组成先计算样本比例再等权平均。
