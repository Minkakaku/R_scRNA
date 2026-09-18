# 方法网页与论文依据

- Seurat PBMC 标准流程：https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
- Seurat v5 integration：https://satijalab.org/seurat/articles/seurat5_integration
- Seurat Stdev/ElbowPlot：https://satijalab.org/seurat/reference/elbowplot
- PCAtools findElbowPoint：https://bioconductor.org/packages/PCAtools
- nFactors nScree：https://cran.r-project.org/package=nFactors
- scran denoisePCA：https://bioconductor.org/packages/scran
- Seurat multimodal/integration 论文：PMID 34062119

代码决策：

1. QC 阈值保留在 config，不写死在复用函数。
2. PCAtools 和 nFactors 同时给出几何拐点，仅作为 dims 参考，不自动覆盖项目配置。
3. marker 和注释使用完整生物学标签，避免技术 cluster 直接成为结论。
4. 细胞比例及 CON/EXP 统计以样本为重复，不以细胞为重复。
5. 最终图成对保存 PNG/PDF，支撑表和中文说明与图同目录。