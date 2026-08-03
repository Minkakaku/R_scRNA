# ==================== 拟时序分析 - 方法3：Seurat 差异表达基因 ====================
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)
library(tibble)

# 创建输出目录
monocle_outdir <- "Monocle_Epi_method3"
dir.create(monocle_outdir, showWarnings = FALSE)

# 读取 Seurat 对象
scEpi <- readRDS("Epi_subset.rds")
Idents(scEpi) <- "celltype"
cat("细胞类型：", unique(scEpi$celltype), "\n")

# 提取 counts 矩阵
count <- GetAssayData(scEpi, layer = "counts")
if (is.null(count)) stop("counts 层不存在，请确保 Seurat 对象包含原始计数矩阵。")

# 准备注释数据框
fdata <- data.frame(gene_short_name = rownames(count), row.names = rownames(count))
pdata <- scEpi@meta.data

# 构建 Monocle 对象
cds <- newCellDataSet(
  count,
  featureData = new("AnnotatedDataFrame", data = fdata),
  phenoData = new("AnnotatedDataFrame", data = pdata),
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

# 标准化和离散度估计
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# 添加 UMI 总数列
pData(cds)$UMI <- Matrix::colSums(exprs(cds))

# 定义颜色向量
n_types <- length(unique(scEpi$celltype))
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")
if (n_types > length(plan)) {
  plan <- colorRampPalette(plan)(n_types)
}

# ---------- 方法3：使用 Seurat 差异表达基因 ----------
# 计算各细胞类型间的差异基因（与 Seurat 聚类一致）
deg.cluster <- FindAllMarkers(scEpi, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
express_genes <- unique(subset(deg.cluster, p_val_adj < 0.05)$gene)
cat("筛选出的差异基因数量：", length(express_genes), "\n")

cds <- setOrderingFilter(cds, express_genes)

# 绘制排序基因图
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

# 降维与排序
cds <- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)

# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "celltype")
pdf(file.path(monocle_outdir, "monocle_cluster-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "celltype") + 
  facet_wrap(~celltype, nrow = 2) +
  scale_color_manual(values = plan)
pdf(file.path(monocle_outdir, "monocle_facet-celltype-3.pdf"), height = 7, width = 12)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

# 保存 cds 对象
saveRDS(cds, file.path(monocle_outdir, "cds_deg-3.rds"))

cat("方法3（差异表达基因）运行完成，结果保存在：", monocle_outdir, "\n")
