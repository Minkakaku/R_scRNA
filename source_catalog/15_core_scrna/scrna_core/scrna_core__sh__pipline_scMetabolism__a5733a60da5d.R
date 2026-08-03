library(Seurat)
# library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)
# library(clusterProfiler)
library(msigdbr)
library(UCell)
library(AUCell)

setwd("/mnt/NexData_file/home/sh/Projects/SC2510235895/18.Mission/Main/")
S_subf<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/cellchat_sub_cell_type.rds")
# DimPlot(S_subf,reduction='umap',label = TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt",sep="\t",header = T)
lines <- readLines("/mnt/NexData_file/home/sh/Projects/SC2510235895/18.Mission/script/all.gmt")
geneset <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset[[set_name]] <- genes
  }
}
# geneset<-geneset$gene
# gene<-list(geneset)
adata<-AddModuleScore_UCell(S_subf,features=geneset)
colnames(adata@meta.data)

# 假设基因集名称存储在向量中
colnames(adata@meta.data) <- gsub(" / ", "_", colnames(adata@meta.data))
gene_set_names <- colnames(adata@meta.data)[16:27]  # 85个基因集的名称

# 创建输出目录
dir.create("./violin_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./umap_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./dot_plots", recursive = TRUE, showWarnings = FALSE)
# 循环处理每个基因集
adata$group <- factor(adata$group,levels = c("asham","aday1","aday3","aday7","aday14"))
for (i in seq_along(gene_set_names)) {
  score_col <- gene_set_names[i]
  if (!score_col %in% colnames(adata@meta.data)) {
    warning(paste("Column not found:", score_col, "- skipping"))
    next
  }
  cat("Processing", i, "/", length(gene_set_names), ":", score_col, "\n")
  # 绘制小提琴图
  p <- VlnPlot(adata, features = score_col, pt.size = 0, group.by = 'group') + NoLegend()+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.7),outlier.shape = NA)
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.pdf"), 
         width = 10, height = 10, bg = 'white')
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.png"), 
         width = 10, height = 10, bg = 'white')
  # 绘制 UMAP 图
  p <- FeaturePlot(adata, features = score_col, reduction = 'umap', 
                   label = TRUE, max.cutoff = 0.99)
  ggsave(paste0("./umap_plots/", score_col, "_umap.png"), 
         width = 7, height = 6, bg = 'white')
  ggsave(paste0("./umap_plots/", score_col, "_umap.pdf"), 
         width = 7, height = 6, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'group', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.png"), 
         width = 7, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.pdf"), 
         width = 7, height = 4, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'cell_type', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.png"), 
         width = 30, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.pdf"), 
         width = 30, height = 4, bg = 'white')
}