
library(Seurat)
library(ggplot2)
library(cols4all)
library(tidydr)
library(openxlsx)
library(tidyr)
script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)

my_theme <- theme_void()+theme(plot.margin = margin(5.5, 15, 5.5, 5.5))+
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, units = "inches"), ends = "last", type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 15, color = "black"))

col_sample <- c("191598A1" = "#49B1DD", "191598A2" = "#E98DAF", "191598A3" = "#E9CF91", "191598B1" = "#95C6C6",
                "191598B2" = "#FFF4F2", "191598B3" = "#2A9D8E", "191598C1" = "#FEB705", "191598C2" = "#F66F69",
                "191598C3" = "#E73F74","191598D1" = "#FFBFFD", "191598D2" = "#EFFDAF", "191598D3" = "#EF5FFF")
col_clinic <- c("PT" = "#FFC24B", "LT" = "#1597A5")
#####################################################################
seurat_dr <- readRDS("../5-harmony移除批次效应/scRNA_harmony.rds")
dim(seurat_dr)
#23291 83229

#特征选择筛选高变基因
head(VariableFeatures(seurat_dr), 20)
#####################################################################
#PCA降维
#seurat_dr <- RunPCA(seurat_dr, npcs = 80)

#量化拐点识别
pct <- seurat_dr[["pca"]]@stdev / sum(seurat_dr[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], 
              sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

#绘制ElbowPlot并可视化拐点
p <- ElbowPlot(seurat_dr, ndims = 80)$data %>% ggplot() + 
  geom_point(aes(x = dims, y = stdev)) + 
  geom_vline(xintercept = pc.use, color = "#E73F74") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach") + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.ticks = element_line(linewidth = 0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave(plot = p, filename = "ElbowPlot.png", width = 6, height = 6, dpi = 300)
ggsave(plot = p, filename = "ElbowPlot.pdf", width = 6, height = 6, dpi = 300)

pc.num <- 1:pc.use
#-------------------------------------------------------------------
#t-SNE降维可视化
seurat_dr <- RunTSNE(seurat_dr, reduction = "harmony", dims = pc.num,
                     perplexity = 30, max_iter = 1000,
                     theta = 0.5, eta = 200, num_threads = 4)
p <- DimPlot(seurat_dr, reduction = "tsne", group.by = "Samples", shuffle = T,
             cols = col_sample, raster = F, pt.size = 0.5) + my_theme
ggsave(plot = p, filename = "t-SNE.pdf", width = 9.5, height = 8, dpi = 300)
ggsave(plot = p, filename = "t-SNE.png", width = 9.5, height = 8, dpi = 300)

#UMAP降维可视化
seurat_dr <- RunUMAP(seurat_dr, reduction = "harmony", dims = pc.num,
                     n.components = 2, n.neighbors = 30, n.epochs = 200,
                     min.dist = 0.3, learning.rate = 1, spread = 1)
q <- DimPlot(seurat_dr, reduction = "umap", group.by = "Samples", shuffle = T,
             cols = col_sample, raster = F, pt.size = 0.5) + my_theme
ggsave(plot = q, filename = "UMAP.pdf", width = 9.5, height = 8, dpi = 300)
ggsave(plot = q, filename = "UMAP.png", width = 9.5, height = 8, dpi = 300)
################################################################
#将分群进行聚类
resolution <- c(0.2, 0.4, 0.6, 0.8, 1)
reduction <- "harmony"
assay <- "RNA"

for (res in resolution) {
  tmp_dr <- FindNeighbors(seurat_dr, dims = pc.num, reduction = reduction, features = VariableFeatures(seurat_dr),
                          assay = assay, verbose = T, graph.name = paste(reduction, assay, sep = "."))
  tmp_dr <- FindClusters(tmp_dr, resolution = res, verbose = T, algorithm = 2,
                         graph.name = paste(reduction, assay, sep = "."))
  print(table(tmp_dr@meta.data$seurat_clusters))
  q <- DimPlot(tmp_dr, reduction = "umap", group.by = "seurat_clusters",
               label.box = T, label = T, shuffle = T, raster = F,
               pt.size = 0.3) + my_theme
  ggsave(plot = q, filename = paste0(res, '_umap_seurat_res.png'), width = 9, height = 8, dpi = 300)
  rm(tmp_dr)
}
#----------------------------------
#确定好resolution值后最终执行
res <- 0.4
seurat_dr <- FindNeighbors(seurat_dr, dims = pc.num, reduction = reduction, features = VariableFeatures(seurat_dr),
                           assay = assay, verbose = T, graph.name = paste(reduction, assay, sep = "."))
seurat_dr <- FindClusters(seurat_dr, resolution = res, verbose = T, algorithm = 2,
                          graph.name = paste(reduction, assay, sep = "."))
table(seurat_dr$seurat_clusters)
#------------------------------------------------------
table(Idents(seurat_dr))

cols2 <- c4a("classic20", 17)

p <- DimPlot(seurat_dr, reduction = "umap", group.by = "seurat_clusters",
             label.box = T, label = T, shuffle = T, raster = F,
             cols = cols2, pt.size = 0.3) + my_theme
ggsave(plot = p, filename = paste0('umap_seurat_res.png'), width = 9, height = 8, dpi = 300)
ggsave(plot = p, filename = paste0('umap_seurat_res.pdf'), width = 9, height = 8, dpi = 300)

q <- DimPlot(seurat_dr, reduction = "tsne", group.by = "seurat_clusters",
             label.box = T, label = T, shuffle = T, raster = F,
             cols = cols2, pt.size = 0.3) + my_theme
ggsave(plot = q, filename = paste0('tsne_seurat_res.png'), width = 9, height = 8, dpi = 300)
ggsave(plot = q, filename = paste0('tsne_seurat_res.pdf'), width = 9, height = 8, dpi = 300)
#------------------------------------------------------
saveRDS(seurat_dr, "seurat_dr.rds")

#计算每个cluster中的特异表达基因
markers <- FindAllMarkers(seurat_dr, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(markers, file = "markers_genes.xlsx", overwrite = T, asTable = T, rowNames = F)

