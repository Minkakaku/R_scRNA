library(Seurat)
library(tidyverse)
library(ConsensusClusterPlus)
source("sc_analysis_tools.R")

out_dir <- 'Epi_consensus_clustering_4slides_sample_filtered'
dir.create(out_dir)

all_data1 <- readRDS("~/R_OUTPUT/data_for_analysis/all_data_filtered2.rds")
all_data2 <- readRDS("~/R_OUTPUT/all_data_lubin.rds")
all_data2$annotation <- all_data2$annotation_re
all_data <- merge(all_data1, all_data2)
all_data$annotation %>% table()
all_data$sample %>% table()




# Tumor annotation -------------
all_data1@meta.data %>% 
  dplyr::filter(res_0.6 %in% c(1, 2, 5, 6, 12, 13, 14, 15, 16)) %>% 
  group_by(sample, res_0.6) %>% 
  summarize(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  ggplot(aes(x = sample, y = freq, fill = res_0.6)) + 
  geom_col()

## spatial plot -----------------------------------------------------------------
spatial_theme <- 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
        #axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.title = element_text(hjust = 0.5, size = 20))


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')
all_data_meta <- all_data1[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!sample %in% all_data_meta$sample) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_cellType.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  




## highlight -----------
### Tumor ----------
all_data1$annotation4 <- fct_collapse(
  all_data1$res_0.6, 
  'Others' = unique(all_data1$res_0.6)[!unique(all_data1$res_0.6) %in% c(1, 2, 5, 6, 12, 13, 14, 15, 16)]
)
all_data1$annotation4 %>% table()
all_data_meta <- all_data1[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!sample %in% all_data_meta$sample) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation4, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey", allcolour)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_epi highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  




Epithelial <- subset(all_data, annotation == 'Epithelial')

Major_marker_ls <- c(
  "Epithelial" = "EPCAM", 
  "Epithelial" = "KRT5", 
  "Epithelial" = "KRT14", 
  "Epithelial" = "KRT20", 
  
  "Endothelial" = "PECAM1", 
  "Endothelial" = "VWF", 
  "Endothelial" = "AQP1", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2", 
  
  "B cell" = "CD19", 
  "B cell" = "CD79A", 
  "B cell" = "MS4A1", 
  "Plasma" = "TNFRSF17", 
  "Plasma" = "JCHAIN", 
  "Plasma" = "IGJ", 
  "Plasma" = "MZB1", 
  "Plasma" = "XBP1", 
  
  "T cell" = "CD3D", 
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "TRAC", 
  "T cell" = "BCL11B", 
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  "Proliferating" = "MKI67", 
  "Proliferating" = "TOP2A", 
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  # "Pericyte" = "PLVAP", 
  
  "Monocyte" = "CD14", 
  "Monocyte" = "FCGR3A", 
  "Monocyte" = "ITGAM", 
  
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  "Macrophage" = "C1QA", 
  "Macrophage" = "C1QC", 
  "Macrophage" = "CSF1R", 
  "Macrophage" = "IL1B", 
  
  "Dendritic cell" = "CD83", 
  "Dendritic cell" = "LILRA4", 
  "Dendritic cell" = "CLEC9A", 
  "Dendritic cell" = "IDO1", 
  "Dendritic cell" = "CD1C", 
  "Dendritic cell" = "FCER1A", 
  "Dendritic cell" = "CST3", 
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  "Lamp3_DC" = "LAMP3", 
  
  "Netrophil" = "CSF3R", 
  "Netrophil" = "FCGR3B", 
  "Netrophil" = "S100A8", 
  "Netrophil" = "S100A9", 
  "Netrophil" = "G0S2", 
  "Netrophil" = "RETNLG", 
  
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = "HDC", 
  "Mast cell" = "TPSAB1", 
  "Mast cell" = "CTSG", 
  "Mast cell" = "CPA3", 
  
  "Smooth muscle cell" = "PLN", 
  "Smooth muscle cell" = "MYLK", 
  "Smooth muscle cell" = "MYH11", 
  "Smooth muscle cell" = "ACTG2", 
  "Smooth muscle cell" = "TAGLN", 
  
  "Schwann cell" = "MPZ", 
  "Schwann cell" = "SOX10", 
  
  "Neuron" = "TUBB3", 
  "Neuron" = "SNAP25", 
  
  "Mesothelial" = "MSLN", 
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2"
)


DefaultAssay(Epithelial) <- 'RNA'

Epithelial <- sc_normalize(Epithelial, scale.factor = "median")

## sketch data --------
# DefaultAssay(Epithelial) <- 'RNA'
# Epithelial[['sketch']] <- NULL
# Epithelial$leverage.score <- NULL

Epithelial <- 
  SketchData(
    Epithelial,
    ncells         =  10000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    over.write     = FALSE,
    seed = 123L
  )
# normalize  again
Epithelial <- sc_normalize(Epithelial, scale.factor = "default", nfeatures = 1500)


Epithelial <- sc_pca(
  Epithelial, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_epi <- 10

Epithelial <- sc_dr(
  Epithelial, 
  ndims = ndims_epi,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr"),
  dir.create.force = TRUE
)


Epithelial <- sc_cluster(
  Epithelial, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

epi_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  Enterocyte  = "SLC26A3",
  
  Goblet = "MUC2",
  Goblet = "TFF3",
  
  
  # Enteroendocrine = "CHGA",
  Enteroendocrine = "TPH1",
  Enteroendocrine = "NEUROD1",
  
  'Paneth cell' = 'DEFA6',
  'Paneth cell' = 'IL9',
  
  'Transit amplifying cell' = 'CCNA2',
  'Transit amplifying cell' = 'HMGB3',
  'Transit amplifying cell' = 'MCM5',
  'Transit amplifying cell' = 'PCNA',
  'Transit amplifying cell' = 'TOP2A',
  
  'Proximal Mature Enterocytes' = 'FABP2',
  'Proximal Mature Enterocytes' = 'ACE',
  'Proximal Mature Enterocytes' = 'FABP6',
  'Proximal Mature Enterocytes' = 'APOC3',
  
  'Proximal Enterocytes' = 'HEBP1',
  'Proximal Enterocytes' = 'CCL25'
)

Epithelial <- sc_dotplot(
  Epithelial, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  dir.create.force = TRUE
)
Epithelial <- sc_dotplot(
  Epithelial, 
  ndims = ndims_epi,
  marker = Major_marker_ls,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_epi <- 
  sc_findmarker(
    Epithelial,
    idents = "sketch_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(Epithelial, all_markers = all_markers_epi, sort_by = "pct")



dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
  height = 6,
  width  = 12
)
sc_mark_heatmap(
  Epithelial,
  markers = all_markers_epi,
  idents = "sketch_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  Epithelial,
  markers = all_markers_epi,
  idents = "sketch_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)




# project  data
DefaultAssay(Epithelial) <- "sketch"
Epithelial[['full.harmonyumap']] <- NULL
Epithelial <-
  Epithelial %>%
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",

    sketched.assay = "sketch",
    sketched.reduction = "harmony", # "pca",

    umap.model = "harmonyumap",
    dims = 1:ndims_epi,

    refdata = list(
      annotation_full = "annotation",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(Epithelial) <- 'annotation_full'
DefaultAssay(Epithelial) <- "RNA"

pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

saveRDS(Epithelial, file.path(out_dir, 'epi', 'Epithelial.rds'))



Epithelial@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(RNA_snn_res.0.4, sample) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = sample, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )
ggsave(file.path(out_dir, 'epi', 'Epi_cell_prop_res_0.4.pdf'), width = 8, height = 6)




Epithelial$sample %>% table()
Epithelial2 <- Epithelial[, !(Epithelial$sample %in% c('ATC12_af_T_1', 'ATC12_af_T_2'))]
Epithelial2$sample %>% table()
Epithelial <- Epithelial2




Epithelial_ls <- SplitObject(JoinLayers(Epithelial), split.by = 'sample')

expr_df <- 
  Epithelial_ls %>% 
  imap(function(Epithelial, sample){ # browser()
    expr_matrix <- LayerData(JoinLayers(Epithelial), 'data')
    expr <- rowMeans(expr_matrix) %>% as.data.frame()
    colnames(expr) <- sample
    expr
  }) %>% 
  bind_cols()


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb


out_dir <- 'Epi_consensus_clustering_4slides_sample_filtered'
dir.create(out_dir)
setwd(out_dir)

set.seed(12345)
results <- ConsensusClusterPlus(
  as.matrix(expr_df), 
  maxK = 10, 
  reps = 1000, 
  plot = "png", 
  seed = 12345
)
setwd('~/')

saveRDS(results, 'Epi_consensus_clustering_4slides_sample_filtered/result.rds')


# 中位数中心化 -----------
#sweep函数减去中位数进行标准化
exprSet = sweep(expr_df,1, apply(expr_df,1,median,na.rm=T))
par(mfrow = c(1,2))
boxplot(expr_df[,1:20],main = "before")
boxplot(exprSet[,1:20],main = "after")

setwd(out_dir)
set.seed(12345)
results <- ConsensusClusterPlus(
  as.matrix(exprSet), 
  maxK = 10, 
  reps = 1000, 
  plot = "png", 
  seed = 12345
)
setwd('~/')




# 使用PAC标准去筛选K值 ------------
maxK = 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK## [1] 4




#输出K=5时的一致性矩阵
results[[3]][["consensusMatrix"]][1:5,1:5]

# colnames(consensusMatrix_select) <- colnames(exprSet)
# rownames(consensusMatrix_select) <- colnames(exprSet)

pheatmap(results[[3]][["consensusMatrix"]],
         cluster_rows = F,
         cluster_cols = TRUE,
         # annotation_row = annotation,
         # scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = T,
         show_rownames = T)

#hclust聚类树
results[[3]][["consensusTree"]]
#Call:
#hclust(d = as.dist(1 - fm), method = finalLinkage)
#Cluster method   : average 
#Number of objects: 128 

#查看样本所属的聚类群
table(results[[3]]$consensusClass)
results[[3]][["consensusClass"]]




# 计算聚类一致性 (cluster-consensus，CLC) 和样品一致性 (item-consensus，IC) -----------
icl <- calcICL(results, title = "./ConsensusCluster",plot = "png", writeTable = TRUE)
## 返回了具有两个元素的list
dim(icl[["clusterConsensus"]])#[1] 20  3
icl[["clusterConsensus"]] 




# plot ------------
## 热图 ------------
message('heatmap')
library(pheatmap)
library(RColorBrewer)

# 准备热图数据
heatmap_data <- expr_df[order(results[[3]][["consensusClass"]])] %>% t()
cluster_order <- results[[3]][["consensusClass"]][order(results[[3]][["consensusClass"]])]

# 创建注释
annotation <- data.frame(Cluster = factor(cluster_order))
rownames(annotation) <- rownames(heatmap_data)

# 绘制热图
pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k3.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = TRUE,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()


# 绘制热图
pdf(file.path(out_dir, 'Clustering Heatmap_k3.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Clustering Heatmap",
         fontsize_row = 8,
         treeheight_row = 0,
         show_colnames = F,
         show_rownames = T)
dev.off()




## 簇间距离矩阵可视化 --------------
# 计算簇中心距离
center_dist <- dist(results[[3]][["consensusClass"]])
center_dist_matrix <- as.matrix(center_dist)

# 绘制簇间距离热图
pdf(file.path(out_dir, 'Distance Between Cluster Centers_k3.pdf'), height = 8, width  = 11)
pheatmap(center_dist_matrix,
         main = "Distance Between Cluster Centers",
         color = colorRampPalette(c("white", "red"))(50),
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()




## PCA ------------
# 如果特征维度较高，先用PCA降维
pca_result <- prcomp(t(expr_df), scale. = TRUE)

# PCA结果可视化
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$cluster <- factor(results[[3]][["consensusClass"]])

ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(title = "Consensus Clustering on PCA-reduced Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, 'Consensus Clustering on PCA-reduced Data_k3.pdf'), height = 6, width  = 11)


# 标记离群点
pct = 0.8
pca_data2 <- pca_data %>%
  group_by(cluster) %>%
  mutate(
    # 计算聚类中心
    center_x = mean(PC1),
    center_y = mean(PC2),
    # 计算到聚类中心的欧氏距离
    dist_to_center = sqrt((PC1 - center_x)^2 + (PC2 - center_y)^2)
  ) %>%
  ungroup() %>% 
  as.data.frame()
rownames(pca_data2) <- rownames(pca_data)

# 定义离群点阈值（例如：距离大于第95百分位数）
pca_data3 <- pca_data2 %>%
  group_by(cluster) %>%
  mutate(
    is_outlier = dist_to_center > quantile(dist_to_center, pct)
  ) %>%
  ungroup() %>% 
  as.data.frame()
rownames(pca_data3) <- rownames(pca_data2)


# 绘制带离群点标记的图
ggplot(pca_data3, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = pct) +
  # # 突出显示离群点
  # geom_point(data = subset(pca_data3, is_outlier), 
  #            aes(x = PC1, y = PC2, color = cluster),
  #            size = 5, shape = 1, stroke = 1.5) +
  # 添加文本标签
  geom_text(data = subset(pca_data3, is_outlier),
            aes(x = PC1, y = PC2, label = rownames(pca_data3[pca_data3$is_outlier, ])),
            vjust = -0.5, size = 3) +
  labs(title = "Consensus Clustering on PCA-reduced Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, 'Consensus Clustering on PCA-reduced Data_k3_outlier.pdf'), height = 6, width  = 11)




# DEGs ----------- 
consensus_cluster <- str_c('cluster_', results[[3]][["consensusClass"]])
names(consensus_cluster) <- names(results[[3]][["consensusClass"]])
Epithelial$consensus_cluster <- str_replace_all(Epithelial$sample, consensus_cluster)
Idents(Epithelial)<-'consensus_cluster'
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Epithelial_markers, file.path(out_dir, 'Epithelial_markers_consensus_clusters_k3.csv'), row.names = F)

Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

## heatmap ---------
# 准备热图数据
heatmap_data <- expr_df[unique(unlist(Epithelial_de_genes_for_enrichmenet)), order(results[[3]][["consensusClass"]])] %>% t()
cluster_order <- results[[3]][["consensusClass"]][order(results[[3]][["consensusClass"]])]

# 创建注释
annotation <- data.frame(Cluster = factor(cluster_order))
rownames(annotation) <- rownames(heatmap_data)

# 绘制热图
pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k3_DEGs.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = TRUE,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()

row_order <- 
  annotation %>%
  rownames_to_column("row") %>%
  with(split(row, Cluster)) %>%
  map(function(x){
    dat <- heatmap_data[x, ]  
    cl <-
      dat %>%
      scale %>%
      dist %>%
      hclust()
    cl$labels[cl$order]
  }) %>%
  unlist()
  

pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k3_DEGs_rowCluster.pdf'), height = 6, width  = 11)
p <- pheatmap(heatmap_data[row_order, ],
         cluster_rows = F,
         cluster_cols = TRUE,
         annotation_row = annotation,
         # annotation_col = annotation_gene,
         scale = 'column', 
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()

heatmap_data %>% 
  scale %>%
  ComplexHeatmap::Heatmap(
    name = " ",
    row_split = annotation$Cluster,
    left_annotation = ComplexHeatmap::rowAnnotation(df = annotation),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2"
  )
    
annotation_gene <- 
  p$tree_col %>%
  cutree(k = 2) %>%
  enframe() %>% 
  mutate(value = as.character(value)) %>% 
  column_to_rownames("name") 

pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k3_DEGs_rowCluster_colCluster.pdf'), height = 6, width  = 11)
p <- pheatmap(heatmap_data[row_order, ],
              cluster_rows = F,
              cluster_cols = TRUE,
              annotation_row = annotation,
              annotation_col = annotation_gene,
              scale = 'column', 
              clustering_method = "ward.D2",
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              main = "Consensus Clustering Heatmap",
              fontsize_row = 8,
              show_colnames = F,
              show_rownames = T)
dev.off()

# id trans 
annotation_gene_id <- 
  annotation_gene %>%
  with(split(rownames(.), value)) %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Hs.eg.db"
    ) %>%
    .$ENTREZID
  })

# enrich 
enrich_cluster_gene_ls <-
  annotation_gene %>%
  with(split(rownames(.), value))%>%
  map(function(x){
  clusterProfiler::enrichGO(
    x,
    OrgDb = 'org.Hs.eg.db', 
    ont='BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05,
    keyType = 'SYMBOL'
  )
})

enrich_cluster_gene_ls %>%
  imap(~ enrichplot::dotplot(.x, title = .y)) %>%
  patchwork::wrap_plots()
ggsave(file.path(out_dir, 'DEGs_cluster2_enrichment.pdf'), width = 12, height = 6)


# all_data cell prop -------
all_data[[]] %>% 
  dplyr::select(sample, annotation) %>% 
  left_join(enframe(cluster_order, "sample", "cluster") ,by = "sample") %>%
  group_by(cluster, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = cluster, y = freq, fill = annotation)) +
  geom_col() +
  ggsci::scale_fill_igv()
ggsave(file.path(out_dir, 'cell_prop_consensus_clustering.pdf'), width = 6, height = 6)


# CN prop ---------

# enrichment ---------------------------------------
Epithelial_er_output <- file.path(out_dir, 'enrichment') 
dir.create(Epithelial_er_output)

enrich_ls <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Epithelial_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = TRUE
    )
  })

# saveRDS(Epithelial_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })





















