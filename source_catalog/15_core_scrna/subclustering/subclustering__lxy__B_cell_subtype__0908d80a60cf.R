library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('result_0622/DC_subtype/all_data_T_DC_annotation.rds')

all_data$annotation2 %>% table()


B_cell_out_dir <- "result_0622/B_cell_subtype" 
dir.create(B_cell_out_dir)


all_data$annotation2 <- str_replace(all_data$annotation2, 'Plasma', 'B cell')

B_cell <- subset(all_data, annotation2 == 'B cell')

## PCA --------
B_cell <- FindVariableFeatures(B_cell)

B_cell <- ScaleData(B_cell, vars.to.regress = 'nCount_RNA')
B_cell <- RunPCA(B_cell)

### pca inspector ------
pdf(file.path(B_cell_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B_cell, ndims = 50)
dev.off()


## dimension reduction --------
# dimplot
pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(B_cell) <- "RNA_snn_res.0.8"

### dotplot -------
B_cell_major_marker_ls <- c(
  "B cell" = "CD19",
  "B cell" = "CD79A",
  "B cell" = "MS4A1", 
  
  "naive B" = "CD20",
  "naive B" = "CD24",
  "naive B" = "IGHM",
  "naive B" = "BANK1",
  "naive B" = "FCRL1",
  "naive B" = "HLA-DQB1",
  "naive B" = "HLA-DQA1",
  "naive B" = "HVCN1",
  "naive B" = "IGHD",
  "naive B" = "FCER2",
  "naive B" = "IL4R", 
  
  "Follicular B" = "CD40", 
  "Follicular B" = "BCL6", 
  
  "memory B"= "CD27",
  "memory B"= "FCRL4",
  "memory B"= "TACI",
  "memory B"= "AIM2",
  "memory B"= "TXNIP",
  "memory B"= "GPR183",
  "memory B"= "TNFRSF13B",
  "memory B"= "CD22",
  "memory B"= "FCRL2",
  "memory B"= "ADAM28",
  "memory B"= "IGLL5",
  "memory B"= "P2RX5",
  "memory B"= "BLK", 
  
  # "Marginal zone B" = "CD27", 
  #"Marginal zone B" = "CR2",
  #"Marginal zone B" = "CD1C", 
  
  "cycling B"="MKI67",
  "cycling B"="STMN1",
  "cycling B"="HMGB2",
  "cycling B"="BIRC5", 
  
  "Plasma" = "CD138",
  "Plasma" = "PRDM1",
  "Plasma" = "BLIMP1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHA2",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHG2",
  "Plasma" = "IGHG3",
  "Plasma" = "IGHG4",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  "Plasma" = "SDC1", 
  
  #"Plasmablast" = "MKI67", 
  "Plasmablast" = "CD38", 
  "Plasmablast" = "IRF4", 
  
  "Breg" = "IL10", 
  "Breg" = "CD1D", 
  "Breg" = "TGFB1", 
  
  "GC B"= "S1PI2",
  # "GC B"="BCL6",
  "GC B"="RGS13",
  "GC B"="TCL1A",
  "GC B"= "LRMP",
  "GC B"= "SUGCT",
  "GC B"= "MME",
  "GC B"="AICDA",
  
  "FDC" = "CR2",
  "FDC" = "CR1",
  #"FDC" = "FCER2",
  "FDC" = "FDC"
)

B_cell_major_marker_ls_1 <- c(
  # 活化B细胞 Activate B
  "Activate B" = "IL2RA",
  "Activate B" = "TNFRSF8",
  # 未成熟B Immature B
  "Immature B" = "MS4A1",
  # 初始B naive B
  "naive B" = "CD19",
  "naive B" = "TCL1A",
  "naive B" = "IL4R",
  # 记忆B memory B
  "memory B" = "CD27",
  "memory B" = "AIM2",
  "memory B" = "TNFRSF13B",
  "memory B" = "CD80",
  "memory B" = "PDCD1LG2",
  "memory B" = "CXCR3",
  # 滤泡B Follicular B
  "Follicular B" = "CD22",
  "Follicular B" = "IGHD",
  "Follicular B" = "FCER2",
  #"Follicular B" = "CR2",
  # 调节性B B reg
  "B reg" = "CD5",
  "B reg" = "CD24",
  "B reg" = "TLR4",
  "B reg" = "TGFB1",
  # 边缘区B Margine zone B
  #"Margine zone B" = "CD27",
  #"Margine zone B" = "NOTCH2",
  #"Margine zone B" = "CR2",
  # 生发中心B Germinal center B
  "Germinal center B" = "MME",
  "Germinal center B" = "AICDA",
  "Germinal center B" = "SUGCT",
  "Germinal center B" = "MKI67",
  # 浆母细胞 Plasmablast
  "Plasmablast" = "IRF4",
  "Plasmablast" = "PRDM1",
  "Plasmablast" = "CD38",
  "FDC" = "CR2",
  "FDC" = "CR1",
  #"FDC" = "FCER2",
  "FDC" = "FDC"
)

Idents(B_cell)<-'RNA_snn_res.0.8'
B_cell_markers <-
  FindAllMarkers(
    JoinLayers(B_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(B_cell_markers, file = "result_0622/B_cell_subtype/B_cell_all_markers.csv", row.names = TRUE)

B_cell_top5_markers <- 
  B_cell_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(B_cell_top5_markers) <- sub('.*\\.','',names(B_cell_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(B_cell_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(B_cell, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(B_cell, group.by = "RNA_snn_res.0.8", features = B_cell_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(B_cell, group.by = "RNA_snn_res.0.8", features = B_cell_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(B_cell, group.by = "RNA_snn_res.0.8", features = B_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(B_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(B_cell ,file = "result_0622/B_cell_subtype/B_cell_cluster.rds")
B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Naive_B", # 
  '1'  = "T_cell-like_B_cell", #
  '2'  = "Memory_B", # 
  '3'  = "Immature_B", # 
  '4'  = "Plasmablast", # 
  '5'  = "Immature_B",#  
  '6'  = "GC_B", #  
  '7'  = "T_cell-like_B_cell", #  
  '8'  = "Immature_B", #  
  '9'  = "Immature_B", #  
  '10'  = "T_cell-like_B_cell", #  
  '11'  = "Immature_B", # 
  '12'  = "Naive_B", # 
  '13'  = "GC_B", #  
  '14'  = "Plasmablast", # 
  '15'  = "Plasmablast"
)

Idents(B_cell) <- "RNA_snn_res.0.8"




B_cell <- RenameIdents(B_cell, cell_cluster)
B_cell$subtype <- Idents(B_cell)


pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(B_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "Plasmablast",
  "Immature_B",
  "Naive_B",
  "Memory_B",
  "GC_B"
)

B_cell_clean <- subset(B_cell, subtype %in% keep)

B_cell_clean$subtype <- droplevels(B_cell_clean$subtype)
Idents(B_cell_clean) <- B_cell_clean$subtype

table(B_cell_clean$subtype)
saveRDS(B_cell_clean, "result_0622/B_cell_subtype/B_cell_clean.rds")

B_cell <- B_cell_clean

B_cell <- RunUMAP(B_cell, dims = 1:{dims},n.neighbors = 50,return.model = TRUE)

pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} B_cell annotation.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} B_cell annotation split.pdf")), height = 6, width = 15)
DimPlot(B_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()



# 保留：0,1,7
B_cell <- subset(B_cell, RNA_snn_res.0.8 %in% c(0,2,3,4,5,6,8,9,11,12,13,14,15))
B_cell$RNA_snn_res.0.8 %>% table()


## PCA --------
B_cell <- FindVariableFeatures(B_cell)

B_cell <- ScaleData(B_cell, vars.to.regress = 'nCount_RNA')
B_cell <- RunPCA(B_cell)

### pca inspector ------
pdf(file.path(B_cell_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(B_cell, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 16 # 18

B_cell <- RunUMAP(B_cell, dims = 1:dims,return.model = TRUE)


B_cell <- B_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

B_cell <- RunUMAP(
  B_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
B_cell  <-
  FindNeighbors(
    B_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


B_cell <-
  FindClusters(
    B_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(B_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(B_cell)<-'RNA_snn_res.0.6'
B_cell_markers <-
  FindAllMarkers(
    JoinLayers(B_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(B_cell_markers, file = "result_0622/B_cell_subtype/B_cell_all_markers_1.csv", row.names = TRUE)

B_cell_top5_markers <- 
  B_cell_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(B_cell_top5_markers) <- sub('.*\\.','',names(B_cell_top5_markers))


pdf(file.path(B_cell_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(B_cell, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(B_cell, group.by = "RNA_snn_res.0.6", features = B_cell_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(B_cell, group.by = "RNA_snn_res.0.4", features = B_cell_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(B_cell, group.by = "RNA_snn_res.0.6", features = B_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(B_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

saveRDS(B_cell ,file = "result_0622/B_cell_subtype/B_cell_cluster_1.rds")
B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Memory_B", # 
  '1'  = "Naive_B", #
  '2'  = "Plasmablast", # 
  '3'  = "Immature_B", # 
  '4'  = "GC_B", # 
  '5'  = "Immature_B",#  
  '6'  = "Immature_B", #  
  '7'  = "Immature_B", #  
  '8'  = "Immature_B", #  
  '9'  = "FDC", #  
  '10'  = "GC_B", #  
  '11'  = "Plasmablast"
)

Idents(B_cell) <- "RNA_snn_res.0.6"




B_cell <- RenameIdents(B_cell, cell_cluster)
B_cell$subtype <- Idents(B_cell)


pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_cell_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(B_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



B_cell_1 <- readRDS("result_0622/B_cell_subtype/B_cell_cluster.rds")
cell_cluster <- c(
  '0'  = "Naive_B", # 
  '1'  = "T_cell-like_B_cell", #
  '2'  = "Memory_B", # 
  '3'  = "Immature_B", # 
  '4'  = "Plasmablast", # 
  '5'  = "Immature_B",#  
  '6'  = "GC_B", #  
  '7'  = "T_cell-like_B_cell", #  
  '8'  = "Immature_B", #  
  '9'  = "Immature_B", #  
  '10'  = "T_cell-like_B_cell", #  
  '11'  = "Immature_B", # 
  '12'  = "Naive_B", # 
  '13'  = "GC_B", #  
  '14'  = "Plasmablast", # 
  '15'  = "Plasmablast"
)

Idents(B_cell_1) <- "RNA_snn_res.0.8"
B_cell_1 <- RenameIdents(B_cell_1, cell_cluster)
B_cell_1$subtype <- Idents(B_cell_1)

B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Memory_B", # 
  '1'  = "Naive_B", #
  '2'  = "Plasmablast", # 
  '3'  = "Immature_B", # 
  '4'  = "GC_B", # 
  '5'  = "Immature_B",#  
  '6'  = "Immature_B", #  
  '7'  = "Immature_B", #  
  '8'  = "Immature_B", #  
  '9'  = "FDC", #  
  '10'  = "GC_B", #  
  '11'  = "Plasmablast"
)
Idents(B_cell) <- "RNA_snn_res.0.6"
B_cell <- RenameIdents(B_cell, cell_cluster)
B_cell$subtype <- Idents(B_cell)


B_cell_1$subtype2 <- as.character(B_cell_1$subtype)
cell_subtype <- B_cell$subtype
names(cell_subtype) <- colnames(B_cell)
B_cell_1$subtype2[colnames(B_cell)] <- as.character(cell_subtype)
table(B_cell_1$subtype2)
table(B_cell_1$subtype)
table(B_cell$subtype)
saveRDS(B_cell_1, "result_0622/B_cell_subtype/B_cell_allannotation.rds")


######下面合并到alldata的注释先不运行
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- B_cell_1$subtype2
names(cell_subtype) <- colnames(B_cell_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(B_cell_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(B_cell_1$subtype2)
saveRDS(all_data, "result_0622/B_cell_subtype/all_data_T_annotation.rds")
saveRDS(B_cell_1 ,file = "result_0622/B_cell_subtype/B_cell_annotation.rds")



B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_clean.rds")
#画一个注释完的点图
B_cell_major_marker_ls_2 <- c(
  # 未成熟B Immature B
  "Immature B" = "CD19",
  "Immature B" = "MS4A1",
  # 初始B naive B
  "naive B" = "BANK1",
  "naive B" = "FCRL1",
  "naive B" = "FCER2",
  "naive B" = "IL4R", 
  # 记忆B memory B
  "memory B"= "AIM2",
  "memory B"= "TNFRSF13B",
  "memory B"= "CD22",
  "memory B"= "FCRL2",
  "memory B"= "ADAM28",
  "memory B"= "BLK", 
  
  "GC B"="TCL1A",
  "GC B"= "MME",
  "GC B"="AICDA",
  "GC B"="MKI67",
  
  "Plasmablast" = "IRF4",
  "Plasmablast" = "PRDM1",
  "Plasmablast" = "CD38",
  
  "FDC" = "CR1",
  "FDC" = "CR2"
)

pdf(file.path(B_cell_out_dir, str_glue("dotplot B_cell subtype major_marker_1.pdf")), width = 12, height = 6)
DotPlot(B_cell, group.by = "subtype", features = B_cell_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


B_cell$subtype %>% table()
B_cell$sample %>% table()
B_cell$group %>% table()
## cell props barplot  ----------------------------------------------------------
B_cell_meta <- B_cell[[]]

B_cell_meta$sample <- factor(B_cell_meta$sample,levels=c(B_cell_meta$sample %>% unique()))

p_col_proj <-
  B_cell_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "B_cell Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(B_cell_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  B_cell_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(B_cell_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    B_cell[['sample']],
    Embeddings(B_cell[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    B_cell, 
    vars = c(subtype_col)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  B_cell$subtype %>% unique %>% as.character() 
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)

cell_colors <- c(
  "Naive_B" = "#F0E685FF",      # 红色系
  "Memory_B" = "#CE3D32FF",     # 青色
  "Immature_B" = "#5050FFFF",   # 蓝色
  "Plasmablast" = "#749B58FF",  # 浅绿色
  "GC_B" = "#466983FF",          # 黄色
  "FDC" = "#BA6338FF"
)

# 确保颜色顺序和因子水平一致
cell_colors <- cell_colors[levels(B_cell$subtype)]


p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[subtype_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path("result_0622/B_cell_subtype/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()


# plot each cell in all sample -----------
p_sp_ct_ls_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    
    p_ls <-
      all_cells %>%
      setNames(., .) %>%
      map(function(x){
        cell_colors[names(cell_colors) != x] <- NA
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[subtype_col]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_colors) +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(sam)
        p
      })
    p_ls
  })

# export by sample
p_sp_ct_ls_ls %>%
  iwalk(function(p_ls, sample){
    pdf(
      file.path(str_glue("result_0622/B_cell_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })
