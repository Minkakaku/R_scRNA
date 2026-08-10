
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

T_cell <- readRDS("T_cell_subtype/T_cell_clean.rds")


out_dir <- "result_0622"
dir.create(out_dir)
T_cell_out_dir <- "result_0622/T_cell_subtype_2" 
dir.create(T_cell_out_dir)

keep <- c(
  "Naive_T",
  "TCM",
  #"Proliferative_T",
  "CD4_Treg",
  #"Gamma-delta_T",
  "CD8_TEM",
  "Tfh"
)

T_cell_clean <- subset(T_cell, subtype %in% keep)

T_cell_clean$subtype <- droplevels(T_cell_clean$subtype)
Idents(T_cell_clean) <- T_cell_clean$subtype

table(T_cell_clean$subtype)
saveRDS(T_cell_clean, "result_0622/T_cell_subtype_2/T_cell_clean.rds")

T_cell <- T_cell_clean

pdf(file.path(T_cell_out_dir, str_glue("Dimplot T annotation.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot T annotation split.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## PCA --------第一次
T_cell <- FindVariableFeatures(T_cell)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 15 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)


T_cell <- T_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

T_cell <- RunUMAP(
  T_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
T_cell  <-
  FindNeighbors(
    T_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


T_cell <-
  FindClusters(
    T_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(T_cell)<-'RNA_snn_res.0.6'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(T_cell_markers, file = "result_0622/T_cell_subtype_2/T_cell_all_markers_1.csv", row.names = TRUE)

T_cell_top5_markers <- 
  T_cell_markers %>% 
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

names(T_cell_top5_markers) <- sub('.*\\.','',names(T_cell_top5_markers))

#Idents(T_cell) <- "RNA_snn_res.0.4"
pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(T_cell ,file = "result_0622/T_cell_subtype_2/T_cell_cluster_1.rds")
T_cell <- readRDS("result_0622/T_cell_subtype_2/T_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Plasma-like_T", # 
  '1'  = "TCM", #
  '2'  = "Naive_T", # 
  '3'  = "Plasma-like_T", # 
  '4'  = "cTEC-like_T", # 
  '5'  = "CD4_Treg",#  
  '6'  = "CD8_TEM", #  
  '7'  = "B cell", #  
  '8'  = "Proliferative_T",# 
  '9'  = "Tfh", #  
  '10' = "Plasma-like_T", # 
  '11' = "Plasma-like_T",#  
  '12' = "TCM"
)

Idents(T_cell) <- "RNA_snn_res.0.6"




T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

keep <- c(
  "Naive_T",
  "TCM",
  #"Proliferative_T",
  "CD4_Treg",
  #"Gamma-delta_T",
  "CD8_TEM",
  "Tfh"
)

T_cell_clean <- subset(T_cell, subtype %in% keep)

T_cell_clean$subtype <- droplevels(T_cell_clean$subtype)
Idents(T_cell_clean) <- T_cell_clean$subtype

table(T_cell_clean$subtype)
saveRDS(T_cell_clean, "result_0622/T_cell_subtype_2/T_cell_clean_2.rds")

T_cell <- T_cell_clean

#T_cell <- RunUMAP(T_cell, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation_1.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation split_1.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




## PCA --------第二次
T_cell <- FindVariableFeatures(T_cell)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 17 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)


T_cell <- T_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

T_cell <- RunUMAP(
  T_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
T_cell  <-
  FindNeighbors(
    T_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


T_cell <-
  FindClusters(
    T_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} cluster_2.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(T_cell)<-'RNA_snn_res.0.8'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(T_cell_markers, file = "result_0622/T_cell_subtype_2/T_cell_all_markers_2.csv", row.names = TRUE)

T_cell_top5_markers <- 
  T_cell_markers %>% 
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

names(T_cell_top5_markers) <- sub('.*\\.','',names(T_cell_top5_markers))

#Idents(T_cell) <- "RNA_snn_res.0.4"
pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker_2.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(T_cell ,file = "result_0622/T_cell_subtype_2/T_cell_cluster_2.rds")
T_cell <- readRDS("result_0622/T_cell_subtype_2/T_cell_cluster_2.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "TCM", # 
  '1'  = "Naive_T", #
  '2'  = "Naive_T", # 
  '3'  = "CD4_Treg", # 
  '4'  = "CD8_TEM", # 
  '5'  = "Tfh",#  
  '6'  = "Macrophage-like_T", #  
  '7'  = "T_other", #  
  '8'  = "Gamma-delta_T"
)

Idents(T_cell) <- "RNA_snn_res.0.8"




T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation_2.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split_2.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

keep <- c(
  "Naive_T",
  "TCM",
  #"Proliferative_T",
  "CD4_Treg",
  #"Gamma-delta_T",
  "CD8_TEM",
  "Tfh"
)

T_cell_clean <- subset(T_cell, subtype %in% keep)

T_cell_clean$subtype <- droplevels(T_cell_clean$subtype)
Idents(T_cell_clean) <- T_cell_clean$subtype

table(T_cell_clean$subtype)
saveRDS(T_cell_clean, "result_0622/T_cell_subtype_2/T_cell_clean_3.rds")

T_cell <- T_cell_clean
table(T_cell$subtype)
#T_cell <- RunUMAP(T_cell, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation_2.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation split_2.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




## PCA -------- 第三次
T_cell <- FindVariableFeatures(T_cell)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot_3.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 14 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)


T_cell <- T_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

T_cell <- RunUMAP(
  T_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
T_cell  <-
  FindNeighbors(
    T_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


T_cell <-
  FindClusters(
    T_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} cluster_3.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(T_cell)<-'RNA_snn_res.0.8'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(T_cell_markers, file = "result_0622/T_cell_subtype_2/T_cell_all_markers_3.csv", row.names = TRUE)

T_cell_top5_markers <- 
  T_cell_markers %>% 
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

names(T_cell_top5_markers) <- sub('.*\\.','',names(T_cell_top5_markers))

#Idents(T_cell) <- "RNA_snn_res.0.4"
pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker_3.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(T_cell ,file = "result_0622/T_cell_subtype_2/T_cell_cluster_3.rds")
T_cell <- readRDS("result_0622/T_cell_subtype_2/T_cell_cluster_3.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Naive_T", # 
  '1'  = "CD4_TCM", #
  '2'  = "Naive_T", # 
  '3'  = "CD4_Treg", # 
  '4'  = "CD8_TEM", # 
  '5'  = "Tfh",#  
  '6'  = "CD8_TCM"
)

Idents(T_cell) <- "RNA_snn_res.0.8"




T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation_3.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split_3.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


saveRDS(T_cell, "result_0622/T_cell_subtype_2/T_cell_annotation.rds")




#逐层映射回原注释
T_cell_1 <- readRDS("result_0622/T_cell_subtype_2/T_cell_cluster_2.rds")
cell_cluster <- c(
  '0'  = "TCM", # 
  '1'  = "Naive_T", #
  '2'  = "Naive_T", # 
  '3'  = "CD4_Treg", # 
  '4'  = "CD8_TEM", # 
  '5'  = "Tfh",#  
  '6'  = "Macrophage-like_T", #  
  '7'  = "T_other", #  
  '8'  = "Gamma-delta_T"
)

Idents(T_cell_1) <- "RNA_snn_res.0.8"
T_cell_1 <- RenameIdents(T_cell_1, cell_cluster)
T_cell_1$subtype <- Idents(T_cell_1)

T_cell_1$subtype2 <- as.character(T_cell_1$subtype)
cell_subtype <- T_cell$subtype
names(cell_subtype) <- colnames(T_cell)
T_cell_1$subtype2[colnames(T_cell)] <- as.character(cell_subtype)
table(T_cell_1$subtype2)
table(T_cell_1$subtype)
table(T_cell$subtype)



T_cell_2 <- readRDS("result_0622/T_cell_subtype_2/T_cell_cluster_1.rds")
cell_cluster <- c(
  '0'  = "Plasma-like_T", # 
  '1'  = "TCM", #
  '2'  = "Naive_T", # 
  '3'  = "Plasma-like_T", # 
  '4'  = "cTEC-like_T", # 
  '5'  = "CD4_Treg",#  
  '6'  = "CD8_TEM", #  
  '7'  = "B cell", #  
  '8'  = "Proliferative_T",# 
  '9'  = "Tfh", #  
  '10' = "Plasma-like_T", # 
  '11' = "Plasma-like_T",#  
  '12' = "TCM"
)

Idents(T_cell_2) <- "RNA_snn_res.0.6"
T_cell_2 <- RenameIdents(T_cell_2, cell_cluster)
T_cell_2$subtype <- Idents(T_cell_2)

T_cell_2$subtype2 <- as.character(T_cell_2$subtype)
cell_subtype <- T_cell_1$subtype2
names(cell_subtype) <- colnames(T_cell_1)
T_cell_2$subtype2[colnames(T_cell_1)] <- as.character(cell_subtype)
table(T_cell_2$subtype2)
table(T_cell_2$subtype)
table(T_cell_1$subtype2)


T_cell_3 <- readRDS("T_cell_subtype/T_cell_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Naive_T", # 
  '1'  = "Naive_T", #
  '2'  = "TCM", #不知道怎么注释，直接根据top基因看的
  '3'  = "Proliferative_T", # 
  '4'  = "CD4_Treg", # 
  '5'  = "Proliferative_T", # 
  '6'  = "Proliferative_T", # 
  '7'  = "mTEC", # CDH1和EGFR像上皮细胞的基因
  '8'  = "Gamma-delta_T",
  '9'  = "CD8_TEM", # 高表达GZMK和GZMA，低表达LEF1和TCF7，是TEM
  '10' = "LAMP3_DC", # 
  '11' = "Naive_T", # 根据位置判断
  '12' = "Tfh", # 
  '13' = "B cell"
)

Idents(T_cell_3) <- "RNA_snn_res.0.6"
T_cell_3 <- RenameIdents(T_cell_3, cell_cluster)
T_cell_3$subtype <- Idents(T_cell_3)

T_cell_3$subtype2 <- as.character(T_cell_3$subtype)
cell_subtype <- T_cell_2$subtype2
names(cell_subtype) <- colnames(T_cell_2)
T_cell_3$subtype2[colnames(T_cell_2)] <- as.character(cell_subtype)
table(T_cell_3$subtype2)
table(T_cell_3$subtype)
table(T_cell_2$subtype2)


######
all_data <- readRDS("output/data/all_data.rds")
all_data$annotation %>% table()

#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- T_cell_3$subtype2
names(cell_subtype) <- colnames(T_cell_3)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(T_cell_3)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(T_cell_3$subtype2)
saveRDS(all_data, "result_0622/T_cell_subtype_2/all_data_T_annotation.rds")
saveRDS(T_cell_3 ,file = "result_0622/T_cell_subtype_2/T_cell_allannotation.rds")

#到这里映射结束，得到alldata和T cell全部注释(包括双胞)




T_cell <- readRDS("result_0622/T_cell_subtype_2/T_cell_annotation.rds")
#画一个注释完的点图
T_major_marker_ls_2 <- c(
  "CD8"="CD8A",
  "CD8"="CD8B",
  
  "CD4"="CD4",
  
  # "Naive"="CCR6",
  "Naive"="LEF1",
  "Naive"="TCF7",
  #"TCM"="SELL",
  
  "TCM"="SELL",
  "TCM"="CD27", 
  "TCM"="CD44",
  #"Effector"="GZMB",
  #"Effector"="CX3CR1", #https://pubmed.ncbi.nlm.nih.gov/35361683/
  
  "TEM"="GZMK", 
  "TEM"="GZMA", 
  
  # regulatory T cells (Treg)
  "Treg" = "FOXP3",
  "Treg" = "IL2RA",
  
  "Tfh"="CD200",
  "Tfh"="ICOS",
  "Tfh"="PDCD1"
)

pdf(file.path(T_cell_out_dir, str_glue("dotplot T subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(T_cell, group.by = "subtype", features = T_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


T_cell$subtype %>% table()
T_cell$sample %>% table()
T_cell$group %>% table()
## cell props barplot  ----------------------------------------------------------
T_cell_meta <- T_cell[[]]

T_cell_meta$sample <- factor(T_cell_meta$sample,levels=c(T_cell_meta$sample %>% unique()))

p_col_proj <-
  T_cell_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  #ggsci::scale_fill_igv() +
  scale_fill_brewer(palette = "Set2", name = "T_cell Subtype") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(T_cell_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  T_cell_meta %>%
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
  filename = file.path(T_cell_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)




# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    T_cell[['sample']],
    Embeddings(T_cell[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    T_cell, 
    vars = c(subtype_col)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  T_cell$subtype %>% unique %>% as.character()
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)



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
  file.path("result_0622/T_cell_subtype_2/spatial plot of celltype in each sample.pdf"),
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
      file.path(str_glue("result_0622/T_cell_subtype_2/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })