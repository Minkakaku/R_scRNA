library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/output/data/all_data.rds')

all_data$annotation %>% table()


Proximal_tubule_out_dir <- "Proximal_tubule_subtype" 
dir.create(Proximal_tubule_out_dir)
Proximal_tubule <- subset(all_data, annotation == 'Proximal tubule')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Proximal_tubule <- FindVariableFeatures(Proximal_tubule)

Proximal_tubule <- ScaleData(Proximal_tubule, vars.to.regress = 'nCount_RNA')
Proximal_tubule <- RunPCA(Proximal_tubule)

### pca inspector ------
pdf(file.path(Proximal_tubule_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Proximal_tubule, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 16 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

Proximal_tubule <- RunUMAP(Proximal_tubule, dims = 1:dims,return.model = TRUE)

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Proximal_tubule <- Proximal_tubule %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Proximal_tubule <- RunUMAP(
  Proximal_tubule, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Proximal_tubule  <-
  FindNeighbors(
    Proximal_tubule ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Proximal_tubule <-
  FindClusters(
    Proximal_tubule,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Proximal_tubule) <- "RNA_snn_res.0.6"

### dotplot -------
PT_major_marker_ls <- c(
  # 全段PT
  "PT" = "LRP2",
  "PT" = "CUBN",
  "PT" = "SLC13A1",
  "PT" = "AQP1",
  "PT" = "SLC27A2",
  
  # 近端小管损伤/损伤活化marker（全PT损伤通用）
  "PT_injury" = "KIM1",
  "PT_injury" = "HAVCR1",
  "PT_injury" = "SLC3A1",
  "PT_injury" = "LCN2",
  "PT_injury" = "NGAL",
  "PT_injury" = "SERPINA3",
  "PT_injury" = "COL1A1",
  "PT_injury" = "FN1",
  "PT_injury" = "CXCL8",
  "PT_injury" = "CCL2",
  
  # 近端小管S1段 (PT-S1)
  "PT-S1" = "PRODH2",
  "PT-S1" = "SLC5A2",
  #"PT-S1" = "SLC22A8",
  "PT-S1" = "SLC16A1",
  "PT-S1" = "SLC34A3",
  "PT-S1" = "SPP2",
  
  #S1段/S2段共同marker
  "PT-S1-2" = "SLC5A12",
  "PT-S1-2" = "SLC13A3",
  "PT-S1-2" = "SLC22A6",
  #"PT-S1-2" = "SLC22A8",
  
  # 近端小管S2段 (PT-S2)
  "PT-S2" = "SLC34A1",
  "PT-S2" = "SLC22A8",
  "PT-S2" = "SLC22A12",
  "PT-S2" = "ACAA1B",
  "PT-S2" = "CNDP2",
  "PT-S2" = "NUDT19",
  
  
  #S2段/S3段共同marker
  "PT-S2-3" = "SLC22A7",
  
  # 近端小管S3段 (PT-S3)
  "PT-S3" = "MOGAT1",
  "PT-S3" = "SLC5A11",
  "PT-S3" = "SLC22A24",
  "PT-S3" = "SLC7A13",
  "PT-S3" = "SLC5A8",
  "PT-S3" = "ABCC3",
  "PT-S3" = "SATB2",
  "PT-S3" = "ATP11A",
  "PT-S3" = "MALAT1",
  "PT-S3" = "TRIM7",
  "PT-S3" = "SLC22A30",
  "PT-S3" = "SLC6A18"
)



PT_major_marker_ls_1 <- c(
  # 近端小管S1段 (PT-S1)
  #"PT-S1" = "PRODH2",
  "PT-S1" = "FTCD",
  #"PT-S1" = "SLC22A6",
  "PT-S1" = "RNF212B",
  "PT-S1" = "SLC7A8",
  "PT-S1" = "PCDH15",
  #"PT-S1" = "PAH",
  "PT-S1" = "SORCS1",
  "PT-S1" = "AGXT2",
  
  # 近端小管S2段 (PT-S2)
  "PT-S2" = "SLC5A10",
  "PT-S2" = "SLC34A1",
  #"PT-S2" = "SLC13A3",
  "PT-S2" = "SLC5A12",
  "PT-S2" = "MIOX",
  
  # 近端小管S3段 (PT-S3)
  "PT-S3" = "SLC5A11",
  "PT-S3" = "SLC22A7",
  "PT-S3" = "LINC00671",
  "PT-S3" = "C5orf67",
  "PT-S3" = "GPT",
  "PT-S3" = "CYP4A11",
  "PT-S3" = "SLC30A8",
  "PT-S3" = "CDH9",
  "PT-S3" = "SLC6A13",
  "PT-S3" = "SLC28A1",
  "PT-S3" = "EBNA1BP2",
  "PT-S3" = "ANKS1B",
  "PT-S3" = "CLDN10",
  "PT-S3" = "SLC22A12",
  "PT-S3" = "GLYATL1",
  "PT-S3" = "ACE2",
  "PT-S3" = "MT1H",
  "PT-S3" = "SLC13A1",
  "PT-S3" = "CUBN",
  "PT-S3" = "SLC16A9",
  "PT-S3" = "AK4",
  
  # S1/S2共用近端小管 (PT-S1-2)
  "PT-S1-2" = "ALB",
  "PT-S1-2" = "FAM151A",
  "PT-S1-2" = "SLC22A8",
  "PT-S1-2" = "GSTA2",
  "PT-S1-2" = "GSTA1",
  "PT-S1-2" = "UGT2B7",
  "PT-S1-2" = "GPX3",
  "PT-S1-2" = "PAH",
  "PT-S1-2" = "SLC22A6",
  "PT-S1-2" = "SLC13A3",
  "PT-S1-2" = "PRODH2"
)

Idents(Proximal_tubule)<-'RNA_snn_res.0.4'
Proximal_tubule_markers <-
  FindAllMarkers(
    JoinLayers(Proximal_tubule),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Proximal_tubule_markers, file = "Proximal_tubule_subtype/Proximal_tubule_all_markers.csv", row.names = TRUE)

Proximal_tubule_top5_markers <- 
  Proximal_tubule_markers %>% 
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

names(Proximal_tubule_top5_markers) <- sub('.*\\.','',names(Proximal_tubule_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Proximal_tubule_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = Proximal_tubule_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Proximal_tubule ,file = "Proximal_tubule_subtype/Proximal_tubule_cluster.rds")
Proximal_tubule <- readRDS("Proximal_tubule_subtype/Proximal_tubule_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "PT_Others", # 
  '1'  = "PT_Others", #
  '2'  = "PTS3", # 
  '3'  = "PTS1", # CD8+TEM
  '4'  = "PTS2", # 
  '5'  = "PTS2",#  
  '6'  = "Endothelial-like_PT", #  
  '7'  = "Mesangial cell-like_PT", #  
  '8'  = "LOH-like_PT"
)

Idents(Proximal_tubule) <- "RNA_snn_res.0.4"




Proximal_tubule <- RenameIdents(Proximal_tubule, cell_cluster)
Proximal_tubule$subtype <- Idents(Proximal_tubule)


pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Proximal_tubule, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "PTS1",
  "PTS2",
  "PTS3"
)

Proximal_tubule_clean <- subset(Proximal_tubule, subtype %in% keep)

Proximal_tubule_clean$subtype <- droplevels(Proximal_tubule_clean$subtype)
Idents(Proximal_tubule_clean) <- Proximal_tubule_clean$subtype

table(Proximal_tubule_clean$subtype)
saveRDS(Proximal_tubule_clean, "Proximal_tubule_subtype/Proximal_tubule_clean.rds")

Proximal_tubule <- Proximal_tubule_clean

Proximal_tubule <- RunUMAP(Proximal_tubule, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} Proximal_tubule annotation.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} Proximal_tubule annotation split.pdf")), height = 6, width = 15)
DimPlot(Proximal_tubule, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




# 保留：0,1,7
Proximal_tubule <- subset(Proximal_tubule, RNA_snn_res.0.4 %in% c(2,3,4,5))
Proximal_tubule$RNA_snn_res.0.4 %>% table()


## PCA --------
Proximal_tubule <- FindVariableFeatures(Proximal_tubule)

Proximal_tubule <- ScaleData(Proximal_tubule, vars.to.regress = 'nCount_RNA')
Proximal_tubule <- RunPCA(Proximal_tubule)

### pca inspector ------
pdf(file.path(Proximal_tubule_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(Proximal_tubule, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 14 # 18

Proximal_tubule <- RunUMAP(Proximal_tubule, dims = 1:dims,return.model = TRUE)


Proximal_tubule <- Proximal_tubule %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Proximal_tubule <- RunUMAP(
  Proximal_tubule, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Proximal_tubule  <-
  FindNeighbors(
    Proximal_tubule ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Proximal_tubule <-
  FindClusters(
    Proximal_tubule,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Proximal_tubule, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Proximal_tubule, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Proximal_tubule)<-'RNA_snn_res.0.4'
Proximal_tubule_markers <-
  FindAllMarkers(
    JoinLayers(Proximal_tubule),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Proximal_tubule_markers, file = "Proximal_tubule_subtype/Proximal_tubule_all_markers_1.csv", row.names = TRUE)

Proximal_tubule_top5_markers <- 
  Proximal_tubule_markers %>% 
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

names(Proximal_tubule_top5_markers) <- sub('.*\\.','',names(Proximal_tubule_top5_markers))

#Idents(Proximal_tubule) <- "RNA_snn_res.0.4"
pdf(file.path(Proximal_tubule_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", features = PT_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.2", features = Proximal_tubule_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(Proximal_tubule_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.4", features = Proximal_tubule_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Proximal_tubule, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

saveRDS(Proximal_tubule ,file = "Proximal_tubule_subtype/Proximal_tubule_cluster_1.rds")
Proximal_tubule <- readRDS("Proximal_tubule_subtype/Proximal_tubule_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "PTS2", # 
  '1'  = "PTS2", #
  '2'  = "PTS1", # 
  '3'  = "PTS3", # 
  '4'  = "PTS2", # 
  '5'  = "Macrophage-like_PT"
)

Idents(Proximal_tubule) <- "RNA_snn_res.0.4"




Proximal_tubule <- RenameIdents(Proximal_tubule, cell_cluster)
Proximal_tubule$subtype <- Idents(Proximal_tubule)


pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Proximal_tubule, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



keep <- c(
  "PTS1",
  "PTS2",
  "PTS3"
)

Proximal_tubule_clean <- subset(Proximal_tubule, subtype %in% keep)

Proximal_tubule_clean$subtype <- droplevels(Proximal_tubule_clean$subtype)
Idents(Proximal_tubule_clean) <- Proximal_tubule_clean$subtype

table(Proximal_tubule_clean$subtype)
saveRDS(Proximal_tubule_clean, "Proximal_tubule_subtype/Proximal_tubule_clean.rds")

Proximal_tubule <- Proximal_tubule_clean



#Proximal_tubule <- RunUMAP(Proximal_tubule, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} Proximal_tubule annotation_1.pdf")), height = 6, width = 8)
DimPlot(Proximal_tubule, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Proximal_tubule_out_dir, str_glue("Dimplot dim{dims} Proximal_tubule annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Proximal_tubule, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()




Proximal_tubule_1 <- readRDS("Proximal_tubule_subtype/Proximal_tubule_cluster.rds")
cell_cluster <- c(
  '0'  = "PT_Others", # 
  '1'  = "PT_Others", #
  '2'  = "PTS3", # 
  '3'  = "PTS1", # 
  '4'  = "PTS2", # 
  '5'  = "PTS2",#  
  '6'  = "Endothelial-like_PT", #  
  '7'  = "Mesangial cell-like_PT", #  
  '8'  = "LOH-like_PT"
)

Idents(Proximal_tubule_1) <- "RNA_snn_res.0.4"
Proximal_tubule_1 <- RenameIdents(Proximal_tubule_1, cell_cluster)
Proximal_tubule_1$subtype <- Idents(Proximal_tubule_1)

Proximal_tubule <- readRDS("Proximal_tubule_subtype/Proximal_tubule_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "PTS2", # 
  '1'  = "PTS2", #
  '2'  = "PTS1", # 
  '3'  = "PTS3", # 
  '4'  = "PTS2", # 
  '5'  = "Macrophage-like_PT"
)

Idents(Proximal_tubule) <- "RNA_snn_res.0.4"
Proximal_tubule <- RenameIdents(Proximal_tubule, cell_cluster)
Proximal_tubule$subtype <- Idents(Proximal_tubule)


Proximal_tubule_1$subtype2 <- as.character(Proximal_tubule_1$subtype)
cell_subtype <- Proximal_tubule$subtype
names(cell_subtype) <- colnames(Proximal_tubule)
Proximal_tubule_1$subtype2[colnames(Proximal_tubule)] <- as.character(cell_subtype)
table(Proximal_tubule_1$subtype2)
table(Proximal_tubule_1$subtype)
table(Proximal_tubule$subtype)
saveRDS(Proximal_tubule_1, "Proximal_tubule_subtype/Proximal_tubule_allannotation.rds")
######下面合并到alldata的注释先不运行
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Proximal_tubule_1$subtype2
names(cell_subtype) <- colnames(Proximal_tubule_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Proximal_tubule_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Proximal_tubule_1$subtype2)
saveRDS(all_data, "Proximal_tubule_subtype/all_data_T_annotation.rds")
saveRDS(Proximal_tubule_1 ,file = "Proximal_tubule_subtype/Proximal_tubule_annotation.rds")



Proximal_tubule <- readRDS("Proximal_tubule_subtype/Proximal_tubule_clean.rds")
#画一个注释完的点图
PT_major_marker_ls_2 <- c(
  "PT" = "LRP2",
  "PT" = "CUBN",
  "PTS1" = "SLC5A2",
  "PTS1" = "SLC16A1",
  "PTS2" = "SLC22A12",
  "PTS2" = "SLC22A8",
  "PTS3" = "ABCC3",
  "PTS3" = "SATB2",
  "PTS3" = "GPT"
)

pdf(file.path(Proximal_tubule_out_dir, str_glue("dotplot Proximal_tubule subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Proximal_tubule, group.by = "subtype", features = PT_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


Proximal_tubule$subtype %>% table()
Proximal_tubule$sample %>% table()
Proximal_tubule$group %>% table()
## cell props barplot  ----------------------------------------------------------
Proximal_tubule_meta <- Proximal_tubule[[]]

Proximal_tubule_meta$sample <- factor(Proximal_tubule_meta$sample,levels=c(Proximal_tubule_meta$sample %>% unique()))

p_col_proj <-
  Proximal_tubule_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Proximal_tubule Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Proximal_tubule_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Proximal_tubule_meta %>%
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
  filename = file.path(Proximal_tubule_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#改图
# 提取meta
Proximal_tubule_meta <- Proximal_tubule[[]]

group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("Normal", "Mild", "Severe")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
Proximal_tubule_meta$Group <- recode(Proximal_tubule_meta$group, !!!group_map)
Proximal_tubule_meta$Group <- factor(Proximal_tubule_meta$Group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Proximal_tubule_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("Proximal_tubule_subtype", "Proximal_tubule_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(Group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("Proximal_tubule_subtype", "Proximal_tubule_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = Group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Proximal_tubule Subtype") +
  labs(title = "Proximal_tubule Subtype Composition by Group",
       x = "Group", y = "Percentage (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(file.path("Proximal_tubule_subtype", "Proximal_tubule_StackedBar.pdf"), 
       p_category_stack, width = 8, height = 6, dpi = 300)
print(p_category_stack)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("Normal","Mild"), c("Normal","Severe"), c("Mild","Severe"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = Group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Proximal_tubule - Cell Type Percentages",
       x = "Group", y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  stat_compare_means(
    comparisons = comp_list,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  )

# 保存+展示
ggsave(file.path("Proximal_tubule_subtype","Proximal_tubule_Boxplot.pdf"), p, width = 12, height = 5, dpi = 300)
print(p)


# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    Proximal_tubule[['sample']],
    Embeddings(Proximal_tubule[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    Proximal_tubule, 
    vars = c(subtype_col)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  Proximal_tubule$subtype %>% unique %>% as.character()
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
  file.path("Proximal_tubule_subtype/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()

