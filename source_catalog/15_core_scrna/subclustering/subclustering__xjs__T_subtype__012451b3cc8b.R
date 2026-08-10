library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
  'epithelials'  = 'EPCAM',
  'epithelials'  = 'KRT5',
  'epithelials'  = 'KRT14',
  
  "Foveolar cells" = "GKN1", #ACT
  "Foveolar cells" = "MUC5AC",
  "Foveolar cells" = "TFF1",
  
  'endothelials' = 'PECAM1', 
  'endothelials' = 'VWF', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  'Macrophage'  = 'CSF1R',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4", 
  'pDCs' = "IRF4", 
  'pDCs' = "IRF7", 
  'pDCs' = "IRF8", 
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  'Lamp3_DCs' = "LAMP3",
  
  'netrophil' = "CSF3R",
  'netrophil' = "G0S2",
  'netrophil' = "RETNLG",
  
  "Myocytes" = "TNNC1", 
  "Myocytes" = "MYH7", 
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "KIT"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )




all_data$annotation %>% table()
#T cells------------------------------------------------------------------------------------
All_T <- subset(all_data, annotation == 'T cell')


## PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')
All_T <- RunPCA(All_T)

T_out_dir<-'T_subtype' #dir.create(T_out_dir)


### pca inspector ------
pdf(file.path(T_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

pdf(file.path(T_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 1:5)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 6:10)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 11:15)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T <- RunUMAP(All_T, dims = 1:dims,return.model = TRUE)



pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


All_T <- All_T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T <- RunUMAP(
  All_T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T  <-
  FindNeighbors(
    All_T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T <-
  FindClusters(
    All_T,
    resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
  )



# dimplot
pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
T_major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  #'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

CD4_major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'T cells' = 'CD3D', 
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'T cells' = 'TRAC',
  'CD4'='CD4',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2'
)


CD8_major_marker_ls <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAYN',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)

pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} multi_res.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "RNA_snn_res.0.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.1", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.2", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.1", features = T_top5_markers)+RotatedAxis()
dev.off()

# find markers
Idents(All_T) <- "RNA_snn_res.0.6"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
T_markers %>% 
  openxlsx::write.xlsx(file.path(T_out_dir, 'diff_gene.xlsx'))

T_top5_markers <- 
  T_markers %>% 
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

names(T_top5_markers) <- sub('.*\\.','',names(T_top5_markers))



# dotplot
pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = T_top5_markers)+RotatedAxis()
dev.off()

DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()

# Idents(All_T) <- "RNA_snn_res.0.8"
# All_T <- FindSubCluster(All_T, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.6)
# Idents(All_T) <- "sub.cluster"
# # find markers
# T_markers_sub <-
#   FindAllMarkers(
#     JoinLayers(All_T),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# 
# T_top5_markers_sub <- 
#   T_markers_sub %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 30) %>%
#   arrange(desc(pct.1)) %>%
#   slice_head(n = 5) %>%
#   nest() %>%
#   deframe() %>% 
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]
# 
# names(T_top5_markers_sub) <- sub('.*\\.','',names(T_top5_markers_sub))
# 
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# 
# pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
# DotPlot(All_T, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_top5_markers_sub)+theme_90angle
# dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "T cell", 
  '1' = "T cell", 
  '2' = "Myeloid cell", 
  '3' = "T cell", 
  '4' = "Fibroblast",  
  '5' = "Fibroblast", 
  '6' = "T cell",  
  '7' = "T cell", 
  '8' = "Ductal cell", # CFTR
  '9' = "Endothelial cell", 
  '10' = "Exocrine glandular cell", # ERP27
  '11' = "T cell"
)


Idents(All_T) <- "RNA_snn_res.0.6"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)

## 映射-------------------------------
library(tibble)
library(data.table)
subtype <- All_T[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('T cell')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])

# saveRDS(All_T, file.path(T_out_dir, 'All_T.rds'))




# refiltering ------------------------------------------------------------------
All_T2 <- subset(All_T, subtype == 'T cell')


hai## PCA --------
DefaultAssay(All_T2) <- "RNA"
All_T2 <- FindVariableFeatures(All_T2)
All_T2 <- ScaleData(All_T2, vars.to.regress = 'nCount_RNA')
All_T2 <- RunPCA(All_T2)

T_out_dir2 <- 'T_subtype/refiltering' #dir.create(T_out_dir2)


### pca inspector ------
pdf(file.path(T_out_dir2, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T2, ndims = 50)
dev.off()

pdf(file.path(T_out_dir2, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(All_T2, dims = 1:5)
dev.off()
pdf(file.path(T_out_dir2, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(All_T2, dims = 6:10)
dev.off()
pdf(file.path(T_out_dir2, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(All_T2, dims = 11:15)
dev.off()
pdf(file.path(T_out_dir2, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(All_T2, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T2 <- RunUMAP(All_T2, dims = 1:dims,return.model = TRUE)



pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(All_T2, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


All_T2 <- All_T2 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T2 <- RunUMAP(
  All_T2, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T2  <-
  FindNeighbors(
    All_T2 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T2 <-
  FindClusters(
    All_T2,
    resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
  )



# dimplot
pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T2, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T2, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
pdf(file.path(T_out_dir2, str_glue("dotplot dims{dims} multi_res.pdf")), width = 35, height = 6)
DotPlot(All_T2, group.by = "RNA_snn_res.0.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.1", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.2", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "RNA_snn_res.0.4", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "RNA_snn_res.1", features = T_top5_markers)+RotatedAxis()
dev.off()

# find markers
Idents(All_T2) <- "RNA_snn_res.1"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T2),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

T_top5_markers <- 
  T_markers %>% 
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

names(T_top5_markers) <- sub('.*\\.','',names(T_top5_markers))



# dotplot
pdf(file.path(T_out_dir2, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T2, group.by = "RNA_snn_res.1", features = T_top5_markers)+RotatedAxis()
dev.off()

DotPlot(All_T2, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()

# Idents(All_T2) <- "RNA_snn_res.0.8"
# All_T2 <- FindSubCluster(All_T2, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.6)
# Idents(All_T2) <- "sub.cluster"
# # find markers
# T_markers_sub <-
#   FindAllMarkers(
#     JoinLayers(All_T2),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# 
# T_top5_markers_sub <- 
#   T_markers_sub %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 30) %>%
#   arrange(desc(pct.1)) %>%
#   slice_head(n = 5) %>%
#   nest() %>%
#   deframe() %>% 
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]
# 
# names(T_top5_markers_sub) <- sub('.*\\.','',names(T_top5_markers_sub))
# 
# DimPlot(All_T2, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(All_T2, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(All_T2, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# 
# pdf(file.path(T_out_dir2, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
# DotPlot(All_T2, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "sub.cluster", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "sub.cluster", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "sub.cluster", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "sub.cluster", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(All_T2, group.by = "sub.cluster", features = T_top5_markers_sub)+theme_90angle
# dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "T cell", 
  '1' = "T cell", 
  '2' = "T cell", 
  '3' = "T cell", 
  '4' = "T cell",   
  '5' = "T cell", 
  '6' = "T cell",  
  '7' = "T cell", 
  '8' = "T cell", 
  '9' = "T cell",  
  '10' = "Mast cell", #KIT
  '11' = "T cell", 
  '12' = "T cell", 
  '13' ="Exocrine glandular cell", # ERP27
  '14' = "T cell",  
  '15' = "Myeloid cell"
)


Idents(All_T2) <- "RNA_snn_res.1"

All_T2 <- RenameIdents(All_T2, cell_cluster)
All_T2$subtype <- Idents(All_T2)

saveRDS(All_T2, file.path(T_out_dir, 'All_T2.rds'))


# 映射-------------------------------
library(tibble)
library(data.table)
subtype <- All_T2[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation3 = ifelse((as.character(annotation2) %in% c('T cell')),
                         as.character(subtype), 
                         as.character(annotation2)))%>%
  dplyr::select(cells, annotation3) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data,meta.data.t['annotation3'])




# refiltering ------------------------------------------------------------------
All_T3 <- subset(All_T2, subtype == 'T cell')

pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T3, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T4  <-
  FindNeighbors(
    All_T4 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T4 <-
  FindClusters(
    All_T4,
    resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
  )



# dimplot
pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T4, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T4, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T4, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T4, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T4, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T4, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
pdf(file.path(T_out_dir2, str_glue("dotplot dims{dims} multi_res.pdf")), width = 35, height = 6)
DotPlot(All_T4, group.by = "RNA_snn_res.0.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.1", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.2", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.1", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T4, group.by = "RNA_snn_res.0.4", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T4, group.by = "RNA_snn_res.1", features = T_top5_markers)+RotatedAxis()
dev.off()

# find markers
Idents(All_T4) <- "RNA_snn_res.0.6"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T4),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

T_top5_markers <- 
  T_markers %>% 
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

names(T_top5_markers) <- sub('.*\\.','',names(T_top5_markers))

# T_out_dir2 <- "T_subtype/refiltering3/"

# dotplot
pdf(file.path(T_out_dir2, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = T_top5_markers)+RotatedAxis()
dev.off()

DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()


## annotation -----------
cell_cluster <- c(
  '0' = "T cell", 
  '1' = "T cell", 
  '2' = "T cell",  
  '3' = "T cell", 
  '4' = "T cell",  
  '5' = "T cell", 
  '6' = "T cell",  
  '7' = "T cell", 
  '8' = "T cell", 
  '9' = "Exocrine glandular cell"
)


Idents(All_T4) <- "RNA_snn_res.0.6"

All_T4 <- RenameIdents(All_T4, cell_cluster)
All_T4$subtype <- Idents(All_T4)

saveRDS(All_T4, file.path(T_out_dir, 'All_T4.rds'))




# 映射-------------------------------
library(tibble)
library(data.table)
subtype <- All_T4[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation4']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation5 = ifelse((as.character(annotation4) %in% c('T cell')),
                         as.character(subtype), 
                         as.character(annotation4)))%>%
  dplyr::select(cells, annotation5) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation5'])




# final -------------------------------------------------------------------------
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = c("CD4", "CD8A", "TRDC", "CD3D", "CD3E", "EPCAM", "ERP27", "ANPEP", "PECAM1", "CCDC80", "CD68", "CD163", "LILR4", "CD83", "FCGR3A", "CLEC9A", "CLEC10A", "LILRA4"), dot.scale = 10) + theme_90angle
DotPlot(All_T4, group.by = "RNA_snn_res.0.6", features = c("CCR7", "SELL", "IL7R", "GZMK", "CXCR5", "CXCR6", "ZNF683", "KLRG1", "PRF1", "FOXP3", "TOX", "TIGIT", "CD6"), dot.scale = 10) + theme_90angle
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = c("CD4", "CD8A", "TRDC", "CD3D", "CD3E", "CCR7", "SELL", "IL7R", "GZMK", "CXCR5", "CXCR6", "ZNF683", "KLRG1", "PRF1", "FOXP3", "TOX", "TIGIT", "CD6"), dot.scale = 10) + theme_90angle

T_out_dir2 <- "T_subtype/annotation" # dir.create(T_out_dir2)

All_T4 <- readRDS(file.path(T_out_dir, 'All_T4.rds'))
All_T_anno <- subset(All_T4, RNA_snn_res.0.6 %in% c(0:8))
All_T_anno$RNA_snn_res.0.6 %>% table()
All_T_anno$RNA_snn_res.0.6 <- as.character(All_T_anno$RNA_snn_res.0.6)
levels(All_T_anno$RNA_snn_res.0.6) <- factor(c(0:8))
All_T_anno$RNA_snn_res.0.6 %>% table()

# # Th17高表达IL-17、IL-22、IL-21，Tcm高表达IL-2（高增殖能力）
# pdf(file.path(T_out_dir2, "dotplot Th17 Tcm.pdf"), width = 8, height = 6)
# DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = c("IL17A", "IL17F", "IL22", "IL21", "IL2"), dot.scale = 10) + theme_90angle
# dev.off()


## annotation -----------
cell_cluster <- c(
  '0' = "CD8 Tem", 
  '1' = "CD4 Trm", 
  '2' = "CD8 T", 
  '3' = "CD8 Tnaive", 
  '4' = "gdT",  
  '5' = "CD4 Treg", 
  '6' = "CD8 Tem",  
  '7' = "Cycling T", 
  '8' = "CD8 Tem"
)

# VlnPlot(All_T, group.by = "RNA_snn_res.0.6", features = 'nCount_RNA', pt.size = 0)

Idents(All_T_anno) <- "RNA_snn_res.0.6"

All_T_anno <- RenameIdents(All_T_anno, cell_cluster)
All_T_anno$subtype <- Idents(All_T_anno)




cellType_marker_ls <- c(
  'T' = 'CD3E', 
  
  'CD8' = 'CD8A',
  
  'CD4' = 'CD4',
  
  'CD8 Tem' = 'GZMK', 
  'CD8 Tem' = 'GZMA',
  'CD8 Tem' = 'KLRG1',
  
  'CD4 Trm'='CXCR6',
  'CD4 Trm'='ITGAE',
  
  "CD8 Tnaive" = "SELL",
  "CD8 Tnaive" = "LEF1",
  "CD8 Tnaive" = "TCF7",
  "CD8 Tnaive" = "CCR7",
  
  "gdT" = "TRDC", # https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/cd8-t-cells-gdt-cells-and-mait-cells/
  
  'CD4 Treg'='FOXP3',
  'CD4 Treg'='CTLA4',
  'CD4 Treg'='TIGIT',
  
  'Cycling T'='MKI67',
  'Cycling T'='TOP2A'
)

pdf(file.path(T_out_dir2, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
DotPlot(All_T_anno, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(T_out_dir2,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(All_T_anno, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 20)
DimPlot(All_T_anno, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(All_T_anno, 'T_subtype/All_T_anno.rds')




# remove Tumor3 -----------------------------------------------------------------
All_T_anno <- readRDS('T_subtype/All_T_anno.rds')
All_T_anno$sample %>% table()

All_T_anno_remove_Tumor3 <- subset(All_T_anno, sample != 'Tumor3')
saveRDS(All_T_anno_remove_Tumor3, 'T_subtype/All_T_anno_remove_Tumor3.rds')

All_T_anno_remove_Tumor3$group <- factor(All_T_anno_remove_Tumor3$group, levels = c('N', 'CP', 'T'))
pdf(file.path(T_out_dir2, str_glue("Dimplot dim{dims} annotation split_rermove_Tumor3.pdf")), height = 7, width = 20)
DimPlot(All_T_anno_remove_Tumor3, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()




# 映射-------------------------------
library(tibble)
library(data.table)
subtype <- All_T_anno[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation5']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    T_annotation = ifelse((as.character(annotation5) %in% c('T cell')),
                         as.character(subtype), 
                         as.character(annotation5)))%>%
  dplyr::select(cells, T_annotation) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['T_annotation'])

all_data$T_annotation <- str_replace_all(all_data$T_annotation, c('Pancreatic β-cell' = 'Pancreatic-beta cell'))
pdf(file.path(T_out_dir2, str_glue("Dimplot all data.pdf")), height = 6, width = 10)
DimPlot(all_data, group.by = "T_annotation", label = T, repel = TRUE)
dev.off()




# spatial plot -----------------------------------------------------------------
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
all_data$T_annotation %>% table()
all_data$T_annotation_T <- fct_collapse(
  all_data$T_annotation, 
  'Others' = unique(all_data$T_annotation)
  [!unique(all_data$T_annotation) %in% 
      c('CD4 Treg', 'CD4 Trm', 'CD8 T', 'CD8 Tem', 'CD8 Tnaive', 'Cycling T', 'gdT')]
)
all_data$T_annotation_T %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')

cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')

p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y, color = T_annotation_T)) + 
      scattermore::geom_scattermore(
        pointsize = 1 , # 1
        alpha = 0.1, # 0.1
        # show.legend = c(color = TRUE, alpha = FALSE), 
        pixels = c(800, 800)
      ) +
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_color_manual(values = c('lightgrey', allcolour)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(T_out_dir2, 'spatial',  
                  str_glue("spatial dim plot {grp}_T highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  




# saveRDS(all_data, 'all_data_T_annotation.rds')
