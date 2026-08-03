library(Seurat)
library(tidyverse)

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

Epi_major_marker_ls <- c(
  "Tumor cell" = 'MIAT', 
  # "Tumor cell" = 'STRIP2', 
  # "Tumor cell" = 'CALB2', 
  # "Tumor cell" = 'MAN1A1', 
  # "Tumor cell" = 'EGLN3', 
  # "Tumor cell" = 'ASPM', 
  # "Tumor cell" = 'SPHK1', 
  # "Tumor cell" = 'SLC2A1', 
  # "Tumor cell" = 'SLC7A7', 
  
  "Anti-tumor cell" = 'FAM107A', 
  "Anti-tumor cell" = 'HSD17B2', 
  "Anti-tumor cell" = 'SERPINA3', 
  
  "Mucous cell" = 'MUC5B', 
  "Basal cell" = 'TP63', 
  "AT1 cell" = 'AGER', 
  
  # 'Mix' = 'POSTN', 
  # 'Mix' = 'EPCAM', 
  
  "Ciliated cell" = 'DLEC1', 
  "Ciliated cell" = 'DNAI1'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$annotation %>% table()

#T cells------------------------------------------------------------------------------------
table(all_data$group, all_data$annotation)
All_T <- subset(all_data, annotation == 'T cell' & group == 'SQ')


## PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')
All_T <- RunPCA(All_T)

T_out_dir<-'T_SQ_subtype' 
dir.create(T_out_dir)


### pca inspector ------
pdf(file.path(T_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 20

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
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )



# dimplot
pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
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

# find markers
Idents(All_T) <- "RNA_snn_res.1.5"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
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
pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = Epi_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.1.5", features = T_top5_markers)+RotatedAxis()
dev.off()


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
VlnPlot(All_T, features = 'nCount_RNA', pt.size = 0)
FeaturePlot(All_T, features = 'EPCAM', reduction = 'harmonyumap')
# cell_cluster <- c(
#   '0' = "",
#   '1' = "",
#   '2' = "",
#   '3' = "",
#   '4' = "",
#   '5' = "",
#   '6' = "",
#   '7' = "",
#   '8' = "",
#   '9' = "",
#   '10' = "",
#   '11' = "",
#   '12' = "",
#   '13' = "",
#   '14' = "",
#   '15' = "",
#   '16' = "",
#   '17' = ""
# )

cell_cluster <- c(
  '0' = "Tumor cell", 
  '1' = "CD4 T cell",
  '2' = "mix", # Fibro
  '3' = "CD8 Tem",
  '4' = "mix", # Macro
  '5' = "mix", # Fibro
  '6' = "CD4 Treg",
  '7' = "CD8 Tem",
  '8' = "mix", # Fibro
  '9' = "B cell",
  '10' = "mix", # Macro
  '11' = "Proliferating T cell",
  '12' = "CD8 Tnaive",
  '13' = "LAMP3 DC",
  '14' = "pDC",
  '15' = "Mast cell",
  '16' = "Neutrophil",
  '17' = "mix" # Tumor
)

Idents(All_T) <- "RNA_snn_res.1.5"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)

pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result.pdf")), height = 7, width = 18)
DimPlot(All_T, group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(subset(All_T, subtype %in% c("Tumor cell", "CD4 T cell", "CD8 Tem", "CD4 Treg", "B cell", "Proliferating T cell", "CD8 Tnaive", 'LAMP3 DC', 'Mast cell', 'pDC', "Neutrophil")), group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result - pure.pdf")), height = 7, width = 18)
DimPlot(subset(All_T, subtype %in% c("Tumor cell", "CD4 T cell", "CD8 Tem", "CD4 Treg", "B cell", "Proliferating T cell", "CD8 Tnaive", 'LAMP3 DC', 'Mast cell', 'pDC', "Neutrophil")), group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


cellType_marker_ls <- c(
  # "Epithelials" = "EPCAM", 
  # "Fibroblast" = 'POSTN', 
  # "Macrophage" = 'CD68', 
  'T' = 'CD3E',
  
  'CD8' = 'CD8A',
  
  'CD4' = 'CD4',
  
  'Tumor cell' = 'MIAT', 
  
  'CD8 Tem'='GZMK',
  'CD8 Tem'='GZMA',
  'CD8 Tem'='GZMB', 
  
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
 
  "B cell" = "CD19", 
  "B cell" = "MS4A1",

  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='STMN1',
  'proliferation'='ANAPC11',
  'proliferation'='TOP2A', 
  
  "CD8_Tnaive" = 'LEF1', 
  "CD8_Tnaive" = 'TCF7', 
  
  "LAMP3 DC" = 'CD83', 
  "LAMP3 DC" = 'LAMP3',  
  
  'pDC' = "CLEC4C", 
  'pDC' = "LILRA4",  
  
  "Mast cell" = 'KIT', 
  
  'Neutrophil' = 'CSF3R'
)

pdf(file.path(T_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
DotPlot(All_T, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(T_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 20, height = 8)
DotPlot(subset(All_T, subtype %in% c("Tumor cell", "CD4 T cell", "CD8 Tem", "CD4 Treg", "B cell", "Proliferating T cell", "CD8 Tnaive", 'LAMP3 DC', 'Mast cell', 'pDC', "Neutrophil")), group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(All_T, file.path(T_out_dir, 'All_T.rds'))
