result_dir <- 'result' 
dir.create(result_dir)

#read RDS----
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
# major_marker_ls <- c(
#   'epithelials'  = 'EPCAM',
#   'epithelials'  = 'KRT5',
#   'epithelials'  = 'KRT14',
#   
#   "Foveolar cells" = "GKN1", #ACT
#   "Foveolar cells" = "MUC5AC",
#   "Foveolar cells" = "TFF1",
#   
#   'endothelials' = 'PECAM1', 
#   'endothelials' = 'VWF', 
#   
#   "Fibroblasts" = "PDGFRA", 
#   "Fibroblasts" = "CCDC80",
#   "Fibroblasts" = "SFRP4",
#   "Fibroblasts" = "COL1A1",
#   "Fibroblasts" = "COL1A2",
#   "Fibroblasts" = "DCN",
#   "Fibroblasts" = "ACTA2",
#   
#   "B cells" = "CD19",
#   "B cells" = "CD79A",
#   "B cells" = "MS4A1",
#   "Plasma" = "TNFRSF17",
#   "Plasma" = "JCHAIN",
#   "Plasma" = "MZB1",
#   "Plasma" = "XBP1",
#   
#   'T cells'           = 'CD3D', 
#   'T cells'           = 'CD3E', 
#   'T cells'           = 'CD3G', 
#   'T cells'           = 'TRAC', 
#   
#   "NK" = "GNLY",
#   "NK" = "NKG7",
#   "NK" = "NCAM1",
#   
#   'Proliferating' = 'MKI67',
#   'Proliferating' = 'TOP2A',
#   
#   "Pericytes" = "PDGFRB",
#   "Pericytes" = "RGS5",
#   "Pericytes" = "MCAM",
#   
#   "Monocytes" = "CD14", 
#   "Monocytes" = "FCGR3A",
#   
#   "Macrophage" = "CD68",
#   "Macrophage" = "CD163",
#   "Macrophage" = "C1QA",
#   "Macrophage" = "C1QC", 
#   'Macrophage'  = 'CSF1R',
#   
#   "Dendritic cell" = "CD83",
#   'Dendritic cell' = "LILRA4", 
#   'pDCs' = "IRF4", 
#   'pDCs' = "IRF7", 
#   'pDCs' = "IRF8", 
#   'Dendritic cell' = "CLEC9A",
#   'Dendritic cell' = "IDO1",
#   'Dendritic cell' = "CD1C",
#   'Dendritic cell' = "FCER1A",
#   'Lamp3_DCs' = "LAMP3",
#   
#   'netrophil' = "CSF3R",
#   'netrophil' = "G0S2",
#   'netrophil' = "RETNLG",
#   
#   "Myocytes" = "TNNC1", 
#   "Myocytes" = "MYH7", 
#   
#   "Mast cell" = "CPA3",
#   "Mast cell" = "CTSG",
#   "Mast cell" = "TPSAB1",
#   "Mast cell" = "KIT"
# )

major_marker_ls <- c(
  # 系膜细胞
  "MES" = "Fhl2",
  "MES" = "Myl9",
  
  # 内皮细胞
  "EC" = "Pecam1",
  "EC" = "Vwf",
  "EC" = "Cldn5",
  "EC" = "Ackr1",
  
  # 足细胞
  "POD" = "Nphs2",
  "POD" = "Podxl",
  "POD" = "Ptpro",
  
  # 近端小管
  "PT" = "Akr1c21",
  "PT" = "Cubn",
  "PT" = "Slc13a1",
  "PT" = "Slc5a2",
  "PT" = "Lrp2",
  "PT" = "Aldob",
  "PT" = "Cryab",
  
  # 远曲小管
  "DT" = "Slc12a3",
  "DT" = "Calb1",
  "DT" = "Slc8a1",
  
  # 髓袢
  "LOH" = "Umod",
  "LOH" = "Slc12a1",
  "LOH" = "Cldn16",
  
  # 主细胞
  "PC" = "Aqp2",
  "PC" = "Aqp3",
  
  # 闰细胞
  "IC" = "Slc26a7",
  "IC" = "Slc4a1",
  "IC" = "Atp6v1g3",
  
  # 巨噬细胞
  "MC" = "Cd68",
  "MC" = "Mrc1",
  "MC" = "C1qa",
  "MC" = "C1qb",
  
  # 单核细胞
  "MON" = "Vcan",
  "MON" = "Cd14",
  
  # 树突状细胞
  "DC" = "Fcer1a",
  "DC" = "Clec10a",
  
  # 成纤维细胞
  "FIB" = "Lum",
  "FIB" = "Dcn",
  
  # 平滑肌 & 肌成纤维细胞（已合并去重）
  "SMC" = "Tagln",
  "SMC" = "Actg2",
  "SMC" = "Sostdc1",
  "SMC" = "Myh11",
  "SMC" = "Col1a1",
  "SMC" = "Col1a2",
  "SMC" = "Fn1",
  
  # 增殖细胞
  "CC" = "Top2a",
  "CC" = "Mki67"
)



valid_genes <- major_marker_ls[major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_genes), "/", length(major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(major_marker_ls, valid_genes)
missing_gene

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

#读取数据for映射----
xenium <- readRDS("xenium_FOV_image.rds")

All_T <- readRDS(file.path(T_out_dir, "all_T.rds"))
Fibroblast <- readRDS(file.path(Fibro_out_dir, "Fibroblast.rds"))
Macrophage <- readRDS(file.path(Macrophage_out_dir , "Macrophage.rds"))

All_T@meta.data %>% head()
All_T@meta.data$cell_id <- gsub("TMA_","",rownames(All_T@meta.data))
rownames(All_T@meta.data)<- paste0(All_T@meta.data$sample,"_",All_T@meta.data$cell_id)

Fibroblast@meta.data %>% head()
Fibroblast@meta.data$cell_id <- gsub("TMA_","",rownames(Fibroblast@meta.data))
rownames(Fibroblast@meta.data)<- paste0(Fibroblast@meta.data$sample,"_",Fibroblast@meta.data$cell_id)

Macrophage@meta.data %>% head()
Macrophage@meta.data$cell_id <- gsub("TMA_","",rownames(Macrophage@meta.data))
rownames(Macrophage@meta.data)<- paste0(Macrophage@meta.data$sample,"_",Macrophage@meta.data$cell_id)

###改行名
# rownames(xenium@meta.data) <- gsub("TMA_","",rownames(xenium@meta.data))
# xenium@meta.data$cell_id <- rownames(xenium@meta.data)
# rownames(xenium@meta.data) <- paste0(xenium@meta.data$sample,"_",xenium@meta.data$cell_id)

# xenium@meta.data$cell_id1 <- rownames(xenium@meta.data)
# xenium@meta.data$cell_id2 <- paste0(xenium@meta.data$orig.ident,"_",xenium@meta.data$cell_id)
# rownames(xenium@meta.data) <- xenium@meta.data$cell_id2

#映射——add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(All_T[['subtype']] %>% rownames_to_column('cells'), 
                 Fibroblast[['subtype']] %>% rownames_to_column('cells'),
                 Macrophage[['subtype']]%>% rownames_to_column('cells'))

meta.data.t <-
  subtype%>%
  right_join(
    xenium[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    Subtype = ifelse((as.character(annotation) %in% c('T cells','Fibroblasts',"Myeloids")),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, Subtype) %>%
  column_to_rownames('cells')


xenium<-AddMetaData(xenium,meta.data.t['Subtype'])
xenium$annotation %>% table()
xenium$Subtype %>% table()

# xenium$annotation2 <- NULL

# 确认是否加入
table(xenium$Subtype, xenium$annotation, useNA = "ifany")

saveRDS(xenium, 'result/xenium_Subtype.rds')


#T cells------------------------------------------------------------------------------------
All_T <- subset(all_data, celltype == 'T cells')

## PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')
All_T <- RunPCA(All_T)

T_out_dir<-'result/T_subtype' 
dir.create(T_out_dir)


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
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
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

valid_Tgenes <- T_major_marker_ls[T_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_Tgenes), "/", length(T_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_Tgene <- setdiff(T_major_marker_ls, valid_Tgenes)
missing_Tgene


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
valid_CD4Tgenes <- CD4_major_marker_ls[CD4_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_CD4Tgenes), "/", length(CD4_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_CD4Tgene <- setdiff(CD4_major_marker_ls, valid_CD4Tgenes)
missing_CD4Tgene


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
valid_CD8Tgenes <- CD8_major_marker_ls[CD8_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_CD8Tgenes), "/", length(CD8_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_CD8Tgene <- setdiff(CD8_major_marker_ls, valid_CD8Tgenes)
missing_CD8Tgene



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

T_top10_markers <- 
  T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(T_top10_markers) <- sub('.*\\.','',names(T_top10_markers))



#提取deepseek
T_top10_markers1 <-
  T_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)



# dotplot
pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = T_top10_markers)+RotatedAxis()
dev.off()

NK_major_marker_ls <- c('CD8'='FCGR3A',
                       'CD8'='NCAM1'
                      )

valid_NKgenes <- NK_major_marker_ls[NK_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_NKgenes), "/", length(NK_major_marker_ls))

DotPlot(All_T, group.by = "RNA_snn_res.0.6", features = NK_major_marker_ls)+RotatedAxis()


valid_Tgenes <- T_major_marker_ls[T_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_Tgenes), "/", length(T_major_marker_ls))

for (gene in valid_Tgenes) {
    p <- FeaturePlot(All_T, features = gene, reduction = 'harmonyumap')+RotatedAxis()
    ggsave(
      filename = file.path(result_dir, paste0(gene, "_feature.png")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
}

#DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls2)+RotatedAxis()

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
  '0' = "CD8_T", # CD8A GZMA
  '1' = "CD4_T", # CD4 TCF7 LEF1
  '2' = "Fib-like_T", #COL5A1
  '3' = "Mo-T", # CD4 CD8A CD68
  '4' = "Epi-like_T", # CD8 EPCAM
  '5' = "Pro_T", # MKI67 TOP2A
  '6' = "Treg", # CD4 FOXP3
  '7' = "CD4_T", # CD4 CXCL13
  '8' = "CD8_T", # CD8A
  '9' = "CD8_T", # CD8A
  '10' = "B-like_T" # SELL naive
)

Idents(All_T) <- "RNA_snn_res.0.6"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)

cellType_marker_ls <- c(
  'T' = 'CD3E',
  
  'CD8' = 'CD8A',
  'CD8' = 'CD8B',
  'CD8' = 'GZMH',
  
  'CD4' = 'CD4',
  'CD4' = 'STAT1',
  
  'Fiblike_T' = 'COL5A1',
  'Fiblike_T' = 'CCDC80',
  
  'MoT' = 'CD14',
  'MoT' = 'CD68',
  
  "Epi" = 'EPCAM',
  
  'Pro_T' = 'MKI67',
  'Pro_T' = 'TOP2A',
  
  'Treg'='FOXP3',
  'Treg'='CTLA4',
  
  'Blike_T'='CD19',
  'Blike_T'='CD79A'
)

pdf(file.path(T_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
DotPlot(All_T, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation split3.pdf")), height = 7, width = 18)
DimPlot(All_T, group.by = "subtype",split.by = 'group3', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
All_T_meta <- All_T[[]]

All_T_meta$sample <- factor(All_T_meta$sample, 
                            levels=c(All_T_meta$sample %>% unique())
)
p_col_proj <-
  All_T_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/T_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  All_T_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path("result/T_subtype/cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)


## volcano-------
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(T_markers,'result/T_subtype/T_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(T_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- T_markers[T_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

T_markers_ls <- list(T_markers)
names(T_markers_ls) <- 'T cells'

# scatter
p_volcanoplot_major_ls <- 
  T_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
T_de_genes_for_enrichmenet <- 
  T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  T_de_genes_for_enrichmenet %>%
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


er_output<-'result/T_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
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
      filename = file.path(er_output, name, "go dotplot.pdf"),
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
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

T_de_genes_id_for_enrichmenet <-
  T_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
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
  T_de_genes_id_for_enrichmenet %>%
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
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(All_T, 'result/T_subtype/all_T.rds')


# B cells ---------------------------------------------------------
b_out_dir <- "B_subtype" # dir.create(b_out_dir)
B <- subset(all_data, annotation %in% c('B cells', 'Plasma cells'))

## PCA ----------------------------------
DefaultAssay(B) <- "RNA"
B <- FindVariableFeatures(B)
B <- ScaleData(B, vars.to.regress = 'nCount_RNA')
B <- RunPCA(B)

### pca inspector ------
pdf(file.path(b_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B, ndims = 50)
dev.off()

pdf(file.path(b_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(B, dims = 1:5)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(B, dims = 6:10)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(B, dims = 11:15)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(B, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 10 

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

B <- RunUMAP(B, dims = 1:dims,return.model = TRUE)

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

B <- B %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

B <- RunUMAP(
  B, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
B <-
  FindNeighbors(
    B ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

B <-
  FindClusters(
    B,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
  )

# dimplot
pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# Idents(B)<-'RNA_snn_res.0.8'
# tmp <- FindSubCluster(B, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(tmp) <- "sub.cluster"
# DimPlot(tmp, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap')
# DotPlot(tmp, group.by = "sub.cluster", features = B_major_marker_ls) + RotatedAxis()
# B <- tmp

### dotplot -------
#https://zhuanlan.zhihu.com/p/25031149689
B_major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'naive_B' = 'CD20',
  'naive_B' = 'CD24',
  'naive_B' = 'IGHM',
  'naive_B' = 'BANK1',
  'naive_B' = 'FCRL1',
  'naive_B' = 'HLA-DQB1',
  'naive_B' = 'HLA-DQA1',
  'naive_B' = 'HVCN1',
  'naive_B' = 'IGHD',
  'naive_B' = 'FCER2',
  'naive_B' = 'IL4R',
  'Follicular B cells' = "CD40", 
  'Follicular B cells' = "BCL6", 
  'memory_B'= 'CD27',
  'memory_B'= 'FCRL4',
  'memory_B'= 'TACI',
  'memory_B'= 'AIM2',
  'memory_B'= 'TXNIP',
  'memory_B'= 'GPR183',
  'memory_B'= 'TNFRSF13B',
  'memory_B'= 'CD22',
  'memory_B'= 'FCRL2',
  'memory_B'= 'ADAM28',
  'memory_B'= 'IGLL5',
  'memory_B'= 'P2RX5',
  'memory_B'= 'BLK',
  # 'Marginal zone B cells' = "CD27", 
  'Marginal zone B cells' = "CR2",
  'Marginal zone B cells' = "CD1C",
  'cyclingB'='MKI67',
  'cyclingB'='STMN1',
  'cyclingB'='HMGB2',
  'cyclingB'='BIRC5',
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
  #'Plasmablasts' = "MKI67", 
  'Plasmablasts' = "CD38", 
  'Plasmablasts' = "IRF4", 
  'Bregs' = "IL10", 
  'Bregs' = "CD1D", 
  'Bregs' = "TGFB1", 
  'GC B'= 'S1PI2',
  # 'GC B'='BCL6',
  'GC B'='RGS13',
  'GC B'='TCL1A',
  'GC B'= 'LRMP',
  'GC B'= 'SUGCT',
  'GC B'= 'MME',
  'GC B'='AICDA'
)

Idents(B)<-'RNA_snn_res.1.2'
B_markers <-
  FindAllMarkers(
    JoinLayers(B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

B_top5_markers <- 
  B_markers %>% 
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

names(B_top5_markers) <- sub('.*\\.','',names(B_top5_markers))


# dotplot
pdf(file.path(b_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(B, group.by = "RNA_snn_res.1.2", features = major_marker_ls) + theme_90angle
DotPlot(B, group.by = "RNA_snn_res.1.2", features = B_major_marker_ls) + RotatedAxis()
DotPlot(B, group.by = "RNA_snn_res.1.2", features = B_top5_markers) + theme_90angle
DotPlot(B, group.by = "RNA_snn_res.1.2", features = Fibro_major_marker_ls2) + RotatedAxis()
DotPlot(B, group.by = "RNA_snn_res.1.2", features = CD4_major_marker_ls) + RotatedAxis()
DotPlot(B, group.by = "RNA_snn_res.1.2", features = CD8_major_marker_ls) + RotatedAxis()
dev.off()


## annotation -------------------------------------------------------------------------
# cell_cluster <- c(
#   '0_0'  = "Memory B cells", 
#   '0_1'  = "Naive B cells",
#   '1'  = "Plasma cells", # MZB1, XBP1
#   '2'  = "B cells", # no other markers
#   '3'  = "Plasma cells", # MZB1, XBP1
#   '4'  = "Plasma cells", # MZB1, XBP1
#   '5'  = "Epithelials", # Epcam
#   '6'  = "Marginal zone B cells", # CR2 https://www.science.org/doi/10.1126/science.abf7470
#   '7' = "rCAFs", # CCDC80, CCL19
#   '8' = "Plasma cells", 
#   '9' = "Plasma cells",  
#   '10' = "Macrophage-like B cells", # CD68, CD19
#   '11' = "Cycling B cells"
# )

#Idents(B) <- "RNA_snn_res.0.8"

cell_cluster <- c(
  '0' = "Plasma cells", #MZB1 ,XBP1
  '1' = "Memory B cells", 
  '2' = "Epithelials",
  '3' = "rCAFs",
  '4' = "Plasma cells",
  '5' = "Plasma cells",
  '6' = "Marginal zone B cells",
  '7' = "Plasma cells",
  '8' = "T-like B cells",
  '9' = "Naive B cells",
  '10' = "Plasma cells",
  '11' = "Macrophage-like B cells",
  '12' = "Cycling B cells",
  '13' = "Memory B cells", # B cells
  '14' = "Plasma cells",
  '15' = "vCAFs",
  '16' = "Marginal zone B cells",
  '17' = "Epithelials"
)

Idents(B) <- "RNA_snn_res.1.2"

B <- RenameIdents(B, cell_cluster)
B$subtype <- Idents(B)


pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 5, width = 11)
DimPlot(B, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

cellTye_markers <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'memory_B'= 'CD27',
  'memory_B'= 'AIM2',
  'memory_B'= 'TNFRSF13B',
  'memory_B'= 'CD22',
  'memory_B'= 'FCRL2',
  'memory_B'= 'ADAM28',
  'memory_B'= 'BLK',
  
  "Epithelials" = 'EPCAM',
  
  "rCAFs" = 'CCDC80', 
  "rCAFs" = 'CCL19', 
  
  'Marginal zone B cells' = "CR2",
  'Marginal zone B cells' = "CD1C",
  
  'T cell' = "CD3E", 
  
  'naive_B' = 'CD20',
  'naive_B' = 'CD24',
  'naive_B' = 'IGHM',
  'naive_B' = 'BANK1',
  'naive_B' = 'FCRL1',
  'naive_B' = 'HLA-DQB1',
  'naive_B' = 'HLA-DQA1',
  'naive_B' = 'HVCN1',
  'naive_B' = 'IGHD',
  'naive_B' = 'FCER2',
  'naive_B' = 'IL4R',
  
  'Macrophages' = "CD68", 
  'Macrophages' = "CD163", 
  'Macrophages' = "CSF1R", 
  
  'cyclingB'='MKI67',
  'cyclingB'='STMN1',
  'cyclingB'='HMGB2',
  'cyclingB'='BIRC5', 
  
  'vCAFs' = "NOTCH3", 
  'vCAFs' = "COL18A1", 
  'vCAFs' = "MCAM"
)

pdf(file.path(b_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(B, group.by = "subtype", features = cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


## cell props barplot  --------
B_meta <- B[[]]

B$sample <- factor(B$sample, levels = c(B_meta$sample %>% unique()))

p_col_proj <-
  B_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("B_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  B_meta %>%
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
  filename = file.path("B_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## volcano -------
B_markers <-
  FindAllMarkers(
    JoinLayers(B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(B_markers,'B_subtype/B_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(B_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- B_markers[B_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('B_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


B_markers_major <-
  B$subtype %>%
  unique %>%
  set_names(., .) %>%
  map(function(x){
    sc <- B[, B$subtype == x]
    sc <- sc %>%
      JoinLayers() %>%
      {Idents(.) <- "group"; .}
    
    print(str_glue('group: {x} cells: {ncol(sc)}'))
    FindAllMarkers(
      sc,
    )
  })

B_markers_ls <- list(B_markers)
names(B_markers_ls) <- 'B cells'

# scatter
p_volcanoplot_major_ls <- 
  B_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('B_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 18
    )
  })


## enrichment ------------------------------------------------------
B_de_genes_for_enrichmenet <- 
  B_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  B_de_genes_for_enrichmenet %>%
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

B_er_output<-'B_subtype/enrichment' # dir.create(B_er_output)
# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(B_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(B_er_output, name, "Go enrichment table.csv"),
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
      filename = file.path(B_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(B_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


B_de_genes_id_for_enrichmenet <-
  B_de_genes_for_enrichmenet %>%
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
  B_de_genes_id_for_enrichmenet %>%
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
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(B_er_output, 'enrich_ls.rds'))


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(B_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(B_er_output, name, "KEGG enriment table.csv"),
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
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(B_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(B_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(B, 'B_subtype/all_B.rds')




# Fibroblasts --------------------------------------------------------------
all_data$annotation %>% table()
all_data$res_0.8 %>% table()
table(all_data$res_0.8,all_data$annotation)

Fibro_out_dir<-'result/Fibroblast_subtype_fib_smc_proli' 
dir.create(Fibro_out_dir)
Fibroblast<-subset(all_data, subset = annotation %in% c("Fibroblast","SMC","Proliferate"))
table(Fibroblast$res_0.8,Fibroblast$annotation)

## PCA --------
DefaultAssay(Fibroblast) <- "RNA"
Fibroblast <- FindVariableFeatures(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)


### pca inspector ------
pdf(file.path(Fibro_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()

# pdf(file.path(Fibro_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 1:5)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 6:10)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 11:15)
# dev.off()
# pdf(file.path(Fibro_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

table(all_data2@meta.data$celltype, all_data2@meta.data$project)
  
Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path('result/Fibroblast_subtype_fib_smc_proli/', str_glue("Dimplot dim{dims} cluster harmonyumap.pdf")), height = 6, width = 8)
DimPlot(FIB1, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path('result/Fibroblast_subtype_fib_smc_proli/', str_glue("Dimplot dim{dims} cluster umap.pdf")), height = 6, width = 8)
DimPlot(FIB1, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(FIB1, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# saveRDS(Fibroblast, "result/Fibroblast_subtype_fib_smc_proli/Fibroblast.rds")
Fibroblast <- readRDS("result/Fibroblast_subtype_fib_smc_proli/Fibroblast.rds")

### dotplot -------
celltype_marker_ls <- c(
  # Epithelial
  'Epithelial' = "Epcam",
  'Epithelial' = "Krt19",
  'Epithelial' = "Cryab",
  
  # POD
  'POD' = "Nphs2",
  'POD' = "Podxl",
  
  # HealthyS1
  'PT' = "Slc34a1",
  'PT' = "Aldob",
  'PT' = "Lrp2",
  
  # LOH
  'LOH' = "Umod",
  'LOH' = "Slc12a1",
  
  # DT
  'DT' = "Slc12a3",
  'DT' = "Calb1",
  
  # PC
  'PC' = "Aqp2",
  'PC' = "Hsd11b2",
  
  # Immune
  'Immune' = "Ptprc",
  'Immune' = "Cd74",
  'Immune' = "Ccl5",
  
  # Fibroblast
  'Fibroblast' = "Dcn",
  'Fibroblast' = "Col1a1",
  
  # SMC
  'SMC' = "Acta2",
  'SMC' = "Tagln",
  
  # Adipocyte
  'Adipocyte' = "Fabp4",
  'Adipocyte' = "Plin4"
)



Fibro_major_marker_ls <- c(
  # S1 fibroblasts
  'S1_fibroblasts' = "Adamdec1", 
  'S1_fibroblasts' = "Cp", 
  'S1_fibroblasts' = "Fabp4",
  
  # S3 fibroblasts
  'S3_fibroblasts' = "Ogn", 
  'S3_fibroblasts' = "Ccdc80", 
  'S3_fibroblasts' = "Grem1",
  
  # Peri-cryptal S2 fibroblasts
  'S2_fibroblasts' = "Sox6",
  'S2_fibroblasts' = "F3",
  
  # S2a fibroblasts
  'S2a_fibroblasts' = "Vstm2",
  
  # S2b fibroblasts
  'S2b_fibroblasts' = "Npy",
  'S2b_fibroblasts' = "Nrg1",
  
  # Ribhi fibroblasts
  'Ribhi_fibroblasts' = "Rpl7a",
  'Ribhi_fibroblasts' = "Rpl28",
  
  # MT fibroblasts
  'MT_fibroblasts' = "Mt-genes",  # 线粒体基因通用标记
  
  # Immediate early-response (IER) fibroblasts
  'IER_fibroblasts' = "Fos", 
  'IER_fibroblasts' = "Fosb", 
  'IER_fibroblasts' = "Irf1",
  
  # Fibroblastic reticular cells (FRCs)
  'FRCs' = "Ccl19", 
  'FRCs' = "Ccl21",
  
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11",
  
  
  # Inflammatory fibroblasts
  'Inflammatory_fibroblasts' = "Il11",
  'Inflammatory_fibroblasts' = "Cxcl5", 
  'Inflammatory_fibroblasts' = "Fap", 
  'Inflammatory_fibroblasts' = "Inhba",
  'Inflammatory_fibroblasts' = "Il24",
  
  # Endothelial cells
  'Endothelial_cells' = "Vwf", 
  'Endothelial_cells' = "Pecam1", 
  'Endothelial_cells' = "Plvap",
  
  # Activated endothelium
  'Activated_endothelium' = "Ackr1", 
  'Activated_endothelium' = "Sele", 
  'Activated_endothelium' = "Selp",
  
  # Lymphatic endothelium
  'Lymphatic_endothelium' = "Mrc1", 
  'Lymphatic_endothelium' = "Mmrn1", 
  'Lymphatic_endothelium' = "Lyve1",
  
  # Pericytes
  'Pericytes' = "Notch3", 
  'Pericytes' = "Rgs5",
  
  # Glial cells
  'Glial_cells' = "Nrxn1", 
  'Glial_cells' = "S100b", 
  'Glial_cells' = "Cdh19"
)


Fibro_major_marker_ls1 <- c(
  # Smooth muscle cells (SMCs)
  'SMC' = "Tagln",
  'SMC' = "Mylk",
  'SMC' = "Cnn1",#'Acta2'
  
  # Fibroblasts (FIB)
  'FIB' = "Lum",
  'FIB' = "Dcn",
  'FIB' = "Col1a1",
  'FIB' = "Col1a2",
  'FIB' = "Fn1",
  'FIB' = "C7",
  
  # Myofibroblasts
  'Myofibroblasts' = "Actg2",
  'Myofibroblasts' = "Sostdc1",
  'Myofibroblasts' = "Myh11", #'Acta2',
  'Myofibroblasts' = "Notch3"
)


Idents(Fibroblast)<-'RNA_snn_res.0.6'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Fibroblast_top10_markers <- 
  Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Fibroblastrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Fibroblast_top10_markers) <- sub('.*\\.','',names(Fibroblast_top10_markers))

Fibroblast_top10_markers1 <-
  Fibroblast_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(Fibro_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = celltype_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_top10_markers)+RotatedAxis()
dev.off()

table(Fibroblast$RNA_snn_res.0.6,Fibroblast$group)

##修改分组--------
Fibroblast$group_raw <- Fibroblast$group

replace_dict <- c(
  "A1" = "H1",
  "A3" = "H12",
  "B1" = "Con",
  "C6" = "D7",
  "D2" = "H3",
  "F1" = "H24",
  "F2" = "H72"
)

Fibroblast$group <- plyr::revalue(Fibroblast$group, replace = replace_dict)
group_order <- c("Con", "H1", "H3", "H12", "H24", "H72", "D7")
Fibroblast$group <- factor(Fibroblast$group, levels = group_order)

table(Fibroblast$group_raw ,Fibroblast$group)

Fibroblast$group %>% table()
Fibroblast$annotation %>% table()

pdf(file.path(Fibro_out_dir, str_glue("Dotplot dims{dims} myofib_marker.pdf")), width = 5, height =5)
DotPlot(Fibroblast, group.by = "group", features = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'))+RotatedAxis()
dev.off()


pdf(file.path(Fibro_out_dir, str_glue("Vlnplot dims{dims} myofib_marker.pdf")), width = 20, height = 10)
VlnPlot(Fibroblast,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'),
        pt.size = 0,
        ncol = 4,
        )
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Vlnplot dims{dims} myofib_marker pt.pdf")), width = 20, height = 10)
VlnPlot(Fibroblast,group.by = "group",
        feature = c('Tagln', 'Actg2', 'Sostdc1', 'Myh11', 'Col1a1', 'Col1a2', 'Fn1'),
        pt.size = 0.1,
        ncol = 4,
)
dev.off()

saveRDS(Fibroblast, "result/Fibroblast_subtype_fib_smc_proli/Fibroblast.rds")
#Fibroblast <- readRDS("result/Fibroblast_subtype_fib_smc_proli/Fibroblast.rds")

p_col_proj <-
  Fibroblast[[]] %>%
  group_by(group, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         RNA_snn_res.0.6 = factor(RNA_snn_res.0.6,levels = c(0:13) )) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.6)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

table(Fibroblast$RNA_snn_res.0.6)

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} big anno.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "res_0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} big anno111.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "res_0.8", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dimplot
pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} oghhfdh.pdf")), height = 6, width = 20)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap',split.by = "group") + NoLegend()
#DimPlot(Fibroblast, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()


# Idents(Fibroblast) <- "RNA_snn_res.0.6"
# Fibroblast <- FindSubCluster(Fibroblast, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(Fibroblast) <- "sub.cluster"
# DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(Fibroblast, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(Fibroblast, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###score--------------------
# CAF <- c('FAP','PDGFRA','PDGFRB','HSF1','STAT3','MYC','CD10','GPR77')
# iCAF <- c('DPT','IGF1','ITM2A','CXCL12','CCL11','CD34','IL6')
# mCAF <-c('MMP11','COL11A1','POSTN','COL10A1')
# myCAF <- c('ACTA2','TAGLN','MMP11', 'MYL9','HOPX','POSTN','TPM1','TPM2')
# meCAF <- c('SLC2A1','BNIP3','ADM','ANGPTL4')
# pCAF <- c('NUSAP1','TOP2A','CENPF','PTTG1','SCG2','STMN1','PLAU','IGFBP2')
# vCAF <- c('NOTCH3','COL18A1','MCAM','RGS5','EPAS1')
# tCAF <- c('PDPN','MME','NDRG1','VEGFA')
# infCAF <- c('IDO1','CXCL9','CXCL10','CXCL11')
# rCAF <- c('CCL19')
# dCAF <- c('MKI67','SOX9','SOX10')
# 
# Macrophages <- c("CD68","CD163","C1QA","C1QC",'CSF1R')
# Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
# Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
# Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
# Fibroblasts <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
# Bcells <- c("CD19","CD79A","MS4A1")
# DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")
# 
# pathway<-list(CAF,iCAF,pCAF,mCAF,myCAF,meCAF,vCAF,tCAF,infCAF,rCAF,dCAF,Macrophages,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
# names(pathway)<-c('CAF','iCAF','pCAF','mCAF','myCAF','meCAF','vCAF','tCAF','infCAF','rCAF','dCAF','Macrophages','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')
# 
# # M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
# n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
# pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])
# 
# for ( i in seq_along(names(pathway))){
#   Fibroblast <- AddModuleScore(JoinLayers(Fibroblast), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
# }
# 
# p<-VlnPlot(Fibroblast,
#            features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
#                         'Macrophage1','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
#            group.by='RNA_snn_res.0.6',
#            combine = T,
#            pt.size = 0)
# 
# FeaturePlot(Fibroblast,
#             features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
#                          'Macrophage','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))
# 
# ggsave(
#   p,
#   filename = 'Fibroblast_subtype/CAF score.pdf',
#   height = 20,
#   width  = 15
# ) 

table(Fibroblast$RNA_snn_res.0.6,Fibroblast$group)

## annotation ----
cell_cluster <- c(
  '0' = "iCAF", #CXCL12 
  '1' = "mCAF", #POSTN COL10A1 COL11A1 CTHRC1 MMP11
  '2'  = "mCAF", #POSTN COL10A1 COL11A1 CTHRC1 MMP11
  '3'  = "mCAF", #POSTN CTHRC1 MMMP11
  '4'  = "Epi_like CAF", #肿瘤上皮细胞
  '5'  = "iCAF", # CXCL9 CXCL10 CXCL11 
  '6'  = "Fib", #成纤维CCDC80 SFRP4
  '7'  = "Mac_like CAF", # 巨噬细胞
  '8'  = "mCAF", #POSTN CTHRC1 MMMP11 
  '9'  = "Plasma_like CAF",#浆细胞
  '10' = "meCAF" # BNIP3 SLC2A1 ADM ANGPTL4
)

Idents(Fibroblast) <- "RNA_snn_res.0.6"
#Idents(Fibroblast) <- "sub.cluster"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)

Fib_marker_ls <- c(
  "iCAF" = 'CXCL12',
  "iCAF" = 'CXCL9',
  "iCAF" = 'CXCL10',
  "iCAF" = 'CXCL11',
  
  "mCAF" = 'POSTN',
  "mCAF" = 'COL10A1',
  "mCAF" = 'COL11A1',
  "mCAF" = 'CTHRC1',
  "mCAF" = 'MMP11',
  
  "Epi_like CAF" = 'EPCAM',
  "Epi_like CAF" = 'CDH1',
  "Epi_like CAF" = 'STMN1',
  
  "Fib" = 'CCDC80',
  "Fib" = 'SFRP4',
  
  "Mac_like CAF" = 'CD68',
  
  "Plasma_like CAF" = 'MZB1',
  "Plasma_like CAF" = 'XBP1',
  
  "meCAF" = "BNIP3" ,
  "meCAF" ='SLC2A1' ,
  "meCAF" ='ADM',
  "meCAF" ='ANGPTL4'
)

pdf(file.path(Fibro_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
DotPlot(Fibroblast, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Fibro_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibro_out_dir, str_glue("Dimplot dim{dims} annotation split3.pdf")), height = 6, width = 14)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group3', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Fibroblast, "result/Fibroblast_subtype/Fibroblast.rds")

## cell props barplot  --------
Fibroblast_meta <- Fibroblast[[]]

Fibroblast_meta$sample <- factor(Fibroblast_meta$sample, 
                            levels=c(Fibroblast_meta$sample %>% unique())
)
p_col_proj <-
  Fibroblast_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Fibroblast_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Fibroblast_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path("result/Fibroblast_subtype/cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)


# all CAFs subset ---------------------------------------------------------------
all_CAFs <- subset(Fibroblast, subtype %in% c('mCAF', 'iCAF', 'meCAF'))
all_CAFs$subtype <- as.character(all_CAFs$subtype)
all_CAFs$subtype <- factor(all_CAFs$subtype, levels = c(all_CAFs$subtype %>% unique()))

## PCA --------
DefaultAssay(all_CAFs) <- "RNA"
all_CAFs <- FindVariableFeatures(all_CAFs)
all_CAFs <- ScaleData(all_CAFs, vars.to.regress = 'nCount_RNA')
all_CAFs <- RunPCA(all_CAFs)


CAFs_out_dir <- "result/Fibroblast_subtype/CAFs_result" # 90
dir.create(CAFs_out_dir )

### pca inspector ------
pdf(file.path(CAFs_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_CAFs, ndims = 50)
dev.off()

pdf(file.path(CAFs_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_CAFs, dims = 1:5)
dev.off()
pdf(file.path(CAFs_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_CAFs, dims = 6:10)
dev.off()
pdf(file.path(CAFs_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_CAFs, dims = 11:15)
dev.off()
pdf(file.path(CAFs_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_CAFs, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

all_CAFs <- RunUMAP(all_CAFs, dims = 1:dims,return.model = TRUE)

pdf(file.path(CAFs_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(all_CAFs, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

all_CAFs <- all_CAFs %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_CAFs <- RunUMAP(
  all_CAFs, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(CAFs_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(all_CAFs, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
all_CAFs  <-
  FindNeighbors(
    all_CAFs ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_CAFs <-
  FindClusters(
    all_CAFs,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(CAFs_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(all_CAFs, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_CAFs, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_CAFs, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_CAFs, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(all_CAFs, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(all_CAFs, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(all_CAFs, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(all_CAFs)<-'RNA_snn_res.0.6'
all_CAFs_markers <-
  FindAllMarkers(
    JoinLayers(all_CAFs),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_CAFs_top5_markers <- 
  all_CAFs_markers %>% 
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
        # all_CAFsrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(all_CAFs_top5_markers) <- sub('.*\\.','',names(all_CAFs_top5_markers))

fibro_markers <- c(
  Epithelial  = "EPCAM", 
  # Epithelial  = "KRT5",
  Epithelial  = "KRT8",
  Epithelial  = "S100A8",
  Epithelial  = "KRT18",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  'mCAF'='MMP11',
  'mCAF'='COL1A1',
  'mCAF'='POSTN',
  'mCAF'='COMP',
  'mCAF'='COL10A1',
  'mCAF'='COL11A1',
  'mCAF'='DCN',
  'mCAF'='CDH11',
  'mCAF'='LUM',
  'iCAF'='PLA2G2A',
  'iCAF'='CFD',
  'iCAF'='CD34',
  'iCAF'='CXCL14',
  'iCAF'='PI16',
  'iCAF'='PRG4',
  'iCAF'='APOD',
  'iCAF'='DPT',
  # 'iCAF'='CXCL12',
  'iCAF'='CXCL1',
  'iCAF'='IL6',
  'iCAF'='CCL2',
  'apCAF'='HLA-DRA',
  'apCAF'='HLA-DRB1',
  'apCAF'='CD74',
  'apCAF'='IGLC2',
  'apCAF'='IGKC',
  'apCAF'='CXCR4',
  'vCAF'='ACTA2',
  'vCAF'='MCAM',
  'vCAF'='RGS5',
  'vCAF'='MYH11',
  'vCAF'='ADIRF',
  'tCAF'='ENO1',
  'tCAF'='GAPDH',
  'tCAF'='NDRG1',
  # 'tCAF'='HSPH1',
  # 'tCAF'='VEGFA',
  'tCAF'='MME',
  'dCAF'='MKI67',
  'dCAF'='STMN1',
  'dCAF'='PTTG1',
  'dCAF'='H2AFZ',
  'dCAF'='CENPF',
  'infCAF'='IDO1',
  'infCAF'='CXCL9',
  'infCAF'='CXCL10',
  'infCAF'='CXCL11',
  'rCAF'='CCL19',
  'rCAF'='CCL21',
  'pericyte'='COL4A1',
  'pericyte'='COL4A2',
  'hsptcaf'='HSPH1',
  'hsptcaf'='VEGFA'
)

# dotplot
pdf(file.path(CAFs_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls1)+RotatedAxis()
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls3)+RotatedAxis()
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = fibro_markers)+RotatedAxis()
DotPlot(all_CAFs, group.by = "RNA_snn_res.0.6", features = all_CAFs_top5_markers)+RotatedAxis()
dev.off()

# Idents(all_CAFs) <- "RNA_snn_res.0.6"
# all_CAFs <- FindSubCluster(all_CAFs, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(all_CAFs) <- "sub.cluster"
# DimPlot(all_CAFs, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(all_CAFs, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(all_CAFs, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(all_CAFs, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###score--------------------
CAF <- c('FAP','PDGFRA','PDGFRB','HSF1','STAT3','MYC','CD10','GPR77')
iCAF <- c('DPT','IGF1','ITM2A','CXCL12','CCL11','CD34','IL6')
mCAF <-c('MMP11','COL11A1','POSTN','COL10A1')
myCAF <- c('ACTA2','TAGLN','MMP11', 'MYL9','HOPX','POSTN','TPM1','TPM2')
meCAF <- c('SLC2A1','BNIP3','ADM','ANGPTL4')
pCAF <- c('NUSAP1','TOP2A','CENPF','PTTG1','SCG2','STMN1','PLAU','IGFBP2')
vCAF <- c('NOTCH3','COL18A1','MCAM','RGS5','EPAS1')
tCAF <- c('PDPN','MME','NDRG1','VEGFA')
infCAF <- c('IDO1','CXCL9','CXCL10','CXCL11')
rCAF <- c('CCL19')
dCAF <- c('MKI67','SOX9','SOX10')

Macrophages <- c("CD68","CD163","C1QA","C1QC",'CSF1R')
Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
all_CAFss <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(CAF,iCAF,pCAF,mCAF,myCAF,meCAF,vCAF,tCAF,infCAF,rCAF,dCAF,Macrophages,Neutrophils,Plasma,Tcells,all_CAFss,Bcells,DC)
names(pathway)<-c('CAF','iCAF','pCAF','mCAF','myCAF','meCAF','vCAF','tCAF','infCAF','rCAF','dCAF','Macrophages','Neutrophils','Plasma','Tcells','all_CAFss','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  all_CAFs <- AddModuleScore(JoinLayers(all_CAFs), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(all_CAFs,
           features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                        'Macrophage1','Neutrophils1','Plasma1','Tcells1','all_CAFss1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           combine = T,
           pt.size = 0)

FeaturePlot(all_CAFs,
            features = c('CAF1','iCAF1','pCAF1','mCAF1','myCAF1','meCAF1','vCAF1','tCAF1','infCAF1','rCAF1','dCAF1',
                         'Macrophage','Neutrophils1','Plasma1','Tcells1','all_CAFss1','Bcells1','DC1'))

ggsave(
  p,
  filename = file.path(CAFs_out_dir, 'CAF score.pdf'),
  height = 20,
  width  = 15
) 


## annotation ----
cell_cluster <- c(
  '0'  = "myCAFs",
  '1' = "iCAFs", #
  '2'  = "CAFs", # https://pubmed.ncbi.nlm.nih.gov/39863019/
  '3'  = "CAFs2", # 
  '4'  = "vCAFs",
  '5'  = "vCAFs", # 
  '6'  = "tCAFs",
  '7'  = "rCAFs",
  '8'  = "CAFs", 
  '9'  = "CAFs2",
  '10' = "dCAFs", 
  '11'  = "CAFs2"
)

Idents(all_CAFs) <- "RNA_snn_res.0.6"
#Idents(all_CAFs) <- "sub.cluster"

all_CAFs <- RenameIdents(all_CAFs, cell_cluster)
all_CAFs$subtype <- Idents(all_CAFs)

pdf(file.path(CAFs_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(all_CAFs, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(all_CAFs, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(CAFs_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 14)
DimPlot(all_CAFs, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



Fibroblast_cellTye_markers <- c(  
  'myCAF'='MMP11',
  'myCAF'='CTHRC1',
  'myCAF'='COL11A1',
  'myCAF'='POSTN',
  'myCAF'='COL10A1',
  'myCAF'='TNC',
  
  'iCAF'='DPT',
  'iCAF'='IGF1',
  'iCAF'='ITM2A',
  'iCAF'='CXCL12',
  'iCAF'='CCL11',
  'iCAF'='CD34',
  'iCAF'='IL6',
  
  'CAF' = 'FAP',
  'CAF' = 'PDGFRA',
  'CAF' = 'PDGFRB',
  'CAF' = 'HSF1',
  'CAF' = 'STAT3',
  'CAF' = 'MYC',
  
  'vCAF' = 'NOTCH3',
  'vCAF' = 'COL18A1',
  'vCAF' = 'MCAM',
  'vCAF' = 'RGS5',
  'vCAF' = 'EPAS1',
  
  'tCAF' = 'PDPN',
  'tCAF' = 'MME',
  'tCAF' = 'NDRG1',
  'tCAF' = 'VEGFA',
  
  'rCAF'='CCL19',
  
  'dCAF' = 'MKI67',
  'dCAF' = 'SOX9',
  'dCAF' = 'SOX10'
)

pdf(file.path(CAFs_out_dir, "dotplot celltype main_marker.pdf"), width = 28, height = 6)
DotPlot(all_CAFs, group.by = "subtype", features = Fibroblast_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()




## cell props barplot  --------
all_CAFs_meta <- all_CAFs[[]]

all_CAFs_meta$sample <- factor(all_CAFs_meta$sample,levels=all_CAFs_meta$sample %>% unique())
p_col_proj <-
  all_CAFs_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(CAFs_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  all_CAFs_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(CAFs_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)




## enrichment ------------------
all_CAFs_markers <-
  FindAllMarkers(
    JoinLayers(all_CAFs),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(all_CAFs_markers,file.path(CAFs_out_dir, 'all_CAFs_diff_gene.csv'))

library(EnhancedVolcano)
type<-unique(all_CAFs_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- all_CAFs_markers[all_CAFs_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(CAFs_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(all_CAFs = all_CAFs_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(CAFs_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })


all_CAFs_de_genes_for_enrichmenet <- 
  all_CAFs_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  all_CAFs_de_genes_for_enrichmenet %>%
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


all_CAFs_er_output <- file.path(CAFs_out_dir, 'enrichment') # dir.create(all_CAFs_er_output)
# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(all_CAFs_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(all_CAFs_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(all_CAFs_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(all_CAFs_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

all_CAFs_de_genes_id_for_enrichmenet <-
  all_CAFs_de_genes_for_enrichmenet %>%
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
  all_CAFs_de_genes_id_for_enrichmenet %>%
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
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(all_CAFs_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(all_CAFs_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(all_CAFs_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(all_CAFs_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(all_CAFs_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })
# saveRDS(all_CAFs, "Fibroblast_subtype/all_CAFs.rds")




## 映射-------------------------------
library(tibble)
library(data.table)
sub <- all_CAFs[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    Fibroblast[['subtype']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype2 = ifelse((as.character(subtype.x) %in% c('myCAFs', 'iCAFs', 'CAFs', 'CAFs2', 'vCAFs', 'tCAFs', 'rCAFs', 'dCAFs')),
                      as.character(subtype.x), 
                      as.character(subtype.y)))%>%
  dplyr::select(cells, subtype2) %>%
  column_to_rownames('cells')


Fibroblast <- AddMetaData(Fibroblast,meta.data.t['subtype2'])
Fibroblast$subtype_raw <- Fibroblast$subtype
Fibroblast$subtype <- Fibroblast$subtype2
saveRDS(Fibroblast, 'Fibroblast_subtype/all_Fibroblast.rds')




# Myeloids/Macrophages -------------------------------------------------------------------
all_data$annotation <- str_replace(all_data$annotation, 'Monocytes', 'Macrophages')

Macrophage_out_dir <- "result/Macrophage_subtype" 
dir.create(Macrophage_out_dir)
Macrophage <- subset(all_data, annotation == 'Myeloids')

## PCA --------
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

pdf(file.path(Macrophage_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 1:5)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 6:10)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 11:15)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Macrophage) <- "RNA_snn_res.0.6"


### dotplot -------
macrophage_major_marker_ls <- c(
  M1 = 'NOS2',
  M1 = 'IL12',
  M1 = 'FCGR1A',
  M1 = 'FCGR1B',
  M1 = 'FCGR1C',
  M1 = 'CD80',
  M1 = 'CXCR10',
  M1 = 'IL23',
  M1 = 'CXCL9',
  M1 = 'CXCL10',
  M1 = 'CXCL11',
  M1 = 'CD86',
  M1 = 'IL1A',
  M1 = 'IL1B',
  M1 =  'IL6',
  M1 =  'TNF',
  M1 = 'MHCII',
  M1 = 'CCL5',
  M1 = 'IRF5',
  M1 = 'IRF1',
  M1 = 'CD40',
  M1 = 'KYNU',
  M1 = 'CCR7',
  M1 =  'VCAN',
  M1 =  'CD300E',
  M1 =  'INHBA',
  M1 =  'TNIP3',
  M2 = 'ARG1',
  M2 = 'IL10',
  M2 = 'CD32',
  M2 = 'CD163',
  M2 = 'FCER2',
  M2 = 'CD200R1',
  M2 = 'PDCD1LG2',
  M2 = 'CD274',
  M2 = 'MARCO',
  M2 = 'CSF1R',
  M2 = 'MRC1',
  M2 = 'IL1RN',
  M2 = 'Il1R2',
  M2 = 'IL4R',
  M2 = 'CCL4',
  M2 = 'CCL13',
  M2 = 'CCL20',
  M2 = 'CCL17',
  M2 = 'CCL18',
  M2 = 'CCL22',
  M2 = 'CCL24',
  M2 = 'LYVE1',
  M2 = 'VEGFA',
  M2 = 'VEGFB',
  M2 = 'VEGFC',
  M2 = 'VEGFD',
  M2 = 'EGF',
  M2 = 'CTSA',
  M2 = 'CTSB',
  M2 = 'CSTC',
  M2 = 'CTSD',
  M2 = 'TGFB1',
  M2 = 'CD209',
  M2 = 'CD163L1',
  M2 = 'FOLR2',
  M0 = 'CD68',
  M0 = 'C1QA',
  M0 = 'C1QB',
  M0 = 'SELENOP',
  M0 = 'LYZ',
  M0 = 'HLA-DPB1',
  M0 = 'HLA-DPA1',
  M0 = 'AIF1', 
  TAM = 'SPP1', 
  TAM = 'CQC'
)


Myeloids_major_marker_ls <- c(
  # 单核细胞亚群
  'Classical Monocyte' = "CD14",
  'Classical Monocyte' = "S100A8",
  'Non-classical Monocyte' = "FCGR3A",  # CD16
  
  # 巨噬细胞亚群
  'M1-lMacrophage' = "CD86",
  'M1-Macrophage' = "CXCL9",
  'M2-Macrophage' = "CD163",
  'M2-Macrophage' = "MSR1",
  'Pro-angiogenic TAM' = "SPP1",       # 促血管生成型[5](@ref)
  'Anti-inflammatory TAM' = "C1QC",    # 抗炎型[5](@ref)
  'Anti-inflammatory TAM' = "APOE",
  
  # 树突状细胞亚群
  "Dendritic" = "FCER1A",
  "Dendritic" = "JCHAIN",
  'cDC1' = "CLEC9A",                   # 抗原呈递型[2,5](@ref)
  'cDC1' = "XCR1",
  'cDC2' = "CD1C",
  'LAMP3+ DC' = "LAMP3",               # 迁移型[5](@ref)
  
  # 粒细胞亚群
  'Neutrophil' = "CSF3R",
  'Neutrophil' = "CXCR2",
  'Eosinophil' = "SIGLEC8",
  'Basophil' = "CD123",
  
  # 肥大细胞
  'Mast' = "KIT",                 # 干细胞因子受体[2](@ref)
  'Mast' = "TPSAB1",              # 类胰蛋白酶[2](@ref)
  "Mast" = "MS4A2",
  "Mast" = "CPA3",
  
  # 髓系祖细胞
  'CMP' = "CD34",                      # 共同髓系祖细胞[1](@ref)
  'GMP' = "CD33"                       # 粒-单核祖细胞[4](@ref)
)

myeloids_cell_markers <- list(
  "Mast_KIT" = c("KIT", "TPSAB1", "CPA3"),
  "pDC_LILRA4" = c("LILRA4", "GZMB", "IL3RA"),
  "cDC1_CLEC9A" = c("CLEC9A", "FLT3", "IDO1"),
  "cDC2_CD1C" = c("CD1C", "FCER1A", "HLA-DQA1"),
  "cDC3_LAMP3" = c("LAMP3", "CCFR7", "FSCN1"),
  "Mono_CD14" = c("FCN1", "S100A8", "S100A9"),
  "Mono_CD16" = c("FCGR3A", "LST1", "LILRB2"),
  "Macro_INHBA" = c("INHBA","IL1RN", "CCL4"),
  "Macro_NLRP3" = c("NLRP3", "EREG", "IL1B"),
  "Macro_LYVE1" = c("LYVE1", "PLTP", "SEPP1"),
  "Macro_C1QC" = c("C1QC", "C1QA", "APOE")
)


Idents(Macrophage)<-'RNA_snn_res.0.6'
macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

macrophage_top10_markers <- 
  macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(macrophage_top10_markers) <- sub('.*\\.','',names(macrophage_top10_markers))


#提取deepseek
macrophage_top10_markers1 <- 
  macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) 


# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 27, height = 6)
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Myeloids_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = macrophage_top10_markers)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = myeloids_cell_markers)+RotatedAxis()
dev.off()


###score--------------------
M1<-c( 'NOS2',
       'IL12',
       'FCGR1A',
       'FCGR1B',
       'FCGR1C',
       'CD80',
       'CXCR10',
       'IL23',
       'CXCL9',
       'CXCL10',
       'CXCL11',
       'CD86',
       'IL1A',
       'IL1B',
       'IL6',
       'TNF',
       'MHCII',
       'CCL5',
       'IRF5',
       'IRF1',
       'CD40',
       'KYNU',
       'CCR7',
       'VCAN',
       'CD300E',
       'INHBA',
       'TNIP3')
M2<-c( 'ARG1',
       'IL10',
       'CD32',
       'CD163',
       'FCER2',
       'CD200R1',
       'PDCD1LG2',
       'CD274',
       'MARCO',
       'CSF1R',
       'MRC1',
       'IL1RN',
       'Il1R2',
       'IL4R',
       'CCL4',
       'CCL13',
       'CCL20',
       'CCL17',
       'CCL18',
       'CCL22',
       'CCL24',
       'LYVE1',
       'VEGFA',
       'VEGFB',
       'VEGFC',
       'VEGFD',
       'EGF',
       'CTSA',
       'CTSB',
       'CSTC',
       'CTSD',
       'TGFB1',
       'CD209',
       'CD163L1',
       'FOLR2')
Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
Fibroblasts <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(M1,M2,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('M1','M2','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Macrophage,
           features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           pt.size = 0)

FeaturePlot(Macrophage,
            features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))

ggsave(
  p,
  filename = 'result/Macrophage_subtype/M1_M2 score.pdf',
  height = 10,
  width  = 15
) 


## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0' = 'M0-Mac', #CD68
  '1' = 'M2-Mac', #CD163 MRC1
  '2' = 'Caf-like_MoMa', #POSTN COL5A1
  '3' = 'cDC2', # CD1C, FCER1A CLEC10A
  '4' = 'Epi-like_MoMa',#EPCAM 
  '5' = 'M1-Mac', #CXCL9 CXCL10
  '6' = 'Pro_MoMa',# MKI67 TOP2A
  '7' = 'Mono', # FCGBP, CDC42
  '8' = 'Mast',# KIT MS4A2 CTSG
  '9' = 'Plasma-like_MoMa',#MZB1 XBP1
  '10' = 'LAMP3_DC',#LAMP3, FSCN1
  '11' = 'pDC',# LILRA4, CLEC4C, IL3RA
  '12' = 'cDC1' # CLEC9A, IDO1, XCR1
)

Idents(Macrophage) <- "RNA_snn_res.0.6"

Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split3.pdf")), height = 6, width = 20)
DimPlot(Macrophage, group.by = "subtype",split.by = 'group3', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# Macrophage_cellTye_markers <- c( 
#   "Macrophage" = "CD68",
#   # "Macrophage" = "CD163",
#   "Macrophage" = "C1QA",
#   "Macrophage" = "C1QC", 
#   # 'Macrophage'  = 'CSF1R',
#   
#   "Fibroblasts" = "PDGFRA", 
#   "Fibroblasts" = "CCDC80",
#   "Fibroblasts" = "POSTN",
#   # "Fibroblasts" = "SFRP4",
#   "Fibroblasts" = "COL1A1",
#   "Fibroblasts" = "DCN",
#   "Fibroblasts" = "ACTA2",
#   
#   'vCAFs' = "MCAM", 
#   'vCAFs' = "NOTCH3", 
#   'vCAFs' = "COL18A1",
#   
#   'M2 macrophages' = "CD163",
#   'M2 macrophages' = "CSF1R",
#   'M2 macrophages' = "MRC1",
#   'M2 macrophages' = "CCL24",
#   'M2 macrophages' = "VEGFB",
#   'M2 macrophages' = "CTSA",
#   'M2 macrophages' = "CD209",
#   'M2 macrophages' = "FOLR2",
#   
#   'Epithelials' = "EPCAM",
#   
#   'Epithelial-like macrophages' = 'PGA3',
#   
#   'cDC2s' = 'CD1C',
#   'cDC2s' = 'CLEC10A',
#   
#   'Proliferating' = 'MKI67',
#   'Proliferating' = 'TOP2A',
#   
#   'LAMP3 DCs' = "LAMP3",
#   'LAMP3 DCs' = 'CD83', 
#   
#   'IFN-responsive macrophages' = "IFIT3", 
#   'IFN-responsive macrophages' = "IFIT2",
#   'IFN-responsive macrophages' = "IFIT1",
#   
#   'cDC1s' = 'CLEC9A',
#   'cDC1s' = 'IDO1'
# )

Macrophage_cellTye_markers <- c( 
  "M0-Mac" = "CD68",
  
  # "Macrophage" = "CD163",
  "M2-Mac" = "CD163",
  "M2-Mac" = "MRC1",
  # 'Macrophage'  = 'CSF1R',
  
  "Caf-like_MoMa" = "POSTN ",
  "Caf-like_MoMa" = "COL5A1 ",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A",
  
  'Epi-like_MoMa' = "EPCAM",
  
  'M1-Mac' = "CXCL9",
  'M1-Mac' = "CXCL10",
 
  'Pro_MoMa' = 'MKI67',
  'Pro_MoMa' = 'TOP2A',
  
  'Mono' = "FCGBP",
  'Mono' = "CDC42",
  
  "MAST" = "KIT",
  "MAST" = "MS4A2",
  
  'Plasma-like_MoMa' = 'MZB1',
  'Plasma-like_MoMa' = 'XBP1',

  'LAMP3_DC' = "LAMP3",
  'LAMP3_DC' = 'FSCN1', 
  
  'pDC' = 'CLEC4C',
  'pDC' = 'LILRA4',
  
  'cDC1s' = 'CLEC9A',
  'cDC1s' = 'IDO1'
)

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(Macrophage, group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) 
#+ theme_90angle
dev.off()


## cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage[[]]

Macrophage_meta$sample <- factor(Macrophage_meta$sample,levels=c(Macrophage_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Macrophage_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Macrophage_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)

#Only_Macrophage----
# 生成子集对象（仅保留三种目标亚型）
target_subtypes <- c('M0-Mac', 'M1-Mac', 'M2-Mac')
Macrophage_subset <- subset(Macrophage, subtype %in% target_subtypes)

# 绘制分组UMAP图（三组对比）
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation filtered.pdf")), 
    height = 6, width = 8)  # 宽度从20调整为15（更紧凑）
DimPlot(
  Macrophage_subset, 
  group.by = "subtype",
  #split.by = 'group3', 
  reduction = 'harmonyumap',
  label = TRUE,
  repel = TRUE
) + NoLegend()
dev.off()

# 绘制分组UMAP图
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split3_filtered.pdf")), 
    height = 6, width = 15)  # 宽度从20调整为15（更紧凑）
DimPlot(
  Macrophage_subset, 
  group.by = "subtype",
  split.by = 'group3', 
  reduction = 'harmonyumap',
  label = TRUE,
  repel = TRUE
) + NoLegend()
dev.off()


Macrophage_subset_meta <- Macrophage_subset[[]]

Macrophage_subset_meta$sample <- factor(Macrophage_subset_meta$sample,levels=c(Macrophage_subset_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_subset_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Macrophage_subtype", "Macro_cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Macrophage_subset_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "Macro_cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)





## vocanol -------------------
Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Macrophage_markers,'Macrophage_subtype/Macrophage_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(Macrophage_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Macrophage_markers[Macrophage_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })

volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


Macrophage_markers_major<-list(Macrophage = Macrophage_markers)
# scatter
p_volcanoplot_major_ls <- 
  Macrophage_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 25
    )
  })


## enrrichment ---------------------------------------
Macrophage_de_genes_for_enrichmenet <- 
  Macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

Macrophage_er_output<-'Macrophage_subtype/enrichment' # dir.create(Macrophage_er_output)
enrich_ls <-
  Macrophage_de_genes_for_enrichmenet %>%
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
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(Macrophage_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Macrophage_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Macrophage_de_genes_id_for_enrichmenet <-
  Macrophage_de_genes_for_enrichmenet %>%
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
  Macrophage_de_genes_id_for_enrichmenet %>%
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
      use_internal_data = FALSE
    )
  })

# saveRDS(Macrophage_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Macrophage_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Macrophage_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(Macrophage, 'result/Macrophage_subtype/Macrophage.rds')




# 映射-------------------------------
library(tibble)
library(data.table)
subtype <- rbind(All_T[['subtype']] %>% rownames_to_column('cells'), 
                 B[['subtype']] %>% rownames_to_column('cells'), 
                 Macrophage[['subtype']] %>% rownames_to_column('cells'), 
                 Fibroblast[['subtype']] %>% rownames_to_column('cells'))
meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('T cells','B cells', 'Plasma cells', 'Macrophages','Fibroblasts')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
all_data$annotation %>% table()
all_data$annotation2 %>% table()




# Epithelials -------------------------------------------------------------------
Epithelial <- subset(all_data, annotation2 %in% c('Epithelials', 'Foveolar cells'))

## PCA --------
DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'Epithelial_subtype' # dir.create(Epithelial_out_dir)


### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()

pdf(file.path(Epithelial_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Epithelial, dims = 1:5)
dev.off()
pdf(file.path(Epithelial_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Epithelial, dims = 6:10)
dev.off()
pdf(file.path(Epithelial_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Epithelial, dims = 11:15)
dev.off()
pdf(file.path(Epithelial_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Epithelial, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)



pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


Epithelial <- Epithelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Epithelial <- RunUMAP(
  Epithelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Epithelial  <-
  FindNeighbors(
    Epithelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
# 1. Epi2: Xiao S, Zhou L. Gastric Stem Cells: Physiological and Pathological Perspectives. Front Cell Dev Biol. 2020
# 2. Epi3: Kim TH et al., Stomach development, stem cells and disease. Development. 2016
# 3. Yaser Greish et al., Biodegradable Scaffolds for Gastric Tissue Regeneration. 2017
# 4. Epi4: Busslinger GA et al., Human gastrointestinal epithelial of the esophagus, stomach, and duodenum resolved at single-cell resolution. Cell Rep. 2021
# 5. Epi5: Engevik AC et al., The Physiology of the Gastric Parietal Cell. Physiol Rev. 2020
# 6. Epi6: Zhang P et al.,Dissecting the Single-Cell Transcriptome Network Underlying Gastric Premalignant Lesions and Early Gastric Cancer. Cell Rep. 2019
# 7. Epi7: Sathe A et al., Single-Cell Genomic Characterization Reveals the Cellular Reprogramming of the Gastric Tumor Microenvironment. Clin Cancer Res. 2020
# 8. Epi8: Zhang M et al., Dissecting transcriptional heterogeneity in primary gastric adenocarcinoma by single cell RNA sequencing. Gut. 2021
Epithelial_major_marker_ls <- c(
  'Surface mucous cells' = "MUC5AC", # 1
  'Surface mucous cells' = "TFF1", # 1 
  'Surface mucous cells' = "PAS", # 1 
  'Surface mucous cells' = "GKN1", # 9 
  'Surface mucous cells' = "GKN2", # 9 
  
  'Mucous gland cells' = "MUC6", # 1 # 5:TFF2 
  
  'Mucous neck cells' = "TFF2", # 1 
  
  'Parietal cells' = "GIF", # 1 
  'Parietal cells' = "ESRRG", # 2
  'Parietal cells' = "ATP4B", # 2
  
  'Enteroendocrine cells' = "CHGA", # 3,6
  'Enteroendocrine cells' = "CHGB", # 6
  'Enteroendocrine cells' = "GAST", # 9
  'Enteroendocrine cells' = "PROX1", # 9
  
  'Acid-secreting parietal cells' = "ATP4A", # 4 
  #'Acid-secreting parietal cells' = "ATP4B", # 4 
  
  # 'Proliferating isthmus cells' = "MKI67", # 4 
  'Proliferating isthmus cells' = "STMN1", # 4
  
  'Proliferative cells' = "MKI67", # 6
  'Proliferative cells' = "BRIC5", # 6
  'Proliferative cells' = "CDKN2A", # 9
  'Proliferative cells' = "RBP4", # 9
  
  'Chief cells' = "LIPF", # 4
  'Chief cells' = "PGA5", # 4
  'Chief cells' = "MIST1", # 5
  'Chief cells' = "PGC", # 5
  'Chief cells' = "PGA3", # 6
  'Chief cells' = "PGA4", # 6
  
  # 'Endocrine cells' = "GAST", # 5
  'Endocrine cells' = "GHRL", # 5
  'Endocrine cells' = "SST", # 5
  'Endocrine cells' = "HDC", # 5
  
  'Tuft cells' = "TRPM5", # 5
  'Tuft cells' = "SH2D6", # 5
  
  'Tumor cells' = "CEACAM5", # 6
  'Tumor cells' = "CEACAM6", # 6
  'Tumor cells' = "EPCAM", # 9
  'Tumor cells' = "CDH17", # 9
  'Tumor cells' = "MDM2", # P53 pathway
  
  'Enterocytes' = "FABP1", # 6
  'Enterocytes' = "APOA1", # 6
  'Enterocytes' = "FABP2", # 8
  
  'Primary cancer cells' = "NOTCH2", # 7
  'Primary cancer cells' = "NOTCH2NL", # 7
  'Primary cancer cells' = "KIF5B", # 7
  'Primary cancer cells' = "ERBB4", # 7
  
  'Metastatic cancer cells' = "CDK12", # 7
  'Metastatic cancer cells' = "ERBB2", # 7
  'Metastatic cancer cells' = "CLDN11", # 7 
  
  'Driver gene' = "FOS", # 7
  'Driver gene' = "JUN", # 7
  
  'Metaplasia genes' = "CD44", # 8 癌前
  # 'Metaplasia genes' = "TFF2", # 8
  # 'Metaplasia genes' = "SOX9", # 8
  'Metaplasia genes' = "HLA-DPA1", # 8
  'Metaplasia genes' = "HLA-DRA", # 8
  
  'Intestinal metaplasia cells' = "TFF3", # 9
  'Intestinal metaplasia cells' = "CDX1", # 9
  'Intestinal metaplasia cells' = "CDX2", # 9
  
  'Goblet cells' = "MUC2", # 6, 9
  'Goblet cells' = "ITLN1", # 6, 9
  'Goblet cells' = "HES6", # 9
  
  'Metaplastic-stem cells' = "EPHB2", # 6
  'Metaplastic-stem cells' = "SOX9", # 6  
  'Metaplastic-stem cells' = "OLFM4", # 6, 9
  'Metaplastic-stem cells' = "REG1A", # 9
  'Metaplastic-stem cells' = "CLDN3" # 9
)

# find markers
Idents(Epithelial) <- "RNA_snn_res.0.6"
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Epithelial_top5_markers <- 
  Epithelial_markers %>% 
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

names(Epithelial_top5_markers) <- sub('.*\\.','',names(Epithelial_top5_markers))



# dotplot
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_top5_markers)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)


Idents(Epithelial) <- "RNA_snn_res.0.6"
Epithelial <- FindSubCluster(Epithelial, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.6)
Idents(Epithelial) <- "sub.cluster"
#find markers
Epithelial_markers_sub <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Epithelial_top5_markers_sub <- 
  Epithelial_markers_sub %>% 
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

names(Epithelial_top5_markers_sub) <- sub('.*\\.','',names(Epithelial_top5_markers_sub))

DimPlot(Epithelial, label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
DimPlot(Epithelial, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()

dims <- 20
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
DotPlot(Epithelial, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "sub.cluster", features = Epithelial_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "sub.cluster", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "sub.cluster", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "sub.cluster", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "sub.cluster", features = Epithelial_top5_markers_sub)+theme_90angle
dev.off()




## inferCNV ---------------------------------------------------------------------
### set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 50Gb

library(infercnv)
infercnv_out_dir = "Epithelial_subtype/inferCNV_episub" #dir.create(infercnv_out_dir) 


### 映射-------------------------------
library(tibble)
library(data.table)
subtype <- 
  Epithelial[['RNA_snn_res.0.6']] %>% 
  rownames_to_column('cells') %>% 
  mutate(cluster = str_replace(.$RNA_snn_res.0.6, '^', 'cluster_'))

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation5 = ifelse((as.character(annotation) %in% c('Plasma cells', 'Mast cells', 'T cells', 'Monocytes', 'Fibroblasts', 'B cells', 'Pericytes', 'Endothelials', 'Smooth muscle cells')),
                         as.character(annotation), 
                         as.character(RNA_snn_res.0.6)))%>%
  dplyr::select(cells, annotation5) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation5'])
all_data$annotation %>% table()
all_data$annotation5 %>% table()

meta_data_infer <- 
  all_data[[c('annotation5', 'sample')]] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(cluster = str_replace(.$annotation5, '^', 'cluster_')) %>% 
  mutate(annotation6 = ifelse(
    annotation5 %in% c('Plasma cells', 'Mast cells', 'T cells', 'Monocytes', 'Fibroblasts', 'B cells', 'Pericytes', 'Endothelials', 'Smooth muscle cells'),
    annotation5,
    cluster
  )) %>% 
  dplyr::select(cell_id, sample, annotation6)

meta_data_infer_ls <-
  meta_data_infer %>% 
  group_by(sample) %>%
  nest() %>% 
  deframe()


# cell is too much
# less than 50000
set.seed(1234)
meta_data_infer_subset_ls <-
  meta_data_infer_ls %>% 
  imap(function(meta_data_infer, sample){
    meta_data_infer2 <- 
      meta_data_infer %>% 
      group_by(annotation6) %>%
      slice_sample(prop = 1)
  })



# export as file
meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){
    if(!file.exists(file.path(infercnv_out_dir,sample))){
      dir.create(file.path(infercnv_out_dir,sample))
    }
    data.table::fwrite(
      meta_data_infer_subset, 
      file = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  })



infercnv_obj_ls <- 
  meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){ # browser()
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(Epithelial)[1:1000], meta_data_infer_subset$cell_id],
      #raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id],
      annotations_file  = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      gene_order_file   = file.path(infercnv_out_dir, "gene_pos47.txt"),
      delim             = "\t",
      min_max_counts_per_cell = c(30,Inf),
      ref_group_names = c('Plasma cells', 'Mast cells', 'T cells', 'Monocytes', 'Fibroblasts', 'B cells', 'Pericytes', 'Endothelials', 'Smooth muscle cells')
    )
    infercnv_obj
  })



infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. > 0.01) } # 788 <--- 
infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. < 0.01) } # 176

infercnv_obj_ls2 <- 
  infercnv_obj_ls %>% 
  imap(function(infercnv_obj, sample){
    infercnv_obj2 <- infercnv::run(
      infercnv_obj,
      cutoff              =    0.01,
      out_dir    =   file.path(infercnv_out_dir, sample),
      cluster_by_groups   =    TRUE,
      plot_steps          =    FALSE,
      denoise             =    TRUE,
      HMM                 =    FALSE,
      no_prelim_plot      =    TRUE,
      png_res             =    300,
      # useRaster   = FALSE,
      output_format = 'pdf',
      no_plot = T #  figure margin too large 
    )
    infercnv_obj2
  })

infercnv_obj_ls2 %>% 
  imap(function(infercnv_obj2, sample){
    infercnv::plot_cnv(
      infercnv_obj2,
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv',
      output_format     =  "pdf",
      png_res = 300,
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })



infercnv_obj2_default_median_filtering_ls <-
  infercnv_obj_ls2 %>%
  map(function(infercnv_obj2){
    infercnv::apply_median_filtering(infercnv_obj2)
  })


infercnv_obj2_default_median_filtering_ls %>% 
  imap(function(infercnv_obj2_default_median_filtering, sample){
    infercnv::plot_cnv(
      infercnv_obj2_default_median_filtering, 
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv_median_filtering',
      output_format="pdf",
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })

saveRDS(infercnv_obj2_default_median_filtering_ls, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering_ls.rds'))

saveRDS(infercnv_obj_ls2, file.path(infercnv_out_dir, 'infercnv_obj_ls2.rds'))
saveRDS(infercnv_obj_ls, file.path(infercnv_out_dir, 'infercnv_obj_ls.rds'))

### plot -------------------------------------------------------------------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_ls <- 
  cnv_score_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_ls <- 
  temp_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation6,log(CNV_score)))+
      geom_violin(aes(fill= annotation6),color="NA",drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_ls2 <-
  temp_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation6) %>%
      
      ggplot(aes(x = annotation6, y = log(CNV_score), fill = annotation6)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot.pdf"),
      height = 6,
      width  = 8
    )
  })


cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_noQuantile_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_noQuantile_ls <- 
  cnv_score_noQuantile_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_noQuantile_ls <- 
  temp_noQuantile_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation6,log(CNV_score)))+
      geom_violin(aes(fill= annotation6),color="NA",drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_noQuantile_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_noQuantile.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_noQuantile_ls2 <-
  temp_noQuantile_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation6) %>%
      
      ggplot(aes(x = annotation6, y = log(CNV_score), fill = annotation6)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_noQuantile_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot_noQuantile.pdf"),
      height = 6,
      width  = 8
    )
  })




## annotation -----------
cell_cluster <- c(
  '0' = "Intestinal metaplasia cells", 
  '1' = "Surface mucous cells", # MUC5AC, GKN1
  '2' = "Tumor cells",
  '3' = "Tumor cells",
  '4' = "Metastatic cancer cells",
  '5' = "low-nCounts cells", 
  '6' = "myCAFs", 
  '7' = "Chief cells",
  '8' = "Intestinal metaplasia cells", # 
  '9' = "Cd8 texs", # 
  '10' = "Goblet cells", # 
  '11' = "Epi-like cells",
  '12' = "Enteroendocrine cells", # tumor
  '13' = "Chief cells",
  '14' = "Plasma cells", # MZB1, XBP1
  '15' = "low-nCounts cells", # tumor
  '16' = "Primary cancer cells",
  '17' = "Surface mucous cells",
  '18' = "low-nCounts cells", 
  '19' = "low-nCounts cells"
)


Idents(Epithelial) <- "RNA_snn_res.0.6"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, 
                                 levels=c(Epithelial_meta$sample %>% unique())
)
p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Epithelial_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
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
  filename = file.path("Epithelial_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  'Intestinal metaplasia cells' = "CDX1", # 9
  'Intestinal metaplasia cells' = "CDX2", # 9
  
  'Surface mucous cells' = "MUC5AC", # 1
  'Surface mucous cells' = "GKN1", # 9 
  
  'Proliferating genes' = "STMN1", # 4
  'Proliferating genes' = "MKI67", # 6
  
  'Metastatic cancer cells' = "ERBB2", # 7
  
  'Tumor cells' = "EPCAM", # 9
  'Tumor cells' = "MDM2", # P53 pathway
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  
  'myCAF' = 'MMP11',
  'myCAF' = 'CTHRC1',
  'myCAF' = 'COL11A1',
  'myCAF' = 'POSTN',
  'myCAF' = 'COL10A1',
  'myCAF' = 'TNC',
  
  'Chief cells' = "LIPF", # 4
  'Chief cells' = "PGA3", # 6
  
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  
  'Goblet cells' = "MUC2", # 6, 9
  'Goblet cells' = "ITLN1", # 6, 9
  'Goblet cells' = "HES6", # 9
  
  'Enteroendocrine cells' = "PROX1", # 9
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'Primary cancer cells' = "NOTCH2", # 7
  'Primary cancer cells' = "KIF5B", # 7
  
  'Metaplasia genes' = "CD44", # 8 癌前
  
  'Metaplastic-stem cells' = "EPHB2", # 6
  'Metaplastic-stem cells' = "SOX9" # 6  
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()




## volcano-------
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_markers,'Epithelial_subtype/Epithelial_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(Epithelial_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_markers[Epithelial_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('Epithelial_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

Epithelial_markers_ls <- list(Epithelial_markers)
names(Epithelial_markers_ls) <- 'Epithelial cells'

# scatter
p_volcanoplot_major_ls <- 
  Epithelial_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('Epithelial_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })




## enrichment ----------------------------
Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


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


er_output<-'Epithelial_subtype/enrichment' # dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
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
      filename = file.path(er_output, name, "go dotplot.pdf"),
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
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
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
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(Epithelial, 'Epithelial_subtype/Epithelial.rds')




## 映射-------------------------------
library(tibble)
library(data.table)
sub <- Epithelial[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation3 = ifelse((as.character(annotation2) %in% c('Epithelials', 'Foveolar cells')),
                         as.character(subtype), 
                         as.character(annotation2)))%>%
  dplyr::select(cells, annotation3) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data,meta.data.t['annotation3'])
saveRDS(all_data, 'all_data_annotation.rds')




# spatial plot ------------------------------------------------------------------
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
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y, color = annotation)) + 
      scattermore::geom_scattermore(
        pointsize = 1 ,
        alpha = 0.1,
        # show.legend = c(color = TRUE, alpha = FALSE), 
        pixels = c(800, 800)
      ) +
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_color_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_cellType.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y, color = annotation3)) + 
      scattermore::geom_scattermore(
        pointsize = 1 ,
        alpha = 0.1,
        # show.legend = c(color = TRUE, alpha = FALSE), 
        pixels = c(800, 800)
      ) +
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_color_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_subType.pdf")), 
        height = 10, 
        width = 15)
    print(p)
    dev.off()
  })  


all_data$annotation3 %>% table()
all_data$annotation4 <- fct_collapse(
  all_data$annotation3, 
  'Others' = unique(all_data$annotation3)[!unique(all_data$annotation3) %in% c('Intestinal metaplasia cells', 'Metastatic cancer cells', 'Tumor cells', 'Primary cancer cells')]
)
all_data$annotation4 %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y, color = annotation4)) + 
      scattermore::geom_scattermore(
        pointsize = 1 ,
        alpha = 0.1,
        # show.legend = c(color = TRUE, alpha = FALSE), 
        pixels = c(800, 800)
      ) +
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_color_manual(values = c('lightgrey',"#DC143C","#0000FF","#20B2AA","#FFA500")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_tumor highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  

# Annotation info output -------------------------------------------------------
all_data$annotation %>% table()
cellType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

cellType_anno_ls %>% 
  iwalk(function(cellType, project){
    write.csv(cellType, str_glue('anno_info/cell type {project}.csv'), row.names = F)
  })


all_data$annotation3 %>% table()
subType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation3')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation3
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

subType_anno_ls %>% 
  iwalk(function(subType, project){
    write.csv(subType, str_glue('anno_info/sub cell type {project}.csv'), row.names = F)
  })

all_data$annotation4 %>% table()
tumor_anno_ls <- 
  all_data@meta.data[c('project', 'annotation4')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation4
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

tumor_anno_ls %>% 
  iwalk(function(tumorType, project){
    write.csv(tumorType, str_glue('anno_info/tumor highlight {project}.csv'), row.names = F)
  })




# Epi reannotation --------------------------------------------------------------
Epithelial_meta <- Epithelial[[]] %>% rownames_to_column('cell_id')
p_dimPlot <- 
  cells_boundary_ls[2] %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(Epithelial_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y, color = seurat_clusters)) + 
      scattermore::geom_scattermore(
        pointsize = 1 ,
        alpha = 0.1,
        # show.legend = c(color = TRUE, alpha = FALSE), 
        pixels = c(800, 800)
      ) +
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_color_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })


Epithelial$subtype %>% table()
epi_sub_anno_ls <- 
  Epithelial@meta.data[c('project', 'subtype')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = subtype
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

epi_sub_anno_ls %>% 
  iwalk(function(epi_subType, project){
    write.csv(epi_subType, str_glue('anno_info/Epithelial subtype {project}.csv'), row.names = F)
  })

Epithelial$RNA_snn_res.0.6 %>% table()
epi_clusters_anno_ls <- 
  Epithelial@meta.data[c('project', 'RNA_snn_res.0.6')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = RNA_snn_res.0.6
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

epi_clusters_anno_ls %>% 
  iwalk(function(epi_seurat_clusters, project){
    write.csv(epi_seurat_clusters, str_glue('anno_info/Epithelial seurat_clusters {project}.csv'), row.names = F)
  })



## annotation -------------------------------------------------------------------
cell_cluster <- c(
  '0' = "Tumor cells", 
  '1' = "Surface mucous cells", # MUC5AC, GKN1
  '2' = "Tumor cells",
  '3' = "Tumor cells",
  '4' = "Tumor cells",
  '5' = "low-nCounts cells", 
  '6' = "myCAFs", 
  '7' = "Chief cells",
  '8' = "Tumor cells", # 
  '9' = "Cd8 texs", # 
  '10' = "Goblet cells", # 
  '11' = "Epi-like cells",
  '12' = "Tumor cells", # tumor
  '13' = "Chief cells",
  '14' = "Plasma cells", # MZB1, XBP1
  '15' = "low-nCounts cells", # tumor
  '16' = "Tumor cells",
  '17' = "Surface mucous cells",
  '18' = "low-nCounts cells", 
  '19' = "low-nCounts cells"
)


Idents(Epithelial) <- "RNA_snn_res.0.6"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation_2.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split_2.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
Epithelial$clinical_type <- str_replace_all(Epithelial$sample, c(
  'P1$' = 'T1', 
  'P2' = 'T1', 
  'P3' = 'T1', 
  'P4' = 'T2', 
  'P5' = 'T2', 
  'P6' = 'T3', 
  'P7' = 'T3', 
  'P8' = 'T4', 
  'P9' = 'T4',
  'P10' = 'T4'
))

Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, 
                                 levels=c(Epithelial_meta$sample %>% unique())
)
p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Epithelial_subtype/cell_freq_sample_2.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
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
  filename = file.path("Epithelial_subtype/cell_freq_group_2.pdf"),
  height = 6,
  width  = 4
)


p_col_proj <-
  Epithelial_meta %>%
  group_by(clinical_type, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = clinical_type, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path("Epithelial_subtype/cell_freq_clinical_type_2.pdf"),
  height = 6,
  width  = 4
)


cellType_marker_ls <- c(
  'Tumor cells' = "ERBB2", # 7
  
  'Tumor cells' = "EPCAM", # 9
  'Tumor cells' = "MDM2", # P53 pathway
  
  'Proliferating genes' = "STMN1", # 4
  'Proliferating genes' = "MKI67", # 6
  
  'Surface mucous cells' = "MUC5AC", # 1
  'Surface mucous cells' = "GKN1", # 9 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  'myCAF' = 'MMP11',
  'myCAF' = 'CTHRC1',
  'myCAF' = 'COL11A1',
  'myCAF' = 'POSTN',
  'myCAF' = 'COL10A1',
  'myCAF' = 'TNC',  
  
  'Chief cells' = "LIPF", # 4
  'Chief cells' = "PGA3", # 6
  
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  
  'Goblet cells' = "MUC2", # 6, 9
  'Goblet cells' = "ITLN1", # 6, 9
  'Goblet cells' = "HES6", # 9
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1"
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker_2.pdf"), width = 20, height = 8)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()




## volcano-------
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_markers,'Epithelial_subtype/Epithelial_diff_gene_2.csv')

library(EnhancedVolcano)
type<-unique(Epithelial_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_markers[Epithelial_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('Epithelial_subtype', str_glue("{name} diffgene valcano_2.pdf")),
      height = 6,
      width = 8
    )
  })

Epithelial_markers_ls <- list(Epithelial_markers)
names(Epithelial_markers_ls) <- 'Epithelial cells'

# scatter
p_volcanoplot_major_ls <- 
  Epithelial_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('Epithelial_subtype', str_glue("{name} diffgene dot_2.pdf")),
      height = 10,
      width = 20
    )
  })




## enrichment ----------------------------
Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


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


er_output<-'Epithelial_subtype/enrichment' # dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table_2.csv"),
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
      filename = file.path(er_output, name, "go dotplot_2.pdf"),
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
      filename = file.path(er_output, name, "go barplot_2.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
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
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table_2.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot_2.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot_2.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(Epithelial, 'Epithelial_subtype/Epithelial2.rds')




## 映射-------------------------------
library(tibble)
library(data.table)
sub <- Epithelial[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation6 = ifelse((as.character(annotation2) %in% c('Epithelials', 'Foveolar cells')),
                         as.character(subtype), 
                         as.character(annotation2)))%>%
  dplyr::select(cells, annotation6) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data,meta.data.t['annotation6'])
saveRDS(all_data, 'all_data_annotation.rds')