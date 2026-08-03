setwd("/mnt/NexData_file/home/lxy/case/XE2601050026-01/stRNA_3")

all_data <- readRDS('output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")

# normalize
all_data <- all_data %>% NormalizeData()

saveRDS(all_data ,file = "output/data/all_data.rds")


library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
  "Adipocyte" = "PLIN1",
  "Adipocyte" = "PLIN4",
  "Adipocyte" = "ADIPOQ",
  "Adipocyte" = "FABP4",
  "Adipocyte" = "LPL",
          
  "cTEC" = "LY75", 
  "cTEC" = "CCL25", 
  "cTEC" = "CXCL12", 
  
  "mTEC" = "EPCAM", 
  "mTEC" = "CCL19", 
  "mTEC" = "CCL21", 
  
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
  "T cell" = "CD4", 
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  "Proliferating T cell" = "MKI67", 
  "Proliferating T cell" = "TOP2A", 
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  "Pericyte" = "PLVAP", 
  
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
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
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
  
  "Mesenchyme" = "PDGFRA", 
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )




# T_cells -------------------------------------------------------------------
all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')

T_cell_out_dir <- "T_cell_subtype" 
dir.create(T_cell_out_dir)
T_cell <- subset(all_data, annotation == 'T cell')

## PCA --------
T_cell <- FindVariableFeatures(T_cell)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

pdf(file.path(T_cell_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(T_cell, dims = 1:5)
dev.off()
pdf(file.path(T_cell_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(T_cell, dims = 6:10)
dev.off()
pdf(file.path(T_cell_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(T_cell, dims = 11:15)
dev.off()
pdf(file.path(T_cell_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(T_cell, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 15 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

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

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

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
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(T_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(T_cell) <- "RNA_snn_res.0.6"

### dotplot -------
T_major_marker_ls <- c(
  "CD8"="CD8A",
  "CD8"="CD8B",
  
  "CD4"="CD4",
  "CD4"="TRDC",
  "CD4"="IL2",
  
  # "Naive"="CCR6",
  "Naive"="LEF1",
  "Naive"="SELL",
  # "Naive"="CCR7",
  "Naive"="YBX3",
  "Naive"="TGFB1",
  # "Naive"="CCL5",
  "Naive"="CXCR3", 
  
  # "Memory"="TCF7",
  "Memory"="GPR183",
  "Memory"="IL7R",
  "Memory"="FOS",
  "Memory"="JUN",
  "Memory"="KLRB1", 
  
  "TRM"="CD89",
  "TRM"="CXCR8", 
  
  # "TEM"="IFNG",
  "TEM"="OCL4",
  "TEM"="OCL5",
  "TEM"="JAK3", 
  "TEM"="GZMK", 
  "TEM"="GZMA", 
  "TEM"="NKG7", 
  
  "TEMRA"="KLRG1",  
  
  "TCM"="TCF7", 
  "TCM"="IDO1",
  "TCM"="SKIV2L",
  "TCM"="TRAF2",
  "TCM"="SLC2A1",
  "TCM"="CD40LG",    
  "TCM"="DUSP2",
  "TCM"="LYAR",
  "TCM"="IL32", 
  "TCM"="CCR7", 
  
  "TRM"="ZNF683", 
  "TRM"="CXCR6", 
  
  "Effector"="FGFBP2",
  "Effector"="TYROBP",
  "Effector"="KLRC2",
  "Effector"="CX3CR1",
  "Effector"="FCGR3A",
  "Effector"="GNLY", 
  
  "cytotoxic"="XCL1",
  "cytotoxic"="IFNG",
  "cytotoxic"="CCL4",
  "cytotoxic"="CCL5", 
  
  "Exhausted"="LAYN",
  "Exhausted"="LAG3",
  "Exhausted"="HAVCR2",
  "Exhausted"="TNFRSF9",
  "Exhausted"="CTLA4",
  "Exhausted"="TIGIT",
  "Exhausted"="TOX",
  #"Exhausted"="PDCD1",
  
  "Tfh"="CD200",
  "Tfh"="CD83",
  "Tfh"="CXCL3",
  "Tfh"="ICOS",
  "Tfh"="BATF",
  "Tfh"="PDCD1", 
  "Tfh"="CXCR5",
  
  "proliferation"="MKI67",
  "proliferation"="H4C3",
  "proliferation"="TUBB",
  "proliferation"="STMN1",
  "proliferation"="TOP1",
  "proliferation"="ANAPC11",
  "proliferation"="MCM7",
  "proliferation"="DNMT1",
  "proliferation"="TOP2A", 
  
  "CD4_Treg"="TNFRSF4",
  # "CD4_Treg"="TNFRSF9",
  "CD4_Treg"="FOXP3",
  "CD4_Treg"="IL2RA",
  # "CD4_Treg"="LAYN",
  # "CD4_Treg"="CTLA4",
  # "CD4_Treg"="TIGIT", 
  "CD4_Treg"="CD25",
  "CD4_Treg"="STAT1",
  
  # "CD4_naive"="LEF1", 
  
  # "CD4_Tem"="GZMK",
  # "CD4_Tem"="GZMA",
  # "CD4_Tem"="NKG7",
  # "CD4_Tem"="CCL4",
  # "CD4_Tem"="CCL5",
  "CD4_Tem"="CAPG", 
  
  "CD4_Th1"="ANXA1",
  "CD4_Th1"="CXCL13",
  "CD4_Th1"="BHLHE40",
  # "CD4_Th1"="IFNG", 
  
  "CD4_Tfh"="TOX2",
  # "CD4_Tfh"="PDCD1",
  # "CD4_Tfh"="CXCR5", 
  
  "CD4_Th17"="IL17A",
  "CD4_Th17"="IL17F",
  "CD4_Th17"="IL23R",
  "CD4_Th17"="RORC",
  "CD4_Th17"="CCR6", 
  
  # "CD4_Tcm"="TCF7",
  # "CD4_Tcm"="CCR7",
  "CD4_Tcm"="S1PR1",
  "CD4_Tcm"="TXK",
  "CD4_Tcm"="ICAM2",
  "CD4_Tcm"="CDC25B",
  "CD4_Tcm"="PTGER2", 
  
  # "CD8_Tex"="HAVCR2",
  # "CD8_Tex"="CTLA4",
  # "CD8_Tex"="LAYN",
  # "CD8_Tex"="LAG3",
  # "CD8_Tex"="TIGIT",
  "CD8_Tex"="ENPD1", 
  # "CD8_Tex"="TOX",
  # "CD8_Tex"="PDCD1", 
  
  # "CD8_Tem"="GZMK",
  # "CD8_Tem"="GZMA",
  # "CD8_Tem"="IFNG",
  # "CD8_Tem"="KLRG1",
  "CD8_Tem"="GZMB",
  "CD8_Tem"="PRF1",
  # "CD8_Tem"="NKG7", 
  
  # "CD8_Tnaive"="SELL",
  # "CD8_Tnaive"="LEF1",
  # "CD8_Tnaive"="CCR7", 
  
  "CD8_MAIT"="NCR3",
  "CD8_MAIT"="SLC4A10",
  # "CD8_MAIT"="KLRB1", 
  
  "CD8_Trm"="MYADM",
  "CD8_Trm"="CD69",
  "CD8_Trm"="ITGAE",
  # "CD8_Trm"="CXCR6",
  "CD8_Trm"="ITGA1", 
  
  # "CD8_Tcm"="IL7R",
  "CD8_Tcm"="CD28" 
  # "CD8_Tcm"="GPR183"
)


Idents(T_cell)<-'RNA_snn_res.0.6'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

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

# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_top5_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


###score--------------------
# 打分这里可做可不做，如果要研究巨噬细胞的亚型，需要打分
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
  filename = 'Macrophage_subtype/M1_M2 score.pdf',
  height = 10,
  width  = 15
) 



## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "vCAFs", # POSTN, PDGFRA, COL18A1
  '1'  = "CD4", #
  '2'  = "CD8",
  '3'  = "Epithelials", # https://zhuanlan.zhihu.com/p/549337063
  '4'  = "TEM", # PGA3
  '5'  = "TCM", # MMP9,10.1007/s11033-022-07184-9, MMP12,https://doi.org/10.1016/j.jhazmat.2023.132733
  '6'  = "cDC2s", # 
  '7'  = "Proliferating T_cells",
  '8'  = "LAMP3 DCs",
  '9'  = "Thf", # 10.1002/jcsm.12584
  '10' = "Epithelials", # 
  '11' = "cDC1s"
)

Idents(T_cell) <- "RNA_snn_res.0.6"

T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# 把这个髓系亚群替换成免疫细胞的，比如T和B细胞
T_cellTye_markers <- c( 
  "T_cell" = "CD68",
  # "T_cell" = "CD163",
  "T_cell" = "C1QA",
  "T_cell" = "C1QC", 
  # 'T_cell'  = 'CSF1R',
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "POSTN",
  # "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  'vCAFs' = "MCAM", 
  'vCAFs' = "NOTCH3", 
  'vCAFs' = "COL18A1",
  
  'M2 T_cells' = "CD163",
  'M2 T_cells' = "CSF1R",
  'M2 T_cells' = "MRC1",
  'M2 T_cells' = "CCL24",
  'M2 T_cells' = "VEGFB",
  'M2 T_cells' = "CTSA",
  'M2 T_cells' = "CD209",
  'M2 T_cells' = "FOLR2",
  
  'Epithelials' = "EPCAM",
  
  'Epithelial-like T_cells' = 'PGA3',
  
  'cDC2s' = 'CD1C',
  'cDC2s' = 'CLEC10A',
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  'LAMP3 DCs' = "LAMP3",
  'LAMP3 DCs' = 'CD83', 
  
  'IFN-responsive T_cells' = "IFIT3", 
  'IFN-responsive T_cells' = "IFIT2",
  'IFN-responsive T_cells' = "IFIT1",
  
  'cDC1s' = 'CLEC9A',
  'cDC1s' = 'IDO1'
)

pdf(file.path(T_cell_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(T_cell, group.by = "subtype", features = T_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()






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
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("T_cell_subtype", "cell_freq_sample.pdf"),
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


## volcano -------------------
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(T_cell_markers,'T_cell_subtype/T_cell_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(T_cell_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- T_cell_markers[T_cell_markers$cluster == type, ]
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
      filename = file.path(T_cell_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


T_cell_markers_major<-list(T_cell = T_cell_markers)
# scatter
p_volcanoplot_major_ls <- 
  T_cell_markers_major %>%
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
      filename = file.path(T_cell_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 25
    )
  })


## enrichment ---------------------------------------
T_cell_de_genes_for_enrichmenet <- 
  T_cell_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

T_cell_er_output<-'T_cell_subtype/enrichment' # dir.create(T_cell_er_output)
enrich_ls <-
  T_cell_de_genes_for_enrichmenet %>%
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
    dir.create(file.path(T_cell_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(T_cell_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(T_cell_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(T_cell_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


T_cell_de_genes_id_for_enrichmenet <-
  T_cell_de_genes_for_enrichmenet %>%
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
  T_cell_de_genes_id_for_enrichmenet %>%
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

# saveRDS(T_cell_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(T_cell_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(T_cell_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(T_cell_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(T_cell_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(T_cell, 'T_cell_subtype/all_T_cell.rds')

