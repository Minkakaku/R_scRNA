library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/output/data/all_data.rds')


# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
  "Adipocyte" = "PLIN1",
  "Adipocyte" = "PLIN4",
  "Adipocyte" = "ADIPOQ",
  "Adipocyte" = "FABP4",
  "Adipocyte" = "LPL",
  
  "Endothelial" = "PECAM1",
  "Endothelial" = "VWF", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2", 
  
  "Epithelial"= "EPCAM",
  
  "Mesangial cell" = "PDGFRB",
  "Mesangial cell" = "CFH",
  "Mesangial cell" = "ZEB2",
  "Mesangial cell" = "COL12A1",
  "Mesangial cell" = "MGP",
  
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
  
  #"Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  #"Pericyte" = "PLVAP", 
  
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
  #"Dendritic cell" = "CLEC9A", 
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
  
  "Neutrophil" = "CSF3R", 
  "Neutrophil" = "FCGR3B", 
  "Neutrophil" = "S100A8", 
  "Neutrophil" = "S100A9", 
  "Neutrophil" = "G0S2", 
  "Neutrophil" = "RETNLG", 
  
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
  
  #"Mesenchyme" = "PDGFRA", 
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2",
  
  # 肾小球
  "podocyte" = "CRYAB",
  "podocyte" = "PDPN",
  "podocyte" = "SYNPO",
  "podocyte" = "ROBO2",
  "podocyte" = "NPHS1",
  "podocyte" = "NPHS2",
  "podocyte" = "TCF21",
  "podocyte" = "PODXL",
  "podocyte" = "VEGFA",
  "podocyte" = "WT1",
  
  # 近端小管
  #"proximal tubule" = "AQP1", # 重复基因：在 Endothelial 出现
  "PT" = "LRP2",
  "PT" = "SLC13A1",
  "PT" = "SLC27A2",
  
  #"proximal convoluted tubule (PCT)" = "SLC34A1",
  #"proximal convoluted tubule (PCT)" = "SLC5A2", # 重复基因
  
  "PST" = "KAP",
  #"proximal straight tubule (PST)" = "SLC22A30", # 重复基因
  
  "PTS1" = "SLC5A2", # 重复基因
  "PTS1" = "SLC16A1",
  "PTS1" = "SLC34A3",
  "PTS1" = "SPP2",
  
  "PTS2" = "CNDP2",
  "PTS2" = "NUDT19",
  "PTS2" = "SLC22A12",
  "PTS2" = "SLC22A6",
  
  "PTS3" = "SLC22A30", # 重复基因
  "PTS3" = "SLC7A13",
  "PTS3" = "SLC6A18",
  "PTS3" = "SLC5A10",
  
  # 远端小管 & 亨利袢
  "DCT" = "CALB1",
  "DCT" = "SLC12A3",
  "DCT" = "PVALB",
  "DCT" = "SLC8A1",
  "DCT" = "TRPM6",
  
  "LOH" = "UMOD",
  "LOH" = "SLC12A1",
  
  # 集合管
  #"CD" = "SLC26A4", # 重复基因
  #"CD" = "INSRR", # 重复基因
  #"CD" = "AQP2", # 重复基因
  #"CD" = "ATP6V1G3", # 重复基因
  #"CD" = "ATP6V1B1", # 重复基因
  "CD" = "TRPV5",
  
  "CD-IC" = "SLC26A4", # 重复基因
  "CD-IC" = "SCNN1A",
  "CD-IC" = "ATP6V0D2",
  "CD-IC" = "ATP6V1G3", # 重复基因
  "CD-IC" = "ATP6V1B1", # 重复基因
  "CD-IC" = "SLC4A1",
  
  "CD-PC" = "HSD11B2",
  "CD-PC" = "AQP2", # 重复基因
  "CD-PC" = "FXYD4",
  "CD-PC" = "STC1",
  "CD-PC" = "AVPR2",
  
  "CD-Trans" = "RHBG",
  "CD-Trans" = "INSRR", # 重复基因
  
  "MCD" = "UPK1B",
  "MCD" = "UPK3A",
  
  # 肾单位祖细胞
  "NP" = "UNCX",
  "NP" = "CITED1"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )




# Fibroblasts -------------------------------------------------------------------
#all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')
all_data$annotation %>% table()
Fibroblast_out_dir <- "Fibroblast_subtype" 
dir.create(Fibroblast_out_dir)
Fibroblast <- subset(all_data, annotation == 'Fibroblast')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Fibroblast <- FindVariableFeatures(Fibroblast)

Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)

### pca inspector ------
pdf(file.path(Fibroblast_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 16 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
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
pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Fibroblast) <- "RNA_snn_res.0.6"

### dotplot -------
Fibro_major_marker_ls <- c(
  # 静息成纤维细胞 (FIB)
  "FIB" = "PDGFRA",
  "FIB" = "CDH11",
  "FIB" = "NEGR1",
  "FIB" = "C7",
  "FIB" = "FBLN5",
  "FIB" = "LAMA2",
  
  # 分化/发育相关成纤维细胞 (dFIB)
  #"dFIB" = "MEG3",
  "dFIB" = "DCN",
  "dFIB" = "B2M",
  "dFIB" = "COL1A2",
  "dFIB" = "VIM",
  "dFIB" = "S100A6",
  "dFIB" = "AEBP1",
  "dFIB" = "COL1A1",
  "dFIB" = "SFRP1",
  "dFIB" = "CXCL12",
  #"dFIB" = "PDGFRA",
  "dFIB" = "COL14A1",
  #"dFIB" = "POSTN",
  "dFIB" = "IGFBP7",
  "dFIB" = "TPT1",
  
  # 活化成纤维细胞 (aFIB)
  #"aFIB" = "POSTN",
  "aFIB" = "FLRT2",
  "aFIB" = "FGF14",
  "aFIB" = "PRRX1",
  "aFIB" = "NAV3",
  "aFIB" = "IGF1",
  "aFIB" = "KRT19",
  "aFIB" = "ABI3BP",
  "aFIB" = "TWIST2",
  
  # 肌成纤维细胞 (MYOF)
  #"MYOF" = "TWIST2",
  "MYOF" = "FAP",
  "MYOF" = "SULF1",
  "MYOF" = "ACTA2",
  "MYOF" = "POSTN",
  "MYOF" = "COL6A3",
  "MYOF" = "COL16A1",
  "MYOF" = "GLI2",
  "MYOF" = "COL5A1",
  "MYOF" = "NTM",
  "MYOF" = "NOTCH3",
  
  # 循环/增殖型肌成纤维细胞 (cycMYOF)
  "cycMYOF" = "MKI67",
  "cycMYOF" = "TOP2A"
)

Fibro_major_marker_ls_1 <- c(
  # 肌成纤维细胞 MYOF
  "MYOF" = "DMD",
  "MYOF" = "COL6A3",
  "MYOF" = "LMOD1",
  "MYOF" = "COL1A1",
  "MYOF" = "SULF1",
  "MYOF" = "SORBS1",
  "MYOF" = "SYNPO2",
  "MYOF" = "PALLD",
  "MYOF" = "PRRX1",
  "MYOF" = "ECRG4",
  "MYOF" = "ADRA2A",
  "MYOF" = "ADGRL3",
  "MYOF" = "INPP4B",
  "MYOF" = "HSPB6",
  "MYOF" = "STEAP4",
  "MYOF" = "MUSTN1",
  "MYOF" = "VGLL3",
  "MYOF" = "AOC3",
  "MYOF" = "ADAM33",
  "MYOF" = "LHFPL6",
  "MYOF" = "KCNMA1",
  "MYOF" = "NTM",
  "MYOF" = "SLC24A3",
  "MYOF" = "RXFP1",
  "MYOF" = "SMOC2",
  
  # 增殖型肌成纤维细胞 cycMYOF
  "cycMYOF" = "GLI2",
  "cycMYOF" = "COL5A1",
  "cycMYOF" = "ADAMTS12",
  "cycMYOF" = "FN1",
  "cycMYOF" = "TNC",
  "cycMYOF" = "MELK",
  "cycMYOF" = "BRIP1",
  
  # 静息成纤维细胞 FIB
  "FIB" = "NEGR1",
  "FIB" = "C7",
  "FIB" = "LPAR1",
  "FIB" = "SVEP1",
  "FIB" = "GNA14",
  "FIB" = "NFASC",
  "FIB" = "ITGA1",
  "FIB" = "ADH1B",
  "FIB" = "NR2F1",
  "FIB" = "ZEB2",
  "FIB" = "RPS23",
  "FIB" = "FBLN5",
  "FIB" = "LUM",
  
  # 基质型肌成纤维细胞 M-FIB
  "M-FIB" = "SYT1",
  "M-FIB" = "GABRG3",
  "M-FIB" = "LSAMP",
  "M-FIB" = "CRISPLD2",
  "M-FIB" = "CLMP",
  "M-FIB" = "NID1",
  "M-FIB" = "MOXD1",
  "M-FIB" = "DACH2",
  "M-FIB" = "PLCXD3",
  "M-FIB" = "HPSE2",
  "M-FIB" = "ADAMTSL1",
  "M-FIB" = "MGP",
  "M-FIB" = "BGN",
  
  # 分化型肌成纤维细胞 dM-FIB
  "dM-FIB" = "AKR1B1",
  "dM-FIB" = "SLPI",
  "dM-FIB" = "LTBP1",
  
  # 活化成纤维细胞 aFIB
  "aFIB" = "DCN",
  "aFIB" = "PDZRN4",
  "aFIB" = "CCDC80",
  "aFIB" = "NAV3",
  "aFIB" = "FGF14",
  "aFIB" = "ESR1",
  "aFIB" = "FBLN1",
  "aFIB" = "FLRT2",
  "aFIB" = "MFAP5",
  "aFIB" = "SFRP4",
  "aFIB" = "C3",
  "aFIB" = "PLA2G2A",
  "aFIB" = "OGN",
  "aFIB" = "PROCR",
  "aFIB" = "SFRP2",
  "aFIB" = "SERPINF1",
  "aFIB" = "PCOLCE",
  
  # 分化型成纤维细胞 dFIB
  #"dFIB" = "MGP",       # 重复 → 已注释
  #"dFIB" = "BGN",       # 重复 → 已注释
  #"dFIB" = "DCN",       # 重复 → 已注释
  "dFIB" = "AEBP1",
  "dFIB" = "COL1A2",
  #"dFIB" = "LUM",       # 重复 → 已注释
  #"dFIB" = "COL1A1",    # 重复 → 已注释
  "dFIB" = "C1R",
  "dFIB" = "TIMP1",
  "dFIB" = "COL6A2"
)

Fibro_major_marker_ls_2 <- c(
  "Pericyte"       = "NOTCH3",
  "Pericyte"       = "RGS5",
  "Pericyte"       = "PDGFRB",
  
  "Fibroblast"   = "PDGFRA", 
  "Fibroblast_1"   = "MEG3",
  "Fibroblast_1"   = "SCARA5",
  
  "Fibroblast_2"   = "CXCL12",
  "Fibroblast_2"   = "COLEC11",
  
  "Myofibroblast_1"= "POSTN",
  "Myofibroblast_1"= "NRP3",
  "Myofibroblast_1"= "NKD2",
  "Myofibroblast_2"= "OGN",
  "Myofibroblast_2"= "COL14A1",
  "Myofibroblast_3"= "CCL19",
  "Myofibroblast_3"= "CCL21",
  "Myofibroblast_1a"= "GREM2",
  "Myofibroblast_1a"= "FRZB"
)


Idents(Fibroblast)<-'RNA_snn_res.0.6'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Fibroblast_markers, file = "Fibroblast_subtype/Fibroblast_all_markers.csv", row.names = TRUE)

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
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

names(Fibroblast_top5_markers) <- sub('.*\\.','',names(Fibroblast_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Fibroblast_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_2)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_markers)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Fibroblast ,file = "Fibroblast_subtype/Fibroblast_cluster.rds")
Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "CXCL12+Fibroblast", # 
  '1'  = "MEG3+Fibroblsat", #
  '2'  = "Macrophage-like_Fibroblsat", # 
  '3'  = "Pericyte", # CD8+TEM
  '4'  = "LOH-like_Fibroblsat", # 
  '5'  = "PT-like_Fibroblsat",#  
  '6'  = "Endothelial-like_Fibroblsat", #  
  '7'  = "Myofibroblast", #  
  '8'  = "Schwann cell"
)

Idents(Fibroblast) <- "RNA_snn_res.0.6"




Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "CXCL12+Fibroblast",
  "MEG3+Fibroblsat",
  "Myofibroblast"
)

Fibroblast_clean <- subset(Fibroblast, subtype %in% keep)

Fibroblast_clean$subtype <- droplevels(Fibroblast_clean$subtype)
Idents(Fibroblast_clean) <- Fibroblast_clean$subtype

table(Fibroblast_clean$subtype)
saveRDS(Fibroblast_clean, "Fibroblast_subtype/Fibroblast_clean.rds")

Fibroblast <- Fibroblast_clean

Fibroblast <- RunUMAP(Fibroblast, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} Fibroblast annotation split.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()




# 保留：0,1,7
Fibroblast <- subset(Fibroblast, RNA_snn_res.0.6 %in% c(0,1,7))
Fibroblast$RNA_snn_res.0.6 %>% table()


## PCA --------
Fibroblast <- FindVariableFeatures(Fibroblast)

Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)

### pca inspector ------
pdf(file.path(Fibroblast_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 16 # 18

library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)


Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



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
pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} cluster_1.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(Fibroblast)<-'RNA_snn_res.0.2'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Fibroblast_markers, file = "Fibroblast_subtype/Fibroblast_all_markers_1.csv", row.names = TRUE)

Fibroblast_top5_markers <- 
  Fibroblast_markers %>% 
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

names(Fibroblast_top5_markers) <- sub('.*\\.','',names(Fibroblast_top5_markers))

#Idents(Fibroblast) <- "RNA_snn_res.0.4"
pdf(file.path(Fibroblast_out_dir, str_glue("dotplot dims{dims} major_marker_1.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.2", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.2", features = Fibro_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_markers)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.2", features = Fibroblast_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("dotplot dims{dims} major_marker_1_0.6.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls_3)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_markers)+RotatedAxis()
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibroblast_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

saveRDS(Fibroblast ,file = "Fibroblast_subtype/Fibroblast_cluster_1.rds")
Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "CXCL12+Fibroblast", # 
  '1'  = "MEG3+Fibroblsat", #
  '2'  = "Myofibroblast"
)


Idents(Fibroblast) <- "RNA_snn_res.0.2"




Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

table(Fibroblast$subtype)
saveRDS(Fibroblast, "Fibroblast_subtype/Fibroblast_clean.rds")

#keep <- c(
#  "CD4+",
#  "CD8+TEM",
#  "CD8+Effector",
#  "Treg",
#  "Proliferating T cell"
#)

#Fibroblast_clean <- subset(Fibroblast, subtype %in% keep)

#Fibroblast_clean$subtype <- droplevels(Fibroblast_clean$subtype)
#Idents(Fibroblast_clean) <- Fibroblast_clean$subtype

#table(Fibroblast_clean$subtype)
#saveRDS(Fibroblast_clean, "Fibroblast_subtype/Fibroblast_clean.rds")

#Fibroblast <- Fibroblast_clean



#Fibroblast <- RunUMAP(Fibroblast, dims = 1:16,n.neighbors = 50,return.model = TRUE)

#pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} T annotation.pdf")), height = 6, width = 8)
#DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
#dev.off()

#pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} T annotation split.pdf")), height = 6, width = 15)
#DimPlot(Fibroblast, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
#dev.off()




Fibroblast_1 <- readRDS("Fibroblast_subtype/Fibroblast_cluster.rds")
cell_cluster <- c(
  '0'  = "CXCL12+Fibroblast", # 
  '1'  = "MEG3+Fibroblsat", #
  '2'  = "Macrophage-like_Fibroblsat", # 
  '3'  = "Pericyte", # CD8+TEM
  '4'  = "LOH-like_Fibroblsat", # 
  '5'  = "PT-like_Fibroblsat",#  
  '6'  = "Endothelial-like_Fibroblsat", #  
  '7'  = "Myofibroblast", #  
  '8'  = "Schwann cell"
)

Idents(Fibroblast_1) <- "RNA_snn_res.0.6"
Fibroblast_1 <- RenameIdents(Fibroblast_1, cell_cluster)
Fibroblast_1$subtype <- Idents(Fibroblast_1)

Fibroblast_1$subtype2 <- as.character(Fibroblast_1$subtype)
cell_subtype <- Fibroblast$subtype
names(cell_subtype) <- colnames(Fibroblast)
Fibroblast_1$subtype2[colnames(Fibroblast)] <- as.character(cell_subtype)
table(Fibroblast_1$subtype2)
table(Fibroblast_1$subtype)
table(Fibroblast$subtype)
saveRDS(Fibroblast_1, "Fibroblast_subtype/Fibroblast_allannotation.rds")
######下面合并到alldata的注释先不运行
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Fibroblast_1$subtype2
names(cell_subtype) <- colnames(Fibroblast_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Fibroblast_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Fibroblast_1$subtype2)
saveRDS(all_data, "Fibroblast_subtype/all_data_T_annotation.rds")
saveRDS(Fibroblast_1 ,file = "Fibroblast_subtype/Fibroblast_annotation.rds")



Fibroblast <- readRDS("Fibroblast_subtype/Fibroblast_clean.rds")
#画一个注释完的点图
Fibro_major_marker_ls_3 <- c(
  "FIB" = "PDGFRA",
  "FIB" = "PDGFRB",
  "FIB" = "CCDC80",
  "CXCL12+Fibroblast" = "CXCL12",
  "CXCL12+Fibroblast" = "COLEC11",
  "MEG3+Fibroblsat" = "MEG3",
  #"TCM"="SELL",
  "Myofibroblast" = "FAP",
  "Myofibroblast" = "SULF1",
  "Myofibroblast" = "ACTA2",
  "Myofibroblast" = "POSTN",
  "Myofibroblast" = "COL5A1"

)

pdf(file.path(Fibroblast_out_dir, str_glue("dotplot Fibroblast subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Fibroblast, group.by = "subtype", features = Fibro_major_marker_ls_3)+RotatedAxis()+theme_90angle
dev.off()


Fibroblast$subtype %>% table()
Fibroblast$sample %>% table()
Fibroblast$group %>% table()
## cell props barplot  ----------------------------------------------------------
Fibroblast_meta <- Fibroblast[[]]

Fibroblast_meta$sample <- factor(Fibroblast_meta$sample,levels=c(Fibroblast_meta$sample %>% unique()))

p_col_proj <-
  Fibroblast_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Fibroblast Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Fibroblast_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Fibroblast_meta %>%
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
  filename = file.path(Fibroblast_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#改图
# 提取meta
Fibroblast_meta <- Fibroblast[[]]

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
Fibroblast_meta$Group <- recode(Fibroblast_meta$group, !!!group_map)
Fibroblast_meta$Group <- factor(Fibroblast_meta$Group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Fibroblast_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("Fibroblast_subtype", "Fibroblast_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(Group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("Fibroblast_subtype", "Fibroblast_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = Group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Fibroblast Subtype") +
  labs(title = "Fibroblast Subtype Composition by Group",
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

ggsave(file.path("Fibroblast_subtype", "Fibroblast_StackedBar.pdf"), 
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
  labs(title = "Fibroblast - Cell Type Percentages",
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
ggsave(file.path("Fibroblast_subtype","Fibroblast_Boxplot.pdf"), p, width = 12, height = 5, dpi = 300)
print(p)



















