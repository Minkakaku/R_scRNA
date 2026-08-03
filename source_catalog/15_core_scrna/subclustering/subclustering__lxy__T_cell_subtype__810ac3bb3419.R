library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")

# normalize
all_data <- all_data %>% NormalizeData()

saveRDS(all_data ,file = "output/data/all_data.rds")



#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

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
  
  "Epi"= "EPCAM",
  
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
  "Pericyte" = "NOTCH3", 
  
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
  "Dendritic cell" = "CST3", 
  
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A", 
  
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




# T_cells -------------------------------------------------------------------
#all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')

T_cell_out_dir <- "T_cell_subtype" 
dir.create(T_cell_out_dir)
T_cell <- subset(all_data, annotation == 'T cell')

all_data$annotation %>% table()
all_data$celltype %>% table()
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
dims <- 16 # 18

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


T_major_marker_ls_1 <- c(
  # CD4+ naive T cells
  "CD4+naïveT" = "CD4",
  "CD4+naïveT" = "CCR7",
  "CD4+naïveT" = "TCF7",
  "CD4+naïveT" = "SELL",
  
  # CD8+ naive T cells
  "CD8+naïveT" = "CD8A",
  "CD8+naïveT" = "CD8B",
  #"CD8+naïveT" = "CCR7",
  #"CD8+naïveT" = "TCF7",
  #"CD8+naïveT" = "SELL",
  
  # CD4+ central memory T cells (CD4+TCM)
  #"CD4+TCM" = "CD4",
  "CD4+TCM" = "CCR6",
  #"CD4+TCM" = "CCR7",
  #"CD4+TCM" = "SELL",
  "CD4+TCM" = "ANXA1",
  "CD4+TCM" = "ANXA2",
  
  # CD8+ central memory T cells (CD8+TCM)
  #"CD8+TCM" = "CD8A",
  #"CD8+TCM" = "CD8B",
  "CD8+TCM" = "GZMK",
  #"CD8+TCM" = "CCR7",
  #"CD8+TCM" = "SELL",
  
  # CD4+ effector memory T cells (CD4+TEM)
  #"CD4+TEM" = "CD4",
  "CD4+TEM" = "KLRG1",
  "CD4+TEM" = "NKG7",
  "CD4+TEM" = "GNLY",
  
  # CD8+ effector memory T cells (CD8+TEM)
  #"CD8+TEM" = "CD8A",
  #"CD8+TEM" = "CD8B",
  #"CD8+TEM" = "GNLY",
  #"CD8+TEM" = "NKG7",
  "CD8+TEM" = "GZMB",
  
  # regulatory T cells (Treg)
  "Treg" = "FOXP3",
  "Treg" = "IL2RA",
  
  # exhausted T cells (exhaustedT)
  "exhaustedT" = "PDCD1",
  "exhaustedT" = "HAVCR2",
  "exhaustedT" = "LAG3",
  
  # Gamma-delta T cells (Gamma-deltaT)
  "Gamma-deltaT" = "KLRB1",
  "Gamma-deltaT" = "TRDV2",
  "Gamma-deltaT" = "TRGV9",
  "Gamma-deltaT" = "TRDC",
  
  # mucosal-associated invariant T cells (MAIT)
  "MAIT" = "SLC4A10",
  #"MAIT" = "KLRB1",
  "MAIT" = "ZBTB16",
  "MAIT" = "NCR3",
  "MAIT" = "RORC",
  "MAIT" = "RORA", 
  
  # Proliferative T cells (ProliferativeT)
  "ProliferativeT" = "MKI67",
  
  # CD8+ Tissue-resident memory T cells (CD8+TRM)
  "CD8+TRM" = "CD6",
  "CD8+TRM" = "XCL1",
  "CD8+TRM" = "XCL2",
  "CD8+TRM" = "MYADM",
  "CD8+TRM" = "CAPG",
  #"CD8+TRM" = "RORA", # 已在MAIT中出现
  "CD8+TRM" = "NR4A1",
  "CD8+TRM" = "NR4A2",
  "CD8+TRM" = "NR4A3",
  "CD8+TRM" = "CD69",
  "CD8+TRM" = "ITGAE",
  
  # CD4+ Tissue-resident memory T cells (CD4+TRM)
  #"CD4+TRM" = "CD69",
  #"CD4+TRM" = "KLRB1", # 已在Gamma-deltaT中出现
  "CD4+TRM" = "PTGER4",
  "CD4+TRM" = "IL7R",
  "CD4+TRM" = "CXCR6"
  #"CD4+TRM" = "NR4A1/2/3", # 已在CD8+TRM中出现
  #"CD4+TRM" = "MYADM"
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
write.csv(T_cell_markers, file = "T_cell_subtype/T_cell_all_markers.csv", row.names = TRUE)

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
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(T_cell ,file = "T_cell_subtype/T_cell_cluster.rds")
T_cell <- readRDS("T_cell_subtype/T_cell_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Fibroblsat-like_T", # 
  '1'  = "CD4+", #偏向CD4+TEM
  '2'  = "Macrophage-like_T", # 
  '3'  = "CD8+", # CD8+TEM
  '4'  = "PT-like_T", # 
  '5'  = "LOH_CD-like_T",#  
  '6'  = "B cell", #  
  '7'  = "Plasma-like_T", #  
  '8'  = "Endothelial-like_T",# 
  '9'  = "LAMP3_DC-like_T", #  
  '10' = "pDC-like_T", # 
  '11' = "Podocyte-like_T",#  
  '12' = "Proliferating T cell"
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


# 保留：0,1,2,3,7,9,12
T_cell <- subset(T_cell, RNA_snn_res.0.6 %in% c(0,1,2,3,7,9,12))
T_cell$RNA_snn_res.0.6 %>% table()


## PCA --------
T_cell <- FindVariableFeatures(T_cell)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot_1.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 16 # 18

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


Idents(T_cell)<-'RNA_snn_res.0.8'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(T_cell_markers, file = "T_cell_subtype/T_cell_all_markers_1.csv", row.names = TRUE)

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
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()



pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker_2.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.8", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()



saveRDS(T_cell ,file = "T_cell_subtype/T_cell_cluster_1.rds")
T_cell <- readRDS("T_cell_subtype/T_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "CD4+", # 
  '1'  = "Fibroblsat-like_T", #
  '2'  = "Macrophage-like_T", # 
  '3'  = "CD8+TEM", # CD8+TEM
  '4'  = "Fibroblsat-like_T", # 
  '5'  = "CD8+Effector",#  
  '6'  = "Macrophage-like_T", #  
  '7'  = "LAMP3_DC-like_T", #  
  '8'  = "Treg",# 
  '9'  = "Macrophage-like_T", #  
  '10' = "Proliferating T cell", # 
  '11' = "Mast cell",#  
  '12' = "Proliferating T cell"
)

Idents(T_cell) <- "RNA_snn_res.0.8"




T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation_1.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split_1.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

keep <- c(
  "CD4+",
  "CD8+TEM",
  "CD8+Effector",
  "Treg",
  "Proliferating T cell"
)

T_cell_clean <- subset(T_cell, subtype %in% keep)

T_cell_clean$subtype <- droplevels(T_cell_clean$subtype)
Idents(T_cell_clean) <- T_cell_clean$subtype

table(T_cell_clean$subtype)
saveRDS(T_cell_clean, "T_cell_subtype/T_cell_clean.rds")

T_cell <- T_cell_clean

T_cell <- RunUMAP(T_cell, dims = 1:16,n.neighbors = 50,return.model = TRUE)

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation_1.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation split_1.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()


T_cell_1 <- readRDS("T_cell_subtype/T_cell_cluster.rds")
cell_cluster <- c(
  '0'  = "Fibroblsat-like_T", # 
  '1'  = "CD4+", #偏向CD4+TEM
  '2'  = "Macrophage-like_T", # 
  '3'  = "CD8+", # CD8+TEM
  '4'  = "PT-like_T", # 
  '5'  = "LOH_CD-like_T",#  
  '6'  = "B cell", #  
  '7'  = "Plasma-like_T", #  
  '8'  = "Endothelial-like_T",# 
  '9'  = "LAMP3_DC-like_T", #  
  '10' = "pDC-like_T", # 
  '11' = "Podocyte-like_T",#  
  '12' = "Proliferating T cell"
)

Idents(T_cell_1) <- "RNA_snn_res.0.6"
T_cell_1 <- RenameIdents(T_cell_1, cell_cluster)
T_cell_1$subtype <- Idents(T_cell_1)

T_cell_1$subtype2 <- as.character(T_cell_1$subtype)
cell_subtype <- T_cell$subtype
names(cell_subtype) <- colnames(T_cell)
T_cell_1$subtype2[colnames(T_cell)] <- as.character(cell_subtype)
table(T_cell_1$subtype2)
table(T_cell_1$subtype)
table(T_cell$subtype)
######
all_data$annotation %>% table()

#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- T_cell_1$subtype2
names(cell_subtype) <- colnames(T_cell_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(T_cell_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(T_cell_1$subtype2)
saveRDS(all_data, "T_cell_subtype/all_data_T_annotation.rds")
saveRDS(T_cell_1 ,file = "T_cell_subtype/T_cell_annotation.rds")



T_cell <- readRDS("T_cell_subtype/T_cell_clean.rds")
#画一个注释完的点图
T_major_marker_ls_2 <- c(
  "CD8"="CD8A",
  "CD8"="CD8B",
  
  "CD4"="CD4",
  
  # "Naive"="CCR6",
  "Naive"="LEF1",
  "Naive"="TCF7",
  #"TCM"="SELL",
  
  "Effector"="GZMB",
  "Effector"="CX3CR1", #https://pubmed.ncbi.nlm.nih.gov/35361683/
  
  "TEM"="GZMK", 
  "TEM"="GZMA", 
  
  # regulatory T cells (Treg)
  "Treg" = "FOXP3",
  "Treg" = "IL2RA",
  
  #"CD4_Tfh"="TOX2",
  
  "Proliferative_T"="MKI67",
  "Proliferative_T"="TOP2A"
)

pdf(file.path(T_cell_out_dir, str_glue("dotplot T subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(T_cell, group.by = "subtype", features = T_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()

# 基础版小提琴图（直接出图）
VlnPlot(object = T_cell,
        features = "nCount_RNA",  # 你要画的指标
        group.by = "seurat_clusters",  # 按cluster分组
        pt.size = 0,  # 去掉散点，更干净
        fill = "seurat_clusters"  # 每个cluster自动配色
)

VlnPlot(object = T_cell,
        features = "total_counts",  # 你要画的指标
        group.by = "seurat_clusters",  # 按cluster分组
        pt.size = 0,  # 去掉散点，更干净
        fill = "seurat_clusters"  # 每个cluster自动配色
)

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

#改图
# 提取meta
T_cell_meta <- T_cell[[]]

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
T_cell_meta$Group <- recode(T_cell_meta$group, !!!group_map)
T_cell_meta$Group <- factor(T_cell_meta$Group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- T_cell_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("T_cell_subtype", "Tcell_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(Group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("T_cell_subtype", "Tcell_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = Group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "T Cell Subtype") +
  labs(title = "T Cell Subtype Composition by Group",
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

ggsave(file.path("T_cell_subtype", "T_cell_StackedBar.pdf"), 
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
  labs(title = "T Cell - Cell Type Percentages",
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
ggsave(file.path("T_cell_subtype","T_cell_Boxplot.pdf"), p, width = 12, height = 8, dpi = 300)
print(p)




















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
    vars = c(subtype_col, some_genes)
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
  file.path("T_cell_subtype/spatial plot of celltype in each sample.pdf"),
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
      file.path(str_glue("T_cell_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })


## volcano -------------------(下面的是差异基因和富集分析，可以不用跑)
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