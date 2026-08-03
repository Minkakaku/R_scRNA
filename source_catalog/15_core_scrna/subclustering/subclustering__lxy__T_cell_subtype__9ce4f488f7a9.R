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

#all_data <- readRDS("output/data/all_data.rds")

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
  
  #"Mesenchyme" = "PDGFRA", 
  
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
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
#options(future.globals.maxSize = 50 * 1024^3)
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

saveRDS(T_cell ,file = "T_cell_subtype/T_cell.rds")
#T_cell <- readRDS("T_cell_subtype/T_cell.rds")
#
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
  "CD4+TRM" = "CXCR6",
  #"CD4+TRM" = "NR4A1/2/3", # 已在CD8+TRM中出现
  #"CD4+TRM" = "MYADM"
  "Tfh" = "CXCR5",
  "Tfh" = "BCL6"
)

Idents(T_cell)<-'RNA_snn_res.0.2'
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(T_cell_markers, file = "T_cell_subtype/T_cell_0.2_all_markers.csv", row.names = TRUE)

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
DotPlot(T_cell, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.2", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(T_cell, group.by = "RNA_snn_res.0.2", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = T_cell_markers)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.2", features = T_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(T_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

saveRDS(T_cell ,file = "T_cell_subtype/T_cell_cluster_0.2.rds")
#T_cell <- readRDS("T_cell_subtype/T_cell_cluster_0.2.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Naive_T", # 
  '1'  = "Proliferative_T", #
  '2'  = "cTEC", #
  '3'  = "Naive_T", # 
  '4'  = "Gamma-delta_T", # 
  '5'  = "Proliferative_T"
)

Idents(T_cell) <- "RNA_snn_res.0.2"

T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$T_annotation <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- T_cell$subtype
names(cell_subtype) <- colnames(T_cell)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$T_annotation[colnames(T_cell)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$T_annotation)
table(all_data$annotation)
table(T_cell$subtype)
saveRDS(all_data, "T_cell_subtype/all_data_T_annotation.rds")


keep <- c(
  "Naive_T",
  "Proliferative_T",
  "Gamma-delta_T"
)

T_cell_clean <- subset(T_cell, subtype %in% keep)

T_cell_clean$subtype <- droplevels(T_cell_clean$subtype)
Idents(T_cell_clean) <- T_cell_clean$subtype

table(T_cell_clean$subtype)
saveRDS(T_cell_clean, "T_cell_subtype/T_cell_clean.rds")
#T_cell <- readRDS("T_cell_subtype/T_cell_clean.rds")

T_cell <- T_cell_clean
#重新分布一下散点
T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} T annotation split.pdf")), height = 6, width = 15)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

#新加入一个点图
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
T_cell_out_dir <- "T_cell_subtype"
T_major_marker_ls_2 <- c(
  "CD8"="CD8A",
  "CD8"="CD8B",
  
  "CD4"="CD4",

  # "Naive"="CCR6",
  "Naive"="LEF1",
  "Naive"="TCF7",
  "Naive"="SELL",
  # "Naive"="CCR7",

  
  "Proliferative_T"="MKI67",
  "Proliferative_T"="H4C3",
  "Proliferative_T"="TUBB",
  "Proliferative_T"="STMN1",
  "Proliferative_T"="TOP1",
  "Proliferative_T"="ANAPC11",
  "Proliferative_T"="MCM7",
  "Proliferative_T"="DNMT1",
  "Proliferative_T"="TOP2A", 
  
  # Gamma-delta T cells (Gamma-deltaT)
  "Gamma-deltaT" = "TRDC"
)
pdf(file.path(T_cell_out_dir, str_glue("dotplot T subtype major_marker.pdf")), width = 25, height = 6)
DotPlot(T_cell, group.by = "subtype", features = T_major_marker_ls_2)+RotatedAxis()+theme_90angle
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
    vars = c(subtype_col
    #, some_genes
    )
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
