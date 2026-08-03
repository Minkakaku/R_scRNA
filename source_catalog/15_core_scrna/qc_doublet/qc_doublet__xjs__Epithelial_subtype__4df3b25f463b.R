library(Seurat)
library(tidyverse)

Major_marker_ls <- c(
  "Epithelial" = "EPCAM", 
  "Epithelial" = "KRT5", 
  "Epithelial" = "KRT14", 
  "Epithelial" = "KRT20", 
  
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
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  "Proliferating" = "MKI67", 
  "Proliferating" = "TOP2A", 
  
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
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$annotation %>% table()


Epithelial <- subset(all_data, annotation %in% c('AT2 cell', 'Club cell', 'Mucous cell', 'Serous cell', 'Tumor cell'))
Epithelial$nFeature_Spatial.008um %>% quantile()
Epithelial$nCount_Spatial.008um %>% quantile()
Epithelial$pct_counts_mt %>% quantile()
Epithelial$pct_counts_mt %>% hist()
Epithelial_raw <- Epithelial
saveRDS(Epithelial_raw, file.path(Epithelial_out_dir, 'Epithelial_raw_harmony.rds'))
Epithelial <- Epithelial_raw
Epithelial <- subset(Epithelial, nFeature_Spatial.008um > 45)
FeaturePlot(Epithelial, features = 'nFeature_Spatial.008um')
FeaturePlot(Epithelial, features = 'nCount_Spatial.008um')
FeaturePlot(Epithelial, features = 'pct_counts_mt')

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

# DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_Spatial.008um')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'Epithelial_subtype' 
dir.create(Epithelial_out_dir)


### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 10

Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE) + NoLegend()
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

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 

reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
Epithelial <- 
  FindNeighbors(
    Epithelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster_raw.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "Spatial.008um_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




Epithelial_major_marker_ls <- c(
  "AT1" = "AGER",
  "AT1" = "CAV1",
  "AT1" = "NTM",
  "AT1" = "TNNC1",
  "AT1" = "CAV2",
  "AT1" = "MYL9",
  "AT1" = "KRT7",
  "AT1" = "EMP2",
  "AT1" = "RTKN2",
  "AT1" = "ADIRF",
  "AT2" = "SFTPC",
  "AT2" = "SFTPA1",
  "AT2" = "SFTPA2",
  "AT2" = "SFTPD",
  "AT2" = "NAPSA",
  "AT2" = "SFTPB",
  "AT2" = "SFTA2",
  "AT2" = "NPC2",
  # "AT2" = "PGC",
  "AT2" = "SERPINA1",
  "Basal" = "KRT17",
  "Basal" = "S100A2",
  "Basal" = "KRT15",
  "Basal" = "KRT5",
  "Basal" = "MMP10",
  "Basal" = "CXCL1",
  "Basal" = "CXCL8",
  "Basal" = "PHLDA1",
  "Basal" = "DST",
  "Basal" = "ALDH3A1",
  # "Basaloid" = "SERPINB3",
  "Basaloid" = "CSTA",
  "Basaloid" = "AQP5",
  "Basaloid" = "PRSS23",
  "Basaloid" = "DUSP5",
  # "Basaloid" = "MDK",
  "Basaloid" = "FAM3D",
  # "Basaloid" = "S100P",
  "Basaloid" = "KRT19",
  "Basaloid" = "TACSTD2",
  "Ciliated" = "CAPS",
  "Ciliated" = "TMEM190",
  "Ciliated" = "C9orf24",
  "Ciliated" = "C20orf85",
  "Ciliated" = "TPPP3",
  "Ciliated" = "RSPH1",
  "Ciliated" = "FAM183A",
  "Ciliated" = "C1orf194",
  "Ciliated" = "OMG",
  "Ciliated" = "TSPAN1",
  "Deuterosomal" = "CDC20B",
  "Deuterosomal" = "HES6",
  "Deuterosomal" = "CCNO",
  "Deuterosomal" = "CARHSP1",
  # "Deuterosomal" = "BPIFA1",
  # "Deuterosomal" = "BPIFB1",
  # "Deuterosomal" = "LCN2",
  "Deuterosomal" = "SOX4",
  "Deuterosomal" = "MT2A",
  "Deuterosomal" = "SCGB1A1",
  "Hillocks" = "KRT6A",
  "Hillocks" = "KIAA0101",
  "Hillocks" = "STMN1",
  "Hillocks" = "S100A8",
  "Hillocks" = "H2AFZ",
  "Hillocks" = "TUBB",
  # "Hillocks" = "KRT5",
  # "Hillocks" = "S100A2",
  "Hillocks" = "HIST1H4C",
  "Hillocks" = "S100A9",
  # "Mucous" = "BPIFB1",
  "Mucous" = "MUC5B",
  # "Mucous" = "LCN2",
  "Mucous" = "SCGB3A1",
  "Mucous" = "WFDC2",
  # "Mucous" = "BPIFA1",
  # "Mucous" = "TFF3",
  # "Mucous" = "MSMB",
  # "Mucous" = "SCGB1A1",
  "Mucous" = "BPIFB2",
  "RAS" = "TMEM45A",
  "RAS" = "MGP",
  "RAS" = "KLK11",
  # "RAS" = "SCGB1A1",
  # "RAS" = "SCGB3A1",
  # "RAS" = "WFDC2",
  # "RAS" = "SCGB3A2",
  # "RAS" = "SOX4",
  "RAS" = "HES1",
  "RAS" = "BPIFA1",
  "SCGB3A2+" = "CALCA",
  "SCGB3A2+" = "TESC",
  "SCGB3A2+" = "APOD",
  "SCGB3A2+" = "SPINK1",
  "SCGB3A2+" = "MDK",
  "SCGB3A2+" = "CEACAM6",
  "SCGB3A2+" = "SPINK5",
  # "SCGB3A2+" = "MGP",
  "SCGB3A2+" = "SCGB3A2",
  "SCGB3A2+" = "PGC",
  "Secretory" = "MSMB",
  "Secretory" = "TFF3",
  "Secretory" = "BPIFB1",
  # "Secretory" = "BPIFA1",
  "Secretory" = "LCN2",
  "Secretory" = "S100P",
  "Secretory" = "SERPINB3",
  "Secretory" = "MUC5AC",
  "Secretory" = "LYPD2",
  # "Secretory" = "SCGB1A1",
  "Serous" = "LTF",
  "Serous" = "ZG16B",
  "Serous" = "PIP",
  "Serous" = "PRB1",
  "Serous" = "PRB3",
  "Serous" = "RP11-1143G9.4",
  "Serous" = "PRR4",
  "Serous" = "LYZ",
  "Serous" = "PRH2",
  "Serous" = "PRB4",
  # "Translational AT2" = "CAV1",
  # "Translational AT2" = "TNNC1",
  # "Translational AT2" = "AGER",
  # "Translational AT2" = "EMP2",
  "Translational AT2" = "HOPX",
  # "Translational AT2" = "KRT7",
  "Translational AT2" = "PEG10",
  "Translational AT2" = "CLDN18",
  # "Translational AT2" = "NTM",
  # "Translational AT2" = "SCGB3A2",
  "Adventitial" = "PI16",
  "Adventitial" = "LEPR",
  # "Adventitial" = "CFD",
  "Adventitial" = "SCARA5",
  # "Adventitial" = "APOD",
  "Adventitial" = "CCDC80",
  "Adventitial" = "DCN",
  "Adventitial" = "SFRP2",
  "Adventitial" = "SERPINF1",
  "Pathologic" = "CTHRC1",
  "Pathologic" = "POSTN",
  "Pathologic" = "THBS2",
  "Pathologic" = "INHBA",
  "Pathologic" = "COL1A1",
  "Pathologic" = "COL3A1",
  "Pericyte" = "MCAM",
  "Pericyte" = "NOTCH3",
  "Pericyte" = "RGS5",
  "Pericyte" = "COX4I2",
  "Pericyte" = "NDUFA4L2",
  "Pericyte" = "IGFBP7",
  "Alveolar-like IEG-high" = "TCF21",
  "Alveolar-like IEG-high" = "ITGA8",
  "Alveolar-like IEG-high" = "FOSB",
  "Alveolar-like IEG-high" = "EGR1",
  "Alveolar-like IEG-high" = "FOS",
  "Alveolar-like IEG-high" = "JUN",
  "Alveolar-like IEG-high" = "JUND",
  "Matrix-remodeling" = "LTBP2",
  "Matrix-remodeling" = "COL8A1",
  "Matrix-remodeling" = "MMP2",
  "Matrix-remodeling" = "HMCN1",
  "Matrix-remodeling" = "IGFBP4",
  "Matrix-remodeling" = "COL14A1",
  "Alveolar-like homeostatic" = "INMT",
  "Alveolar-like homeostatic" = "FGFR4",
  "Alveolar-like homeostatic" = "MFAP4",
  "Alveolar-like homeostatic" = "ADH1B",
  "Alveolar-like homeostatic" = "SEPP1",
  "Alveolar-like homeostatic" = "GPC3",
  "Alveolar-like homeostatic" = "A2M",
  "Alveolar-like homeostatic" = "LIMCH1",
  "Alveolar-like homeostatic" = "PTGDS"
  # "SMC" = "MYH11",
  # "SMC" = "CNN1",
  # "SMC" = "LMOD1",
  # "SMC" = "ACTA2",
  # "SMC" = "TAGLN",
  # "SMC" = "DES",
  # "SMC" = "MYLK",
  # "AMs" = "FABP4",
  # "AMs" = "SERPING1",
  # "AMs" = "APOC1",
  # "AMs" = "CD52",
  # "AMs" = "PCOLCE2",
  # # "AMs" = "INHBA",
  # "AMs" = "MME",
  # "AMs" = "CD10",
  # "AMs" = "CES1",
  # "AMs" = "GPD1",
  # "AMs" = "CFD",
  # "AMs" = "RBP4",
  # "AMs" = "ITIH5",
  # "AMs" = "SCD",
  # "IMs" = "LGMN",
  # "IMs" = "MARCKS",
  # "IMs" = "FOLR2",
  # "IMs" = "SELENOP",
  # "IMs" = "F13A1",
  # "IMs" = "SLC40A1",
  # "IMs" = "SPP1",
  # "IMs" = "PLA2G7",
  # "IMs" = "MMP9",
  # "IMs" = "HAMP",
  # # "recMac" = "S100A8",
  # "recMac" = "FN1",
  # "recMac" = "VCAN",
  # "recMac" = "RNASE1",
  # "recMac" = "FCN1",
  # "recMac" = "EMP1",
  # "recMac" = "CD14",
  # "cDC1" = "CLEC9A",
  # "cDC1" = "XCR1",
  # "cDC1" = "CADM1",
  # "cDC1" = "BTLA",
  # "cDC1" = "IDO1",
  # "cDC1" = "DNASE1L3",
  # "cDC1" = "CCND1",
  # "cDC1" = "S100B",
  # "cDC1" = "LGALS2",
  # "cDC2" = "CD1C",
  # "cDC2" = "CD1E",
  # "cDC2" = "CLEC10A",
  # "cDC2" = "FCER1A",
  # "cDC2" = "FCGR2B",
  # "cDC2" = "CCL17",
  # "cDC2" = "PKIB",
  # # "cDC2" = "LGALS2",
  # "cDC2" = "BIRC3",
  # "pDCs" = "GZMB",
  # "pDCs" = "JCHAIN",
  # "pDCs" = "MZB1",
  # "pDCs" = "CLEC4C",
  # "pDCs" = "CLIC3",
  # "pDCs" = "SPIB",
  # "pDCs" = "TCF4"
)
DotPlot(Epithelial, group.by = "res_1.0", features = Epithelial_major_marker_ls)+RotatedAxis()
ggsave(file.path(Epithelial_out_dir, 'dotplot_all_data_epi_cluster.pdf'), width = 35, height = 8)
DotPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.2", features = Epithelial_major_marker_ls)+RotatedAxis()
ggsave(file.path(Epithelial_out_dir, 'dotplot_epi_sub_cluster_res0.2.pdf'), width = 35, height = 8)
DotPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.4", features = Epithelial_major_marker_ls)+RotatedAxis()
ggsave(file.path(Epithelial_out_dir, 'dotplot_epi_sub_cluster_res0.4.pdf'), width = 35, height = 8)


DotPlot(Epithelial, group.by = "Spatial.008um_snn_res.0.4", features = Major_marker_ls)+RotatedAxis()
ggsave(file.path(Epithelial_out_dir, 'dotplot_epi_sub_cluster_res0.4_maintype.pdf'), width = 25, height = 8)


saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial_gene_over_45.rds'))

























# find markers
Idents(Epithelial) <- "RNA_snn_res.0.8"
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
  dplyr::filter(avg_log2FC > 0.75) %>% 
  dplyr::filter(pct.1 - pct.2 > 0.1) %>% 
  arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 30) %>%
  # arrange(desc(pct.1)) %>%
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
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Epithelial_top5_markers)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)
pdf(file.path(Epithelial_out_dir, 'for REST tumor cell.pdf'), width = 8, height = 6)
VlnPlot(subset(Epithelial, REST > 0), features = 'REST', group.by = 'RNA_snn_res.0.8', pt.size = 0)
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)
FeaturePlot(Epithelial, features = 'REST', reduction = 'harmonyumap')
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




## annotation -----------
Epithelial$RNA_snn_res.0.8 %>% table()
table(Epithelial$RNA_snn_res.0.8, Epithelial$sample)
cell_cluster <- c(
  '0' = "Tumor cell", # MIAT
  '1' = "Tumor cell", # MIAT
  '2' = "Tumor cell", # CCNB1 https://www.x-mol.com/paper/1707085704578486272/t#:~:text=%E5%88%A9%E7%94%A8%E8%82%BA%E8%85%BA%E7%99%8C%E6%95%B0%E6%8D%AE%EF%BC%8C%E6%88%91%E4%BB%AC%E7%A0%94%E7%A9%B6%E4%BA%86%E5%8F%82%E4%B8%8E%20CCNB1%20%E8%B0%83%E8%8A%82%E7%9A%84%E6%BD%9C%E5%9C%A8%E4%B8%8A%E6%B8%B8%E9%9D%9E%E7%BC%96%E7%A0%81%20RNA%EF%BC%8C%E5%B9%B6%E9%AA%8C%E8%AF%81%E4%BA%86%20CCNB1%20%E5%AF%B9%E4%BA%8E%E8%AF%A5%E7%96%BE%E7%97%85%E9%83%A8%E4%BD%8D%E7%9A%84%E8%9B%8B%E7%99%BD%E6%B0%B4%E5%B9%B3%E5%92%8C%E9%A2%84%E5%90%8E%E4%BB%B7%E5%80%BC%E3%80%82,CCNB1%20%E5%9C%A8%E5%A4%A7%E5%A4%9A%E6%95%B0%E7%99%8C%E7%97%87%E4%B8%AD%E9%AB%98%E8%A1%A8%E8%BE%BE%EF%BC%8C%E5%90%AF%E5%8A%A8%E5%AD%90%E7%94%B2%E5%9F%BA%E5%8C%96%E9%99%8D%E4%BD%8E%E3%80%82%20CCNB1%E7%9A%84%E5%9F%BA%E5%9B%A0%E8%A1%A8%E8%BE%BE%E4%B8%8E%E8%82%BF%E7%98%A4%E6%82%A3%E8%80%85%E7%9A%84%E4%B8%8D%E8%89%AF%E9%A2%84%E5%90%8E%E5%91%88%E6%AD%A3%E7%9B%B8%E5%85%B3%EF%BC%8C%E8%BF%99%E4%BA%9B%E7%BB%93%E6%9E%9C%E5%9C%A8%E8%82%BA%E8%85%BA%E7%99%8C%E7%9A%84%E8%9B%8B%E7%99%BD%E8%B4%A8%E6%B0%B4%E5%B9%B3%E4%B8%8A%E5%BE%97%E5%88%B0%E4%BA%86%E8%AF%81%E5%AE%9E%E3%80%82%20CCNB1%20%E8%A1%A8%E8%BE%BE%E4%B8%8E%20T%20%E8%BE%85%E5%8A%A9%E7%BB%86%E8%83%9E%E6%B5%B8%E6%B6%A6%E7%9B%B8%E5%85%B3%EF%BC%8C%E8%BF%99%E8%BF%9B%E4%B8%80%E6%AD%A5%E4%B8%8E%E6%9F%90%E4%BA%9B%E7%99%8C%E7%97%87%EF%BC%88%E5%8C%85%E6%8B%AC%E8%82%BE%E9%80%8F%E6%98%8E%E7%BB%86%E8%83%9E%E7%99%8C%E5%92%8C%E8%82%BA%E8%85%BA%E7%99%8C%EF%BC%89%E7%9A%84%E4%B8%8D%E8%89%AF%E9%A2%84%E5%90%8E%E7%9B%B8%E5%85%B3%E3%80%82
  '3' = "Tumor cell", # MIAT
  '4' = "Tumor cell", # MIAT
  '5' = "Tumor cell", # MIAT
  '6' = "low-nCounts",
  '7' = "Tumor cell", # MIAT
  '8' = "Tumor cell", # MIAT 
  '9' = "Tumor cell", # MIAT
  '10' = "Tumor cell", # MIAT
  '11' = "Tumor cell", # MIAT
  '12' = "Tumor cell", # MIAT
  '13' = "Anti-tumor cell", # SERPINA3 https://pmc.ncbi.nlm.nih.gov/articles/PMC10552721/
  '14' = "Tumor cell", # MIAT
  '15' = "Mucous cell", # MUC5B
  '16' = "Mix", 
  '17' = "low-nCounts",
  '18' = "low-nCounts",
  '19' = "Ciliated cell", # DLEC1, DNAI1
  '20' = "Tumor cell", # MIAT
  '21' = "low-nCounts",
  '22' = "low-nCounts",
  '23' = "Mix" # EPCAM, CD83
)


Idents(Epithelial) <- "RNA_snn_res.0.8"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# Epithelial$group <- str_replace_all(Epithelial$sample, c(
#   'CW_SH' = 'non-PCR', 
#   'CW_SQ' = 'non-PCR', 
#   'HJD_SH' = 'non-PCR', 
#   'HJD_SQ' = 'non-PCR', 
#   'LGRO_SH' = 'non-PCR', 
#   'LGRO_SQ' = 'non-PCR', 
#   'WBS_SH' = 'non-PCR', 
#   'WBS_SQ' = 'non-PCR', 
#   'WLJ_SH' = 'non-PCR', 
#   'WLJ_SQ' = 'non-PCR', 
#   'JJW_SH' = 'PCR', 
#   'JJW_SQ' = 'PCR', 
#   'YBK_SH' = 'PCR', 
#   'YBK_SQ' = 'PCR'
# ))
# 
# pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 7, width = 18)
# DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - PCR -SH.pdf")), height = 6, width = 8)
# DimPlot(subset(Epithelial, sample %in% c('JJW_SH', 'YBK_SH')), group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# dev.off()




## cell props barplot  --------
Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, levels=c(
  'CW_SQ', 
  'CW_SH', 
  'HJD_SQ',
  'HJD_SH',
  'LGRO_SQ',
  'LGRO_SH',
  'WBS_SQ',
  'WBS_SH',
  'WLJ_SQ',
  'WLJ_SH',
  'JJW_SQ',
  'JJW_SH',
  'YBK_SQ',
  'YBK_SH'
)
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
  filename = file.path(Epithelial_out_dir, "cell_freq_sample.pdf"),
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
  filename = file.path(Epithelial_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  'Tumor cell' = 'MIAT', 
  'Tumor cell' = 'DSG3', 
  'Tumor cell' = 'MUC1', 
  
  "Anti-tumor cell" = 'FAM107A', 
  "Anti-tumor cell" = 'HSD17B2', 
  "Anti-tumor cell" = 'SERPINA3', 
  
  "Mucous cell" = 'MUC5B', 
  
  # 'Mix' = 'AGER',
  # 'Mix' = 'ABCA3',
  
  "Ciliated cell" = 'DLEC1', 
  "Ciliated cell" = 'DNAI1'
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 10, height = 5)
DotPlot(subset(Epithelial, subtype %in% c("Tumor cell", "Anti-tumor cell", "Mucous cell", "Ciliated cell")), group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial.rds'))



# spatial plot -------------------
## 映射-------------------------------
library(tibble)
library(data.table)
all_data$annotation %>% table()
all_data$group %>% table()
subtype <- 
  Epithelial[['subtype']] %>% 
  rownames_to_column('cells')
all_data_SH <- subset(all_data, group == 'SH')

meta.data.t <- 
  subtype %>%
  right_join(
    all_data_SH[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Endothelial cell', 'Fibroblast', 'Macrophage', 'Plasma cell', 'T cell')),
                         as.character(annotation), 
                         as.character(subtype)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')

all_data_SH <- AddMetaData(all_data_SH,meta.data.t['annotation2'])
all_data_SH$annotation %>% table()
all_data_SH$annotation2 %>% table()


## Spatial Dimplot -----------------
all_data_SH_meta <- 
  all_data_SH[[]] %>% 
  rownames_to_column('cell_id')
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


ct_sp_out_dir <- file.path(Epithelial_out_dir, 'Spatial dimplot')
dir.create(ct_sp_out_dir)
ct_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    if(sample %in% unique(all_data_SH$sample)){
      anno_data   <- all_data_SH_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_SH_meta, by = "cell_id")
      
      p <- dat %>%
        ggplot(aes(x = vertex_x, y = vertex_y)) +
        geom_polygon(aes(fill = annotation2, group = cell_id))  +
        scale_y_reverse() +
        coord_equal() + 
        ggsci::scale_fill_igv() +
        theme_void() +
        labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    }
    else{
      return(NULL)
    }
    
    p
  }) %>%
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir, str_glue("{name} - spatial dimplot.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


ct_sp_out_dir2 <- file.path(Epithelial_out_dir, 'Spatial dimplot_split')
dir.create(ct_sp_out_dir2)
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(all_data_SH$sample)){
      anno_data   <- all_data_SH_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_SH_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$annotation2 %>% unique(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$annotation2 %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$annotation2 <- as.character(dat2$annotation2)
          dat2$annotation2[dat2$annotation2 != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = annotation2, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{name} - spatial dimplot - split by annotation2.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

all_data_SH_meta$annotation2 %>% 
  unique() %>% 
  setNames(., .) %>% 
  # .[1] %>% 
  iwalk(function(ct, ct_name){
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{ct_name} - spatial dimplot - split by sample.pdf")),
      width  = 20,
      height = 20
    )
    ct_sp_p_ls2 %>% 
      walk(function(p){ # browser()
        if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
        if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
        print(p[[ct_name]])
      })
    dev.off()
  })




# REST tumor annotation ----------
Epithelial_out_dir <-'Epithelial_SH_subtype'
Epithelial <- readRDS("~/Epithelial_SH_subtype/Epithelial.rds")
Epithelial$group %>% table()
