library(Seurat)
library(tidyverse)

Major_marker_ls <- c(
  "Macrophage" = "EPCAM", 
  "Macrophage" = "KRT5", 
  "Macrophage" = "KRT14", 
  "Macrophage" = "KRT20", 
  
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
  # "T cell" = "TRAC", 
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
  # "Pericyte" = "PLVAP", 
  
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


Macrophage <- subset(all_data, annotation == 'Macrophage')
Macrophage$nFeature_Spatial.008um %>% quantile()
Macrophage$nCount_Spatial.008um %>% quantile()
Macrophage$pct_counts_mt %>% quantile()
Macrophage$pct_counts_mt %>% hist()
Macrophage_raw <- Macrophage
Macrophage <- Macrophage_raw
Macrophage <- subset(Macrophage, nFeature_Spatial.008um > 45)
FeaturePlot(Macrophage, features = 'nFeature_Spatial.008um')
FeaturePlot(Macrophage, features = 'nCount_Spatial.008um')
FeaturePlot(Macrophage, features = 'pct_counts_mt')

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

# DefaultAssay(Macrophage) <- "RNA"
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_Spatial.008um')
Macrophage <- RunPCA(Macrophage)

Macrophage_out_dir <-'Macrophage_subtype' 
dir.create(Macrophage_out_dir)


### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 10

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "group", label = TRUE) + NoLegend()
dev.off()


Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
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

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 

reduction.use <- 'pca' # 'harmony' 'pca'
### cluster -------
Macrophage <- 
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Spatial.008um_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




Macrophage_major_marker_ls <- c(
  # "AT1" = "AGER",
  # "AT1" = "CAV1",
  # "AT1" = "NTM",
  # "AT1" = "TNNC1",
  # "AT1" = "CAV2",
  # "AT1" = "MYL9",
  # "AT1" = "KRT7",
  # "AT1" = "EMP2",
  # "AT1" = "RTKN2",
  # "AT1" = "ADIRF",
  # "AT2" = "SFTPC",
  # "AT2" = "SFTPA1",
  # "AT2" = "SFTPA2",
  # "AT2" = "SFTPD",
  # "AT2" = "NAPSA",
  # "AT2" = "SFTPB",
  # "AT2" = "SFTA2",
  # "AT2" = "NPC2",
  # # "AT2" = "PGC",
  # "AT2" = "SERPINA1",
  # "Basal" = "KRT17",
  # "Basal" = "S100A2",
  # "Basal" = "KRT15",
  # "Basal" = "KRT5",
  # "Basal" = "MMP10",
  # "Basal" = "CXCL1",
  # "Basal" = "CXCL8",
  # "Basal" = "PHLDA1",
  # "Basal" = "DST",
  # "Basal" = "ALDH3A1",
  # # "Basaloid" = "SERPINB3",
  # "Basaloid" = "CSTA",
  # "Basaloid" = "AQP5",
  # "Basaloid" = "PRSS23",
  # "Basaloid" = "DUSP5",
  # # "Basaloid" = "MDK",
  # "Basaloid" = "FAM3D",
  # # "Basaloid" = "S100P",
  # "Basaloid" = "KRT19",
  # "Basaloid" = "TACSTD2",
  # "Ciliated" = "CAPS",
  # "Ciliated" = "TMEM190",
  # "Ciliated" = "C9orf24",
  # "Ciliated" = "C20orf85",
  # "Ciliated" = "TPPP3",
  # "Ciliated" = "RSPH1",
  # "Ciliated" = "FAM183A",
  # "Ciliated" = "C1orf194",
  # "Ciliated" = "OMG",
  # "Ciliated" = "TSPAN1",
  # "Deuterosomal" = "CDC20B",
  # "Deuterosomal" = "HES6",
  # "Deuterosomal" = "CCNO",
  # "Deuterosomal" = "CARHSP1",
  # # "Deuterosomal" = "BPIFA1",
  # # "Deuterosomal" = "BPIFB1",
  # # "Deuterosomal" = "LCN2",
  # "Deuterosomal" = "SOX4",
  # "Deuterosomal" = "MT2A",
  # "Deuterosomal" = "SCGB1A1",
  # "Hillocks" = "KRT6A",
  # "Hillocks" = "KIAA0101",
  # "Hillocks" = "STMN1",
  # "Hillocks" = "S100A8",
  # "Hillocks" = "H2AFZ",
  # "Hillocks" = "TUBB",
  # # "Hillocks" = "KRT5",
  # # "Hillocks" = "S100A2",
  # "Hillocks" = "HIST1H4C",
  # "Hillocks" = "S100A9",
  # # "Mucous" = "BPIFB1",
  # "Mucous" = "MUC5B",
  # # "Mucous" = "LCN2",
  # "Mucous" = "SCGB3A1",
  # "Mucous" = "WFDC2",
  # # "Mucous" = "BPIFA1",
  # # "Mucous" = "TFF3",
  # # "Mucous" = "MSMB",
  # # "Mucous" = "SCGB1A1",
  # "Mucous" = "BPIFB2",
  # "RAS" = "TMEM45A",
  # "RAS" = "MGP",
  # "RAS" = "KLK11",
  # # "RAS" = "SCGB1A1",
  # # "RAS" = "SCGB3A1",
  # # "RAS" = "WFDC2",
  # # "RAS" = "SCGB3A2",
  # # "RAS" = "SOX4",
  # "RAS" = "HES1",
  # "RAS" = "BPIFA1",
  # "SCGB3A2+" = "CALCA",
  # "SCGB3A2+" = "TESC",
  # "SCGB3A2+" = "APOD",
  # "SCGB3A2+" = "SPINK1",
  # "SCGB3A2+" = "MDK",
  # "SCGB3A2+" = "CEACAM6",
  # "SCGB3A2+" = "SPINK5",
  # # "SCGB3A2+" = "MGP",
  # "SCGB3A2+" = "SCGB3A2",
  # "SCGB3A2+" = "PGC",
  # "Secretory" = "MSMB",
  # "Secretory" = "TFF3",
  # "Secretory" = "BPIFB1",
  # # "Secretory" = "BPIFA1",
  # "Secretory" = "LCN2",
  # "Secretory" = "S100P",
  # "Secretory" = "SERPINB3",
  # "Secretory" = "MUC5AC",
  # "Secretory" = "LYPD2",
  # # "Secretory" = "SCGB1A1",
  # "Serous" = "LTF",
  # "Serous" = "ZG16B",
  # "Serous" = "PIP",
  # "Serous" = "PRB1",
  # "Serous" = "PRB3",
  # "Serous" = "RP11-1143G9.4",
  # "Serous" = "PRR4",
  # "Serous" = "LYZ",
  # "Serous" = "PRH2",
  # "Serous" = "PRB4",
  # # "Translational AT2" = "CAV1",
  # # "Translational AT2" = "TNNC1",
  # # "Translational AT2" = "AGER",
  # # "Translational AT2" = "EMP2",
  # "Translational AT2" = "HOPX",
  # # "Translational AT2" = "KRT7",
  # "Translational AT2" = "PEG10",
  # "Translational AT2" = "CLDN18",
  # # "Translational AT2" = "NTM",
  # # "Translational AT2" = "SCGB3A2",
  # "Adventitial" = "PI16",
  # "Adventitial" = "LEPR",
  # # "Adventitial" = "CFD",
  # "Adventitial" = "SCARA5",
  # # "Adventitial" = "APOD",
  # "Adventitial" = "CCDC80",
  # "Adventitial" = "DCN",
  # "Adventitial" = "SFRP2",
  # "Adventitial" = "SERPINF1",
  # "Pathologic" = "CTHRC1",
  # "Pathologic" = "POSTN",
  # "Pathologic" = "THBS2",
  # "Pathologic" = "INHBA",
  # "Pathologic" = "COL1A1",
  # "Pathologic" = "COL3A1",
  # "Pericyte" = "MCAM",
  # "Pericyte" = "NOTCH3",
  # "Pericyte" = "RGS5",
  # "Pericyte" = "COX4I2",
  # "Pericyte" = "NDUFA4L2",
  # "Pericyte" = "IGFBP7",
  # "Alveolar-like IEG-high" = "TCF21",
  # "Alveolar-like IEG-high" = "ITGA8",
  # "Alveolar-like IEG-high" = "FOSB",
  # "Alveolar-like IEG-high" = "EGR1",
  # "Alveolar-like IEG-high" = "FOS",
  # "Alveolar-like IEG-high" = "JUN",
  # "Alveolar-like IEG-high" = "JUND",
  # "Matrix-remodeling" = "LTBP2",
  # "Matrix-remodeling" = "COL8A1",
  # "Matrix-remodeling" = "MMP2",
  # "Matrix-remodeling" = "HMCN1",
  # "Matrix-remodeling" = "IGFBP4",
  # "Matrix-remodeling" = "COL14A1",
  # "Alveolar-like homeostatic" = "INMT",
  # "Alveolar-like homeostatic" = "FGFR4",
  # "Alveolar-like homeostatic" = "MFAP4",
  # "Alveolar-like homeostatic" = "ADH1B",
  # "Alveolar-like homeostatic" = "SEPP1",
  # "Alveolar-like homeostatic" = "GPC3",
  # "Alveolar-like homeostatic" = "A2M",
  # "Alveolar-like homeostatic" = "LIMCH1",
  # "Alveolar-like homeostatic" = "PTGDS"
  # "SMC" = "MYH11",
  # "SMC" = "CNN1",
  # "SMC" = "LMOD1",
  # "SMC" = "ACTA2",
  # "SMC" = "TAGLN",
  # "SMC" = "DES",
  # "SMC" = "MYLK",
  "AMs" = "FABP4",
  "AMs" = "SERPING1",
  "AMs" = "APOC1",
  "AMs" = "CD52",
  "AMs" = "PCOLCE2",
  # "AMs" = "INHBA",
  "AMs" = "MME",
  "AMs" = "CD10",
  "AMs" = "CES1",
  "AMs" = "GPD1",
  "AMs" = "CFD",
  "AMs" = "RBP4",
  "AMs" = "ITIH5",
  "AMs" = "SCD",
  "IMs" = "LGMN",
  "IMs" = "MARCKS",
  "IMs" = "FOLR2",
  "IMs" = "SELENOP",
  "IMs" = "F13A1",
  "IMs" = "SLC40A1",
  "IMs" = "SPP1",
  "IMs" = "PLA2G7",
  "IMs" = "MMP9",
  "IMs" = "HAMP",
  # "recMac" = "S100A8",
  "recMac" = "FN1",
  "recMac" = "VCAN",
  "recMac" = "RNASE1",
  "recMac" = "FCN1",
  "recMac" = "EMP1",
  "recMac" = "CD14"
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
DotPlot(Macrophage, group.by = "res_1.0", features = Macrophage_major_marker_ls)+RotatedAxis()
ggsave(file.path(Macrophage_out_dir, 'dotplot_all_data_macro_cluster.pdf'), width = 15, height = 8)
DotPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.2", features = Macrophage_major_marker_ls)+RotatedAxis()
ggsave(file.path(Macrophage_out_dir, 'dotplot_macro_sub_cluster_res0.2.pdf'), width = 15, height = 8)
DotPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()
ggsave(file.path(Macrophage_out_dir, 'dotplot_macro_sub_cluster_res0.4.pdf'), width = 15, height = 8)

DotPlot(Macrophage, group.by = "Spatial.008um_snn_res.0.4", features = Major_marker_ls)+RotatedAxis()
ggsave(file.path(Macrophage_out_dir, 'dotplot_macro_sub_cluster_res0.4_maintype.pdf'), width = 25, height = 8)


saveRDS(Macrophage, file.path(Macrophage_out_dir, 'Macrophage.rds'))
