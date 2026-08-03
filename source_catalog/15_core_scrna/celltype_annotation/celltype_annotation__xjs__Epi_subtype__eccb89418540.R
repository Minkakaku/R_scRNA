library(Seurat)
library(tidyverse)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")

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


Epithelial <- readRDS("~/CN_analysis_epi/all_data_CN_epis.rds")
Epithelial$annotation %>% table()
Epithelial$CN2 %>% table()




## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'Epithelial_subtype' 
dir.create(Epithelial_out_dir)


### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb


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

DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = "harmonyumap", cols = allcolour, raster=FALSE)

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
    k.param = 10
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# CN - res prop -----------
p_col_proj <- 
  Epithelial[[]] %>% 
  ungroup() %>% 
  group_by(CN2, RNA_snn_res.0.4) %>% 
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(RNA_snn_res.0.4 = factor(RNA_snn_res.0.4, levels = seq(0, 30))) %>% 
  mutate(CN2 = factor(CN2, levels = c(
    "Microenvironmrnt-poor tumor", 
    "Angio-associated tumor", 
    "Myeloid infiltration tumor", 
    "Lymphoid infiltration tumor", 
    "CAF-associated tumor"
  ))) %>% 
  
  ggplot(aes(x = CN2, y = freq, fill = RNA_snn_res.0.4)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_CN2_RNA_snn_res.0.4.pdf"),
  height = 6,
  width  = 6
)


saveRDS(Epithelial, 'R_OUTPUT/Epithelial_CN_reUMAP.rds')



