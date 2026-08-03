library(Seurat)
library(tidyverse)

all_data <- readRDS("~/R_OUTPUT/all_data_4slides_sample_filtered.rds")
all_data$annotation %>% table()
all_data$sample %>% table()
# all_data2 <- all_data[all_data$annotation != 'necrosis-area tumor cell', ]
# all_data2$annotation %>% table()
# all_data2$sample %>% table()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

marker_ls <- c(
  'B cell' = 'CD19', 
  'B cell' = 'CD79A', 
  'B cell' = 'MS4A1', 
  
  'Endothelial' = 'PECAM1', 
  'Endothelial' = 'PLVAP', 
  
  'Epithelial' = 'EPCAM', 
  
  'Fibroblast' = 'COL5A1', 
  'Fibroblast' = 'COL5A2', 
  'Fibroblast' = 'POSTN', 
  'Fibroblast' = 'PDGFRA', 
  'Fibroblast' = 'CCDC80', 
  'Fibroblast' = 'SFRP4', 
  
  'Macrophage' = 'CD68', 
  'Macrophage' = 'CD163',
  'Macrophage' = 'CSF1R', 
  # 'Macrophage' = 'IL1B', 
  
  'Neutrophil' = 'CSF3R', 
  
  'Pericyte' = 'PDGFRB', 
  'Pericyte' = 'RGS5', 
  'Pericyte' = 'MCAM', 
  
  'Plasma' = 'XBP1', 
  'Plasma' = 'MZB1', 
  
  'Skeletal muscle cell' = 'S100A1', 
  'Skeletal muscle cell' = 'RYR1', 
  'Skeletal muscle cell' = 'MYH1', 
  
  'T cell' = 'CD3E', 
  'T cell' = 'CD3G'
)

DotPlot(all_data, group.by = 'annotation', features = marker_ls, dot.scale = 10) + theme_90angle
ggsave('Dotplot_annotation.pdf', width = 15, height = 6)




# re annotation -------------------------------------------------------------------
out_dir <- "annotation_4slide" 
dir.create(out_dir)
table(all_data$group, all_data$annotation)

## PCA --------
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')
all_data <- RunPCA(all_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 20

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

all_data <- RunUMAP(all_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE) + NoLegend()
dev.off()

all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()




reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
all_data  <-
  FindNeighbors(
    all_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


all_data <- 
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


saveRDS(all_data, "~/R_OUTPUT/all_data_4slides_sample_filtered_reUMAP.rds")