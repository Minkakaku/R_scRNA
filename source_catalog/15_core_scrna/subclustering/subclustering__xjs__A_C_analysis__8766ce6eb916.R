library(Seurat)
library(tidyverse)
source('cellchat_fix.R')

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#FA8072",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")

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

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

out_dir <- 'A_C_analysis'
dir.create(out_dir)

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$group %>% table()
all_data$group2 %>% table()
all_data$annotation %>% table()


# Epithelial annotation ----------
Epithelial_out_dir <- file.path(out_dir, 'Epithelial_annotation')
dir.create(Epithelial_out_dir)
Epithelial <- subset(all_data, annotation == 'Epithelial' & sample != 'D2' & group2 %in% c('A', 'C'))
Epithelial$sample %>% table()
Epithelial$group2 %>% table()
Epithelial$annotation %>% table()


## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)


### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group2.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group2", label = TRUE) + NoLegend()
dev.off()


Epithelial <- Epithelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample", 
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony" 
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

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 7)
DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group2 - harmony.pdf")), height = 6, width = 7.5)
DimPlot(Epithelial, group.by = "group2", label = TRUE, reduction = "harmonyumap")
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
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 7)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




# https://zhuanlan.zhihu.com/p/536775605
Epithelial_major_marker_ls <- c(
  'Basal cells/Myoepithelial' = "KRT5", # 
  'Basal cells/Myoepithelial' = "ACTA2", # 
  'Basal cells/Myoepithelial' = "KRT14", # 
  'Basal cells/Myoepithelial' = "KRT17 ", # 
  'Basal cells/Myoepithelial' = "MYL19", # 
  'Basal cells/Myoepithelial' = "MYLK", # 
  'Basal cells/Myoepithelial' = "MYH11", # 
  
  'Luminal cells' = 'KRT8', 
  'Luminal cells' = 'KRT18', 
  'Luminal cells' = 'KRT19', 
  
  'Luminal-HS cells' = "PRLR", 
  'Luminal-HS cells' = "CITED1", 
  'Luminal-HS cells' = "PGR", 
  'Luminal-HS cells' = "PROM1", 
  'Luminal-HS cells' = "ESR1", 
  
  'Luminal-AV cells' = "MFGE8", 
  'Luminal-AV cells' = "TRF", 
  'Luminal-AV cells' = "CSN3", 
  'Luminal-AV cells' = "WFDC18", 
  'Luminal-AV cells' = "ELF5", 
  'Luminal-AV cells' = "LTF", 
  
  'Luminal progenitor' = "KIT", 
  'Luminal progenitor' = "ALDH1A3", 
  'Luminal progenitor' = "CD14", 
  
  'Tumor cell' = 'MDM2', 
  
  'HER2-tumor cell' = 'ERBB2'
)

# find markers
Idents(Epithelial) <- "RNA_snn_res.1.2"
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
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 8)
DotPlot(Epithelial, group.by = "RNA_snn_res.1.2", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.1.2", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.1.2", features = Epithelial_major_marker_ls2)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.1.2", features = Epithelial_top5_markers)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', group.by = "RNA_snn_res.0.6", pt.size = 0)
DotPlot(Epithelial, group.by = "RNA_snn_res.1.5", features = major_marker_ls)+RotatedAxis()

# pdf(file.path(Epithelial_out_dir, 'dimplot_18_19_21.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), res_1.2 %in% c(18, 19, 21)))) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir, 'dimplot_C3_D2.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), sample %in% c('C3', 'D2')))) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir, 'dimplot_D2.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), sample %in% c('D2')))) + NoLegend()
# dev.off()





## annotation -----------
Epithelial$RNA_snn_res.1.2 %>% table()
table(Epithelial$RNA_snn_res.1.2, Epithelial$sample)

Epithelial[[]] %>%
  group_by(RNA_snn_res.1.2, sample) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.1.2, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(sample, RNA_snn_res.1.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.1.2, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample.pdf'), width = 5, height = 6)

# Epithelial[[]] %>%
#   group_by(sample, group2, RNA_snn_res.1.2) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = RNA_snn_res.1.2, y = freq, fill = group2)) + 
#   geom_boxplot() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# ggsave(file.path(Epithelial_out_dir, 'boxplot_cluster prop by sample.pdf'), width = 8, height = 6)


Epithelial[[]] %>%
  group_by(RNA_snn_res.1.2, group2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.1.2, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(group2, RNA_snn_res.1.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.1.2, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2.pdf'), width = 5, height = 6)




cell_cluster <- c(
  '0' = "Epi_0", 
  '1' = "Epi_1", 
  '2' = "Epi_2", 
  '3' = "Epi_3", 
  '4' = "Epi_4", 
  '5' = "Epi_5", 
  '6' = "Epi_6", 
  '7' = "Epi_7", 
  '8' = "Epi_8", 
  '9' = "Epi_9", 
  '10' = "noise", 
  '11' = "mix", 
  '12' = "Epi_12", 
  '13' = "Epi_13", 
  '14' = "noise", 
  '15' = "Epi_15", 
  '16' = "Epi_16", 
  '17' = "Epi_17", 
  '18' = "Epi_18", # C3
  '19' = "Epi_19", 
  '20' = "Epi_20", 
  '21' = "mix", 
  '22' = "Epi_22", # C2
  '23' = "Epi_23", 
  '24' = "noise", 
  '25' = "noise"
)


Idents(Epithelial) <- "RNA_snn_res.1.2"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Epithelial$subtype[str_starts(Epithelial$subtype, 'Epi')] %>% unique()
Epithelial_pure <- subset(JoinLayers(Epithelial), subtype %in% (Epithelial$subtype[str_starts(Epithelial$subtype, 'Epi')] %>% unique()))
Epithelial_pure <- FindVariableFeatures(Epithelial_pure)
Epithelial_pure <- ScaleData(Epithelial_pure, vars.to.regress = 'nCount_RNA')
Epithelial_pure <- RunPCA(Epithelial_pure)

pdf(file.path(Epithelial_out_dir, 'pca plot - pure.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial_pure, ndims = 50)
dev.off()

dims <- 15

Epithelial_pure <- RunUMAP(Epithelial_pure, dims = 1:dims, return.model = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

Epithelial_pure <- Epithelial_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Epithelial_pure <- RunUMAP(
  Epithelial_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

DimPlot(Epithelial_pure, group.by = "sample", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 10)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split - pure.pdf")), height = 7, width = 16)
DimPlot(Epithelial_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial.rds'))




Epithelial <- Epithelial_pure
Epithelial_out_dir <- file.path(Epithelial_out_dir, 'filter')
dir.create(Epithelial_out_dir)
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
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 7)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




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
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 8)
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls2)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_top5_markers)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', group.by = "RNA_snn_res.0.6", pt.size = 0)
DotPlot(Epithelial, group.by = "RNA_snn_res.1.5", features = major_marker_ls)+RotatedAxis()

# pdf(file.path(Epithelial_out_dir, 'dimplot_18_19_21.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), res_1.2 %in% c(18, 19, 21)))) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir, 'dimplot_C3_D2.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), sample %in% c('C3', 'D2')))) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir, 'dimplot_D2.pdf'), width = 8, height = 6)
# DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = 'harmonyumap', cells.highlight = Cells(subset(JoinLayers(Epithelial), sample %in% c('D2')))) + NoLegend()
# dev.off()





## annotation -----------
Epithelial$RNA_snn_res.0.6 %>% table()
table(Epithelial$RNA_snn_res.0.6, Epithelial$sample)

Epithelial[[]] %>% 
  group_by(RNA_snn_res.0.6, sample) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(sample, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = sample)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by sample.pdf'), width = 5, height = 6)

# Epithelial[[]] %>%
#   group_by(sample, group2, RNA_snn_res.0.6) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
#   geom_boxplot() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# ggsave(file.path(Epithelial_out_dir, 'boxplot_cluster prop by sample.pdf'), width = 8, height = 6)


Epithelial[[]] %>% 
  group_by(RNA_snn_res.0.6, group2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2_scale.pdf'), width = 5, height = 6)

Epithelial[[]] %>%
  group_by(group2, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.6, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Epithelial_out_dir, 'cluster prop by group2.pdf'), width = 5, height = 6)




cell_cluster <- c(
  '0' = "Epi_0", 
  '1' = "Epi_1", 
  '2' = "Epi_2", 
  '3' = "low-nCounts", 
  '4' = "Epi_4", 
  '5' = "Epi_5", 
  '6' = "Epi_6", 
  '7' = "Epi_7", 
  '8' = "Epi_8", 
  '9' = "Epi_9", 
  '10' = "mix", 
  '11' = "Epi_11", 
  '12' = "low-nCounts", 
  '13' = "mix", 
  '14' = "mix"
)


Idents(Epithelial) <- "RNA_snn_res.0.6"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 11)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Epithelial$subtype[str_starts(Epithelial$subtype, 'Epi')] %>% unique()
Epithelial_pure <- subset(JoinLayers(Epithelial), subtype %in% (Epithelial$subtype[str_starts(Epithelial$subtype, 'Epi')] %>% unique()))
Epithelial_pure <- FindVariableFeatures(Epithelial_pure)
Epithelial_pure <- ScaleData(Epithelial_pure, vars.to.regress = 'nCount_RNA')
Epithelial_pure <- RunPCA(Epithelial_pure)

pdf(file.path(Epithelial_out_dir, 'pca plot - pure.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial_pure, ndims = 50)
dev.off()

dims <- 15

Epithelial_pure <- RunUMAP(Epithelial_pure, dims = 1:dims, return.model = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

Epithelial_pure <- Epithelial_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    theta = 3, # default 2, bigger value -- more diversity per cluster
    lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Epithelial_pure <- RunUMAP(
  Epithelial_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

DimPlot(Epithelial_pure, group.by = "sample", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(Epithelial_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split - pure.pdf")), height = 6, width = 12)
DimPlot(Epithelial_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial.rds'))
saveRDS(Epithelial_pure, file.path(Epithelial_out_dir, 'Epithelial_pure.rds'))




# Epi 0 2 11 DEGs & enrichment ------------
Epithelial_pure$subtype %>% table()
Idents(Epithelial_pure) <- 'subtype'

Epithelial_pure_markers <- 
  FindAllMarkers(
    JoinLayers(Epithelial_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_pure_markers,file.path(Epithelial_out_dir, 'Epithelial_pure_DEGs_subtype.csv'))

library(EnhancedVolcano)
type<-unique(Epithelial_pure_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_pure_markers[Epithelial_pure_markers$cluster == type, ]
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
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Epithelial_pure = Epithelial_pure_markers)
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
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })


Epithelial_pure_de_genes_for_enrichmenet <- 
  Epithelial_pure_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1-pct.2 > 0.05) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  Epithelial_pure_de_genes_for_enrichmenet %>%
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


Epithelial_pure_er_output <- file.path(Epithelial_out_dir, 'enrichment') 
dir.create(Epithelial_pure_er_output)
# expot enrichment table
enrich_ls %>% 
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(Epithelial_pure_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Epithelial_pure_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_pure_de_genes_id_for_enrichmenet <-
  Epithelial_pure_de_genes_for_enrichmenet %>%
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
  Epithelial_pure_de_genes_id_for_enrichmenet %>%
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

# saveRDS(enrich_ls,file.path(Epithelial_pure_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_pure_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_pure_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Epithelial_pure_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Epithelial_pure_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })





# Fibroblast plot ---------
Fibroblast_out_dir <- file.path(out_dir, 'Fibroblast_annotation')
dir.create(Fibroblast_out_dir)
Fibroblast <- readRDS("Fibroblast_subtype_remove_D2/Fibroblast.rds")
Fibroblast_A_C <- subset(Fibroblast, group2 %in% c('A', 'C'))
Fibroblast_A_C$subtype %>% table()
Fibroblast <- subset(Fibroblast_A_C, subtype %in% c('LENG8-mCAF', 'THBS1-mCAF', 'iCAF', 'tCAF', 'ifnCAF', 'vCAF'))

Fibroblast[[]] %>% 
  group_by(group2, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Fibroblast_out_dir, 'cluster prop by group2.pdf'), width = 5, height = 6)

Fibroblast[[]] %>% 
  group_by(subtype, group2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq, fill = group2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(file.path(Fibroblast_out_dir, 'cluster prop by group2_scaled.pdf'), width = 5, height = 6)

Fibroblast$group2 %>% table()




# THBS1_mCAF DEGs -------------
Fibroblast$subtype %>% table()
Idents(Fibroblast) <- 'subtype'

THBS1_mCAF_markers <-  
  FindMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox", 
    ident.1 = 'THBS1-mCAF', 
    # logfc.threshold = 0.25,
    # min.pct = 0.1,
    only.pos = TRUE
  ) %>% 
  dplyr::filter(p_val_adj < 0.05)
write.csv(THBS1_mCAF_markers, file.path(Fibroblast_out_dir, 'THBS1_mCAF_DEGs.csv'))

THBS1_mCAF_markers %>% rownames()




# 映射------------
all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$group2 %>% table()
all_data$annotation %>% table()
all_data_A_C <- subset(all_data, group2 %in% c('A', 'C'))
all_data_A_C$group2 %>% table()

Fibroblast <- readRDS("Fibroblast_subtype_remove_D2/Fibroblast.rds")
Fibroblast_A_C <- subset(Fibroblast, group2 %in% c('A', 'C'))
Fibroblast_A_C$subtype %>% table()

Epithelial_pure <- readRDS("A_C_analysis/Epithelial_annotation/filter/Epithelial_pure.rds")
Epithelial <- readRDS("A_C_analysis/Epithelial_annotation/Epithelial.rds")
Epithelial_pure$subtype <- as.character(Epithelial_pure$subtype)
Epithelial$annotation <- as.character(Epithelial$annotation)

subtype <- 
  Epithelial_pure[['subtype']] %>% 
  rownames_to_column('cells')
subtype$subtype %>% table()
Epithelial$annotation %>% table()

meta.data.t <- 
  subtype %>% 
  right_join(
    Epithelial[['annotation']] %>% rownames_to_column('cells'), 
    by = 'cells'
    ) %>% 
  mutate(
    subtype_all = ifelse(
      is.na(subtype), 
      annotation, 
      subtype
    )
  ) %>% 
  dplyr::select(cells, subtype_all) %>% 
  column_to_rownames('cells')

meta.data.t$subtype_all <- factor(meta.data.t$subtype_all, levels = meta.data.t$subtype_all %>% unique())
meta.data.t$subtype_all %>% table()

Epithelial <- AddMetaData(Epithelial, meta.data.t['subtype_all'])
Epithelial$subtype_all %>% table()
Epithelial$group2 %>% table()

Fibroblast_A_C$subtype_all <- Fibroblast_A_C$subtype

subtype_all <- rbind(
  Fibroblast_A_C[['subtype_all']] %>% rownames_to_column('cells'), 
  Epithelial[['subtype_all']] %>% rownames_to_column('cells')
)
all_data_A_C$annotation %>% table()

meta.data.t <- 
  subtype_all %>% 
  right_join(
    all_data_A_C[['annotation']] %>% rownames_to_column('cells'), 
    by = 'cells'
  ) %>% 
  mutate(
    subtype_all = ifelse(
      (as.character(annotation) %in% c('Epithelial', 'Fibroblast')), 
      as.character(subtype_all), 
      as.character(annotation)
    )
  ) %>% 
  dplyr::select(cells, subtype_all) %>% 
  column_to_rownames('cells')
meta.data.t$subtype_all %>% table()

all_data_A_C<-AddMetaData(all_data_A_C,meta.data.t['subtype_all'])
all_data_A_C$subtype_all %>% table()
all_data_A_C_pure <- subset(all_data_A_C, subtype_all != 'low-nCounts' & subtype_all != 'mix' & subtype_all != 'Epithelial')
all_data_A_C_pure$subtype_all %>% table()
saveRDS(all_data_A_C_pure, file.path(Epithelial_out_dir, 'all_data_A_C_pure.rds'))


## Annotation info output -------------------------------------------------------
dir.create('anno_info')
all_data_A_C_pure$annotation %>% table()
cellType_anno_ls <- 
  all_data_A_C_pure@meta.data[c('project', 'annotation')] %>% 
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
    write.csv(cellType, str_glue('anno_info/AC_cell type {project}.csv'), row.names = F)
  })

all_data_A_C_pure$subtype_all %>% table()
subType_anno_ls <- 
  all_data_A_C_pure@meta.data[c('project', 'subtype_all')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = subtype_all
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

subType_anno_ls %>% 
  iwalk(function(subType, project){
    write.csv(subType, str_glue('anno_info/AC_sub cell type {project}.csv'), row.names = F)
  })




# THBS1 expr --------------
## spatial expr ----------
### target gene ------------
target_gene_ls <- c(
  'THBS1'
)
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data_A_C_pure)


all_data_A_C_pure_meta <- 
  all_data_A_C_pure[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data_A_C_pure %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_A_C_pure_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_A_C_pure_meta, data_expr, by = 'cell_id')

p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        if (!sample %in% unique(all_data_A_C_pure_meta$sample)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], group = cell_id)
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          # scale_color_manual(values = c("lightgrey", allcolour)) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.1,
              barheight = 10
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Epithelial_out_dir,  
                  str_glue("spatial dim plot {grp}_middle0.1.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 




### target gene + cell type ------------
p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        if (!sample %in% unique(all_data_A_C_pure_meta$sample)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = subtype_all, group = cell_id), 
            size = 0.05,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", allcolour)) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.1,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Epithelial_out_dir,  
                  str_glue("spatial dim plot {grp}_subtype_middle0.1.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 




# THBS1 SDC4 CD47 expr --------------
VlnPlot(all_data_A_C_pure, features = c('THBS1', 'SDC4', 'CD47'), group.by = 'subtype_all', pt.size = 0, ncol = 1) + NoLegend()
ggsave(file.path(Epithelial_out_dir, 'Vlnplot_THBS1-SDC4-CD47_expr.pdf'), width = 8, height = 15)


## spatial expr ----------
### target gene + celltype ------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

target_gene_ls <- c(
  'THBS1', 'SDC4', 'CD47'
)
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data)


all_data_A_C_pure$subtype_all %>% table()
all_data_A_C_pure$annotation1 <- fct_collapse(
  all_data_A_C_pure$subtype_all, 
  'Others' = unique(all_data_A_C_pure$subtype_all)[!unique(all_data_A_C_pure$subtype_all) %in% c('Epi_0', 'Epi_2', 'Epi_11', 'THBS1-mCAF')]
)
all_data_A_C_pure$annotation1 %>% table()

all_data_A_C_pure_meta <- 
  all_data_A_C_pure[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data_A_C_pure %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_A_C_pure_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_A_C_pure_meta, data_expr, by = 'cell_id')


p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        if (!sample %in% unique(all_data_A_C_pure_meta$sample)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = annotation1, group = cell_id), 
            size = 0.1,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", "#40E0D0","#FA8072","#FFFF00","#0000FF")) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.1,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Epithelial_out_dir,  
                  str_glue("spatial dim plot {grp}_target gene - Epi0 - Epi2 - Epi11 - THBS1_mCAF_highlight_middle0.1.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 




# SP1 expr --------------
VlnPlot(all_data_A_C_pure, features = c('SP1'), group.by = 'subtype_all', pt.size = 0, ncol = 1) + NoLegend()
ggsave(file.path(Epithelial_out_dir, 'Vlnplot_SP1_expr.pdf'), width = 8, height = 5)


## spatial expr ----------
### target gene + celltype ------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

target_gene_ls <- c(
  'SP1'
)
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data)


all_data_A_C_pure$subtype_all %>% table()
all_data_A_C_pure$annotation1 <- fct_collapse(
  all_data_A_C_pure$subtype_all, 
  'Others' = unique(all_data_A_C_pure$subtype_all)[!unique(all_data_A_C_pure$subtype_all) %in% c('Epi_0', 'Epi_2', 'Epi_11', 'THBS1-mCAF')]
)
all_data_A_C_pure$annotation1 %>% table()

all_data_A_C_pure_meta <- 
  all_data_A_C_pure[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data_A_C_pure %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_A_C_pure_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_A_C_pure_meta, data_expr, by = 'cell_id')


p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        if (!sample %in% unique(all_data_A_C_pure_meta$sample)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = annotation1, group = cell_id), 
            size = 0.1,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", "#40E0D0","#FA8072","#FFFF00","#0000FF")) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.1,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(Epithelial_out_dir,  
                  str_glue("spatial dim plot {grp}_target gene - Epi0 - Epi2 - Epi11 - THBS1_mCAF_highlight_middle0.1.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 








# cellchat --------------
library(CellChat)

## define default-------------
out_dir_cellchat <- "A_C_analysis" 
out_dir_cellchat <- file.path(out_dir_cellchat, 'cellchat_result')
dir.create(out_dir_cellchat)

all_data <- all_data_A_C_pure
all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

all_data_meta <- all_data_meta_all
all_data_norm <- all_data_norm_all

all_grp_ls <-
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, group2) %>%
  group_by(group2) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.name = FALSE)




## using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    # table
    cells <- all_data_meta_all$subtype_all %>% unique()
    if(!all(cells %in% meta$subtype_all)){
      missing_cells <- cells %>% setdiff(unique(meta$subtype_all))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(subtype_all %in% missing_cells) %>%
        group_by(subtype_all) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$group2)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "subtype_all")
    
    # cellchat db
    CellChatDB.use <- 
      CellChatDB.human %>% 
      subsetDB(
        search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), 
        key    = "annotation"
      )
    
    cellchat@DB <- CellChatDB.use
    
    
    # preprocess
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # optinal - project PPI data
    # cellchat = projectData(cellchat, PPI.human)
    cellchat <- smoothData(cellchat, adj = PPI.human)
    
    # Compute the communication probability
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    # cellchat <- netAnalysis_computeCentrality(cellchat)
    
    cellchat
  })


saveRDS(all_cellchat_ls, file.path(out_dir_cellchat, "all_cellchat_ls2.rds"))




# export ligand receptor table
all_cellchat_ls  %>%
  iwalk(function(cellchat, name){ 
    df.net <- tryCatch(
      cellchat %>% subsetCommunication() %>% as.data.frame(),
      error = function(e) NULL
    )
    if(!is.null(df.net)){
      df.net$relation <-  paste0(df.net[[1]], "->", df.net[[2]])
      write.table(
        df.net,
        file      = file.path(out_dir_cellchat, str_glue("net_lr - {name}.txt")),
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
    
    # network
    df.netp <- tryCatch(
      subsetCommunication(cellchat, slot.name="netP"),
      error = function(e) NULL
    )
    if(!is.null(df.netp)){
      write.table(
        df.netp,
        file      = file.path(out_dir_cellchat, str_glue("netPathways_lr - {name}.txt")),
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
  })



# merge cellchat
idents.list <- map(all_cellchat_ls, ~ levels(.x@idents))
idents.list.len <- idents.list %>% map_dbl(length)

if(any(idents.list.len[1] != idents.list.len)){
  
  idents.max      <- idents.list[[which.max(idents.list.len)]]
  all_cellchat_ls <- all_cellchat_ls %>% map(function(x) liftCellChat(x, idents.max))
}



all_cellchat_new_ls <-
  all_cellchat_ls %>%
  imap(function(cellchat, name){ # browser()
    cellchat <- tryCatch(netAnalysis_computeCentrality(cellchat), error = function(e) NULL)
    if(is.null(cellchat)){
      cat(str_glue("group {name} cannot run netAnalysis_computeCentrality and removed!\n\n"))
    }
    cellchat
  }) %>% 
  purrr::discard(is.null)  

groups   <- names(all_cellchat_new_ls)
cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)


saveRDS(all_cellchat_new_ls, file.path(out_dir_cellchat, "all_cellchat_new_ls.rds")) 
saveRDS(cellchat, file.path(out_dir_cellchat, "cellchat.rds")) 


all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)


## plot ---------------
### comprare interation ------
pdf(
  file.path(out_dir_cellchat, str_glue("comprare interation plot.pdf")),
  height = 4,
  width  = 6
)
compareInteractions(cellchat, show.legend = FALSE, measure = "count")
compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
dev.off()


### diff network plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("diff network plot - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    dev.off()
  })




# THBS1-mCAF
all_combns_ls %>% 
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("THBS1-mCAF - diff network plot - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      targets.use = 'THBS1-mCAF', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      targets.use = 'THBS1-mCAF', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      sources.use = 'THBS1-mCAF', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      sources.use = 'THBS1-mCAF', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    dev.off()
  })


### diff heatmap plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    
    pdf(
      file.path(out_dir_cellchat, str_glue("diff heatmap - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    print(netVisual_heatmap(cellchat, comparison = x))
    print(netVisual_heatmap(cellchat, comparison = x, measure = "weight"))
    dev.off()
  })


### diff rank plot --------
all_combns_ls %>% 
  walk(function(x){
    x <- c(x[2], x[1])
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked     = T,
      do.stat     = TRUE
    ) +
      ggtitle("Relative information flow - count")
    print(p1)
    
    p2 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked       = T,
      do.stat       = TRUE
    ) +
      ggtitle("Relative information flow - weight")
    print(p2)
    
    p3 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked     = F,
      do.stat     = TRUE
    ) +
      ggtitle("Information flow - count")
    print(p3)
    
    p4 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p4)
    
    dev.off()
  })


# THBS1-mCAF-source
all_combns_ls %>% 
  walk(function(x){
    x <- c(x[2], x[1])
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("THBS1-mCAF-source diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = 'THBS1-mCAF',
      targets.use = NULL, 
      comparison  = x,
      stacked     = T,
      do.stat     = TRUE
    ) +
      ggtitle("Relative information flow - count")
    print(p1)
    
    p2 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = 'THBS1-mCAF',
      targets.use = NULL, 
      comparison  = x,
      stacked       = T,
      do.stat       = TRUE
    ) +
      ggtitle("Relative information flow - weight")
    print(p2)
    
    p3 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = 'THBS1-mCAF',
      targets.use = NULL, 
      comparison  = x,
      stacked     = F,
      do.stat     = TRUE
    ) +
      ggtitle("Information flow - count")
    print(p3)
    
    p4 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = 'THBS1-mCAF',
      targets.use = NULL, 
      comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p4)
    
    dev.off()
  })

# THBS1-mCAF-target
all_combns_ls %>% 
  walk(function(x){
    x <- c(x[2], x[1])
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("THBS1-mCAF-target diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = 'THBS1-mCAF', 
      comparison  = x,
      stacked     = T,
      do.stat     = TRUE
    ) +
      ggtitle("Relative information flow - count")
    print(p1)
    
    p2 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = 'THBS1-mCAF', 
      comparison  = x,
      stacked       = T,
      do.stat       = TRUE
    ) +
      ggtitle("Relative information flow - weight")
    print(p2)
    
    p3 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = 'THBS1-mCAF', 
      comparison  = x,
      stacked     = F,
      do.stat     = TRUE
    ) +
      ggtitle("Information flow - count")
    print(p3)
    
    p4 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = 'THBS1-mCAF', 
      comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p4)
    
    dev.off()
  })




### diff incoming and outgoing cells --------
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - outgoing.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){ # browser()
    pathways <- cellchat@netP$pathways
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "outgoing", 
        signaling = pathways,
        title     = name,
        width     = 10,
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

# incoming 
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - incoming.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    
    tryCatch({
      ## YOUR CODE
      p <- cellchat %>%
        netAnalysis_computeCentrality() %>%
        netAnalysis_signalingRole_heatmap(
          pattern   = "incoming", 
          signaling = pathways,
          title     = name,
          width     = 10, 
          height    = length(pathways)/3
        )
      print(p)
    },
    error = function(x) {
      print(x$message)
      NULL
    }
    )
    
  })
dev.off()

# all 
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - all.pdf")),
  height = 10, 
  width  = 6
)
all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "all", 
        signaling = pathways,
        title     = name,
        width     = 10, 
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()


### diff interaction strength -----
pdf(
  file.path(out_dir_cellchat, str_glue("diff interaction strength - all groups.pdf")),
  height = 5,
  width  = 15
)
num.link <- sapply(all_cellchat_new_ls, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
all_cellchat_new_ls %>%
  imap(function(cellchat, name){
    cellchat %>% 
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_scatter(
        title = name,
        weight.MinMax = range(num.link)
      )
  }) %>%
  patchwork::wrap_plots()
dev.off()




### THBS1-mCAF - C diff incoming and outgoing cells --------
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("THBS1-mCAF - C_epi contribution of signals to cell groups - outgoing.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){ # browser()
    pathways <- c('CXCL', 'COLLAGEN', 'LAMININ', 'THBS', 'MPZ', 'HSPG', 'CD99')
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "outgoing", 
        signaling = pathways,
        title     = name,
        width     = 10,
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

# incoming 
pdf(
  file.path(out_dir_cellchat, str_glue("THBS1-mCAF - C_epi contribution of signals to cell groups - incoming.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- c('CXCL', 'COLLAGEN', 'LAMININ', 'THBS', 'MPZ', 'HSPG', 'CD99')
    
    tryCatch({
      ## YOUR CODE
      p <- cellchat %>%
        netAnalysis_computeCentrality() %>%
        netAnalysis_signalingRole_heatmap(
          pattern   = "incoming", 
          signaling = pathways,
          title     = name,
          width     = 10, 
          height    = length(pathways)/3
        )
      print(p)
    },
    error = function(x) {
      print(x$message)
      NULL
    }
    )
    
  })
dev.off()




### THBS1-mCAF - all diff incoming and outgoing cells --------
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("THBS1-mCAF - all_epi contribution of signals to cell groups - outgoing.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){ # browser()
    pathways <- c('CXCL', 'COLLAGEN', 'LAMININ', 'THBS', 'MPZ', 'HSPG', 'CD99', 'GRN')
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "outgoing", 
        signaling = pathways,
        title     = name,
        width     = 10,
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

# incoming 
pdf(
  file.path(out_dir_cellchat, str_glue("THBS1-mCAF - all_epi contribution of signals to cell groups - incoming.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- c('CXCL', 'COLLAGEN', 'LAMININ', 'THBS', 'MPZ', 'HSPG', 'CD99', 'GRN')
    
    tryCatch({
      ## YOUR CODE
      p <- cellchat %>%
        netAnalysis_computeCentrality() %>%
        netAnalysis_signalingRole_heatmap(
          pattern   = "incoming", 
          signaling = pathways,
          title     = name,
          width     = 10, 
          height    = length(pathways)/3
        )
      print(p)
    },
    error = function(x) {
      print(x$message)
      NULL
    }
    )
    
  })
dev.off()




### chordDiagram 4 `THBS1-mCAF 2 Epi` -----------
# df_ct <-
#   data.frame(
#     source = rep(LETTERS[1:3], 4),
#     target = rep(letters[1:6], 2),
#     prob = 1:12,
#     pathway = c(
#       'CXCL',
#       'CXCL',
#       'COLLAGEN',
#       'COLLAGEN',
#       'THBS',
#       'THBS',
#       'THBS',
#       'THBS',
#       'THBS',
#       'MPZ',
#       'MPZ',
#       'MPZ'
#     )
#   )

#### A -----------------
##### 3Epi -----------
df_ct <- 
  read.table('A_C_analysis/cellchat_result/net_lr - A.txt', header = T, sep = '\t') %>% 
  dplyr::filter(source == 'THBS1-mCAF') %>% 
  dplyr::filter(target %in% c('Epi_0', 'Epi_2', 'Epi_11')) %>% 
  select(source, target, prob, interaction_name_2)

# color alpha
alpha_val <-  0.5

# pathway and pathway color
col_map <- structure(
  # ggsci::pal_igv()(length(unique(df_ct$interaction_name_2))), 
  allcolour[1:length(unique(df_ct$interaction_name_2))],
  names = unique(df_ct$interaction_name_2)
)

# linke color
link.colors <-
  df_ct$interaction_name_2 %>%
  fct_recode(!!! structure(names(col_map), names = col_map)) %>%
  alpha(alpha = alpha_val)

# plot
pdf(file.path(out_dir_cellchat, 'chordDiagram - THBS1-mCAF 2 3Epi - A.pdf'), height = 8, width = 12)
# par(mar = c(5, 4, 4, 2))
par(mar = c(2, 2, 2, 10))
circlize::circos.par()

circlize::chordDiagram(
  df_ct,
  # linke color
  col = link.colors,
  
  # transparency
  transparency = 1-alpha_val,
  
  # link is a big arrow
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

# add legend
legend(
  # 'topright',
  x = 1.2,
  y = 0.5,
  legend = names(col_map),
  fill = alpha(col_map, alpha_val),
  bty = "o",
  xpd = TRUE
)

circlize::circos.clear()
dev.off()


##### Epi_0 -----------
df_ct <- 
  read.table('A_C_analysis/cellchat_result/net_lr - A.txt', header = T, sep = '\t') %>% 
  dplyr::filter(source == 'THBS1-mCAF') %>% 
  dplyr::filter(target %in% c('Epi_0')) %>% 
  select(source, target, prob, interaction_name_2)

# color alpha
alpha_val <-  0.5

# pathway and pathway color
col_map <- structure(
  # ggsci::pal_igv()(length(unique(df_ct$interaction_name_2))), 
  allcolour[1:length(unique(df_ct$interaction_name_2))],
  names = unique(df_ct$interaction_name_2)
)

# linke color
link.colors <-
  df_ct$interaction_name_2 %>%
  fct_recode(!!! structure(names(col_map), names = col_map)) %>%
  alpha(alpha = alpha_val)

# plot
pdf(file.path(out_dir_cellchat, 'chordDiagram - THBS1-mCAF 2 Epi_0 - A.pdf'), height = 8, width = 12)
# par(mar = c(5, 4, 4, 2))
par(mar = c(2, 2, 2, 10))
circlize::circos.par()

circlize::chordDiagram(
  df_ct,
  # linke color
  col = link.colors,
  
  # transparency
  transparency = 1-alpha_val,
  
  # link is a big arrow
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

# add legend
legend(
  # 'topright',
  x = 1.2,
  y = 0.5,
  legend = names(col_map),
  fill = alpha(col_map, alpha_val),
  bty = "o",
  xpd = TRUE
)

circlize::circos.clear()
dev.off()


#### C -----------------
##### 3Epi ----------
df_ct <- 
  read.table('A_C_analysis/cellchat_result/net_lr - C.txt', header = T, sep = '\t') %>% 
  dplyr::filter(source == 'THBS1-mCAF') %>% 
  dplyr::filter(target %in% c('Epi_0', 'Epi_2', 'Epi_11')) %>% 
  select(source, target, prob, interaction_name_2)

# color alpha
alpha_val <-  0.5

# pathway and pathway color
col_map <- structure(
  # ggsci::pal_igv()(length(unique(df_ct$interaction_name_2))), 
  allcolour[1:length(unique(df_ct$interaction_name_2))],
  names = unique(df_ct$interaction_name_2)
)

# linke color
link.colors <-
  df_ct$interaction_name_2 %>%
  fct_recode(!!! structure(names(col_map), names = col_map)) %>%
  alpha(alpha = alpha_val)

# plot
pdf(file.path(out_dir_cellchat, 'chordDiagram - THBS1-mCAF 2 3Epi - C.pdf'), height = 8, width = 12)
# par(mar = c(5, 4, 4, 2))
par(mar = c(2, 2, 2, 10))
circlize::circos.par()

circlize::chordDiagram(
  df_ct,
  # linke color
  col = link.colors,
  
  # transparency
  transparency = 1-alpha_val,
  
  # link is a big arrow
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

# add legend
legend(
  # 'topright',
  x = 1.2,
  y = 0.5,
  legend = names(col_map),
  fill = alpha(col_map, alpha_val),
  bty = "o",
  xpd = TRUE
)

circlize::circos.clear()
dev.off()


##### Epi_0 ----------
df_ct <- 
  read.table('A_C_analysis/cellchat_result/net_lr - C.txt', header = T, sep = '\t') %>% 
  dplyr::filter(source == 'THBS1-mCAF') %>% 
  dplyr::filter(target %in% c('Epi_0')) %>% 
  select(source, target, prob, interaction_name_2)

# color alpha
alpha_val <-  0.5

# pathway and pathway color
col_map <- structure(
  # ggsci::pal_igv()(length(unique(df_ct$interaction_name_2))), 
  allcolour[1:length(unique(df_ct$interaction_name_2))],
  names = unique(df_ct$interaction_name_2)
)

# linke color
link.colors <-
  df_ct$interaction_name_2 %>%
  fct_recode(!!! structure(names(col_map), names = col_map)) %>%
  alpha(alpha = alpha_val)

# plot
pdf(file.path(out_dir_cellchat, 'chordDiagram - THBS1-mCAF 2 Epi_0 - C.pdf'), height = 8, width = 12)
# par(mar = c(5, 4, 4, 2))
par(mar = c(2, 2, 2, 10))
circlize::circos.par()

circlize::chordDiagram(
  df_ct,
  # linke color
  col = link.colors,
  
  # transparency
  transparency = 1-alpha_val,
  
  # link is a big arrow
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

# add legend
legend(
  # 'topright',
  x = 1.2,
  y = 0.5,
  legend = names(col_map),
  fill = alpha(col_map, alpha_val),
  bty = "o",
  xpd = TRUE
)

circlize::circos.clear()
dev.off()




# paper figure ------------
out_dir <- 'A_C_analysis/paper_fig'
dir.create(out_dir)

## all_data ----------
all_data <- readRDS("all_data_annotation_removeD2_filtered.rds")
all_data$annotation_full %>% table()
all_data$group2 %>% table()
A_C <- subset(all_data, group2 %in% c('A', 'C'))


### dimplot ------------
DimPlot(A_C, cols = allcolour, label = T, alpha = 0.4, label.size = 3, group.by = 'annotation') + NoLegend()
ggsave(file.path(out_dir, 'umap_celltype_noLegend.pdf'), width = 6, height = 5.5)
DimPlot(A_C, cols = allcolour, label = F, alpha = 0.4, label.size = 3, group.by = 'annotation')
ggsave(file.path(out_dir, 'umap_celltype.pdf'), width = 7, height = 5.5)

DimPlot(A_C, cols = allcolour, label = T, alpha = 0.4, label.size = 3, split.by = 'group2', ncol = 2, group.by = 'annotation') + NoLegend()
ggsave(file.path(out_dir, 'umap_celltype_noLegend_split.pdf'), width = 10, height = 5.5)
DimPlot(A_C, cols = allcolour, label = F, alpha = 0.4, label.size = 3, split.by = 'group2', ncol = 2, group.by = 'annotation')
ggsave(file.path(out_dir, 'umap_celltype_split.pdf'), width = 11, height = 5.5)

### dotplot -------------
major_marker_ls <- c(
  "B cell" = "CD19",
  "B cell" = "CD79A",
  "B cell" = "MS4A1",
  
  'Endothelial' = 'PECAM1', 
  'Endothelial' = 'VWF', 
  
  'Epithelial'  = 'EPCAM',
  'Epithelial'  = 'KRT5',
  'Epithelial'  = 'KRT14',

  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80",
  "Fibroblast" = "SFRP4",
  "Fibroblast" = "COL5A1",
  "Fibroblast" = "COL5A2",
  
  'Myeloid cell' = 'ITGAM', 
  'Myeloid cell' = 'CD14', 
  'Myeloid cell' = 'CD68', 
  'Myeloid cell' = 'CD83', 
  
  "Pericyte" = "PDGFRB",
  "Pericyte" = "RGS5",
  "Pericyte" = "MCAM",
  
  "Plasma cell" = "TNFRSF17",
  "Plasma cell" = "MZB1",
  "Plasma cell" = "XBP1"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

A_C$annotation_re <- 
  str_replace_all(A_C$annotation, c(
    'Adipocyte' = 'Epithelial'
  ))
A_C$annotation_re %>% table()

DotPlot(A_C, features = major_marker_ls, group.by = 'annotation_re') + theme_90angle
ggsave(file.path(out_dir, 'dotplot_celltype.pdf'), width = 12, height = 5)


### col plot ---------------
A_C_meta <- 
  A_C@meta.data

p_col_proj <-
  A_C_meta %>%
  group_by(group2, annotation_re) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group2, y = freq, fill = annotation_re)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(out_dir, "cell_freq_group2.pdf"),
  height = 5,
  width  = 4
)




## Fibroblast ----------
Fibroblast <- readRDS("Fibroblast_subtype_remove_D2/Fibroblast.rds")
Fibroblast_A_C <- subset(Fibroblast, group2 %in% c('A', 'C'))
Fibroblast_A_C$subtype %>% table()
Fibroblast <- subset(Fibroblast_A_C, subtype %in% c('LENG8-mCAF', 'THBS1-mCAF', 'iCAF', 'tCAF', 'ifnCAF', 'vCAF'))
Fibroblast$subtype %>% table()
Fibroblast$group2 %>% table()
saveRDS(Fibroblast, file.path(out_dir, 'Fibroblast_subtype_A_C_filtered.rds'))

### dimplot ------------
### PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

DefaultAssay(Fibroblast) <- "RNA"
Fibroblast <- FindVariableFeatures(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')
Fibroblast <- RunPCA(Fibroblast)


#### pca inspector ------
pdf(file.path(Fibroblast_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()


### dimension reduction --------
dims <- 10

Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "group", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast_out_dir, str_glue("Dimplot dim{dims} group2.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "group2", label = TRUE) + NoLegend()
dev.off()


Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample", 
    # theta = 3, # default 2, bigger value -- more diversity per cluster
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


DimPlot(Fibroblast, cols = allcolour, label = T, alpha = 0.4, label.size = 3, group.by = 'subtype', reduction = 'harmonyumap') + NoLegend()
ggsave(file.path(out_dir, 'Fibroblast_umap_celltype_noLegend.pdf'), width = 6, height = 5.5)
DimPlot(Fibroblast, cols = allcolour, label = F, alpha = 0.4, label.size = 3, group.by = 'subtype', reduction = 'harmonyumap')
ggsave(file.path(out_dir, 'Fibroblast_umap_celltype.pdf'), width = 7.5, height = 5.5)

DimPlot(Fibroblast, cols = allcolour, label = T, alpha = 0.4, label.size = 3, split.by = 'group2', ncol = 2, group.by = 'subtype', reduction = 'harmonyumap') + NoLegend()
ggsave(file.path(out_dir, 'Fibroblast_umap_celltype_noLegend_split.pdf'), width = 10, height = 5.5)
DimPlot(Fibroblast, cols = allcolour, label = F, alpha = 0.4, label.size = 3, split.by = 'group2', ncol = 2, group.by = 'subtype', reduction = 'harmonyumap')
ggsave(file.path(out_dir, 'Fibroblast_umap_celltype_split.pdf'), width = 11, height = 5.5)


### dotplot -------------
Fibroblast_marker_ls <- c(
  'mCAF' = 'MMP11', 
  'mCAF' = 'POSTN', 
  "LENG8-mCAF" = 'LENG8', 
  "THBS1-mCAF" = 'THBS1', 
  
  "iCAF" = 'IGF1', 
  "iCAF" = 'ITM2A', 
  "iCAF" = 'CXCL12', 
  "iCAF" = 'CD34', 
  
  "tCAF" = 'NDRG1', 
  "tCAF" = 'MME', 
  "tCAF" = 'SLC2A1', 
  "tCAF" = 'BNIP3', 
  "tCAF" = 'ADM', 
  "tCAF" = 'ANGPTL4', 
  
  "ifnCAF" = 'IDO1', 
  "ifnCAF" = 'CXCL9', 
  "ifnCAF" = 'CXCL10', 
  "ifnCAF" = 'CXCL11', 
  
  "vCAF" = 'NOTCH3', 
  "vCAF" = 'COL18A1', 
  "vCAF" = 'MCAM', 
  "vCAF" = 'RGS5'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

DotPlot(Fibroblast, features = Fibroblast_marker_ls, group.by = 'subtype') + theme_90angle
ggsave(file.path(out_dir, 'Fibroblast_dotplot_celltype.pdf'), width = 12, height = 5)


### col plot ---------------
Fibroblast_meta <- Fibroblast@meta.data

p_col_proj <-
  Fibroblast_meta %>%
  group_by(group2, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group2, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(out_dir, "Fibroblast_cell_freq_group2.pdf"),
  height = 5,
  width  = 4
)


### enrichment --------------- 
Idents(Fibroblast) <- "subtype"
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Fibroblast_markers,file.path(out_dir, 'Fibroblast_DEGs_subtype.csv'))


Fibroblast_de_genes_for_enrichmenet <- 
  Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Fibroblast_de_genes_for_enrichmenet %>%
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


Fibroblast_er_output <- file.path(out_dir, 'Fibroblast_enrichment') 
dir.create(Fibroblast_er_output)


# expot enrichment table
enrich_ls %>% 
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Fibroblast_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Fibroblast_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(Fibroblast_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Fibroblast_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Fibroblast_de_genes_id_for_enrichmenet <-
  Fibroblast_de_genes_for_enrichmenet %>%
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
  Fibroblast_de_genes_id_for_enrichmenet %>%
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

# saveRDS(enrich_ls,file.path(Fibroblast_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Fibroblast_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Fibroblast_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Fibroblast_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Fibroblast_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })




## monocle ----------
### CAF trajactory -------
library(monocle)
source("monocle fix bugs.R")
source('color_pallete.R')

Fibroblast <- readRDS("Fibroblast_subtype_remove_D2/Fibroblast.rds")
Fibroblast_A_C <- subset(Fibroblast, group2 %in% c('A', 'C'))
Fibroblast_A_C$subtype %>% table()
CAF <- subset(Fibroblast_A_C, subtype %in% c('LENG8-mCAF', 'THBS1-mCAF', 'iCAF', 'tCAF', 'ifnCAF', 'vCAF'))
CAF$subtype %>% table()

CAF <- CAF[, CAF$nCount_RNA > 100]


### set default ------- 
monocle_out_dir <- file.path(out_dir, "monocle_result_caf")
dir.create(monocle_out_dir)

VariableFeatures(CAF)

# subset
# set.seed(1234)
# meta_CAF_subset <-
#   CAF[['group2']] %>%
#   rownames_to_column("cell_id") %>%
#   group_by(group2) %>%
#   slice_sample(n = 5000)
meta_CAF_subset <- CAF[['group']] %>% rownames_to_column("cell_id")

### cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(CAF), "counts"), "sparseMatrix")
p_data_monocle <- CAF[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(CAF),
  row.names       = rownames(CAF)
)


cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_CAF_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_CAF_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


# cds <- clusterCells(cds)


# ### 使用seurat选择的高变基因⚠️
# CAF <- FindVariableFeatures(CAF, nfeatures = 2000)
# express_genes_monocle <- VariableFeatures(CAF)
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

# ### 使用clusters差异表达基因
# Idents(CAF) <- "group2"
# deg.cluster <- FindAllMarkers(JoinLayers(CAF))
# express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

#### 使用monocle选择的高变基因⚠️ --------
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.001 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

# #这一步输入的expressed_genes来自于步骤4。
# fData(cds)$num_cells_expressed <- rowSums(exprs(cds) > 0)
# expressed_genes2 <- row.names(subset(fData(cds), num_cells_expressed >= 10))
# diff <- differentialGeneTest(cds[expressed_genes2, ],fullModelFormul aStr="~annotation",cores=1) 
# 
# head(diff)
# 
# ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
# deg <- subset(diff, qval < 0.01) #选出2829个基因
# deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##差异基因的结果文件保存
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
# 
# ## 轨迹构建基因可视化
# ordergene <- rownames(deg) 
# cds <- setOrderingFilter(cds, ordergene)  
# plot_ordering_genes(cds)


cds <- reduceDimension(cds)
# undebug(orderCells)
cds <- orderCells(cds)

# Error:
#   ! `nei()` was deprecated in igraph 2.1.0 and is now defunct.
# ℹ Please use `.nei()` instead.
# install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz', type = "source", repos = NULL)


# plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_seuratHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="subtype", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_seuratHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="group2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_group2_seuratHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_seuratHVGs.pdf")), height = 6, width = 8)
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_seuratHVGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, subtype ) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = subtype ))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - group2_seuratHVGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, group2) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = group2))
# dev.off()
# 
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_seuratHVGs.pdf"), width = 12, height = 6)
# local({
#   Biobase::pData(cds)$subtype <- CAF[['subtype']][rownames(Biobase:::pData(cds)), ]
#   
#   # cds2 <- cds[, !is.na(pData(cds)$subtype)]
#   
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ subtype) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by group2_seuratHVGs.pdf"), width = 12, height = 6)
# local({
#   Biobase::pData(cds)$group2 <- CAF[['group2']][rownames(Biobase:::pData(cds)), ]
#   
#   # cds2 <- cds[, !is.na(pData(cds)$group2)]
#   
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ group2) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()




# plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="subtype", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="group2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_group2_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_clusterDEGs.pdf")), height = 6, width = 8)
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_clusterDEGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, subtype ) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = subtype ))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - group2_clusterDEGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, group2) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = group2))
# dev.off()
# 
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_clusterDEGs.pdf"), width = 12, height = 6)
# local({
#   Biobase::pData(cds)$subtype <- CAF[['subtype']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$subtype)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ subtype) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by group2_clusterDEGs.pdf"), width = 8, height = 3)
# local({
#   Biobase::pData(cds)$group2 <- CAF[['group2']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$group2)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ group2) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()




plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="subtype", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="group2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_group2_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_monocleHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_monocleHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, subtype ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = subtype ))
dev.off()
pdf(file.path(monocle_out_dir, "barplot of cell props group by State - group2_monocleHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, group2) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group2))
dev.off()

pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_monocleHVGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$subtype <- CAF[['subtype']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$subtype)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ subtype) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir, "cell trajectory - State - split by group2_monocleHVGs.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$group2 <- CAF[['group2']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$group2)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ group2) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()




p_col_state <- cds %>% 
  Biobase::pData() %>%
  dplyr::select(group2, State) %>%
  group_by(group2, State) %>%
  mutate(group2 = as.character(group2)) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group2, y = freq)) +
  geom_col(aes(fill = State)) +
  # scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
p_col_annotation <- cds %>% 
  Biobase::pData() %>%
  dplyr::select(group2, State) %>%
  group_by(State, group2) %>%
  mutate(group2 = as.character(group2)) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group2)) +
  # scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
p_col_state + p_col_annotation





cds <- orderCells(cds, root_state = 3)

# saveRDS(cds, file = file.path(monocle_out_dir, "cds_seuratHVGs.rds"))
# saveRDS(cds, file = file.path(monocle_out_dir, "cds_clusterDEGs.rds"))
saveRDS(cds, file = file.path(monocle_out_dir, "cds_monocleHVGs.rds"))




# 