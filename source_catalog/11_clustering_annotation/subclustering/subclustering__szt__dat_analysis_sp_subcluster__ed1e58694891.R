library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")

library(future)
plan(multicore, workers = 6)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 10Gb



all_data <- readRDS("output/data/adata_celltype.rds")
all_data_raw <- readRDS("output/data/adata_raw.rds")
# sp_epi <- readRDS("output/data/adata_epi.rds")

out_dir <- "sp_subcluster"
dir.create(out_dir)

all_data$annotation_full <- all_data$celltype
all_data$annotation_full %>% table
# AT2 cell                    B          Endothelial          Epithelials           Fibroblast 
# 44199                26739                70207               462241               139348 
# Mast cell               T cell marcophage_monocytes               plasma 
# 9849               116559               122606                71003 

major_marker_ls <- c(
  
  "B cell" = "CD19",
  "B cell" = "CD79A",
  "B cell" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "MZB1",
  "T cell" = "CD3E",
  "T cell" = "CD4",
  "T cell" = "FOXP3",
  "T cell" = "CD8A",
  "T cell" = "GZMB",
  "NK" = "NCAM1",
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Mast cell" = "CTSG",
  "Mast cell" = "MS4A2",
  "Mast cell" = "KIT",
  "Mast cell" = "HDC",
  "Monocytes" = "CD14",
  "Monocytes" = "FCGR3A",
  "Dendritic" = "FCER1A",
  "Dendritic" = "CD1C",
  "Dendritic" = "CD83",
  "Dendritic" = "LILRA4",
  "Dendritic" = "CLEC9A",
  "Dendritic" = "IDO1",
  "Neutrophi" = "CSF3R",
  "Neutrophi" = "FCGR3B",
  "Endothelial" = "PECAM1",
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT20",
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "CCDC80",
  "Smooth muscle cells" = "PDGFRB",
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Smooth muscle cells" = "TAGLN",
  "Pericyte" = "RGS5",
  "Cholangiocyte" = "KRT7",
  "Cholangiocyte" = "KRT19",
  "Cholangiocyte" = "SOX9",
  "Cholangiocyte" = "CFTR",
  "Hepatocyte" = "ALB",
  "Hepatocyte" = "CYP3A4",
  "Hepatocyte" = "CYP2E1",
  "Hepatocyte" = "APOB",
  "Hepatocyte" = "APOE",
  "Hepatocyte" = "HP",
  "Liver sinusoidal endothelial cells" = "CLEC4G",
  "Liver sinusoidal endothelial cells" = "CLEC4M",
  "Hepatic stellate cells" = "ACTA2",
  # "Kupffer cell" = "CD68",
  "Kupffer cell" = "MAFB",
  "Kupffer cell" = "VSIG4"
)


# T subcluster ----------------
sp_t <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "T"] ]
sp_t$nCount_RNA <- sp_t$total_counts
sp_t$project <- sp_t$sample

DefaultAssay(sp_t) <- 'RNA'

sp_t <- sc_normalize(sp_t, scale.factor = "def")
sp_t <- sc_pca(
  sp_t, 
  sc_out_dir = file.path(out_dir, "t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_t <- 15

sp_t <- sc_dr(
  sp_t, 
  ndims = ndims_t,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "t/sc.2 dr"),
  dir.create.force = TRUE
)


sp_t <- sc_cluster(
  sp_t, 
  ndims = ndims_t,
  k.param = 30,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "t/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

t_markers <- c(
  Plasmablast = 'XBP1',
  Plasmablast = 'MZB1',
  
  'Stromal cell' = 'COL5A1',
  'Stromal cell' = 'THY1',
  
  Macrophage = 'CD14',
  Macrophage = 'MRC1',
  
  
  'T cell' = 'CD3D',
  'T cell' = 'CD3E',
  'T cell' = 'CD4',
  'T cell' = 'CD8A',
  
  'Smooth muscle cell' = 'HSPA6',
  'Smooth muscle cell' = 'HSPA8',
  'Smooth muscle cell' = 'HSPD1'
  
)

sp_t <- sc_dotplot(
  sp_t, 
  ndims = ndims_t,
  marker = c(
    "CD3D", "CD3E",
    "CD4",
    "CD40LG",
    "CD8A",
    "FCGR3A",
    "NCAM1",
    "MKI67",
    "TOP2A"
  ),
  width = 8, 
  height = 6,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot-2"), 
  dir.create.force = TRUE
)

sp_t <- sc_dotplot(
  sp_t, 
  ndims = ndims_t,
  marker = t_markers,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_t <- sc_dotplot(
  sp_t, 
  ndims = ndims_t,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

# dotplot of t cluster
sp_t <- sc_dotplot(
  ,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot - final"),
  dir.create.force = TRUE
)


p_dot_anno <-
  DotPlot(
    sp_t,
    group.by = "annotation_t",
    features = c(
      'T' = "CD3D", 
      'T' = "CD3E",
      'CD4T' = "CD4",
      # 'T' = "CD40LG",
      'CD8T' = "CD8A",
      'Tregs' = "FOXP3",
      'gamma delta T' = "TRDC",
      Cycling = "MKI67",
      Cycling = "TOP2A",
      
      Hepatocyte = "APOB",
      Cholangiocyte = "SOX9",
      Endothelial = "PECAM1",
      Fibroblast = "CCDC80"
    )
  ) +
  theme(strip.clip = 'off')
    
pdf(
  file.path("sp_subcluster//t/sc.7 dotplot - final.pdf"),
  width  = 15,
  height = 6
)
p_dot_anno %>% print()
dev.off()






all_markers_t <- 
  sc_findmarker(
    sp_t,
    idents =  "RNA_snn_res.0.4", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


sc_findmarker(sp_t, all_markers = all_markers_t, sort_by = "pct")

dir.create(file.path(out_dir, "t/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "t/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.4", # "sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_t[[c('RNA_snn_res.0.4', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col()

sc_mark_heatmap(
  sp_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

t_cell_cluster <- c(
  '0'  = "CD8T",
  '1'  = "Hepatocyte-like T cell",
  '2'  = "Cholangiocyte-like T cell",
  '3'  = "Fibroblast-like T cell",
  '4'  = "CD4T",
  '5'  = "Endothelial-like T cell",
  '6'  = "Treg",
  '7'  = "Gamma delta T cell",
  '8'  = "Cycling T cell",
  '9'  = "Cycling T cell",
  '10'  = "Cholangiocyte-like T cell"
)



Idents(sp_t) <- "RNA_snn_res.0.4"
sp_t <- RenameIdents(sp_t, t_cell_cluster)
sp_t$annotation_t <- Idents(sp_t)

dir.create(file.path(out_dir, "t/sc.5 annotation"))
pdf(file.path(out_dir, 't/sc.5 annotation', str_glue("Dimplot dim{ndims_t} t annotation.pdf")), height = 6, width = 8)
DimPlot(sp_t, group.by = "annotation_t", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_t, group.by = "annotation_t", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_t, file.path(out_dir, 'sp_t.rds'))

# cd4T subcluster ----------------
sp_cd4t <- sp_t[, sp_t$annotation_t == "CD4T" ]
# sp_cd4t <- sp_t[, sp_t$RNA_snn_res.0.2 %in% c('1', '2') ]
# sp_cd4t$nCount_RNA <- sp_cd4t$total_counts
# sp_cd4t$project <- sp_cd4t$sample

DefaultAssay(sp_cd4t) <- 'RNA'

sp_cd4t <- sc_normalize(sp_cd4t, scale.factor = "def")
sp_cd4t <- sc_pca(
  sp_cd4t, 
  sc_out_dir = file.path(out_dir, "cd4t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_cd4t <- 10

sp_cd4t <- sc_dr(
  sp_cd4t, 
  ndims = ndims_cd4t,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "cd4t/sc.2 dr"),
  dir.create.force = TRUE
)


sp_cd4t <- sc_cluster(
  sp_cd4t, 
  ndims = ndims_cd4t,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "cd4t/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sp_cd4t <- sc_dotplot(
  sp_cd4t, 
  ndims = ndims_cd4t,
  marker = c(
    "CCDC80", "POSTN",
    
    "CD3D", "CD3E",
    "CD4",
    "CD8A",
    "NCAM1",
    "MKI67",
    "TOP2A",
    "CXCL13",
    "CD200",
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "SELL",
    "TCF7",
    
    "GZMA", "GZMH", "KLRG1",
    "TRDC", "TRAC",
    "CD14","CD68","FCN1",
    "HAVCR2"
  ),
  width = 15, 
  height = 6,
  sc_out_dir = file.path(out_dir, "cd4t/sc.4 dotplot"), 
  dir.create.force = TRUE
)

sp_cd4t <- sc_dotplot(
  sp_cd4t, 
  ndims = ndims_cd4t,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "cd4t/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_cd4t <- 
  sc_findmarker(
    sp_cd4t,
    idents =  "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


sc_findmarker(sp_cd4t, all_markers = all_markers_cd4t, sort_by = "pct")

dir.create(file.path(out_dir, "cd4t/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "cd4t/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_cd4t,
  markers = all_markers_cd4t,
  idents = "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_cd4t[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sp_cd4t,
  markers = all_markers_cd4t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

cd4t_cell_cluster <- c(
  '0'  = "CD8-Temra",
  '1'  = "Fibroblast-like T cell",
  '2'  = "CD4-Treg"
)



Idents(sp_cd4t) <- "RNA_snn_res.0.2"
sp_cd4t <- RenameIdents(sp_cd4t, cd4t_cell_cluster)
sp_cd4t$annotation_cd4t <- Idents(sp_cd4t)

dir.create(file.path(out_dir, "cd4t/sc.5 annotation"))
pdf(file.path(out_dir, 'cd4t/sc.5 annotation', str_glue("Dimplot dim{ndims_cd4t} t annotation.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "annotation_cd4t", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_cd4t, group.by = "annotation_cd4t", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_cd4t, file.path(out_dir, 'sp_cd4t.rds'))


# cd8T subcluster ----------------
sp_cd8t <- sp_t[, sp_t$annotation_t == "CD8T" ]
# sp_cd8t <- sp_t[, sp_t$RNA_snn_res.0.2 == "0" ]


# sp_cd8t$nCount_RNA <- sp_cd8t$total_counts
# sp_cd8t$project <- sp_cd8t$sample

DefaultAssay(sp_cd8t) <- 'RNA'

sp_cd8t <- sc_normalize(sp_cd8t, scale.factor = "def")
sp_cd8t <- sc_pca(
  sp_cd8t, 
  sc_out_dir = file.path(out_dir, "cd8t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_cd8t <- 7

sp_cd8t <- sc_dr(
  sp_cd8t, 
  ndims = ndims_cd8t,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "cd8t/sc.2 dr"),
  dir.create.force = TRUE
)


sp_cd8t <- sc_cluster(
  sp_cd8t, 
  ndims = ndims_cd8t,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "cd8t/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  k.param = 15,
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sp_cd8t <- sc_dotplot(
  sp_cd8t, 
  ndims = ndims_cd8t,
  marker = c(
    "CCDC80", "POSTN",
    
    "CD3D", "CD3E",
    "CD4",
    "CD8A",
    "NCAM1",
    "MKI67",
    "TOP2A",
    "CXCL13",
    "CD200",
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "SELL",
    "TCF7",
    "MZB1","XBP1",
    "GZMA", "GZMH", "KLRG1",
    "CD14","CD68","FCN1",
    "LAMP3",
    "XCR1",
    "CD1C",
    "CLEC9A",
    "TRDC", "TRAC",  
    "HAVCR2"
  ),
  width = 18, 
  height = 6,
  sc_out_dir = file.path(out_dir, "cd8t/sc.4 dotplot"), 
  dir.create.force = TRUE
)

sp_cd8t <- sc_dotplot(
  sp_cd8t, 
  ndims = ndims_cd8t,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "cd8t/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_cd8t <- 
  sc_findmarker(
    sp_cd8t,
    idents =  "RNA_snn_res.0.4", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


sc_findmarker(sp_cd8t, all_markers = all_markers_cd8t, sort_by = "pct")

dir.create(file.path(out_dir, "cd8t/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "cd8t/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_cd8t,
  markers = all_markers_cd8t,
  idents = "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_cd8t[[c('RNA_snn_res.0.4', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sp_cd8t,
  markers = all_markers_cd8t,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

cd8t_cell_cluster <- c(
  '0'  = "CD8-Tem",
  '1'  = "CD4-Tfh",
  '2'  = "CD8-Temra",
  '3'  = "Cycling T cell",
  '4'  = "Fibrobalst-like T cell",
  '5'  = "Plasma-like T cell",
  '6'  = "DC-like T cell"
)



Idents(sp_cd8t) <- "RNA_snn_res.0.4"
sp_cd8t <- RenameIdents(sp_cd8t, cd8t_cell_cluster)
sp_cd8t$annotation_cd8t <- Idents(sp_cd8t)

dir.create(file.path(out_dir, "cd8t/sc.5 annotation"))
pdf(file.path(out_dir, 'cd8t/sc.5 annotation', str_glue("Dimplot dim{ndims_cd8t} t annotation.pdf")), height = 6, width = 8)
DimPlot(sp_cd8t, group.by = "annotation_cd8t", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_cd8t, group.by = "annotation_cd8t", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_cd8t, file.path(out_dir, 'sp_cd8t.rds'))


# macrophage subcluster --------------
sp_mac <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "marcophage_monocytes"] ]
sp_mac$nCount_RNA <- sp_mac$total_counts
sp_mac$project <- sp_mac$sample

DefaultAssay(sp_mac) <- 'RNA'

sp_mac <- sc_normalize(sp_mac, scale.factor = "def")
sp_mac <- sc_pca(
  sp_mac, 
  sc_out_dir = file.path(out_dir, "mac/sc.1 pca"),
  ndims = 50,
  dir.create.force = TRUE
)

ndims_mac <- 20

sp_mac <- sc_dr(
  sp_mac, 
  ndims = ndims_mac,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "mac/sc.2 dr"),
  dir.create.force = TRUE
)


sp_mac <- sc_cluster(
  sp_mac, 
  ndims = ndims_mac,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "mac/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

mac_markers <- c(
  'PTPRC',
  'CD3D',
  'CD3G',
  'CD4',
  'CD8A',
  'NCAM1',
  'CSF3R',
  'CSF1R',
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "FCN1", 
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  # pDC
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  # cDC1
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  #cDC2
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'TOP2A',
  'MKI67',
  'TREM2',
  
  # 'CXCR4',
  'CYP27A1',
  'MARCO',
  'CXCL10'
)

sp_mac <- sc_dotplot(
  sp_mac, 
  ndims = ndims_mac,
  marker = mac_markers,
  sc_out_dir = file.path(out_dir, "mac/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_mac <- sc_dotplot(
  sp_mac, 
  ndims = ndims_mac,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "mac/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_mac <- 
  sc_findmarker(
    sp_mac,
    idents =  "RNA_snn_res.0.4", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )

VlnPlot(
  sp_mac,
  features = "total_counts"
)
sc_findmarker(sp_mac, all_markers = all_markers_mac, sort_by = "pct")
all_markers_mac %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(pct.1), .by_group = TRUE) %>%
  slice_head(n = 5)


dir.create(file.path(out_dir, "mac/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "mac/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_mac,
  markers = all_markers_mac,
  idents = "RNA_snn_res.0.4", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_mac[[c('RNA_snn_res.0.4', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.4)) +
  geom_col()

sc_mark_heatmap(
  sp_mac,
  markers = all_markers_mac,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

mac_cell_cluster <- c(
  '0'  = "low quality cells", # 无差异基因，nCount显著低于其他群
  '1'  = "CYP27A1 Macrophage",
  '2'  = "Fibroblast-like Macrophage",
  '3'  = "Epithelial-like Macrophage",
  '4'  = "Macro Macrophage", # 组织驻留型
  '5'  = "Proliferating Macrophage",
  '6'  = "CXCL10 Macrophage", # INF response
  '7'  = "cDC2",
  '8'  = "cDC1"
)



Idents(sp_mac) <- "RNA_snn_res.0.4"
sp_mac <- RenameIdents(sp_mac, mac_cell_cluster)
sp_mac$annotation_mac <- Idents(sp_mac)

dir.create(file.path(out_dir, "mac/sc.5 annotation"))
pdf(file.path(out_dir, 'mac/sc.5 annotation', str_glue("Dimplot dim{ndims_t} mac annotation.pdf")), height = 6, width = 8)
DimPlot(sp_mac, group.by = "annotation_mac", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mac, group.by = "annotation_mac", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_mac, file.path(out_dir, 'sp_mac.rds'))
# Epi subcluster ----------------
sp_epi <- all_data_raw[, colnames(all_data)[all_data$celltype %in% c('Hepatocyte', 'Tumor cell')] ]
sp_epi$nCount_RNA <- sp_epi$total_counts
sp_epi$project <- sp_epi$sample

DefaultAssay(sp_epi) <- 'RNA'
sp_epi[['RNA']]  <- sp_epi[['RNA']] %>% split(sp_epi$project)

# sp_epi <- sc_normalize(sp_epi, scale.factor = "median")

sp_epi <- sc_normalize(sp_epi, scale.factor = "median", method = "SCT", nfeatures = 2000, ncells = 1e4)

sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

## group2: LR与PH的对比  
sp_epi$group2 <- sp_epi$group %>% str_extract("^[A-Z]{2}")

# others group
grp_info <- readxl::read_excel("annotation.xlsx")

grp_info_full <- 
  sp_epi[['sample']] %>%
  rownames_to_column("id") %>%
  left_join(grp_info, by = c('sample')) %>%
  column_to_rownames("id")


sp_epi <- 
  sp_epi %>%
  AddMetaData(grp_info_full)
sp_epi$paired_sample <- sp_epi$`paired sample`


## project harmony ----
ndims_epi <- 20
# ndims_epi <- 10


sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  group.by = "project",
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr - project"),
  dir.create.force = TRUE
)

sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster - project"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8, 1)
)


## group2 harmony ----
# ndims_epi <- 20
ndims_epi <- 20


sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  group.by = "group2",
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr - group2"),
  dir.create.force = TRUE
)

sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster - group2"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8, 1)
)


## paried sample harmony ----
# 20
ndims_epi <- 20 # 10 # 20 # 30
sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  group.by = "paired_sample",
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr - paired sample"),
  dir.create.force = TRUE
)
sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster - paired sample"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8, 1)
)
# saveRDS(sp_epi, "sp_subcluster/sp_epi-pariedsample-dim20.rds")


# dimplot
pdf(file.path(file.path(out_dir, "epi/sc.3 cluster"), str_glue("Dimplot dim{ndims_epi} cluster - harmony.pdf")), height = 6, width = 8)
walk(
  str_glue("RNA_snn_res.{c(0.2, 0.4, 0.6, 0.8, 1)}"),
  function(group.by){
    p <- DimPlot(sp_epi, group.by = group.by, label = TRUE, reduction = "umap") + NoLegend()
    print(p)
  })
# harmony
walk(
  str_glue("RNA_snn_res.{c(0.2, 0.4, 0.6, 0.8, 1)}"),
  function(group.by){
    p <- DimPlot(sp_epi, group.by = group.by, label = TRUE, reduction = 'harmonyumap') + NoLegend()
    print(p)
  })
  
dev.off()


# subcluster 
Idents(sp_epi) <- "RNA_snn_res.0.6"
sp_epi <- 
  FindSubCluster(
    sp_epi,
    cluster = '0',
    graph.name = 'RNA_snn',
    subcluster.name = "sub.cluster",
    resolution = 0.3,
    algorithm  = 1
  )

pdf(file.path(file.path(out_dir, "epi"), str_glue("Dimplot subcluster0 - harmony.pdf")), height = 6, width = 8)
  DimPlot(
    sp_epi,
    group.by = 'sub.cluster',
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()
dev.off()

pdf(file.path(file.path(out_dir, "epi"), str_glue("Dimplot subcluster0 - split - harmony.pdf")), height = 26, width = 16)
DimPlot(
  sp_epi,
  group.by = 'sub.cluster',
  split.by = 'sample',
  reduction = "harmonyumap",
  ncol = 4
) +
  ggsci::scale_color_igv()
dev.off()




# 
epi_markers <- c(
  # 'EPCAM',
  "CD3D",
  # "CD3E",
  # "CD8A",
  
  'B cell' = 'CD19',
  'B cell' = 'CD79A',
  'B cell' = 'MS4A1',
  'Plasma' = 'TNFRSF17',
  'Plasma' = 'MZB1',
  'T cell' = 'CD3E',
  'T cell' = 'CD4',
  'T cell' = 'FOXP3',
  'T cell' = 'CD8A',
  'T cell' = 'GZMB',
  'NK' = 'NCAM1',
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  'Myeloid' = 'ITGAM',
  'Macrophage' = 'CD68',
  'Macrophage' = 'CD163',
  'Mast cell' = 'CTSG',
  'Mast cell' = 'MS4A2',
  'Mast cell' = 'KIT',
  'Mast cell' = 'HDC',
  'Monocytes' = 'CD14',
  'Monocytes' = 'FCGR3A',
  'Dendritic' = 'FCER1A',
  'Dendritic' = 'CD1C',
  'Dendritic' = 'CD83',
  'Dendritic' = 'LILRA4',
  'Dendritic' = 'CLEC9A',
  'Dendritic' = 'IDO1',
  'Neutrophi' = 'CSF3R',
  'Neutrophi' = 'FCGR3B',
  'Endothelial' = 'PECAM1',
  'Epithelial' = 'EPCAM',
  'Epithelial' = 'KRT20',
  'Fibroblasts' = 'PDGFRA',
  'Fibroblasts' = 'CCDC80',
  'Smooth muscle cells' = 'PDGFRB',
  'Smooth muscle cells' = 'PLN',
  'Smooth muscle cells' = 'MYLK',
  'Smooth muscle cells' = 'RGS5',
  'Cholangiocyte' = 'KRT7',
  'Cholangiocyte' = 'KRT19',
  'Cholangiocyte' = 'SOX9',
  'Hepatocyte' = 'ALB',
  'Hepatocyte' = 'CYP3A4',
  'Hepatocyte' = 'CYP2E1',
  'Hepatocyte' = 'APOB',
  'Hepatocyte' = 'APOE',
  'Hepatocyte' = 'HP',
  'Liver sinusoidal endothelial cells' = 'CLEC4G',
  'Liver sinusoidal endothelial cells' = 'CLEC4M',
  'Hepatic stellate cells' = 'ACTA2',
  # 'Kupffer cell' = 'CD68',
  'Kupffer cell' = 'MAFB',
  'Kupffer cell' = 'VSIG4'
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  width = 30,
  dir.create.force = TRUE
)
sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

# sp_epi2 <- sp_epi[, !sp_epi$RNA_snn_res.0.2 %in% c('13', '14', '15')]

all_markers_epi <- 
  sc_findmarker(
    sp_epi2,
    idents = "RNA_snn_res.0.6",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sp_epi, all_markers = all_markers_epi, sort_by = "pct")


all_markers_epi %>%
  write.csv(
    file = file.path(out_dir, "epi/de analysis - RNA_snn_res.0.6.csv")
  )


# epi top 5 markers
all_markers_epi_top5 <-
  all_markers_epi %>%
  group_by(cluster) %>% 
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 1
  ) %>%
  dplyr::slice_max(order_by = pct.1 - pct.2, n = 5)

p_dot_epi_2 <- 
  DotPlot(
    sp_epi,
    features = unique(all_markers_epi_top5$gene),
    group.by = "RNA_snn_res.0.6"
  ) +
  RotatedAxis()

ggsave(
  p_dot_epi_2,
  filename = file.path(out_dir, "epi/Dotplot - top 5 epi markers - RNA_snn_res.0.6.pdf"),
  height = 8,
  width  = 16
)


# dimplot of sample 
p_dim_sample <- 
  DimPlot(
    sp_epi,
    group.by = "sample",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()


ggsave(
  p_dim_sample,
  filename = file.path(out_dir, str_glue("epi/dimplot dim{ndims_epi} - sample.pdf")),
  height = 6,
  width  = 8
)


p_dim_sample_split <- 
  DimPlot(
    sp_epi,
    group.by = "RNA_snn_res.0.6",
    split.by = "sample",
    reduction = "harmonyumap",
    label = TRUE,
    ncol = 4
  ) &
  ggsci::scale_color_igv()

ggsave(
  p_dim_sample_split,
  filename = file.path(out_dir, str_glue("epi/dimplot dim{ndims_epi} - sample - split.pdf")),
  height = 26,
  width  = 16
)


# dimplot of group 
p_dim_group <- 
  DimPlot(
    sp_epi,
    group.by = "group",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()


ggsave(
  p_dim_group,
  filename = file.path(out_dir, str_glue("epi/dimplot dim{ndims_epi} - group.pdf")),
  height = 6,
  width  = 8
)


# dimplot of group 
p_dim_some_sample <- 
  DimPlot(
    sp_epi[, sp_epi$sample %in% c("7320N", "7140N", "2484N", "5547N", "5785", "2484T")],
    group.by = c("sample", "RNA_snn_res.0.6"),
    reduction = "harmonyumap"
  ) &
  ggsci::scale_color_igv()


ggsave(
  p_dim_some_sample,
  filename = file.path(out_dir, str_glue("epi/dimplot dim{ndims_epi} - some samples.pdf")),
  height = 6,
  width  = 12
)
##  only subcluster
pdf(
  file = file.path(out_dir, str_glue("epi/dimplot dim{ndims_epi} - some samples - subcluster.pdf")),
  height = 6,
  width  = 12
)
  DimPlot(
    sp_epi[
      , 
      sp_epi$sample %in% c("7320N", "7140N", "2484N", "5547N", "5785", "2484T") &
        str_detect(sp_epi$sub.cluster, "_" )
    ],
    group.by = c("sample", "sub.cluster"),
    reduction = "harmonyumap"
  ) &
  ggsci::scale_color_igv()
dev.off()



# vlnplot dotplot of some gene
p_vln_ubd <-
  VlnPlot(
    sp_epi,
    features = "UBD",
    # pt.size = 0,
    group.by = "RNA_snn_res.0.6"
  )

ggsave(
  p_vln_ubd,
  filename = file.path(out_dir, "epi/vlnplot of UBD - RNA_snn_res.0.6.pdf"),
  height = 6,
  width  = 8
)


p_dot_epi_1 <- 
  DotPlot(
    sp_epi,
    features = c(
      'APOB', 'SOX9', 'EPCAM', 'MKI67', 'TOP2A', 'UBD', 'MET',
      'CD14', "CD68", "CD3E", "CD19", "MS4A1", "PECAM1", "CLEC4M", "CCDC80"
    ),
    group.by = "RNA_snn_res.0.6"
  )

ggsave(
  p_dot_epi_1,
  filename = file.path(out_dir, "epi/Dotplot - RNA_snn_res.0.6.pdf"),
  height = 8,
  width  = 16
)

# cell props 
p_cell_props_ls <-
  str_glue("RNA_snn_res.0.{seq(2,8, by = 2)}") %>%
  c("sub.cluster") %>%
  setNames(., .) %>%
  map(function(x){
    sp_epi[[]] %>% 
      group_by(group, .data[[x]]) %>%
      summarise(count = n()) %>%
      
      mutate(freq = count / sum(count)) %>%
      
      ggplot(aes(x = group, y = freq, fill = .data[[x]])) +
      geom_col() +
      ggsci::scale_fill_igv()
  })

p_cell_props_ls2 <-
  str_glue("RNA_snn_res.0.{seq(2,8, by = 2)}") %>%
  c("sub.cluster") %>%
  setNames(., .) %>%
  map(function(x){
    sp_epi[[]] %>% 
      group_by(.data[[x]], group) %>%
      summarise(count = n()) %>%
      
      mutate(freq = count / sum(count)) %>%
      
      ggplot(aes(x = .data[[x]], y = freq, fill = group)) +
      geom_col() +
      ggsci::scale_fill_igv()
  })

p_cell_props_sample_ls <-
  str_glue("RNA_snn_res.0.{seq(2,8, by = 2)}") %>%
  c("sub.cluster") %>%
  setNames(., .) %>%
  map(function(x){
    sp_epi[[]] %>% 
      group_by(sample, .data[[x]]) %>%
      summarise(count = n()) %>%
      
      mutate(freq = count / sum(count)) %>%
      
      ggplot(aes(x = sample, y = freq, fill = .data[[x]])) +
      geom_col() +
      ggsci::scale_fill_igv() +
      RotatedAxis() +
      ggtitle(x)
  })

# export
pdf(
  file = file.path("sp_subcluster/epi/", str_glue("cell props - dim{ndims_epi} - group.pdf")),
  height = 6,
  width = 8
)
p_cell_props_ls %>% walk(print)
dev.off()


pdf(
  file = file.path("sp_subcluster/epi/", str_glue("cell props - dim{ndims_epi} - group - celltype(x).pdf")),
  height = 6,
  width = 8
)
p_cell_props_ls2 %>% walk(print)
dev.off()

pdf(
  file = file.path("sp_subcluster/epi/", str_glue("cell props - dim{ndims_epi} - sample.pdf")),
  height = 6,
  width  = 8
)
p_cell_props_sample_ls %>% walk(print)
dev.off()



# heatmap
dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_epi2,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()



sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

epi_cell_cluster <- c(
  # '0'  =  "Hepatocyte",
  "0_0"=  "Hepatocyte", 
  "0_1"=  "Hepatocyte", 
  "0_2"=  "Tumor", 
  "0_3"=  "Tumor LR", 
  '1'  =  "Tumor",
  '2'  =  "Tumor LR",
  '3'  =  "Tumor",
  '4'  =  "Tumor PH",
  '5'  =  "Tumor",
  '6'  =  "Tumor",
  '7'  =  "Tumor PH",
  '8'  =  "Hepatocyte",
  '9'  =  "Tumor PH", # 7298
  '10'  = "Tumor LR",
  '11'  = "Proliferating tumor",
  '12'  = "Cholangiocyte",
  '13'  = "Others",
  '14'  = "Others",
  '15'  = "Others"
)



# Idents(sp_epi) <- "RNA_snn_res.0.6"
Idents(sp_epi) <- "sub.cluster"
sp_epi <- RenameIdents(sp_epi, epi_cell_cluster)
sp_epi$annotation_epi <- Idents(sp_epi)

dir.create(file.path(out_dir, "epi - v3 - pair_sample/sc.5 annotation"))
pdf(file.path(out_dir, 'epi - v3 - pair_sample/sc.5 annotation', str_glue("Dimplot dim{ndims_epi} epi annotation .pdf")), height = 6, width = 8)
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# spatial plot 
cell_position <-
  cbind(
    all_data_raw[[c('sample', 'group')]],
    Embeddings(all_data_raw[['spatial']])
  ) %>%
  rownames_to_column('cell_id')

cell_position_ls <-
  cell_position %>%
  group_by(sample) %>%
  nest() %>%
  deframe()  %>%
  imap(function(dat, name){
    dat %>%
      left_join(
        rownames_to_column(sp_epi[['sub.cluster']], "cell_id"),
        by = "cell_id"
      ) 
  })


sp_p_ls <- 
  cell_position_ls %>% # .[1] %>%
  imap(function(df, name) {
    df %>%
      ggplot(aes(x = SPATIAL_1, y = SPATIAL_2, fill = sub.cluster)) +
      geom_tile(height = 10, width = 10) +
      coord_equal() +
      ggsci::scale_fill_igv(na.value = "grey50") +
      scale_y_reverse() +
      theme_void() +
      ggtitle(name)
  })

# export 
pdf(
  file.path(file.path('sp_subcluster/epi - v3 - pair_sample/', "sp_p.pdf")),
  height = 8,
  width = 8
)
sp_p_ls %>% walk(print)
dev.off()
## 
sp_p_ls2 <- 
  cell_position %>%
  group_by(sample) %>%
  nest() %>%
  deframe()  %>%
  imap(function(dat, name){
    dat %>%
      left_join(
        rownames_to_column(sp_epi[['annotation_epi']], "cell_id"),
        by = "cell_id"
      ) 
  }) %>% # .[1] %>%
  imap(function(df, name) {
    df %>%
      ggplot(aes(x = SPATIAL_1, y = SPATIAL_2, fill = annotation_epi)) +
      geom_tile(height = 10, width = 10) +
      coord_equal() +
      ggsci::scale_fill_igv(na.value = "grey50") +
      scale_y_reverse() +
      theme_void() +
      ggtitle(name)
  })

# export 
pdf(
  file.path(file.path('sp_subcluster/epi - v3 - pair_sample/', "sp_p2.pdf")),
  height = 8,
  width = 8
)
sp_p_ls2 %>% walk(print)
dev.off()


# saveRDS(sp_epi, file.path(out_dir, 'sp_epi.rds'))




# fibro subcluster ----------------
sp_fibro <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "Fibroblast"] ]
sp_fibro$nCount_RNA <- sp_fibro$total_counts
sp_fibro$project <- sp_fibro$sample

DefaultAssay(sp_fibro) <- 'RNA'

sp_fibro <- sc_normalize(sp_fibro, scale.factor = "def")
sp_fibro <- sc_pca(
  sp_fibro, 
  sc_out_dir = file.path(out_dir, "fibro/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_fibro <- 20

sp_fibro <- sc_dr(
  sp_fibro, 
  ndims = ndims_fibro,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "fibro/sc.2 dr"),
  dir.create.force = TRUE
)


sp_fibro <- sc_cluster(
  sp_fibro, 
  ndims = ndims_fibro,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "fibro/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

fibro_markers <-  c(
  Epithelial  = "EPCAM", 
  # Epithelial  = "KRT5",
  Epithelial  = "KRT8",
  Epithelial  = "S100A8",
  Epithelial  = "KRT18",
  
  B = "CD79A",
  plasma = "MZB1",
  plasma = "XBP1",
  
  'Fibrobalst' = 'CCDC80',
  'Fibrobalst' = 'FAP',
  'Fibrobalst' = 'TNXB',

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
  # 'vCAF'='RGS5',
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
  
  'pericyte'='RGS5',
  'pericyte'='COL4A1',
  'pericyte'='COL4A2',
  
  'hsptcaf'='HSPH1',
  'hsptcaf'='VEGFA'
)

sp_fibro <- sc_dotplot(
  sp_fibro, 
  ndims = ndims_fibro,
  marker = fibro_markers,
  sc_out_dir = file.path(out_dir, "fibro/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_fibro <- sc_dotplot(
  sp_fibro, 
  ndims = ndims_fibro,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "fibro/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_fibro <- 
  sc_findmarker(
    sp_fibro,
    idents =  "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )

sc_findmarker(sp_fibro, all_markers = all_markers_fibro, sort_by = "avg")

dir.create(file.path(out_dir, "fibro/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "fibro/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_fibro,
  markers = all_markers_fibro,
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_fibro[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sp_fibro,
  markers = all_markers_fibro,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

fibro_cell_cluster <- c(
  '0'  = "mCAF",
  '1'  = "infCAF",
  '2'  = "dCAF",
  '3'  = "NF", # normal fibroblast
  '4'  = "Plasma-like fibroblast",
  '5'  = "Pericyte",
  '6'  = "CDH1 tCAF",
  '7'  = "SOX2 tCAF",
  '8'  = "Others"
)



Idents(sp_fibro) <- "RNA_snn_res.0.2"
sp_fibro <- RenameIdents(sp_fibro, fibro_cell_cluster)
sp_fibro$annotation_fibro <- Idents(sp_fibro)

dir.create(file.path(out_dir, "fibro/sc.5 annotation"))
pdf(file.path(out_dir, 'fibro/sc.5 annotation', str_glue("Dimplot dim{ndims_fibro} fibro annotation.pdf")), height = 6, width = 8)
DimPlot(sp_fibro, group.by = "annotation_fibro", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_fibro, group.by = "annotation_fibro", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_fibro, file.path(out_dir, 'sp_fibro.rds'))


# combined annotation ---------
# meta_t_anno <- 
#   sp_t[['annotation_t']] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(
#     sp_cd4t[['annotation_cd4t']] %>% rownames_to_column('cell_id'),
#     by = "cell_id"
#   ) %>%
#   left_join(
#     sp_cd8t[['annotation_cd8t']] %>% rownames_to_column('cell_id'),
#     by = "cell_id"
#   ) %>%
#   mutate(annotation_t = case_when(
#     !is.na(annotation_cd4t) ~ annotation_cd4t,
#     !is.na(annotation_cd8t) ~ annotation_cd8t,
#     TRUE ~ annotation_t
#   )) %>%
#   column_to_rownames('cell_id') %>%
#   .['annotation_t']

# sub_cluster_df <-
#   rbind(
#     sp_fibro[['annotation_fibro']] %>% dplyr::rename(annotation = 1),
#     meta_t_anno                    %>% dplyr::rename(annotation = 1),
#     sp_mac[['annotation_mac']]     %>% dplyr::rename(annotation = 1)
#   ) %>% 
#   rownames_to_column("cell_id")

sub_cluster_df <- sp_epi[['annotation_epi']] %>% dplyr::rename(annotation = 1) %>% rownames_to_column("cell_id")

sub_cluster_df_tidy <-
  all_data[['celltype']] %>% 
  rownames_to_column("cell_id") %>% 
  left_join(sub_cluster_df, by = 'cell_id') %>%
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(is.na(annotation), celltype,  annotation)) %>% 
  dplyr::select(cell_id, annotation2) %>% 
  column_to_rownames("cell_id")

all_data_subcluster <- 
  all_data %>%
  AddMetaData(sub_cluster_df_tidy)

# factor
all_data_subcluster$annotation <-
  all_data_subcluster$celltype %>%
  fct_recode(
    Epithelial = 'Hepatocyte',
    Epithelial = 'Tumor cell'
  ) %>%
  fct_relevel(c(
    "Epithelial",
    "T", "B/Plasma cell", "Kupffer cell",
    "Hepatic stellate cell", "Liver sinusoidal endothelial cell", "Pericyte"
  ))

# remove others
others_cells <- all_data_subcluster$annotation2 %>% unique %>% str_subset("like|Others|low")
all_data_subcluster2 <- all_data_subcluster[, !all_data_subcluster$annotation2 %in% others_cells]


# order cels - factor
all_data_subcluster2$annotation <-
  all_data_subcluster2$celltype %>%
  fct_recode(
    Epithelial = 'Hepatocyte',
    Epithelial = 'Tumor cell'
  ) %>%
  fct_relevel(c(
    "Epithelial",
    "T", "B/Plasma cell", "Kupffer cell",
    "Hepatic stellate cell", "Liver sinusoidal endothelial cell", "Pericyte"
  ))


all_data_subcluster2$annotation2  <-
  all_data_subcluster2$annotation2 %>%
  factor(
    levels = c(
      "Hepatocyte", "Cholangiocyte",
      "Tumor", "Tumor PH", "Tumor LR", "Proliferating tumor",
      "T", "B/Plasma cell", "Kupffer cell",
      "Hepatic stellate cell", "Liver sinusoidal endothelial cell", "Pericyte"
    )
  )
                          

saveRDS(all_data_subcluster,  'all_data_subcluster.rds')
saveRDS(all_data_subcluster2, 'all_data_subcluster2.rds')

pdf(
  'dimplot of  all annotation.pdf',
  height = 6,
  width  = 8
)
all_data_subcluster2 %>% 
  DimPlot(group.by = "annotation", reduction = "umap") +
  ggsci::scale_color_igv()

all_data_subcluster2 %>% 
  DimPlot(group.by = "annotation2", reduction = "umap") +
  ggsci::scale_color_igv()
dev.off()

pdf(
  'dimplot of  all annotation - split annotation.pdf',
  height = 16,
  width  = 20
)
all_data_subcluster2 %>% 
  DimPlot(group.by = "annotation2", split.by = "annotation2", ncol = 4, reduction = "umap") +
  ggsci::scale_color_igv()
dev.off()



# combined annotation2 ---------
# meta_t_anno <-
#   sp_t[['annotation_t']] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(
#     sp_cd4t[['annotation_cd4t']] %>% rownames_to_column('cell_id'),
#     by = "cell_id"
#   ) %>%
#   left_join(
#     sp_cd8t[['annotation_cd8t']] %>% rownames_to_column('cell_id'),
#     by = "cell_id"
#   ) %>%
#   mutate(annotation_t = case_when(
#     !is.na(annotation_cd4t) ~ annotation_cd4t,
#     !is.na(annotation_cd8t) ~ annotation_cd8t,
#     TRUE ~ annotation_t
#   )) %>%
#   column_to_rownames('cell_id') %>%
#   .['annotation_t']

sub_cluster_df <-
  rbind(
    sp_t[[c('cell_id_sample', 'annotation_t')]]         %>% dplyr::rename(annotation = 2),
    sp_epi[[c('cell_id_sample', 'annotation_epi')]]     %>% dplyr::rename(annotation = 2)
  )  

sub_cluster_df_sorted <-
  all_data[[c('cell_id_sample', "celltype")]] %>% 
  left_join(sub_cluster_df, by = "cell_id_sample") %>%
  mutate(annotation = ifelse(is.na(annotation), as.character(celltype), as.character(annotation)))

all_data$annotation <- sub_cluster_df_sorted$annotation

table( all_data$celltype, all_data$annotation )

saveRDS(all_data, 'all_data.epi.t.rds')



