library(tidyverse)
library(Seurat)
library(patchwork)
source("subcluster/sc_analysis_tools.R")

# multi workers  -----
library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 50 * 1024 ^ 3)
nbrOfFreeWorkers()



all_data <- readRDS("all_data.rds")

out_dir <- "sp_subcluster_LUSC"
dir.create(out_dir)

all_data$annotation_full %>% table
# AT2 cell                    B          Endothelial          Epithelials           Fibroblast 
# 44199                26739                70207               462241               139348 
# Mast cell               T cell marcophage_monocytes               plasma 
# 9849               116559               122606                71003 

major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  "Dendritic cell" = "FCER1A",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A"
)


# T subcluster ----------------
# sp_t <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "T cell"] ]
# sp_t$nCount_RNA <- sp_t$total_counts
# sp_t$project <- sp_t$sample

sp_t <- all_data[, all_data$annotation_full == "T"]

DefaultAssay(sp_t) <- 'RNA'

sp_t <- sc_normalize(sp_t, scale.factor = "median")
sp_t <- sc_pca(
  sp_t, 
  sc_out_dir = file.path(out_dir, "t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_t <- 10

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

all_markers_t <- 
  sc_findmarker(
    sp_t,
    idents =  "RNA_snn_res.0.6", # "sketch_snn_res.0.2",
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
  idents = "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_t[[c('RNA_snn_res.0.6', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.6)) +
  geom_col()

sc_mark_heatmap(
  sp_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

t_cell_cluster <- c(
  '0'   = "CD4T",
  '1'   = "CD8T",
  '2'   = "Macrophage-like T cell",
  '3'   = "CD4T",
  '4'   = "CD8T",
  '5'   = "Fibroblast-like T cell",
  '6'   = "Epithelial-like T cell",
  '7'   = "Cycling T cell",
  '8'   = "Endothelial-like T cell",
  '9'   = "B-like T cell",
  '10'  = "DC-like T cell"
)



Idents(sp_t) <- "RNA_snn_res.0.6"
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

sp_cd4t <- sc_normalize(sp_cd4t, scale.factor = "med")
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
    idents =  "RNA_snn_res.0.4", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


marker_de <-
  all_markers_cd4t %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_max(order_by = pct.1 - pct.2, n = 30) %>%
  slice_max(order_by = avg_log2FC, n = 5)

sp_cd4t <- sc_dotplot(
  sp_cd4t, 
  ndims = ndims_cd4t,
  marker = marker_de$gene,
  sc_out_dir = file.path(out_dir, "cd4t/sc.4 dotplot - de marker"), 
  dir.create.force = TRUE
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


sp_cd4t[[c('RNA_snn_res.0.6', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.6) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.6)) +
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
  '0'  = "CXCL13 CD8T",
  '1'  = "CD4 Treg",
  '2'  = "Others",
  '3'  = "Fibroblast-like T cell"
)



Idents(sp_cd4t) <- "RNA_snn_res.0.4"
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

sp_cd8t <- sc_normalize(sp_cd8t, scale.factor = "med")
sp_cd8t <- sc_pca(
  sp_cd8t, 
  sc_out_dir = file.path(out_dir, "cd8t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_cd8t <- 5

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
    "SELL", "CCR7", "IL7R", "TCF7",
    "MZB1","XBP1",
    "GZMA", "GZMH", "KLRG1", 
    "CD14","CD68","FCN1",
    "LAMP3",
    "XCR1",
    "CD1C",
    "CLEC9A",
    "TRDC", "TRAC",  
    "HAVCR2", "TIGIT", "LAG3"
  ),
  width = 22, 
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
  '0'  = "CD8 Temra",
  '1'  = "CD8 Tn",
  '2'  = "CD8 Tn",
  '3'  = "CD8 Tex",
  '4'  = "Others"
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
# sp_mac <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "marcophage_monocytes"] ]
# sp_mac$nCount_RNA <- sp_mac$total_counts
# sp_mac$project <- sp_mac$sample

sp_mac <- all_data[, all_data$annotation_full == "Myeloid cell"]

DefaultAssay(sp_mac) <- 'RNA'

sp_mac <- sc_normalize(sp_mac, scale.factor = "med")
sp_mac <- sc_pca(
  sp_mac, 
  sc_out_dir = file.path(out_dir, "mac/sc.1 pca"),
  ndims = 50,
  dir.create.force = TRUE
)

ndims_mac <- 30

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
  'CD19', 'MS4A1',
  'MZB1', 'XBP1',
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
  "Mast cell" = "KIT",
  "Mast cell" = "MS4A2",
  
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
  'Dendritic cell' = "CLEC10A",
  'Dendritic cell' = "FCER1A",
  
  'TOP2A',
  'MKI67',
  'TREM2',
  
  # 'CXCR4',
  'MMP12',
  'MMP9',
  'CCL18', 'CYP27A1',
  'MRC1',
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
    idents =  "RNA_snn_res.0.6", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )

write.csv(
  all_markers_mac,
  file = file.path(out_dir, 'mac', 'de analysis of RNA_snn_res.0.6.csv')
)

VlnPlot(
  sp_mac,
  group.by = "RNA_snn_res.0.6",
  features = "nCount_RNA",
  pt.size = 0
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

DotPlot(
  sp_mac,
  group.by = "RNA_snn_res.0.6",
  features = c(  'CD68',  'CD14', 
                 'TOP2A','MKI67','TREM2',
                 'MMP12','MMP9',
                 'CCL18', 'CYP27A1',
                 'MRC1','CXCL10'),
  
)


# mac_cell_cluster <- c(
#   '0'  = "low quality cells", # 无差异基因，nCount显著低于其他群
#   '1'  = "CYP27A1 Macrophage",
#   '2'  = "Fibroblast-like Macrophage",
#   '3'  = "Epithelial-like Macrophage",
#   '4'  = "Macro Macrophage", # 组织驻留型
#   '5'  = "Proliferating Macrophage",
#   '6'  = "CXCL10 Macrophage", # INF response
#   '7'  = "cDC2",
#   '8'  = "cDC1"
# )

mac_cell_cluster <- c(
  '0'  = "CCL18 Macrophage", # M2
  '1'  = "Fibroblast-like myeloid cell",
  '2'  = "T-like myeloid cell",
  '3'  = "Epithelial-like myeloid cell",
  '4'  = "MMP12 Macrophage", # M1
  '5'  = "low quality cells", # 无差异基因，nCount显著低于其他群
  '6'  = "CXCL10 Macrophage", # IFN response
  '7'  = "cDC2",
  '8'  = "Plasma-like myeloid cell",
  '9'  = "B-like myeloid cell",
  '10' = "Mast cell",
  '11' = "LAMP3 DC",
  '12' = "cDC1",
  '13' = "MMP9 Macrophage"
)


Idents(sp_mac) <- "RNA_snn_res.0.6"
sp_mac <- RenameIdents(sp_mac, mac_cell_cluster)
sp_mac$annotation_mac <- Idents(sp_mac)

dir.create(file.path(out_dir, "mac/sc.5 annotation"))
pdf(file.path(out_dir, 'mac/sc.5 annotation', str_glue("Dimplot dim{ndims_t} mac annotation.pdf")), height = 6, width = 8)
DimPlot(sp_mac, group.by = "annotation_mac", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mac, group.by = "annotation_mac", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_mac, file.path(out_dir, 'sp_mac.rds'))

# Epi subcluster ----------------
# sp_epi <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "Epithelials"] ]
# sp_epi$nCount_RNA <- sp_epi$total_counts
# sp_epi$project <- sp_epi$sample

sp_epi <- all_data[, all_data$annotation_full == "Epithelial"]

DefaultAssay(sp_epi) <- 'RNA'

sp_epi <- sc_normalize(sp_epi, scale.factor = "median")
sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

ElbowPlot(sp_epi, ndims = 50)

ndims_epi <- 15

sp_epi <- sc_dr(
  sp_epi, 
  ndims = ndims_epi,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "epi/sc.2 dr"),
  dir.create.force = TRUE
)


sp_epi <- sc_cluster(
  sp_epi, 
  ndims = ndims_epi,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "epi/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

saveRDS(sp_epi, file.path(out_dir, 'sp_epi.rds'))

# saveRDS(
#   sp_epi[[c('RNA_snn_res.0.2', 'RNA_snn_res.0.4', 'RNA_snn_res.0.6', 'RNA_snn_res.0.8')]],
#   file = file.path("infercnv", "sp_epi_meta.rds")
# )


epi_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  AT1 = 'HOPX',
  AT1 = 'AGER',
  AT1 = 'RTKN2',
  
  AT2 = 'SFTPC',
  AT2 = 'LAMP3',
  
  Goblet = 'MUC5AC',
  # Goblet = 'SPDEF',
  
  'Mucous cells' = 'SPDEF',
  'Mucous cells' = 'MUC5B',
  
  Tuft = 'DCLK1',
  Tuft = 'ASCL2',
  Tuft = 'POU2F3',
  
  'Ciliated cells' = 'FOXJ1',
  'Ciliated cells' = 'RSPH1'
  
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
  dir.create.force = TRUE
)


epi_markers2 <- c(
  'EPCAM',
  'TP63',
  'SOX2',
  'IGFBP2',
  'WNT5A',
  'CD44',
  
  'NKX2-1',
  'MUC1',
  'HOPX',
  'ERBB3',
  'SDC4',
  
  'TOP2A',
  'MKI67',
  'CCNB1',
  'CDK1',
  
  'SOX2-OT'
)


sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers2,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - epi - de"), 
  dir.create.force = TRUE
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)



all_markers_epi <- 
  sc_findmarker(
    sp_epi,
    idents = "RNA_snn_res.0.2",
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
    file = file.path(out_dir, "epi/de analysis - RNA_snn_res.0.2.csv")
  ) 


sp_epi[[c("group", "RNA_snn_res.0.2")]] %>%
  group_by(group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x  = group, y = freq, fill  = RNA_snn_res.0.2)) +
  geom_col()


VlnPlot(
  sp_epi,
  features = "nCount_RNA",
  group.by = "RNA_snn_res.0.2",
  pt.size = 0
)


DotPlot(
  sp_epi,
  features = c(
    "EPCAM", "CD68", "CD3E", "CD18", "PECAM1", "CCDC80",
    'APOL1','SLC2A1','TYMP','ALG11','APOL2', # 0
    'BIRC5','RRM1','SMC4','TUBB','NTS', # 1
    'CEBPA','CES1','CLDN1','NDRG4','GJA1', # 2 
    'CDKN1A', 'PLAU', 'SDC4', 'A2ML1', 'FUT3'  #  3
  ),
  group.by = "RNA_snn_res.0.2",
  scale = T
) +
  RotatedAxis()



dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


pdf(
  file.path(out_dir, "epi/dimplot of nCount.pdf"),
  height = 4,
  width  = 8
)
VlnPlot(
  sp_epi,
  group.by = 'RNA_snn_res.0.2',
  features = "nCount_RNA",
  log = TRUE,
  pt.size = 0
)
dev.off()

sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

epi_cell_cluster <- c(
  '0'  =  "ALG11 Epithelial",
  '1'  =  "BIRC5 Epithelial",
  '2'  =  "CEBPA Epithelial",
  '3'  =  "CDKN1A Epithelial"
)

DotPlot(
  sp_epi,
  features = c(
    'APOL1', 'SLC2A1', 'TYMP', 'APOL2',
    'ALG11', 'BIRC5', 'CEBPA', 'CDKN1A'
  ),
  group.by = 'RNA_snn_res.0.2'
)

Idents(sp_epi) <- "RNA_snn_res.0.2"
sp_epi <- RenameIdents(sp_epi, epi_cell_cluster)
sp_epi$annotation_epi <- Idents(sp_epi)

dir.create(file.path(out_dir, "epi/sc.5 annotation"))
pdf(file.path(out_dir, 'epi/sc.5 annotation', str_glue("Dimplot dim{ndims_epi} epi annotation .pdf")), height = 6, width = 8)
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_epi, group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # project  data
# DefaultAssay(sp_epi) <- "sketch"
# sp_epi[['full.harmonyumap']] <- NULL
# sp_epi <- 
#   sp_epi %>% 
#   ProjectData(
#     assay = "RNA",
#     full.reduction = "harmony.full",
#     
#     sketched.assay = "sketch",
#     sketched.reduction = "harmony", # "pca",
#     
#     umap.model = "harmonyumap", 
#     dims = 1:ndims_epi, 
#     
#     refdata = list(
#       annotation_full = "annotation",
#       RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
#       RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
#       RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
#       RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
#     )
#   )
# 
# Idents(sp_epi) <- 'annotation_full'
# DefaultAssay(sp_epi) <- "RNA"
# 
# pdf(file.path(out_dir, 'epi', str_glue("Dimplot dim{ndims_epi} annotation - projected.pdf")), height = 6, width = 8)
# DimPlot(sp_epi, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
# dev.off()

# saveRDS(sp_epi, file.path(out_dir, 'sp_epi.rds'))




# fibro subcluster ----------------
# sp_fibro <- all_data_raw[, colnames(all_data)[all_data$annotation_full == "Fibroblast"] ]
# sp_fibro$nCount_RNA <- sp_fibro$total_counts
# sp_fibro$project <- sp_fibro$sample
sp_fibro <- all_data[, all_data$annotation_full == "Fibroblast"]


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
    idents =  "RNA_snn_res.0.6", # "sketch_snn_res.0.2",
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

DotPlot(
  sp_fibro,
  features = c('CCDC80', 'POSTN', 'PDGFRA'),
  group.by = "RNA_snn_res.0.6",
  scale = T
)

# fibro_cell_cluster <- c(
#   '0'  = "mCAF",
#   '1'  = "infCAF",
#   '2'  = "dCAF",
#   '3'  = "NF", # normal fibroblast
#   '4'  = "Plasma-like fibroblast",
#   '5'  = "Pericyte",
#   '6'  = "CDH1 tCAF",
#   '7'  = "SOX2 tCAF",
#   '8'  = "Others"
# )

fibro_cell_cluster <- c(
  '0'  = "mCAF",
  '1'  = "mCAF",
  '2'  = "mCAF",
  '3'  = "T-like fibroblast",
  '4'  = "mCAF",
  '5'  = "Plasma-like fibroblast",
  '6'  = "Myeloid-like fibroblast",
  '7'  = "tCAF",
  '8'  = "infCAF",
  '9'  = "NF", # normal fibroblast
  '10'  = "B-like fibroblast",
  '11'  = "dCAF"
)


Idents(sp_fibro) <- "RNA_snn_res.0.6"
sp_fibro <- RenameIdents(sp_fibro, fibro_cell_cluster)
sp_fibro$annotation_fibro <- Idents(sp_fibro)

dir.create(file.path(out_dir, "fibro/sc.5 annotation"))
pdf(file.path(out_dir, 'fibro/sc.5 annotation', str_glue("Dimplot dim{ndims_fibro} fibro annotation.pdf")), height = 6, width = 8)
DimPlot(sp_fibro, group.by = "annotation_fibro", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_fibro, group.by = "annotation_fibro", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_fibro, file.path(out_dir, 'sp_fibro.rds'))


# B/Plasma -----
sp_b <- all_data[, all_data$annotation_full == "B/Plasma"]

DefaultAssay(sp_b) <- 'RNA'

sp_b <- sc_normalize(sp_b, scale.factor = "median")
sp_b <- sc_pca(
  sp_b, 
  sc_out_dir = file.path(out_dir, "b/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_b <- 10

sp_b <- sc_dr(
  sp_b, 
  ndims = ndims_b,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "b/sc.2 dr"),
  dir.create.force = TRUE
)


sp_b <- sc_cluster(
  sp_b, 
  ndims = ndims_b,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "b/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

b_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CCDC80",
  "PECAM1",
  "CD68",
  "CD14",
  "LILRA4",
  "CD1C",
  "IDO1",
  "XCR1",
  "CLEC9A",
  "CLEC10A",
  
  "CD19",
  "MS4A1",
  "CD79A",
  
  'CD27',
  
  
  "TNFRSF17",
  "JCHAIN",
  "MZB1",
  "XBP1"
)

sp_b <- sc_dotplot(
  sp_b, 
  ndims = ndims_b,
  marker = b_markers,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot"), 
  dir.create.force = TRUE
)

sp_b <- sc_dotplot(
  sp_b, 
  ndims = ndims_b,
  marker = major_marker_ls,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)



DotPlot(
  sp_b,
  features = c("CCDC80", "EPCAM", "CD3E", "CD68", "CD14", "MKI67", "TOP2A", "MS4A1", "CD27", "MZB1", "XBP1"),
  group.by = "RNA_snn_res.0.2"
)



all_markers_b <- 
  sc_findmarker(
    sp_b,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sp_b, all_markers = all_markers_b, sort_by = "avg")
sc_findmarker(sp_b, all_markers = all_markers_b, sort_by = "pct")



dir.create(file.path(out_dir, "b/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "b/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_b,
  markers = all_markers_b,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sp_b,
  markers = all_markers_b,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

b_cell_cluster <- c(
  '0'  =  "Plasma",
  '1'  =  "Plasma",
  '2'  =  "B cell"
)



Idents(sp_b) <- "RNA_snn_res.0.2"
sp_b <- RenameIdents(sp_b, b_cell_cluster)
sp_b$annotation_b <- Idents(sp_b)

dir.create(file.path(out_dir, "b/sc.5 annotation"))
pdf(file.path(out_dir, 'b/sc.5 annotation', str_glue("Dimplot dim{ndims_b} b annotation .pdf")), height = 6, width = 8)
DimPlot(sp_b, group.by = "annotation_b", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_b, group.by = "annotation_b", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sp_b, file.path(out_dir, "sp_b.rds"))


# combined annotation ---------
meta_t_anno <- 
  sp_t[['annotation_t']] %>%
  rownames_to_column('cell_id') %>%
  left_join(
    sp_cd4t[['annotation_cd4t']] %>% rownames_to_column('cell_id'),
    by = "cell_id"
  ) %>%
  left_join(
    sp_cd8t[['annotation_cd8t']] %>% rownames_to_column('cell_id'),
    by = "cell_id"
  ) %>%
  mutate(annotation_t = case_when(
    !is.na(annotation_cd4t) ~ annotation_cd4t,
    !is.na(annotation_cd8t) ~ annotation_cd8t,
    TRUE ~ annotation_t
  )) %>%
  column_to_rownames('cell_id') %>%
  .['annotation_t']

sub_cluster_df <-
  rbind(
    sp_fibro[['annotation_fibro']] %>% dplyr::rename(annotation = 1),
    meta_t_anno                    %>% dplyr::rename(annotation = 1),
    sp_mac[['annotation_mac']]     %>% dplyr::rename(annotation = 1),
    sp_epi[['annotation_epi']]     %>% dplyr::rename(annotation = 1),
    sp_b[['annotation_b']]         %>% dplyr::rename(annotation = 1)
    
  ) %>% 
  rownames_to_column("cell_id")

sub_cluster_df_tidy <-
  all_data[['annotation_full']] %>% 
  rownames_to_column("cell_id") %>% 
  left_join(sub_cluster_df, by = 'cell_id') %>%
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(is.na(annotation), annotation_full,  annotation)) %>% 
  dplyr::select(cell_id, annotation2) %>% 
  column_to_rownames("cell_id")

all_data_subcluster <- 
  all_data %>%
  AddMetaData(sub_cluster_df_tidy)

others_cells <- all_data_subcluster$annotation2 %>% unique %>% str_subset("like|Others|low")
all_data_subcluster2 <- all_data_subcluster[, !all_data_subcluster$annotation2 %in% others_cells]

saveRDS(all_data_subcluster,  'all_data_subcluster.rds')
saveRDS(all_data_subcluster2, 'all_data_subcluster2.rds')

all_data_subcluster2 %>% 
  DimPlot(group.by = "annotation2", reduction = "harmonyumap") +
  ggsci::scale_color_igv()









# spatial plot used for tumor cells --------------
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
parse_data_path <- readRDS("data_for_analysis/parse_data_path.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

sp_out_dir_epi <- file.path(out_dir, 'epi/spatial plot') 
dir.create(sp_out_dir_epi, recursive = TRUE)

meta_sp_epi <- 
  sp_epi[['annotation_epi']] %>%
  # dplyr::rename(annotation = annotation_epi) %>%
  rownames_to_column('cell_id') 

all_data_meta_epi <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%

  left_join(meta_sp_epi, by = "cell_id") %>%
  mutate(annotation2 = ifelse(is.na(annotation_epi), as.character(annotation_full), as.character(annotation_epi)))


sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){#  browser()
    cell_cols <- c(
      "ALG11 Epithelial"  = "#E41A1C",
      "BIRC5 Epithelial"  = "#377EB8",
      "CEBPA Epithelial"  = "#4DAF4A",
      "CDKN1A Epithelial" = "#984EA3",
      
      "T cell"         = "gray",
      "B"              = "gray",
      "plasma"         = "gray",
      "Mast cell"      = "gray",
      "marcophage_monocytes" = "gray",
      
      "AT2 cell"       = "gray",
      "Endothelial"    = "gray",
      "Fibroblast"     = "gray"
    )
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta_epi, by = c("cell_id"))
    
    group_info <- all_data_meta_epi$group[all_data_meta_epi$project == project] %>%  .[1]
    
    p <- dat %>%
      mutate(annotation2 = factor(annotation2, levels = names(cell_cols))) %>%
      
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_cols) +
      theme_void() +
      labs(fill = NULL , title = paste(project, group_info, sep = "_") )
    p
  })

# export image
sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir_epi, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })






