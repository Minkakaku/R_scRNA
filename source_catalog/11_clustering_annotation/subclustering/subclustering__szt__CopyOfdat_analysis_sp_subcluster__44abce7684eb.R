library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")

all_data <- readRDS('output/data/adata_filtered_add_celltype.rds')

out_dir <- "sp_subcluster"
dir.create(out_dir)

all_data$annotation_full <- all_data$annotation
all_data$annotation_full %>% table
# B              Endothelial Epithelial proliferative             Epithelial_1 
# 89515                   119845                   328204                    61475 
# Epithelial_2               Fibroblast             Fibroblast-T        Immune-Epithelial 
# 347381                   538602                    48274                    30859 
# Macrophage                 Pericyte                   Plasma                        T 
# 169641                   123042                    80867                   212926 


# Epi subcluster ----------------
sp_epi <- all_data[, all_data$annotation_full %in% c("Epithelial_1", "Epithelial_2", "Epithelial proliferative")]

DefaultAssay(sp_epi) <- 'RNA'

sp_epi <- sc_normalize(sp_epi, scale.factor = "median")

# ## sketch data --------
# sp_epi <-
#   SketchData(
#     sp_epi,
#     ncells         = 50000L,
#     sketched.assay = "sketch",
#     method         = "LeverageScore",
#     var.name       = "leverage.score",
#     seed = 123L
#   )


sp_epi <- sc_pca(
  sp_epi, 
  sc_out_dir = file.path(out_dir, "epi/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_epi <- 20

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

sp_epi2 <- sp_epi[, sample(colnames(sp_epi), 5000)]
library(clustree)
clustree::clustree(
  sp_epi,
  prefix = "RNA_snn_res."
)

epi_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CD8A",
  
  Enterocyte  = "SLC26A3",
  
  Goblet = "MUC2",
  Goblet = "TFF3",
  
  
  # Enteroendocrine = "CHGA",
  Enteroendocrine = "TPH1",
  Enteroendocrine = "NEUROD1",
  
  'Paneth cell' = 'DEFA6',
  'Paneth cell' = 'IL9',
  
  'Transit amplifying cell' = 'CCNA2',
  'Transit amplifying cell' = 'HMGB3',
  'Transit amplifying cell' = 'MCM5',
  'Transit amplifying cell' = 'PCNA',
  'Transit amplifying cell' = 'TOP2A',
  
  'Proximal Mature Enterocytes' = 'FABP2',
  'Proximal Mature Enterocytes' = 'ACE',
  'Proximal Mature Enterocytes' = 'FABP6',
  'Proximal Mature Enterocytes' = 'APOC3',
  
  'Proximal Enterocytes' = 'HEBP1',
  'Proximal Enterocytes' = 'CCL25'
)

sp_epi <- sc_dotplot(
  sp_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = file.path(out_dir, "epi/sc.4 dotplot"), 
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
    idents = "RNA_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sp_epi, all_markers = all_markers_epi, sort_by = "pct")



dir.create(file.path(out_dir, "epi/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "epi/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
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
  '0'  =  "Proximal mature enterocytes",
  '1'  =  "Peneth cell",
  '2'  =  "Transit amplifying cell",
  '3'  =  "Others",
  '4'  =  "Enterocytes",
  '5'  =  "Others",
  '6'  =  "Others",
  '7'  =  "CD8T",
  '8'  =  "Transit amplifying cell",
  '9'  =  "Others",
  '10'  = "Others",
  '11'  = "Others",
  '12'  = "Others",
  '13'  = "Proximal mature enterocytes",
  '14'  = "Others"

)



Idents(sp_epi) <- "RNA_snn_res.0.4"
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
sp_fibro <- all_data[, all_data$annotation_full == "Fibroblast"]

DefaultAssay(sp_fibro) <- 'RNA'

sp_fibro <- sc_normalize(sp_fibro, scale.factor = "def")
sp_fibro <- sc_pca(
  sp_fibro, 
  sc_out_dir = file.path(out_dir, "fibro/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_fibro <- 15

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
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "fibro/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

fibro_markers <- c(
  # Enteroendocrine = "CHGA",
  Enteroendocrine = "TPH1",
  Enteroendocrine = "NEUROD1",
  
  'Paneth cell' = 'DEFA6',
  'Paneth cell' = 'IL9',
  
  'Transit amplifying cell' = 'CCNA2',
  'Transit amplifying cell' = 'HMGB3',
  'Transit amplifying cell' = 'MCM5',
  'Transit amplifying cell' = 'PCNA',
  'Transit amplifying cell' = 'TOP2A',
  
  'Proximal Mature Enterocytes' = 'FABP2',
  'Proximal Mature Enterocytes' = 'ACE',
  'Proximal Mature Enterocytes' = 'FABP6',
  'Proximal Mature Enterocytes' = 'APOC3',
  
  'Proximal Enterocytes' = 'HEBP1',
  'Proximal Enterocytes' = 'CCL25'
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
  marker = major_marker_for_subcluster,
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
  '0'  = "COL5A1 Fibroblast",
  '1'  = "TNXB Fibroblast",
  '2'  = "ADAMDEC1 Fibroblast",
  '3'  = "THBS1 Fibroblast",
  '4'  = "EEF1G Fibroblast"
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

# T subcluster ----------------
sp_t <- all_data[, all_data$annotation_full == "T"]

DefaultAssay(sp_t) <- 'RNA'

sp_t <- sc_normalize(sp_t, scale.factor = "def")
sp_t <- sc_pca(
  sp_t, 
  sc_out_dir = file.path(out_dir, "t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_t <- 20

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
  marker = t_markers,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_t <- sc_dotplot(
  sp_t, 
  ndims = ndims_t,
  marker = major_marker_for_subcluster,
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
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_t[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
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
  '0'  = "Others",
  '1'  = "CD4T naive",
  '2'  = "GZMK CD8T",
  '3'  = "CD4T",
  '4'  = "Tregs",
  '5'  = "Plasma",
  '6'  = "CD4T",
  '7'  = "Enterocyte-like", # Eppithelial
  '8'  = "CD8T",
  '9'  = "MYLK CD8T",
  '10' = "Others",
  '11' = "CD8T",
  '12' = "CD8T",
  '13' = "CD8T",
  '14' = "Others"
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

# macrophage subcluster --------------

sp_mac <- all_data[, all_data$annotation_full == "Macrophage"]

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
  'Dendritic cell' = "FCER1A"
  
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
  marker = major_marker_for_subcluster,
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
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_mac[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
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
  '0'  = "CCDC80 Macrophage",
  '1'  = "IL1B Monocyte",
  '2'  = "CD163 Macrophage",
  '3'  = "Others",
  '4'  = "Macrophage",
  '5'  = "IDO1 DC",
  '6'  = "CD14 Monocyte",
  '7'  = "Plasma",
  '8'  = "CD8T",
  '9'  = "CD68 Macrophage",
  '10' = "CLEC9A DC",
  '11' = "CD83 DC",
  '12' = "Others",
  '13' = "Others"
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


# combined annotation ---------

sub_cluster_df <-
  rbind(
    sp_epi[['annotation_epi']]     %>% dplyr::rename(annotation = 1),
    sp_fibro[['annotation_fibro']] %>% dplyr::rename(annotation = 1),
    sp_t[['annotation_t']]         %>% dplyr::rename(annotation = 1),
    sp_mac[['annotation_mac']]     %>% dplyr::rename(annotation = 1)
  ) %>% 
  rownames_to_column("cell_id")

sub_cluster_df_tidy <-
  all_data[['annotation_full2']] %>% 
  rownames_to_column("cell_id") %>% 
  left_join(sub_cluster_df, by = 'cell_id') %>%
  map_dfc(as.character) %>%
  mutate(annotation3 = ifelse(is.na(annotation), annotation_full2,  annotation)) %>% 
  dplyr::select(cell_id, annotation3) %>% 
  column_to_rownames("cell_id")

all_data <- 
  all_data %>%
  AddMetaData(sub_cluster_df_tidy)


# saveRDS(all_data, 'all_data_subcluster.rds')

all_data %>% 
  DimPlot(group.by = "annotation3", reduction = "harmonyumap") +
  ggsci::scale_color_igv()
















