library(tidyverse)
library(Seurat)
library(patchwork)

source("subcluster/sc_analysis_tools.R")

all_data_sc <- readRDS("all_data_sc_dim20.rds")

out_dir <- "sc_subcluster"
dir.create(out_dir)

all_data_sc$annotation %>% table
# B      Endothelial    Epithelial  Fibroblast  Gial  Macrophage  Mast cell  Myefibroblast  NK    Others 
# 65393  108071         410131      200715      26781 182161      13356      102490         2102  21 

# Pericyte Plasma  T 
# 48124    201339  209290 



# peircyte/smc subcluster ----------------
sc_smc <- all_data_sc[, all_data_sc$annotation == "Pericyte/SMC"]

DefaultAssay(sc_smc) <- 'RNA'

sc_smc <- sc_normalize(sc_smc, scale.factor = "def")
sc_smc <- sc_pca(
  sc_smc, 
  sc_out_dir = file.path(out_dir, "smc/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_fibro <- 10

sc_smc <- sc_dr(
  sc_smc, 
  ndims = ndims_fibro,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "smc/sc.2 dr"),
  dir.create.force = TRUE
)


sc_smc <- sc_cluster(
  sc_smc, 
  ndims = ndims_fibro,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "smc/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

fibro_markers <- c(
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  # "Fibroblasts" = "DCN",
  # "Fibroblasts" = "ACTA2",
  "Fibroblasts" = "Postn",
  # "Mural cells" = "Pdgfrb",
  SMC = "Myh11",
  SMC = "Tagln",
  # Pdgfrb # Pericytes
  
  # "Pericytes" = "LFHP",
  "Pericytes" = "Rgs5",
  "Pericytes" = "Mcam"
)

sc_smc <- sc_dotplot(
  sc_smc, 
  ndims = ndims_fibro,
  marker = fibro_markers,
  sc_out_dir = file.path(out_dir, "smc/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sc_smc <- sc_dotplot(
  sc_smc, 
  ndims = ndims_fibro,
  marker = marker_ls,
  sc_out_dir = file.path(out_dir, "smc/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_fibro <- 
  sc_findmarker(
    sc_smc,
    idents =  "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )

sc_findmarker(sc_smc, all_markers = all_markers_fibro, sort_by = "avg")

dir.create(file.path(out_dir, "smc/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "smc/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sc_smc,
  markers = all_markers_fibro,
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_smc[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sc_smc,
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
  '0'  = "Pericyte",
  '1'  = "Others",
  '2'  = "Pericyte-like CM",
  '3'  = "Immune-like cell",
  '4'  = "EC-like cell",
  '5'  = "Fibroblast-like cell",
  '6'  = "SMC",
  '7'  = "Fibroblast-like proliferative"
)



Idents(sc_smc) <- "RNA_snn_res.0.2"
sc_smc <- RenameIdents(sc_smc, fibro_cell_cluster)
sc_smc$annotation_fibro <- Idents(sc_smc)

dir.create(file.path(out_dir, "smc/sc.5 annotation"))
pdf(file.path(out_dir, 'smc/sc.5 annotation', str_glue("Dimplot dim{ndims_fibro} fibro annotation.pdf")), height = 6, width = 8)
DimPlot(sc_smc, group.by = "annotation_fibro", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_smc, group.by = "annotation_fibro", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sc_smc, file.path(out_dir, 'sc_smc.rds'))


# macrophage subcluster --------------

sc_immune <- all_data_sc[, all_data_sc$annotation == "Immune cell"]

DefaultAssay(sc_immune) <- 'RNA'

sc_immune <- sc_normalize(sc_immune, scale.factor = "def")
sc_immune <- sc_pca(
  sc_immune, 
  sc_out_dir = file.path(out_dir, "mac/sc.1 pca"),
  ndims = 50,
  dir.create.force = TRUE
)

ndims_mac <- 10

sc_immune <- sc_dr(
  sc_immune, 
  ndims = ndims_mac,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "mac/sc.2 dr"),
  dir.create.force = TRUE
)


sc_immune <- sc_cluster(
  sc_immune, 
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
  'Dendritic cell' = "FCER1A",
  
  Neutrophil = "CSF3R",
  Neutrophil = "S100A8",
  Neutrophil = "S100A9"
)


DotPlot(
  sc_immune,
  features = c(
    'Fcn1', 'Adgre1', 'Spp1', 'Vcan', "Cd163", "Cd68",
    'Csf3r', 'S100a8', 'Fcgr3',
    'Cxcl10', 'Ifit1', 'Ifit2', 'Ifit3',
    'Siglech', 'Ly6c1', 'Ly6c2',
    'Cd3d', 'Cd3e', 'Cd3g', 'Cd4', 'Cd8a',
    'Ncam1'
  ),
  group.by = 'RNA_snn_res.0.2'
)
sc_immune <- sc_dotplot(
  sc_immune, 
  ndims = ndims_mac,
  marker = mac_markers %>% {structure(names = names(.), str_to_title(.))},
  sc_out_dir = file.path(out_dir, "mac/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sc_immune <- sc_dotplot(
  sc_immune, 
  ndims = ndims_mac,
  marker = c(
    'Fcn1', 'Adgre1', 'Spp1', 'Vcan', "Cd163", "Cd68",
    'Csf3r', 'S100a8', 'Fcgr3',
    'Cxcl10', 'Ifit1', 'Ifit2', 'Ifit3',
    'Siglech', 'Ly6c1', 'Ly6c2',
    'Cd3d', 'Cd3e', 'Cd3g', 'Cd4', 'Cd8a',
    'Ncam1'
  ),
  sc_out_dir = file.path(out_dir, "mac/sc.4 dotplot - selected marker"), 
  dir.create.force = TRUE
)
sc_immune <- sc_dotplot(
  sc_immune, 
  ndims = ndims_mac,
  marker = marker_ls,
  sc_out_dir = file.path(out_dir, "mac/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)


all_markers_mac <- 
  sc_findmarker(
    sc_immune,
    idents =  "RNA_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


sc_findmarker(sc_immune, all_markers = all_markers_mac, sort_by = "pct")
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
  sc_immune,
  markers = all_markers_mac,
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_immune[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sc_immune,
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
  '0'  = "Spp1 Macrophage",
  '1'  = "Cd163 Macrophage",
  '2'  = "Neutrophil",
  '3'  = "Others",
  '4'  = "Cxcl10 Macrophage",
  '5'  = "Others",
  '6'  = "Others"
)



Idents(sc_immune) <- "RNA_snn_res.0.2"
sc_immune <- RenameIdents(sc_immune, mac_cell_cluster)
sc_immune$annotation_mac <- Idents(sc_immune)

dir.create(file.path(out_dir, "mac/sc.5 annotation"))
pdf(file.path(out_dir, 'mac/sc.5 annotation', str_glue("Dimplot dim{ndims_mac} mac annotation.pdf")), height = 6, width = 8)
DimPlot(sc_immune, group.by = "annotation_mac", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "annotation_mac", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sc_immune, file.path(out_dir, 'sc_immune.rds'))


# combined annotation ---------

sub_cluster_df <-
  rbind(
    sc_smc[['annotation_fibro']] %>% dplyr::rename(annotation = 1),
    sc_immune[['annotation_mac']]     %>% dplyr::rename(annotation = 1)
  ) %>% 
  rownames_to_column("cell_id")

sub_cluster_df_tidy <-
  all_data_sc[['annotation']] %>% 
  rownames_to_column("cell_id") %>%
  dplyr::rename(annotation_full = annotation) %>%
  left_join(sub_cluster_df, by = 'cell_id') %>%
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(is.na(annotation), annotation_full,  annotation)) %>% 
  dplyr::select(cell_id, annotation2) %>% 
  column_to_rownames("cell_id")


others_cell <- 
  all_data_sc$annotation2 %>%
  unique() %>%
  str_subset(
    "like|Others|high-mt"
  )
all_data_sc2 <- all_data_sc[, !all_data_sc$annotation2 %in% others_cell]

# saveRDS(all_data_sc, 'all_data_sc_subcluster.rds')
# saveRDS(all_data_sc2, 'all_data_sc_subcluster2.rds')


# plot dimplot
pdf(
  file.path(out_dir, "dimploto of annotation - all.pdf"),
  height = 6,
  width = 12
)

all_data_sc %>% 
  DimPlot(group.by = "annotation2", reduction = "harmonyumap") +
  ggsci::scale_color_igv()
dev.off()

pdf(
  file.path(out_dir, "dimploto of annotation - remove others.pdf"),
  height = 6,
  width = 8
)
all_data_sc2 %>% 
  DimPlot(group.by = "annotation2", reduction = "harmonyumap") +
  ggsci::scale_color_igv()

dev.off()













