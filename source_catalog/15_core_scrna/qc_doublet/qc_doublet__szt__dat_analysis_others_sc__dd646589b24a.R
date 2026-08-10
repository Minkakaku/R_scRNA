library(Seurat)
library(tidyverse)
source("code base/sc_analysis_tools.R")

all_data_other <- readRDS("others_sc/GSE209664_Sampath_finalobject-3.rds")


all_data_other %>% 
  DimPlot(
    group.by = c("lineage", "cell_type"),
    label = TRUE
  ) +
  NoLegend()

DotPlot(
  all_data_other,
  features = c(
    "Cd3d", "Cd3e", "Cd4", "Cd8a",
    "Trac", "Il7r", "Trdc", "Nkg7"
  ),
  group.by = "cell_type"
)


#### NK/T cells subcluster ----------------
sc_t_other <- all_data_other[, all_data_other$cell_type == "T cells"]

DefaultAssay(sc_t_other) <- "RNA"
sc_t_other[['SCT']] <- NULL
sc_t_other[['integrated']] <- NULL

sc_t_other$project <- sc_t_other$individual

sc_t_other <- sc_normalize(sc_t_other)
sc_t_other <- sc_pca(
  sc_t_other, 
  sc_out_dir = "others_sc/GSE209664/subcluster/t/sc.1 pca",
  dir.create.force = TRUE
)

ndims_t_other <- 20

sc_t_other <- sc_dr(
  sc_t_other, 
  ndims = ndims_t_other,
  sc_out_dir = "others_sc/GSE209664/subcluster/t/sc.2 dr",
  dir.create.force = TRUE
)


sc_t_other <- sc_cluster(
  sc_t_other, 
  ndims = ndims_t_other,
  reduction =  "pca", # "harmony",
  sc_out_dir = "others_sc/GSE209664/subcluster/t/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

t_markers_other <- c(
  "Cd3d",
  "Cd3e",
  "Cd4",
  "Cd8a",
  "Foxp3",
  "Il7r",
  "Trac",
  "Trdc",
  "Ncam1",
  "Nkg7",
  
  'Mki67',
  'Top2a'
  
  # 'Malat1','H3f3b','Trbc2','Eif1','Satb1','H3f3a','Rbm39','Ubb','Sftpc',
  # 'Mier1',
  # 
  # 
  # 'Hbb-bt', 
  # 'Hbb-bs'
)

DotPlot(
  sc_t_other,
  group.by = "RNA_snn_res.0.2",
  features = c(
    "Cd19", "Cd20", "Ighd",
    "Epcam", "Pecam1",
    "Col1a1", "Dcn",
    "Acta2", "Vim",
    "Ager",
    "Sftpa1",
    "Sftpd",
    ""
  )
)






sc_t_other <- sc_dotplot(
  sc_t_other, 
  ndims = ndims_t_other,
  marker = t_markers_other,
  sc_out_dir = 'others_sc/GSE209664/subcluster/t/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 25
)

all_markers_t <- 
  sc_findmarker(
    sc_t_other,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

sc_findmarker(sc_t_other, all_markers = all_markers_t, sort_by = "pct")

write.csv(all_markers_t, file = file.path("others_sc/GSE209664/subcluster/t", "all_markers.csv"))

sc_mark_heatmap(
  sc_t_other,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_t_other,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


t_cell_cluster_other <- c(
  '0'   = "CD8T",
  '1'   = "gd T cells",
  '2'   = "ILC2",
  '3'   = "CD4T cells",
  '4'   = "NK cells",
  '5'   = "Treg cells",
  '6'   = "CD4-gd T cells",
  '7'   = "NK",
  '8'   = "B naive"
)


Idents(sc_t_other) <- "RNA_snn_res.0.2"
sc_t_other <- RenameIdents(sc_t_other, t_cell_cluster_other)
sc_t_other$annotation <- Idents(sc_t_other)

dir.create('others_sc/GSE209664/subcluster/t/sc.5 annotation')
pdf(file.path('others_sc/GSE209664/subcluster/t/sc.5 annotation', str_glue("Dimplot dim{ndims_t_other} harmony t annotation.pdf")), height = 6, width = 8)
DimPlot(sc_t_other, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t_other, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_t_other[, sc_t_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_t_other[, sc_t_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_t_other, 'others_sc/GSE209664/subcluster/sc_t_other.rds')




#### AT1 cells subcluster ----------------
sc_at1_other <- all_data_other[, all_data_other$cell_type == "AT1"]

DefaultAssay(sc_at1_other) <- "RNA"
sc_at1_other[['SCT']] <- NULL
sc_at1_other[['integrated']] <- NULL

sc_at1_other$project <- sc_at1_other$individual

sc_at1_other <- sc_normalize(sc_at1_other)
sc_at1_other <- sc_pca(
  sc_at1_other, 
  sc_out_dir = "others_sc/GSE209664/subcluster/at1/sc.1 pca",
  dir.create.force = TRUE
)

ndims_at1_other <- 20

sc_at1_other <- sc_dr(
  sc_at1_other, 
  ndims = ndims_at1_other,
  sc_out_dir = "others_sc/GSE209664/subcluster/at1/sc.2 dr",
  dir.create.force = TRUE
)


sc_at1_other <- sc_cluster(
  sc_at1_other, 
  ndims = ndims_at1_other,
  reduction =  "pca", # "harmony",
  sc_out_dir = "others_sc/GSE209664/subcluster/at1/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

at1_markers_other <- c(
  "Cd3d",
  "Cd3e",
  "Cd4",
  "Cd8a",
  "Foxp3",
  "Il7r",
  "Trac",
  "Trdc",
  "Ncam1",
  "Nkg7",
  
  'Mki67',
  'Top2a',
  
  'Cldn3', 'S100a14', 
  'Lars2',
  'Mettl7a1',
  'Adgrg1', 'Nupr1' # 'man1c1
  
)

DotPlot(
  sc_at1_other,
  group.by = "RNA_snn_res.0.2",
  features = c(
    "Cd19", "Cd20", "Ighd",
    "Epcam", "Pecam1",
    "Col1a1", "Dcn",
    "Acta2", "Vim",
    "Ager",
    "Sftpa1",
    "Sftpd",
    ""
  )
)






sc_at1_other <- sc_dotplot(
  sc_at1_other, 
  ndims = ndims_at1_other,
  marker = at1_markers_other,
  sc_out_dir = 'others_sc/GSE209664/subcluster/at1/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 25
)

all_markers_at1 <- 
  sc_findmarker(
    sc_at1_other,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

sc_findmarker(sc_at1_other, all_markers = all_markers_at1, sort_by = "pct")

write.csv(all_markers_at1, file = file.path("others_sc/GSE209664/subcluster/at1", "all_markers.csv"))

sc_mark_heatmap(
  sc_at1_other,
  markers = all_markers_at1,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_at1_other,
  markers = all_markers_at1,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


at1_cell_cluster_other <- c(
  '0'   = "CD8T",
  '1'   = "gd T cells",
  '2'   = "ILC2",
  '3'   = "CD4T cells",
  '4'   = "NK cells",
  '5'   = "Treg cells",
  '6'   = "CD4-gd T cells",
  '7'   = "NK",
  '8'   = "B naive"
)

Idents(sc_at1_other) <- "RNA_snn_res.0.2"
sc_at1_other <- RenameIdents(sc_at1_other, t_cell_cluster)
sc_at1_other$annotation <- Idents(sc_at1_other)

dir.create('others_sc/GSE209664/subcluster/at1/sc.5 annotation')
pdf(file.path('others_sc/GSE209664/subcluster/at1/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony at1 annotation.pdf")), height = 6, width = 8)
DimPlot(sc_at1_other, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at1_other, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_at1_other[, sc_at1_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at1_other[, sc_at1_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_at1_other, 'others_sc/GSE209664/subcluster/sc_at1_other.rds')


#### AT2 cells subcluster ----------------
sc_at2_other <- all_data_other[, all_data_other$cell_type == "AT2"]

DefaultAssay(sc_at2_other) <- "RNA"
sc_at2_other[['SCT']] <- NULL
sc_at2_other[['integrated']] <- NULL

sc_at2_other$project <- sc_at2_other$individual

sc_at2_other <- sc_normalize(sc_at2_other)
sc_at2_other <- sc_pca(
  sc_at2_other, 
  sc_out_dir = "others_sc/GSE209664/subcluster/at2/sc.1 pca",
  dir.create.force = TRUE
)

ndims_at2_other <- 20

sc_at2_other <- sc_dr(
  sc_at2_other, 
  ndims = ndims_at2_other,
  sc_out_dir = "others_sc/GSE209664/subcluster/at2/sc.2 dr",
  dir.create.force = TRUE
)


sc_at2_other <- sc_cluster(
  sc_at2_other, 
  ndims = ndims_at2_other,
  reduction =  "pca", # "harmony",
  sc_out_dir = "others_sc/GSE209664/subcluster/at2/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

at2_markers_other <- c(
  "Cd3d",
  "Cd3e",
  "Cd4",
  "Cd8a",
  "Foxp3",
  "Il7r",
  "Trac",
  "Trdc",
  "Ncam1",
  "Nkg7",
  
  'Mki67',
  'Top2a',
  
  'Dlk1', 'Adamts1', 'Areg', 'Acoxl', 'Phgdh'
  
)

DotPlot(
  sc_at2_other,
  group.by = "RNA_snn_res.0.2",
  features = c(
    "Cd19", "Cd20", "Ighd",
    "Epcam", "Pecam1",
    "Col1a1", "Dcn",
    "Acta2", "Vim",
    "Ager",
    "Sftpa1",
    "Sftpd",
    ""
  )
)






sc_at2_other <- sc_dotplot(
  sc_at2_other, 
  ndims = ndims_at2_other,
  marker = at2_markers_other,
  sc_out_dir = 'others_sc/GSE209664/subcluster/at2/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 25
)

all_markers_at2 <- 
  sc_findmarker(
    sc_at2_other,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

sc_findmarker(sc_at2_other, all_markers = all_markers_at2, sort_by = "pct")

write.csv(all_markers_at2, file = file.path("others_sc/GSE209664/subcluster/at2", "all_markers.csv"))

sc_mark_heatmap(
  sc_at2_other,
  markers = all_markers_at2,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_at2_other,
  markers = all_markers_at2,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


at2_cell_cluster_other <- c(
  '0'   = "CD8T",
  '1'   = "gd T cells",
  '2'   = "ILC2",
  '3'   = "CD4T cells",
  '4'   = "NK cells",
  '5'   = "Treg cells",
  '6'   = "CD4-gd T cells",
  '7'   = "NK",
  '8'   = "B naive"
)

Idents(sc_at2_other) <- "RNA_snn_res.0.2"
sc_at2_other <- RenameIdents(sc_at2_other, t_cell_cluster)
sc_at2_other$annotation <- Idents(sc_at2_other)

dir.create('others_sc/GSE209664/subcluster/at2/sc.5 annotation')
pdf(file.path('others_sc/GSE209664/subcluster/at2/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony at2 annotation.pdf")), height = 6, width = 8)
DimPlot(sc_at2_other, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at2_other, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_at2_other[, sc_at2_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at2_other[, sc_at2_other$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_at2_other, 'others_sc/GSE209664/subcluster/sc_at2_other.rds')


### merge annotation
table(all_data_other$cell_type, all_data_other$CellType_finename)

meta_data_for_merge_other <-
  all_data_other[[c('CellType_finename', 'cell_type')]] %>%
  rownames_to_column("ID") %>%
  left_join(rownames_to_column(sc_t_other[['annotation']], "ID"), by = "ID") %>%
  mutate(annotation2 = ifelse(is.na(annotation), as.character(CellType_finename), as.character(annotation))) %>%
  dplyr::select(ID, annotation2)
  
  
all_data_other <- all_data_other %>% AddMetaData(meta_data_for_merge_other)

saveRDS(all_data_other, 'others_sc/GSE209664/all_data_other.rds')

### validation -----------
p_dim_other <-
  DimPlot(
    all_data_other,
    group.by = "annotation2",
    label = TRUE,
    repel = TRUE
  ) +
  NoLegend()

p_dim_other


p_dot_light_other  <-
  DotPlot(
    all_data_other,
    assay = "RNA",
    features = c('Tnfsf14', 'Ltbr'),
    group.by = "annotation2"
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  )


p_dot_light2_other <-
  DotPlot(
    all_data_other[, all_data_other$annotation2 == "ILC2"],
    assay = "RNA",
    features = c('Tnfsf14', 'Ltbr'),
    group.by = "condition" # "Oxygen" # Sample
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('ILC2')


# Vwf	ITGAV_ITGB3
p_dot_vwf_other <-
  DotPlot(
    all_data_other,
    assay = "RNA",
    features = c('Vwf', 'Itgav', 'Itgb3'),
    group.by = "annotation2"
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  )

p_dot_vwf2_other <-
  DotPlot(
    all_data_other[, all_data_other$annotation2 == "ILC2"],
    assay = "RNA",
    features = c('Vwf', 'Itgav', 'Itgb3'),
    group.by = "condition" # "Oxygen" # Sample
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('ILC2')


DotPlot(
  all_data_other,
  assay = "RNA",
  features = c('Vwf', 'Itgav', 'Itgb3'),
  group.by = "condition" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  )

# pdf(
#   file.path(out_dir_cellchat, str_glue("dotplot - LIGHT pathway.pdf")),
#   height = 8,
#   width  = 6
# )
# p_dot_light %>% print()
# p_dot_light2 %>% print()
# dev.off()
# 
# pdf(
#   file.path(out_dir_cellchat, str_glue("dotplot - VWF pathway.pdf")),
#   height = 8,
#   width  = 6
# )
# p_dot_vwf %>% print()
# p_dot_vwf2 %>% print()
# dev.off()

#
DotPlot(
  sc_t_other,
  group.by = "RNA_snn_res.0.2",
  features = c(
    "Tnfsf14",
    "Itgav",
    "Itgb3"
  )
)

DotPlot(
  all_data_other,
  group.by = "cell_type",
  features = c(
    "Tnfsf14",
    "Itgav",
    "Itgb3",
    'Ltbr'
  )
)


# hallmarkers ----
options(timeout = 10^3)
h_geneset <- msigdbr::msigdbr(collection = "H", species = "Mus musculus")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

sc_at1_score <-
  sc_at1_2 %>%
  JoinLayers() %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

sc_at1_score %>% head

# msigdb anno
h_at1_score_anno <- 
  sc_at1_score[[]] %>%
  dplyr::select(annotation, Cluster1:Cluster50) %>%
  group_by(annotation) %>%
  summarise(across(everything(), mean))

p_at1_hallmark_1 <-
  pheatmap::pheatmap(
    h_at1_score_anno %>% dplyr::filter(annotation != "doublet") %>% column_to_rownames("annotation") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )



