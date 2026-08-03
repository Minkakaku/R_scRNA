library(Seurat)
library(tidyverse)



all_data_sample1 <- all_data_ls_pass$A02994A6

# sample1 -----------------

## normalize ------
all_data_sample1 <- NormalizeData(all_data_sample1, scale.factor = median(all_data_sample1$nCount_RNA))
all_data_sample1 <- FindVariableFeatures(all_data_sample1)
all_data_sample1 <- ScaleData(all_data_sample1, vars.to.regress = 'nCount_RNA')

## ## sketch data --------
# all_data_sample1 <-
#   SketchData(
#     all_data_sample1,
#     ncells         = 50000L,
#     sketched.assay = "sketch",
#     method         = "LeverageScore",
#     var.name       = "leverage.score",
#     seed = 123L
#   )


## dimension reduction - pca ------------
dr_dir_sample1 <- "dimension reduction - sample1"
dir.create(dr_dir_sample1)

# DefaultAssay(all_data_sample1) <- "sketch"
# all_data_sample1 <- FindVariableFeatures(all_data_sample1, nfeatures = 3000)
# all_data_sample1 <- ScaleData(all_data_sample1, vars.to.regress = 'nCount_RNA')

# all_data_sample1 <- SCTransform(all_data_sample1, assay  = "sketch")

all_data_sample1 <- RunPCA(all_data_sample1)

pdf(file.path(dr_dir_sample1, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data_sample1, ndims = 30)
dev.off()


set.seed(1234)
subset_cells_sample1 <- 
  all_data_sample1[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data_sample1) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir_sample1, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data_sample1[, subset_cells_sample1], dims = 1:5)
dev.off()

pdf(file.path(dr_dir_sample1, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data_sample1[, subset_cells_sample1], dims = 6:10)
dev.off()

pdf(file.path(dr_dir_sample1, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data_sample1[, subset_cells_sample1], dims = 11:15)
dev.off()

pdf(file.path(dr_dir_sample1, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data_sample1[, subset_cells_sample1], dims = 16:20)
dev.off()


## dimension reduction - umap --------
dims <- 10


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


all_data_sample1 <- RunUMAP(all_data_sample1, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_dir_sample1, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data_sample1, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data_sample1 <- all_data_sample1 %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "roi", # "project",
    reduction.save = "harmony"
  )

all_data_sample1 <- RunUMAP(
  all_data_sample1,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf( file.path(dr_dir_sample1, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data_sample1, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'

## cluster ------
all_data_sample1 <-
  FindNeighbors(
    all_data_sample1,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data_sample1 <-
  FindClusters(
    all_data_sample1,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(dr_dir_sample1, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_sample1, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

##  









## # cell marker ----
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  Enterocyte  = "SLC26A3",
  Enterocyte  = "CA1",
  Enterocyte  = "GUCA2A",
  Enterocyte  = "GUCA2B",
  
  Goblet = "MUC2",
  Goblet = "TFF3",
  
  Enteroendocrine = "CHGA",
  
  # "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  # "Fibroblasts" = 'COL1A2',
  # "Fibroblasts" = 'COL3A1',
  # "Fibroblasts" = 'LUM',
  # "Fibroblasts" = 'FN1',
  # "Fibroblasts" = 'COL6A1',
  # "Fibroblasts" = 'FBLN1',
  # 
  # "Fibroblasts" = 'VIM',
  # "Fibroblasts" = 'THY',
  # "Fibroblasts" = 'S100A4',
  # "Fibroblasts" = 'PDPN',
  # "Fibroblasts" = 'FAP',
  
  
  # 
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  
  # "Smooth muscle cells" = "PDGFRB", 
  # "Smooth muscle cells" = "PLN", 
  # "Smooth muscle cells" = "MYLK", 
  # "Smooth muscle cells" = "MYH11", 
  # "Smooth muscle cells" = "ACTA2",
  
  
  # "Myefibroblast" = "PDGFRB", 
  "Myefibroblast" = "PLN", 
  "Myefibroblast" = "MYLK", 
  "Myefibroblast" = "MYH11", 
  # "Myefibroblast" = "ACTA2",
  
  # "Breast glandular cells" = "SERPINA3", 
  # "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  # Plasmablast = "MZB1",
  # Plasmablast = "XBP1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "MS4A2",
  "Mast cell" = "KIT",
  
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Neutrophi' = "S100A8",
  'Neutrophi' = "S100A9",
  'Neutrophi' = "CSF3R",
  'Neutrophi' = "FCGR3B",
  
  Gial = "MPZ"
)


fibro_markers <- c(
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "Fibroblasts" = 'COL1A2',
  "Fibroblasts" = 'COL3A1',
  "Fibroblasts" = 'LUM',
  "Fibroblasts" = 'FN1',
  "Fibroblasts" = 'COL6A1',
  "Fibroblasts" = 'FBLN1',
  
  "Fibroblasts" = 'VIM',
  "Fibroblasts" = 'THY',
  "Fibroblasts" = 'S100A4',
  "Fibroblasts" = 'PDPN',
  "Fibroblasts" = 'FAP'
)

fibro_markers %in% rownames(all_data_sample1)

# dotplot
pdf(file.path(dr_dir_sample1, str_glue("dotplot dims{dims} {reduction.use} marker.pdf")), width = 25, height = 6)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()

pdf(file.path(dr_dir_sample1, str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data_sample1, group.by = "RNA_snn_res.0.8", features = major_marker_ls)
dev.off()

## annotation ----
Idents(all_data_sample1) <- "sketch_snn_res.0.2"



cell_cluster <- c(
  '0'  = "Epithelial",
  '1'  = "Fibroblast",
  '2'  = "Epithelial",
  '3'  = "Macrophage",
  '4'  = "Plasma",
  '5'  = "T",
  '6'  = "Endothelial",
  '7'  = "Myefibroblast",
  '8'  = "Pericyte",
  '9'  = "B",
  '10' = "Gial",
  '11' = "Mast cell",
  '12' = "NK",
  '13' = "Others",
  '14' = "Others"
)


all_data_sample1 <- RenameIdents(all_data_sample1, cell_cluster)
all_data_sample1$annotation <- Idents(all_data_sample1)


pdf(file.path(dr_dir_sample1, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data_sample1, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data_sample1, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_dir_sample1, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data_sample1, group.by = "annotation", features = major_marker_ls)
dev.off()


## find markers -------
Idents(all_data_sample1) <- "sketch_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data_sample1),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data_sample1))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data_sample1), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )


## project data ------
DefaultAssay(all_data_sample1) <- "sketch"
all_data_sample1 <- 
  all_data_sample1 %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",# "pca",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(all_data_sample1) <- 'annotation_full'
DefaultAssay(all_data_sample1) <- "RNA"

pdf(file.path(dr_dir_sample1, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data_sample1, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data_sample1, "all_data_sample1_projected.rds")
# all_data_sample1 <- readRDS("all_data_sample1_projected.rds")
