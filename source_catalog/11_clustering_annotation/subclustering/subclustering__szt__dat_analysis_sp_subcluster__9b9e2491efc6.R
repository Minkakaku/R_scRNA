library(tidyverse)
library(Seurat)
library(patchwork)

source("sc_analysis_tools.R")

# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb



all_data <- readRDS("all_data_projected - project batch.rds")

out_dir <- "sp_subcluster"
dir.create(out_dir)

all_data$annotation_full %>% table
# CEACAM6 Epithelial      CEACAM6 Goblet          Enterocyte    EPCAM Epithelial         Fibroblasts         MUC2 Goblet             Myeloid 
# 95997              120849               37931               40551               72181               87837                5448 
# Others           Pericytes              Plasma Proliferating cells Smooth muscle cells         TFF3 Goblet 
# 21156                7228               33760                3550               13420               22175


# fibro subcluster ----------------
sp_fibro <- all_data[, all_data$annotation_full == "Fibroblasts"]

DefaultAssay(sp_fibro) <- 'Spatial.008um'

# sp_fibro <- sp_fibro[, sp_fibro$nCount_Spatial.008um > 50]

sp_fibro <- sc_normalize(sp_fibro, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
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


sp_fibro <- sc_dotplot(
  sp_fibro, 
  ndims = ndims_fibro2,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "fibro/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

# saveRDS(sp_fibro, file.path(out_dir, 'sp_fibro.rds'))



# fibro2 subcluster ----------------
sp_fibro2 <- sp_fibro[, ! sp_fibro$Spatial.008um_snn_res.0.2 %in% c(1, 7:11)]

DefaultAssay(sp_fibro2) <- 'Spatial.008um'

# sp_fibro2 <- sp_fibro2[, sp_fibro2$nCount_Spatial.008um > 50]

sp_fibro2 <- sc_normalize(sp_fibro2, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
sp_fibro2 <- sc_pca(
  sp_fibro2, 
  sc_out_dir = file.path(out_dir, "fibro2/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_fibro2 <- 10

sp_fibro2 <- sc_dr(
  sp_fibro2, 
  ndims = ndims_fibro2,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "fibro2/sc.2 dr"),
  dir.create.force = TRUE
)


sp_fibro2 <- sc_cluster(
  sp_fibro2, 
  ndims = ndims_fibro2,
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "fibro2/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
)


caf.signature_ls <- 
  readr::read_tsv("caf.signature.txt", comment = "#") %>%
  map(function(x){
    x[! is.na(x)]
  })

sp_fibro2_score <- sp_fibro2 %>% JoinLayers() %>% AddModuleScore(features = caf.signature_ls)

sp_fibro2_score[[c("Spatial.008um_snn_res.0.2", "Cluster1", "Cluster2", "Cluster3")]] %>%
  group_by(Spatial.008um_snn_res.0.2) %>%
  pivot_longer(Cluster1:Cluster3, names_to = "Caf", values_to = "Score") %>%
  ggplot(aes(x = Spatial.008um_snn_res.0.2, y = Score, fill = Caf)) +
  geom_violin()
  
  
major_marker_for_subcluster <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "RGS5",
  
  "Smooth muscle cells" = "TAGLN", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  
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
  'Dendritic cell' = "CLEC9A",
  
  Neutrophil = "G0S2",
  Neutrophil = "CSF3R",
  Neutrophil = "S100A8",
  Neutrophil = "S100A9"
)

fibro_markers <- c(
  Epithelial  = "EPCAM", 
  # Epithelial  = "KRT5",
  Epithelial  = "KRT8",
  Epithelial  = "S100A8",
  Epithelial  = "KRT18",
  
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
  'vCAF'='RGS5',
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
  'pericyte'='COL4A1',
  'pericyte'='COL4A2',
  'hsptcaf'='HSPH1',
  'hsptcaf'='VEGFA'
)

sp_fibro2 <- sc_dotplot(
  sp_fibro2, 
  ndims = ndims_fibro2,
  marker = fibro_markers,
  sc_out_dir = file.path(out_dir, "fibro2/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_fibro2 <- sc_dotplot(
  sp_fibro2, 
  ndims = ndims_fibro2,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "fibro2/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_fibro <- 
  sc_findmarker(
    sp_fibro2,
    idents =  "Spatial.008um_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.15,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )

sc_findmarker(sp_fibro2, all_markers = all_markers_fibro, sort_by = "avg")

dir.create(file.path(out_dir, "fibro2/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "fibro2/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_fibro2,
  markers = all_markers_fibro,
  idents = "Spatial.008um_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_fibro2[[c('Spatial.008um_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

sc_mark_heatmap(
  sp_fibro2,
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
  '0'  = "apCAF",
  '1'  = "Epi-like Fibroblast",
  '2'  = "mCAF",
  '3'  = "infCAF",
  '4'  = "vCAF"
)



Idents(sp_fibro2) <- "Spatial.008um_snn_res.0.2"
sp_fibro2 <- RenameIdents(sp_fibro2, fibro_cell_cluster)
sp_fibro2$annotation_fibro <- Idents(sp_fibro2)

dir.create(file.path(out_dir, "fibro2/sc.5 annotation"))
pdf(file.path(out_dir, 'fibro2/sc.5 annotation', str_glue("Dimplot dim{ndims_fibro2} fibro annotation.pdf")), height = 6, width = 8)
DimPlot(sp_fibro2, group.by = "annotation_fibro", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_fibro2, group.by = "annotation_fibro", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

sp_fibro2_2 <- sp_fibro2[, sp_fibro2$annotation_fibro != "Epi-like Fibroblast"]

dir.create(file.path(out_dir, "fibro2/sc.5 annotation"))
pdf(file.path(out_dir, 'fibro2/sc.5 annotation', str_glue("Dimplot dim{ndims_fibro2} fibro annotation - remove epi-like cellls.pdf")), height = 6, width = 8)
DimPlot(sp_fibro2_2, group.by = "annotation_fibro", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_fibro2_2, group.by = "annotation_fibro", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

sp_fibro2_2 <- sc_normalize(sp_fibro2_2, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
sp_fibro2_2 <- sc_pca(
  sp_fibro2_2, 
  sc_out_dir = file.path(out_dir, "fibro2_2/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_fibro2_2 <- 10

sp_fibro2_2 <- sc_dr(
  sp_fibro2_2, 
  ndims = ndims_fibro2_2,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "fibro2_2/sc.2 dr"),
  dir.create.force = TRUE
)

saveRDS(sp_fibro2, file.path(out_dir, 'sp_fibro2.rds'))
saveRDS(sp_fibro2_2, file.path(out_dir, 'sp_fibro2_2.rds'))




# B subcluster ----------------
sp_b <- all_data[, all_data$annotation_full == "Plasma"]

DefaultAssay(sp_b) <- 'Spatial.008um'

sp_b <- sc_normalize(sp_b, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
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
  resolution = 0.1 # c(0.2, 0.4, 0.6, 0.8)
)



b_markers2 <- c(
  'CD19',
  'MS4A1',
  'CD79A',
  'JCHAIN',
  'IGHA1',
  'IGHG1',
  
  B_naive  = "IGHD",
  B_naive  = "TCL1A",
  B_naive  = "IGHM",
  
  B_memory = "CD27",
  B_memory = "TNFRSF13B",
  
  B_gc    = "BCL6",
  B_gc    = "AICDA",
  B_gc    = "RGS13",
  
  ASC     = "MZB1",
  ASC     = "XBP1",
  
  cycling = "MKI67"
)

sp_b <- sc_dotplot(
  sp_b, 
  ndims = ndims_b,
  marker = b_markers2,
  resolution = 0.1,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_b <- sc_dotplot(
  sp_b, 
  ndims = ndims_b,
  marker = major_marker_for_subcluster,
  resolution = 0.1,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_b <- 
  sc_findmarker(
    sp_b,
    idents =  "Spatial.008um_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.10,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


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
  idents = "RNA_snn_res.0.2", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  sort_by = "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_b[[c('RNA_snn_res.0.2', 'project', 'group')]]  %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = RNA_snn_res.0.2)) +
  geom_col()

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



Idents(sp_b) <- "RNA_snn_res.0.6"
sp_b <- RenameIdents(sp_b, b_cell_cluster)
sp_b$annotation_b <- Idents(sp_b)

dir.create(file.path(out_dir, "b/sc.5 annotation"))
pdf(file.path(out_dir, 'b/sc.5 annotation', str_glue("Dimplot dim{ndims_b} b annotation.pdf")), height = 6, width = 8)
DimPlot(sp_b, group.by = "annotation_b", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_b, group.by = "annotation_b", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_b, file.path(out_dir, 'sp_b.rds'))

# myeloid subcluster --------------
sp_mye <- all_data[, all_data$annotation_full == "Myeloid"]

DefaultAssay(sp_mye) <- 'Spatial.008um'

sp_mye <- sc_normalize(sp_mye, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
sp_mye <- sc_pca(
  sp_mye, 
  sc_out_dir = file.path(out_dir, "mye/sc.1 pca"),
  ndims = 50,
  dir.create.force = TRUE
)

ndims_mye <- 10

sp_mye <- sc_dr(
  sp_mye, 
  ndims = ndims_mye,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "mye/sc.2 dr"),
  dir.create.force = TRUE
)


sp_mye <- sc_cluster(
  sp_mye, 
  ndims = ndims_mye,
  reduction =  "harmony", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "mye/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

mye_markers <- c(
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
  
  #Neutrophil
  Neutrophil = "G0S2",
  Neutrophil = "CSF3R",
  Neutrophil = "S100A8",
  Neutrophil = "S100A9"
  
)

sp_mye <- sc_dotplot(
  sp_mye, 
  ndims = ndims_mye,
  marker = mye_markers,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_mye <- sc_dotplot(
  sp_mye, 
  ndims = ndims_mye,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)
sp_mye <- sc_dotplot(
  sp_mye, 
  ndims = ndims_mye,
  marker = c(
    "SPP1", "C1QA", "C1QC", "SELENOP", 
    "CSF1R", "IL1B", "NLRP3", "PLTP", "MARCO"
  ),
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot - papers"), 
  dir.create.force = TRUE
)


all_markers_mye <- 
  sc_findmarker(
    sp_mye,
    idents =  "Spatial.008um_snn_res.0.2", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )


sc_findmarker(sp_mye, all_markers = all_markers_mye, sort_by = "pct")
all_markers_mye %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(pct.1), .by_group = TRUE) %>%
  slice_head(n = 5)

openxlsx::write.xlsx(
  all_markers_mye,
  file = "sp_subcluster/mye/markers.xlsx",
  asTable = TRUE
)


mac_cell_cluster <- c(
  '0'  = "SELENOP Macrophage",
  '1'  = "C1QC Macrophage"
)



Idents(sp_mye) <- "Spatial.008um_snn_res.0.2"
sp_mye <- RenameIdents(sp_mye, mac_cell_cluster)
sp_mye$annotation_mac <- Idents(sp_mye)

dir.create(file.path(out_dir, "mye/sc.5 annotation"))
pdf(file.path(out_dir, 'mye/sc.5 annotation', str_glue("Dimplot dim{ndims_mye} mac annotation.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "annotation_mac", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mye, group.by = "annotation_mac", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

## marker gene enrichment  -------
df_all_marker_mye <-
  FindMarkers(
    JoinLayers(sp_mye),
    ident.1 = "C1QC Macrophage",
    ident.2 = "SELENOP Macrophage",
    group.by = "annotation_mac",
    logfc.threshold = 0,
    min.pct = 0
  )

diff_gene_ls <- 
  df_all_marker_mye %>%
  rownames_to_column("gene") %>%
  mutate(Group = case_when(
    avg_log2FC >  0.25 ~ "C1QC Macrophage",
    avg_log2FC < -0.25 ~ "SELENOP Macrophage",
    TRUE ~ "Others"
  )) %>%
  dplyr::filter(Group != "Others") %>%
  group_by(Group) %>%
  nest() %>%
  deframe() %>%
  map(function(df){
    df %>%
      dplyr::filter(
        # pct.1 > 0.05,
        p_val_adj < 0.05
      )
  }) %>%
  map("gene")

diff_gene_id_ls <-
  diff_gene_ls %>%
  map(function(gene){
    clusterProfiler::bitr(
      gene,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$ENTREZID
  })

all_enrich_mye_ls <- 
  diff_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# add gene
all_enrich_mye_ls <-
  all_enrich_mye_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    xx2 <-
      xx %>%
      map_chr(function(g){
        clusterProfiler::bitr(
          g,
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$SYMBOL %>%
          str_c(collapse = "/")
      })
    x@result$geneSymbol <- xx2
    x
  })

# enrichment plot
p_all_enrich_mye_ls <-
  all_enrich_mye_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })

# export
all_enrich_mye_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path("sp_subcluster/mye", "gene enrichment analysis table.xlsx"),
    asTable = TRUE
  )

pdf(
  file.path("sp_subcluster/mye", "enrichment analysis dotplot.pdf"), 
  height = 6, 
  width = 8
)
p_all_enrich_mye_ls %>%  walk(print)
dev.off()





## macrophage pathway score -------
mac_pathway_ls <- list(
  M1 = c("IL12A", "IL23A", "TNF", "IL6", "CD86", "HLA-DMA", "HLA-DMB", 
        "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", 
        "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "IL1B", 
        "MARCO", "NOS2", "FCGR1A", "CD80", "CXCL9", "CXCL10", "CXCL11", 
        "CD86", "IL1A", "IL1B", "IL6", "TNF", "CCL5", "IRF5", "IRF1", 
        "CD40", "IDO1", "KYNU", "CCR7", "PTPRC", "CD68", "CSF1R", "LY75", 
        "CD14"),
  M2 = c("ARG1", "ARG2", "IL10", "FCGR2A", "CD163", "FCER2", "CD200R1", 
         "PDCD1LG2", "CD274", "MARCO", "CSF1R", "MRC1", "IL1RN", "IL1R2", 
         "IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", 
         "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", 
         "CTSA", "CTSB", "CST3", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", 
         "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASLG", "TNFSF12", "TNFSF8", 
         "CD276", "VTCN1", "MSR1", "FN1", "IRF4", "PTPRC", "CD68", "CSF1R", 
         "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "LY75", "CD14"),
  Phagocytosis = c("MRC1", "CD163", "MERTK", "C1QB"),
  Angiogenesis = c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", "EZH2", 
                   "FGF18", "FGFR1", "FYN", "HEY1", "ITGAV", "JAG1", "JAG2", "MMP9", 
                   "NOTCH1", "PDGFA", "PTK2", "SPP1", "STC1", "TNFAIP6", "TYMP", 
                   "VAV2", "VCAN", "VEGFA")
  
)

sp_mye_score <- 
  sp_mye %>%
  JoinLayers() %>%
  AddModuleScore(features = mac_pathway_ls) 
colnames(sp_mye_score@meta.data)[match(paste0("Cluster", 1:4), colnames(sp_mye_score@meta.data))] <- names(mac_pathway_ls)


pdf(
  file.path("sp_subcluster/mye", "macrophage pathway score feature plot.pdf"), 
  height = 4, 
  width = 14
)
  dim_path_ls <- 
    FeaturePlot(sp_mye_score, features = c("M1", "M2"), order = T, combine = FALSE) %>%
    map(function(p) p +  scale_color_gradientn(colours = Seurat:::SpatialColors(50), transform = "log1p"))
  dim_path_ls[['CellType']] <- DimPlot(sp_mye_score, group.by = "annotation_mac") & ggsci::scale_color_aaas()
  patchwork::wrap_plots(dim_path_ls, nrow = 1)
dev.off()


pdf(
  file.path("sp_subcluster/mye", "macrophage pathway score Vlnplot.pdf"), 
  height = 5, 
  width = 10
)
VlnPlot(sp_mye_score, group.by = "annotation_mac", features = c('M1', 'M2'), pt.size = 0) &
  geom_boxplot(width = 0.15, fill = "white", outliers = F) &
  ggsci::scale_fill_aaas() &
  # scale_x_discrete(labels = function(x) str_remove(x, "\\s.*$")) &
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) &
  labs(
    x = NULL
  )

VlnPlot(sp_mye_score, group.by = "annotation_mac", features = c('Phagocytosis', 'Angiogenesis'), pt.size = 0) &
  geom_boxplot(width = 0.15, fill = "white", outliers = F) &
  ggsci::scale_fill_aaas() &
  # scale_x_discrete(labels = function(x) str_remove(x, "\\s.*$")) &
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) &
  labs(
    x = NULL
  )

dev.off()


pdf(
  file.path("sp_subcluster/mye", "macrophage pathway score - group.pdf"), 
  height = 5, 
  width = 10
)

VlnPlot(
  sp_mye_score, 
  group.by = "annotation_mac",
  features = c('Phagocytosis', 'Angiogenesis'),
  split.by = "group",
  pt.size = 0
)  

p_selenop_dot <- 
  DotPlot(
    sp_mye_score[, sp_mye_score$annotation_mac == "SELENOP Macrophage"],
    features = c('Phagocytosis', 'Angiogenesis'),
    scale = F,
    group.by = "group"
  ) +
  ggtitle("SELENOP Macrophage")
p_c1qc_dot <-
  DotPlot(
    sp_mye_score[, sp_mye_score$annotation_mac == "C1QC Macrophage"],
    features = c('Phagocytosis', 'Angiogenesis'),
    scale = F,
    group.by = "group"
  )  +
  ggtitle("C1QC Macrophage")
p_selenop_dot + p_c1qc_dot
dev.off()  


sp_mye[[]] %>%
  # rownames_to_column("cell_id") %>%
  # left_join(all_data_meta_cn[c("cell_id", "CN")], by = "cell_id") %>% 
  # 
  
  group_by(group, annotation_mac) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_mac)) +
  geom_col() +
  ggthemes::theme_clean() +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )

DotPlot(
  sp_mye_score,
  features = mac_pathway_ls$Phagocytosis,
  scale = F,
  split.by = "annotation_mac",
  group.by = "group"
)



sp_mye_score[[]] %>% 
  dplyr::select(all_of(c(
    # names(mac_pathway_ls),
    c('Phagocytosis', 'Angiogenesis'),
    'group',
    'annotation_mac'
  ))) %>%
  pivot_longer(
    # all_of(names(mac_pathway_ls)), 
    c('Phagocytosis', 'Angiogenesis'), 
    names_to = "pathway",
    values_to = "score"
  ) %>%
  ggplot(aes(x = annotation_mac, y = score)) +
  geom_violin(aes(fill = group), width = 1) +
  facet_wrap(~ pathway, scales = "free") +
  ggthemes::theme_clean()






## macrophage de analysis -----------
FindMarkers(
  JoinLayers(sp_mye),
  ident.1 = "Tumor",
  group.by = "group",
  only.pos = TRUE,
  subset.ident = "C1QC Macrophage"
) %>%
  DT::datatable(filter = "bottom")



all_data



# saveRDS(sp_mye, file.path(out_dir, 'sp_mye.rds'))


# myeloid subcluster 2 --------------
sp_mye2 <- sp_mye[, sp_mye$Spatial.008um_snn_res.0.4 %in% c(2, 3, 4)]

DefaultAssay(sp_mye2) <- 'Spatial.008um'

sp_mye2 <- sc_normalize(sp_mye2, scale.factor = "def", vars.to.regress = "nCount_Spatial.008um")
sp_mye2 <- sc_pca(
  sp_mye2, 
  sc_out_dir = file.path(out_dir, "mye2/sc.1 pca"),
  ndims = 50,
  dir.create.force = TRUE
)

ndims_mye2 <- 10

sp_mye2 <- sc_dr(
  sp_mye2, 
  ndims = ndims_mye2,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "mye2/sc.2 dr"),
  dir.create.force = TRUE
)


sp_mye2 <- sc_cluster(
  sp_mye2, 
  ndims = ndims_mye2,
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "mye2/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

mye_markers <- c(
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
  
  #Neutrophil
  Neutrophil = "G0S2",
  Neutrophil = "CSF3R",
  Neutrophil = "S100A8",
  Neutrophil = "S100A9"
  
)

sp_mye2 <- sc_dotplot(
  sp_mye2, 
  ndims = ndims_mye2,
  marker = mye_markers,
  sc_out_dir = file.path(out_dir, "mye2/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sp_mye2 <- sc_dotplot(
  sp_mye2, 
  ndims = ndims_mye2,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "mye2/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_mye <- 
  sc_findmarker(
    sp_mye2,
    idents =  "Spatial.008um_snn_res.0.4", # "sketch_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "project"
  )



sc_findmarker(sp_mye2, all_markers = all_markers_mye, sort_by = "pct")

pdf(
  file.path(out_dir, "mye2", 'annotation.pdf'),
  height = 6,
  width  = 9
)
FeaturePlot(
  sp_mye2,
  features = c("CD14", "CD68", "C1QA"),
  order = TRUE
)

DimPlot(sp_mye2, group.by = "Spatial.008um_snn_res.0.4")

all_markers_mye %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(pct.1), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  dplyr::pull("gene") %>%
  # c("CD14", "CD68", "C1QC", "C1QA") %>%
  
  DotPlot(
    sp_mye2,
    features = .,
    group.by = "Spatial.008um_snn_res.0.4"
  )
dev.off()



dir.create(file.path(out_dir, "mye2/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "mye2/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sp_mye2,
  markers = all_markers_mye,
  assay = "Spatial.008um",
  idents = "Spatial.008um_snn_res.0.4", #"sketch_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.1,
  sort_by = "pct", # "avg",
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sp_mye2[[c('Spatial.008um_snn_res.0.4', 'project', 'group')]]  %>%
  group_by(project, group, Spatial.008um_snn_res.0.4) %>%
  summarise(count = n()) %>%
  # group_by(project) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = Spatial.008um_snn_res.0.4)) +
  geom_col()

sc_mark_heatmap(
  sp_mye2,
  markers = all_markers_mye,
  idents = "Spatial.008um_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

mac_cell_cluster2 <- c(
  '0'  = "C1QC Macrophage",
  '1'  = "SELENOP Monocyte"
)



Idents(sp_mye2) <- "Spatial.008um_snn_res.0.4"
sp_mye2 <- RenameIdents(sp_mye2, mac_cell_cluster2)
sp_mye2$annotation_mac <- Idents(sp_mye2)

dir.create(file.path(out_dir, "mye2/sc.5 annotation"))
pdf(file.path(out_dir, 'mye2/sc.5 annotation', str_glue("Dimplot dim{ndims_mye2} mac annotation.pdf")), height = 6, width = 8)
DimPlot(sp_mye2, group.by = "annotation_mac", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mye2, group.by = "annotation_mac", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sp_mye2, file.path(out_dir, 'sp_mye2.rds'))

## macrophage pathway score -------
sp_mye_score2 <- 
  sp_mye2 %>%
  JoinLayers() %>%
  AddModuleScore(features = mac_pathway_ls) 
colnames(sp_mye_score2@meta.data)[match(paste0("Cluster", 1:4), colnames(sp_mye_score2@meta.data))] <- names(mac_pathway_ls)


pdf(
  file.path("sp_subcluster/mye2", "macrophage pathway score feature plot.pdf"), 
  height = 4, 
  width = 14
)
dim_path_ls <- 
  FeaturePlot(sp_mye_score2, features = c("M1", "M2"), order = T, combine = FALSE) %>%
  map(function(p) p +  scale_color_gradientn(colours = Seurat:::SpatialColors(50), transform = "log1p"))
dim_path_ls[['CellType']] <- DimPlot(sp_mye_score2, group.by = "annotation_mac") & ggsci::scale_color_aaas()
patchwork::wrap_plots(dim_path_ls, nrow = 1)
dev.off()


pdf(
  file.path("sp_subcluster/mye2", "macrophage pathway score Vlnplot.pdf"), 
  height = 5, 
  width = 10
)
VlnPlot(sp_mye_score2, group.by = "annotation_mac", features = c('M1', 'M2'), pt.size = 0) &
  geom_boxplot(width = 0.15, fill = "white", outliers = F) &
  ggsci::scale_fill_aaas() &
  # scale_x_discrete(labels = function(x) str_remove(x, "\\s.*$")) &
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) &
  labs(
    x = NULL
  )

VlnPlot(sp_mye_score2, group.by = "annotation_mac", features = c('Phagocytosis', 'Angiogenesis'), pt.size = 0) &
  geom_boxplot(width = 0.15, fill = "white", outliers = F) &
  ggsci::scale_fill_aaas() &
  # scale_x_discrete(labels = function(x) str_remove(x, "\\s.*$")) &
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) &
  labs(
    x = NULL
  )
dev.off()






# combined annotation ---------

sub_cluster_df <-
  rbind(
    sp_fibro2[['annotation_fibro']] %>% dplyr::rename(annotation = 1),
    sp_mye[['annotation_mac']]      %>% dplyr::rename(annotation = 1)
  ) %>% 
  rownames_to_column("cell_id")

sub_cluster_df_tidy <-
  all_data[['annotation_full']] %>% 
  rownames_to_column("cell_id") %>% 
  left_join(sub_cluster_df, by = 'cell_id') %>%
  map_dfc(as.character) %>%
  mutate(annotation3 = ifelse(is.na(annotation), annotation_full,  annotation)) %>%
  dplyr::select(cell_id, annotation3) %>% 
  column_to_rownames("cell_id")

all_data <-  
  all_data %>% 
  AddMetaData(sub_cluster_df_tidy) 


# saveRDS(all_data, 'all_data_subcluster.rds')

all_data %>% 
  DimPlot(group.by = "annotation3", reduction = "harmonyumap") +
  ggsci::scale_color_igv()



all_data[[c("group", "annotation3")]] %>%
  group_by(group, annotation3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation3)) +
  geom_col() +
  ggsci::scale_fill_igv()
  












