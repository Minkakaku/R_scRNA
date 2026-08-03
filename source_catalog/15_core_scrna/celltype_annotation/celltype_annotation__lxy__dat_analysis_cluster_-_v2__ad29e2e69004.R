library(Seurat)
library(tidyverse)

options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
future::plan("multisession", workers = 4, gc = TRUE)
# future::plan("sequential")
# future:::ClusterRegistry("stop")


out_dir_cluster <- "dat_analysis_cluster"
dir.create(out_dir_cluster)

# # running base report pipeline and get cd_final
# cd_final %>% dim
# cd_final %>% LayerData() %>% dim
#
# saveRDS(cd_final, file.path(out_dir_cluster, "cd_final.rds"))
#
# # anno data    ------
# anno_data <- readxl::read_excel("data/Annotation.xlsx") %>% dplyr::slice(-c(1:2))
#
#
# # as seurat    ------
# sc <- cd_final %>% celldive2Seurat()
#
# meta_data_for_add_seurat <-
#   sc[['ROI_ID']] %>%
#   rownames_to_column("ID") %>%
#   left_join(anno_data[c("ROI_ID", "Group")], by = "ROI_ID") %>%
#   column_to_rownames("ID")
#
# sc <- sc %>% AddMetaData(meta_data_for_add_seurat)


sc <- readRDS("pipeline/exportDat/sc.rds")
# sc <- readRDS("sc_raw.rds")


# QC    ------
sc$log10GenesPerUMI <- log10(sc$nFeature_RNA + 1) / log10(sc$nCount_RNA + 1)

pdf(
  file.path(out_dir_cluster, "QC.pdf"),
  height = 6,
  width  = 8
)
VlnPlot(sc, features = "nCount_RNA", pt.size = 0, group.by = 'Slide_ID')
VlnPlot(sc, features = "nCount_RNA", pt.size = 0, log = TRUE, group.by = 'Slide_ID')
VlnPlot(sc, features = "log10GenesPerUMI", pt.size = 0, group.by = 'Slide_ID')
dev.off()


ncount_quantile_rng <- quantile(sc$nCount_RNA, c(0.001, 0.999))
log10GeneperUMI_quantile_rng <- quantile(sc$log10GenesPerUMI, c(0.001, 0.999))


sc <- sc[ ,
          sc$nCount_RNA > ncount_quantile_rng[1] &
            sc$nCount_RNA < ncount_quantile_rng[2] &
            sc$log10GenesPerUMI > log10GeneperUMI_quantile_rng[1] &
            sc$log10GenesPerUMI < log10GeneperUMI_quantile_rng[2]
]


# hvg    ------
sc <- sc %>% FindVariableFeatures()

p_hvg <- VariableFeaturePlot(sc)
p_hvg$data$gene <- rownames(p_hvg$data)
p_hvg2 <- p_hvg + ggrepel::geom_text_repel(aes(label = gene))
p_hvg2

VariableFeatures(sc)

sc@assays$RNA@layers$data

sc %>% head

pdf(
  file.path(out_dir_cluster, "hvg.pdf"),
  height = 6,
  width  = 8
)
p_hvg2 %>% print()
dev.off()

# norm   ------
normalize_q3 <- function(dat, log_transform = TRUE, min_val = 1e-5){
  stopifnot(length(dim(dat)) == 2)
  # stopifnot(is.numeric(dat))

  # dat <- dat + min_val

  q3 <- apply(dat, 2, quantile, 0.75)
  mean_q3  <- mean(q3)

  scale_fators <- mean_q3/q3
  dat2 <- sweep(dat, 2, scale_fators, '*')

  if(log_transform) dat2 <- log1p(dat2)
  dat2
}


norm_assay <-
  SetAssayData(
    sc@assays$RNA,
    slot = 'data',
    new.data = normalize_q3(LayerData(sc, "counts"))
  )

pdf(
  file.path(out_dir_cluster, "normalize.pdf"),
  height = 6,
  width  = 8
)

boxplot(as.matrix(norm_assay$counts[,1:50]))
boxplot(as.matrix(norm_assay$data[,1:50]))
boxplot(log1p(as.matrix(norm_assay$counts[,1:50])))
dev.off()

# norm_assay$data[, "CWMDA_PCF_Slide_1...TS1_X2_Y1...TotalSegment...6233"]

sc[['RNA']] <- norm_assay



# # 20251107 only using typing marker
# typing_feature <- c(
#   "CD31",
#   "CD3e", "CD4", "CD8", "FoxP3",
#   "Ly6G", "S100A9",
#   "CD20", "CD45R/B220",
#   "F4/80", "CD206", "CD68", "CD11c", "Iba1",
#   "COL1A1", "Vimentin",
#   "PanCK",
#   "CD45", "Ki67", "Ter119"
# )
# sc[['type']] <- sc[typing_feature, ][['RNA']]
#
# DefaultAssay(sc) <- "type"


# # scaled by DAPI
# sc$DAPI <- FetchData(sc, vars = "DAPI")
# VariableFeatures(sc) <- VariableFeatures(sc) %>% setdiff("DAPI")
#
# sc <- sc %>% ScaleData( vars.to.regress = "DAPI" )

sc <- sc %>% ScaleData( )



# pca    ------
sc <- sc %>% RunPCA()
ElbowPlot(sc, ndims = 30)



pdf(
  file.path(out_dir_cluster, "pca.pdf"),
  height = 6,
  width  = 8
)
ElbowPlot(sc, ndims = 30)
dev.off()



ndims <- 7

sc <- RunUMAP(sc, dims = 1:ndims)

# harmony    ------
sc <- sc %>%
  harmony::RunHarmony(
    dims.use = 1:ndims,
    group.by.vars = 'ROI_ID', # 'ROI_ID',
    reduction.save = "harmony",
    assay.use = "RNA"
  )

sc <- sc %>%
  RunUMAP(
    dims = 1:ndims,
    # min.dist = 0.1,
    # spread   = 1.5,
    # n.neighbors = 50,
    # local.connectivity = 3,
    reduction = "harmony",
    reduction.name = "harmonyumap",
    reduction.key  = "harmonyUMAP_"
  )


pdf(
  file.path(out_dir_cluster, str_glue("umap - sample - dim{ndims}.pdf")),
  height = 6,
  width  = 10
)
DimPlot(sc, reduction = "umap", group.by  = c("ROI_ID", "Group"))
DimPlot(sc,reduction = "harmonyumap", group.by  = c("ROI_ID", "Group") )
dev.off()


# group
DimPlot(
  sc,
  reduction = "umap",
  group.by  = c("ROI_ID", "Group")
)

# group
DimPlot(
  sc,
  reduction = "harmonyumap",
  group.by  = c("ROI_ID", "Group")
)



FeaturePlot(
  sc,
  features = "CD3e",
  reduction = "harmonyumap"
) +
  scale_color_gradientn(
    colours = c("white", "white", RColorBrewer::brewer.pal(8, "Blues"))
  )

# cluster    ------
sc <-
  FindNeighbors(
    sc,
    reduction = 'pca', # 'harmony',
    dims = 1:ndims,
    k.param = 15
  )

sc <-
  FindClusters(
    sc,
    resolution = c(0.2, 0.3, 0.4, 0.5, 0.6)
  )


pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims}.pdf")),
  height = 11,
  width  = 20
)
  DimPlot(
    sc,
    group.by = str_glue('RNA_snn_res.0.{c(2:6)}'),
    reduction = 'umap', # "harmonyumap",
    label = TRUE,
    raster = TRUE,
    # raster.dpi = 100
  ) &
  NoLegend()
dev.off()



# sc$RNA_snn_res.0.3 %>% table()  %>%  .[str_sort(names(.), numeric = TRUE)]

# cell_gt_1000 <- sc$RNA_snn_res.0.3 %>% table() %>% .[. > 1000] %>% names
# sc <- sc[, sc$RNA_snn_res.0.3 %in% cell_gt_1000]
# sc$RNA_snn_res.0.3 <- sc$RNA_snn_res.0.3 %>% factor(., levels=str_sort(unique(.), numeric = T))

# find markers -------
select_res <- "RNA_snn_res.0.6"


Idents(sc) <- select_res
DefaultAssay(sc) <- "RNA"
all_markers <-
  FindAllMarkers(
    sc,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# # export marker gene
# all_markers %>%
#   dplyr::select(gene, everything()) %>%
#   write.csv(
#     file = file.path(mk_output, "marker.csv"),
#     row.names = FALSE
#   )

all_de_genes <-
  all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")


all_top5_gene <-
  all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  # dplyr::filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))


# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc,
    features = all_top5_gene %>% {structure(unique(.), names = names(.)[!duplicated(.)])}
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue', 'blue'))

ggsave(
  p_dot_top5,
  filename = file.path(out_dir_cluster, str_glue("top 5 gene dotplot - {select_res} - dim{ndims}.pdf")),
  height = 8,
  width  = 18
)


# DimPlot(
#   sc[, !sc$RNA_snn_res.0.3 %in% 12:14],
#   group.by = "RNA_snn_res.0.3",
#   reduction = "harmonyumap",
#   label = TRUE
# ) +
#   NoLegend()

# 分组内皮细胞(cd31)
# T cell (CD4, CD44, CD8, cd3e, foxp3)
# 粒细胞(Ly6g, s100a9)
# B 细胞 （B220, cd20）
# 巨噬细胞 （F4/80， cd206, Iba1, cd68, cd11c）
# 成纤维细胞 （col1a1, vimentin）
# 上皮细胞（pan-Cytokeratin）

pdf(
  file.path(out_dir_cluster, str_glue("dotplot - {select_res} - dim{ndims}.pdf")),
  height = 8,
  width  = 18
)
DotPlot(
  sc, # [, !sc$RNA_snn_res.0.3 %in% 12:14],
  group.by = select_res,
  features = c(
    "CD31",
    "CD3e", "CD4", "CD8", "FoxP3", "CD44",
    "Ly6G", "S100A9",
    "CD20", "CD45R/B220",
    "F4/80", "CD206", "CD68", "CD11c", "Iba1",
    "COL1A1", "Vimentin",
    "PanCK",
    "Caveolin", "CD36", "CD45", "FCRY", "H2AX","Ki67", "Ter119"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue')) +
  ggtitle(select_res)


DotPlot(
  sc, # [, !sc$RNA_snn_res.0.3 %in% 12:14],
  group.by = select_res,
  features = c(
    "CD31",
    "CD3e", "CD4", "CD8", "FoxP3", "CD44",
    "Ly6G", "S100A9",
    "CD20", "CD45R/B220",
    "F4/80", "CD206", "CD68", "CD11c", "Iba1",
    "COL1A1", "Vimentin",
    "PanCK",
    "Caveolin", "CD36", "CD45", "FCRY", "H2AX","Ki67", "Ter119"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  ggtitle(select_res)
dev.off()



#
# DoHeatmap(
#   sc,
#   features = c(
#     "CD14" , "CD20"       , "CD31"    , "CD34" ,  "CD3e"        , "CD4"             ,
#     "CD45" , "CD45RO"     , "CD68"    , "CD8"  ,  "Collagen IV" , "E-cadherin"  ,
#     "FOXP3", "Ki67"       , "Pan-Cytokeratin" ,
#     "SMA"        , "Vimentin"
#   ),
#   group.by = "RNA_snn_res.0.3"
# )

# annotation --------
# dim7 res0.6 all markers
anno_cluster <- c(
  '0'  = "T",
  '1'  = "Macrophage",
  '2'  = "Macrophage",
  '3'  = "T",
  '4'  = "Fibroblast",
  '5'  = "Epithelial",
  '6'  = "Endothelial",
  '7'  = "Fibroblast",
  '8'  = "T",
  '9'  = "B",
  '10' = "Others",
  '11' = "Neutrophil",
  '12' = "CD36 cell",
  '13' = "Fibroblast",
  '14' = "Neutrophil",
  '15' = "Erythroid cell"
)


#  res0.6 all markers
# anno_cluster <- c(
#   '0'  = "T",
#   '1'  = "Macrophage",
#   '2'  = "Fibroblast",
#   '3'  = "T",
#   '4'  = "T-like Macrophage",
#   '5'  = "Macrophage",
#   '6'  = "Epithelial",
#   '7'  = "Tregs",
#   '8'  = "CD8T",
#   '9'  = "Endothelial",
#   '10' = "Endothelial",
#   '11' = "Neutrophil",
#   '12' = "B",
#   '13' = "Epithelial",
#   '14' = "CD8T",
#   '15' = "Fibroblast",
#   '16' = "CD36 cells",
#   '17' = "Neutrophil",
#   '18' = "Erythroid cell",
#   '19' = "Proliferating T"
# )

# # sc_raw: dim15 res0.6 all markers
# anno_cluster <- c(
#   '0'  = "CD4T",
#   '1'  = "CD8T",
#   '2'  = "Fibroblast",
#   '3'  = "Endothelial",
#   '4'  = "Macrophage",
#   '5'  = "Macrophage",
#   '6'  = "T-like Macrophage",
#   '7'  = "Proliferating Macrophage",
#   '8'  = "Epithelial",
#   '9'  = "Tregs",
#   '10' = "Epithelial",
#   '11' = "Endothelial",
#   '12' = "Fibroblast",
#   '13' = "Neutrophil",
#   '14' = "B",
#   '15' = "CD36 cell",
#   '16' = "Epithelial",
#   '17' = "Epi-like Neutrophil",
#   '18' = "Fibro-like Epithelial",
#   '19' = "Erythroid cell",
#   '20' = "Macrophage",
#   '21' = "Proliferating T"
# )

Idents(sc) <- "RNA_snn_res.0.6"
sc <- sc %>% RenameIdents(anno_cluster)
sc$annotation <- Idents(sc)


# export umap of annotation
pdf(
  file.path(out_dir_cluster,  str_glue("umap - annotation - dim{ndims} - v2.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "annotation",
  reduction = "umap", # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()

# export umap of group
pdf(
  file.path(out_dir_cluster,  str_glue("umap - group - dim{ndims} - v2.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "Group",
  reduction = "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()


# dotplot
pdf(
  file.path(out_dir_cluster, "dotplot - annotation.pdf"),
  height = 6,
  width  = 18
)
DotPlot(
  sc,
  group.by = 'annotation',
  features = c(
    "CD14" , "CD20"       , "CD31"    , "CD34" ,  "CD3e"        , "CD4"         ,
    "CD44" , "CD45RO"     , "CD68"    , "CD8"  ,  "Collagen IV" , "E-cadherin"  , "b-Catenin1",
    "FOXP3", "Ki67"       , "Pan-Cytokeratin"  , "ICOS"         ,   "IDO1"      ,
    "SMA",
    "HLA-A", "HLA-DR", "PD-1", "PD-L1", "TIGIT", "Granzyme B"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue'))

dev.off()

# saveRDS(sc, file.path(out_dir_cluster, "sc.rds"))



# fator
cell_order <- c(
  'CD4T', 'CD8T', 'Tregs', 'Epithelials', 'Proliferating epithelials',
  'CollagenIV cells',
  'Myofibroblasts', 'Fibroblasts', 'Endothelials', 'B cells',
  'Monocytes', 'Macrophage'
)

protein_order <- c(
  'Immune cells' = "CD45" ,
  'Immune cells' = "CD45RO",
  'T' = "CD3e",
  'T' = "CD4",
  'T' = "CD8",
  Tregs = "FOXP3",

  Epithelials = "Pan-Cytokeratin",
  Epithelials = "E-cadherin",
  Proliferating = "Ki67",

  "Collagen IV",
  "CD34",

  Fibroblasts = "Vimentin",
  Myofibroblasts = "SMA",

  Endothelials = "CD31",
  B = "CD20",

  Monocytes = "CD14",
  Macrophage = "CD68"
)

sc$annotation  <- sc$annotation %>% factor(levels = cell_order)


## plot umap annotation ------
p_anno <-
  DimPlot(
    sc,
    label = TRUE,
    reduction = "harmonyumap",
    group.by = "annotation"
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

ggsave(
  p_anno,
  filename = file.path(out_dir_cluster, "umap - annotation.pdf"),
  width  = 8,
  height = 6
)

## plot dotplot 1 ------
p_dot_annotation <-
  DotPlot(
    sc,
    group.by = "annotation",
    features = protein_order,
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.clip = 'off'
  ) +
  # scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  coord_cartesian(clip = "off")

p_dot_annotation2 <-
  p_dot_annotation +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues")))



pdf(
  file.path(out_dir_cluster, "dot - annotation.pdf"),
  height = 6,
  width  = 20
)
p_dot_annotation2 %>% print()
p_dot_annotation %>% print()
dev.off()


## plot dotplot 2 ------
p_dot_annotation3 <-
  DotPlot(
    sc,
    group.by = "annotation",
    features = protein_order %>% c(., setdiff(rownames(sc), c(., "DAPI")) ),
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.clip = 'off'
  ) +
  # scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  coord_cartesian(clip = "off")


p_dot_annotation4 <-
  p_dot_annotation3 +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues")))


pdf(
  file.path(out_dir_cluster, "dot - annotation - all markers.pdf"),
  height = 6,
  width  = 25
)
p_dot_annotation4 %>% print()
p_dot_annotation3 %>% print()
dev.off()


## plot avg_expr --------
avg_expr_rna <-
  sc %>%
  LayerData("data") %>%
  as.matrix() %>%
  apply(1, function(x) tapply(x, sc$annotation, function(y) log1p(mean(expm1(y))) ) )

pheatmap::pheatmap(
  avg_expr_rna[cell_order, protein_order],
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
  scale = "column",
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)


# spatial plot -----------
p_sp_ls <-
  sc$ROI_ID %>%
  unique() %>%
  setNames(., .) %>%
  map(function(roi){ # browser()
    meta <- sc[[]] %>% dplyr::filter(ROI_ID == roi)


    meta %>%
      ggplot() +
      geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = annotation)) +
      # annotate(
      #   "text",
      #   label = roi,
      #   x = mean((meta$XMin + meta$XMax) * 0.5),
      #   y = max((meta$YMin + meta$YMax) * 0.5),
      #   vjust = 2
      # ) +
      ggtitle(roi) +
      ggsci::scale_fill_igv() +
      scale_y_reverse() +
      coord_equal(clip = 'off') +
      theme_void() +
      theme(plot.margin = margin(5, 5, 5, 5, 'mm'))
  })

pdf(
  file.path(out_dir_cluster, "spatial plot of annotation - saplit by sample.pdf"),
  height = 4,
  width  = 6
)
p_sp_ls %>% walk(print)
dev.off()


p_sp <-
  sc[[]] %>% #  head(10000) %>%
  ggplot() +
  geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = annotation)) +
  geom_text(
    data = ~ summarise(
      group_by(., ROI_ID),
      x = mean((XMin + XMax) * 0.5),
      y = max((YMin + YMax) * 0.5)
    ),
    aes(x = x, y = y, label = ROI_ID),
    show.legend = FALSE,
    vjust = 1
  ) +
  # annotate(
  #   "text",
  #   label = roi,
  #   x = mean((meta$XMin + meta$XMax) * 0.5),
  #   y = max((meta$YMin + meta$YMax) * 0.5),
  #   vjust = 2
  # ) +
  ggsci::scale_fill_igv() +
  scale_y_reverse() +
  coord_equal(clip = 'off') +
  theme_void()

ggsave(
  p_sp,
  filename = file.path(out_dir_cluster, "spatial plot of annotation.pdf"),
  height = 22,
  width  = 10,
  limitsize = FALSE
)


# export
saveRDS(sc, str_glue("sc-dim{ndims}.rds"))


# export cell porps table ----------------

sc[[]] %>%
  group_by(ID, Group, ROI_ID, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  dplyr::select(-count) %>%
  pivot_wider(
    names_from = "annotation",
    values_from = "freq"
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_cluster, "cell props.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

sc[[]] %>%
  group_by(ID, Group, ROI_ID, annotation) %>%
  summarise(count = n()) %>%

  pivot_wider(
    names_from = "annotation",
    values_from = "count"
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_cluster, "cell count.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )














