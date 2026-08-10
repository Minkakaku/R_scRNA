library(tidyverse)
library(Seurat)
library(patchwork)

# AT2 subcluster ---------
AT2_out_dir <- "AT2_subcluster" # dir.create(AT2_out_dir)
dir.create(AT2_out_dir)

sc_AT2 <- all_data[, all_data$annotation_full == "AT2"]
DefaultAssay(sc_AT2) <- "RNA"
sc_AT2[['sketch']] <- NULL


sc_AT2 <-
  sc_AT2 %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_AT2 <- sc_AT2 %>% RunPCA()

pdf(file.path(AT2_out_dir, 'pca plot-AT2.pdf'), height = 6, width = 8)
ElbowPlot(sc_AT2, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_AT2 <- 
  sc_AT2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_AT2) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(AT2_out_dir, 'AT2 dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_AT2[, suset_cells_AT2], dims = 1:5)
dev.off()
pdf(file.path(AT2_out_dir, 'AT2 dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_AT2[, suset_cells_AT2], dims = 6:10)
dev.off()
pdf(file.path(AT2_out_dir, 'AT2 dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_AT2[, suset_cells_AT2], dims = 11:15)
dev.off()
pdf(file.path(AT2_out_dir, 'AT2 dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_AT2[, suset_cells_AT2], dims = 16:20)
dev.off()


## dimension reduction --------
AT2_dims <- 12 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_AT2 <- RunUMAP(sc_AT2, dims = 1:AT2_dims)


pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_AT2, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

# pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} project - split.pdf")), height = 8, width = 16)
# DimPlot(sc_AT2, group.by = "project", label = TRUE, ncol = 4, split.by = "project") + NoLegend()
# dev.off()


sc_AT2 <- sc_AT2 %>% 
  harmony::RunHarmony(
    dims.use = 1:AT2_dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_AT2 <- RunUMAP(
  sc_AT2, 
  dims = 1:AT2_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_AT2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} project - harmony - split.pdf")), height = 8, width = 16)
# DimPlot(sc_AT2, group.by = "project", label = TRUE, reduction = "harmonyumap", ncol = 4, split.by = "project") + NoLegend()
# dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_AT2 <-
  FindNeighbors(
    sc_AT2,
    reduction = reduction.use,
    dims = 1:AT2_dims,
    k.param = 10
  )

sc_AT2 <- 
  FindClusters(
    sc_AT2,
    resolution = c(0.05, 0.1, 0.2, 0.3, 0.4) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_AT2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_AT2, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_AT2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_AT2, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_AT2) <- "RNA_snn_res.0.3"
DefaultAssay(sc_AT2) <- "RNA"

# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_AT2),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
all_markers <-
  FindAllConservedMarkers(
    JoinLayers(sc_AT2),
    grouping.var = "sample",
    test.use = "wilcox",
    # logfc.threshold = 0.25,
    # min.pct = 0.1,
    only.pos = TRUE
  ) %>%
  map(rownames_to_column, "gene") %>%
  enframe("cluster") %>%
  unnest(value) %>%
  rowwise() %>%
  mutate(
    avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
    pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
  ) %>%
  mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_AT2))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(AT2_out_dir, "all_marker_table.xlsx"),
    asTable = TRUE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc_AT2, 
    features = all_top5_gene
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(AT2_out_dir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 20
)



# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  # arrange(p_val_adj) %>% 
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 20) %>% 
  nest() %>% 
  deframe() %>% 
  map(pull, 'gene') %>% 
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>% 
  walk(cat, "\n")

# avg_expr[all_de_genes, ] %>% 
#   as.matrix() %>% 
#   as.data.frame() %>% 
#   setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t() %>%
#   pheatmap::pheatmap(
#     filename = "AT2_subcluster/all_markers_heatmap.pdf", width = 20, height = 5,
#     scale = "column",
#     cluster_row = avg_expr[VariableFeatures(sc_AT2), ] %>% t %>% dist %>% hclust,
#     cluster_col = FALSE,
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_col = NA
#   )




AT2_marker_ls1 <- c(
  'mCAF'='Col1a1',
  'mCAF'='Postn',
  'mCAF'='Col10a1',
  'mCAF'='Cdh11',
  'mCAF'='Lum',
  'iCAF'='Cd34',
  'iCAF'='Cxcl14',
  'iCAF'='Cxcl1',
  'iCAF'='Ccl2',
  'apCAF'='Igkc',
  'apCAF'='Cxcr4',
  'vCAF'='Mcam',
  'vCAF'='Myh11',
  'tCAF'='Eno1',
  'tCAF'='Mme',
  'dCAF'='Mki67',
  'dCAF'='Pttg1',
  'dCAF'='H2afz',
  'dCAF'='Cenpf',
  'infCAF'='Ido1',
  'infCAF'='Cxcl9',
  'infCAF'='Cxcl10',
  'rCAF'='Ccl19',
  'pericyte'='Col4a1',
  'pericyte'='Col4a2',
  'hsptcaf'='Hsph1',
  'hsptcaf'='Vegfa'
)
# AT2_marker_ls1 <- AT2_marker_ls1 %>% str_to_title


major_marker_ls <- c(
  "Epithelial"  = "Epcam",
  AT1 = "Ager",
  AT1 = "Pdpn",
  AT1 = "Csf1r",
  AT1 = "Mafb",
  AT1 = "Batf3",
  AT1 = "Krt79",
  AT1 = "Msr1",
  AT2 = "Abca3",
  AT2 = "Wif1",
  AT2 = "Etv5",
  "Basel cell" = "Ngfr",
  "B cell" = "Cd79a",
  "B cell" = "Cd79b",
  "B cell" = "Cd19",
  "B cell" = "Ms4a1",
  "T cells" = "Cd3d",
  "T cells" = "Cd3e",
  "T cells" = "Cd3g",
  "CD4T cell" = "Cd4",
  "CD4T cell" = "Cxcl13",
  "CD4T cell" = "Il21",
  "CD8T cell" = "Cd2",
  "CD8T cell" = "Cd8a",
  "CD8T cell" = "Cd8b1",
  "CD8T cell" = "Il2ra",
  "CD8T cell" = "Ptprc",
  Cilia = "Gsta1",
  Cilia = "Pifo",
  Cilia = "Foxj1",
  # "Club cell" = "PECAM1",
  Cilliated = "Muc5ac",
  "Dendritic cell" = "Ccl17",
  "Dendritic cell" = "Cd83",
  "Dendritic cell" = "Clec9a",
  "Dendritic cell" = "Ido1",
  Endothelial = "Pecam1",
  Endothelial = "Cdh5",
  Endothelial = "Vwf",
  Endothelial = "Plvap",
  Fibroblast = "Col1a1", 
  Fibroblast = "Col1a2",
  Fibroblast = "Pdgfra",
  Fibroblast = "Lum",
  Granulocyte = "Csf3r",
  "Mono/Mac" = "Cd68",
  "Mono/Mac" = "Fcgr1",
  "Mono/Mac" = "Cd163",
  "Mono/Mac" = "Cd14",
  "Mono/Mac" = "Fcgr3",
  "Mono/Mac" = "Vcan",
  "Mono/Mac" = "Mrc1",
  "Mast cell" = "Gata2",
  "Mast cell" = "Kit",
  "Mast cell" = "Ms4a2",
  "Mast cell" = "Tpsab1",
  "Mast cell" = "Tpsb2",
  "NK cell" = "Ncam1",
  "NK cell" = "Klrd1",
  Pericytes = "Pdgfrb",
  Pericytes = "Mcam",
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Tagln"
)


avg_expr[rownames(avg_expr) %in% major_marker_ls, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )


# dotplot
pdf(file.path(AT2_out_dir, str_glue("dotplot dims{AT2_dims} {reduction.use} major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(sc_AT2, group.by = "RNA_snn_res.0.05", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_AT2, group.by = "RNA_snn_res.0.1", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_AT2, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_AT2, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_AT2, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


saveRDS(sc_AT2, file = file.path(AT2_out_dir, "sc_AT2.rds"))

## remove 6 7 8--------
sc_AT2 <- sc_AT2[, !sc_AT2$RNA_snn_res.0.3 %in% c(2, 3, 6, 8)]

pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_AT2, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_AT2, file = file.path(AT2_out_dir, "sc_AT2_trimmed.rds"))


sample_map_group <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "")
sc_AT2@meta.data$sample <- factor(sc_AT2@meta.data$sample, levels = sample_map_group[,1])


# cell props 
sc_AT2_cell_props_df <- 
  sc_AT2[[]] %>% 
  dplyr::select(sample, group, RNA_snn_res.0.3) %>% 
  group_by(sample, group, RNA_snn_res.0.3) %>% 
  summarise(count = n()) %>% 
  group_by(sample) %>% 
  mutate(freq_project = count/sum(count)) %>% 
  group_by(group) %>% 
  mutate(freq_group = count/sum(count))

p_sc_AT2_props_project <-
  sc_AT2_cell_props_df %>%
  ggplot(aes(x = sample, y = freq_project, fill = RNA_snn_res.0.3)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_igv()

p_sc_AT2_props_group <-
  sc_AT2_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.3)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_igv()

pdf(file.path(AT2_out_dir, str_glue("AT2 - cell props baplot.pdf")), height = 4, width = 6)
p_sc_AT2_props_project %>% print()
p_sc_AT2_props_group   %>% print()
dev.off()

# saveRDS(sc_AT2, file = file.path(AT2_out_dir, "sc_AT2.rds"))

## dotplot ----------
# FeaturePlot(
#   sc_AT2,
#   features = str_to_title(AT2_marker_ls1),
#   reduction = "harmonyumap"
# )


# DotPlot(
#   sc_AT2,
#   features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
#   group.by = "RNA_snn_res.0.6"
# )




# spatial plot --------------
AT2_meta <-
  sc_AT2[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(AT2_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = RNA_snn_res.0.3, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(AT2_out_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })







# # 剔除 cluster 2,7,8,10
# sc_AT2_1 <- sc_AT2[, sc_AT2@meta.data$RNA_snn_res.0.4 %in% c("0", "1", "3", "4", "5", "6", "9")]
# 
# 
# # annotation ----
# # annotation ---- res.0.4细胞注释
# Idents(sc_AT2_1) <- "RNA_snn_res.0.3"
# cell_cluster <- c(
#   '0'  = "mCAF",
#   '1'  = "iCAF",
#   '3'  = "iCAF",
#   '4'  = "iCAF",
#   '5'  = "mCAF",
#   '6'  = "dCAF",
#   '9' = "infCAF"
# )
# 
# sc_AT2_1 <- RenameIdents(sc_AT2_1, cell_cluster)
# sc_AT2_1$annotation_fb <- Idents(sc_AT2_1)
# 
# 
# pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} annotation_fb {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_AT2_1, group.by = "annotation_fb", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(sc_AT2_1, group.by = "annotation_fb", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# dev.off()
# 
# DefaultAssay(sc_AT2_1) <- "RNA"
# Idents(sc_AT2_1) <- "annotation_fb"
# 
# pdf(file.path(AT2_out_dir, str_glue("AT2 dotplot dims{AT2_dims} {reduction.use} annotation_fb major_marker.pdf")), width = 25, height = 6)
# DotPlot(sc_AT2_1, group.by = "annotation_fb", features = major_marker_ls, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()
# 
# pdf(file.path(AT2_out_dir, str_glue("AT2 dotplot dims{AT2_dims} {reduction.use} annotation_fb AT2_marker_ls1.pdf")), width = 15, height = 6)
# DotPlot(sc_AT2_1, group.by = "annotation_fb", features = AT2_marker_ls1, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()
# 
# pdf(file.path(AT2_out_dir, str_glue("AT2 dotplot dims{AT2_dims} {reduction.use} annotation_fb AT2_marker_ls2.pdf")), width = 8, height = 6)
# DotPlot(sc_AT2_1, group.by = "annotation_fb", features = AT2_marker_ls2, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()
# 
# 
# ### 保存data
# saveRDS(sc_AT2_1, file = file.path(AT2_out_dir, "sc_AT2_1.rds"))
# 
# 
# 
# ### RNA_snn_res.0.3 和 RNA_snn_res.0.4 FindAllMarkers
# Idents(sc_AT2) <- "RNA_snn_res.0.3"
# DefaultAssay(sc_AT2) <- "RNA"
# 
# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_AT2),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# # all_markers <-
# #   FindAllConservedMarkers(
# #     JoinLayers(sc_AT2),
# #     grouping.var = "project",
# #     test.use = "wilcox",
# #     # logfc.threshold = 0.25,
# #     # min.pct = 0.1,
# #     only.pos = TRUE
# #   ) %>%
# #   map(rownames_to_column, "gene") %>%
# #   enframe("cluster") %>%
# #   unnest(value) %>%
# #   rowwise() %>%
# #   mutate(
# #     avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
# #     pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
# #   ) %>%
# #   mutate(p_val_adj  = minimump_p_val)
# 
# all_de_genes <-
#   all_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   # dplyr::filter(pct.1 > 0.25) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 5) %>%
#   pull("gene")
# 
# # DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
# if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_AT2))[[1]]
# 
# # export all_markers
# all_markers %>% 
#   as.data.frame() %>%
#   openxlsx::write.xlsx(
#     file = file.path(AT2_out_dir, "all_marker_table_res.0.3.xlsx"),
#     asTable = TRUE
#   )
# 
# all_top5_gene <- 
#   all_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   dplyr::filter(pct.1 > 0.25) %>%
#   dplyr::filter(abs(avg_log2FC) > 1) %>% ##
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   slice_head(n = 5) %>%
#   with(structure(
#     gene,
#     names = as.character(cluster)
#   ))
# 
# # top5 dotplot
# p_dot_top5 <-
#   DotPlot(
#     sc_AT2, 
#     features = unique(all_top5_gene)
#   ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggsave(
#   p_dot_top5,
#   filename = file.path(AT2_out_dir, "top 5 gene dotplot res.0.3.pdf"),
#   height = 6,
#   width  = 20
# )
# 
# 
# 
# Idents(sc_AT2) <- "RNA_snn_res.0.4"
# DefaultAssay(sc_AT2) <- "RNA"
# 
# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_AT2),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# # all_markers <-
# #   FindAllConservedMarkers(
# #     JoinLayers(sc_AT2),
# #     grouping.var = "project",
# #     test.use = "wilcox",
# #     # logfc.threshold = 0.25,
# #     # min.pct = 0.1,
# #     only.pos = TRUE
# #   ) %>%
# #   map(rownames_to_column, "gene") %>%
# #   enframe("cluster") %>%
# #   unnest(value) %>%
# #   rowwise() %>%
# #   mutate(
# #     avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
# #     pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
# #   ) %>%
# #   mutate(p_val_adj  = minimump_p_val)
# 
# all_de_genes <-
#   all_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   # dplyr::filter(pct.1 > 0.25) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 5) %>%
#   pull("gene")
# 
# # DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
# if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_AT2))[[1]]
# 
# # export all_markers
# all_markers %>% 
#   as.data.frame() %>%
#   openxlsx::write.xlsx(
#     file = file.path(AT2_out_dir, "all_marker_table_res.0.4.xlsx"),
#     asTable = TRUE
#   )
# 
# all_top5_gene <- 
#   all_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   dplyr::filter(pct.1 > 0.25) %>%
#   dplyr::filter(abs(avg_log2FC) > 1) %>% ##
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   slice_head(n = 5) %>%
#   with(structure(
#     gene,
#     names = as.character(cluster)
#   ))
# 
# # top5 dotplot
# p_dot_top5 <-
#   DotPlot(
#     sc_AT2, 
#     features = unique(all_top5_gene)
#   ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggsave(
#   p_dot_top5,
#   filename = file.path(AT2_out_dir, "top 5 gene dotplot res.0.4.pdf"),
#   height = 6,
#   width  = 20
# )
# 
# 
# 
# 
# 
# ### cell proportion ------
# library(reshape2)
# library(patchwork)
# library(dplyr)
# library(ggplot2)
# 
# 
# plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
#           "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
#           "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
#           "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
#           "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
#           "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
#           "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
#           "#f17d80", "#737495")
# 
# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
# 
# 
# 
# sc_AT2_1 <- readRDS(file.path(AT2_out_dir, "sc_AT2_1.rds"))
# 
# fb_ratio <- file.path(AT2_out_dir, "Ratio")
# dir.create(fb_ratio)
# 
# out_gene_umi=as.matrix(sc_AT2_1@meta.data)
# row_name=c("Cell",rownames(out_gene_umi))
# out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
# out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
# write.table(out_gene_umi, file = file.path(fb_ratio, "meta_data.txt"), sep="\t", quote=F, row.names=F, col.names=F)
# 
# 
# data <- read.table(file.path(fb_ratio, "meta_data.txt"), header = T, sep = "\t", check.names = F)
# colnames(data)
# data <- data[, c(1,5,6,42)]
# colnames(data) <- c("Cells","Sample","Group","Annotation")
# 
# # data <- data %>% filter(Annotation != "Unknown")
# 
# 
# mydata<-data %>% count(Sample,Annotation)
# mymtrix<-dcast(mydata,Sample~Annotation)
# mymtrix[is.na(mymtrix)]<-0
# mymtrix_out<-t(mymtrix)
# write.table(as.data.frame(mymtrix_out),file = file.path(fb_ratio, "Cell_type_stat.txt"), quote=F, row.names=T, col.names=F, sep="\t")
# 
# rownames(mymtrix)<-mymtrix[,"Sample"]
# mymtrix<-mymtrix[,-1]
# mymtrix<-as.matrix(mymtrix)
# percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
#   mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
# }))
# rownames(percent)<-rownames(mymtrix)
# data4plot1<-melt(percent)
# percent<-t(percent)
# DATA<-cbind(rownames(percent),percent)
# DATA<-rbind(c("Sample",colnames(percent)),DATA)
# write.table(DATA,file = file.path(fb_ratio, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
# colnames(data4plot1)<-c("Sample","Celltype","Abund")
# 
# # samplemap <- read.table("fibroblast_subcluster/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# # colnames(samplemap)<-c("Sample","Group")
# # final_data <- data %>% left_join(samplemap, by = "Sample")
# 
# final_data <- data
# 
# # write.table(final_data,file="fibroblast_subcluster/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
# mydata<-final_data %>% count(Group,Annotation)
# mymtrix<-dcast(mydata,Group~Annotation)
# mymtrix[is.na(mymtrix)]<-0
# mymtrix_out<-t(mymtrix)
# write.table(as.data.frame(mymtrix_out),file=file.path(fb_ratio, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")
# 
# rownames(mymtrix)<-mymtrix[,"Group"]
# mymtrix<-mymtrix[,-1]
# mymtrix<-as.matrix(mymtrix)
# percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
#   mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
# }))
# rownames(percent)<-rownames(mymtrix)
# data4plot2<-melt(percent)
# percent<-t(percent)
# DATA<-cbind(rownames(percent),percent)
# DATA<-rbind(c("Group",colnames(percent)),DATA)
# write.table(DATA,file=file.path(fb_ratio, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
# colnames(data4plot2)<-c("Group","Celltype","Abund")
# 
# # 限制样本顺序和分组顺序
# sample_map <- read.table("sample_map_group.txt",header=F,sep="\t",quote="",check.names=F)
# data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
# data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))
# 
# plot_list <- list()
# plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
#   geom_bar(stat="identity",size=0.5,width=0.7)+
#   labs(x="Sample",y="Percent(%)",title="")+
#   scale_fill_manual(values = plan)+
#   # ggsci::scale_fill_igv() +
#   theme_classic()+
#   theme(
#     axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
#     axis.text.y = element_text(size=15,colour="black"),
#     axis.title = element_text(size=15,colour="black",face="bold"),
#     legend.title = element_text(size=15,colour="black",face="bold"),
#     plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
#     strip.background=element_rect(fill="white",size=2))+NoLegend()
# 
# plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
#   geom_bar(stat="identity",size=0.5,width=0.7)+
#   labs(x="Group",y="",title="")+
#   scale_fill_manual(values = plan)+
#   # ggsci::scale_fill_igv() +
#   theme_classic()+
#   theme(
#     axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
#     axis.text.y = element_text(size=15,colour="black"),
#     legend.text = element_text(size=15,colour="black"),
#     legend.title = element_text(size=15,colour="black",face="bold"),
#     axis.title = element_text(size=15,colour="black",face="bold"),
#     plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
#     strip.background=element_rect(fill="white",size=2))
# 
# p <- wrap_plots(plot_list,ncol=2)+
#   plot_layout(widths = c(2, 1)) + # 两张图各自的宽度
#   plot_annotation(
#     title = "Cell Type Ratio",
#     theme = theme(plot.title = element_text(
#       size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
#     ))
#   )
# ggsave(p,file=file.path(fb_ratio, "ratio_barplot.pdf"),width=12, height=6)
# ggsave(p,file=file.path(fb_ratio, "ratio_barplot.png"),width=12, height=6,dpi=300)
# 
# 
# 
# 
# 
# # 箱线图（带显著性标记）
# library(ggpubr)
# 
# # 普通版
# data<-read.table(file = file.path(fb_ratio, "Cell_type_stat_percent.txt"), header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
# colnames(data)[1]<-"Gene"
# data[,1]<-factor(data[,1])
# data<-melt(data,value.name = "value",variable.name="Sample")
# 
# rename<-read.table(file="sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
# plot_data<-left_join(data,rename,by=c("Sample" = "SampleID"))
# 
# sample<-read.table(file="sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
# colnames(sample) <- c("SampleID","Group")
# 
# plot_data<-left_join(plot_data,sample,by=c("Sample" = "SampleID"))
# 
# # 限制分组顺序
# plot_data$Group <- factor(plot_data$Group, levels=unique(sample_map[,2]))
# 
# groups<-unique(plot_data[,"Group"])
# 
# Y_Max<-max(plot_data[,"value"])*1.5
# Y_Max_2<-max(plot_data[,"value"])*1.001
# 
# p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
#   geom_boxplot(color="black")+
#   geom_jitter(position = position_dodge(0.75),size=1)+
#   stat_compare_means(label = "p.format",size=3,
#                      aes(group = Group))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
#                                   size=15,lineheight = 0.5),
#         axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
#                                  size=12),
#         axis.text.y=element_text(angle=0,colour="black",
#                                  # hjust=0.5,vjust=0.5,
#                                  size=12),
#         legend.position="right",
#         legend.title=element_text(angle=0,colour="black",face="bold",
#                                   size=15,lineheight = 0.5),
#         legend.text=element_text(angle=0,colour="black",size=12),
#         panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
#         plot.margin=unit(c(6,6,6,25),"mm"))+
#   labs(x="",y="Percent of Celltype (%)")+
#   theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
# ggsave(p,file="fibroblast_subcluster/Ratio/diff_boxplot.pdf",width=10,height=5)
# ggsave(p,file="fibroblast_subcluster/Ratio/diff_boxplot.png",width=10,height=5)




### RNA_snn_res.0.3 --- annotation ------------
Idents(sc_AT2) <- "RNA_snn_res.0.3"
cell_cluster <- c(
  '0'  = "C0",
  '1'  = "C1",
  '4'  = "C4",
  '5'  = "C5",
  '7'  = "C7"
)

sc_AT2 <- RenameIdents(sc_AT2, cell_cluster)
sc_AT2$AT2_subtype <- Idents(sc_AT2)






## 替换AT2注释 -----
out_gene_umi = as.matrix(sc_AT2@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
write.table(out_gene_umi, file = file.path(AT2_out_dir, "AT2_meta_data.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

AT2_meta_data <- read.table(file.path(AT2_out_dir, "AT2_meta_data.txt"), header = T, sep = "\t", check.names = F)
colnames(AT2_meta_data)
AT2_meta_data <- AT2_meta_data[,c(1,31)]
colnames(AT2_meta_data) <- c("barcode","AT2_subtype")


all_data@meta.data$Cells <- rownames(all_data@meta.data)
all_data@meta.data$barcode <- all_data@meta.data$Cells
all_data@meta.data <- all_data@meta.data %>%
  left_join(AT2_meta_data, by = "barcode") %>%
  mutate(
    annotation_updated = case_when(
      # 在亚群对象中有细化的注释
      !is.na(AT2_subtype) ~ AT2_subtype,
      
      # 保持非T细胞类型的原始注释
      !annotation_full %in% "AT2" ~ annotation_full,
      
      # 对于未细化的T细胞给予特殊标记
      TRUE ~ paste0(annotation_full, " (unclassified)")
    )
  ) %>%
  select(-c(barcode, AT2_subtype))  # 移除临时列
unique(all_data$annotation_updated)
all_data@meta.data$AT2_subtype <- all_data@meta.data$annotation_updated
rownames(all_data@meta.data) <- all_data@meta.data$Cells

all_data@meta.data <- all_data@meta.data[, -27]



pdf(file.path(AT2_out_dir, str_glue("AT2 Dimplot dim{AT2_dims} AT2_subtype {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_AT2, group.by = "AT2_subtype", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_AT2, group.by = "AT2_subtype", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


data <- data %>% filter(Annotation != "Unknown")

# spatial plot --------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>% 
  filter(AT2_subtype != "AT2 (unclassified)")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

cells_boundary_ls_new <- cells_boundary_ls

sample_mapping <- c(
  "BD3-1" = "Dex3_1",
  "BD3-2" = "Dex3_2",
  "B3-1" = "IPF3_1",
  "B3-2" = "IPF3_2",
  "B3-3" = "IPF3_3"
)

names(cells_boundary_ls_new) <- sample_mapping[names(cells_boundary_ls_new)]

sp_p_ls <-
  cells_boundary_ls_new %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = AT2_subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(AT2_out_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })








### FindAllMarkers --------
Idents(sc_AT2) <- "AT2_subtype"
DefaultAssay(sc_AT2) <- "RNA"

AT2_enrich_outdir <- file.path(AT2_out_dir, "enrichment analysis")
dir.create(AT2_enrich_outdir)


library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Mm.eg.db")
# library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 



diffmarkers <- FindAllMarkers(
  JoinLayers(sc_AT2), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE
)

write.csv(diffmarkers, file.path(AT2_enrich_outdir, 'cluster diffmarker.csv'))






celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C0")


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(AT2_enrich_outdir, 'C0 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(AT2_enrich_outdir, 'C0 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(AT2_enrich_outdir, 'C0 KEGG barplot.pdf'),
  height = 6,
  width = 8)





celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C1")


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(AT2_enrich_outdir, 'C1 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(AT2_enrich_outdir, 'C1 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(AT2_enrich_outdir, 'C1 KEGG barplot.pdf'),
  height = 6,
  width = 8)






celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C4")


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(AT2_enrich_outdir, 'C4 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(AT2_enrich_outdir, 'C4 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(AT2_enrich_outdir, 'C4 KEGG barplot.pdf'),
  height = 6,
  width = 8)






celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C5")


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(AT2_enrich_outdir, 'C5 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(AT2_enrich_outdir, 'C5 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(AT2_enrich_outdir, 'C5 KEGG barplot.pdf'),
  height = 6,
  width = 8)







celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "C7")


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(celltype_diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(AT2_enrich_outdir, 'C7 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(AT2_enrich_outdir, 'C7 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(AT2_enrich_outdir, 'C7 KEGG barplot.pdf'),
  height = 6,
  width = 8)






