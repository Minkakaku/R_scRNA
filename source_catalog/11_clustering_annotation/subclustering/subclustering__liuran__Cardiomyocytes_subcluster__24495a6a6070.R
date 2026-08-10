library(tidyverse)
library(Seurat)
library(patchwork)

# Cardiomyocytesblast subcluster ---------
cm_out_dir <- "Cardiomyocytes_subcluster" # dir.create(cm_out_dir)
dir.create(cm_out_dir)

sc_Cardiomyocytes <- all_data[, all_data$cell_type_res05 == "Cardiomyocytes"]
DefaultAssay(sc_Cardiomyocytes) <- "RNA"
sc_Cardiomyocytes[['sketch']] <- NULL


sc_Cardiomyocytes <-
  sc_Cardiomyocytes %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_Cardiomyocytes <- sc_Cardiomyocytes %>% RunPCA()

pdf(file.path(cm_out_dir, 'pca plot-Cardiomyocytes.pdf'), height = 6, width = 8)
ElbowPlot(sc_Cardiomyocytes, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_Cardiomyocytes <- 
  sc_Cardiomyocytes[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_Cardiomyocytes) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(cm_out_dir, 'Cardiomyocytes dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_Cardiomyocytes[, suset_cells_Cardiomyocytes], dims = 1:5)
dev.off()
pdf(file.path(cm_out_dir, 'Cardiomyocytes dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_Cardiomyocytes[, suset_cells_Cardiomyocytes], dims = 6:10)
dev.off()
pdf(file.path(cm_out_dir, 'Cardiomyocytes dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_Cardiomyocytes[, suset_cells_Cardiomyocytes], dims = 11:15)
dev.off()
pdf(file.path(cm_out_dir, 'Cardiomyocytes dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_Cardiomyocytes[, suset_cells_Cardiomyocytes], dims = 16:20)
dev.off()


## dimension reduction --------
Cardiomyocytes_dims <- 10 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_Cardiomyocytes <- RunUMAP(sc_Cardiomyocytes, dims = 1:Cardiomyocytes_dims)


pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes, group.by = "project", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} project - split.pdf")), height = 8, width = 16)
DimPlot(sc_Cardiomyocytes, group.by = "project", ncol = 4, split.by = "project") + NoLegend()
dev.off()


sc_Cardiomyocytes <- sc_Cardiomyocytes %>% 
  harmony::RunHarmony(
    dims.use = 1:Cardiomyocytes_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_Cardiomyocytes <- RunUMAP(
  sc_Cardiomyocytes, 
  dims = 1:Cardiomyocytes_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} project - harmony - split.pdf")), height = 8, width = 16)
DimPlot(sc_Cardiomyocytes, group.by = "project", reduction = "harmonyumap", ncol = 4, split.by = "project") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_Cardiomyocytes <-
  FindNeighbors(
    sc_Cardiomyocytes,
    reduction = reduction.use,
    dims = 1:Cardiomyocytes_dims,
    k.param = 10
  )

sc_Cardiomyocytes <- 
  FindClusters(
    sc_Cardiomyocytes,
    resolution = c(0.1, 0.2, 0.3, 0.4) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_Cardiomyocytes, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_Cardiomyocytes) <- "RNA_snn_res.0.4"
DefaultAssay(sc_Cardiomyocytes) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(sc_Cardiomyocytes),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# all_markers <-
#   FindAllConservedMarkers(
#     JoinLayers(sc_Cardiomyocytes),
#     grouping.var = "project",
#     test.use = "wilcox",
#     # logfc.threshold = 0.25,
#     # min.pct = 0.1,
#     only.pos = TRUE
#   ) %>%
#   map(rownames_to_column, "gene") %>%
#   enframe("cluster") %>%
#   unnest(value) %>%
#   rowwise() %>%
#   mutate(
#     avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
#     pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
#   ) %>%
#   mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_Cardiomyocytes))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(cm_out_dir, "all_marker_table.xlsx"),
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
    sc_Cardiomyocytes, 
    features = all_top5_gene
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(cm_out_dir, "top 5 gene dotplot.pdf"),
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

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t() %>%
  pheatmap::pheatmap(
    filename = "Cardiomyocytes_subcluster/all_markers_heatmap.pdf", width = 20, height = 5,
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )





major_marker_ls <- c(
  "EC" = "Cdh5",
  "EC" = "Pecam1",
  "FB" = "Col1a2",
  "FB" = "Col1a1",
  "FB" = "Pdgfra",
  "SMC" = "Tagln",
  "SMC" = "Myh11",
  "TC" = "Cd3e",
  "TC" = "Cd4",
  "TC" = "Cd8a",
  "BC" = "Cd79a",
  "BC" = "Ms4a1",
  "Cardiomyocytes" = "Tnni3k",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Pln",
  "Cardiomyocytes" = "Myh7",
  "Cardiomyocytes" = "Actn2",
  "MAC" = "Adgre1",
  "MAC" = "Fcgr1",
  "MAC" = "Csf1r",
  "Endocardial" = "Plagl1",
  "Endocardial" = "Nrg1",
  "Endocardial" = "Plvap",
  "Endocardial" = "Smoc1",
  "LECs" = "Ccl21a",
  "LECs" = "Ccl21b",
  "LECs" = "Mmrn1",
  "LECs" = "Reln",
  "LECs" = "Gfra3",
  "LECs" = "Prox1",
  "LECs" = "Lyve1",
  "Schwann cells" = "Sox10",
  "Schwann cells" = "Cnp",
  "Schwann cells" = "Cd59a",
  "Epicardial" = "Wt1",
  "Epicardial" = "Tbx18",
  "Epicardial" = "Aldh1a2",
  "Neutrophil" = "Csf3r",
  "Neutrophil" = "Cxcr2",
  "Neutrophil" = "Ccr1",
  "Pericytes" = "Pdgfrb",
  "Pericytes" = "Cspg4",
  "Pericytes" = "Mcam"
)



avg_expr[rownames(avg_expr) %in% Cardiomyocytes_marker_ls2, ] %>% 
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

saveRDS(sc_Cardiomyocytes, file = file.path(cm_out_dir, "sc_Cardiomyocytes.rds"))


# dotplot
pdf(file.path(cm_out_dir, str_glue("dotplot dims{Cardiomyocytes_dims} {reduction.use} Cardiomyocytes_marker_ls1.pdf")), width = 15, height = 6)
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", features = Cardiomyocytes_marker_ls1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", features = Cardiomyocytes_marker_ls1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", features = Cardiomyocytes_marker_ls1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", features = Cardiomyocytes_marker_ls1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


pdf(file.path(cm_out_dir, str_glue("dotplot dims{Cardiomyocytes_dims} {reduction.use} Cardiomyocytes_marker_ls2.pdf")), width = 8, height = 6)
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", features = Cardiomyocytes_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", features = Cardiomyocytes_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", features = Cardiomyocytes_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", features = Cardiomyocytes_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


pdf(file.path(cm_out_dir, str_glue("dotplot dims{Cardiomyocytes_dims} {reduction.use} major_marker_ls.pdf")), width = 15, height = 6)
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

pdf(file.path(cm_out_dir, str_glue("dotplot dims{Cardiomyocytes_dims} {reduction.use} major_marker_ls_1.pdf")), width = 15, height = 6)
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.1", features = major_marker_ls_1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.2", features = major_marker_ls_1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.3", features = major_marker_ls_1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", features = major_marker_ls_1) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()



sample_map_group <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "")
sc_Cardiomyocytes@meta.data$project <- factor(sc_Cardiomyocytes@meta.data$project, levels = sample_map_group[,1])


# cell props 
sc_Cardiomyocytes_cell_props_df <-
  sc_Cardiomyocytes[[]] %>% 
  dplyr::select(project, group, RNA_snn_res.0.4) %>%
  group_by(project, group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_Cardiomyocytes_props_project <-
  sc_Cardiomyocytes_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_igv()

p_sc_Cardiomyocytes_props_group <-
  sc_Cardiomyocytes_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_igv()

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes - cell props baplot.pdf")), height = 4, width = 8)
p_sc_Cardiomyocytes_props_project %>% print()
p_sc_Cardiomyocytes_props_group   %>% print()
dev.off()

saveRDS(sc_Cardiomyocytes, file = file.path(cm_out_dir, "sc_Cardiomyocytes.rds"))

## dotplot ----------
# FeaturePlot(
#   sc_Cardiomyocytes,
#   features = str_to_title(Cardiomyocytes_marker_ls1),
#   reduction = "harmonyumap"
# )


# DotPlot(
#   sc_Cardiomyocytes,
#   features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
#   group.by = "RNA_snn_res.0.6"
# )


# 剔除 cluster 0,4,6,7,8,9
sc_Cardiomyocytes_1 <- sc_Cardiomyocytes[, !sc_Cardiomyocytes@meta.data$RNA_snn_res.0.4 %in% c("0", "4", "6", "7", "8", "9")]


pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_Cardiomyocytes_1, file = file.path(cm_out_dir, "sc_Cardiomyocytes_1.rds"))


# annotation ---- res.0.4细胞注释
Idents(sc_Cardiomyocytes_1) <- "RNA_snn_res.0.4"
cell_cluster <- c(
  '0'  = "mCAF",
  '1'  = "iCAF",
  '3'  = "iCAF",
  '4'  = "iCAF",
  '5'  = "mCAF",
  '6'  = "dCAF",
  '9' = "infCAF"
  )

sc_Cardiomyocytes_1 <- RenameIdents(sc_Cardiomyocytes_1, cell_cluster)
sc_Cardiomyocytes_1$annotation_fb <- Idents(sc_Cardiomyocytes_1)


pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} annotation_fb {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes_1, group.by = "annotation_fb", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_Cardiomyocytes_1, group.by = "annotation_fb", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(sc_Cardiomyocytes_1) <- "RNA"
Idents(sc_Cardiomyocytes_1) <- "annotation_fb"

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes dotplot dims{Cardiomyocytes_dims} {reduction.use} annotation_fb major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_Cardiomyocytes_1, group.by = "annotation_fb", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes dotplot dims{Cardiomyocytes_dims} {reduction.use} annotation_fb Cardiomyocytes_marker_ls1.pdf")), width = 15, height = 6)
DotPlot(sc_Cardiomyocytes_1, group.by = "annotation_fb", features = Cardiomyocytes_marker_ls1, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

pdf(file.path(cm_out_dir, str_glue("Cardiomyocytes dotplot dims{Cardiomyocytes_dims} {reduction.use} annotation_fb Cardiomyocytes_marker_ls2.pdf")), width = 8, height = 6)
DotPlot(sc_Cardiomyocytes_1, group.by = "annotation_fb", features = Cardiomyocytes_marker_ls2, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


### 保存data
saveRDS(sc_Cardiomyocytes_1, file = file.path(cm_out_dir, "sc_Cardiomyocytes_1.rds"))



### RNA_snn_res.0.3 和 RNA_snn_res.0.4 FindAllMarkers
Idents(sc_Cardiomyocytes) <- "RNA_snn_res.0.3"
DefaultAssay(sc_Cardiomyocytes) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(sc_Cardiomyocytes),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
# all_markers <-
#   FindAllConservedMarkers(
#     JoinLayers(sc_Cardiomyocytes),
#     grouping.var = "project",
#     test.use = "wilcox",
#     # logfc.threshold = 0.25,
#     # min.pct = 0.1,
#     only.pos = TRUE
#   ) %>%
#   map(rownames_to_column, "gene") %>%
#   enframe("cluster") %>%
#   unnest(value) %>%
#   rowwise() %>%
#   mutate(
#     avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
#     pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
#   ) %>%
#   mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_Cardiomyocytes))[[1]]

# export all_markers
all_markers %>% 
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(cm_out_dir, "all_marker_table_res.0.3.xlsx"),
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
    sc_Cardiomyocytes, 
    features = unique(all_top5_gene)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(cm_out_dir, "top 5 gene dotplot res.0.3.pdf"),
  height = 6,
  width  = 20
)



Idents(sc_Cardiomyocytes) <- "RNA_snn_res.0.4"
DefaultAssay(sc_Cardiomyocytes) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(sc_Cardiomyocytes),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
# all_markers <-
#   FindAllConservedMarkers(
#     JoinLayers(sc_Cardiomyocytes),
#     grouping.var = "project",
#     test.use = "wilcox",
#     # logfc.threshold = 0.25,
#     # min.pct = 0.1,
#     only.pos = TRUE
#   ) %>%
#   map(rownames_to_column, "gene") %>%
#   enframe("cluster") %>%
#   unnest(value) %>%
#   rowwise() %>%
#   mutate(
#     avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
#     pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
#   ) %>%
#   mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_Cardiomyocytes))[[1]]

# export all_markers
all_markers %>% 
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(cm_out_dir, "all_marker_table_res.0.4.xlsx"),
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
    sc_Cardiomyocytes, 
    features = unique(all_top5_gene)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(cm_out_dir, "top 5 gene dotplot res.0.4.pdf"),
  height = 6,
  width  = 20
)





### cell proportion ------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")



sc_Cardiomyocytes_1 <- readRDS(file.path(cm_out_dir, "sc_Cardiomyocytes_1.rds"))

fb_ratio <- file.path(cm_out_dir, "Ratio")
dir.create(fb_ratio)

out_gene_umi=as.matrix(sc_Cardiomyocytes_1@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi, file = file.path(fb_ratio, "meta_data.txt"), sep="\t", quote=F, row.names=F, col.names=F)


data <- read.table(file.path(fb_ratio, "meta_data.txt"), header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,5,6,42)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

# data <- data %>% filter(Annotation != "Unknown")


mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file = file.path(fb_ratio, "Cell_type_stat.txt"), quote=F, row.names=T, col.names=F, sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file = file.path(fb_ratio, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("Cardiomyocytesblast_subcluster/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="Cardiomyocytesblast_subcluster/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file=file.path(fb_ratio, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file=file.path(fb_ratio, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) + # 两张图各自的宽度
  plot_annotation(
    title = "Cell Type Ratio",
    theme = theme(plot.title = element_text(
      size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
    ))
  )
ggsave(p,file=file.path(fb_ratio, "ratio_barplot.pdf"),width=12, height=6)
ggsave(p,file=file.path(fb_ratio, "ratio_barplot.png"),width=12, height=6,dpi=300)





# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = file.path(fb_ratio, "Cell_type_stat_percent.txt"), header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"
data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "SampleID"))

sample<-read.table(file="sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by=c("Sample" = "SampleID"))

# 限制分组顺序
plot_data$Group <- factor(plot_data$Group, levels=unique(sample_map[,2]))

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="Cardiomyocytesblast_subcluster/Ratio/diff_boxplot.pdf",width=10,height=5)
ggsave(p,file="Cardiomyocytesblast_subcluster/Ratio/diff_boxplot.png",width=10,height=5)










### RNA_snn_res.0.3 --- annotation ------------
Idents(sc_Cardiomyocytes) <- "RNA_snn_res.0.3"
cell_cluster <- c(
  '0'  = "f1",
  '1'  = "f2",
  '2'  = "f3",
  '3'  = "f4",
  '4'  = "f5",
  '5'  = "f6",
  '6'  = "f7",
  '7'  = "f8",
  '8'  = "f9",
  '9'  = "f10"
)

sc_Cardiomyocytes <- RenameIdents(sc_Cardiomyocytes, cell_cluster)
sc_Cardiomyocytes$annotation_fb_0918 <- Idents(sc_Cardiomyocytes)


pdf(file.path(modify_out_dir, str_glue("Cardiomyocytes Dimplot dim{Cardiomyocytes_dims} annotation_fb_0918 {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_Cardiomyocytes, group.by = "annotation_fb_0918", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_Cardiomyocytes, group.by = "annotation_fb_0918", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### 左下角小标志umap -------
modify_out_dir <- "Cardiomyocytesblast_subcluster/modify_0918"
dir.create(modify_out_dir)
# 无边框
p2 <-DimPlot(sc_Cardiomyocytes, group.by = "annotation_fb_0918", pt.size=0.01, reduction="harmonyumap", raster=FALSE, cols=plan) + 
  scale_x_continuous("") + scale_y_continuous("") + 
  theme(
    axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "p2") # 标题

# 左下角小坐标
xstart <- min(sc_Cardiomyocytes@reductions$harmonyumap@cell.embeddings[,1]) - 1
xend <- min(sc_Cardiomyocytes@reductions$harmonyumap@cell.embeddings[,1]) + 3
ystart <- min(sc_Cardiomyocytes@reductions$harmonyumap@cell.embeddings[,2]) - 1
yend <- min(sc_Cardiomyocytes@reductions$harmonyumap@cell.embeddings[,2]) + 3

p3 <- p2 +
  geom_segment(aes(x = xstart, y = ystart , xend = xend, yend = ystart),
               colour = "black", arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = xstart, y = ystart, xend = xstart , yend = yend),
               colour = "black", arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = xstart +1.5, y = ystart -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = xstart -1, y = ystart + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)+
  theme(
    legend.position = "right", #去掉图例
  )+
  labs(title="Cardiomyocytesblast")

ggsave(p3, file = file.path(modify_out_dir, "umap_annot.pdf"), width = 8, height = 6)
ggsave(p3, file = file.path(modify_out_dir, "umap_annot.png"), width = 8, height = 6, dpi=300)






### Fig c
out_gene_umi=as.matrix(sc_Cardiomyocytes@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi, file = file.path(modify_out_dir, "meta_data.txt"), sep="\t", quote=F, row.names=F, col.names=F)


data <- read.table(file.path(modify_out_dir, "meta_data.txt"), header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,5,6,42)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

# data <- data %>% filter(Annotation != "Unknown")


mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file = file.path(modify_out_dir, "Cell_type_stat.txt"), quote=F, row.names=T, col.names=F, sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file = file.path(modify_out_dir, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("Cardiomyocytesblast_subcluster/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="Cardiomyocytesblast_subcluster/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file=file.path(modify_out_dir, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file=file.path(modify_out_dir, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

# 限制细胞顺序
data4plot1$Celltype <- factor(data4plot1$Celltype, levels = paste0("f", 1:10))
data4plot2$Celltype <- factor(data4plot2$Celltype, levels = paste0("f", 1:10))

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) + # 两张图各自的宽度
  plot_annotation(
    title = "Cell Type Ratio",
    theme = theme(plot.title = element_text(
      size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
    ))
  )
ggsave(p,file=file.path(modify_out_dir, "ratio_barplot.pdf"),width=12, height=6)
ggsave(p,file=file.path(modify_out_dir, "ratio_barplot.png"),width=12, height=6,dpi=300)



p2 <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))
ggsave(p2,file=file.path(modify_out_dir, "ratio_barplot_group.pdf"),width=8, height=6)
ggsave(p2,file=file.path(modify_out_dir, "ratio_barplot_group.png"),width=8, height=6,dpi=300)


### 保存data
saveRDS(sc_Cardiomyocytes, file = file.path(modify_out_dir, "sc_Cardiomyocytes_0918.rds"))




### 小提琴图

p9 <- VlnPlot(sc_Cardiomyocytes,
              pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
              stack = F, # 是否堆叠
              flip = F, # 翻转
              #adjust=5, # 调整粗细
              fill.by = "ident",
              cols = plan,
              features = "Thbs4")+
  labs(x="")+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12,color="black",angle=90,vjust=0.5),
        strip.text.y = element_text(hjust=0))
ggsave(p9, file = file.path(modify_out_dir, "vlnplot_annot.pdf"), width = 6, height = 4)
ggsave(p9, file = file.path(modify_out_dir, "vlnplot_annot.png"), width = 6, height = 4, dpi = 300)

p10 <- VlnPlot(sc_Cardiomyocytes,
               pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
               stack = F, # 是否堆叠
               flip = F, # 翻转
               #adjust=5, # 调整粗细
               fill.by = "ident",
               group.by = "group", 
               cols = plan,
               features = "Thbs4")+
  labs(x="")+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12,color="black",angle=90,vjust=0.5),
        strip.text.y = element_text(hjust=0))
ggsave(p10, file = file.path(modify_out_dir, "vlnplot_annot_group.pdf"), width = 5, height = 4)
ggsave(p10, file = file.path(modify_out_dir, "vlnplot_annot_group.png"), width = 5, height = 4, dpi = 300)









### FindAllMarkers --------

Idents(sc_Cardiomyocytes) <- "RNA_snn_res.0.3"
DefaultAssay(sc_Cardiomyocytes) <- "RNA"



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
  JoinLayers(sc_Cardiomyocytes), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE
)

write.csv(diffmarkers, file.path(modify_out_dir, 'sc_Cardiomyocytes/cluster diffmarker.csv'))


celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "f1")


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
write.csv(enrichKEGG_ls@result,file = file.path(modify_out_dir, 'sc_Cardiomyocytes/KEGG enriment table.csv'),row.names = FALSE)

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
  filename = file.path(modify_out_dir, 'sc_Cardiomyocytes/KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(modify_out_dir, 'sc_Cardiomyocytes/KEGG barplot.pdf'),
  height = 6,
  width = 8)





celltype_diffmarkers <- diffmarkers %>%
  filter(cluster == "f2")


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
write.csv(enrichKEGG_ls@result,file = file.path(modify_out_dir, 'sc_Cardiomyocytes/f2 KEGG enriment table.csv'),row.names = FALSE)

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
  filename = file.path(modify_out_dir, 'sc_Cardiomyocytes/f2 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(modify_out_dir, 'sc_Cardiomyocytes/f2 KEGG barplot.pdf'),
  height = 6,
  width = 8)







# feature plot of thbs4
p_thbs4 <- FeaturePlot(
  sc_Cardiomyocytes, 
  raster = F, 
  features = "Thbs4", 
  reduction = "harmonyumap", 
  pt.size = 0.05, 
  label = T
) + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
    legend.position = "none", 
    plot.subtitle = element_text(hjust = 0.05, vjust = -9, size = 15, color = "red"), 
    plot.title = element_text(vjust = -5, hjust = 0.5, size = 17)
  ) + 
  labs(title = "Cardiomyocytesblast", subtitle = "Thbs4") + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))

ggsave(p_thbs4, file = file.path(cm_out_dir, "Thbs4_umap.pdf"), width = 8, height = 7)
ggsave(p_thbs4, file = file.path(cm_out_dir, "Thbs4_umap.png"), width = 8, height = 7, dpi = 300)














