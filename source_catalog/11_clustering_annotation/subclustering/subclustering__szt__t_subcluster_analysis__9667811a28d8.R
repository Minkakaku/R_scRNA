library(tidyverse)
library(Seurat)
library(patchwork)



# multi workers  -----
library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3)
plan("multisession", workers = 10, gc = TRUE)





# t subcluster ---------
gene_discription <- read_tsv("gene_description_biomart.txt")
ig_genes <- 
  gene_discription %>% 
  dplyr::filter(str_detect(`Gene description`, "^immunoglobulin")) %>% 
  dplyr::pull(`Gene name`) %>%
  .[!is.na(.)]

tcr_genes <-
  gene_discription %>% 
  dplyr::filter(str_detect(`Gene description`, "^T cell receptor")) %>% 
  dplyr::pull(`Gene name`) %>%
  .[!is.na(.)]

mt_genes <- str_subset(rownames(all_data), pattern = "^MT-")
ribo_genes <- str_subset(rownames(all_data), pattern = "^RP[SL]")

# table(ig_genes %in% rownames(all_data) )
# table(tcr_genes %in% rownames(all_data) )



t_out_dir <- "t_subcluster"
dir.create(t_out_dir)


sc_t <- all_data %>% JoinLayers() %>% .[, .$annotation_full == "T"]
DefaultAssay(sc_t) <- "RNA"
sc_t[['sketch']] <- NULL

VlnPlot(sc_t, features = "nFeature_RNA", group.by = 'project', pt.size = 0)

sc_t[['RNA']] <- split(sc_t[['RNA']], sc_t$project, layers = "counts")


# cellcycle
sc_t <- 
  sc_t %>%
  CellCycleScoring(
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes
  )

# preprocess
sc_t <-
  sc_t %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 4000) %>%
  NormalizeData(scale.factor = median(sc_t$nCount_RNA)) %>%
  ScaleData(vars.to.regress = c('nCount_RNA')) # , 'S.Score', 'G2M.Score'))

# sc_t <- sc_t %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))

VariableFeaturePlot(sc_t)
# VariableFeatures(sc_t) <- VariableFeatures(sc_t) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))



## pca -----
sc_t <- sc_t %>% RunPCA()

pdf(file.path(t_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sc_t, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_t <- 
  sc_t[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_t) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(t_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_t[, suset_cells_t], dims = 1:5)
dev.off()
pdf(file.path(t_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_t[, suset_cells_t], dims = 6:10)
dev.off()
pdf(file.path(t_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_t[, suset_cells_t], dims = 11:15)
dev.off()
pdf(file.path(t_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_t[, suset_cells_t], dims = 16:20)
dev.off()


## dimension reduction --------
t_dims <-  20

 
sc_t <- RunUMAP(sc_t, dims = 1:t_dims)


pdf(file.path(t_out_dir, str_glue("Dimplot dim{t_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_t <- sc_t %>% 
  harmony::RunHarmony(
    dims.use = 1:t_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_t <- RunUMAP(
  sc_t, 
  dims = 1:t_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(t_out_dir, str_glue("Dimplot dim{t_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# plot QC figure
pdf(file.path(t_out_dir, str_glue("qc - Dimplot dim{t_dims}.pdf")), height = 6, width = 8)
FeaturePlot(sc_t, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI')) + NoLegend()
FeaturePlot(sc_t, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(t_out_dir, str_glue("qc - Vlnplot dim{t_dims}.pdf")), height = 6, width = 8)
VlnPlot(sc_t, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), group.by = 'project', pt.size = 0) + NoLegend()
dev.off()

FeaturePlot(sc_t, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
DimPlot(sc_t, group.by = c('Phase'), reduction = "harmonyumap")


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_t <-
  FindNeighbors(
    sc_t,
    reduction = reduction.use,
    dims = 1:t_dims,
    k.param = 10
  )

sc_t <-
  FindClusters(
    sc_t,
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(t_out_dir, str_glue("Dimplot dim{t_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_t, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(t_out_dir, str_glue("Dimplot dim{t_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_t, group.by = "group2")
# dev.off()

DotPlot(sc_t, group.by = "RNA_snn_res.0.4", features = c(
  "CD3D",
  "CD3E",
  "TRAC",
  "CD4",
  "CD8A",
  "MS4A1",
  "CD79A",
  "JCHAIN",
  "EPCAM",
  "PECAM1",
  "COL1A1",
  "RGS5",
  "TAGLN",
  "LYZ",
  'CD74',
  'CD14',
  'FCN1',
  'CD68',
  'C1QA',
  'KIT'
))



## remove mix cluster --------
sc_t$RNA_snn_res.0.4 %>% table
sc_t2 <- sc_t[, !sc_t$RNA_snn_res.0.4 %in% 9:13]

pdf(file.path(t_out_dir, str_glue("Dimplot dim{t_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_t2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


pdf(file.path(t_out_dir, str_glue("Dotplot dim{t_dims} cluster - {reduction.use} - remove small cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t2,
  features = major_marker_ls,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()


pdf(file.path(t_out_dir, str_glue("Dotplot dim{t_dims} cluster - {reduction.use} - cd4t - remove small cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t2,
  features = cd4t_marker,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

pdf(file.path(t_out_dir, str_glue("Dotplot dim{t_dims} cluster - {reduction.use} - cd8t - remove small cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t2,
  features = cd8t_marker,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

# sc_t3 subcluster --------------------
# 0: fibroblast, 
# 3: myeloid, myeloid cells和T细胞交较多，因此不剔除此群
# 5: B, 
# 7: CD3-NCAM+ -> NK cells

sc_t3 <- sc_t2[, !sc_t2$RNA_snn_res.0.4 %in% c(0, 5)]


sc_t3 <- sc_t3 %>% ScaleData(vars.to.regress = c('nCount_RNA')) # , 'S.Score', 'G2M.Score'))
sc_t3 <- RunPCA(sc_t3)

ElbowPlot(sc_t3, ndims = 30)


sc_t3_dims <- 20

sc_t3 <- sc_t3 %>% 
  harmony::RunHarmony(
    dims.use = 1:sc_t3_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_t3 <- RunUMAP(
  sc_t3, 
  dims = 1:sc_t3_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

sc_t3 <-
  FindNeighbors(
    sc_t3,
    reduction = reduction.use,
    dims = 1:sc_t3_dims,
    k.param = 10
  )

sc_t3 <-
  FindClusters(
    sc_t3,
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
  )




# CD8T marker
cd8t_marker <- c(
  "CD8_Tex" = "PDCD1",
  "CD8_Tex" = "HAVCR2",
  "CD8_Tex" = "CTLA4",
  "CD8_Tex" = "CXCL13",
  "CD8_Tex" = "LAYN",
  "CD8_Tex" = "LAG3",
  "CD8_Tex" = "TIGIT",
  "CD8_Tex" = "ENTPD1",
  
  "CD8_Tem" = "GZMK",
  "CD8_Tem" = "GZMA",
  "CD8_Tem" = "GZMB",
  "CD8_Tem" = "IFNG",
  "CD8_Temra" = "KLRG1",
  "CD8_Temra" = "PRF1",
  "CD8_Temra" = "NKG7",
  
  "CD8_Trm" = "MYADM",
  "CD8_Trm" = "CD69",
  "CD8_Trm" = "ITGAE",
  "CD8_Trm" = "CXCR6",
  "CD8_Trm" = "ITGA1",
  
  "CD8_MAIT" = "NCR3",
  "CD8_MAIT" = "SLC4A10",
  
  "CD8_Naive" = "SELL",
  "CD8_Naive" = "LEF1",
  "CD8_Naive" = "TCF7",
  "CD8_Naive" = "CCR7"
)

cd4t_marker <- c(
  "CD4T_Treg" = "TNFRSF4",
  "CD4T_Treg" = "TNFRS9",
  "CD4T_Treg" = "FOXP3",
  "CD4T_Treg" = "IL2RA",
  "CD4T_Treg" = "LAYN",
  "CD4T_Treg" = "CTLA4",
  "CD4T_Treg" = "TIGIT",
  "CD4T_Treg" = "RTKN2",
  
  "CD4T_Tem" = "GZMK",
  "CD4T_Tem" = "GZMA",
  "CD4T_Tem" = "NKG7",
  "CD4T_Tem" = "CCL4",
  "CD4T_Tem" = "CCL5",
  "CD4T_Tem" = "CAPG",
  
  "CD4T_Th17" = "IL17A",
  "CD4T_Th17" = "IL17F",
  
  "CD4_Th1" = "ANXA1",
  "CD4_Th1" = "CXCL13",
  "CD4_Th1" = "BHLHE40",
  "CD4_Th1" = "IFNG",
  
  "CD4_Tfh" = "TOX2",
  "CD4_Tfh" = "PDCD1",
  
  "CD4_Tcm" = "TCF7",
  "CD4_Tcm" = "CCR7",
  "CD4_Tcm" = "S1PR1",
  "CD4_Tcm" = "TXK"
)

# dimplot
pdf(file.path(t_out_dir, str_glue("Dimplot dim{sc_t3_dims} cluster - remove small and mix  cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_t3, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_t3, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t3, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(t_out_dir, str_glue("Dotplot dim{sc_t3_dims} cluster - {reduction.use} - remove small and mix  cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t3,
  features = major_marker_ls,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

pdf(file.path(t_out_dir, str_glue("Dotplot dim{sc_t3_dims} cluster - {reduction.use} - cd4t - remove small and mix cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t3,
  features = cd4t_marker,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

pdf(file.path(t_out_dir, str_glue("Dotplot dim{sc_t3_dims} cluster - {reduction.use} - cd8t - remove small and mix  cluster.pdf")), height = 6, width = 18)
DotPlot(
  sc_t3,
  features = cd8t_marker,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

# annotation ----
Idents(sc_t2) <- "RNA_snn_res.0.4"

cell_cluster <- c(
  '0' = "Fibroblast",
  '1' = "CD8_Tem",
  '2' = "CD4_Tcm",
  '3' = "CD8_Tem",
  '4' = "CD8_Temra",
  '5' = "B cells",
  '6' = "CD4T_Treg",
  '7' = "NK",
  '8' = "CD4T_Treg"
)

sc_t2 <- RenameIdents(sc_t2, cell_cluster)
sc_t2$annotation <- Idents(sc_t2)


pdf(file.path(t_out_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_t2, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_t2, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(t_out_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(sc_t2, group.by = "annotation", features = major_marker_ls)
dev.off()



# enrichment analysis ------
t_all_markers <- 
  sc_t2[, !as.character(sc_t2$annotation) %in% c('Fibroblast', 'B cells', 'NK')] %>%
  FindAllMarkers(
    logfc.threshold = 0.1,
    min.pct = 0.1,
    return.thresh = 1.1,
    only.pos = FALSE
  )

t_all_markers_sub <- 
  t_all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(
    avg_log2FC > 1,
    p_val_adj < 0.05
  )

t_all_markers_id <-
  t_all_markers %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map('gene') %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = 'org.Hs.eg.db',
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  })

all_enrich_gobp_ls <-
  t_all_markers_id %>%
  map(function(gene){
    clusterProfiler::enrichGO(
      gene,
      OrgDb = 'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.2
    )
  })


p_enrich_gobp_ls <-
  all_enrich_gobp_ls %>%
  imap(function(x, y){
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 60,
      title = y
    )
  })

pdf(
  file.path(t_out_dir, "t subcluster gobp enrichment analysis.pdf"),
  width  = 8,
  height = 6
)
p_enrich_gobp_ls %>% walk(print)
dev.off()

# multi volcano plot ---------

set.seed(1234)
rect_height <- 0.6
dodge_width <- 0.1


# scatter
df_scatter <-
  t_all_markers %>%
  mutate(cluster = factor(cluster), x_base = as.numeric(cluster)) %>%
  slice_sample(prop = 1) %>%
  group_by(cluster) %>%
  
  mutate(
    x_delta = 1:n(),
    x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
    xx = x_base - 1 + x_delta,
    y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
  ) %>%
  mutate(
    group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
    group = factor(group, levels = c("p > 0.01", "p < 0.01"))
  ) %>%
  dplyr::arrange(group)

# df_label 
df_label_up <- df_scatter %>%
  dplyr::filter(as.character(group) == "p < 0.01") %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10)

df_label_down <- df_scatter %>%
  dplyr::filter(as.character(group) == "p < 0.01") %>%
  group_by(cluster) %>% 
  slice_min(avg_log2FC, n = 10)

df_label_center <- df_scatter %>%
  dplyr::select(cluster, x_base) %>%
  distinct() %>%
  mutate(x_base2 = x_base - 0.5)
  


# grp rect  
df_grp_rect <-
  t_all_markers %>%
  group_by(cluster) %>%
  summarise() %>%
  mutate(x_base = as.numeric(factor(cluster))) %>%
  mutate(xmin = x_base - 1, xmax = x_base) %>%
  mutate(ymin = -rect_height, ymax = rect_height)


# background rect
df_bg_rect <-
  t_all_markers  %>%
  group_by(cluster) %>%
  summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
  mutate(x_base = as.numeric(factor(cluster))) %>%
  mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)



p_multi_volcano <-
  df_scatter %>%
  ggplot(aes(x = xx, y = y)) +
  geom_rect(
    data = df_bg_rect, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray90",
    inherit.aes = FALSE
  ) +
  geom_rect(
    data = df_grp_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
    inherit.aes = FALSE
  ) +
  geom_point(aes(color = group), size = 1) +
  
  ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
  # ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
  geom_text(data = df_label_center, aes(x = x_base2, y = 0, label = cluster), inherit.aes = FALSE) +
  scale_color_manual(values = c(
    'p < 0.01' = "red",
    'p > 0.01' = "gray70"
  )) +
  ggsci::scale_fill_igv() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Log2FC"
  )

# export
ggsave(
  p_multi_volcano,
  filename = file.path(t_out_dir, "t subcluster de analysis - multi volcano plot.pdf"),
  width  = 8,
  height = 6
)
