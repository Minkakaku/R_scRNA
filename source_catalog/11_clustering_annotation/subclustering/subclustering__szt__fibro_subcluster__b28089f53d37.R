library(tidyverse)
library(Seurat)
library(patchwork)

# list.files("/data/xenium/20241010__051140__FWSP20241029")

# # library --------------
# options(java.parameters = '-Xmx10g') # options('java.parameters')
# library(rJava)
# library(RBioFormats)
# source("xenium.R")

# multi workers  -----
library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 10 * 1024 ^ 3)

# nbrOfFreeWorkers()


# Find all markers using FindConservedMarkers
FindAllConservedMarkers <-
  function(obj, 
           ident.1 = NULL, 
           grouping.var,
           ...) {
    
  stopifnot(is(obj, "Seurat"))
  
  all_ident <- unique(Idents(obj))
  stopifnot(length(all_ident) > 1)
  
  
  if(is.null(ident.1)){
    message("ident.1 is null, will find all markers of each ident in all idents.")
  }else{
    all_ident <-  ident.1
  }
  
  all_ident %>%
    set_names(., .) %>%
    map(function(ident){
      FindConservedMarkers(
        object = obj,
        ident.1 = ident,
        grouping.var = grouping.var,
        ...
      )
    })
}


# Fibroblast subcluster ---------
fb_out_dir <- "fibroblast_subcluster" # dir.create(fb_out_dir)

all_data <- readRDS("output/data/adata_filtered_add_celltype.rds")

sc_fibro <- all_data[, all_data$annotation_full == "Fibroblast"]
DefaultAssay(sc_fibro) <- "RNA"
sc_fibro[['sketch']] <- NULL


sc_fibro <-
  sc_fibro %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_fibro <- sc_fibro %>% RunPCA()

pdf(file.path(fb_out_dir, 'pca plot-fibro.pdf'), height = 6, width = 8)
ElbowPlot(sc_fibro, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_fibro <- 
  sc_fibro[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_fibro) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(fb_out_dir, 'fibro dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 1:5)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 6:10)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 11:15)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 16:20)
dev.off()


## dimension reduction --------
fibro_dims <- 10 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_fibro <- RunUMAP(sc_fibro, dims = 1:fibro_dims)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_fibro <- sc_fibro %>% 
  harmony::RunHarmony(
    dims.use = 1:fibro_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_fibro <- RunUMAP(
  sc_fibro, 
  dims = 1:fibro_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'pca' # 'harmony' 'pca'

## cluster ------
sc_fibro <-
  FindNeighbors(
    sc_fibro,
    reduction = reduction.use,
    dims = 1:fibro_dims,
    k.param = 10
  )

sc_fibro <-
  FindClusters(
    sc_fibro,
    resolution = c(0.05, 0.1, 0.2) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_fibro, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_fibro) <- "RNA_snn_res.0.2"
DefaultAssay(sc_fibro) <- "RNA"

# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_fibro),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
all_markers <-
  FindAllConservedMarkers(
    JoinLayers(sc_fibro),
    grouping.var = "project",
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
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_fibro))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(fb_out_dir, "all_marker_table.xlsx"),
    asTable = TRUE
  )

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

avg_expr[all_de_genes %>% unique, -(7:9)] %>% 
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

fibro_genes <- c(
  'VIM', 'THY', 'S100A4', 'DCN',
  'COL1A1', 'ACTA2', 'PDPN', 'FAP'
)

avg_expr[rownames(avg_expr) %in% fibro_genes, ] %>% 
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

# saveRDS(sc_fibro, file = file.path(fb_out_dir, "sc_fibro.rds"))
## remove 6 7 8--------
sc_fibro2 <- sc_fibro[, !sc_fibro$RNA_snn_res.0.2 %in% 6:8]

pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_fibro2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


# heatmap 7-9是小群，剔除
pdf(file.path(fb_out_dir, str_glue("fibro marker heatmapdim dim{fibro_dims} - {reduction.use} - trimmed.pdf")), height = 4, width = 8)
avg_expr[all_de_genes %>% unique %>% setdiff(c('CCDC80', 'TCF21')), -(7:9)] %>% 
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
dev.off()

# cell props 
sc_fibro_cell_props_df <-
  sc_fibro2[[]] %>% 
  dplyr::select(project, group, RNA_snn_res.0.2) %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_fibro_props_project <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_bmj()

p_sc_fibro_props_group <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_bmj()

pdf(file.path(fb_out_dir, str_glue("fibro - cell props baplot - trimmed.pdf")), height = 4, width = 8)
p_sc_fibro_props_project %>% print()
p_sc_fibro_props_group   %>% print()
dev.off()


## dotplot ----------
FeaturePlot(
  sc_fibro,
  features = c(
    "ACTA2", "S100A4", "COL1A1", "VIM", 
    "PDGFRA", "PDGFRB", "THY1",
    "DCN", "FAP", "PDPN"
  ),
  reduction = "harmonyumap"
)


DotPlot(
  sc_fibro,
  features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
  group.by = "RNA_snn_res.0.6"
)


DimPlot(sc_fibro, reduction = "harmonyumap", label = TRUE) + NoLegend() 


