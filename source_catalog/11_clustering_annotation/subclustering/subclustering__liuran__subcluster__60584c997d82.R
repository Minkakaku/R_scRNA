
library(tidyverse)
library(Seurat)
library(patchwork)


all_data <- readRDS("data_for_analysis/all_data_h5ad_ROI.rds")

all_data_filt <- all_data[, !all_data@meta.data$sample %in% c("B1", "B4")]

# microgliablast subcluster ---------
Micro_out_dir <- "Mircoglia_subcluster" # dir.create(Micro_out_dir)
dir.create(Micro_out_dir)

sc_micro <- all_data_filt[, all_data_filt$annotation == "Microglial"]
sc_micro_filt <- subset(sc_micro, subset = Tmem119 >0)

saveRDS(sc_micro, "Mircoglia_subcluster/sc_micro.rds")


DefaultAssay(sc_micro_filt) <- "RNA"
# sc_micro_filt[['sketch']] <- NULL


sc_micro_filt <-
  sc_micro_filt %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_micro_filt <- sc_micro_filt %>% RunPCA()

pdf(file.path(Micro_out_dir, 'pca plot-microglia.pdf'), height = 6, width = 8)
ElbowPlot(sc_micro_filt, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_microglia <- 
  sc_micro_filt[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_micro_filt) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(Micro_out_dir, 'microglia dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_micro_filt[, suset_cells_microglia], dims = 1:5)
dev.off()
pdf(file.path(Micro_out_dir, 'microglia dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_micro_filt[, suset_cells_microglia], dims = 6:10)
dev.off()
pdf(file.path(Micro_out_dir, 'microglia dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_micro_filt[, suset_cells_microglia], dims = 11:15)
dev.off()
pdf(file.path(Micro_out_dir, 'microglia dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_micro_filt[, suset_cells_microglia], dims = 16:20)
dev.off()


## dimension reduction --------
microglia_dims <- 12 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_micro_filt <- RunUMAP(sc_micro_filt, dims = 1:microglia_dims)


pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_micro_filt, group.by = "project", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} project - split.pdf")), height = 16, width = 16)
DimPlot(sc_micro_filt, group.by = "project", ncol = 4, split.by = "project") + NoLegend()
dev.off()


sc_micro_filt <- sc_micro_filt %>% 
  harmony::RunHarmony(
    dims.use = 1:microglia_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_micro_filt <- RunUMAP(
  sc_micro_filt, 
  dims = 1:microglia_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_micro_filt, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} project - harmony - split.pdf")), height = 16, width = 16)
DimPlot(sc_micro_filt, group.by = "project", reduction = "harmonyumap", ncol = 4, split.by = "project") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_micro_filt <-
  FindNeighbors(
    sc_micro_filt,
    reduction = reduction.use,
    dims = 1:microglia_dims,
    k.param = 10
  )

sc_micro_filt <- 
  FindClusters(
    sc_micro_filt,
    resolution = c(0.2, 0.3, 0.4, 0.5, 0.6) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_micro_filt, group.by = "group2")
# dev.off()

major_marker_ls <- c(
  "B cell" = "Cd19",
  "B cell" = "Cd79a", 
  "B cell" = "Ms4a1",
  "Plasma" = "Mzb1",
  "Plasma" = "Jchain",
  "T cell" = "Cd3e", 
  "T cell" = "Cd4",
  "T cell" = "Foxp3", 
  "T cell" = "Cd8a", 
  "T cell" = "Gzmb",
  "NK" = "Ncam1",
  "Endothelial" = "Pecam1",
  "Endothelial" = "Vwf",
  "Epithelial" = "Epcam",
  # "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Col1a1",
  "Oligo" = "Mog",
  "Oligo" = "Sox10",
  "Oligo" = "Olig1",
  "OPC" = "Pdgfra",
  "OPC" = "Cspg4",
  "Neuron" = "Tubb3",
  "Neuron" = "Map2",
  "Neuron" = "Gap43",
  "Neuron" = "Snap25",
  "Astrocytes" = "Aqp4",
  "Astrocytes" = "Agt",
  "Astrocytes" = "Rorb",
  "Astrocytes" = "Slc7a10",
  "Astrocytes" = "Gja1",
  "Microglial" = "Tmem119",
  "Microglial" = "Dock8",
  "Macrophage" = "Cd68", 
  "Macrophage" = "Cd163",
  "Dendritic" = "Cd209a", 
  "Dendritic" = "Cd83", 
  "Dendritic" = "Clec9a", 
  "Dendritic" = "Ido1",
  "Neutrophil" = "Csf3r",
  "Epen" = "Foxj1",
  "Epen" = "Dnah11",
  "Glutamatergic neurons" = "Slc17a7",
  "GABAergic neurons" = "Gad1",
  "GABAergic neurons" = "Gad2",
  "MSN" = "Syt6",
  "MSN" = "Otof",
  "MSN" = "Drd1",
  "MSN" = "Adora2a",
  "Dopaminergic neurons" = "Th",
  "Serotonergic neurons" = "Tph1",
  "Serotonergic neurons" = "Tph2",
  "Cholinergic neurons" = "Pdgfrb",
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Myh11",
  "Pericytes" = "Kcnj8",
  "Pericytes" = "Abcc9",
  "Radial glial" = "Nes",
  "Radial glial" = "Cdh2",
  "Radial glial" = "Fabp7"
)

# dotplot
pdf(file.path(Micro_out_dir, str_glue("dotplot dims{microglia_dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.5", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()



pdf(file.path(Micro_out_dir, str_glue("FeaturePlot all_data umap major_marker.pdf")), width = 8, height = 6)
FeaturePlot(sc_micro_filt, features = "Tmem119", reduction = "umap")
dev.off()

pdf(file.path(Micro_out_dir, str_glue("FeaturePlot all_data {reduction.use} major_marker.pdf")), width = 8, height = 6)
FeaturePlot(sc_micro_filt, features = "Tmem119", reduction = "harmonyumap")
dev.off()



## find markers -------
Idents(sc_micro_filt) <- "RNA_snn_res.0.4"
DefaultAssay(sc_micro_filt) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(sc_micro_filt),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
# all_markers <-
#   FindAllConservedMarkers(
#     JoinLayers(sc_micro_filt),
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
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_micro_filt))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(Micro_out_dir, "all_marker_table.xlsx"),
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
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]


# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc_micro_filt, 
    features = clean_top5
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(Micro_out_dir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 12
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
    filename = "Mircoglia_subcluster/all_markers_heatmap.pdf", width = 20, height = 5,
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
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




# saveRDS(sc_micro_filt, file = file.path(Micro_out_dir, "sc_micro_filt.rds"))
## remove 6 7 8--------
# sc_micro_filt2 <- sc_micro_filt[, !sc_micro_filt$RNA_snn_res.0.2 %in% 6:8]

# pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
# DimPlot(sc_micro_filt, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()


# dotplot
pdf(file.path(Micro_out_dir, str_glue("dotplot dims{microglia_dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.5", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_micro_filt, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


# sample_map_group <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "")
# sc_micro_filt@meta.data$project <- factor(sc_micro_filt@meta.data$project, levels = sample_map_group[,1])
# sc_micro_filt@meta.data$group <- factor(sc_micro_filt@meta.data$group, levels = unique(sample_map_group[,2]))

# cell props 
sc_micro_filt_cell_props_df <-
  sc_micro_filt[[]] %>% 
  dplyr::select(project, group, RNA_snn_res.0.4) %>%
  group_by(project, group, RNA_snn_res.0.4) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_micro_filt_props_project <-
  sc_micro_filt_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_igv()

p_sc_micro_filt_props_group <-
  sc_micro_filt_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.4)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_igv()

pdf(file.path(Micro_out_dir, str_glue("microglia - cell props baplot.pdf")), height = 5, width = 8)
p_sc_micro_filt_props_project %>% print()
p_sc_micro_filt_props_group   %>% print()
dev.off()

saveRDS(sc_micro_filt, file = file.path(Micro_out_dir, "sc_micro_filt.rds"))

## dotplot ----------
# FeaturePlot(
#   sc_micro_filt,
#   features = str_to_title(microglia_marker_ls1),
#   reduction = "harmonyumap"
# )


# DotPlot(
#   sc_micro_filt,
#   features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
#   group.by = "RNA_snn_res.0.6"
# )


# 剔除 cluster 1,4,5,6,10
# sc_micro_filt_1 <- sc_micro_filt[, !sc_micro_filt@meta.data$RNA_snn_res.0.4 %in% c("2", "3", "4", "6", "7", "8", "9", "10", "12")]


## annotation ---- res.0.4 细胞注释
Idents(sc_micro_filt) <- "RNA_snn_res.0.4"
cell_cluster <- c(
  '0'  = "Micro_0",
  '1'  = "Micro_1",
  '2'  = "Micro_2",
  '3'  = "Micro_3",
  '4'  = "Micro_4",
  '5'  = "Micro_5",
  '6'  = "Micro_6",
  '7'  = "Micro_7",
  '8'  = "Micro_8",
  '9'  = "Micro_9",
  '10' = "Micro_10",
  '11' = "Micro_11"
)

cell_cluster <- c(
  '0'  = "Neutrophil_like_Microglia",
  '1'  = "Neuron_like_Microglia",
  '2'  = "Astrocyte_like_Microglia",
  '3'  = "Neuron_like_Microglia",
  '4'  = "DC_like_Microglia",
  '5'  = "Micro_5",
  '6'  = "Oligo_like_Microglia",
  '7'  = "Neuron_like_Microglia",
  '8'  = "OPC_like_Microglia",
  '9'  = "Neuron_like_Microglia",
  '10' = "Micro_10",
  '11' = "Glutamatergic_like_Microglia"
)

sc_micro_filt <- RenameIdents(sc_micro_filt, cell_cluster)
sc_micro_filt$annotation_micro <- Idents(sc_micro_filt)


pdf(file.path(Micro_out_dir, str_glue("microglia Dimplot dim{microglia_dims} annotation_micro {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_micro_filt, group.by = "annotation_micro", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_micro_filt, group.by = "annotation_micro", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(sc_micro_filt) <- "RNA"
Idents(sc_micro_filt) <- "annotation_micro"

pdf(file.path(Micro_out_dir, str_glue("microglia dotplot dims{microglia_dims} {reduction.use} annotation_micro major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_micro_filt, group.by = "annotation_micro", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# pdf(file.path(Micro_out_dir, str_glue("microglia dotplot dims{microglia_dims} {reduction.use} annotation_micro microglia_marker_ls1.pdf")), width = 15, height = 6)
# DotPlot(sc_micro_filt, group.by = "annotation_micro", features = microglia_marker_ls1, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()
# 
# pdf(file.path(Micro_out_dir, str_glue("microglia dotplot dims{microglia_dims} {reduction.use} annotation_micro microglia_marker_ls2.pdf")), width = 8, height = 6)
# DotPlot(sc_micro_filt, group.by = "annotation_micro", features = microglia_marker_ls2, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()


### 保存data
saveRDS(sc_micro_filt, file = file.path(Micro_out_dir, "sc_micro_filt.rds"))







## cell proportion ------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


sc_micro_filt <- readRDS(file.path(Micro_out_dir, "sc_micro_filt.rds"))

micro_ratio <- file.path(Micro_out_dir, "Ratio")
dir.create(micro_ratio)

out_gene_umi = as.matrix(sc_micro_filt@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
write.table(out_gene_umi, file = file.path(micro_ratio, "meta_data.txt"), sep = "\t", quote = F, row.names = F, col.names = F)


data <- read.table(file.path(micro_ratio, "meta_data.txt"), header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1, 8, 6, 47)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

# data <- data %>% filter(Annotation != "Unknown")


mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file = file.path(micro_ratio, "Cell_type_stat.txt"), quote=F, row.names=T, col.names=F, sep="\t")

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
write.table(DATA,file = file.path(micro_ratio, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("microgliablast_subcluster/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="microgliablast_subcluster/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file=file.path(micro_ratio, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file=file.path(micro_ratio, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# # 限制样本顺序和分组顺序
# sample_map <- read.table("sample_map_group.txt",header=F,sep="\t",quote="",check.names=F)
# data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
# data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

data4plot1$Celltype <- factor(data4plot1$Celltype, levels = paste0("Micro_", 0:11))
data4plot2$Celltype <- factor(data4plot2$Celltype, levels = paste0("Micro_", 0:11))

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
  plot_layout(widths = c(3, 1)) + # 两张图各自的宽度
  plot_annotation(
    title = "Cell Type Ratio",
    theme = theme(plot.title = element_text(
      size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
    ))
  )
ggsave(p,file=file.path(micro_ratio, "ratio_barplot.pdf"),width=12, height=6)
ggsave(p,file=file.path(micro_ratio, "ratio_barplot.png"),width=12, height=6,dpi=300)





# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = file.path(micro_ratio, "Cell_type_stat_percent.txt"), header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
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
ggsave(p,file="microgliablast_subcluster/Ratio/diff_boxplot.pdf",width=10,height=5)
ggsave(p,file="microgliablast_subcluster/Ratio/diff_boxplot.png",width=10,height=5)














