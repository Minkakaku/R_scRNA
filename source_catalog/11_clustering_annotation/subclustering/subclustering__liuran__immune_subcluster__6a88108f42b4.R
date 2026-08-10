
library(Seurat)
library(tidyverse)
library(patchwork)

# Immune subcluster ---------
immune_out_dir <- "Immune_subcluster" # dir.create(immune_out_dir)
dir.create(immune_out_dir)

# sc_immune <- all_data[, all_data$annotation_full %in% c("Club/Basel", "Epithelial", "Luminal", "NE", "Macrophage")]
sc_immune <- all_data[, all_data$annotation_full %in% c("B and Plasma cells", "Endothelial", "Fibroblasts", "Mast cells", "Pericytes", "Smooth Muscle Cell", "T cells")]
DefaultAssay(sc_immune) <- "RNA"
sc_immune[['sketch']] <- NULL


sc_immune <-
  sc_immune %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_immune <- sc_immune %>% RunPCA()

pdf(file.path(immune_out_dir, 'pca plot-immune.pdf'), height = 6, width = 8)
ElbowPlot(sc_immune, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_immune <- 
  sc_immune[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_immune) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(immune_out_dir, 'immune dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_immune[, suset_cells_immune], dims = 1:5)
dev.off()
pdf(file.path(immune_out_dir, 'immune dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_immune[, suset_cells_immune], dims = 6:10)
dev.off()
pdf(file.path(immune_out_dir, 'immune dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_immune[, suset_cells_immune], dims = 11:15)
dev.off()
pdf(file.path(immune_out_dir, 'immune dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_immune[, suset_cells_immune], dims = 16:20)
dev.off()


## dimension reduction --------
immune_dims <- 10 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_immune <- RunUMAP(sc_immune, dims = 1:immune_dims)


pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_immune, group.by = "project", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} project - split.pdf")), height = 8, width = 16)
DimPlot(sc_immune, group.by = "project", ncol = 3, split.by = "project") + NoLegend()
dev.off()


sc_immune <- sc_immune %>% 
  harmony::RunHarmony(
    dims.use = 1:immune_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_immune <- RunUMAP(
  sc_immune, 
  dims = 1:immune_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)



pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_immune, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} project - harmony - split.pdf")), height = 8, width = 16)
DimPlot(sc_immune, group.by = "project", reduction = "harmonyumap", ncol = 3, split.by = "project") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_immune <-
  FindNeighbors(
    sc_immune,
    reduction = reduction.use,
    dims = 1:immune_dims,
    k.param = 10
  )

sc_immune <- 
  FindClusters(
    sc_immune,
    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_immune, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.7", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_immune, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.7", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_immune, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_immune, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_immune) <- "RNA_snn_res.0.8"
DefaultAssay(sc_immune) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(sc_immune),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
# all_markers <-
#   FindAllConservedMarkers(
#     JoinLayers(sc_immune),
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
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_immune))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(immune_out_dir, "all_marker_table.xlsx"),
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
all_top5_gene_unique <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc_immune, 
    features = all_top5_gene_unique
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(immune_out_dir, "top 5 gene dotplot.pdf"),
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
    filename = "Immune_subcluster/all_markers_heatmap.pdf", width = 20, height = 5,
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )




avg_expr[rownames(avg_expr) %in% immune_marker_ls2, ] %>% 
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

# saveRDS(sc_immune, file = file.path(immune_out_dir, "sc_immune.rds"))
## remove 6 7 8--------
# sc_immune2 <- sc_immune[, !sc_immune$RNA_snn_res.0.2 %in% 6:8]

# pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
# DimPlot(sc_immune, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()


# dotplot
pdf(file.path(immune_out_dir, str_glue("dotplot dims{immune_dims} {reduction.use} major_marker_ls.pdf")), width = 20, height = 6)
DotPlot(sc_immune, group.by = "RNA_snn_res.0.1", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.5", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.7", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(sc_immune, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


# annotation -------
Idents(sc_immune) <- "RNA_snn_res.0.8"
cell_cluster <- c(
  '0'  = "SMC",
  '1'  = "Fibroblast",
  '2'  = "Fibroblast",
  '3'  = "Fibroblast",
  '4'  = "Fibroblast",
  '5'  = "Fibroblast",
  '6'  = "T cell",
  '7'  = "Pericyte",
  '8'  = "Fibroblast",
  '9' = "Endothelial",
  '10'  = "T cell",
  '11'  = "SMC",
  '12'  = "Endothelial",
  '13'  = "B and Plasma",
  '14'  = "Fibroblast",
  '15'  = "Fibroblast",
  '16'  = "Mast",
  '17'  = "Mast", 
  '18'  = "T cell", 
  '19'  = "Fibroblast", 
  '20'  = "Fibroblast", 
  '21'  = "Fibroblast", 
  '22'  = "Fibroblast", 
  '23'  = "T cell"
)

sc_immune <- RenameIdents(sc_immune, cell_cluster)
sc_immune$annotation_sub2 <- Idents(sc_immune)


pdf(file.path(immune_out_dir, str_glue("immune Dimplot dim{immune_dims} annotation_sub2 {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_immune, group.by = "annotation_sub2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_immune, group.by = "annotation_sub2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(sc_immune) <- "RNA"
Idents(sc_immune) <- "annotation_sub2"

pdf(file.path(immune_out_dir, str_glue("immune dotplot dims{immune_dims} {reduction.use} annotation_sub2 major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_immune, group.by = "annotation_sub2", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# pdf(file.path(immune_out_dir, str_glue("immune dotplot dims{immune_dims} {reduction.use} annotation_sub2 immune_marker_ls1.pdf")), width = 15, height = 6)
# DotPlot(sc_immune, group.by = "annotation_sub2", features = immune_marker_ls1, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()
# 
# pdf(file.path(immune_out_dir, str_glue("immune dotplot dims{immune_dims} {reduction.use} annotation_sub2 immune_marker_ls2.pdf")), width = 8, height = 6)
# DotPlot(sc_immune, group.by = "annotation_sub2", features = immune_marker_ls2, dot.scale = 10) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
# dev.off()



sample_map_group <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "")
sc_immune@meta.data$project <- factor(sc_immune@meta.data$project, levels = sample_map_group[,1])
sc_immune@meta.data$group <- factor(sc_immune@meta.data$group, levels = unique(sample_map_group[,2]))


# cell props 
sc_immune_cell_props_df <- 
  sc_immune[[]] %>% 
  dplyr::select(project, group, annotation_sub2) %>%
  group_by(project, group, annotation_sub2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_immune_props_project <-
  sc_immune_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = annotation_sub2)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_igv()

p_sc_immune_props_group <-
  sc_immune_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = annotation_sub2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_igv()

pdf(file.path(immune_out_dir, str_glue("immune - cell props baplot.pdf")), height = 4, width = 8)
p_sc_immune_props_project %>% print()
p_sc_immune_props_group   %>% print()
dev.off()

saveRDS(sc_immune, file = file.path(immune_out_dir, "sc_immune.rds"))









## 替换 immune_subcluster 注释 ----------
out_gene_umi = as.matrix(sc_immune@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
write.table(out_gene_umi, file = "Immune_subcluster/Filter_Basicinfo_nGene_nUMI_mito.txt", sep = "\t", quote = F, row.names = F, col.names = F)

annotation_sub2_metadata <- read.table("Immune_subcluster/Filter_Basicinfo_nGene_nUMI_mito.txt", header = T, sep = "\t", check.names = F)
colnames(annotation_sub2_metadata)
annotation_sub2_metadata <- annotation_sub2_metadata[, c(1,30)]
colnames(annotation_sub2_metadata) <- c("barcode", "annotation_sub2")


all_data@meta.data$Cells <- rownames(all_data@meta.data)
all_data@meta.data$barcode <- all_data@meta.data$Cells
all_data@meta.data <- all_data@meta.data %>%
  left_join(annotation_sub2_metadata, by = "barcode") %>%
  mutate(
    annotation_updated = case_when(
      # 在亚群对象中有细化的注释
      !is.na(annotation_sub2) ~ annotation_sub2,
      
      # 保持非T细胞类型的原始注释
      !annotation_sub %in% 
        c("B and Plasma cells", "Endothelial", "Fibroblasts", "Mast cells", "Pericytes", "Smooth Muscle Cell", "T cells") ~ annotation_sub,
      
      # 对于未细化的T细胞给予特殊标记
      TRUE ~ paste0(annotation_sub, " (unclassified)")
    )
  ) %>%
  select(-c(barcode, annotation_sub2))  # 移除临时列
unique(all_data$annotation_updated)
all_data@meta.data$annotation_sub <- all_data@meta.data$annotation_updated

rownames(all_data@meta.data) <- all_data@meta.data$Cells


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation_sub {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_sub", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation_sub", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation_sub {reduction.use} split full.harmonyumap.pdf")), height = 6, width = 14)
DimPlot(all_data, group.by = "annotation_sub", label = TRUE, reduction = 'umap', ncol = 2, split.by = "group") + NoLegend()
DimPlot(all_data, group.by = "annotation_sub", label = TRUE, reduction = 'full.harmonyumap', ncol = 2, split.by = "group") + NoLegend()
dev.off()








# spatial plot --------------
spatial_outdir <- "Immune_subcluster/Spatial plot"
dir.create(spatial_outdir)

all_data_meta <-
  all_data[[]] %>%
  filter(annotation_sub != "Other") %>% 
  rownames_to_column("cell_id")

plan <- c("#00B4F0", "#7CAE00", "#FF64B0", "#619CFF", "#F8766D", "#00BA38", "#00C08B", "#DE8C00", "#C77CFF", "#005995")
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
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_sub, group = cell_id))  +
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
      file.path(spatial_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })







# 写出新的metadata
out_gene_umi = as.matrix(all_data@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
out_gene_umi %>% colnames()
# out_gene_umi <- out_gene_umi[, -c(2:8, 12:19, 20:23, 25)]
write.table(out_gene_umi, file = "all_data_annotation_sub_metadata.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.csv(out_gene_umi, file = "all_data_annotation_sub_metadata.csv", row.names = F, col.names = F, quote = F)


# 写出csv导入loupe browser
meta_data <- read.table("all_data_annotation_sub_metadata.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data)
meta_data$Cell <- gsub("^[0-9]+_", "", meta_data$Cell)
meta_data$sample_orig <- paste(meta_data$project, meta_data$orig.ident, sep = "_")
# "Pre_2"  "Post_2" "Pre_3"  "Post_3" "Pre_1"  "Post_1"


Pre_1_meta_data <- meta_data[meta_data$project == "Pre_1", ]
Pre_1_meta_data <- Pre_1_meta_data[, c(1, 26)]
colnames(Pre_1_meta_data) <- c("cell_id", "group")

Pre_2_meta_data <- meta_data[meta_data$project == "Pre_2", ]
Pre_2_meta_data <- Pre_2_meta_data[, c(1, 26)]
colnames(Pre_2_meta_data) <- c("cell_id", "group")

Pre_3_meta_data <- meta_data[meta_data$project == "Pre_3", ]
Pre_3_meta_data <- Pre_3_meta_data[, c(1, 26)]
colnames(Pre_3_meta_data) <- c("cell_id", "group")

Post_1_meta_data <- meta_data[meta_data$project == "Post_1", ]
Post_1_meta_data <- Post_1_meta_data[, c(1, 26)]
colnames(Post_1_meta_data) <- c("cell_id", "group")

Post_2_meta_data <- meta_data[meta_data$project == "Post_2", ]
Post_2_meta_data <- Post_2_meta_data[, c(1, 26)]
colnames(Post_2_meta_data) <- c("cell_id", "group")

Post_3_meta_data <- meta_data[meta_data$project == "Post_3", ]
Post_3_meta_data <- Post_3_meta_data[, c(1, 26)]
colnames(Post_3_meta_data) <- c("cell_id", "group")


write.table(Pre_1_meta_data, file = "Pre_1_849384_12_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(Pre_2_meta_data, file = "Pre_2_1140629_15_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(Pre_3_meta_data, file = "Pre_3_1206894_4_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(Post_1_meta_data, file = "Post_1_883994_34_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(Post_2_meta_data, file = "Post_2_1160884_28_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(Post_3_meta_data, file = "Post_3_1227144_35_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)











