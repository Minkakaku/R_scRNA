library(tidyverse)
library(Seurat)
library(scMetabolism)

out_dir <- "20250724 metabolism"
dir.create(out_dir)


all_data <- readRDS('output/data/adata_celltype.rds')
all_data_raw <- readRDS('output/data/adata_raw.rds')

# cell coords ---------
cell_coords <-
  all_data_raw@reductions$spatial %>% 
  Embeddings() %>% 
  as.data.frame() %>%
  dplyr::rename(
    x = SPATIAL_1,
    y = SPATIAL_2
  )


########################
# # norm_data to raw_data
# scale_factor <- all_data[['RNA']]@data[, 1:20] %>% expm1 %>%  colSums()
# stopifnot(all(scale_factor[1] == scale_factor))
# 
# scale_factor_1 <- round(scale_factor[1]) ## round
# stopifnot(scale_factor_1 ==  median(all_data$total_counts))
# 
# scale_factor_series <- all_data$total_counts / scale_factor_1
# 
# raw_data <- sweep(all_data[['RNA']]$data, MARGIN = 2, STATS = scale_factor_series, FUN = "*")


#######################


grp_info <- readxl::read_excel("annotation.xlsx")

grp_info_full <- 
  all_data[['sample']] %>%
  rownames_to_column("id") %>%
  left_join(grp_info, by = c('sample')) %>%
  column_to_rownames("id")


all_data <- 
  all_data %>%
  AddMetaData(grp_info_full)

# metabolisam pathway --------------
panel_gene <- rownames(all_data)

kegg_metabolisam_df <- clusterProfiler::read.gmt(system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism"))
# system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")
kegg_metabolisam_df_tidy <- 
  kegg_metabolisam_df %>%
  dplyr::filter(gene %in% panel_gene)  

kegg_metabolisam_vector <- 
  kegg_metabolisam_df %>% 
  group_by(gene) %>%
  nest() %>%
  deframe() %>%
  map_chr(function(x){
    str_c(unlist(x), collapse = ";")
  })
  
# scMetabolism --------------
all_data_raw[['RNA']] <- all_data_raw[['RNA']]  %>% as('Assay5') %>% split(f = all_data_raw$sample) 

all_data_raw <- NormalizeData(all_data_raw, scale.factor = floor(median(all_data_raw$total_counts)))
all_data_raw <- FindVariableFeatures(all_data_raw)


## sketch data and scMetabolisam --------
all_data_raw <-
  SketchData(
    all_data_raw,
    ncells         = 5000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


all_data_raw2 <- CreateSeuratObject(
  all_data_raw[['sketch']],
  meta.data = all_data_raw[[]]
)
all_data_raw2 <- JoinLayers(all_data_raw2)
all_data_raw2[['RNA']] <- as(all_data_raw2[['RNA']], 'Assay')

all_data_raw2 <- 
  sc.metabolism.Seurat(
    all_data_raw2,
    method = "VISION",
    ncores = 6,
    metabolism.type = "KEGG"
  )

# saveRDS(all_data_raw2, "all_data_raw.metabolisam.rds")
all_data_raw2 <- readRDS("all_data_raw.metabolisam.rds")



## tidy operation of seurat object --------------
all_data_raw2_meta <- 
  all_data_raw2@meta.data %>%
  rownames_to_column("id") %>% 
  left_join(
    rownames_to_column(all_data[[c('celltype', 'group2', 'group3', 'paired sample', 'info')]], "id"),
    by = "id"
  ) %>%
  column_to_rownames("id")
stopifnot(all(rownames(all_data_raw2_meta) == colnames(all_data_raw2)))
all_data_raw2@meta.data <- all_data_raw2_meta

colnames(all_data_raw2@assays$METABOLISM$score) <- colnames(all_data_raw2)

# only used cells which have celltype annotation 
inter_cells <- intersect(
  rownames(all_data_raw2_meta),
  colnames(all_data)
)
df_score <- all_data_raw2@assays$METABOLISM$score[, inter_cells]
meta_score <- all_data_raw2[[]][inter_cells, ]

all_data_raw2_tidy <- all_data_raw2[, inter_cells]
all_data_raw2_tidy@reductions$umap <- all_data[['umap']] %>% subset(cells = inter_cells)
all_data_raw2_tidy@meta.data <- cbind(all_data_raw2_tidy@meta.data, t(df_score))


# saveRDS(all_data_raw2_tidy, "all_data_raw.metabolisam2.rds")
all_data_raw2_tidy <- readRDS("all_data_raw.metabolisam2.rds")


## plot metabolisam  ------------
### heatmap by group -----
heatmap_grp  <-
  apply(df_score, 1, function(x)(
    tapply(x, meta_score$group, mean)
  )) %>%
  scale %>%
  t %>%
  ComplexHeatmap::Heatmap(
    .,
    name = "score",
    row_names_max_width = unit(12, "cm"),
    row_names_gp = grid::gpar(fontsize = 10),
    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",
    # top_annotation = ComplexHeatmap::HeatmapAnnotation(
    #   group = colnames(.)
    # )
  )

# export
pdf(
  file = file.path(out_dir, 'metabolism - heatmap - group.pdf'),
  height = 12,
  width =  6
)
print(heatmap_grp)
dev.off()

### heatmap by celltype ----
heatmap_ct <-
  apply(df_score, 1, function(x)(
    tapply(x, meta_score$celltype, mean)
  )) %>%
  scale %>%
  t %>%
  ComplexHeatmap::Heatmap(
    .,
    name = "score",    
    row_names_gp = grid::gpar(fontsize = 10),
    row_names_max_width = unit(12, "cm"),
    column_names_max_height = unit(6, "cm"),
    height = unit(20, "cm"),
    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",
    # top_annotation = ComplexHeatmap::HeatmapAnnotation(
    #   group = colnames(.)
    # )
  )

# export
pdf(
  file = file.path(out_dir, 'metabolism - heatmap - celltype.pdf'),
  height = 12,
  width =  8
)
print(heatmap_ct)
dev.off()

### dimplot  ---------
metabolism_pathways <- rownames(df_score)
p_dim_ct <- 
  DimPlot(
    all_data_raw2_tidy,
    group.by = 'celltype'
  ) +
  ggsci::scale_color_igv()

p_pathway_ls <-
  metabolism_pathways %>%
  set_names(., .) %>%
  map(function(x){
    FeaturePlot(
      all_data_raw2_tidy,
      feature = x,
      order = TRUE
    )
  })
p_pathway_ls <- c(list(Cell = p_dim_ct), p_pathway_ls)

# export
pdf(
  file = file.path(out_dir, 'metabolism - feature plot.pdf'),
  height = 6,
  width =  8
)
p_pathway_ls %>% walk(print)
dev.off()



# dotplot by group
DotPlot.metabolism(
  obj = all_data_raw2, 
  pathway = rownames(all_data_raw2@assays[["METABOLISM"]][["score"]])[1:20], 
  phenotype = "group", 
  norm = "y"
)


## cluster by metabolisam -----------------
all_data_metabolism <- 
  CreateSeuratObject(
    CreateAssayObject(
      data = as.matrix(all_data_raw2@assays$METABOLISM$score)
    ),
    meta.data = all_data_raw2@meta.data
  )
all_data_metabolism_new_name <- all_data_metabolism %>% colnames() %>% str_remove("^X")
all_data_metabolism <- RenameCells(all_data_metabolism, new.names = all_data_metabolism_new_name)

all_data_metabolism <- 
  all_data_metabolism %>% 
  AddMetaData(cell_coords)

all_data_metabolism <-
  SeuratWrappers::RunBanksy(
    all_data_metabolism,
    lambda  = 0.2,
    verbose = TRUE, 
    assay   = 'RNA', 
    slot    = 'data', 
    assay_name = "BANKSY",
    dimx = 'x',
    dimy = 'y',
    features = 'all',
    k_geom  = 15
  )

# Run PCA and UMAP
all_data_metabolism <- RunPCA(all_data_metabolism, assay = 'BANKSY', features = rownames(all_data_metabolism), npcs = 30)
ElbowPlot(all_data_metabolism)
all_data_metabolism <- RunUMAP(all_data_metabolism, dims = 1:5)

DimPlot(all_data_metabolism, group.by = "sample")

# Clustering
all_data_metabolism <- FindNeighbors(all_data_metabolism, dims = 1:5)
all_data_metabolism <- FindClusters(all_data_metabolism, resolution = c(0.1, 0.2,0.5, 0.8))

# dimplot
pdf(
  file.path(out_dir, 'mebabolisam cluster - umap.pdf'),
  height = 6,
  width  = 8
)
walk(
  c(0.1, 0.2,0.5, 0.8),
  ~ all_data_metabolism %>% DimPlot(group.by = str_glue("BANKSY_snn_res.{.x}")) %>% print
)
dev.off()

# spatial plot
pdf(
  file.path(out_dir, 'mebabolisam cluster - spatial plot.pdf'),
  height = 16,
  width  = 8
)
all_data_metabolism[[c('x', 'y', 'BANKSY_snn_res.0.1')]] %>%
  ggplot(aes(x = x, y = y)) +
  geom_tile(aes(fill = BANKSY_snn_res.0.1), height = 20, width = 20) +
  coord_equal() +
  scale_y_reverse() +
  theme_void()
dev.off()

# dotlot
# DotPlot(
#   all_data_metabolism,
#   features = rownames(all_data_metabolism),
#   group.by = "BANKSY_snn_res.0.1"
# ) +
#   coord_flip()

pdf(
  file.path(out_dir, 'mebabolisam cluster - heatmap of metabolism pathway.pdf'),
  height = 12,
  width  = 6
)
AverageExpression(
  all_data_metabolism,
  group.by = "BANKSY_snn_res.0.1"
) %>%
  .$RNA %>%
  pheatmap::pheatmap(
    .,
    scale = "row",
    border_color = NA,
    labels_col = str_remove(colnames(.), "^g")
  )
dev.off()

# cell props
pdf(
  file.path(out_dir, 'mebabolisam cluster - cell props of celltype and group.pdf'),
  height = 12,
  width  = 6
)
inner_join(
  all_data_metabolism[[c('BANKSY_snn_res.0.1')]] %>% rownames_to_column("cell_id"),
  all_data[[c('celltype')]] %>% rownames_to_column("cell_id"),
  by = "cell_id"
) %>%
  group_by(BANKSY_snn_res.0.1, celltype) %>%
  summarise(count = n()) %>%
  group_by(BANKSY_snn_res.0.1) %>%
  mutate(freq = count/sum(count))  %>%
  
  ggplot(aes(x = BANKSY_snn_res.0.1, y = freq, fill = celltype)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw()

inner_join(
  all_data_metabolism[[c('BANKSY_snn_res.0.1')]] %>% rownames_to_column("cell_id"),
  all_data[[c('group')]] %>% rownames_to_column("cell_id"),
  by = "cell_id"
) %>%
  group_by(BANKSY_snn_res.0.1, group) %>%
  summarise(count = n()) %>%
  group_by(BANKSY_snn_res.0.1) %>%
  mutate(freq = count/sum(count))  %>%
  
  ggplot(aes(x = BANKSY_snn_res.0.1, y = freq, fill = group)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw()
dev.off()





# T  差异基因 代谢酶相关 ----------
t_grp1 <- all_data$group2 == 'LR' & all_data$group3 == 'tumor' & all_data$celltype == 'T'
t_grp2 <- all_data$group2 == 'PH' & all_data$group3 == 'tumor' & all_data$celltype == 'T'

t_grp1_cells <- colnames(all_data)[t_grp1]
t_grp2_cells <- colnames(all_data)[t_grp2]


all_markers_t <- 
  FindMarkers(
    all_data[["RNA"]],
    cells.1 = t_grp1_cells,
    cells.2 = t_grp2_cells
  )

# export de table
all_markers_t %>%
  rownames_to_column("gene") %>%
  mutate(
    pathway = kegg_metabolisam_vector[gene] 
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "T cell - LR vs PH - de table.xlsx"),
    asTable = TRUE
  )

# volcano plot 
p_volcano_t <-
  all_markers_t %>% 
  rownames_to_column("gene") %>% 
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC >  1 ~ "pos",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "neg",
    TRUE ~ "other"
    
  ))  %>%
  mutate(`-log10(p_val_adj)` = -log10(p_val_adj) %>% pmin(320)) %>%
  
  ggplot(aes(x = avg_log2FC,  y = `-log10(p_val_adj)`, color = group)) +
  geom_point(alpha = 0.5) +
  geom_point(
    data = ~ dplyr::filter(., gene %in% kegg_metabolisam_df_tidy$gene),
    shape = 21,
    color = "black",
    fill = NA,
    size = 3,
    show.legend = FALSE,
    alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(., gene %in% kegg_metabolisam_df_tidy$gene, group != "other"),
    aes(label = gene),
    color = "black",
    min.segment.length = 0.5,
    force = 3,
    force_pull = 0.5,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  scale_color_manual(values = c(
    pos = "red",
    neg = "blue",
    other = "gray"
  )) +
  theme_classic() +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "LR vs PH"
  )

# export
ggsave(
  p_volcano_t,
  filename = file.path(out_dir, "T cell - LR vs PH - volcano plot.pdf"),
  height = 6, 
  width = 8
)


#  enrichment 
enrich_metabolism_t <-
  all_markers_t %>% 
  rownames_to_column('gene')  %>% 
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC >  1 ~ "LR",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "PH",
    TRUE ~ "other"
  )) %>%
  dplyr::filter(group != "other") %>%
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  map('gene') %>%
  map(function(x){
    clusterProfiler::enricher(
      x,
      TERM2GENE = kegg_metabolisam_df,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  })

p_enrich_metabolism_t <-
  enrich_metabolism_t  %>%
  imap(function(x, name){
    enrichplot::dotplot(
      x,
      showCategory = 15,
      title = name
    )
  })

# export
pdf(
  file = file.path(out_dir, "T cell - LR vs PH - metabolisam pathway enrichment plot.pdf"),
  height = 6, 
  width = 8
)
p_enrich_metabolism_t %>% walk(print)
dev.off()




gsea_metabolism_t <-
  structure(
    all_markers_t$avg_log2FC,
    names = rownames(all_markers_t)
  ) %>%
  sort(decreasing = TRUE) %>%
  clusterProfiler::GSEA(
    TERM2GENE = kegg_metabolisam_df,
    pvalueCutoff = 0.2
  )

sig_pathway_t <- gsea_metabolism_t@result$ID


p_gsea_t <-
  enrichplot::gseaplot2(
    gsea_metabolism_t,
    geneSetID = sig_pathway_t
  )

# export
ggsave(
  p_gsea_t,
  filename = file.path(out_dir, "T cell - LR vs PH - gsea plot.pdf"),
  height = 6, 
  width = 8
)


 
# tumor 差异基因 代谢酶相关 ----------
tumor_grp1 <- all_data$group2 == 'LR' & all_data$group3 == 'tumor' & all_data$celltype == 'Tumor cell'
tumor_grp2 <- all_data$group2 == 'PH' & all_data$group3 == 'tumor' & all_data$celltype == 'Tumor cell'

tumor_grp1_cells <- colnames(all_data)[tumor_grp1]
tumor_grp2_cells <- colnames(all_data)[tumor_grp2]


all_markers_tumor <- 
  FindMarkers(
    all_data[["RNA"]],
    cells.1 = tumor_grp1_cells,
    cells.2 = tumor_grp2_cells,
    logfc.threshold = 0,
    min.pct = 0
  )

# export de table
all_markers_tumor %>%
  rownames_to_column("gene") %>%
  mutate(
    pathway = kegg_metabolisam_vector[gene] 
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "tumor cell - LR vs PH - de table.xlsx"),
    asTable = TRUE
  )

# volcano plot 
p_volcano_tumor <-
  all_markers_tumor %>% 
  rownames_to_column("gene") %>% 
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC >  1 ~ "pos",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "neg",
    TRUE ~ "other"
    
  ))  %>%
  mutate(`-log10(p_val_adj)` = -log10(p_val_adj) %>% pmin(320)) %>%

  ggplot(aes(x = avg_log2FC,  y = `-log10(p_val_adj)`, color = group)) +
  geom_point(alpha = 0.5) +
  geom_point(
    data = ~ dplyr::filter(., gene %in% kegg_metabolisam_df_tidy$gene),
    shape = 21,
    color = "black",
    fill = NA,
    size = 3,
    show.legend = FALSE,
    alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(., gene %in% kegg_metabolisam_df_tidy$gene, group != "other"),
    aes(label = gene),
    color = "black",
    min.segment.length = 0.5,
    force = 3,
    force_pull = 0.5,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  scale_color_manual(values = c(
    pos = "red",
    neg = "blue",
    other = "gray"
  )) +
  theme_classic() +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "LR vs PH"
  )

# export
ggsave(
  p_volcano_tumor,
  filename = file.path(out_dir, "tumor cell - LR vs PH - volcano plot.pdf"),
  height = 6, 
  width = 8
)


#  enrichment 
enrich_metabolism_tumor <-
  all_markers_tumor %>% 
  rownames_to_column('gene')  %>% 
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC >  1 ~ "LR",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "PH",
    TRUE ~ "other"
  )) %>%
  dplyr::filter(group != "other") %>%
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  map('gene') %>%
  map(function(x){
    clusterProfiler::enricher(
      x,
      TERM2GENE = kegg_metabolisam_df,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  })
  
p_enrich_metabolism_tumor <-
  enrich_metabolism_tumor  %>%
  imap(function(x, name){
    enrichplot::dotplot(
      x,
      showCategory = 15,
      title = name
    )
  })

# export
pdf(
  file = file.path(out_dir, "tumor cell - LR vs PH - metabolisam pathway enrichment plot.pdf"),
  height = 6, 
  width = 8
)
p_enrich_metabolism_tumor %>% walk(print)
dev.off()




gsea_metabolism_tumor <-
  structure(
    all_markers_tumor$avg_log2FC,
    names = rownames(all_markers_tumor)
  ) %>%
  sort(decreasing = TRUE) %>%
  clusterProfiler::GSEA(
    TERM2GENE = kegg_metabolisam_df,
    pvalueCutoff = 0.1
  )

sig_pathway_tumor <- gsea_metabolism_tumor@result$ID


p_gsea_tumor <-
  enrichplot::gseaplot2(
    gsea_metabolism_tumor,
    geneSetID = sig_pathway_tumor
  )

# export
ggsave(
  p_gsea_tumor,
  filename = file.path(out_dir, "tumor cell - LR vs PH - gsea plot.pdf"),
  height = 6, 
  width = 8
)

# hippo-yap pathway -----------
hippo_yap_genes <- c(
  "RASSF1",
  "NF2",
  
  "MST1",
  "STK3",
  "SAV1",
  
  "LATS1",
  "LATS2",
  "CXCL10",
  
  "MAP4K1",
  "MAP4K2",
  "MAP4K3",
  "MAP4K4",
  "MINK1",
  "TNIK",
  
  "NDRG1",
  "STK38L",
  
  "YY1AP1",
  "TAZ",
  "YWHAQ",
  
  "TEAD1",
  "TEAD2",
  "TEAD3",
  "TEAD4"
)

p_tumor_volcano_hippo <-
  all_markers_tumor %>% 
  rownames_to_column("gene") %>%
  mutate(`-log10(p_val_adj)` = -log10(p_val_adj)) %>%
  mutate(`-log10(p_val_adj)` = case_when(
    is.infinite(`-log10(p_val_adj)`) ~ 300,
    TRUE ~ `-log10(p_val_adj)`
  )) %>%
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1 ~ "Pos",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Neg",
    TRUE ~ "non"
  )) %>%
  
  ggplot(aes(x = avg_log2FC, y = `-log10(p_val_adj)`, color = group)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  scale_color_manual(values = c(
    Pos = "red",
    Neg = "blue",
    non = "gray"
  )) +
  ylim(c(-10, 310)) +
  geom_point(
    data = ~ dplyr::filter(., gene %in% hippo_yap_genes),
    size = 3,
    color = "black",
    shape = 21,
    fill = NA
  ) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(., gene %in% hippo_yap_genes),
    aes(label = gene),
    color = "black"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("LR vs PH")

# export
ggsave(
  p_tumor_volcano_hippo,
  filename = file.path(out_dir, "tumor cell - volcano - hippo yap.pdf"),
  height = 6, 
  width = 8
)


# tf gene -----------
# some TF gene
gene_description <- read_tsv("gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

p_tumor_tf_gene <-
  all_markers_tumor %>% 
  rownames_to_column("gene") %>%
  mutate(`-log10(p_val_adj)` = -log10(p_val_adj)) %>%
  mutate(`-log10(p_val_adj)` = case_when(
    is.infinite(`-log10(p_val_adj)`) ~ 300,
    TRUE ~ `-log10(p_val_adj)`
  )) %>%
  mutate(group = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1 ~ "Pos",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Neg",
    TRUE ~ "non"
  )) %>%
  
  ggplot(aes(x = avg_log2FC, y = `-log10(p_val_adj)`, color = group)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  scale_color_manual(values = c(
    Pos = "red",
    Neg = "blue",
    non = "gray"
  )) +
  ylim(c(-10, 310)) +
  geom_point(
    data = ~ dplyr::filter(., gene %in% tf_gene, group != "non"),
    size = 3,
    color = "black",
    shape = 21,
    fill = NA
  ) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(., gene %in% tf_gene, group != "non"),
    aes(label = gene),
    color = "black",
    max.overlaps = 20
  ) +
  geom_label(
    data = ~ dplyr::filter(., gene == "YAP1"),
    aes(label = gene)
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("LR vs PH")

# export
ggsave(
  p_tumor_tf_gene,
  filename = file.path(out_dir, "tumor cell - volcano - tf gene.pdf"),
  height = 6, 
  width = 8
)

# all_data_hippo <-
#   all_data %>% 
#   AddModuleScore(
#     features = hippo_yap_genes,
#     name = hippo_yap,
#     nbin = 1,
#     ctrl = 10
#   )

# umap plot ----------------
p_umap_ct <-
  DimPlot(
    all_data,
    group.by = "celltype"
  ) +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_ct,
  filename = file.path(out_dir, "umap - major annotation.pdf"),
  height = 6, 
  width = 8
)


# spatial plot of YAP1 ------------------
sp_plot_dat <-
  cell_coords %>% 
  rownames_to_column("cell_id") %>%
  left_join(
    FetchData(all_data, vars = c("YAP1", "sample")) %>% rownames_to_column("cell_id"),
    by = "cell_id"
  ) %>%
  dplyr::filter(! is.na(sample))
  
sp_plot_dat_ls <- 
  sp_plot_dat %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

sp_yap1_ls <-
  sp_plot_dat_ls %>%
  imap(function(dat, name){
    dat %>%
      ggplot(aes(x = x, y = y)) +
      geom_tile(
        aes(fill = YAP1, group = cell_id, alpha = YAP1), 
        height = 5, 
        width  = 5,
        show.legend = c(alpha = FALSE, fill = TRUE)
      ) +
      theme_void() +
      coord_equal() +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu")),
        # trans = "log1p"
      ) +
      scale_alpha(range = c(0.3, 1)) +
      ggtitle(name)
  })

# epxort
pdf(
  file = file.path(out_dir, "spatial plot of YAP1.pdf"),
  height = 8,
  width = 8
)
sp_yap1_ls %>% walk(print)
dev.off()



