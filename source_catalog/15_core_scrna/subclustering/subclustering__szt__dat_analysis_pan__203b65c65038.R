library(tidyverse)
library(Seurat)


all_data_sc2 <- readRDS('all_data_sc2.rds')

out_dir <- "20251106 panoptosis"
dir.create(out_dir)

#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}


# read PANoptosis signature --------------
panoptosis_ls <-
  list.files("PANoptosis/", full.names = T) %>%
  setNames(., str_remove(basename(.), "-.*$")) %>%
  map(readxl::read_excel) %>%
  # map(~ dplyr::filter(.x, `Relevance score`))
  map('Gene Symbol')


p_pan_venn <-
  ggvenn::ggvenn(
    panoptosis_ls,
    fill_color = c("#B4FCC2", "#FDBEB4", "#BCB4FF"),
    fill_alpha = 1,
    show_stats = "c",
    stroke_size = 0.2
  )

ggsave(
  p_pan_venn,
  filename = file.path(out_dir, "venn plot of panoptosis.pdf"),
  height = 5,
  width = 6
)

pan_gene <- purrr::reduce(panoptosis_ls, intersect) %>% str_to_title()

# 1 panoptosis signature score ------------
all_data_sc_score <- all_data_sc2 %>% AddModuleScore(features = list(PANoptosis = pan_gene))
all_data_sc_score$PANoptosis <- all_data_sc_score$Cluster1

df_signature <-
  all_data_sc_score %>%
  FetchData(
    vars = c('Age', 'Oxygen', 'Sample', 'annotation', 'PANoptosis')
  )

# Age 
p_heat_pan_age <-
  df_signature %>%
  group_by(Age) %>%
  mutate(Age = factor(Age, levels = c("P3", "P7", "P14"))) %>%
  summarise(PANoptosis = mean(PANoptosis)) %>%
  ggplot(aes(x = Age, y = 'PANoptosis', fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)

p_heat_pan_oxygen <-
  df_signature %>%
  group_by(Oxygen) %>%
  summarise(PANoptosis = mean(PANoptosis)) %>%
  ggplot(aes(x = Oxygen, y = 'PANoptosis', fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)


p_heat_pan_anno <-
  df_signature %>%
  group_by(annotation) %>%
  summarise(PANoptosis = mean(PANoptosis)) %>% 
  mutate(annotation = reorder(annotation, PANoptosis)) %>% 
  ggplot(aes(x = annotation, y = 'PANoptosis', fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  RotatedAxis()


p_heat_pan_sample <-
  df_signature %>%
  group_by(Sample) %>%
  summarise(PANoptosis = mean(PANoptosis)) %>% 
  # mutate(annotation = reorder(annotation, PANoptosis)) %>% 
  ggplot(aes(x = Sample, y = 'PANoptosis', fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  RotatedAxis()

df_pan_sample_anno <- 
  df_signature %>%
  group_by(Sample, annotation) %>%
  summarise(PANoptosis = mean(PANoptosis))

p_heat_pan_sample_anno <-
  df_pan_sample_anno %>% 
  mutate(annotation = reorder(annotation, PANoptosis)) %>%
  ggplot(aes(x = annotation, y = Sample, fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  RotatedAxis()

p_line_pan_sam_anno_ls <-
  df_pan_sample_anno %>%
  # dplyr::filter(annotation %in% c(
  #   "Alveolar Macrophage", 
  #   'Interstitial Macrophage',
  #   'Monocyte',
  #   'DC'
  # )) %>%
  group_by(annotation) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat2 <- dat %>%
      mutate(Time = str_remove(Sample, "_.*$"), Group = str_remove(Sample, "^.*_")) %>%
      mutate(Time = factor(Time, levels = c("P3", "P7", "P14")))
    
    # dat3 <- dat2 %>% 
    #   dplyr::select(-Sample) %>%
    #   arrange(Time) %>%
    #   pivot_wider(
    #     names_from = "Group",
    #     values_from = "PANoptosis"
    #   )
    
    # print(dat3)
    # pval <- t.test(dat3$Normoxia, dat3$Hyperoxia, paired = TRUE)$p.value
    dat2 %>% 
      ggplot(aes(x = Time, y = PANoptosis, color = Group, group = Group)) +
      geom_line() +
      theme_classic() +
      ggtitle(name)  
  })



# export
ggsave(
  p_heat_pan_age,
  filename = file.path(out_dir, "1 panoptosis expression - Age.pdf"),
  height = 2,
  width = 6
)
ggsave(
  p_heat_pan_oxygen,
  filename = file.path(out_dir, "1 panoptosis expression - Oxygen.pdf"),
  height = 2,
  width = 6
)
ggsave(
  p_heat_pan_anno,
  filename = file.path(out_dir, "1 panoptosis expression - annotation.pdf"),
  height = 3,
  width = 6
)
ggsave(
  p_heat_pan_sample,
  filename = file.path(out_dir, "1 panoptosis expression - sample.pdf"),
  height = 3,
  width = 6
)
ggsave(
  p_heat_pan_sample_anno,
  filename = file.path(out_dir, "1 panoptosis expression - sample - annotation.pdf"),
  height = 6,
  width = 6
)
pdf(
  file = file.path(out_dir, "1 panoptosis expression - line - annotation and sample.pdf"),
  height = 4,
  width = 6
)
p_line_pan_sam_anno_ls %>% walk(print)
dev.off()


# 2 de analysis ----------------
all_data_sc2 <- JoinLayers(all_data_sc2)
all_marker_oxy_groupby_age_ls <-
  all_data_sc2$Age %>%
  unique() %>%
  setNames(., .) %>%
  map(function(x){
    Idents(all_data_sc2) <- "Age"
    FindMarkers(
      all_data_sc2,
      ident.1 = "Hyperoxia",
      ident.2 = "Normoxia",
      group.by = "Oxygen",
      subset.ident = x,
      logfc.threshold = 0
    )
  })
  
all_marker_oxy_groupby_age_ls %>%
  map(rownames_to_column, 'gene') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "de analysis of Oxygen groupby Age.xlsx"),
    asTable = TRUE
  )


# scatter
p_volcanoplot_oxy_groupby_age_ls <- 
  all_marker_oxy_groupby_age_ls %>%
  map(rownames_to_column, 'gene') %>%
  enframe("cluster") %>%
  unnest(value) %>%
  list(Hyperoxia_vs_Normoxia = .) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
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
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
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
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      ggsci::scale_fill_aaas() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })




pdf(
  file.path(out_dir, "multi volcano plot of oxygen group by age.pdf"), 
  height = 5, 
  width = 8
)
p_volcanoplot_oxy_groupby_age_ls %>%  walk(print)
dev.off()



# 3 enrichment analysis -----------
all_de_gene_ls <- 
  all_marker_oxy_groupby_age_ls %>%
  map(function(df){
    up_gene <-
      df %>%
      rownames_to_column("gene") %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      dplyr::pull('gene')
    
    down_gene <-
      df %>%
      rownames_to_column("gene") %>%
      dplyr::filter(avg_log2FC < -0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      dplyr::pull('gene')
    
    list(
      Hyperoxia = up_gene,
      Normoxia = down_gene
    )
  }) %>%
  flatDoubleList()


all_de_geneid_ls <-
  all_de_gene_ls %>%
  map(function(gene){
    clusterProfiler::bitr(
      gene,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$ENTREZID

  })

# enrichment 
all_enrich_ls <- 
  all_de_geneid_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })


all_enrich_ls <-
  all_enrich_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Mm.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_ls <-
  all_enrich_ls %>%
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



pdf(
  file.path(out_dir, "2 de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_ls %>%  walk(print)
dev.off()

all_enrich_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "enrichment analysis of Oxygen groupby Age.xlsx"),
    asTable = TRUE
  )

# 4 cell props -----------
df_freq_major <- 
  all_data_sc2[[]] %>%
  group_by(Sample, annotation) %>%
  summarise(count = n()) %>%
  group_by(Sample) %>%
  mutate(freq = count / sum(count)) 

p_col_major_anno <- 
  df_freq_major %>%
  ggplot(aes(x = Sample, y = freq, fill = annotation)) +
  geom_col(width = 0.75) +
  theme_classic() +
  ggsci::scale_fill_igv() +
  RotatedAxis() +
  labs(
    x = NULL
  )
  
p_line_major_anno_ls <-
  df_freq_major %>%
  # dplyr::filter(annotation %in% c(
  #   "Alveolar Macrophage", 
  #   'Interstitial Macrophage',
  #   'Monocyte',
  #   'DC'
  # )) %>%
  group_by(annotation) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat2 <- dat %>%
      mutate(Time = str_remove(Sample, "_.*$"), Group = str_remove(Sample, "^.*_")) %>%
      mutate(Time = factor(Time, levels = c("P3", "P7", "P14")))
    
    dat2 %>% 
      ggplot(aes(x = Time, y = freq, color = Group, group = Group)) +
      geom_line() +
      theme_classic() +
      ggtitle(name)  
  })

# export
ggsave(
  p_col_major_anno,
  filename = file.path(out_dir, "4 cell props - barplot - major annotation.pdf"),
  height = 8,
  width = 6
)


pdf(
  file = file.path(out_dir, "4 cell props - line - major annotation.pdf"),
  height = 4,
  width = 6
)
p_line_major_anno_ls %>% walk(print)
dev.off()

# 5. umap --------------
p_umap_pan <-
  FeaturePlot(
    all_data_sc_score,
    features  = "PANoptosis",
    reduction = 'harmonyumap',
    label = TRUE
  ) +
  # NoLegend() +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds") %>% c("white","white", .))

p_umap_pan_sample <-
  FeaturePlot(
    all_data_sc_score,
    features  = "PANoptosis",
    reduction = 'harmonyumap',
    split.by = "Sample",
    # ncol = 2,
    label = TRUE
  ) &
  # NoLegend() +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds") %>% c("white","white", .))


p_umap_major <-
  DimPlot(
    all_data_sc2,
    group.by = "annotation",
    reduction = 'harmonyumap',
    label = TRUE
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

p_umap_major_oxy <-
  DimPlot(
    all_data_sc2,
    group.by = "annotation",
    reduction = 'harmonyumap',
    split.by = "Oxygen",
    label = TRUE
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

p_umap_major_age <-
  DimPlot(
    all_data_sc2,
    group.by = "annotation",
    reduction = 'harmonyumap',
    split.by = "Age",
    label = TRUE
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

p_umap_major_sample <-
  DimPlot(
    all_data_sc2,
    group.by = "annotation",
    reduction = 'harmonyumap',
    split.by = "Sample",
    label = TRUE,
    ncol = 2
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_major_sample,
  filename = file.path(out_dir, "5 umap - Panoptosis.pdf"),
  height = 6,
  width = 8
)

ggsave(
  p_umap_pan_sample,
  filename = file.path(out_dir, "5 umap - Panoptosis - sample.pdf"),
  height = 6,
  width = 36
)

ggsave(
  p_umap_major,
  filename = file.path(out_dir, "5 umap - major annotation.pdf"),
  height = 6,
  width = 8
)



ggsave(
  p_umap_major_oxy,
  filename = file.path(out_dir, "5 umap - major annotation - oxy.pdf"),
  height = 6,
  width = 12
)



ggsave(
  p_umap_major_age,
  filename = file.path(out_dir, "5 umap - major annotation - age.pdf"),
  height = 6,
  width = 18
)

ggsave(
  p_umap_major_sample,
  filename = file.path(out_dir, "5 umap - major annotation - sample.pdf"),
  height = 18,
  width = 12
)


# 6. GSVA analysis ---------------
all_panoptosis <- c(list(PANoptosis = pan_gene), panoptosis_ls)
all_panoptosis2 <- 
  all_panoptosis %>% 
  map(str_to_title) %>%
  map(intersect, rownames(all_data_sc2))


all_log2fc_oxy_groupby_age_ls <-
  all_marker_oxy_groupby_age_ls %>%
  map(function(x){
    log2fc <- x$avg_log2FC 
    names(log2fc) <- rownames(x)
    log2fc <- sort(log2fc, decreasing = T)
    log2fc
  })


all_gsea_oxy_groupby_age_ls <- 
  
  all_log2fc_oxy_groupby_age_ls %>%
  map(function(x){
    clusterProfiler::GSEA(
      geneList = x,
      TERM2GENE = all_panoptosis2 %>% enframe('term', 'name') %>% unnest(name),
      maxGSSize = 20000,
      pvalueCutoff = 1
    )
  })


# all_gsea_oxy_groupby_age_ls$P3 %>%
#   enrichplot::gseaplot2(
#     geneSetID = all_panoptosis2 %>% names()
#   )

df_all_gsea_oxy_groupby_age <-
  all_gsea_oxy_groupby_age_ls %>%
  map('result') %>%
  enframe('Age') %>%
  unnest(value)

# df_all_gsea_oxy_groupby_age %>%
#   mutate(Age = factor(Age, levels = c('P3', 'P7', 'P14'))) %>% 
#   mutate(sig = as.character(pvalue < 0.05)) %>%
#   
#   ggplot(aes(x = NES, y = Description, color = sig, fill = Age)) +
#   geom_col(position = 'dodge', width = 0.75) +
#   theme_classic() +
#   ggsci::scale_fill_aaas() +
#   scale_color_manual(
#     values = c(
#       'TRUE' = "red",
#       'FALSE' = NA
#     ),
#     na.value = NA
#   ) +
#   guides(color = FALSE)  

p_gsea_panopotosis <-
  df_all_gsea_oxy_groupby_age %>%
  mutate(Age = factor(Age, levels = c('P3', 'P7', 'P14'))) %>% 
  mutate(sig = as.character(pvalue < 0.05)) %>%
  group_by(Age) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat %>% 
      ggplot(aes(x = NES, y = Description)) +
      geom_col(fill = '#F3CBC2', width = 0.75) +
      theme_classic() +
      guides(color = FALSE) +
      geom_text(
        data = ~ dplyr::filter(., NES > 0),
        aes(x = 0, label = Description, color = sig),
        hjust = - 0.1,
        size = 8
      ) +
      geom_text(
        data = ~ dplyr::filter(., NES < 0),
        aes(x = 0, label = Description, color = sig),
        hjust = 1.1,
        size = 8
      ) +
      geom_vline(xintercept = 0) +
      scale_color_manual(values = c(
        'TRUE'  = "red",
        'FALSE' = "black"
      )) +
      theme_void() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.line.x = element_line(color = "black"),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 18)
      ) +
      labs(
        x = "Normalized enrichment score",
        title = name
      ) 
  })


ggsave(
  patchwork::wrap_plots(p_gsea_panopotosis, ncol = 1),
  filename = file.path(out_dir, "6 GSVA - hyper vs normal - split by Age.pdf"),
  height = 12,
  width = 6
)


# 7. de analysis of alveolar mac --------------
## de - mac ------------------
all_marker_al_mac_oxy_groupby_age_ls <-
  all_data_sc2$Age %>%
  unique() %>%
  setNames(., .) %>%
  map(function(x){
    sc <- all_data_sc2[, all_data_sc2$annotation == 'Alveolar Macrophage']
    Idents(all_data_sc2) <- "Age"
    FindMarkers(
      all_data_sc2,
      ident.1 = "Hyperoxia",
      ident.2 = "Normoxia",
      group.by = "Oxygen",
      subset.ident = x,
      logfc.threshold = 0
    )
  })

all_marker_al_mac_oxy_groupby_age_ls %>%
  map(rownames_to_column, 'gene') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "de analysis of Oxygen groupby Age - alveolar mac.xlsx"),
    asTable = TRUE
  )

# scatter
p_volcanoplot_al_mac_oxy_groupby_age_ls <- 
  all_marker_al_mac_oxy_groupby_age_ls %>%
  map(rownames_to_column, 'gene') %>%
  enframe("cluster") %>%
  unnest(value) %>%
  list(Hyperoxia_vs_Normoxia = .) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
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
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
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
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      ggsci::scale_fill_aaas() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })




pdf(
  file.path(out_dir, "7 multi volcano plot of oxygen group by age - mac.pdf"), 
  height = 5, 
  width = 8
)
p_volcanoplot_al_mac_oxy_groupby_age_ls %>%  walk(print)
dev.off()


## enrichment ---------------
all_de_gene_mac_ls <- 
  all_marker_al_mac_oxy_groupby_age_ls %>%
  map(function(df){
    up_gene <-
      df %>%
      rownames_to_column("gene") %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      dplyr::pull('gene')
    
    down_gene <-
      df %>%
      rownames_to_column("gene") %>%
      dplyr::filter(avg_log2FC < -0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      dplyr::pull('gene')
    
    list(
      Hyperoxia = up_gene,
      Normoxia = down_gene
    )
  }) %>%
  flatDoubleList()


all_de_geneid_mac_ls <-
  all_de_gene_mac_ls %>%
  map(function(gene){
    clusterProfiler::bitr(
      gene,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$ENTREZID
    
  })

# enrichment 
all_enrich_mac_ls <- 
  all_de_geneid_mac_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })


all_enrich_mac_ls <-
  all_enrich_mac_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Mm.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_mac_ls <-
  all_enrich_mac_ls %>%
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



pdf(
  file.path(out_dir, "7 de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_mac_ls %>%  walk(print)
dev.off()

all_enrich_mac_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "enrichment analysis of Oxygen groupby Age - mac.xlsx"),
    asTable = TRUE
  )

# 8. msigdb hallmark -----
h_geneset <- msigdbr::msigdbr(category = "H")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

sc_mac_score  <-
  all_data_sc2[, all_data_sc2$annotation == "Alveolar Macrophage"] %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

sc_mac_score  %>% head

# msigdb anno
h_score_sample <- 
  sc_mac_score[[]] %>%
  dplyr::select(Sample, Cluster1:Cluster50) %>%
  group_by(Sample) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_score_sample %>% column_to_rownames("Sample") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path(out_dir, "8 Hallmark pathway score of Alveolar Macrophage - sample.pdf"), 
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()

# 9. PPI analysis --------
all_de_gene_mac_ls['P14---Hyperoxia'] %>%
  map(function(x){
    intersect(x, pan_gene)
  }) %>%
  .[[1]] %>%
  walk(cat, "\n")
  


DotPlot(
  all_data_sc2[, all_data_sc2$annotation == "Alveolar Macrophage"],
  features = c(
    'Casp3',
    'Stat3',
    'Hsp90aa1'
  ),
  group.by = "Sample"
)


# ppi 
ppi_core_gene <- list(
  Betweenness = c("Hsp90aa1", "Myd88", "Casp3", "Stat3", "Tubb6", "Nfkbia"),
  Stress = c("Stat3", "Hsp90aa1", "Myd88", "Casp3", "Casp4", "Calm3"),
  BottleNeck = c("Casp3", "Stat3", "Hsp90aa1", "Myd88", "Nfkbia"),
  Closeness = c("Casp3", "Stat3", "Hsp90aa1", "Myd88", "Nfkbia"),
  Degree = c("Nlrp3", "Sqstm1", "Lmna", "Nfe2", "Nfkbia"),
  DMNC = c("Casp3", "Stat3", "Hsp90aa1", "Myd88", "Nfkbia"),
  EcCentricity = c("Casp3", "Nfkbia", "Stat3", "Hsp90aa1", "Myd88"),
  EPC = c("Casp3", "Stat3", "Hsp90aa1", "Myd88", "Nfkbia"),
  MCC = c("Stat3", "Casp3", "Nfkbia", "Hsp90aa1", "Myd88"),
  MNC = c("Casp3", "Stat3", "Hsp90aa1", "Myd88", "Nfkbia"),
  Radiality = c("Casp3", "Hsp90aa1", "Myd88", "Stat3", "Nfkbia")
)

ppi_core_gene %>% unlist %>% unique %>% str_to_title() %>% .[! . %in% rownames(all_data_sc2)]

ppi_core_gene %>% unlist %>% table %>% sort(decreasing = T)

# top 10
ppi_gene_core_10 <-
  readxl::read_excel("基因节点排名工作簿.xlsx") %>%
  .[-1] %>%
  map(c)

ppi_gene_core_10 %>% unlist %>% unique %>% .[! . %in% rownames(all_data_sc2)]

ppi_gene_core_10 %>% unlist %>% table %>% sort(decreasing = T)


# 10. cellchat ------------
DotPlot(
  all_data_sc2, # [, all_data_sc2$Age == "P14"],
  group.by = "annotation",
  features = c('Lair1',	'Lilrb4a')
)

DotPlot(
  all_data_sc2[, all_data_sc2$Age == "P14"],
  group.by = "annotation",
  features = c('Lair1',	'Lilrb4a')
)

DotPlot(
  all_data_sc2,
  group.by = "Sample",
  features = c('Lair1',	'Lilrb4a')
)

# panoptosis and each optosis signature ------------

set.seed(1234)
all_data_sc_score2 <- 
  all_data_sc2 %>%
  AddModuleScore(
    features = all_panoptosis2
  )

# Age 
h_score_age <- 
  all_data_sc_score2[[]] %>%
  dplyr::select(Age, Cluster1:Cluster4) %>%
  group_by(Age) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_score_age %>% column_to_rownames("Age") %>% t,
    scale = 'row',
    labels_row = names(all_panoptosis),
    cluster_cols = F,
    cluster_rows = F,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )
# export 
pdf(
  file.path(out_dir, "1 panoptosis expression - Age2.pdf"),
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()
    


 
all_data_sc_score2 %>%
  FetchData(
    vars = c('Age', 'Oxygen', 'Sample', 'annotation', 'Cluster1')
  ) %>%
 
  group_by(Age) %>%
  mutate(Age = factor(Age, levels = c("P3", "P7", "P14"))) %>%
  summarise(PANoptosis = mean(Cluster1)) %>%
  ggplot(aes(x = Age, y = 'PANoptosis', fill = PANoptosis)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)



# 20251124 update geneset ---------------
path_kegg <-
  msigdbr::msigdbr(
    db_species = "HS",
    species = "human",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  )

path_kegg2 <-
  msigdbr::msigdbr(
    db_species = "HS",
    species = "human",
    collection = "C2",
    subcollection = "CP:KEGG_MEDICUS"
  )

path_reactome <-
  msigdbr::msigdbr(
    db_species = "HS",
    species = "human",
    collection = "C2",
    subcollection = "CP:REACTOME"
  )

path_gobp <-
  msigdbr::msigdbr(
    db_species = "HS",
    species = "human",
    collection = "C5",
    subcollection = "GO:BP"
  )

# msigdbr::msigdbr_collections('HS')


path_db <- rbind(
  path_kegg,
  # path_kegg2,
  path_reactome,
  path_gobp
)


# APOTOSIS
path_db$gs_name %>% unique %>% str_subset("APOPTOSIS")
"KEGG_APOPTOSIS"                                                                                                                               
"REACTOME_APOPTOSIS" 

# necroptosis
path_db$gs_name %>% unique %>% str_subset("NECROPTO")
"GOBP_NECROPTOTIC_PROCESS"
"GOBP_NECROPTOTIC_SIGNALING_PATHWAY" 

# pyroptosis
path_db$gs_name %>% unique %>% str_subset("PYROPTO")
"REACTOME_PYROPTOSIS" 
"GOBP_PYROPTOTIC_CELL_DEATH" 

PANoptosis_path_ls <-
  path_db %>%
  dplyr::filter(gs_name %in% c(
    "KEGG_APOPTOSIS",
    "REACTOME_APOPTOSIS",
    "GOBP_NECROPTOTIC_PROCESS",
    "GOBP_NECROPTOTIC_SIGNALING_PATHWAY",
    "REACTOME_PYROPTOSIS",
    "GOBP_PYROPTOTIC_CELL_DEATH"
  )) %>%
  with(split(gene_symbol, gs_name)) 

# interset gene, change to mm label
PANoptosis_path_ls2 <-
  PANoptosis_path_ls %>%
  map(function(x){
    xx <- str_to_title(x)
    xx <- xx[xx %in% rownames(all_data_sc2)]
    xx
  })

pan_gene2 <- purrr::reduce(PANoptosis_path_ls2, intersect)

all_data_sc_score2 <- 
  all_data_sc2 %>%
  AddModuleScore(
    features = PANoptosis_path_ls2
  )

## Age --------------
h_score_age <- 
  all_data_sc_score2[[]] %>%
  dplyr::select(Age, Cluster1:Cluster6) %>%
  group_by(Age) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_score_age %>% column_to_rownames("Age") %>% t,
    scale = 'row',
    labels_row = names(PANoptosis_path_ls2),
    cluster_cols = F,
    cluster_rows = F,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )
# export 
pdf(
  file.path(out_dir, "1 panoptosis expression - Age2.pdf"),
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()


