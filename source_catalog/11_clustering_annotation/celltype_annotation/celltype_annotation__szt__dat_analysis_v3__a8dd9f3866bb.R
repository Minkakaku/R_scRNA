# marker gene using findmarkers
# cell vs others

library(tidyverse)
library(Seurat)

library(clusterProfiler)

# out dir
out_dir_analysis <- "out_dir_20250711"
dir.create(out_dir_analysis)

# using new annotation
all_data_full <- readRDS("all_data_full.new.annos.rds")
all_data_full <- JoinLayers(all_data_full)

all_marker <-
  FindAllMarkers(
    all_data_full,
    group.by = "group",
    return.thresh = 1.1
  )


all_de_genes_for_enrichmenet <-
  all_marker %>%
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    # pct.1 - pct.2 > 0.1,
    avg_log2FC > 1
  ) %>%
  nest() %>%
  deframe() %>%
  map('gene')

# enrichemt analysis -----
all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)

## KEGG ------------
enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(out_dir_analysis, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(out_dir_analysis, str_glue("KEGG - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
    ggtitle(name)
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(out_dir_analysis, str_glue("dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })

## GOBP -------
gobp_enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = 'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2
    )
  })

# saveRDS(gobp_enrich_ls, file.path(out_dir_analysis, 'gobp_enrich_ls.rds'))

# expot enrichment table
gobp_enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(out_dir_analysis, str_glue("GOBP - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
gobp_enrich_dot_ls <-
  gobp_enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
gobp_enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(out_dir_analysis, str_glue("BOBP - dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })
# multi volcano -------------
all_marker_ls <- 
  list(
    Total = all_marker
  )
all_marker_ls$Total %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj cluster       gene
# MTRNR2L12  1.259519e-114  2.0699800 0.421 0.139 3.364427e-110    CD-F  MTRNR2L12
# FP236383.3 1.647682e-110  0.5967926 1.000 1.000 4.401288e-106    CD-F FP236383.3
# FP236383.5  4.720774e-79  0.9386205 0.750 0.558  1.261013e-74    CD-F FP236383.5
# FP236383.4  8.018641e-69  0.9385583 0.673 0.449  2.141939e-64    CD-F FP236383.4
# FP671120.7  3.710999e-63  0.6419110 0.977 0.949  9.912819e-59    CD-F FP671120.7
# SLCO5A1     1.146982e-47  0.9227888 0.460 0.258  3.063819e-43    CD-F    SLCO5A1

set.seed(1234)
rect_height <- 1
dodge_width <- 0.1

# scatter
df_scatter <-
  all_marker_ls$Total %>%
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
    group = ifelse(p_val_adj > 0.001, "p > 0.001", "p < 0.001"),
    group = factor(group, levels = c("p > 0.001", "p < 0.001"))
  ) %>%
  dplyr::arrange(group)

# df_label 
df_label_up <- df_scatter %>%
  dplyr::filter(as.character(group) == "p < 0.001") %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10)

df_label_down <- df_scatter %>%
  dplyr::filter(as.character(group) == "p < 0.001") %>%
  group_by(cluster) %>% 
  slice_min(avg_log2FC, n = 10)



# grp rect  
df_grp_rect <-
  all_marker_ls$Total %>%
  group_by(cluster) %>%
  summarise() %>%
  mutate(x_base = as.numeric(factor(cluster))) %>%
  mutate(xmin = x_base - 1, xmax = x_base) %>%
  mutate(ymin = -rect_height, ymax = rect_height)


# background rect
df_bg_rect <-
  all_marker_ls$Total %>%
  group_by(cluster) %>%
  summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
  mutate(x_base = as.numeric(factor(cluster))) %>%
  mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)


# text label
df_text_grp <-
  df_scatter %>%
  group_by(cluster) %>%
  summarise() %>%
  mutate(x_base = as.numeric(factor(cluster))) %>%
  mutate(x = x_base - 1) %>%
  mutate(x = (x_base + x)/2)


p_multi_volcano <-
  df_scatter %>%
  .[sample(1:nrow(.), nrow(.)),] %>%
  ggplot(aes(x = xx, y = y)) +
  geom_rect(
    data = df_bg_rect, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray95",
    inherit.aes = FALSE
  ) +
  geom_point(aes(color = group), size = 0.1) +
  
  geom_rect(
    data = df_grp_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_text_grp,
    aes(x = x, y = 0, label = cluster)
  ) +
  
  ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
  ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
  scale_color_manual(values = c(
    'p < 0.001' = "red",
    'p > 0.001' = "gray70"
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
    y = "Log2FC"
  )

# export
ggsave(
  p_multi_volcano,
  filename = file.path(out_dir_analysis, "multi volcano plot - group.pdf"),
  height = 6,
  width = 15
)


# plot by group --------
all_markers_fib_df <-
  sc_fib %>%
  JoinLayers() %>%
  {Idents(.) <- "annotation"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1.2
  )

all_markers_fib_df %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj        cluster     gene
# EBF1      0.000000e+00   3.184000 0.690 0.179  0.000000e+00 S3 fibroblasts     EBF1
# MIR99AHG 2.867768e-269   2.796024 0.619 0.153 7.660382e-265 S3 fibroblasts MIR99AHG
# KCNN3    5.304018e-180   4.171928 0.391 0.073 1.416809e-175 S3 fibroblasts    KCNN3
# ADGRL3   4.410816e-151  -4.038387 0.083 0.451 1.178217e-146 S3 fibroblasts   ADGRL3
# HPSE2    1.850722e-148  -3.198233 0.107 0.485 4.943647e-144 S3 fibroblasts    HPSE2
# C7       3.038981e-142   2.497229 0.484 0.165 8.117726e-138 S3 fibroblasts       C7


# scatter
p_volcanoplot_fibro_ls <- 
  list("fibroblast" = all_markers_fib_df) %>%
  imap(function(dat, name){browser()
    
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
    
    
    # text label
    df_text_grp <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(x = x_base - 1) %>%
      mutate(x = (x_base + x)/2)
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      
      geom_point(aes(color = group), size = 1) +
      
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_text_grp,
        aes(x = x, y = 0, label = cluster)
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      scale_fill_manual(values =  c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4")) +
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

# marker expre ------------
plot_dat <-
  FetchData(
    all_data_full2,
    vars = "IL10"
  )

# plot_dat <- cbind( all_data_norm_all['IL10',], all_data_meta_all[c('group', "celltype")] )
# colnames(plot_dat)[1] <- "IL10"


plot_dat2 <- cbind( plot_dat, all_data_full2[[c('group', "celltype")]] )



plot_dat2 %>% 
  # dplyr::filter(IL10 != 0) %>%
  ggplot(aes(x = celltype, y = IL10)) +
  geom_violin(trim = FALSE) +
  geom_point(data = ~ dplyr::filter(., IL10 != 0)) +
  RotatedAxis()


