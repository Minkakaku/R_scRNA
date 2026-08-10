library(Seurat)
library(tidyverse)

out_dir <- "dat_analysis"
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


##### nestSingleList
nestSingleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  name <- names(x)
  id_ <- seq_along(name)
  stopifnot(all(str_detect(name, sep)))
  name_1 <- str_split(name, pattern = sep) %>% map_chr(~ .x[1])
  name_2 <- str_split(name, pattern = sep) %>% map(~ .x[2:length(.x)]) %>% map_chr(str_c, collapse = sep)
  
  res_ls <- list()
  for(i in id_){    
    n1 <- name_1[i]
    n2 <- name_2[i]
    res <- x[[i]]
    res_ls[[n1]][[n2]] <- res
  }
  res_ls
}



# read data and add roi info ------
all_data_full <- readRDS("all_data_full.rds")
all_data_meta_cn <- readRDS("CN_analysis/all_data_meta_cn2.rds")


meta_all_data_grp <- 
  all_data_full[['project']] %>% 
  rownames_to_column("cell_id") %>% 
  left_join(distinct(all_data_meta_cn[c("project", "group")]), by = "project") %>%
  dplyr::select(-project) %>% 
  column_to_rownames("cell_id") 
 

all_data_full <- 
  all_data_full %>% 
  AddMetaData(meta_all_data_grp)

# saveRDS(all_data_full, 'all_data_full_add_grp.rds')


file.copy(
  "/mnt/NexData_file/home/ztt/test_data/Visium_HD/Gastric_Cancer/ROIgroup.parquet",
  "data"
)

roi_dat <- arrow::read_parquet("data/ROIgroup.parquet")

roi_dat2 <- 
  roi_dat %>%
  mutate(project = id_sample2 %>% str_remove("^.*-1_")) %>%
  mutate(cell_id = str_c(project, id, sep = "_")) %>%
  dplyr::select(
    cell_id, 
    group2 = group, 
    group_all, 
    group_sample
  ) %>% 
  column_to_rownames("cell_id")




all_data_full <-
  all_data_full %>%
  AddMetaData(roi_dat2) %>%
  AddMetaData(column_to_rownames(all_data_meta_cn[c("cell_id", "CN", "CN2")], "cell_id"))


# Others cells not used to compute CN
all_data_full2 <- all_data_full[, !is.na(all_data_full$CN2)]

all_data_full %>%  dim
all_data_meta_cn %>% dim

table(is.na(all_data_full$group2))
table(is.na(all_data_full$group_all))


# saveRDS(all_data_full2, "all_data_full2.rds")

all_data_full2 <- readRDS("all_data_full2.rds")

# cell props ------
all_cn_count_df_group2 <-
  all_data_full2[[]] %>% 
  dplyr::filter(!is.na(group2)) %>%
  group_by(project, group2, CN2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df_group2 %>%
  # mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  ggplot(aes(x = project, y = freq, fill = CN2)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

all_cn_count_all_zeros_df <-
  expand.grid(
    CN2     = all_cn_count_df$CN2 %>% unique(),
    project = all_cn_count_df$project %>% unique(),
    group2   = all_cn_count_df$group2 %>% unique()
  )

p_cn_props_group2 <-
  all_cn_count_df_group2 %>%
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>%
  # mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  
  ggplot(aes(x = CN2, y = freq, fill = group2)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 0.5) + 
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(out_dir, "CN_boxplot_and_barplot - group2.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
p_cn_props_group2 %>% print
dev.off()

# cell props - each CN ------
p_cell_props_cn_ls_group2 <-
  all_data_full2[[]] %>%
  dplyr::filter(!is.na(group2)) %>% 
  group_by(CN2) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat2 <-
      dat %>%
      group_by(group2, celltype) %>%
      summarise(count = n()) %>%
      group_by(group2) %>%
      mutate(freq = count/sum(count))
    dat2 %>%
      ggplot(aes(x = group2, y = freq, fill = celltype)) +
      geom_col() +
      cowplot::theme_cowplot() +
      scale_y_continuous(labels = scales::percent) +
      ggsci::scale_fill_igv() +
      RotatedAxis() +
      labs(
        x = NULL,
        y = NULL,
        title = name
      )
  })

pdf(
  file.path(out_dir, "cell props barplot - CN - group2.pdf"),
  height = 6,
  width  = 6
)
p_cell_props_cn_ls_group2 %>% walk(print)
dev.off()


# cell props - group_all ------
all_cn_count_df_group_all <-
  all_data_full2[[]] %>% 
  dplyr::filter(!is.na(group_all)) %>%
  group_by(project, group_all, CN2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df_group_all %>%
  # mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  ggplot(aes(x = project, y = freq, fill = CN2)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

all_cn_count_all_zeros_df <-
  expand.grid(
    CN2     = all_cn_count_df$CN2 %>% unique(),
    project = all_cn_count_df$project %>% unique(),
    group_all   = all_cn_count_df$group_all %>% unique()
  )

p_cn_props_group_all <-
  all_cn_count_df_group_all %>%
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>%
  # mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  
  ggplot(aes(x = CN2, y = freq, fill = group_all)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 0.5) + 
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(out_dir, "CN_boxplot_and_barplot - group_all.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
p_cn_props_group_all %>% print
dev.off()

# cell props - each CN - group_all ------
p_cell_props_cn_ls_group_all <-
  all_data_full2[[]] %>%
  dplyr::filter(!is.na(group_all)) %>% 
  group_by(CN2) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat2 <-
      dat %>%
      group_by(group_all, celltype) %>%
      summarise(count = n()) %>%
      group_by(group_all) %>%
      mutate(freq = count/sum(count))
    dat2 %>%
      ggplot(aes(x = group_all, y = freq, fill = celltype)) +
      geom_col() +
      cowplot::theme_cowplot() +
      scale_y_continuous(labels = scales::percent) +
      ggsci::scale_fill_igv() +
      RotatedAxis() +
      labs(
        x = NULL,
        y = NULL,
        title = name
      )
  })

pdf(
  file.path(out_dir, "cell props barplot - CN - group_all.pdf"),
  height = 6,
  width  = 6
)
p_cell_props_cn_ls_group_all %>% walk(print)
dev.off()



# de analysis - CN pwc -----
all_data_full2_remove_na <-
  all_data_full2[, !is.na(all_data_full2$group2)] %>%
  {Idents(.) <- "group2"; .} 

all_markers_cn_ls <-
  all_data_full2_remove_na$CN2 %>%
  unique %>%
  set_names(., .) %>%
  map(function(x){
    all_data_full2_remove_na <- all_data_full2_remove_na[, all_data_full2_remove_na$CN2 == x]
    print(str_glue('CN: {x}: {ncol(all_data_full2_remove_na)}'))
    FindAllMarkers(
      JoinLayers(all_data_full2_remove_na),
    )
  })


all_markers_cn_ls$`CN10_Tumor boundary` %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj cluster       gene
# MTRNR2L12  1.259519e-114  2.0699800 0.421 0.139 3.364427e-110    CD-F  MTRNR2L12
# FP236383.3 1.647682e-110  0.5967926 1.000 1.000 4.401288e-106    CD-F FP236383.3
# FP236383.5  4.720774e-79  0.9386205 0.750 0.558  1.261013e-74    CD-F FP236383.5
# FP236383.4  8.018641e-69  0.9385583 0.673 0.449  2.141939e-64    CD-F FP236383.4
# FP671120.7  3.710999e-63  0.6419110 0.977 0.949  9.912819e-59    CD-F FP671120.7
# SLCO5A1     1.146982e-47  0.9227888 0.460 0.258  3.063819e-43    CD-F    SLCO5A1



# scatter
p_volcanoplot_cn_ls <- 
  all_markers_cn_ls %>%
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
  file.path(out_dir, "de scatter plot - multi volcano plot.pdf"), 
  height = 5, 
  width = 15
)
p_volcanoplot_cn_ls %>%  walk(print)
dev.off()




# enrichment analysis -----------

all_geneid_cn_ls <-
  all_markers_cn_ls %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_cn_ls <- 
  all_geneid_cn_ls %>%
  flatDoubleList

all_enrich_cn_ls <- 
  all_geneid_cn_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment plot
p_all_enrich_cn_ls <-
  all_enrich_cn_ls %>%
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
  file.path(out_dir, "de enrichment analysis dotplot.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_cn_ls %>%  walk(print)
dev.off()


# export 
all_enrich_cn_ls %>%
  nestSingleList() %>%
  map(function(x){
    x %>%
      map(as.data.frame) %>%
      enframe("group2") %>%
      unnest(value)
  }) %>%
  
  openxlsx::write.xlsx(
    file = file.path(out_dir, "enrichment analysis table.xlsx"),
    asTable = TRUE
  )

# de analysis - CN9 - cell pwc -----
all_data_full2_remove_na_cn9 <-
  all_data_full2[, !is.na(all_data_full2$group2)] %>%
  {Idents(.) <- "group2"; .} %>%
  .[, .$CN == "CN9"] 

  

all_markers_cn9_ls <-
  all_data_full2_remove_na_cn9$celltype %>%
  unique %>%
  set_names(., .) %>%
  map(function(x){
    all_data_full2_remove_na <- all_data_full2_remove_na_cn9[, all_data_full2_remove_na_cn9$celltype == x]
    print(str_glue('CELL: {x}: {ncol(all_data_full2_remove_na_cn9)}'))
    FindAllMarkers(
      JoinLayers(all_data_full2_remove_na_cn9),
    )
  })


all_markers_cn9_ls$Enterocytes %>% head
# p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
# GKN1   0.000000e+00   3.727350 0.234 0.029 0.000000e+00  NORMAL   GKN1
# GKN2   8.984759e-63   2.554579 0.074 0.014 1.624894e-58  NORMAL   GKN2
# MT-CO2 1.143933e-38  -1.352239 0.264 0.440 2.068802e-34  NORMAL MT-CO2
# LIPF   3.767994e-36   1.699258 0.050 0.011 6.814417e-32  NORMAL   LIPF
# MT-ND4 3.757597e-30  -1.350001 0.204 0.361 6.795614e-26  NORMAL MT-ND4
# MT-CYB 5.587950e-29  -1.258958 0.259 0.403 1.010581e-24  NORMAL MT-CYB



# scatter
p_volcanoplot_cn9_ls <- 
  all_markers_cn9_ls %>%
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
  file.path(out_dir, "de scatter plot - CN9 - multi volcano plot.pdf"), 
  height = 5, 
  width = 15
)
p_volcanoplot_cn9_ls %>%  walk(print)
dev.off()

# enrichment analysis - CN9 -----------

all_geneid_cn9_ls <-
  all_markers_cn9_ls %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_cn9_ls <- 
  all_geneid_cn9_ls %>% 
  flatDoubleList

all_enrich_cn9_ls <- 
  all_geneid_cn9_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment plot
p_all_enrich_cn9_ls <-
  all_enrich_cn9_ls %>%
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
  file.path(out_dir, "de enrichment analysis dotplot.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_cn_ls %>%  walk(print)
dev.off()


# export 
all_enrich_cn_ls %>%
  nestSingleList() %>%
  map(function(x){
    x %>%
      map(as.data.frame) %>%
      enframe("group2") %>%
      unnest(value)
  }) %>%
  
  openxlsx::write.xlsx(
    file = file.path(out_dir, "enrichment analysis table.xlsx"),
    asTable = TRUE
  )

