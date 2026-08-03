library(Seurat)
library(tidyverse)

source("spatial analysis.R")

cn_out_dir <- "CN_analysis - fibro"
dir.create(cn_out_dir)

# grp_cols_sp <-c(
#   Ctrl   = '#1F77B4FF', # Z 
#   CD     = "#D62728FF", # Y 
#   'CD-F' = "#9467BDFF"  # X 
# )



all_data_full2 <- readRDS("all_data_full2.rds")



all_cellid_not_others <- colnames(all_data_full2)[all_data_full2$celltype != "Others"]
cells_position_ls <- 
  all_data_full2@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  setNames(., str_remove(names(.), "\\..*$")) %>%
  map(function(x){
    x[x$cell_id %in% all_cellid_not_others, ]
  })
cells_position_df <- cells_position_ls %>% bind_rows() 




# implement CN ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  all_data_full2[[]] %>%
  dplyr::filter(celltype != "Others") %>% 
  dplyr::mutate(annotation_full = celltype ) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # important!
    obj_id = cell_id # important!
  )

Sys.time()
res_cn <- CN_compute_CN2(
  all_data_meta_for_cn_analysis,  # %>% head(10000),
  
  query =  c("ICAF", "EMT_CAF", "ECAF", "apCAF") ,

  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 10,
  k.cluster      = 10
)
Sys.time()

# saveRDS(res_cn, file.path(cn_out_dir, "res_cn.rds"))

# heatmap ----------------------
# FoldChange cluster/allCell
cn_cluster_mean <- res_cn$cluster_centers
cn_baseMean     <- res_cn$base_mean

# log2(FC)
cn_cluster_FC <- log2(1 + sweep(cn_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn_baseMean))

# heatmap
p_cn_heat <-
  pheatmap::pheatmap(
    cn_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(cn_out_dir, "CN_heatmap_CN.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_heat$gtable)
dev.off()

# CN cell props ---------------
all_data_meta_cn_fibs <-
  all_data_full2[[]] %>%
  rownames_to_column('cell_id') %>%
  dplyr::filter(
    celltype %in% c("ICAF", "EMT_CAF",  "ECAF", "apCAF")
  ) %>%
  left_join(
    rownames_to_column(res_cn$cn_cluster, "cell_id"),
    suffix = c("", "_fibs"),
    by = "cell_id"
  )

p_fibs_cn_props <-
  all_data_meta_cn_fibs %>% 
  group_by(
    celltype, CN_fibs
  ) %>%
  summarise(count = n()) %>%
  group_by(CN_fibs) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = freq, y = CN_fibs, fill = celltype)) +
  geom_col() +
  theme_classic()

ggsave(
  p_fibs_cn_props, 
  filename = file.path(cn_out_dir, "Fibs CN cell props.pdf"),
  height = 4, 
  width  = 8
)


p_fibs_cn_props_CN2 <-
  all_data_meta_cn_fibs %>% 
  group_by(
    CN2, CN_fibs
  ) %>%
  summarise(count = n()) %>%
  group_by(CN_fibs) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = freq, y = CN_fibs, fill = CN2)) +
  geom_col() +
  theme_classic() +
  ggsci::scale_fill_igv()

ggsave(
  p_fibs_cn_props_CN2, 
  filename = file.path(cn_out_dir, "Fibs CN cell props - all cells CN.pdf"),
  height = 4, 
  width  = 8
)


p_fibs_cn_props_group <-
  all_data_meta_cn_fibs %>% 
  group_by(
    group, CN_fibs
  ) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(group = factor(group, levels = c(
    "Normal", "Precancerous", "EarlyStage", "LateStage" 
  ))) %>% 
  
  ggplot(aes(x = group, y = freq, fill = CN_fibs)) +
  geom_col() +
  theme_classic() +
  ggsci::scale_fill_igv() +
  RotatedAxis()

ggsave(
  p_fibs_cn_props_group, 
  filename = file.path(cn_out_dir, "Fibs CN cell props - group.pdf"),
  height = 6, 
  width  = 4
)




# CN annotaton --------------------
# CN k10


cn_annotation <- c(
  "CN1"  = "apCAF core",
  "CN2"  = "ECAF_1",
  "CN3"  = "ECAF core",
  "CN4"  = "ICAF core",
  "CN5"  = "EMT_CAF_tumor_boundary",
  "CN6"  = "ECAF_2",
  "CN7"  = "ICAF_immune",
  "CN8"  = "apCAF_ECAF",
  "CN9"  = "ICAF_muscle",
  "CN10" = "EMT_CAF_tumor_core"
)

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5",
#   "CN6"  = "CN6",
#   "CN7"  = "CN7",
#   "CN8"  = "CN8",
#   "CN9"  = "CN9",
#   "CN10" = "CN10"
# )


p_cn_heat2 <-
  cn_cluster_FC %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = cn_annotation,
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(cn_out_dir, "CN_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_heat2$gtable)
dev.off()


# # spatial plot -----------------
# sp_out_dir <- file.path(cn_out_dir, 'spatial plot - fibro CN') 
# dir.create(sp_out_dir)
# 
# p_sp_plot_anno <-
#   Images(all_data_full2) %>%
#   setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
#   imap(function(img, img_name){
#     if(img_name == "HC") img_name <- "NML"
#     
#     p <- SpatialDimPlot(
#       all_data_full2,
#       images = img, 
#       group.by = "celltype",
#       pt.size.factor	= 2
#     ) +
#       # ggsci::scale_fill_aaas() +
#       ggsci::scale_fill_igv() +
#       guides(fill = guide_legend(override.aes = list(size = 5))) +
#       ggtitle(img_name)
#     p
#   })
# 
# # export
# 
# p_sp_plot_anno %>%
#   iwalk(function(p, name){
#     on.exit(dev.off())
#     pdf(
#       file.path(sp_out_dir, str_glue("spatial plot annotation - {name}.pdf")),
#       height = 8,
#       width  = 8
#     )
#     print(p)
#     # dev.off()
#   })



# CN spatial plot -----------------
cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta_sp_fibs <- 
  all_data_full2[[]] %>%
  
  dplyr::filter(celltype != "Others") %>%
  
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id", suffix = c("", "_fibs")) %>%
  mutate(CN2_fibs = fct_recode(CN_fibs, !!!structure(names(cn_annotation), names = cn_annotation)))  
  # mutate(CN2_fibs = ifelse(is.na(CN2_fibs), as.character(celltype), as.character(CN2_fibs)))  
  

# all_data_meta_sp_fibs <-
#   all_data_meta_sp_fibs %>% 
#   left_join(sample_df, by = c(project = "sample"))

# all_data_meta$group <- 
#   factor(
#     all_data_meta$group,
#     levels = c("Normal", "Precancerous", "EarlyStage", "LateStage")
#   )

saveRDS(all_data_meta_sp_fibs, file.path(cn_out_dir, 'all_data_meta_sp_fibs2.rds') )

cn_sp_p_ls <-
  cells_position_ls  %>% 
  imap(function(cells_position, project){ # browser()
    
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_position %>%
      right_join(all_data_meta_sp_fibs, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      # geom_tile(aes(fill = CN2), height = 10, width = 10, alpha = 0.5)  +
      geom_tile(fill = "gray", height = 10, width = 10, alpha = 0.5)  +
      geom_tile(
        data = ~ dplyr::filter(., !is.na(CN2_fibs)),
        aes(fill = CN2_fibs), height = 10, width = 10
      )  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn_sp_p_ls %>% 
  iwalk(function(p, name){
    on.exit(dev.off())
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
  })


# CN spatial plot2 -----------------
cn_sp_out_dir2 <- file.path(cn_out_dir, 'cn spatial plot2') 
dir.create(cn_sp_out_dir2)

# all_data_meta <- readRDS( file.path(cn_out_dir, 'all_data_meta_cn.rds') )
cn_sp_p_ls2 <-
  cells_position_ls %>% # .[1] %>%
  imap(function(cells_position, project){ # browser()
    cn_color <- structure(
      names = paste0("CN", 1:10),
      ggsci::pal_igv()(10)
    ) %>%
      c(
        Others = "gray90"
      )
    
    dat <-
      cells_position %>%
      right_join(all_data_meta, by = "cell_id")
    
    p_ls <-
      paste0("CN", 1:10) %>% # .[1] %>%
      set_names(., .) %>%
      map(function(cn){
        dat2 <- dat
        
        dat2$CN2 <- as.character(dat2$CN2)
        dat2$CN2[dat2$CN2 != cn] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = y, y = x)) +
          geom_tile(aes(fill = CN2), height = 10, width = 10)  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cn_color) +
          theme_void() +
          labs(fill = NULL) 
        p
      })
    
    p_ls
  })


# export image
cn_sp_p_ls2 %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



# cell props ------

all_cn_count_df <-
  all_data_meta_sp_fibs %>% 
  dplyr::filter(!is.na(CN2_fibs)) %>%
  group_by(project, group, CN2_fibs) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df %>%
  mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  ggplot(aes(x = project, y = freq, fill = CN2_fibs)) +
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
    group   = all_cn_count_df$group %>% unique()
  )
p_cn_props_group <-
  all_cn_count_df %>%
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>%
  mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
  
  ggplot(aes(x = CN2_fibs, y = freq, fill = group)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) + 
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(cn_out_dir, "CN_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
p_cn_props_group %>% print
dev.off()

# cell props - each CN ------
p_cell_props_cn_ls <-
  all_data_meta_sp_fibs %>% 
  group_by(CN2_fibs) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  imap(function(dat, name){
    dat2 <-
      dat %>%
      group_by(group, celltype) %>%
      summarise(count = n()) %>%
      group_by(group) %>%
      mutate(freq = count/sum(count))
    dat2 %>%
      ggplot(aes(x = group, y = freq, fill = celltype)) +
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
  file.path(cn_out_dir, "cell props barplot - CN.pdf"),
  height = 6,
  width  = 6
)
p_cell_props_cn_ls %>% walk(print)
dev.off()


# de analysis of sub cluster of CN_fibs -------
cell_id_cn_fibs <- all_data_meta_cn_fibs$cell_id[!is.na(all_data_meta_cn_fibs$CN_fibs)] 
all_data_full2_cn_fibs <- all_data_full2[, cell_id_cn_fibs]

all_markers_at1_df <-
  sc_at1_2 %>%
  .[, .$annotation != "doublet"] %>%
  JoinLayers() %>%
  {Idents(.) <- "annotation"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1.2
  )

all_markers_at1_df %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj   cluster   gene
# Actg1  3.671699e-104  0.9636172 0.986 0.805 7.365061e-100 Cldn3 AT1  Actg1
# Malat1  1.330267e-93 -1.1808514 0.990 0.992  2.668383e-89 Cldn3 AT1 Malat1
# Cldn18  1.543076e-91  0.7951087 0.994 0.849  3.095256e-87 Cldn3 AT1 Cldn18
# H3f3b   2.011950e-82  0.7796797 0.994 0.859  4.035771e-78 Cldn3 AT1  H3f3b
# Ubb     5.249673e-78  0.8113941 0.980 0.734  1.053032e-73 Cldn3 AT1    Ubb
# Rplp0   4.980331e-77  0.8401779 0.975 0.736  9.990046e-73 Cldn3 AT1  Rplp0


# scatter
p_volcanoplot_AT1_ls <- 
  list("AT1" = all_markers_at1_df) %>%
  imap(function(dat, name){ # browser()
    
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

pdf(
  file.path("subcluster/AT1/", "de scatter plot - AT1 annotation.pdf"), 
  height = 5, 
  width  = 16
)
p_volcanoplot_AT1_ls %>%  walk(print)
dev.off()



# enrichment analysis of sub cluster of AT1 -------

all_geneid_AT1_ls <-
  list("AT1" = all_markers_at1_df) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        print(nrow(.x))
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Mm.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_AT1_ls <- all_geneid_AT1_ls$AT1

all_enrich_AT1_ls <- 
  all_geneid_AT1_ls %>%
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

# enrichment plot
p_all_enrich_AT1_ls <-
  all_enrich_AT1_ls %>%
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
  file.path("subcluster/AT1/", "de enrichment analysis dotplot - sub cluster of AT1.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_AT1_ls %>%  walk(print)
dev.off()


# export enrichemt table

all_enrich_AT1_df <-
  all_enrich_AT1_ls %>%
  map('result') %>%
  enframe(name = "cluster") %>%
  unnest(value) 

all_enrich_AT1_df$geneSymbol <-
  all_enrich_AT1_df$geneID %>%
  str_split("/") %>%
  map_chr(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "ENTREZID",
      toType   = "SYMBOL",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$SYMBOL %>%
      str_c(collapse = "/")
  })

all_enrich_AT1_df %>%
  openxlsx::write.xlsx(
    file.path("subcluster/AT1/", "enrichment analysis table of sub cluster of AT1.xlsx"),
    asTable = TRUE
  )
