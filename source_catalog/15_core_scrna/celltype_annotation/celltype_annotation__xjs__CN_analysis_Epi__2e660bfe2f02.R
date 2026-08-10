library(Seurat)
library(tidyverse)

source("spatial analysis.R")

cn_out_dir <- "CN_analysis_epi"
dir.create(cn_out_dir)

# grp_cols_sp <-c(
#   Ctrl   = '#1F77B4FF', # Z 
#   CD     = "#D62728FF", # Y 
#   'CD-F' = "#9467BDFF"  # X 
# )



all_data_full2 <- readRDS("~/R_OUTPUT/all_data_4slides_sample_filtered_reUMAP.rds")



# all_cellid_not_others <- colnames(all_data_full2)[all_data_full2$annotation != "Others"]
# cells_position_ls <- 
#   all_data_full2@images %>%
#   map(GetTissueCoordinates) %>%
#   map(dplyr::rename, cell_id = cell) %>% 
#   setNames(., str_remove(names(.), "\\..*$")) %>%
#   map(function(x){
#     x[x$cell_id %in% all_cellid_not_others, ]
#   })
cells_position_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_position_ls_4slides_filtered.rds')
cells_position_df <- cells_position_ls %>% bind_rows() 




# implement CN ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  all_data_full2[[]] %>%
  # dplyr::filter(annotation != "Others") %>% 
  dplyr::mutate(annotation_full = annotation ) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation_full = factor(annotation_full), # important!
    obj_id = cell_id # important!
  )

# saveRDS(all_data_meta_for_cn_analysis, file.path(cn_out_dir, 'all_data_meta_for_cn_analysis.rds') )


Sys.time()
res_cn <- CN_compute_CN2(
  all_data_meta_for_cn_analysis, # %>% head(10000),
  
  query =  c("Epithelial") ,
  
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation_full",
  k.nn           = 10,
  k.cluster      = 5
)
Sys.time()

saveRDS(res_cn, file.path(cn_out_dir, "res_cn5_2.rds"))
all_data_full2[[]][(res_cn$cn_cluster %>% rownames()), 'annotation'] %>% unique


# cell_id_epi <- all_data_meta_for_cn_analysis$cell_id[all_data_meta_for_cn_analysis$annotation == "Epithelial"]




# get NN cells of each EPI-CN cells  ----------------------
red_nn_epi_cn_ls <-
  all_data_meta_for_cn_analysis$sample %>%
  unique %>%
  setNames(., .) %>% 
  map(function(samp){
    idx1 <- all_data_meta_for_cn_analysis$sample == samp
    idx2 <- all_data_meta_for_cn_analysis$annotation_full == "Epithelial"
    
    dat   <- all_data_meta_for_cn_analysis[idx1, ]
    query <- all_data_meta_for_cn_analysis[idx1 & idx2,  ]
    
    sp_cols <- c("X", "Y")
    
    red_nn_epi_cn <- RANN::nn2(
      data  = dat[, sp_cols], 
      query = query[, sp_cols],
      k = 10
    )
    
    # return
    list(
      res_cn = red_nn_epi_cn,
      all_cells = dat$obj_id,
      query_cells = query$obj_id
    )
  })

# cell stats
res_nn_byCN_bySample_ls_ls <-
  all_data_meta_sp_epis2$CN2 %>% # .[1] %>%
  unique() %>%
  .[!is.na(.)] %>%
  setNames(., .) %>%
  map(function(cn){
    print(cn)
    red_nn_epi_cn_ls %>%
      imap(function(res_ls, name){
        print(name)
        cn_name <- all_data_meta_sp_epis2$CN2[match(res_ls$query_cells, all_data_meta_sp_epis2$cell_id)]
        
        nn_idx <- res_ls$res_cn$nn.idx[cn_name == cn, ] %>% c() %>% unique()
        nn_cell_id <- res_ls$all_cells[nn_idx]
        
        # return all nn cell_id
        return(nn_cell_id)
      })
  })
res_nn_byCN_ls <- res_nn_byCN_bySample_ls_ls %>% map(unlist, use.names = FALSE)
res_nn_byCN_ls %>% map(~ sum(duplicated(.)))
saveRDS(res_nn_byCN_ls, file.path(cn_out_dir, 'res_nn_byCN_ls.rds'))


# 
res_nn_byCN_all_dup_cells <- res_nn_byCN_ls %>% unlist(use.names = FALSE)  %>% .[duplicated(.)] %>% unique()
res_nn_byCN_uniq_ls <- res_nn_byCN_ls %>% map(~ .[! . %in% res_nn_byCN_all_dup_cells])
saveRDS(res_nn_byCN_uniq_ls, file.path(cn_out_dir, 'res_nn_byCN_uniq_ls.rds'))




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
  file.path(cn_out_dir, "CN_heatmap_CN5.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_heat$gtable)
dev.off()

# # CN cell props ---------------
# all_data_meta_cn_fibs <-
#   all_data_full2[[]] %>%
#   rownames_to_column('cell_id') %>%
#   dplyr::filter(
#     annotation %in% c("Epithelial")
#   ) %>%
#   left_join(
#     rownames_to_column(res_cn$cn_cluster, "cell_id"),
#     # suffix = c("", "_epis"),
#     by = "cell_id"
#   )
# 
# p_fibs_cn_props <-
#   all_data_meta_cn_fibs %>% 
#   group_by(
#     annotation, CN
#   ) %>%
#   summarise(count = n()) %>%
#   group_by(CN) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = freq, y = CN, fill = annotation)) +
#   geom_col() + 
#   theme_classic()
# 
# ggsave(
#   p_fibs_cn_props, 
#   filename = file.path(cn_out_dir, "Epis CN cell props.pdf"),
#   height = 4, 
#   width  = 8
# )
# 
# 
# p_fibs_cn_props_CN2 <-
#   all_data_meta_cn_fibs %>% 
#   group_by(
#     CN2, CN_fibs
#   ) %>%
#   summarise(count = n()) %>%
#   group_by(CN_fibs) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = freq, y = CN_fibs, fill = CN2)) +
#   geom_col() +
#   theme_classic() +
#   ggsci::scale_fill_igv()
# 
# ggsave(
#   p_fibs_cn_props_CN2, 
#   filename = file.path(cn_out_dir, "Epis CN cell props - all cells CN.pdf"),
#   height = 4, 
#   width  = 8
# )


# p_fibs_cn_props_group <-
#   all_data_meta_cn_fibs %>% 
#   group_by(
#     group, CN_fibs
#   ) %>%
#   summarise(count = n()) %>%
#   group_by(group) %>%
#   mutate(freq = count / sum(count)) %>%
#   mutate(group = factor(group, levels = c(
#     "Normal", "Precancerous", "EarlyStage", "LateStage" 
#   ))) %>% 
#   
#   ggplot(aes(x = group, y = freq, fill = CN_fibs)) +
#   geom_col() +
#   theme_classic() +
#   ggsci::scale_fill_igv() +
#   RotatedAxis()
# 
# ggsave(
#   p_fibs_cn_props_group, 
#   filename = file.path(cn_out_dir, "Epis CN cell props - group.pdf"),
#   height = 6, 
#   width  = 4
# )




# CN annotaton --------------------
# CN k10


cn_annotation <- c(
  "CN1"  = "Myeloid infiltration tumor",
  "CN2"  = "Angio-associated tumor",
  "CN3"  = "Lymphoid infiltration tumor",
  "CN4"  = "Microenvironmrnt-poor tumor",
  "CN5"  = "CAF-associated tumor"
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

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5"
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
#       group.by = "annotation",
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
  
  # dplyr::filter(annotation != "Others") %>%
  
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id", suffix = c("", "_epis")) %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))  
# mutate(CN2_fibs = ifelse(is.na(CN2_fibs), as.character(annotation), as.character(CN2_fibs)))  


# all_data_meta_sp_fibs <-
#   all_data_meta_sp_fibs %>% 
#   left_join(sample_df, by = c(project = "sample"))

# all_data_meta$group <- 
#   factor(
#     all_data_meta$group,
#     levels = c("Normal", "Precancerous", "EarlyStage", "LateStage")
#   )

saveRDS(all_data_meta_sp_fibs, file.path(cn_out_dir, 'all_data_meta_sp_epis2.rds') )

cn_sp_p_ls <-
  cells_position_ls  %>% 
  imap(function(cells_position, project){ # browser()
    
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_position %>%
      right_join(all_data_meta_sp_fibs, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y_centroid, y = x_centroid)) +
      # geom_tile(aes(fill = CN2), height = 10, width = 10, alpha = 0.5)  +
      geom_tile(fill = "gray", height = 10, width = 10, alpha = 0.5)  +
      geom_tile(
        data = ~ dplyr::filter(., !is.na(CN2)),
        aes(fill = CN2), height = 10, width = 10
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
    
    dat <-
      cells_position %>%
      right_join(all_data_meta_sp_fibs, by = "cell_id")
    
    cn_color <- structure(
      names = levels(all_data_meta_sp_fibs$CN2),
      ggsci::pal_igv()(5)
    ) %>%
      c(
        Others = "gray90"
      )
    
    p_ls <-
      levels(all_data_meta_sp_fibs$CN2) %>% # .[1] %>%
      set_names(., .) %>%
      map(function(cn){
        dat2 <- dat
        
        dat2$CN2 <- as.character(dat2$CN2)
        dat2$CN2[dat2$CN2 != cn] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = y_centroid, y = x_centroid)) +
          # geom_tile(aes(fill = CN2), height = 10, width = 10)  +
          geom_tile(fill = "gray", height = 10, width = 10, alpha = 0.5)  +
          geom_tile(
            data = ~ dplyr::filter(., !is.na(CN2)),
            aes(fill = CN2), height = 10, width = 10
          )  +
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
  dplyr::filter(!is.na(CN2)) %>%
  group_by(sample, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df %>%
  mutate(CN2 = factor(CN2, levels = c("Microenvironmrnt-poor tumor", "Angio-associated tumor", "Myeloid infiltration tumor", "Lymphoid infiltration tumor", "CAF-associated tumor"))) %>% 
  ggplot(aes(x = sample, y = freq, fill = CN2)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     sample = all_cn_count_df$sample %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )
# p_cn_props_group <-
#   all_cn_count_df %>%
#   # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'sample', 'group')) %>% 
#   mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>%
#   mutate(group = factor(group, levels = c("Normal", "Precancerous", "EarlyStage", "LateStage"))) %>% 
#   
#   ggplot(aes(x = CN2, y = freq, fill = group)) +
#   geom_boxplot(outliers = FALSE) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) + 
#   ggsci::scale_fill_aaas() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) +
#   labs(
#     x = NULL
#   )

pdf(
  file.path(cn_out_dir, "CN_prop_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
# p_cn_props_group %>% print
dev.off()

# # cell props - each CN ------
# p_cell_props_cn_ls <-
#   all_data_meta_sp_fibs %>% 
#   dplyr::filter(!is.na(CN2)) %>%
#   group_by(CN2) %>%
#   nest() %>%
#   deframe() %>% # .[1] %>%
#   imap(function(dat, name){
#     dat2 <-
#       dat %>%
#       group_by(sample, annotation) %>%
#       summarise(count = n()) %>%
#       group_by(sample) %>%
#       mutate(freq = count/sum(count))
#     dat2 %>%
#       ggplot(aes(x = sample, y = freq, fill = annotation)) +
#       geom_col() +
#       cowplot::theme_cowplot() +
#       scale_y_continuous(labels = scales::percent) +
#       ggsci::scale_fill_igv() +
#       RotatedAxis() +
#       labs(
#         x = NULL,
#         y = NULL,
#         title = name
#       )
#   })
# 
# pdf(
#   file.path(cn_out_dir, "cell props barplot - CN.pdf"),
#   height = 6,
#   width  = 10
# )
# p_cell_props_cn_ls %>% walk(print)
# dev.off()


# de analysis of sub cluster of CN_fibs -------
tmp <- all_data_meta_sp_fibs %>% column_to_rownames('cell_id')
all_data_full2 <- AddMetaData(all_data_full2, tmp['CN2'])
# cell_id_cn_fibs <- all_data_meta_sp_fibs$cell_id[!is.na(all_data_meta_sp_fibs$CN2)] 
# all_data_full2_cn_fibs <- all_data_full2[, cell_id_cn_fibs]
all_data_full2_cn_fibs <- all_data_full2[, !is.na(all_data_full2$CN2)]

all_data_full2_cn_fibs$CN2 %>% table()
all_data_full2_cn_fibs$annotation %>% table()  
?FindAllMarkers  
saveRDS(all_data_full2_cn_fibs, 'CN_analysis_epi/all_data_CN_epis.rds')


all_markers_at1_df <-
  all_data_full2_cn_fibs %>%
  # .[, .$annotation != "doublet"] %>%
  JoinLayers() %>%
  {Idents(.) <- "CN2"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0
  )

all_markers_at1_df %>% head
# p_val avg_log2FC pct.1 pct.2 p_val_adj                    cluster     gene
# TG           0  -2.429483 0.137 0.480         0 Myeloid infiltration tumor       TG
# EPCAM        0  -1.991188 0.201 0.535         0 Myeloid infiltration tumor    EPCAM
# BMP8A        0  -3.210238 0.053 0.328         0 Myeloid infiltration tumor    BMP8A
# SH3BGRL2     0  -2.483685 0.073 0.336         0 Myeloid infiltration tumor SH3BGRL2
# PAX8         0  -1.803535 0.197 0.459         0 Myeloid infiltration tumor     PAX8
# FMOD         0  -2.888344 0.060 0.317         0 Myeloid infiltration tumor     FMOD
write.csv(all_markers_at1_df, file.path(cn_out_dir, 'DEGs_CN2_epi_all.csv'), row.names = F)


# scatter
p_volcanoplot_Tumor_ls <- 
  list("Tumor" = all_markers_at1_df) %>%
  # all_markers_at1_df %>% 
  # group_by(cluster) %>% 
  # nest() %>% 
  # deframe() %>% 
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
  file.path(cn_out_dir, "de scatter plot - CN.pdf"), 
  height = 5, 
  width  = 16
)
p_volcanoplot_Tumor_ls %>%  walk(print)
dev.off()



# enrichment analysis of sub cluster of Tumor -------

all_geneid_AT1_ls <-
  list("Tumor" = all_markers_at1_df) %>% 
  # all_markers_at1_df %>% 
  # group_by(cluster) %>% 
  # nest() %>% 
  # deframe() %>% 
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
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_AT1_ls <- all_geneid_AT1_ls$Tumor

all_enrich_AT1_ls <- 
  all_geneid_AT1_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 0.05, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2 # 0.2
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
  file.path(cn_out_dir, "de enrichment analysis dotplot - CN.pdf"), 
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
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$SYMBOL %>%
      str_c(collapse = "/")
  })

all_enrich_AT1_df %>%
  openxlsx::write.xlsx(
    file.path(cn_out_dir, "enrichment analysis table of CN.xlsx"),
    asTable = TRUE
  )




# Hallmark 50 --------------
all_data <- readRDS("~/R_OUTPUT/all_data_4slides_sample_filtered_reUMAP.rds")
all_data_meta_sp_epis2 <- readRDS("~/CN_analysis_epi/all_data_meta_sp_epis2.rds")
all_data_meta_sp_epis2 %>% head()
all_data_meta <- all_data@meta.data %>% rownames_to_column('cell_id')
all_data_meta %>% head()
tmp <- 
  all_data_meta %>% 
  left_join(all_data_meta_sp_epis2[c('cell_id', 'CN2')], by = 'cell_id') %>% 
  column_to_rownames('cell_id')
all_data <- AddMetaData(all_data, tmp['CN2'])
all_data@meta.data %>% head()


CN2_tumor <- all_data[, !is.na(all_data$CN2)]
CN2_tumor$CN2 %>% table()

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

sc_CN2_tumor_score <-
  JoinLayers(CN2_tumor) %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

sc_CN2_tumor_score %>% head

# msigdb anno
h_score_anno_CN2_tumor <- 
  sc_CN2_tumor_score[[]] %>%
  dplyr::select(CN2, Cluster1:Cluster50) %>%
  group_by(CN2) %>%
  summarise(across(everything(), mean))

p_hallmark_CN2_tumor_1 <-
  pheatmap::pheatmap(
    h_score_anno_CN2_tumor %>% column_to_rownames("CN2") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    cluster_cols = F, 
    silent = TRUE
  )

# export 
pdf(
  file.path(cn_out_dir, "CN_epi_hallmark_score_by_CN.pdf"), 
  height = 8, 
  width  = 4.5
)
p_hallmark_CN2_tumor_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()
