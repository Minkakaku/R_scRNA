library(Seurat)
library(tidyverse)

source("spatial analysis.R")

cn_out_dir <- "CN_analysis"
dir.create(cn_out_dir)

grp_cols_sp <-c(
  Ctrl   = '#1F77B4FF', # Z 
  CD     = "#D62728FF", # Y 
  'CD-F' = "#9467BDFF"  # X 
)


# all_data <- readRDS("all_data.rds")
# all_data <- readRDS("all_data_subcluster.rds")
all_data <- readRDS("all_data_singler.rds")
 

# Idents(all_data) <- "annotation_full"
Idents(all_data) <- "annotation3"

# # group
# all_data$grp <- all_data$group
# all_data$group <- fct_recode(
#   all_data$grp, 
#   Ctrl   = 'Z',
#   CD     = 'Y',
#   'CD-F' = 'X'
# )


cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()



# implement CN ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  # all_data[[c("project", "group", "annotation_full2")]] %>% 
  # dplyr::rename(annotation_full = annotation_full2 ) %>% 

  all_data[[]] %>%
  mutate(annotation4 = ifelse(
    as.character(annotation_full) == "Fibroblast",
    as.character(annotation3),
    as.character(annotation_full)
  )) %>%
  
  dplyr::mutate(annotation_full = annotation4 ) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )

res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)


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

# CN annotaton --------------------
# CN k10


# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

cn_annotation <- c(
  "CN1"  = "CN1",
  "CN2"  = "CN2",
  "CN3"  = "CN3",
  "CN4"  = "CN4",
  "CN5"  = "CN5",
  "CN6"  = "CN6",
  "CN7"  = "CN7",
  "CN8"  = "CN8",
  "CN9"  = "CN9",
  "CN10" = "CN10"
)


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


# CN spatial plot -----------------
cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

saveRDS(all_data_meta, file.path(cn_out_dir, 'all_data_meta_cn.rds') )

cn_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]

    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN2, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })

# CN spatial plot2 -----------------
cn_sp_out_dir2 <- file.path(cn_out_dir, 'cn spatial plot2') 
dir.create(cn_sp_out_dir2)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(cn_out_dir, 'all_data_meta_cn.rds') )

all_data_meta <- readRDS( file.path(cn_out_dir, 'all_data_meta_cn.rds') )
cn_sp_p_ls2 <-
  cells_boundary_ls %>% # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    cn_color <- structure(
        names = paste0("CN", 1:10),
        ggsci::pal_igv()(10)
      ) %>%
      c(
        Others = "gray90"
      )
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p_ls <-
      paste0("CN", 1:10) %>% # .[1] %>%
      set_names(., .) %>%
      map(function(cn){
        dat2 <- dat
        
        dat2$CN2 <- as.character(dat2$CN2)
        dat2$CN2[dat2$CN2 != cn] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = CN2, group = cell_id))  +
          scale_y_reverse() +
          coord_fixed(ratio = 1/image_wh_r) +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cn_color) +
          theme_void() +
          labs(fill = NULL, title = paste(project, group_info, cn, sep = "_") )
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
  all_data_meta %>% 
  group_by(project, group, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df %>%
  
  ggplot(aes(x = project, y = freq, fill = CN2)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn_props_group <-
  all_cn_count_df %>%
  
  ggplot(aes(x = CN2, y = freq, fill = group)) +
  geom_boxplot() +
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



# Proximity analysis -------
# 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    all_data_meta <- 
      all_data[[]] %>%
      mutate(annotation4 = ifelse(
        as.character(annotation_full) == "Fibroblast",
        as.character(annotation3),
        as.character(annotation_full)
      )) %>%
      dplyr::select(-annotation_full) %>%
      dplyr::rename(annotation_full = annotation4)
      
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_full, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_full = factor(annotation_full)) %>%
      dplyr::rename(Classification1 = annotation_full,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


# 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

## 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

# saveRDS(res_pa, file = file.path(cn_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by group and celltype
pa_groups <-
  with(
    all_data_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

# pa_groups2 <-
#   all_data_meta_pos_subset %>%
#   left_join(all_data_meta[c('cell_id', 'CN2')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN2, sep = '..')
#   )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
# res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# 
# res_pa_df2  <-
#   res_pa_df2 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "Group", "project", "CN2"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
    # dplyr::filter(Target == "Epithelial cells") %>%
    ggplot(aes(x = Others, y = cell_props, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.colour = NA) +
    # geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y", ncol = 2)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(cn_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()

# 
# pdf(
#   file.path(cn_out_dir, "proximity analysis - CN.pdf"),
#   height = 50,
#   width  = 25
# )
# p_res_pa_ls %>% print()
# dev.off()

## 2.2 proximity analysis 2 -------
##total data  clean
all_data_meta_pa <-
  all_data[[]] %>%
  mutate(annotation4 = ifelse(
    as.character(annotation_full) == "Fibroblast",
    as.character(annotation3),
    as.character(annotation_full)
  )) %>%
  dplyr::select(-annotation_full) %>%
  dplyr::rename(annotation_full = annotation4) %>%
  mutate(annotation_final = annotation_full %>%
           as.character() %>%
           factor(., levels = unique(.))
  ) %>%
  mutate(
    Group = group
  )



all_data_meta_pos <-
  cells_position_ls  %>%
  imap(function(cells_position, slide){# browser()
    
    if('cell_id' %in% colnames(all_data_meta_pa) ){
      all_data_meta_pa <- all_data_meta_pa %>%
        remove_rownames() %>%
        column_to_rownames("cell_id")
    }
    
    proj_meta_sub <-
      all_data_meta_pa %>%
      dplyr::select(annotation_final, Group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    if(! 'cell_id' %in% colnames(proj_meta_sub) ){
      proj_meta_sub <- rownames_to_column(proj_meta_sub, "cell_id")
    }
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = annotation_final,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")


res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

# saveRDS(res_pa_ls, file.path(cn_out_dir, "res_pa_ls.rds"))


# epithelial proximity
# proximity_cells <- c(
#   "THBS1 Fibroblast",
#   "TNXB Fibroblast",
#   "COL5A1 Fibroblast",
#   "EEF1G Fibroblast",
#   "ADAMDEC1 Fibroblast"
# )
proximity_cells <- c(
  'Immediate early-response (IER) fibroblasts',
  'S3 fibroblasts',
  'S1 fibroblasts',
  'Peri-cryptal S2 fibroblasts' 
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- all_data_meta_pos$Classification1 == cell
    
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- all_data_meta_pos[cell_idx, ]
    
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    
    
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('Group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


# proximity_all_data_ls$Fibroblast$res_pa
# # res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") 


p_proximity_cell_count_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = Group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })



# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
    #   
    #   ggplot(aes(x = distance, fill = Group)) + #, color = Group)) +
    #   geom_density(alpha = 0.6) +
    #   # geom_smooth(aes(y = count), method = "loess") +
    #   facet_wrap(~ Cell, scales = "free_y") +
    #   theme_bw() +
    #   scale_fill_manual(values = grp_cols_sp) +
    #   labs(
    #     x = 'Distance(um)',
    #     y = 'Median count'
    #   )
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, Group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = Group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })




# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - barplot.pdf"), width = 14, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 

p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%

      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(Group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = Group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - radarplot.pdf"), width = 10, height = 8)
p_proximity_radar_ls %>% walk(print)
dev.off()



# NN analysis ----------------

all_data_meta_pos_subset
mpp <- 1 # micron per pixel, xenium 1


all_nn_df_ls <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>%
  nest() %>%
  deframe() %>%
  map(function(dat){ # browser()
      cn_nn <- RANN::nn2(
        dat[c("X", "Y")], 
        k = 10
      )
      
      nn_df <- data.frame(
        target_idx = cn_nn$nn.idx[, 1],
        nn_idx     = cn_nn$nn.idx[, 2],
        nn_dist    = cn_nn$nn.dists[, 2] * mpp
      )
      
      
      nn_df$cell_id      <- dat$cell_id[nn_df$target_idx]
      nn_df$target_cells <- dat$Classification1[nn_df$target_idx]
      nn_df$nn_cells     <- dat$Classification1[nn_df$nn_idx]
      
      nn_df
  })


# add metadata
all_nn_df_full_ls <-
  all_nn_df_ls %>%
  map(function(df){
    df %>% 
      left_join(
        all_data_meta_pos_subset[c("cell_id", "group")], 
        by = "cell_id"
      )
  })

all_nn_df_full <- bind_rows(all_nn_df_full_ls)



# plot
p_nn_density_ls <-
  all_nn_df_full %>%
  group_by(target_cells) %>%
  nest() %>%
  deframe() %>% .[1] %>%
  imap(function(dat, cell){
    
    dat %>%
      ggplot(aes(x = nn_dist, y = group, fill = group)) +
      # geom_density(alpha = 0.3) +
      ggridges::geom_density_ridges(alpha = 0.3) +
      scale_x_continuous(transform = "log2") +
      ggpubr::theme_pubclean() +
      facet_wrap(~ nn_cells) +
      labs(
        title = cell
      )
  })


# export
