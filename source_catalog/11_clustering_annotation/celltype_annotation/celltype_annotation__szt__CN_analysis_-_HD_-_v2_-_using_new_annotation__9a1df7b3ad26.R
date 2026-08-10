library(Seurat)
library(tidyverse)

source("spatial analysis.R")

cn_out_dir <- "CN_analysis - v2"
dir.create(cn_out_dir)


all_data_full <- readRDS("all_data_full.rds")

meta_data <- arrow::read_parquet("data/20250703_group_celltype_columns.parquet")

meta_data %>%
  group_by(celltype) %>%
  # summarise(cnv_score = mean(as.numeric(as.character(cnv_score)), na.omit = T)) %>%
  mutate(cnv_score = as.numeric(as.character(cnv_score)) ) %>% 
  ggplot(aes(x = celltype, y = cnv_score)) +
  geom_boxplot(outliers = F) +
  RotatedAxis()


meta_data_celltype <-
  meta_data %>% 
  dplyr::mutate(
    # sample,
    id = `id_sample`,
    # celltype 
  ) %>%
  mutate(sample = str_replace(sample, "LC", "LC_")) %>%
  mutate(id2 = id %>% str_replace("-1_.*$", "-1")) %>%
  mutate(id = str_c(sample, id2, sep = "_")) %>%
  dplyr::select(
    id, celltype,
    group, group_all, group_stage, group_TLS
  )

meta_data_celltype_full <-
  meta_data_celltype %>%
  mutate(group = id %in% colnames(all_data_full)) %>%
  group_by(group) %>%
  summarise(
    count = list(table(celltype))
  )
  
# many cells are droped by QC
# use all_data instead
# all_data_full <- merge(visim_ls[[1]], visim_ls[-1])

all_data_full$celltype <- NULL
all_data_full <- 
  all_data_full %>% 
  AddMetaData(column_to_rownames(meta_data_celltype, "id"))


table(is.na(all_data_full$celltype))



Idents(all_data_full) <- "celltype"


all_cellid_not_others <- colnames(all_data_full)[(all_data_full$celltype != "Others") & (!is.na(all_data_full$celltype))]
cells_position_ls <- 
  all_data_full@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  setNames(., str_remove(names(.), "\\..*$")) %>%
  map(function(x){
    x[x$cell_id %in% all_cellid_not_others, ]
  })
cells_position_df <- cells_position_ls %>% bind_rows() 


# save rds
# saveRDS(all_data_full, 'all_data_full.new.annos.rds')

# implement CN ---------

# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k.nn 计算的最近邻的细胞数，也就是每一个细胞周围多少个细胞参与计算邻域
# k.cluster 返回多少个CN


# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  all_data_full[[]] %>%
  dplyr::filter(celltype != "Others") %>% 
  dplyr::filter(!is.na(celltype)) %>% 
  dplyr::mutate(annotation_full = celltype ) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # important!
    obj_id = cell_id # important!
  )

stopifnot(all(!is.na(all_data_meta_for_cn_analysis$celltype)))

Sys.time()
res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
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

# CN annotaton --------------------
# CN k10


cn_annotation <- c(
  "CN1"  = "CN1_Tumor core",
  "CN2"  = "CN2_Immune + SMC",
  "CN3"  = "CN3_Intestinal metaplasia",
  "CN4"  = "CN4_SMC",
  "CN5"  = "CN5_PMC",
  "CN6"  = "CN6_Tumor boundary",
  "CN7"  = "CN7_Oxyntic gland",
  "CN8"  = "CN8_Immune core",
  "CN9"  = "CN9_Fibroblast",
  "CN10" = "CN10_GMC"
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

# spatial plot -----------------
sp_out_dir <- file.path(cn_out_dir, 'spatial plot') 
dir.create(sp_out_dir)

p_sp_plot_anno <-
  Images(all_data_full) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    if(img_name == "HC") img_name <- "NML"
    
    p <- SpatialDimPlot(
      all_data_full,
      images = img, 
      group.by = "celltype",
      pt.size.factor	= 2
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export

p_sp_plot_anno %>%
  iwalk(function(p, name){
    on.exit(dev.off())
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot annotation - {name}.pdf")),
      height = 8,
      width  = 8
    )
    print(p)
    # dev.off()
  })



# CN spatial plot -----------------
cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- 
  all_data_full[[]] %>%
  
  dplyr::filter(celltype != "Others") %>%
  
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

# all_data_meta <-
#   all_data_meta %>% 
#   left_join(sample_df, by = c(project = "sample"))

all_data_meta$group <- 
  factor(
    all_data_meta$group,
    levels = c("NORMAL", "IM", "LGD", "HGD", "early_cancer", "late_cancer", "TLS", "other")
  )

saveRDS(all_data_meta, file.path(cn_out_dir, 'all_data_meta_cn2.rds') )

cn_sp_p_ls <-
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]

    dat <-
      cells_position %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN2), height = 10, width = 10)  +
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
  all_data_meta %>% 
  group_by(project, group, CN2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
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

  ggplot(aes(x = CN2, y = freq, fill = group)) +
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
  file.path(cn_out_dir, "CN_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
p_cn_props_group %>% print
dev.off()

# cell props - each CN ------
p_cell_props_cn_ls <-
  all_data_meta %>% 
  group_by(CN2) %>%
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

# Proximity analysis -------
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



# 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

all_data_meta <- 
  readRDS(file.path(cn_out_dir, 'all_data_meta_cn2.rds')) %>% 
  dplyr::mutate(annotation_full = celltype )

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    # all_data_meta <- 
    #   all_data_full[[]] %>%
    #   dplyr::filter(celltype != "Others") %>%
    #   dplyr::mutate(annotation_full = celltype )

    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_full, cell_id, group, project) %>%
      # rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_full = factor(annotation_full)) %>%
      dplyr::rename(Classification1 = annotation_full,
                    X = x,
                    Y = y,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

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
    k = 300, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_full@images[[1]]@scale.factors$spot, # !important, 1 for xenium
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

res_pa <- readRDS(file.path(cn_out_dir, 'res_pa.rds'))


# summarized by group and celltype
pa_groups <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>%
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    with(
      dat,
      str_c(Classification1, group, name, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)
  

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
  mutate(Group = factor(Group, levels = c(
    "NORMAL", "IM",  "LGD", "HGD",  "early_cancer", "late_cancer", "TLS", "other"        
  ))) %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA, linewidth = 0.1) +
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
  width  = 30
)
p_res_pa %>% print()
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
# ##total data  clean
# all_data_meta_pa <-
#   all_data[[]] %>%
#   mutate(annotation4 = ifelse(
#     as.character(annotation_full) == "Fibroblast",
#     as.character(celltype),
#     as.character(annotation_full)
#   )) %>%
#   dplyr::select(-annotation_full) %>%
#   dplyr::rename(annotation_full = annotation4) %>%
#   mutate(annotation_final = annotation_full %>%
#            as.character() %>%
#            factor(., levels = unique(.))
#   ) %>%
#   mutate(
#     Group = group
#   )
# 
# 
# 
# all_data_meta_pos <-
#   cells_position_ls  %>%
#   imap(function(cells_position, slide){# browser()
#     
#     if('cell_id' %in% colnames(all_data_meta_pa) ){
#       all_data_meta_pa <- all_data_meta_pa %>%
#         remove_rownames() %>%
#         column_to_rownames("cell_id")
#     }
#     
#     proj_meta_sub <-
#       all_data_meta_pa %>%
#       dplyr::select(annotation_final, Group, project) %>%
#       rownames_to_column("cell_id") %>%
#       dplyr::filter(project == slide)
#     
#     if(! 'cell_id' %in% colnames(proj_meta_sub) ){
#       proj_meta_sub <- rownames_to_column(proj_meta_sub, "cell_id")
#     }
#     
#     cells_position[, 1:3] %>%
#       inner_join(proj_meta_sub, by = "cell_id") %>%
#       dplyr::rename(Classification1 = annotation_final,
#                     X = x_centroid,
#                     Y = y_centroid,
#                     Slide_ID = project)
#   }) %>%
#   enframe() %>%
#   unnest(cols = value) %>%
#   dplyr::select(-name) %>%
#   column_to_rownames("cell_id")

mpp <- 8 / all_data@images$TSA_233841.008um@scale.factors$spot

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
    mpp = mpp, # !important, 1 for xenium
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
    mpp = mpp, # !important, 1 for xenium
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
    mpp = mpp, # !important, 1 for xenium
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
    mpp = mpp, # !important, 1 for xenium
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
    mpp = mpp, # !important, 1 for xenium
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
proximity_cells <- all_data_meta_pos$Classification1 %>% unique() %>% as.character()

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){ # browser()
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
    
    # median 
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
    
    # turn to freq
    res_pa_sub_diff_ls2 <-
      res_pa_sub_diff_ls %>% 
      map(function(dat){
        col_sum <- colSums(dat)
        dat <- sweep(dat, 2, STATS = col_sum, FUN = "/")
        if(any(col_sum == 0)){
          dat[, col_sum == 0] <- 0
        }
        dat
      })
    
    
    # median freq
    res_pa_sub_freq <-
      res_pa_sub_diff_ls2 %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_freq <-
      res_pa_sub_freq %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count,
      res_pa2  = res_pa_sub_freq
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
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot(outliers = FALSE, linewidth = 0.1) +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })
p_proximity_cell_freq_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa2 %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot(outliers = FALSE, linewidth = 0.1) +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median freq',
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


res_pa_freq_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa2 %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <- res_pa_count_raw
      # 1:nrow(res_pa_count_raw) %>%  
      # map(function(x){
      #   count <- res_pa_count_raw$count[x] 
      #   if(count < 1) return(NULL) # <- count是median count，有可能是0.5
      #   xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
      #   xx
      # }) %>%
      # purrr::discard(is.null) %>%
      # bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    res_pa_count_tidy$freq  <- res_pa_count_tidy$count
    
    res_pa_count_tidy
  })

# # only selecet two macrohpages
# p_proximity_cell_count_distribution_ls2 <-
#   proximity_all_data_ls %>%
#   imap(function(dat_ls, name){
#     
#     all_data_group_sub <- dat_ls$all_data
#     
#     p <-
#       dat_ls$res_pa %>%
#       pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#       left_join(all_data_group_sub, by = "Slide_ID")  %>%
#       # mutate(distance = as.numeric(distance)) %>%
#       mutate(distance = factor(distance, levels = c("10", "50", "100", "150", "200"))) %>%
#       
#       dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
#       
#       ggplot(aes(x = Cell, y = count, fill = distance)) +
#       # geom_violin(trim = FALSE) +
#       geom_boxplot() +
#       # geom_point() +
#       facet_wrap(~ group, scales = "free_y") +
#       theme_bw() +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#       # scale_fill_manual(values = grp_cols_sp) +
#       labs(
#         # x = 'Distance(um)',
#         y = 'Median count',
#         title = name
#       )
#     
#     p 
#   })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){ # browser()

    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_x_continuous(breaks = c(10, 50, 100, 150, 200)) +
      # scale_color_manual(values = grp_cols_sp) +
      ggsci::scale_color_igv(guide = guide_legend(override.aes = list(linewidth = 3, fill = NA))) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })

p_proximity_freq_ls <-
  res_pa_freq_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_x_continuous(breaks = c(10, 50, 100, 150, 200)) +
      # scale_color_manual(values = grp_cols_sp) +
      ggsci::scale_color_igv(guide = guide_legend(override.aes = list(linewidth = 3, fill = NA))) +
      labs(
        x = 'Distance(um)',
        y = 'Median freq',
        title = name
      )
    
    p
  })


# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - barplot.pdf"), width = 30, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(cn_out_dir, "proximity analysis - freq - barplot.pdf"), width = 30, height = 8)
p_proximity_cell_freq_distribution_ls %>% walk(print)
dev.off()

# pdf(file.path(cn_out_dir, "proximity analysis - barplot2.pdf"), width = 14, height = 6)
# p_proximity_cell_count_distribution_ls2 %>% walk(print)
# dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_ls %>% walk(print)
dev.off()

pdf(file.path(cn_out_dir, "proximity analysis - freq - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_freq_ls %>% walk(print)
dev.off()

# radar plot 

p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%

      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
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
      # scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - radarplot.pdf"), width = 10, height = 8)
p_proximity_radar_ls %>% walk(print)
dev.off()


# Nearest neighbor analysis -------
# 由于最近邻只需要最近的那个细胞的信息
# 因此不同于proximity analysis，有以下变化：
# classification调为NULL
# searchtype 邻近度使用radius
# k不需要太大，10即可，最近邻分析只需要第二个细胞的信息
# radius适当调大100um，防止部分细胞比较离散，距离其他细胞较远
# returntype与cell_enough_ratio_cutoff也没有意义了

# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification NULL即可，返回每一个query细胞的细胞编号、距离、query和total细胞的id
# k 计算的最近邻的细胞数
# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325


res_nn_ls <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = NULL, # 若为NULL: 提前返回window的细胞近邻的list
    k = 10, # 'auto',
    searchtype = "radius",
    radius = 100, 
    mpp = 8 / all_data_full@images[[1]]@scale.factors$spot, # !important, 1 for xenium
    cell_enough_ratio_cutoff = 20
  )

# res_nn_ls <- readRDS("data/res_nn_ls.rds")

res_nn_df <-
  res_nn_ls %>%
  map(function(x){ # browser()
    cell_idx  <- x$nn.idx[, 2]
    cell_idx[cell_idx == 0] <- NA # some cells have no nn cells
    cell_id   <- x$total_cell_id[cell_idx]
    cell_dists <- x$nn.dists[, 2]
    
    all_data_used <- all_data_meta_pos_subset %>% column_to_rownames("cell_id")
    cell_anno <- all_data_used[cell_id, 'Classification1']
    
    res <-
      data.frame(
        Others = cell_anno,
        distance = cell_dists
      ) %>%
      cbind(all_data_used[x$query_cell_id, c('Classification1', 'group', 'Slide_ID')]) %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(Target = Classification1) 
    
    res
  }) %>%
  bind_rows()

if(any(is.na(res_nn_df$Others) )){
  message(str_glue("{sum(is.na(res_nn_df$Others))} cells have no nearest neighbor cell"))
  res_nn_df <- drop_na(res_nn_df)
}

# plot
# p_nn_ls <-
#   res_nn_df %>%
#   group_by(Target) %>%
#   nest() %>%
#   deframe() %>%
#   imap(function(x, target){ # browser()
#     x %>%
#       ggplot(aes(x = distance, color = group, fill = group)) + 
#       geom_density(alpha = 0.5) +
#       facet_wrap(~ Others) +
#       ggpubr::theme_pubclean() +
#       scale_y_continuous(trans = "log1p") +
#       labs(
#         title = str_glue("Target cell: {target}"),
#         x = "distance(um)"
#       )
#     
#   })

p_nn_ls2 <-
  res_nn_df %>%
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){ # browser()
    x %>%
      ggplot(aes(x = distance, y = group, fill = group)) + 
      # geom_density(alpha = 0.5) +
      ggridges::geom_density_ridges() +
      facet_wrap(~ Others) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        x = "distance(um)"
      )
    
  })
# export plot
pdf(
  file.path(cn_out_dir, "nearest neighbor analysis.pdf"),
  height = 14,
  width  = 18
)
p_nn_ls2 %>% print()
dev.off()




