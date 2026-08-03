library(Seurat)
library(tidyverse)
library(patchwork)
source('CN_mTEC//spatial analysis.R')


cn_out_dir<- "CN_mTEC" 
dir.create(cn_out_dir)

mTEC <- readRDS("mTEC_subtype_3/mTEC_clean.rds") #重新读取mTEC，之前做完拟时序分析筛掉了一部分细胞
mTEC$subtype %>% table()

all_data$mTEC_3_Plasma_T_annotation %>% table()
# 新建 annotation 列，并按要求批量修改注释
all_data$annotation <- dplyr::recode(
  all_data$mTEC_3_Plasma_T_annotation,
  
  # 改为 TEC
  `CCL19+ mTEClow` = "TEC",
  `Cycling TEC` = "TEC",
  `TEC precursor` = "TEC",
  
  # 改为 T cell
  CD4_Treg = "T cell",
  CD8_TEM = "T cell",
  `Gamma-delta_T` = "T cell",
  Naive_T = "T cell",
  Proliferative_T = "T cell",
  TCM = "T cell",
  Tfh = "T cell",
  
  `Mito-high` = "Plasma",
  `Proliferative plasmablast` = "Plasma",
  
  # 改为 mix
  cTEC_Naive_T = "mix"
)

# 查看结果（验证）
table(all_data$annotation)


# 1. 读取你的CSV文件（细胞ID列表）
nc_cells <- read.csv("mTEC_monocle/NC.csv", stringsAsFactors = FALSE)

# 2. 提取CSV里的所有细胞ID
csv_cell_ids <- nc_cells$Cell.ID

# 3. 提取当前Seurat对象里的所有细胞ID
seurat_cell_ids <- rownames(all_data@meta.data)  # 或 mTEC$cell_id2 都可以

seurat_clean <- gsub("-H[0-9]+$", "", seurat_cell_ids)

# 4. 分组：在CSV里 = NC，不在 = Hyperplasia
all_data$groupfinal <- ifelse(
  seurat_clean %in% csv_cell_ids,
  "NC",
  "Hyperplasia"
)

# 5. 查看结果
table(all_data$groupfinal)

table(all_data$groupfinal, all_data$sample)

#以上都在准备all_data文件，下面的代码annotation_full对应我的annotation，group对应我的groupfinal

#接着准备cells_position_ls.rds文件



names(all_data@images)

















#all_data <- readRDS("mTEC_subtype_3/all_data_mTEC_3_Plasma_T_annotation.rds")
all_data$annotation %>% table()
all_data$groupfinal %>% table()

Idents(all_data) <- "annotation"

# cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
# cells_position_df <- cells_position_ls %>% bind_rows()

cells_position_ls <-
all_data %>%
  FetchData(vars = c('SPATIAL_1', 'SPATIAL_2', 'sample')) %>%
  rownames_to_column('cell_id') %>%
  dplyr::rename(x_centroid = SPATIAL_1, y_centroid = SPATIAL_2) %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

cells_position_df <- cells_position_ls %>% bind_rows()

## implement CN ---------
# options(warn = 0)

# Classification1 must be a factor
# add a column: obj_id
cn_out_dir <- 'CN_result'
dir.create(cn_out_dir)

all_data_meta_for_cn_analysis <-  
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  dplyr::rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation = factor(annotation), # import!
    obj_id = cell_id # import!
  )

res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn, file.path(cn_out_dir, "res_cn_25_10_cellType.rds"))


## heatmap ----------------------
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




## CN spatial plot -----------------
cn_annotation <- c(
  "CN1"  = "CN1",
  "CN2"  = "CN2",
  "CN3"  = "CN3",
  "CN4"  = "CN4",
  "CN5"  = "CN5"
)

cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

saveRDS(all_data_meta, file.path(cn_out_dir, 'all_data_meta_cn.rds'))

cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "treat_result")] %>% distinct()
    group_info  <- anno_data$treat_result[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN2, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })

cn_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })




## CN spatial plot2 -----------------
cn_sp_out_dir2 <- file.path(cn_out_dir, 'cn spatial plot2') 
dir.create(cn_sp_out_dir2)

all_data_meta <- readRDS( file.path(cn_out_dir, 'all_data_meta_cn.rds') )
all_data_meta_cn <- readRDS( file.path(cn_out_dir, 'all_data_meta_cn.rds'))
cn_sp_p_ls2 <- 
  cells_boundary_ls %>% # .[1] %>%
  imap(function(cells_boundary, sample){ # browser()
    # cn_color <- structure(
    #   names = paste0("CN", 1:10),
    #   ggsci::pal_igv()(10)
    # ) %>%
    #   c(
    #     Others = "gray90"
    #   )
    # image_wh_r <- image_wh_r_ls[[sample]]
    
    anno_data   <- all_data_meta[c("sample", "treat_result")] %>% distinct()
    group_info  <- anno_data$treat_result[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    cn_color <- structure(
      names = dat$CN2 %>% levels(),
      ggsci::pal_igv()(10)
    ) %>%
      c(
        Others = "gray90"
      )
    
    p_ls <-
      dat$CN2 %>% levels() %>% # .[1] %>%
      set_names(., .) %>%
      map(function(cn){
        dat2 <- dat
        
        dat2$CN2 <- as.character(dat2$CN2)
        dat2$CN2[dat2$CN2 != cn] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = CN2, group = cell_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cn_color) +
          theme_void() +
          labs(fill = NULL, title = paste(sample, group_info, cn, sep = "_") )
        p
      })
    
    p_ls
  })


# export image
cn_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
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
  group_by(sample, groupfinal, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN2)) +
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
  
  ggplot(aes(x = CN2, y = freq, fill = groupfinal)) +
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


all_ct_count_df <-
  all_data_meta %>% 
  group_by(sample, groupfinal, CN2, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

CN_ls <- all_ct_count_df$CN2 %>% levels()
names(CN_ls) <- CN_ls


p_ls <- 
  CN_ls %>%
  map(function(CN_v){ # browser()
    p_ct_props <-
      all_data_meta_cn %>% 
      dplyr::filter(CN2 == CN_v) %>% 
      group_by(groupfinal, annotation) %>%
      summarise(count = n()) %>% 
      mutate(freq = count/sum(count)) %>% 
      
      ggplot(aes(x = groupfinal, y = freq, fill = annotation)) +
      geom_col() +
      #ggsci::scale_fill_aaas() + 
      ggsci::scale_fill_igv() + 
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      labs(
        x = NULL
      )
    
    
    p_ct_props_group <-
      all_data_meta_cn %>% 
      dplyr::filter(CN2 == CN_v) %>% 
      group_by(sample, groupfinal, annotation) %>%
      summarise(count = n()) %>%
      mutate(freq = count/sum(count)) %>% 
      
      ggplot(aes(x = annotation, y = freq, fill = groupfinal)) +
      geom_boxplot() +
      ggsci::scale_fill_aaas() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      labs(
        x = NULL
      )
    p_ct_props / p_ct_props_group
    # patchwork::wrap_plots(
    #   p_ct_props,
    #   p_ct_props_group,
    #   widths = c(3, 7)
    # )
  })

p_ls %>% 
  iwalk(function(p, annotation){
    ggsave(
      p,
      filename = file.path(cn_out_dir, str_glue("cell_pos/CN_boxplot_and_barplot_{annotation}.pdf")),
      height = 12,
      width  = 15
    )
  })

pdf(
  file.path(cn_out_dir, str_glue("CN_boxplot_and_barplot_2.pdf")),
  height = 8,
  width  = 7
)
p_ls %>% 
  iwalk(function(p, annotation){
    print(p$patches$plots[[1]] + ggtitle(annotation))
  })
dev.off()
