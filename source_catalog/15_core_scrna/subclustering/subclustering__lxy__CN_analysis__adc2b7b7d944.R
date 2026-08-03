library(Seurat)
library(tidyverse)
source('CN_result_10/spatial analysis.R')

all_data <- readRDS("/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2.rds")

all_data$annotation2 %>% table()
all_data$group %>% table()

# 改组名
group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
group_levels <- c("Normal", "Mild", "Severe")
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


all_data@meta.data <- all_data@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map),
    Group = factor(Group, levels = group_levels)
  )

table(all_data$Group)

#all_data$treat_result %>% table()

Idents(all_data) <- "annotation2"

cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

cells_position_df <- cells_position_df %>%
  mutate(core_id = str_remove(cell_id, "^TMA_")) # 移除开头 TMA_

## implement CN ---------
# options(warn = 0)

# Classification1 must be a factor
# add a column: obj_id
cn_out_dir <- 'CN_result_10'
dir.create(cn_out_dir)

all_data$core_id <- all_data$cell_id2
head(all_data)


all_data_meta_for_cn_analysis <-  
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  # 改用 core_id 内连接
  left_join(cells_position_df[c('core_id', 'x_centroid', 'y_centroid')], by = "core_id") %>%
  
  dplyr::rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation2 = factor(annotation2),
    obj_id = cell_id
  )


head(all_data_meta_for_cn_analysis)

res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation2",
  k.nn           = 25,
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
  file.path(cn_out_dir, "CN_heatmap_CN10.pdf"),
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
  "CN5"  = "CN5",
  "CN6"  = "CN6",
  "CN7"  = "CN7",
  "CN8"  = "CN8",
  "CN9"  = "CN9",
  "CN10"  = "CN10"
)

cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
cells_boundary_ls <- map(cells_boundary_ls, ~mutate(.x, core_id = str_remove(cell_id, "^TMA_")))
#cells_boundary_df <- cells_boundary_ls %>% bind_rows()
#cells_boundary_df <- cells_boundary_df %>%
#  mutate(core_id = str_remove(cell_id, "^TMA_")) # 移除开头 TMA_
head(cells_boundary_ls)


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

saveRDS(all_data_meta, file.path(cn_out_dir, 'all_data_meta_cn.rds'))
all_data_meta$sample %>% table()


head(all_data_meta)
cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "Group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "core_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN2, group = core_id))  +
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
      file.path(cn_out_dir, str_glue("CN spatial plot {name}.pdf")),
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
    
    anno_data   <- all_data_meta[c("sample", "Group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "core_id")
    
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
          geom_polygon(aes(fill = CN2, group = core_id))  +
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
  group_by(sample, Group, CN2) %>%
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



library(ggpubr)
library(rstatix)
# 1. 配色和主题（和你之前保持一致）
my_theme <- theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# 2. 画基础箱线图（和示例p3_box完全一样）
p3_box <- ggboxplot(
  data = all_cn_count_df,
  x = "CN2",          # 对应示例的Gene
  y = "freq",         # 对应示例的Expression
  fill = "Group"
) +
  ggsci::scale_fill_aaas() +
  my_theme +
  labs(x = NULL) +
  coord_cartesian(ylim = c(0, 0.55)) # 给显著性留空间，和你之前一致

# 3. 计算每个CN内三组两两比较的p值（wilcox.test，和示例逻辑一致）
stat.test_raw <- all_cn_count_df %>%
  group_by(CN2) %>%
  wilcox_test(freq ~ Group) %>%
  add_significance("p") %>%  # 使用原始p值（p列）
  add_xy_position(x = "CN2", dodge = 0.8,
                  step.increase = 0.05)

# 绘图
p_cn_props_group <- p3_box +
  stat_pvalue_manual(
    stat.test_raw,
    label = "p",  # 显示原始p值
    hide.ns = FALSE,
    tip.length = 0.01,
    label.size = 3,  # 可能需要调小一点以便显示数字
    bracket.size = 0.5
  ) +
  theme(
    plot.title = element_text(size = 20, color = "black", hjust = 0.5)
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
  group_by(sample, Group, CN2, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

CN_ls <- all_ct_count_df$CN2 %>% levels()
names(CN_ls) <- CN_ls



p_ls <- 
  CN_ls %>%
  map(function(CN_v){ #browser()
    p_ct_props <-
      all_data_meta_cn %>% 
      dplyr::filter(CN2 == CN_v) %>% 
      group_by(Group, annotation2) %>%
      summarise(count = n()) %>% 
      mutate(freq = count/sum(count)) %>% 
      
      ggplot(aes(x = Group, y = freq, fill = annotation2)) +
      geom_col() +
      #ggsci::scale_fill_aaas() + 
      ggsci::scale_fill_igv() + 
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      labs(
        x = NULL
      )
    
    barplot_order <- sort(unique(p_ct_props$data$annotation2))
    
    ct_sample_freq <- all_data_meta_cn %>%
      dplyr::filter(CN2 == CN_v) %>%
      group_by(sample, Group,annotation2) %>%
      summarise(count = n()) %>%
      mutate(freq = count / sum(count)) %>%
      mutate(annotation2 = factor(annotation2, levels = barplot_order)) 
    
    p4_box <- ggboxplot(
      data = ct_sample_freq %>% 
        mutate(annotation2 = factor(annotation2, levels = barplot_order)),
      x = "annotation2",          # 对应示例的Gene
      y = "freq",         # 对应示例的Expression
      fill = "Group"
    )+
      ggsci::scale_fill_aaas() +
      my_theme+
      labs(x = NULL) +
      coord_cartesian(ylim = c(0, max(ct_sample_freq$freq, na.rm = TRUE) * 1.2))
    
    # 计算每个 annotation2 内三组两两比较的 p 值
    stat.test_ct <- ct_sample_freq %>%
      group_by(annotation2) %>%
      filter(n_distinct(Group) >= 2) %>% 
      wilcox_test(freq ~ Group) %>%
      add_significance("p") %>%
      add_xy_position(x = "annotation2", dodge = 0.8,
                      step.increase = 0.05)
    
    # 绘图（带上显著性）
    p_ct_props_group <- p4_box +
      stat_pvalue_manual(
        stat.test_ct,
        label = "p",
        hide.ns = FALSE,
        tip.length = 0.01,
        label.size = 2.5,
        bracket.size = 0.5
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


# 转h5ad
# 后续在jupyter中运行NMF:/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_NMF/nmf.ipynb
# 在服务器运行cn tissue：/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_tissue
# 读取文件
# 主Seurat对象（原始注释数据）
all_data <- readRDS("/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2.rds")
all_data$annotation2 %>% table()
all_data$group %>% table()

# 改组名
group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
group_levels <- c("Normal", "Mild", "Severe")
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


all_data@meta.data <- all_data@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map),
    Group = factor(Group, levels = group_levels)
  )

table(all_data$Group)

# CN分析结果数据框
cn_df <- readRDS("CN_result_10/all_data_meta_cn.rds")


# 把Seurat的细胞名转为cell_id列，用于匹配
seurat_meta <- all_data@meta.data %>%
  rownames_to_column("cell_id")  # 细胞名转为cell_id

# 按cell_id拼接（左连接，保留所有Seurat细胞，无CN信息则填充NA）
merge_meta <- left_join(
  x = seurat_meta,
  y = cn_df,
  by = "cell_id"   # 两边细胞ID列名必须完全一致
) %>%
  column_to_rownames("cell_id")  # 还原行名为细胞名

# 追加元数据到Seurat对象
# 只追加【CN新增的列】，避免重复覆盖原有列
new_cn_cols <- setdiff(colnames(merge_meta), colnames(all_data@meta.data))
if (length(new_cn_cols) > 0) {
  all_data <- AddMetaData(all_data, metadata = merge_meta[, new_cn_cols, drop = FALSE])
  cat("成功追加CN相关列：", paste(new_cn_cols, collapse = ", "), "\n")
} else {
  cat("未检测到新列，列名已存在，跳过追加\n")
}

#保存整合后的Seurat RDS
saveRDS(all_data, file = "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_result_10/alldata_CN_10.rds")

table(all_data$annotation2)
immune_cell_list <- c(
  "B cell",
  "CD4+",
  "CD8+Effector",
  "CD8+TEM",
  "cDC1",
  "cDC2",
  "LAMP3_DC",
  "pDC",
  "M1",
  "M2",
  "Mast cell",
  "Neutrophil",  
  "Plasma",
  "Proliferating T cell",
  "Treg"
)

data_immune <- subset(all_data, subset = annotation2 %in% immune_cell_list)


table(data_immune$annotation2)

saveRDS(data_immune, file = "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_result_10/immune_cells_CN_10.rds")


# 绑定Python环境
use_python('/home/szt/miniforge3/envs/rpy/bin/python', required = TRUE)
ad <- import('anndata')
## sceasy转换格式，conda env tools
all_data <- FindVariableFeatures(all_data)
VariableFeaturePlot(all_data)
all_data2 <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

all_data@assays$RNA <- as(all_data@assays$RNA, "Assay")


convertFormat(
  obj = all_data,
  from = "seurat",
  to = "anndata",
  outFile = "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_result_10/alldata_CN_10.h5ad"
)

## sceasy转换格式，conda env tools
data_immune <- FindVariableFeatures(data_immune)
VariableFeaturePlot(data_immune)
data_immune2 <- ScaleData(data_immune, vars.to.regress = 'nCount_RNA')

data_immune@assays$RNA <- as(data_immune@assays$RNA, "Assay")


convertFormat(
  obj = data_immune,
  from = "seurat",
  to = "anndata",
  outFile = "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/CN_result_10/immune_cells_CN_10.h5ad"
)
