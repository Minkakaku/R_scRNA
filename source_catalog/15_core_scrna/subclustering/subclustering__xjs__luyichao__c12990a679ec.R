library(Seurat)
library(tidyverse)

out_dir <- 'luyichao'
dir.create(out_dir)

all_data <- readRDS('R_OUTPUT/all_data_4slides_sample_filtered.rds')
all_data$sample %>% table()

all_data2 <- all_data[, all_data$sample %in% c(
  'A31', 'A32', 
  'P54', 'A131', 'A71'
)]
all_data2$sample %>% table()

all_data2$group <- str_replace_all(all_data2$sample, c(
  'A31' = 'with TLS', 
  'A32' = 'with TLS', 
  'P54' = 'without TLS', 
  'A131' = 'without TLS', 
  'A71' = 'without TLS'
))
all_data2$group %>% table()


# cell prop by group ----------
p_col_proj <- 
  all_data2[[]] %>% 
  ungroup() %>% 
  group_by(group, annotation) %>% 
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(out_dir, "cell_freq_group_annotation.pdf"),
  height = 6,
  width  = 6
)


saveRDS(all_data2, file.path(out_dir, 'all_data2_TLS.rds'))




# inferCNV ---------------------------------------------------------------------
Epithelial <- all_data2[, all_data2$annotation == 'Epithelial']
Epithelial <- FindVariableFeatures(Epithelial)

## set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

library(infercnv)
infercnv_out_dir = "luyichao/inferCNV_episub" 
dir.create(infercnv_out_dir) 


## 映射-------------------------------
library(tibble)
library(data.table)

ref_cells <- c('B cell', 'Plasma', 'T cell') # , 'Endothelial', 'Fibroblast', 'Macrophage', 'Pericyte'
tumor_cells <- "Epithelial"

meta_data_infer_ls <- 
  all_data2[[c('annotation', 'sample')]] %>% 
  dplyr::filter(annotation %in% c(ref_cells, tumor_cells)) %>%
  rownames_to_column('cell_id') %>% 
  group_by(sample) %>%
  nest() %>% 
  deframe()


# cell is too much
# less than 50000
set.seed(1234)
meta_data_infer_subset_ls <-
  meta_data_infer_ls %>% 
  imap(function(meta_data_infer, sample){
    meta_data_infer2 <- 
      meta_data_infer %>% 
      group_by(annotation) %>%
      slice_sample(prop = 1)
  })


# export as file
meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){
    if(!file.exists(file.path(infercnv_out_dir, sample))){
      dir.create(file.path(infercnv_out_dir, sample))
    }
    data.table::fwrite(
      meta_data_infer_subset, 
      file = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  })



infercnv_obj_ls <- 
  meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){ # browser()
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
      #raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:1000], meta_data_infer_subset$cell_id],
      annotations_file  = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      gene_order_file   = "luyichao/gene_pos47.txt",
      delim             = "\t",
      min_max_counts_per_cell = c(50,Inf),
      ref_group_names = ref_cells
    )
    infercnv_obj
  })



infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. > 0.01) } # 3280 <--- 
infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. < 0.01) } # 1477

undebug(infercnv:::.single_tumor_leiden_subclustering)

options(scipen = 100)
infercnv_obj_ls2 <- 
  infercnv_obj_ls %>% 
  imap(function(infercnv_obj, sample){
    infercnv_obj2 <- infercnv::run(
      infercnv_obj,
      cutoff              =    0.01,
      out_dir    =   file.path(infercnv_out_dir, sample),
      cluster_by_groups   =    TRUE, 
      # k_nn                =    20,
      # leiden_method       =    "simple", # PCA 
      plot_steps          =    FALSE,
      denoise             =    TRUE,
      HMM                 =    FALSE,
      no_prelim_plot      =    TRUE,
      png_res             =    300,
      # useRaster   = FALSE,
      output_format = 'pdf',
      no_plot = T #  figure margin too large 
    )
    infercnv_obj2
  })

infercnv_obj_ls2 %>% 
  imap(function(infercnv_obj2, sample){
    infercnv::plot_cnv(
      infercnv_obj2,
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv',
      output_format     =  "pdf",
      png_res = 300,
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  mean(infercnv_obj@expr.data), # 1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })



infercnv_obj2_default_median_filtering_ls <-
  infercnv_obj_ls2 %>%
  map(function(infercnv_obj2){
    infercnv::apply_median_filtering(infercnv_obj2)
  })


infercnv_obj2_default_median_filtering_ls %>% 
  imap(function(infercnv_obj2_default_median_filtering, sample){
    infercnv::plot_cnv(
      infercnv_obj2_default_median_filtering, 
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv_median_filtering',
      output_format="pdf",
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })

saveRDS(infercnv_obj2_default_median_filtering_ls, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering_ls.rds'))

saveRDS(infercnv_obj_ls2, file.path(infercnv_out_dir, 'infercnv_obj_ls2.rds'))
saveRDS(infercnv_obj_ls, file.path(infercnv_out_dir, 'infercnv_obj_ls.rds'))


par(mafr
    )

## plot -------------------------------------------------------------------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_ls <- 
  cnv_score_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_ls <- 
  temp_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation, log(CNV_score)))+
      geom_violin(aes(fill= annotation), color="NA", drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_ls2 <-
  temp_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation) %>%
      
      ggplot(aes(x = annotation, y = log(CNV_score), fill = annotation)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot.pdf"),
      height = 6,
      width  = 8
    )
  })


cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_noQuantile_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_noQuantile_ls <- 
  cnv_score_noQuantile_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(., as.data.frame(meta_data_infer_subset_ls[[sample]]), 'cell_id') -> temp
    temp
  })


plot_noQuantile_ls <- 
  temp_noQuantile_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation, log(CNV_score)))+
      geom_violin(aes(fill= annotation), color="NA", drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_noQuantile_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_noQuantile.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_noQuantile_ls2 <-
  temp_noQuantile_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation) %>%
      
      ggplot(aes(x = annotation, y = log(CNV_score), fill = annotation)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_noQuantile_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot_noQuantile.pdf"),
      height = 6,
      width  = 8
    )
  })




# Epi DEGs --------------
Epithelial$sample %>% table()
Epithelial$group %>% table()
table(Epithelial$sample, Epithelial$group)

TLS_epi_markers <- 
  FindMarkers(
    JoinLayers(Epithelial), 
    group.by = 'group', 
    ident.1 = 'without TLS', 
    ident.2 = 'with TLS', 
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0
  ) %>% 
  rownames_to_column('gene')
# TLS_epi_markers <- rownames_to_column(TLS_epi_markers, 'gene')
TLS_epi_markers %>% head

write.csv(TLS_epi_markers, file.path(out_dir, 'TLS_epi_markers.csv'), row.names = F)




# all_data DEGs -------------
all_data2$sample %>% table()
all_data2$group %>% table()
table(all_data2$sample, all_data2$group)

TLS_epi_markers <- 
  FindMarkers(
    JoinLayers(all_data2), 
    group.by = 'group', 
    ident.1 = 'without TLS', 
    ident.2 = 'with TLS', 
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0
  ) %>% 
  rownames_to_column('gene')
# TLS_epi_markers <- rownames_to_column(TLS_epi_markers, 'gene')
TLS_epi_markers %>% head

write.csv(TLS_epi_markers, file.path(out_dir, 'TLS_all_data_markers.csv'), row.names = F)




# immune DEGs -------------
all_data2$annotation %>% table()
immune <- all_data2[, all_data2$annotation %in% c('B cell', 'Macrophage', 'Plasma', 'T cell')]
immune$annotation %>% table()
table(immune$sample, immune$group)

TLS_immune_markers <- 
  FindMarkers(
    JoinLayers(immune), 
    group.by = 'group', 
    ident.1 = 'without TLS', 
    ident.2 = 'with TLS', 
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0
  ) %>% 
  rownames_to_column('gene')

TLS_immune_markers %>% head

write.csv(TLS_immune_markers, file.path(out_dir, 'TLS_immune_markers.csv'), row.names = F)




# T_cell DEGs -------------
all_data2$annotation %>% table()
T_cell <- all_data2[, all_data2$annotation %in% c('T cell')]
T_cell$annotation %>% table()
table(T_cell$sample, T_cell$group)

TLS_T_cell_markers <- 
  FindMarkers(
    JoinLayers(T_cell), 
    group.by = 'group', 
    ident.1 = 'without TLS', 
    ident.2 = 'with TLS', 
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0
  ) %>% 
  rownames_to_column('gene')

TLS_T_cell_markers %>% head

write.csv(TLS_T_cell_markers, file.path(out_dir, 'TLS_T_cell_markers.csv'), row.names = F)




# Ro/e ----------
roe_out_dir <- file.path(out_dir, 'roe_analysis')
dir.create(roe_out_dir)


## sample -----------
meta <- all_data2@meta.data
meta$sample <- factor(meta$sample, levels = c(
  'A31', 'A32', 
  'P54', 'A131', 'A71')
  )
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_sample'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    roe_out_dir_cn <- str_glue("{roe_out_dir}/{name}/") 
    if(!file.exists(roe_out_dir_cn)){
      dir.create(roe_out_dir_cn)
    }
    
    anno_column <- "annotation" # <---------------------
    grp_column  <- "sample"    # <---------------------
    
    all_data_meta_roe <- meta
    
    
    ## Ro/e -------
    # Cell Preference
    library(ComplexHeatmap)
    # 生成细胞数目分布频数统计表observe.data
    observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
    stopifnot(is.matrix(observe.data))
    
    expected <-
      function(x){
        if(!is.matrix(x)){stop("Must be a matrix")}
        rtot <- margin.table(x, 1)
        ctot <- margin.table(x, 2)
        tot  <- margin.table(x)
        outer(rtot, ctot, "*")/tot
      }
    expected.data <- expected(observe.data) # 生成期望频数表
    plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
    # 绘图数据的列顺序
    col.order <- ggplot2:::`%||%`(
      levels(all_data_meta_roe[[grp_column]]),
      sort(unique(all_data_meta_roe[[grp_column]]))
    ) 
    plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
    plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
    # 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max
    
    # 热图颜色设定
    roe_cols <- structure(
      c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
      names = c("1", "1.5", "2", ">3")
    )
    
    cell_fun = function(j, i, x, y, width, height, fill) { # browser()
      
      # print(c(j, i, x, y, width, height, fill))
      
      if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
        grid::grid.text("\u00B1", x = x, y = y)
      } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
        grid::grid.text("+", x = x, y = y)
      } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
        grid::grid.text("++", x = x, y = y)
      } else {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
        grid::grid.text("+++", x = x, y = y)
      }
      
    }
    
    # 生成热图
    roi_hm <-
      ComplexHeatmap::Heatmap(
        plot.data, 
        cell_fun = cell_fun,
        cluster_rows = F, 
        cluster_columns = F,
        width = unit(5, "cm"), 
        height = unit(17, "cm"), 
        show_heatmap_legend = F
      )
    
    # stat comparasion
    chi2.T <- sum((observe.data - expected.data)^2/expected.data)
    p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
    print(p_chisq)
    
    # 生成图例
    lgd <-
      ComplexHeatmap::Legend(
        labels = names(roe_cols), 
        title = expression(R[o/e]), 
        legend_gp = grid::gpar(fill = roe_cols)
      )
    
    pdf(
      file.path(roe_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })


## group -----------
meta <- all_data2@meta.data
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_group'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    roe_out_dir_cn <- str_glue("{roe_out_dir}/{name}/") 
    if(!file.exists(roe_out_dir_cn)){
      dir.create(roe_out_dir_cn)
    }
    
    anno_column <- "annotation" # <---------------------
    grp_column  <- "group"    # <---------------------
    
    all_data_meta_roe <- meta
    
    
    ## Ro/e -------
    # Cell Preference
    library(ComplexHeatmap)
    # 生成细胞数目分布频数统计表observe.data
    observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
    stopifnot(is.matrix(observe.data))
    
    expected <-
      function(x){
        if(!is.matrix(x)){stop("Must be a matrix")}
        rtot <- margin.table(x, 1)
        ctot <- margin.table(x, 2)
        tot  <- margin.table(x)
        outer(rtot, ctot, "*")/tot
      }
    expected.data <- expected(observe.data) # 生成期望频数表
    plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
    # 绘图数据的列顺序
    col.order <- ggplot2:::`%||%`(
      levels(all_data_meta_roe[[grp_column]]),
      sort(unique(all_data_meta_roe[[grp_column]]))
    ) 
    plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
    plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
    # 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max
    
    # 热图颜色设定
    roe_cols <- structure(
      c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
      names = c("1", "1.5", "2", ">3")
    )
    
    cell_fun = function(j, i, x, y, width, height, fill) { # browser()
      
      # print(c(j, i, x, y, width, height, fill))
      
      if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
        grid::grid.text("\u00B1", x = x, y = y)
      } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
        grid::grid.text("+", x = x, y = y)
      } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
        grid::grid.text("++", x = x, y = y)
      } else {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
        grid::grid.text("+++", x = x, y = y)
      }
      
    }
    
    # 生成热图
    roi_hm <-
      ComplexHeatmap::Heatmap(
        plot.data, 
        cell_fun = cell_fun,
        cluster_rows = F, 
        cluster_columns = F,
        width = unit(5, "cm"), 
        height = unit(17, "cm"), 
        show_heatmap_legend = F
      )
    
    # stat comparasion
    chi2.T <- sum((observe.data - expected.data)^2/expected.data)
    p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
    print(p_chisq)
    
    # 生成图例
    lgd <-
      ComplexHeatmap::Legend(
        labels = names(roe_cols), 
        title = expression(R[o/e]), 
        legend_gp = grid::gpar(fill = roe_cols)
      )
    
    pdf(
      file.path(roe_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })




# cell prop by annotation ----------
# 单样本进行细胞比例统计后，按分组取均值画堆叠图，百分比标在条的右侧
plot_data <- 
  all_data2[[]] %>% 
  ungroup() %>% 
  group_by(sample, annotation, group) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(freq = count / sum(count)) %>% 
  ungroup() %>% 
  group_by(group, annotation) %>% 
  summarise(aver_freq = mean(freq), .groups = "drop") %>% 
  group_by(annotation) %>% 
  mutate(global_mean = mean(aver_freq),
         total = sum(aver_freq)) %>% 
  ungroup() %>% 
  mutate(annotation = forcats::fct_reorder(annotation, global_mean, .desc = TRUE))

# 准备标签数据：每个 annotation 下两个组的百分比，用 '/' 拼接
label_data <- 
  plot_data %>%
  group_by(annotation) %>%
  summarise(label = paste(sprintf("%.1f%%", rev(aver_freq) * 100), collapse = "/"),
            total = sum(aver_freq))  # 总长度用于定位标签x坐标

p_col_proj <- 
  ggplot(plot_data, aes(x = aver_freq, y = annotation, fill = group)) + 
  geom_col(position = "stack", width = 0.7) +  # 堆叠条形
  geom_text(data = label_data, 
            aes(x = total + 0.02, y = annotation, label = label), 
            inherit.aes = FALSE,   # 关键：不继承主图的 aes
            hjust = 0, size = 3.5) +
  ggsci::scale_fill_igv() +
  labs(x = "cell proportion (%)", y = NULL) +
  scale_x_continuous(limits = c(0, 0.9)) +    
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export
ggsave(
  p_col_proj,
  filename = file.path(out_dir, "cell_proportion_annotation.pdf"),
  height = 5,
  width  = 8
)




# T subtype -----------
all_data <- readRDS("luyichao/all_data2_TLS.rds")
all_data$sample %>% table()
all_data$annotation %>% table()

Major_marker_ls <- c(
  "Epithelial" = "EPCAM", 
  "Epithelial" = "KRT5", 
  "Epithelial" = "KRT14", 
  "Epithelial" = "KRT20", 
  
  "Endothelial" = "PECAM1", 
  "Endothelial" = "VWF", 
  "Endothelial" = "AQP1", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2", 
  
  "B cell" = "CD19", 
  "B cell" = "CD79A", 
  "B cell" = "MS4A1", 
  "Plasma" = "TNFRSF17", 
  "Plasma" = "JCHAIN", 
  "Plasma" = "IGJ", 
  "Plasma" = "MZB1", 
  "Plasma" = "XBP1", 
  
  "T cell" = "CD3D", 
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "TRAC", 
  "T cell" = "BCL11B", 
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  "Proliferating" = "MKI67", 
  "Proliferating" = "TOP2A", 
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  #  "Pericyte" = "PLVAP", 
  
  "Monocyte" = "CD14", 
  "Monocyte" = "FCGR3A", 
  "Monocyte" = "ITGAM", 
  
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  "Macrophage" = "C1QA", 
  "Macrophage" = "C1QC", 
  "Macrophage" = "CSF1R", 
  "Macrophage" = "IL1B", 
  
  "Dendritic cell" = "CD83", 
  "Dendritic cell" = "LILRA4", 
  "Dendritic cell" = "CLEC9A", 
  "Dendritic cell" = "IDO1", 
  "Dendritic cell" = "CD1C", 
  "Dendritic cell" = "FCER1A", 
  "Dendritic cell" = "CST3", 
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  "Lamp3_DC" = "LAMP3", 
  
  "Netrophil" = "CSF3R", 
  "Netrophil" = "FCGR3B", 
  "Netrophil" = "S100A8", 
  "Netrophil" = "S100A9", 
  "Netrophil" = "G0S2", 
  "Netrophil" = "RETNLG", 
  
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = "HDC", 
  "Mast cell" = "TPSAB1", 
  "Mast cell" = "CTSG", 
  "Mast cell" = "CPA3", 
  
  "Smooth muscle cell" = "PLN", 
  "Smooth muscle cell" = "MYLK", 
  "Smooth muscle cell" = "MYH11", 
  "Smooth muscle cell" = "ACTG2", 
  "Smooth muscle cell" = "TAGLN", 
  
  "Schwann cell" = "MPZ", 
  "Schwann cell" = "SOX10", 
  
  "Neuron" = "TUBB3", 
  "Neuron" = "SNAP25", 
  
  "Mesothelial" = "MSLN", 
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2", 
  
  "Adipocyte" = "PLIN1", 
  "Adipocyte" = "PLIN4", 
  "Adipocyte" = "ADIPOQ", 
  "Adipocyte" = "CIDEC", 
  "Adipocyte" = "AQP7", 
  "Adipocyte" = "LPL"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


All_T <- subset(all_data, annotation == 'T cell')


## PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')
All_T <- RunPCA(All_T)

T_out_dir <- file.path(out_dir, 'T_subtype') 
dir.create(T_out_dir)


### pca inspector ------
pdf(file.path(T_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T <- RunUMAP(All_T, dims = 1:dims,return.model = TRUE)

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


All_T <- All_T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T <- RunUMAP(
  All_T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T  <-
  FindNeighbors(
    All_T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T <-
  FindClusters(
    All_T,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
T_major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  #'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

CD4_major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'T cells' = 'CD3D', 
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'T cells' = 'TRAC',
  'CD4'='CD4',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2'
)


CD8_major_marker_ls <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAYN',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)

# find markers
Idents(All_T) <- "RNA_snn_res.0.4"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(T_markers, file.path(T_out_dir, 'T_markers.csv'))

T_top5_markers <- 
  T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(T_top5_markers) <- sub('.*\\.','',names(T_top5_markers))



# dotplot
pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = Major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = Epi_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = T_top5_markers)+RotatedAxis()
dev.off()


# Idents(All_T) <- "RNA_snn_res.0.8"
# All_T <- FindSubCluster(All_T, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.6)
# Idents(All_T) <- "sub.cluster"
# # find markers
# T_markers_sub <-
#   FindAllMarkers(
#     JoinLayers(All_T),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# 
# T_top5_markers_sub <- 
#   T_markers_sub %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 30) %>%
#   arrange(desc(pct.1)) %>%
#   slice_head(n = 5) %>%
#   nest() %>%
#   deframe() %>% 
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]
# 
# names(T_top5_markers_sub) <- sub('.*\\.','',names(T_top5_markers_sub))
# 
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# 
# pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
# DotPlot(All_T, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_top5_markers_sub)+theme_90angle
# dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "mix", # T, Epi
  '1' = "Macrophage",
  '2' = "CD4 naive T", # 
  '3' = "Fibroblast", # 
  '4' = "CD8 Tem", # 
  '5' = "CD4 Treg", # 
  '6' = "Proliferating T", # 
  '7' = "mix", # 
  '8' = "Plasma", # 
  '9' = "Endothelial", # 
  '10' = "Epithelial", # 
  '11' = "Lamp3 DC", # 
  '12' = "Epithelial" # 
)


Idents(All_T) <- "RNA_snn_res.0.4"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)

pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 7, width = 18)
DimPlot(All_T, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



saveRDS(All_T, file.path(T_out_dir, 'All_T.rds'))


## pure -----------
All_T_pure <- All_T[, All_T$subtype %in% c('CD4 naive T', 'CD8 Tem', 'CD4 Treg', 'Proliferating T')]
All_T_pure <- RunUMAP(All_T_pure, dims = 1:dims, return.model = TRUE)


pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(All_T_pure, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure.pdf")), height = 7, width = 18)
DimPlot(All_T_pure, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()


cellType_marker_ls <- c(
  'T' = 'CD3E',
  
  'CD8' = 'CD8A',
  
  'CD4' = 'CD4',

  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='GZMB',
  
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='STMN1',
  'proliferation'='ANAPC11',
  'proliferation'='TOP2A'
)

pdf(file.path(T_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 20, height = 8)
DotPlot(All_T_pure, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(All_T_pure, file.path(T_out_dir, 'All_T_pure.rds'))
