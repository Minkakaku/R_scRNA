#1、亚类注释-----------
##Myeloid 亚群空间解析，将 Myeloid 拆解为 M1（促炎）和 M2（修复）亚分类注释。------
# 1、dat_analysis_h5ad_transfer_v3.R的  4、髓系亚类注释
##Fibroblast 的修复性亚群鉴定，区分"修复性成纤维细胞"与"纤维化性肌成纤维细胞"。------
# 1、dat_analysis_h5ad_transfer_v3.R的 5、成纤维、肌成纤维亚类注释  


#2、细胞通讯-----
all_data2$annotation_Lgr5 %>% table()

cellchat_data <- all_data2
cellchat_data <- JoinLayers(cellchat_data)

library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 50 * 1024^3) 

cellchat_data <- NormalizeData(cellchat_data)

cellchat_data <- cellchat_data[ ,!cellchat_data$annotation_Lgr5 %in% c("delete")]

cellchat_data$Celltype <- cellchat_data$annotation_Lgr5
cellchat_data$Celltype %>% unique()

cellchat_data$Celltype <- gsub(" ","_",cellchat_data$Celltype)

cellchat_data$group%>% table()

saveRDS(cellchat_data, file = "/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/Cellchat_all/cellchat_data.rds")

## 投bash------
#/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/Cellchat_all




#2、MSC-Myeloid-ISC 三元空间互作网络-----
# ①BM_MSCs          Myeloid
# ②Myeloid        Lgr5 Intestinal Stem Cell
# ③BM_MSCs          Lgr5 Intestinal Stem Cell

all_data2$annotation_Lgr5 %>% table()

cellchat_data <- all_data2
cellchat_data <- JoinLayers(cellchat_data)

library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 50 * 1024^3) 

cellchat_data <- NormalizeData(cellchat_data)

cellchat_3types_data <- cellchat_data[ ,cellchat_data$annotation_Lgr5 %in% c("BM_MSCs","Myeloid","Lgr5 Intestinal Stem Cell")]

cellchat_3types_data$Celltype <- cellchat_3types_data$annotation_Lgr5
cellchat_3types_data$Celltype %>% unique()

cellchat_3types_data$Celltype <- gsub(" ","_",cellchat_3types_data$Celltype)

cellchat_3types_data$group%>% table()

saveRDS(cellchat_3types_data, file = "/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/Cellchat_3types/cellchat_3types_data.rds")

## 投bash------
#/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/Cellchat_3types








#3、MSC拟时序，组间趋势------
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)
library(tibble)

#BM-MSC内
BMMSC <- readRDS(file.path(BMMSC_out_dir,'BMMSC.rds'))
table(BMMSC@meta.data$RNA_snn_res.0.2)

BMMSC_monocle_outdir <- "result/BMMSC_monocle_result/"
dir.create(BMMSC_monocle_outdir)

unique(BMMSC@meta.data$RNA_snn_res.0.2)
Idents(BMMSC)
saveRDS(BMMSC,file.path(BMMSC_monocle_outdir,"BMMSC_mono.rds"))

## 投shell------
#/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/BMMSC_monocle_result

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

cds_BMMSC <- readRDS("/mnt/NexData_file/home/hyj/shell/FWSP20254307_yurongbing/BMMSC_monocle_result/cds_VariableFeatures-2.rds")

pData(cds_BMMSC)$RNA_snn_res.0.2 <- factor(pData(cds_BMMSC)$RNA_snn_res.0.2#, levels = c(paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9))
)

pdf(file.path(BMMSC_monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds_BMMSC,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds_BMMSC,color_by = "RNA_snn_res.0.2", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds_BMMSC, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds_BMMSC, color_by = "group", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(BMMSC_monocle_outdir, "cell trajectory - - split by RNA_snn_res.0.2.pdf"), width = 10, height = 12)
plot_cell_trajectory(cds_BMMSC, color_by = "RNA_snn_res.0.2") + facet_wrap("~RNA_snn_res.0.2", nrow = 5) + theme(legend.position = "right") + scale_color_manual(values = plan) 
dev.off()

pdf(file.path(BMMSC_monocle_outdir, "cell trajectory - - split by state.pdf"), width = 10, height = 12)
plot_cell_trajectory(cds_BMMSC, color_by = "State") + facet_wrap("~State", nrow = 3) + scale_color_manual(values = plan)
dev.off()

pdf(file.path(BMMSC_monocle_outdir, "cell trajectory - cluster - split by group.pdf"), width = 12, height = 10)
local({
  #pData(cds_BMMSC)$group <- pData(cds_BMMSC)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds_BMMSC, color_by = "RNA_snn_res.0.2", size = 1, show_backbone = TRUE) + 
    scale_color_manual(values = plan) + 
    facet_wrap(~ group) + 
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(BMMSC_monocle_outdir, "cell trajectory - tree.pdf"), width = 30, height = 30)
plot_complex_cell_trajectory(cds_BMMSC, x = 1, y = 2, color_by = "RNA_snn_res.0.2") + theme(legend.position = "right") + 
  scale_color_manual(values = plan) +
  facet_wrap(~ RNA_snn_res.0.2)
dev.off()

pdf(file.path(BMMSC_monocle_outdir, "cell trajectory  group - tree.pdf"), width = 30, height = 30)
plot_complex_cell_trajectory(cds_BMMSC, x = 1, y = 2, color_by = "RNA_snn_res.0.2") + theme(legend.position = "right") + 
  scale_color_manual(values = plan) +
  facet_wrap(~ group)
dev.off()



# cell barplot 
p_col_cell_state <- 
  cds_BMMSC %>% 
  pData() %>% 
  dplyr::select(RNA_snn_res.0.2, State) %>% 
  group_by(RNA_snn_res.0.2, State) %>% 
  mutate(RNA_snn_res.0.2 = as.character(RNA_snn_res.0.2)) %>% 
  mutate(RNA_snn_res.0.2 = factor(RNA_snn_res.0.2#, levels = c(
                          #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                          # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = RNA_snn_res.0.2, y = freq)) + 
  geom_col(aes(fill = State)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = NULL, y = NULL) + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(BMMSC_monocle_outdir, "cell stat barplot.pdf"), height = 12, width = 8)


# cell barplot 
p_col_cell_state <- 
  cds_BMMSC %>% 
  pData() %>% 
  dplyr::select(RNA_snn_res.0.2, group) %>% 
  group_by(RNA_snn_res.0.2, group) %>% 
  mutate(RNA_snn_res.0.2 = as.character(RNA_snn_res.0.2)) %>% 
  mutate(RNA_snn_res.0.2 = factor(RNA_snn_res.0.2#, levels = c(
                                  #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                                  # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = RNA_snn_res.0.2, y = freq)) + 
  geom_col(aes(fill = group)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "cluster", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(BMMSC_monocle_outdir, "cell stat cluster group barplot.pdf"), height = 12, width = 8)


p_col_cell_state <- 
  cds_BMMSC %>% 
  pData() %>% 
  dplyr::select(State, group) %>% 
  group_by(State, group) %>% 
  mutate(State = as.character(State)) %>% 
  mutate(State = factor(State#, levels = c(
                                  #   paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
                                  # )
  )
  ) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) + 
  geom_col(aes(fill = group)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = "State", y = "") + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(BMMSC_monocle_outdir, "cell stat State group barplot.pdf"), height = 12, width = 8)



# 沿时间轴的细胞密度图 -----------
colour = c('#C6DADF', '#E26EB8', '#918FBE', '#A2C5A6', '#EED784', '#D4B796', 
           '#F5BDC8', '#F0BDD8', '#77A1EC', '#F8D39C', '#F17591', '#CCC9E6')

df <- pData(cds_BMMSC)[, c("Pseudotime", "RNA_snn_res.0.2")]
p_density <- ggplot(df, aes(Pseudotime, colour = RNA_snn_res.0.2, fill = RNA_snn_res.0.2)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(BMMSC_monocle_outdir, "cell density cluster plot.pdf"), height = 6, width = 12)


df <- pData(cds_BMMSC)[, c("Pseudotime", "group")]
p_density <- ggplot(df, aes(Pseudotime, colour = group, fill = group)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(BMMSC_monocle_outdir, "cell density group plot.pdf"), height = 6, width = 12)

#4、邻域分析------
source("spatial analysis.R")


#使用annotation_Lgr5列做加入M1M2，剔除delete
all_data2$annotation_Lgr5 %>% table()

Myeloids <- readRDS("result/Myeloids_subtype/Myeloids.rds")
Myeloids$subtype %>% table()

## 映射-------------------------------
library(tibble)
library(data.table)
sub <- Myeloids[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    all_data2[['annotation_Lgr5']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_Lgr5_M1M2 = ifelse((as.character(annotation_Lgr5) %in% c('Myeloid')),
                      as.character(subtype), 
                      as.character(annotation_Lgr5)))%>%
  dplyr::select(cells, annotation_Lgr5_M1M2) %>%
  column_to_rownames('cells')

meta.data.t$annotation_Lgr5_M1M2 %>% table()

all_data2 <- AddMetaData(all_data2,meta.data.t['annotation_Lgr5_M1M2'])

all_data2@meta.data %>% head()
table(all_data2$annotation_Lgr5_M1M2,all_data2$annotation_Lgr5)

saveRDS(all_data2,"result/all_data_Lgr5_M1M2.rds")

all_data2 <- readRDS("result/all_data_Lgr5_M1M2.rds")

##空间meta构建--------
cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

# options(warn = 0)
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  all_data2[[c("sample", "group", "annotation_Lgr5_M1M2")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation_Lgr5_M1M2 = factor(annotation_Lgr5_M1M2), # import!
    obj_id = cell_id # import!
  )

write.csv(all_data_submeta_for_cn_analysis,file.path(CN_out_dir,"all_data_submeta_for_cn_analysis.csv"))

## CN analysis ---------
CN_out_dir <- "result/CN analysis"
dir.create(CN_out_dir)

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation_Lgr5_M1M2",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))


### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  all_data2[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds'))
All_type_meta_cn10$project <- All_type_meta_cn10$sample

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### CN k15 ---------
res_cn_15 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation_Lgr5_M1M2",
  k.nn           = 10,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 9,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()


### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data2$cell_id_raw <- all_data2$cell_id
all_data2$cell_id %>% head()
all_data2$cell_id <- NULL

All_type_meta_cn15 <- 
  all_data2[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(CN_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn15[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn15_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
All_type_cn15_count_df <-
  All_type_meta_cn15 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn15_props_group <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()


#5、邻近度分析------
cells_boundary_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("/mnt/NexData_file/home/hyj/spatial analysis.R")

proximity_out_dir <- "result/proximity_out_dir_all"
dir.create(proximity_out_dir)

all_data2@meta.data %>% head()
all_data2@meta.data$cell_id %>% head()
rownames(all_data2@meta.data) %>% head()
#xenium_RCTD_noEPI@meta.data$cell_id <- paste0("1_",xenium_RCTD_noEPI@meta.data$cell_id1)

table(all_data2$group,all_data2$sample)

#saveRDS(all_data2,"RCTD/result_noEPI/all_data2.rds")

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
## Proximity analysis -------
## 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- all_data2[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_Lgr5_M1M2, sample, "group") %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_Lgr5_M1M2 = factor(annotation_Lgr5_M1M2)) %>%
      dplyr::rename(Classification1 = annotation_Lgr5_M1M2,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos.csv"))
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)

## 2. proximity analysis -------
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

### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
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
#pa_groups <-
#  with(
#    all_data_meta_pos_subset,
#    str_c(Classification1, Group, Slide_ID, sep = '..')
#  )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)



# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

res_pa_df$Target %>% unique()
p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Others != "delete") %>% 
  #dplyr::filter(Target == "BM_MSCs") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  #geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y")+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis.pdf"),
  height = 14,
  width  = 18
)
p_res_pa %>% print()
dev.off()

####以各亚群为核心循环-------
cell_list <- res_pa_df$Target %>% unique()
# res_pa_df$Group <- factor(res_pa_df$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    #dplyr::filter(Others %in% cell_list) %>% 
    dplyr::filter(Others != "other") %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
    )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
    height = 6,
    width  = 18
  )
  p_res_pa %>% print()
  dev.off()
})



### 2.2 proximity analysis - series analysis -----
init_cells <- c(
  '10'  = 750,
  '50'  = 800,
  '100' = 1000,
  '150' = 1200,
  '200' = 1500
)

res_pa_ls <-
  c(10, 50, 100, 150, 200) %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos,
        query = all_data_meta_pos,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = k, # 'auto',
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })

# res_pa_ls <- readRDS("data/res_pa_ls.rds")

# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  c(10, 50, 100, 150, 200) %>% 
  as.character() %>%
  rev() %>%
  {matrix(c(.[-length(.)], .[-1]), ncol = 2)} %>%
  apply(1, c, simplify = FALSE)

for(i in diff_idx_ls){
  first  <- i[1]
  second <- i[2]
  res_pa_ls2[[ first ]] <- res_pa_ls2[[ first ]] - res_pa_ls2[[ second ]] 
}

# count to freq
res_pa_ls2 <- 
  res_pa_ls2 %>%
  map(function(x){# browser()
    x_colsum <- colSums(x)
    xx <- sweep(x, 2, x_colsum, FUN = "/")
    if(any(x_colsum == 0)){
      xx[, x_colsum == 0] <- 0
    }
    xx
  })

# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy <- 
  res_pa_ls2 %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    # 分组内取中位数
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "median")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()


# plot 
# only plot "Epithelial cells" ....
res_pa_count_tidy_epi <-
  res_pa_count_tidy %>%
  mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
  dplyr::filter(Target == "BM_MSCs") %>% 
  dplyr::filter(Others != "delete")  

res_pa_count_tidy_epi_median <-
  res_pa_count_tidy_epi %>%
  group_by(Others, Group, distance) %>% 
  summarise(cell_props = median(cell_props))

p_proximity_series_epi <-
  res_pa_count_tidy_epi %>%
  ggplot(aes(x = distance, y = cell_props, color = Group)) +
  geom_boxplot( ) +
  ggsci::scale_color_aaas() +
  # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
  geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
  geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
  facet_wrap(~ Others, scales = 'free_y') +
  ggpubr::theme_pubclean()


# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()


####循环各亚群-------
cell_list <- res_pa_count_tidy$Target %>% unique()
#res_pa_count_tidy$Group <- factor(res_pa_count_tidy$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  
  res_pa_count_tidy_epi <-
    res_pa_count_tidy %>%
    mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
    dplyr::filter(Target == cell) %>% 
    dplyr::filter(Others != "other") 
  
  res_pa_count_tidy_epi_median <-
    res_pa_count_tidy_epi %>%
    group_by(Others, Group, distance) %>% 
    summarise(cell_props = median(cell_props))
  
  p_proximity_series_epi <-
    res_pa_count_tidy_epi %>%
    ggplot(aes(x = distance, y = cell_props, color = Group)) +
    geom_boxplot( ) +
    ggsci::scale_color_aaas() +
    # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
    geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
    geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
    facet_wrap(~ Others, scales = 'free_y') +
    ggpubr::theme_pubclean()
  
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis_series distance.pdf")),
    height = 14,
    width  = 18
  )
  p_proximity_series_epi %>% print()
  dev.off()
})


