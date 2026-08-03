library(Seurat)
library(tidyverse)

#setwd('~/')
cn_out_dir <- 'CN_analysis_epi/kMeans_cluster3'
dir.create(cn_out_dir)
# all_data <- readRDS("R_OUTPUT/data_for_analysis/all_data_filtered2.rds")
# all_data$annotation %>% table()
# all_data$group %>% table()

all_data_meta <- 
  readRDS('CN_analysis_epi/all_data_meta_sp_epis2.rds') %>% 
  filter(!is.na(CN2))
all_data_meta$CN2 <- factor(all_data_meta$CN2, levels = c(
  "Microenvironmrnt-poor tumor", 
  "Angio-associated tumor", 
  "Myeloid infiltration tumor", 
  "Lymphoid infiltration tumor", 
  "CAF-associated tumor"
))

all_cn_count_df <- 
  all_data_meta %>% 
  group_by(sample, group, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))
all_cn_count_df %>% head()


all_cn_count_matrix <- 
  all_cn_count_df %>% 
  as.data.frame() %>% 
  dplyr::select(sample, freq, CN2) %>% 
  pivot_wider(names_from = CN2, values_from = freq) %>% 
  column_to_rownames('sample')

na_idx <- is.na(all_cn_count_matrix)
all_cn_count_matrix[na_idx] <- 0
all_cn_count_matrix


set.seed(1234)
cn_cluster <- kmeans(all_cn_count_matrix, centers = 3, iter.max = 100)
cn_cluster$cluster




# plot ------------
## 热图 ------------
message('heatmap')
library(pheatmap)
library(RColorBrewer)

# 准备热图数据
heatmap_data <- all_cn_count_matrix[order(cn_cluster$cluster), ]
cluster_order <- cn_cluster$cluster[order(cn_cluster$cluster)]

# 创建注释
annotation <- data.frame(Cluster = factor(cluster_order))
rownames(annotation) <- rownames(heatmap_data)

# 绘制热图
pdf(file.path(cn_out_dir, 'Kmeans Clustering Heatmap.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = F,
         annotation_row = annotation,
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "K-means Clustering Heatmap",
         fontsize_row = 8,
         show_rownames = FALSE)
dev.off()




## 平行坐标图 ----------
message('parallel')
library(GGally)

# 添加聚类标签
data_with_clusters <- all_cn_count_matrix
data_with_clusters$cluster <- factor(cn_cluster$cluster)

# 平行坐标图
ggparcoord(data_with_clusters, 
           columns = 1:5,
           groupColumn = "cluster",
           scale = "globalminmax") +
  labs(title = "Parallel Coordinates Plot",
       x = "Features", 
       y = "Standardized Values") +
  theme_minimal()
ggsave(file.path(cn_out_dir, 'Parallel Coordinates Plot.pdf'), height = 6, width  = 11)




## 簇间距离矩阵可视化 --------------
# 计算簇中心距离
center_dist <- dist(cn_cluster$centers)
center_dist_matrix <- as.matrix(center_dist)

# 绘制簇间距离热图
pdf(file.path(cn_out_dir, 'Distance Between Cluster Centers.pdf'), height = 6, width  = 11)
pheatmap(center_dist_matrix,
         main = "Distance Between Cluster Centers",
         color = colorRampPalette(c("white", "red"))(50),
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()




## PCA ------------
# 如果特征维度较高，先用PCA降维
pca_result <- prcomp(all_cn_count_matrix, scale. = TRUE)

# PCA结果可视化
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$cluster <- factor(cn_cluster$cluster)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(title = "K-means Clustering on PCA-reduced Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
ggsave(file.path(cn_out_dir, 'Kmeans Clustering on PCA-reduced Data.pdf'), height = 6, width  = 11)




# arrange samples --------------
all_data_meta %>% head()
sample_cluster_dist <- cn_cluster$cluster %>% as.character() %>% str_c('cluster_', .)
names(sample_cluster_dist) <- names(cn_cluster$cluster)
saveRDS(sample_cluster_dist, file.path(cn_out_dir, 'sample_cluster_dist_kMeans.rds'))


all_cn_count_df <- 
  all_data_meta %>% 
  group_by(sample, group, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# 3
all_cn_count_df$sample <- factor(all_cn_count_df$sample, levels = c(
  sample_cluster_dist[sample_cluster_dist == 'cluster_1'] %>% names(), 
  sample_cluster_dist[sample_cluster_dist == 'cluster_2'] %>% names(), 
  sample_cluster_dist[sample_cluster_dist == 'cluster_3'] %>% names()
))
# "A11"          "A12"          "A121"         "A131"         "A141"         "A21"          "A31"          "A32"  "A41"          "A51"          "A61"          "ATC12_be_T_1" "ATC14_be_T"   "ATC18_be_T_1"
# "A101"         "A111"         "A112"         "A132"         "A71"          "ATC11_be_T"   "ATC16_be_T_1" "ATC16_be_T_2" "ATC19_be_T"   "P11"          "P12"          "P21"          "P31"          "P32"          "P41"          "P51"  "P52"          "P53"          "P54"
# "ATC10_be_T_1" "ATC10_be_T_3" "ATC12_be_T_2" "ATC17_be_T"   "ATC18_be_T_2"




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

ggsave(file.path(cn_out_dir, 'cell_prop_CN_kmeans_clustering.pdf'), height = 6, width  = 11)




## consensus clustering --------------
results <- readRDS('Epi_consensus_clustering_4slides_sample_filtered/result.rds')
results[[3]][["consensusClass"]]
all_cn_count_df$sample <- as.character(all_cn_count_df$sample)
all_cn_count_df$sample <- factor(all_cn_count_df$sample, levels = c(
  results[[3]][["consensusClass"]][results[[3]][["consensusClass"]] == '1'] %>% names(), 
  results[[3]][["consensusClass"]][results[[3]][["consensusClass"]] == '2'] %>% names(), 
  results[[3]][["consensusClass"]][results[[3]][["consensusClass"]] == '3'] %>% names()
))
# "A121" "A101" "A21" "A141" "ATC18_be_T_2" "ATC10_be_T_1" "ATC11_be_T" "ATC14_be_T" "ATC10_be_T_3" "ATC12_be_T_2" "ATC12_be_T_1"
# "A112" "A71" "P31" "P32" "P21" "A111" "A61" "A32" "A41" "P12" "P11" "A12" "A11"  "A51" "A131" "A31" "ATC17_be_T" "ATC18_be_T_1"
# "P54" "P53" "A132" "P41" "P52" "P51" "ATC16_be_T_1" "ATC16_be_T_2""ATC19_be_T" 

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

ggsave(file.path(cn_out_dir, 'cell_prop_CN_consensus_clustering.pdf'), height = 6, width  = 11)




# DEGs -------------
## CN kMeans ----------
all_data_CN_epis <- readRDS('CN_analysis_epi/all_data_CN_epis.rds')
names(sample_cluster_dist) <- str_c(names(sample_cluster_dist), '$')
all_data_CN_epis$group_CN_kMeans <- str_replace_all(all_data_CN_epis$sample, c(
  sample_cluster_dist
))
table(all_data_CN_epis$group_CN_kMeans, all_data_CN_epis$sample)

all_markers_at1_df <-
  all_data_CN_epis %>%
  JoinLayers() %>%
  {Idents(.) <- "group_CN_kMeans"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0
  )

all_markers_at1_df %>% head
# p_val avg_log2FC pct.1 pct.2 p_val_adj   cluster   gene
# TG         0  -4.788768 0.081 0.561         0 cluster_1     TG
# PAX8       0  -3.993365 0.079 0.545         0 cluster_1   PAX8
# EPCAM      0  -2.800793 0.187 0.601         0 cluster_1  EPCAM
# FOXE1      0  -4.430288 0.046 0.452         0 cluster_1  FOXE1
# NKX2-1     0  -3.823196 0.036 0.430         0 cluster_1 NKX2-1
# BMP8A      0  -5.928775 0.021 0.389         0 cluster_1  BMP8A
write.csv(all_markers_at1_df, file.path(cn_out_dir, 'DEGs_group_CN_kMeans.csv'), row.names = F)


### consensus clustering ------------
sample_cluster_dist <- results[[3]][["consensusClass"]] %>% as.character() %>% str_c('cluster_', .)
names(sample_cluster_dist) <- names(results[[3]][["consensusClass"]])
names(sample_cluster_dist) <- str_c(names(sample_cluster_dist), '$')
all_data_CN_epis$group_consensus_clusters <- str_replace_all(all_data_CN_epis$sample, c(
  sample_cluster_dist
))
table(all_data_CN_epis$group_consensus_clusters, all_data_CN_epis$sample)

all_markers_at1_df <-
  all_data_CN_epis %>%
  JoinLayers() %>%
  {Idents(.) <- "group_consensus_clusters"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0
  )

all_markers_at1_df %>% head
# p_val avg_log2FC pct.1 pct.2 p_val_adj   cluster   gene
# TG         0  -4.788768 0.081 0.561         0 cluster_1     TG
# PAX8       0  -3.993365 0.079 0.545         0 cluster_1   PAX8
# EPCAM      0  -2.800793 0.187 0.601         0 cluster_1  EPCAM
# FOXE1      0  -4.430288 0.046 0.452         0 cluster_1  FOXE1
# NKX2-1     0  -3.823196 0.036 0.430         0 cluster_1 NKX2-1
# BMP8A      0  -5.928775 0.021 0.389         0 cluster_1  BMP8A
write.csv(all_markers_at1_df, file.path(cn_out_dir, 'DEGs_group_consensus_clusters.csv'), row.names = F)
