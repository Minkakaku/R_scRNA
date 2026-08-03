library(Seurat)
library(tidyverse)
library(ConsensusClusterPlus)

out_dir <- 'Epi_consensus_clustering'
dir.create(out_dir)

all_data <- readRDS("~/R_OUTPUT/data_for_analysis/all_data_filtered2.rds")
all_data$annotation %>% table()
all_data$group %>% table()

Epithelial <- subset(all_data, annotation == 'Epithelial')

Epithelial_ls <- SplitObject(Epithelial, split.by = 'sample')

expr_df <- 
  Epithelial_ls %>% 
  imap(function(Epithelial, sample){ # browser()
    expr_matrix <- LayerData(JoinLayers(Epithelial), 'data')
    expr <- rowMeans(expr_matrix) %>% as.data.frame()
    colnames(expr) <- sample
    expr
  }) %>% 
  bind_cols()


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

setwd(out_dir)
set.seed(12345)
results <- ConsensusClusterPlus(
  as.matrix(expr_df), 
  maxK = 20, 
  reps = 1000, 
  plot = "png", 
  seed = 12345
)
setwd('~/')




# 中位数中心化 -----------
#sweep函数减去中位数进行标准化
exprSet = sweep(expr_df,1, apply(expr_df,1,median,na.rm=T))
par(mfrow = c(1,2))
boxplot(expr_df[,1:20],main = "before")
boxplot(exprSet[,1:20],main = "after")

setwd(out_dir)
set.seed(12345)
results <- ConsensusClusterPlus(
  as.matrix(exprSet), 
  maxK = 20, 
  reps = 1000, 
  plot = "png", 
  seed = 12345
)
setwd('~/')




# 使用PAC标准去筛选K值 ------------
maxK = 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK## [1] 5




#输出K=5时的一致性矩阵
results[[5]][["consensusMatrix"]][1:5,1:5]

# colnames(consensusMatrix_select) <- colnames(exprSet)
# rownames(consensusMatrix_select) <- colnames(exprSet)

pheatmap(results[[5]][["consensusMatrix"]],
         cluster_rows = F,
         cluster_cols = TRUE,
         # annotation_row = annotation,
         # scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)

#hclust聚类树
results[[5]][["consensusTree"]]
#Call:
#hclust(d = as.dist(1 - fm), method = finalLinkage)
#Cluster method   : average 
#Number of objects: 128 

#查看样本所属的聚类群
table(results[[5]]$consensusClass)
results[[5]][["consensusClass"]]




# 计算聚类一致性 (cluster-consensus，CLC) 和样品一致性 (item-consensus，IC) -----------
icl <- calcICL(results, title = "./ConsensusCluster",plot = "png", writeTable = TRUE)
## 返回了具有两个元素的list
dim(icl[["clusterConsensus"]])#[1] 20  3
icl[["clusterConsensus"]] 




# plot ------------
## 热图 ------------
message('heatmap')
library(pheatmap)
library(RColorBrewer)

# 准备热图数据
heatmap_data <- expr_df[order(results[[2]][["consensusClass"]])] %>% t()
cluster_order <- results[[2]][["consensusClass"]][order(results[[2]][["consensusClass"]])]

# 创建注释
annotation <- data.frame(Cluster = factor(cluster_order))
rownames(annotation) <- rownames(heatmap_data)

# 绘制热图
pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k2_2.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = TRUE,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()


# 绘制热图
pdf(file.path(out_dir, 'Clustering Heatmap_k2_2.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Clustering Heatmap",
         fontsize_row = 8,
         treeheight_row = 0,
         show_colnames = F,
         show_rownames = T)
dev.off()




## 簇间距离矩阵可视化 --------------
# 计算簇中心距离
center_dist <- dist(results[[2]][["consensusClass"]])
center_dist_matrix <- as.matrix(center_dist)

# 绘制簇间距离热图
pdf(file.path(out_dir, 'Distance Between Cluster Centers_k2_2.pdf'), height = 6, width  = 11)
pheatmap(center_dist_matrix,
         main = "Distance Between Cluster Centers",
         color = colorRampPalette(c("white", "red"))(50),
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()




## PCA ------------
# 如果特征维度较高，先用PCA降维
pca_result <- prcomp(t(expr_df), scale. = TRUE)

# PCA结果可视化
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$cluster <- factor(results[[2]][["consensusClass"]])

ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(title = "Consensus Clustering on PCA-reduced Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, 'Consensus Clustering on PCA-reduced Data_k2_2.pdf'), height = 6, width  = 11)


# 标记离群点
pct = 0.8
pca_data2 <- pca_data %>%
  group_by(cluster) %>%
  mutate(
    # 计算聚类中心
    center_x = mean(PC1),
    center_y = mean(PC2),
    # 计算到聚类中心的欧氏距离
    dist_to_center = sqrt((PC1 - center_x)^2 + (PC2 - center_y)^2)
  ) %>%
  ungroup() %>% 
  as.data.frame()
rownames(pca_data2) <- rownames(pca_data)

# 定义离群点阈值（例如：距离大于第95百分位数）
pca_data3 <- pca_data2 %>%
  group_by(cluster) %>%
  mutate(
    is_outlier = dist_to_center > quantile(dist_to_center, pct)
  ) %>%
  ungroup() %>% 
  as.data.frame()
rownames(pca_data3) <- rownames(pca_data2)


# 绘制带离群点标记的图
ggplot(pca_data3, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = pct) +
  # # 突出显示离群点
  # geom_point(data = subset(pca_data3, is_outlier), 
  #            aes(x = PC1, y = PC2, color = cluster),
  #            size = 5, shape = 1, stroke = 1.5) +
  # 添加文本标签
  geom_text(data = subset(pca_data3, is_outlier),
            aes(x = PC1, y = PC2, label = rownames(pca_data3[pca_data3$is_outlier, ])),
            vjust = -0.5, size = 3) +
  labs(title = "Consensus Clustering on PCA-reduced Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, 'Consensus Clustering on PCA-reduced Data_k2_outlier_2.pdf'), height = 6, width  = 11)




# DEGs ----------- 
consensus_cluster <- str_c('cluster_', results[[2]][["consensusClass"]])
names(consensus_cluster) <- names(results[[2]][["consensusClass"]])
Epithelial$consensus_cluster <- str_replace_all(Epithelial$sample, consensus_cluster)
Idents(Epithelial)<-'consensus_cluster'
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Epithelial_markers, file.path(out_dir, 'Epithelial_markers_consensus_clusters_k2.csv'), row.names = F)

Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

## heatmap ---------
# 准备热图数据
heatmap_data <- expr_df[unlist(Epithelial_de_genes_for_enrichmenet), ] %>% t()
cluster_order <- results[[2]][["consensusClass"]][order(results[[2]][["consensusClass"]])]

# 创建注释
annotation <- data.frame(Cluster = factor(cluster_order))
rownames(annotation) <- rownames(heatmap_data)

# 绘制热图
pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k2_DEGs.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = TRUE,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()

pdf(file.path(out_dir, 'Consensus Clustering Heatmap_k2_DEGs_rowCluster.pdf'), height = 6, width  = 11)
pheatmap(heatmap_data,
         cluster_rows = T,
         cluster_cols = TRUE,
         annotation_row = annotation,
         scale = 'column', 
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Consensus Clustering Heatmap",
         fontsize_row = 8,
         show_colnames = F,
         show_rownames = T)
dev.off()


# enrichment ---------------------------------------
Epithelial_er_output <- file.path(out_dir, 'enrichment') 
dir.create(Epithelial_er_output)

enrich_ls <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Epithelial_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(Epithelial_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Epithelial_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Epithelial_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Epithelial_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

