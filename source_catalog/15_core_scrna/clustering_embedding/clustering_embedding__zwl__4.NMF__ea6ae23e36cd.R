setwd('/home/zwl/Project/FW2024-428-02/6.NMF/')
rm(list = ls())
gc()
seurat = readRDS("./12.2_seurat_raw_count_nonormal.rds")
seurat <- NormalizeData(seurat) %>% FindVariableFeatures() %>% ScaleData(do.center = F)#promise all value >0

# NMF clustering
library(NMF)
library(dplyr)
library(tidyverse)

#data.matrix = as.matrix(seurat@assays$RNA$counts)
data.matrix = readRDS('./feature2000_martix.rds') %>% as.matrix()
# Check for and remove rows with null or NA values
data.matrix <- data.matrix[rowSums(is.na(data.matrix)) == 0, ]
# Also remove rows with zero variance (all zeros)
data.matrix <- data.matrix[rowSums(data.matrix != 0) > 0, ]

nmf_result <- nmf(data.matrix, rank = 20, method = "brunet", seed = 123)
## 储存NMF的结果
embeddings_matrix <- t(coef(nmf_result))
colnames(embeddings_matrix) <- paste0("nmf_", 1:ncol(embeddings_matrix))
seurat@reductions$nmf <- CreateDimReducObject(
  embeddings = embeddings_matrix,
  loadings = basis(nmf_result),
  assay = DefaultAssay(seurat),
  key = "nmf_")

library(Seurat)
set.seed(123)
data.nmf <- RunUMAP(seurat, 
                    reduction='nmf', 
                    dims = 1:20) %>% FindNeighbors(reduction = "nmf", dims = 1:20)%>%
  FindClusters(resolution = 0.5)
data.nmf = FindClusters(data.nmf,resolution = 0.3)

## 聚类结果可视化  
DimPlot(data.nmf, label = T) + ggsci::scale_color_igv()
p0 = DimPlot(data.nmf,group.by = c('RNA_snn_res.0.3','cell_type','celltype2'), label = F) + ggsci::scale_color_igv()
ggsave(filename = '0.NMF_umap.png',w=12,h=4,p0)

## 查看细胞在NMF因子上的权重
coef_df <- data.frame(t(coef(nmf_result)), check.names = F)
colnames(coef_df) <- paste0("factor", 1:dim(coef_df)[2])
data.nmf <- AddMetaData(data.nmf, metadata = coef_df)

## 细胞在NMF因子上的权重可视化
p <- FeaturePlot(data.nmf, features = paste0("factor", 1:dim(coef_df)[2]), ncol = 4)
p
ggsave(filename = '1.NMF_featureplot.png',w=12,h=12,p)

## 提取所有基因模块中特征基因，每个模块30个特征基因
library(reshape2)
signatureFeature <- extractFeatures(nmf_result, 30L)
signatureFeature <- lapply(signatureFeature, function(x) rownames(nmf_result)[x])
signatureFeature <- do.call("rbind", signatureFeature)
rownames(signatureFeature) <- paste0("factor", 1:(dim(coef_df)[2]))
colnames(signatureFeature) <- paste('gene', 1:30, sep = '_')
signatureFeature.out <- as.data.frame(t(signatureFeature)) %>% rownames_to_column(var = "genes") %>% melt(id.vars ="genes")
write.csv(file = './signatureFeature_out.csv',signatureFeature.out)

library(pheatmap)
library(data.table)
load("FW2024-428-01/NMF_res.RData")
meta = data.nmf@meta.data %>% dplyr::select(seurat_clusters) %>% merge(as.data.frame(embeddings_matrix), by="row.names") %>% melt(id.vars = c("Row.names","seurat_clusters")) %>% 
  group_by(seurat_clusters, variable) %>% summarise(value = mean(value)) %>% dcast(seurat_clusters ~ variable)
meta$seurat_clusters = paste0("C", meta$seurat_clusters)
colnames(meta)[-1] <- gsub("nmf_","Module", colnames(meta)[-1])
meta = meta %>% column_to_rownames(var ="seurat_clusters")
p = pheatmap(t(scale(meta)), treeheight_row = 0,color=colorRampPalette(c("darkblue","white","#EE0000FF"))(100))
ggsave("./NMF_heatmap.pdf",plot = p, width = 6, height = 8)

nmf_matrix = as.data.frame(t(scale(meta))) %>% rownames_to_column(var ="Module")
fwrite(nmf_matrix,"./NMF_heatmap data.xls",row.names = F, sep = "\t")

# save.image("FW2024-428-01/NMF_res.RData")
# 寻找每列最大值对应的行名和列名
max_rows <- apply(meta, 2, function(x) rownames(meta)[which.max(x)])
max_cols <- colnames(meta)
c2m <- data.frame(variable = gsub('Module','factor',max_cols), cluster = max_rows,row.names = NULL)

signatureFeature.out = merge(signatureFeature.out, c2m, by="variable")

library(clusterProfiler)
library(org.Hs.eg.db)
entrezid = bitr(signatureFeature.out$value, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
signatureFeature.out = merge(signatureFeature.out, entrezid, by.x = "value", by.y = "SYMBOL", all.x =TRUE)
go_res_nmf <- lapply(unique(signatureFeature.out$cluster), function(x){
  go.res = enrichGO(gene = unique(signatureFeature.out$ENTREZID[signatureFeature.out$cluster == x]),
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, 
                    readable = TRUE)
  return(go.res)
})
names(go_res_nmf) = unique(signatureFeature.out$cluster)
temp = lapply(names(go_res_nmf), function(x) {
  t = go_res_nmf[[x]]@result %>% filter(nchar(Description) < 60) %>%  top_n(Count, n=6) %>% arrange(desc(Count)) 
  t = t[1:6,]
})
names(temp) = unique(signatureFeature.out$cluster)
temp = rbindlist(temp, idcol = "NMF_cluster" )
fwrite(temp,file = "./GO enrichment results by NMF_cluster.xls",sep="\t",row.names = F)


library(Hmisc)
go.plot = lapply(names(go_res_nmf), function(x){
  temp = go_res_nmf[[x]]@result %>% filter(nchar(Description) < 60) %>%  top_n(Count, n=6) %>% arrange(desc(Count)) 
  p = temp[1:6,] %>% ggplot(aes(x =Count, y = capitalize(Description))) + geom_bar(stat= "identity", color = "#F2F2F2") + theme_classic(base_size = 14) + ylab("")+ xlab("") + ggtitle(x) +theme(title = element_text(size = 16, hjust = 0.5))
  return(p)
})

library(egg)
library(ggpubr)
p = ggarrange(plotlist= go.plot, ncol= 4, nrow=4)
ggsave("./GO_enrichment_by_NMF_defined_clusters.pdf", plot = p,width = 24, height = 18)
