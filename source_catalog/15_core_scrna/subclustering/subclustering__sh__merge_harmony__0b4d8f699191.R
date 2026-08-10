.libPaths()
library(Seurat)

a<-readRDS('Main_mapping.rds')
b<-readRDS("D:/项目/FWAE20252866/FWAE20252866/FWAE20252866.rds")

a<-subset(a,group=='AS')
b
head(b@meta.data)

library(Seurat)
library(Matrix)
# 假设 adata 是 Seurat v5 对象
# 1. 提取 counts 和 data 层
counts_matrix <- as(adata[["RNA"]]$counts, "dgCMatrix")
adata<-NormalizeData(adata)
data_matrix   <- as(adata[["RNA"]]$data, "dgCMatrix")
# 2. 提取 feature 名称（行名）
features <- rownames(adata)
# 3. 提取 cell 名称（列名）
cells <- colnames(adata)
# 4. 提取 meta.data
meta <- adata@meta.data
# 5. 构建 Seurat v4 对象
seu_v4 <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = meta
)
# 6. 添加 normalized data（可选）
seu_v4[["RNA"]]@data <- data_matrix
# 7. 添加降维结果（如果有的话）
# 例如 PCA
if ("pca" %in% Reductions(adata)) {
  seu_v4[["pca"]] <- adata[["pca"]]
}
# 例如 UMAP
if ("umap" %in% Reductions(adata)) {
  seu_v4[["umap"]] <- adata[["umap"]]
}
# 8. 保存为 v4 兼容的 rds
saveRDS(seu_v4, file = "seurat_v4_compat.rds")
# 检查 assay 类型
seu_v4 <- JoinLayers(seu_v4)
seu_v4[["RNA"]] <- as(object = seu_v4[["RNA"]], Class = "Assay")

b<-readRDS("D:/项目/FWAE20252866/02.Mission/seurat_v4_compat.rds")
b$sample<-b$orig.ident
b$group<-'control'
combined <- merge(a, y = b, add.cell.ids = c("AS", "control"))
saveRDS(combined,'D:/项目/FWAE20252866/02.Mission/combined.rds')

S$sample<-as.factor(S$sample)
levels(S$sample)[5]<-'AS1'
levels(S$sample)[6]<-'AS1'
levels(S$sample)[6]<-'AS2'
levels(S$sample)[7]<-'AS2'
levels(S$sample)[7]<-'AS3'
levels(S$sample)[8]<-'AS3'
S<-combined
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
library(harmony)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50)
S <- FindNeighbors(S,reduction = "harmony",dims = 1:15)
S <- FindClusters(S,resolution = 0.5)
S <- RunUMAP(S,reduction = "harmony",dims = 1:15)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap')
DimPlot(S,reduction='umap',split.by='sample')

#levels(S$cell_type)[2]<-'Mast cells'
S$cell_type<-as.character(S$cell_type)
S$sample<-as.character(S$sample)
S$group<-as.character(S$group)
saveRDS(S,'D:/项目/FWAE20252866/02.Mission/result_harmony.rds')

library(ggplot2)
setwd('D:/项目/FWAE20252866/02.Mission/')
getwd()
S_sub<-readRDS('result_harmony.rds')

DimPlot(S_sub,reduction='umap',label = TRUE)
ggsave('D:/项目/FWAE20252866/02.Mission/Main/sub_cluster_umap.png',width =7,height =6)
ggsave('D:/项目/FWAE20252866/02.Mission/Main/sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S_sub,reduction='umap',label = TRUE,split.by = 'group')
ggsave('D:/项目/FWAE20252866/02.Mission/Main/umap_split_groups.png',width =14,height =6)
ggsave('D:/项目/FWAE20252866/02.Mission/Main/umap_split_groups.pdf',width =14,height =6)

S_diffgene <- FindAllMarkers(object = S_sub, only.pos = TRUE)
write.table(S_diffgene,"./Main/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
S_diffgene <- read.table("./Main/diffgene.xls",header = TRUE)
library(magrittr)
library(dplyr)
top10 <- S_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
all.genes <- rownames(top10$gene)
S_sub <- ScaleData(S_sub, features = all.genes)
DoHeatmap(S_sub, features = top10$gene)+NoLegend()
ggsave(file="./Main/top10_diffgene_heatmap.png", width=16, height=21, bg = 'white')
ggsave(file="./Main/top10_diffgene_heatmap.pdf", width=16, height=21, bg = 'white')

FeaturePlot(S_sub,features = c('BPIFA2','TFF1'))
FeaturePlot(S_sub,features = c('CLDN4', 'TM4SF1','EPCAM','KRT14'))
DimPlot(S_sub,reduction='umap')
levels(Idents(S_sub)) <- c('Fibroblast','Fibroblast','Epithelial','Epithelial','Myofibroblast','Epithelial','Epithelial','Endothelial cell','Ciliated cell','NK cell','DC','Monocyte','T cell','Fibroblast','Proliferating cell','Epithelial','Epithelial')
DimPlot(S_sub,reduction='umap')
names(S_sub@meta.data)
S_sub@meta.data <- S_sub@meta.data[, !names(S_sub@meta.data) %in% c("RNA_snn_res.0.5","orig.ident","cluster","sub_cell_type")]
S_sub$cell_type <- Idents(S_sub)
metadata <- S_sub@meta.data
write.table(metadata,"./Main/all_cells_features.xls",sep="\t")
DimPlot(S_sub,reduction='umap',group.by = 'cell_type')
ggsave(file="./Main/dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="./Main/dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')
DimPlot(S_sub,reduction='umap',split.by = 'group')
ggsave(file="./Main/dimplot_umap_group_celltype.png", width=14, height=6, bg = 'white')
ggsave(file="./Main/dimplot_umap_group_celltype.pdf", width=14, height=6, bg = 'white')



saveRDS(S_sub,'./Main/Main.rds')
