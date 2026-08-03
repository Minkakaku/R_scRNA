# conda activate scMetabolism
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/08.Mission/annotation2/all_data_sc_doublet.rds")
adata
levels(adata$annotation)
head(adata@meta.data)
adata[["RNA"]] <- adata[["RNA3"]]
DefaultAssay(adata) <- 'RNA'
adata[["RNA3"]] <- NULL
adata <- subset(adata,subset=(group == 'LAHS'))
adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.4","RNA_snn_res.0.6","RNA_snn_res.0.8")]
adata$sample<-as.factor(adata$sample)
adata$group<-adata$sample
levels(adata$group) <- c('Before','Post','Before','Before','Post','Before')
levels(adata$sample) <- c('LAHS1','LAHS1-Post','LAHS4','LAHS2','LAHS2-Post','LAHS3')

p1<-DimPlot(adata, reduction = "umap", label = FALSE) + ggsci::scale_color_igv() +  theme(title=element_blank())
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/Dimplot_annotation.pdf", height = 7, width = 8)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/Dimplot_annotation.png", height = 480, width = 560)
print(p1)
dev.off()

p1<-DimPlot(adata, reduction = "umap", label = FALSE,split.by='sample') + ggsci::scale_color_igv() +  theme(title=element_blank())
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/Dimplot_annotation_sample.pdf", height = 7, width = 30)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/Dimplot_annotation_sample.png", height = 480, width = 2400)
print(p1)
dev.off()

adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata)
adata@reductions$pca <- NULL
adata@reductions$umap <- NULL
adata@reductions$harmony <- NULL
adata@reductions$harmonyumap <- NULL
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
set.seed(666)
adata <- RunHarmony(adata, group.by.vars = "sample")

p1<-ElbowPlot(adata,ndims = 50)
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

# adata@graphs$RNA_nn <- NULL
# adata@graphs$RNA_snn <- NULL
# adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.2","RNA_snn_res.0.3")]
adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:20)
adata <- FindClusters(adata,resolution = 0.3)
adata <- RunUMAP(adata,reduction = "harmony",dims = 1:20)
# adata <- RunTSNE(adata,reduction = "harmony",dims = 1:30)
DimPlot(adata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
FeaturePlot(adata,features = c('CD14', 'FCGR3A','PTTG1','ZFP36L1','S100A8'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')

DimPlot(adata,reduction='umap',split.by='sample')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap_sample.png", width=21, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap_sample.pdf", width=21, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap_group.png", width=11, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/cluster_umap_group.pdf", width=11, height=7, dpi = 300, bg = 'white')

saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_sc_doublet.rds")
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

levels(Idents(adata))[1] <- 'Monocyte'
levels(Idents(adata))[2] <- 'CD4T'
levels(Idents(adata))[3] <- 'Neutrophil'
levels(Idents(adata))[4] <- 'CD8T'
levels(Idents(adata))[5] <- 'T proliferating'
levels(Idents(adata))[6] <- 'Monocyte'
levels(Idents(adata))[7] <- 'NK'
levels(Idents(adata))[8] <- 'B'
levels(Idents(adata))[9] <- 'NK'
levels(Idents(adata))[10] <- 'Monocyte'
levels(Idents(adata))[11] <- 'Platelet'
levels(Idents(adata))[12] <- 'Erythrocyte'
levels(Idents(adata))[13] <- 'Monocyte'
levels(Idents(adata))[14] <- 'B'
levels(Idents(adata))[15] <- 'Neutrophil'
levels(Idents(adata))[16] <- 'Monocyte'
levels(Idents(adata))[17] <- 'Other'
levels(Idents(adata)) <- c('Monocyte','CD4T','Neutrophil','CD8T','T proliferating','Monocyte','NK','B','NK','Monocyte','Platelet','Erythrocyte','Monocyte','B','Neutrophil','Monocyte','Other')
Idents(adata)<-factor(Idents(adata),levels=c("Monocyte","Neutrophil","CD4T","CD8T","T proliferating","NK","B","Other","Erythrocyte","Platelet"))
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_raw.rds")

adata <- subset(adata,idents=c('Other','Platelet','Erythrocyte'), invert = TRUE)
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_anno.rds")

my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
count_df <- table(adata$group, adata$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

count_df <- table(adata$sample, adata$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("sample", "cell_type", "proportion")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')



adata<-subset(adata,idents = 'Monocyte')
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
set.seed(666)
adata <- RunHarmony(adata, group.by.vars = "sample")

p1<-ElbowPlot(adata,ndims = 50)
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:15)
adata <- FindClusters(adata,resolution = 0.5)
adata <- RunUMAP(adata,reduction = "harmony",dims = 1:15)
# adata <- RunTSNE(adata,reduction = "harmony",dims = 1:30)
DimPlot(adata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',label = TRUE,group.by='annotation2')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='sample')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap_sample.png", width=21, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap_sample.pdf", width=21, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap_group.png", width=11, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/cluster_umap_group.pdf", width=11, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_sc_doublet.rds")
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

FeaturePlot(adata,features = c('CD14', 'FCGR3A','MKI67','TOP2A','STMN1','TPX2','PTTG1'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')
FeaturePlot(adata,features = c('CD14', 'FCGR3A','PTTG1','ZFP36L1','S100A8'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')

bdata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_anno.rds")
adata_cell_names <- rownames(adata@meta.data)
bdata_cell_names <- rownames(bdata@meta.data)
adata_cell_types <- as.character(adata@meta.data$RNA_snn_res.0.5)
bdata_cell_types <- as.character(bdata@meta.data$cell_type)
# 匹配行名并替换对应列的值，保留未匹配行的原始值
mapped_indices <- match(bdata_cell_names, adata_cell_names)
matched_cell_types <- adata_cell_types[mapped_indices]
bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
# 将新的细胞类型注释赋值给bdata，并转换回factor类型
bdata@meta.data$cell_type <- as.factor(bdata_cell_types)
Idents(bdata) <- bdata@meta.data$cell_type
DimPlot(bdata,reduction='umap',label=TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')

levels(Idents(adata))[1] <- 'S100A8 classical monocyte'
levels(Idents(adata))[2] <- 'ZFP36L1 classical monocyte'
levels(Idents(adata))[3] <- 'Non-classical monocyte'
levels(Idents(adata))[4] <- 'Classical monocyte'
levels(Idents(adata))[5] <- 'ZFP36L1 classical monocyte'
levels(Idents(adata))[6] <- 'Dendritic cell'
levels(Idents(adata))[7] <- 'Other'
levels(Idents(adata))[8] <- 'PTTG1 classical monocyte'
levels(Idents(adata))[9] <- 'Other'
levels(Idents(adata))[10] <- 'Dendritic cell'
levels(Idents(adata)) <- c('S100A8 classical monocyte','ZFP36L1 classical monocyte','Non-classical monocyte','Classical monocyte','ZFP36L1 classical monocyte','Dendritic cell','Other','PTTG1 classical monocyte','Other','Dendritic cell')
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')
adata <- subset(adata,idents=c('Other'), invert = TRUE)
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_raw.rds")

adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata)
adata@reductions$pca <- NULL
adata@reductions$umap <- NULL
adata@reductions$harmony <- NULL
adata@reductions$harmonyumap <- NULL
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
set.seed(666)
adata <- RunHarmony(adata, group.by.vars = "sample")

p1<-ElbowPlot(adata,ndims = 50)
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

# adata@graphs$RNA_nn <- NULL
# adata@graphs$RNA_snn <- NULL
adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.2","RNA_snn_res.0.3")]
adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:25)
adata <- FindClusters(adata,resolution = 0.5)
adata <- RunUMAP(adata,reduction = "harmony",dims = 1:25)
FeaturePlot(adata,features = c('CD14', 'FCGR3A','PTTG1','ZFP36L1','S100A8'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

levels(Idents(adata)) <- c('Classical monocyte','ZFP36L1 classical monocyte','Non-classical monocyte','S100A8 classical monocyte','Dendritic cell','S100A8 classical monocyte','Classical monocyte','PTTG1 classical monocyte','Dendritic cell','Classical monocyte')
Idents(adata)<-factor(Idents(adata),levels=c('Non-classical monocyte','PTTG1 classical monocyte','ZFP36L1 classical monocyte','S100A8 classical monocyte','Classical monocyte','Dendritic cell'))
adata$cell_type <- Idents(adata)
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_anno.rds")

DimPlot(adata,reduction='umap',group.by='cell_type')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap.png", width=9, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap.pdf", width=9, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')

my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
count_df <- table(adata$group, adata$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

count_df <- table(adata$sample, adata$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("sample", "cell_type", "proportion")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')

metadata <- adata@meta.data
write.table(metadata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/all_cells_features.xls",sep="\t")