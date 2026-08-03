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
adata <- subset(adata,subset=(sample2 %in% c("HLH01-LAHS-B","HLH01-LAHS-B-Post","HLH03-LAHS-T","HLH03-LAHS-T-Post")))
adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.4","RNA_snn_res.0.6","RNA_snn_res.0.8")]
adata$sample<-as.factor(adata$sample)
adata$group<-adata$sample
levels(adata$group) <- c('Before','Post','Before','Post')
levels(adata$sample) <- c('LAHS1','LAHS1-Post','LAHS2','LAHS2-Post')

p1<-DimPlot(adata, reduction = "umap", label = FALSE) + ggsci::scale_color_igv() +  theme(title=element_blank())
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/Dimplot_annotation.pdf", height = 7, width = 10)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/Dimplot_annotation.png", height = 480, width = 700)
print(p1)
dev.off()

p1<-DimPlot(adata, reduction = "umap", label = FALSE,split.by='sample') + ggsci::scale_color_igv() +  theme(title=element_blank())
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/Dimplot_annotation_sample.pdf", height = 7, width = 30)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/Dimplot_annotation_sample.png", height = 480, width = 2400)
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

# adata@graphs$RNA_nn <- NULL
# adata@graphs$RNA_snn <- NULL
# adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.2","RNA_snn_res.0.3")]
adata <- RunUMAP(adata,reduction = "harmony",dims = 1:10)
adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:10)
adata <- FindClusters(adata,resolution = 0.2)
# adata <- RunTSNE(adata,reduction = "harmony",dims = 1:30)
DimPlot(adata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')

DimPlot(adata,reduction='umap',split.by='sample')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_sample.png", width=21, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_sample.pdf", width=21, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_group.png", width=11, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_group.pdf", width=11, height=7, dpi = 300, bg = 'white')

saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_sc_doublet.rds")
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

levels(Idents(adata))[1] <- 'Monocyte'
levels(Idents(adata))[2] <- 'CD8T'
levels(Idents(adata))[3] <- 'CD4T'
levels(Idents(adata))[4] <- 'Monocyte'
levels(Idents(adata))[4] <- 'Neutrophil'
levels(Idents(adata))[5] <- 'T proliferating'
levels(Idents(adata))[6] <- 'B'
levels(Idents(adata))[7] <- 'Other'
levels(Idents(adata))[8] <- 'Platelet'
levels(Idents(adata))[9] <- 'Erythrocyte'
Idents(adata)<-factor(Idents(adata),levels=c("Monocyte","Neutrophil","CD4T","CD8T","T proliferating","B","Other","Erythrocyte","Platelet"))
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_raw.rds")

adata <- subset(adata,idents=c('Other','Platelet','Erythrocyte'), invert = TRUE)
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')

bdata<-subset(adata,idents = 'CD8T')
p1<-DimPlot(bdata, reduction = "umap", label = FALSE) + ggsci::scale_color_igv() +  theme(title=element_blank())
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/Dimplot_annotation.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/Dimplot_annotation.png", height = 480, width = 480)
print(p1)
dev.off()

bdata <- FindClusters(bdata,resolution = 0.4)
# bdata <- RunTSNE(bdata,reduction = "harmony",dims = 1:30)
DimPlot(bdata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')

bdata_diffgene <- FindAllMarkers(object = bdata, only.pos = TRUE)
write.table(bdata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)

FeaturePlot(bdata,features = c('GNLY','NCAM1','NCR1','NKG7','KLRD1','GZMA','FGFBP2','CD8A','CD8B'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/marker_umap.png", bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/marker_umap.pdf", bg = 'white')

levels(Idents(bdata))[1] <- 'NK'
levels(Idents(bdata))[2] <- 'CD8T'
levels(Idents(bdata))[3] <- 'CD8T'
levels(Idents(bdata))[3] <- 'CD8T'
levels(Idents(bdata))[3] <- 'NK'
bdata$cell_type<-Idents(bdata)
DimPlot(bdata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/CD8T-NK/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')

adata_cell_names <- rownames(adata@meta.data)
bdata_cell_names <- rownames(bdata@meta.data)
adata_cell_types <- as.character(adata@meta.data$cell_type)
bdata_cell_types <- as.character(bdata@meta.data$cell_type)
mapped_indices <- match(adata_cell_names, bdata_cell_names)
matched_cell_types <- bdata_cell_types[mapped_indices]
adata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
adata@meta.data$cell_type <- as.factor(adata_cell_types)
adata$cell_type<-factor(adata$cell_type,levels=c("Monocyte","Neutrophil","CD4T","CD8T","T proliferating",'NK',"B"))
Idents(adata)<-adata$cell_type
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

# adata@graphs$RNA_nn <- NULL
# adata@graphs$RNA_snn <- NULL
# adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("RNA_snn_res.0.2","RNA_snn_res.0.3")]
adata <- RunUMAP(adata,reduction = "harmony",dims = 1:10)
adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:10)
adata <- FindClusters(adata,resolution = 0.2)
# adata <- RunTSNE(adata,reduction = "harmony",dims = 1:30)
DimPlot(adata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='sample')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_sample.png", width=21, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_sample.pdf", width=21, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_group.png", width=11, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/cluster_umap_group.pdf", width=11, height=7, dpi = 300, bg = 'white')

saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_sc_doublet.rds")
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

levels(Idents(adata))[1] <- 'Monocyte'
levels(Idents(adata))[2] <- 'CD8T'
levels(Idents(adata))[3] <- 'CD4T'
levels(Idents(adata))[4] <- 'Monocyte'
levels(Idents(adata))[4] <- 'Neutrophil'
levels(Idents(adata))[5] <- 'T proliferating'
levels(Idents(adata))[6] <- 'B'
levels(Idents(adata))[7] <- 'Other'
levels(Idents(adata))[8] <- 'Platelet'
levels(Idents(adata))[9] <- 'Erythrocyte'
Idents(adata)<-factor(Idents(adata),levels=c("Monocyte","Neutrophil","CD4T","CD8T","T proliferating","B","Other","Erythrocyte","Platelet"))
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_raw.rds")

adata <- subset(adata,idents=c('Other','Platelet','Erythrocyte'), invert = TRUE)
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')



adata<-subset(adata,idents = 'Monocyte')
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
set.seed(666)
adata <- RunHarmony(adata, group.by.vars = "sample")

p1<-ElbowPlot(adata,ndims = 50)
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/ElbowPlot.pdf", height = 7, width = 7)
print(p1)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/ElbowPlot.png", height = 480, width = 480)
print(p1)
dev.off()

adata <- RunUMAP(adata,reduction = "harmony",dims = 1:12)
adata <- FindNeighbors(adata,reduction = "harmony",dims = 1:12)
adata <- FindClusters(adata,resolution = 0.2)
# adata <- RunTSNE(adata,reduction = "harmony",dims = 1:30)
DimPlot(adata,reduction='umap',label = TRUE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',label = TRUE,group.by='annotation2')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap.png", width=7, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap.pdf", width=7, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='sample')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap_sample.png", width=21, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap_sample.pdf", width=21, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap_group.png", width=11, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/cluster_umap_group.pdf", width=11, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/Monocyte_sc_doublet.rds")
adata_diffgene <- FindAllMarkers(object = adata, only.pos = TRUE)
write.table(adata_diffgene,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/Monocyte_diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
# adata_diffgene<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/Monocyte_diffgene.xls",header=TRUE)
top10 <- adata_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# all.genes <- rownames(adata_diffgene)
adata <- ScaleData(adata, features = top10$gene)
DoHeatmap(adata, features = top10$gene)+NoLegend()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/top10_diffgene_heatmap.png", width=7, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/top10_diffgene_heatmap.pdf", width=7, height=10, dpi = 300, bg = 'white')

FeaturePlot(adata,features = c('CD14', 'FCGR3A','MKI67','TOP2A','STMN1','TPX2','PTTG1'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')
FeaturePlot(adata,features = c('CD14', 'FCGR3A','PTTG1','ZFP36L1','S100A8'))
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/marker_umap.png", width=11, height=10, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/marker_umap.pdf", width=11, height=10, dpi = 300, bg = 'white')

levels(Idents(adata))[1] <- 'Classical monocyte'
levels(Idents(adata))[2] <- 'Classical monocyte'
levels(Idents(adata))[2] <- 'Non-classical monocyte'
levels(Idents(adata))[3] <- 'Classical monocyte'
levels(Idents(adata))[3] <- 'Other'
levels(Idents(adata))[4] <- 'Classical monocyte'
adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap_raw.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap_raw.pdf", width=8, height=7, dpi = 300, bg = 'white')
adata <- subset(adata,idents=c('Other'), invert = TRUE)
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/Monocyte_raw.rds")

adata$cell_type <- Idents(adata)
DimPlot(adata,reduction='umap')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap.png", width=8, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap.pdf", width=8, height=7, dpi = 300, bg = 'white')
DimPlot(adata,reduction='umap',split.by='group')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap_group.png", width=15, height=7, dpi = 300, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_umap_group.pdf", width=15, height=7, dpi = 300, bg = 'white')
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/Monocyte_anno.rds")

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_proportion_group.pdf", width=7, height=7,bg = 'white')

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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_proportion_sample.png", width=12, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte/anno_proportion_sample.pdf", width=12, height=7,bg = 'white')