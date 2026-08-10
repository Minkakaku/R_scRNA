.libPaths()
library(Seurat)
library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)

S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
# "Fibroblast" ,"Epithelial","Myofibroblast","Endothelial cell","Ciliated cell","NK cell","DC","Monocyte","T cell","Proliferating cell"

# 巨噬细胞细分，基本全是M1
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Monocyte/")
getwd()
S<-subset(S,cell_type=='Monocyte')
DimPlot(S,reduction='umap')
ggsave("umap_mapping.pdf",width = 8,height = 7)
ggsave("umap_mapping.png",width = 8,height = 7)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(999)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50)
ggsave("ElbowPlot.pdf",width = 7,height = 7,bg='white')
ggsave("ElbowPlot.png",width = 7,height = 7,bg='white')
S <- FindNeighbors(S,reduction = "harmony",dims = 1:10)
S <- FindClusters(S,resolution = 0.2)
S <- RunUMAP(S,reduction = "harmony",dims = 1:10)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap')
ggsave("sub_cluster_umap.pdf",width = 7,height = 7,bg='white')
ggsave("sub_cluster_umap.png",width = 7,height = 7,bg='white')
DimPlot(S,reduction='umap',split.by='group')
ggsave("umap_split_groups.pdf",width = 7,height = 7,bg='white')
ggsave("umap_split_groups.png",width = 7,height = 7,bg='white')

FeaturePlot(S,features = c('CD14','LYZ','IL1B','TNFA','MRC1','CD163','CLEC10A','CD1C'))
ggsave("umap_marker.pdf",width = 21,height = 21,bg='white')
ggsave("umap_marker.png",width = 21,height = 21,bg='white')
saveRDS(S,'Macrophage.rds')

S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
# 成纤维细分
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Fibroblast/")
getwd()
S<-subset(S,cell_type=='Fibroblast')
DimPlot(S,reduction='umap',label = TRUE)
ggsave('sub_cluster_umap.png',width =7,height =6)
ggsave('sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S,reduction='umap',label = TRUE,split.by = 'group')
ggsave('umap_split_groups.png',width =14,height =6)
ggsave('umap_split_groups.pdf',width =14,height =6)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(999)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50)
ggsave("ElbowPlot.pdf",width = 7,height = 7,bg='white')
ggsave("ElbowPlot.png",width = 7,height = 7,bg='white')
S <- FindNeighbors(S,reduction = "harmony",dims = 1:20)
S <- FindClusters(S,resolution = 0.3)
S <- RunUMAP(S,reduction = "harmony",dims = 1:20)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap',label='TRUE')
ggsave("sub_cluster_umap.pdf",width = 7,height = 7,bg='white')
ggsave("sub_cluster_umap.png",width = 7,height = 7,bg='white')
DimPlot(S,reduction='umap',split.by='group')
ggsave("umap_split_groups.pdf",width = 14,height = 7,bg='white')
ggsave("umap_split_groups.png",width = 14,height = 7,bg='white')
saveRDS(S,"Fibroblast.rds")

FeaturePlot(S,features = c('COL6A1', 'COL6A3','FAP','THY1'))
ggsave("CAF_marker.pdf",width = 7,height = 7,bg='white')
ggsave("CAF_marker.png",width = 7,height = 7,bg='white')
S_diffgene <- FindAllMarkers(object = S, only.pos = TRUE)
write.table(S_diffgene,"diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
S_diffgene <- read.table("diffgene.xls",header = TRUE)
S_diffgene <- S_diffgene[which(S_diffgene$cluster!='5'),]
top10 <- S_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
all.genes <- rownames(top10$gene)
S <- ScaleData(S, features = all.genes)
DoHeatmap(S, features = top10$gene)+NoLegend()
ggsave(file="top10_diffgene_heatmap.png", width=12, height=16, bg = 'white')
ggsave(file="top10_diffgene_heatmap.pdf", width=12, height=16, bg = 'white')

S<-subset(S,idents=5,invert=TRUE)
levels(Idents(S)) <- c('Fibroblast C1','Fibroblast C2','Fibroblast C2','Fibroblast C3','Fibroblast C4')
DimPlot(S,reduction='umap')
ggsave(file="dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')
DimPlot(S,reduction='umap',split.by='group')
ggsave(file="dimplot_umap_celltype_group.png", width=14, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype_group.pdf", width=14, height=6, bg = 'white')
S$cell_type<-Idents(S)

count_df <- table(S$sample, S$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("sample", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_sample.png", width=14, height=7,bg = 'white')
ggsave(file="celltype_proportion_sample.pdf", width=14, height=7,bg = 'white')

count_df <- table(S$group, S$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

metadata <- S@meta.data
write.table(metadata,"all_cells_features.xls",sep="\t")
saveRDS(S,"Fibroblast_anno.rds")


S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
# "NK cell","T cell"细分
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/TandNK/")
getwd()
S<-subset(S,cell_type %in% c("NK cell","T cell"))
DimPlot(S,reduction='umap',label = TRUE)
ggsave('sub_cluster_umap.png',width =7,height =6)
ggsave('sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S,reduction='umap',label = TRUE,split.by = 'group')
ggsave('umap_split_groups.png',width =14,height =6)
ggsave('umap_split_groups.pdf',width =14,height =6)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(999)
S <- RunHarmony(S, group.by.vars = "group")

ElbowPlot(S,ndims = 50)
ggsave("ElbowPlot.pdf",width = 7,height = 7,bg='white')
ggsave("ElbowPlot.png",width = 7,height = 7,bg='white')
S <- FindNeighbors(S,reduction = "harmony",dims = 1:15)
S <- FindClusters(S,resolution = 0.3)
S <- RunUMAP(S,reduction = "harmony",dims = 1:15)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap',label='TRUE')
ggsave("sub_cluster_umap.pdf",width = 7,height = 7,bg='white')
ggsave("sub_cluster_umap.png",width = 7,height = 7,bg='white')
DimPlot(S,reduction='umap',split.by='group')
ggsave("umap_split_groups.pdf",width = 14,height = 7,bg='white')
ggsave("umap_split_groups.png",width = 14,height = 7,bg='white')
saveRDS(S,"TandNK.rds")

FeaturePlot(S,features = c('CD4', 'CD8A','CD8B','TRAC','TRDC','GNLY','NKG7','GZMA'))
ggsave("TandNK_marker.pdf",width = 21,height = 21,bg='white')
ggsave("TandNK_marker.png",width = 21,height = 21,bg='white')


S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
# 上皮细分
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Epithelial/")
getwd()
S<-subset(S,cell_type=='Epithelial')
DimPlot(S,reduction='umap',label = TRUE)
ggsave('sub_cluster_umap.png',width =7,height =6)
ggsave('sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S,reduction='umap',label = TRUE,split.by = 'group')
ggsave('umap_split_groups.png',width =14,height =6)
ggsave('umap_split_groups.pdf',width =14,height =6)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(999)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50)
ggsave("ElbowPlot.pdf",width = 7,height = 7,bg='white')
ggsave("ElbowPlot.png",width = 7,height = 7,bg='white')
S <- FindNeighbors(S,reduction = "harmony",dims = 1:15)
S <- FindClusters(S,resolution = 0.3)
S <- RunUMAP(S,reduction = "harmony",dims = 1:15)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap',label='TRUE')
ggsave("sub_cluster_umap.pdf",width = 7,height = 7,bg='white')
ggsave("sub_cluster_umap.png",width = 7,height = 7,bg='white')
DimPlot(S,reduction='umap',split.by='group')
ggsave("umap_split_groups.pdf",width = 14,height = 7,bg='white')
ggsave("umap_split_groups.png",width = 14,height = 7,bg='white')
saveRDS(S,"Epithelial.rds")

FeaturePlot(S,features = c('SNHG9', 'MMP7','HSPA1A','SOCS3','SLPI'))
ggsave("Epithelial_marker.pdf",width = 7,height = 7,bg='white')
ggsave("Epithelial_marker.png",width = 7,height = 7,bg='white')
VlnPlot(S,features = c('SNHG9', 'MMP7','HSPA1A','SOCS3','SLPI'))
ggsave("Epithelial_marker_vln.pdf",width = 21,height = 14,bg='white')
ggsave("Epithelial_marker_vln.png",width = 21,height = 14,bg='white')

levels(Idents(S)) <- c('Ciliated Epithelium','Epithelium','Epithelium','Epithelium','Epithelium','Ciliated Epithelium','Ciliated Epithelium','AS Epithelium','Epithelium','Epithelium')
Idents(S)<-factor(Idents(S),levels=c('AS Epithelium','Epithelium','Ciliated Epithelium'))
DimPlot(S,reduction='umap')
ggsave(file="dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')
DimPlot(S,reduction='umap',split.by='group')
ggsave(file="dimplot_umap_celltype_group.png", width=14, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype_group.pdf", width=14, height=6, bg = 'white')
S$cell_type<-Idents(S)

VlnPlot(S,features = c('SNHG9', 'MMP7','HSPA1A','SOCS3','SLPI'),pt.size = 0)
ggsave("AS_marker_vln.pdf",width = 21,height = 14,bg='white')
ggsave("AS_marker_vln.png",width = 21,height = 14,bg='white')

count_df <- table(S$sample, S$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("sample", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_sample.png", width=14, height=7,bg = 'white')
ggsave(file="celltype_proportion_sample.pdf", width=14, height=7,bg = 'white')

count_df <- table(S$group, S$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

metadata <- S@meta.data
write.table(metadata,"all_cells_features.xls",sep="\t")
saveRDS(S,"Epithelium_anno.rds")


S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
levels(S$cell_type)[8]<-"Macrophage"
S@meta.data$sub_cell_type<-S@meta.data$cell_type
S_sub<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Fibroblast/Fibroblast_anno.rds")
S_sub<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Epithelial/Epithelium_anno.rds")
adata_cell_names <- rownames(S_sub@meta.data)
bdata_cell_names <- rownames(S@meta.data)
adata_cell_types <- as.character(S_sub@meta.data$cell_type)
bdata_cell_types <- as.character(S@meta.data$sub_cell_type)
# 匹配行名并替换对应列的值，保留未匹配行的原始值
mapped_indices <- match(bdata_cell_names, adata_cell_names)
matched_cell_types <- adata_cell_types[mapped_indices]
bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
# 将新的细胞类型注释赋值给bdata，并转换回factor类型
S@meta.data$sub_cell_type <- as.factor(bdata_cell_types)
DimPlot(S,reduction='umap',group.by = 'cell_type',raster=FALSE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/dimplot_umap_celltype.png", width=8, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/dimplot_umap_celltype.pdf", width=8, height=7,bg = 'white')
DimPlot(S,reduction='umap',group.by = 'sub_cell_type',raster=FALSE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/dimplot_umap_sub_celltype.png", width=9, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/dimplot_umap_sub_celltype.pdf", width=9, height=7,bg = 'white')

metadata <- S@meta.data
write.table(metadata,"/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/all_cells_features.xls",sep="\t")
levels(S@meta.data$sub_cell_type)
Idents(S)<-S$sub_cell_type
saveRDS(S,"Main_mapping.rds")