.libPaths()
library(Seurat)
library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)

setwd('D:/项目/FWAE20254246/')
getwd()
S<-readRDS('result.rds')
DimPlot(S,reduction='umap',raster=FALSE)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
set.seed(666)
S <- RunHarmony(S, group.by.vars = "sample")

ElbowPlot(S,ndims = 50,reduction = "harmony")
S <- FindNeighbors(S,reduction = "harmony",dims = 1:30)
S <- FindClusters(S,resolution = 0.2)
S <- RunUMAP(S,reduction = "harmony",dims = 1:30)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
bdata_cell_names <- rownames(S@meta.data)
adata_cell_names <- all_cells_features$cell_id
adata_cell_types <- as.character(all_cells_features$cell_type)
bdata_cell_types <- as.character(S@meta.data$cell_type)
mapped_indices <- match(bdata_cell_names, adata_cell_names)
matched_cell_types <- adata_cell_types[mapped_indices]
bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
S@meta.data$cell_type <- as.factor(bdata_cell_types)

DimPlot(S,reduction='umap',raster=FALSE,group.by = 'cell_type')
DimPlot(S,reduction='umap',split.by='group',raster=FALSE)

saveRDS(S,'D:/项目/FWAE20254246/result_harmony_pc30.rds')

S_sub<-S
DimPlot(S_sub,reduction='umap',label = TRUE,raster=FALSE)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/sub_cluster_umap.png',width =7,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S_sub,reduction='umap',label = TRUE,split.by = 'group',raster=FALSE)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/umap_split_groups.png',width =28,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/umap_split_groups.pdf',width =28,height =6)
DimPlot(S_sub,reduction='umap',group.by = 'cell_type',label = TRUE,raster=FALSE)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/cell_type_umap.png',width =7,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/cell_type_umap.pdf',width =7,height =6)

S_diffgene <- FindAllMarkers(object = S_sub, only.pos = TRUE)
write.table(S_diffgene,"./02.Mission/Main_test/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
S_diffgene <- read.table("./02.Mission/Main_test/diffgene.xls",header = TRUE)

top10 <- S_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
all.genes <- rownames(top10$gene)
S_sub <- ScaleData(S_sub, features = all.genes)
DoHeatmap(S_sub, features = top10$gene)+NoLegend()
ggsave(file="./02.Mission/Main_test/top10_diffgene_heatmap.png", width=12, height=16, bg = 'white')
ggsave(file="./02.Mission/Main_test/top10_diffgene_heatmap.pdf", width=12, height=16, bg = 'white')

S_sub$raw_cell_type<-S_sub$cell_type
levels(Idents(S_sub)) <- c('T cells','Myeloid cells','Neutrophil','Hepatocyte','Myeloid cells','NK',
                           'Endothelial cells','Myeloid cells','Plasma cellS','Mast cells','Myeloid cells',
                           'Fibroblast','Hepatocyte','T cells','doublet','T cells','T cells','Hepatocyte',
                           'Myeloid cells','Hepatocyte','Myeloid cells','T cells','doublet','doublet','lowQ')
DimPlot(S_sub,reduction='umap',raster=FALSE)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/annonation_umap.png',width =7,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/annonation_umap.pdf',width =7,height =6)

count_df <- table(S_sub$group, S_sub$cell_type)
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
ggsave(file="./02.Mission/Main_test/celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="./02.Mission/Main_test/celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

metadata <- S_sub@meta.data
write.table(metadata,"./02.Mission/Main_test/all_cells_features.xls",sep="\t")
saveRDS(S_sub,'./02.Mission/Main_test/anno.rds')


#### 去低质量重新注释 ####
S_subf<-subset(S_sub,idents=c('doublet','lowQ'),invert=TRUE)
S_subf <- NormalizeData(S_subf, normalization.method = "LogNormalize", scale.factor = 10000)
S_subf<- FindVariableFeatures(S_subf, selection.method = "vst", nfeatures = 2000)
S_subf <- ScaleData(S_subf)
S_subf <- RunPCA(S_subf, features = VariableFeatures(object = S_subf))
set.seed(666)
S_subf <- RunHarmony(S_subf, group.by.vars = "sample")

ElbowPlot(S_subf,ndims = 50,reduction = "harmony")
S_subf <- FindNeighbors(S_subf,reduction = "harmony",dims = 1:30)
S_subf <- FindClusters(S_subf,resolution = 0.2)
DimPlot(S_subf,reduction='umap',raster=FALSE,group.by = 'cell_type')
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/annonation_umap_filter.png',width =7,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/annonation_umap_filter.pdf',width =7,height =6)
DimPlot(S_subf,reduction='umap',label = TRUE,group.by = 'cell_type',split.by = 'group',raster=FALSE)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/umap_split_groups.png',width =28,height =6)
ggsave('D:/项目/FWAE20254246/02.Mission/Main_test/umap_split_groups.pdf',width =28,height =6)

count_df <- table(S_subf$group, S_subf$cell_type)
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
ggsave(file="./02.Mission/Main_test/celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="./02.Mission/Main_test/celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

metadata <- S_subf@meta.data
write.table(metadata,"./02.Mission/Main_test/all_cells_features.xls",sep="\t")

S_subf@meta.data <- S_subf@meta.data[, !names(S_subf@meta.data) %in% c("RNA_snn_res.0.5","cluster","leiden","n_genes_by_counts")]
S_sub@meta.data$cell_type <- Idents(S_sub)
adata_cell_names <- row.names(S_sub@meta.data)
bdata_cell_names <- row.names(S_subf@meta.data)
adata_cell_types <- as.character(S_sub@meta.data$cell_type)
bdata_cell_types <- as.character(S_subf@meta.data$cell_type)
mapped_indices <- match(bdata_cell_names, adata_cell_names)
matched_cell_types <- adata_cell_types[mapped_indices]
bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
S_subf@meta.data$cell_type <- bdata_cell_types
S_subf@meta.data$cell_type <- as.factor(S_subf@meta.data$cell_type)
Idents(S_subf)<-S_subf@meta.data$cell_type
saveRDS(S_subf,'./02.Mission/Main_test/annof.rds')

#### 基因集评分 ####
geneset<-read.gmt("02.Mission/genset-9.gmt")
geneset<-geneset$gene
gene<-list(geneset)
adata<-AddModuleScore_UCell(S_subf,features=gene)
colnames(adata@meta.data)
colnames(adata@meta.data)[23]<-'select_geneset_Score'
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)
ggsave("./02.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("./02.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

#### 定义HCC ####
S_subf<-readRDS("02.Mission/Main_test/annof.rds")
S_sub<-subset(S_subf,idents=c('Hepatocyte'))

S_sub <- NormalizeData(S_sub, normalization.method = "LogNormalize", scale.factor = 10000)
S_sub <- FindVariableFeatures(S_sub, selection.method = "vst", nfeatures = 2000)
S_sub <- ScaleData(S_sub)
S_sub <- RunPCA(S_sub, features = VariableFeatures(object = S_sub))
set.seed(999)
S_sub <- RunHarmony(S_sub, group.by.vars = "sample")

ElbowPlot(S_sub,ndims = 50,reduction = 'harmony')
S_sub <- FindNeighbors(S_sub,reduction = "harmony",dims = 1:5)
S_sub <- FindClusters(S_sub,resolution = 0.3)
S_sub <- RunUMAP(S_sub,reduction = "harmony",dims = 1:5,metric="euclidean")
# S <- RunTSNE(S_sub,reduction = "harmony",dims = 1:30)
DimPlot(S_sub,reduction='umap',label = TRUE)
ggsave('./02.Mission/HCC/sub_cluster_umap.png',width =7,height =6)
ggsave('./02.Mission/HCC/sub_cluster_umap.pdf',width =7,height =6)
DimPlot(S_sub,reduction='umap',label = TRUE,split.by = 'group')
ggsave('./02.Mission/HCC/umap_split_groups.png',width =28,height =6)
ggsave('./02.Mission/HCC/umap_split_groups.pdf',width =28,height =6)

S_sub_test<-subset(S_sub,group %in% c("LRT","PHT"))
DimPlot(S_sub_test,reduction='umap',label = TRUE,split.by = 'group')
ggsave('./02.Mission/HCC/umap_split_groups.png',width =14,height =6)
ggsave('./02.Mission/HCC/umap_split_groups.pdf',width =14,height =6)

#### HCC基因集评分 ####
geneset9<-read.gmt("02.Mission/genset-9.gmt")
geneset<-geneset9$gene
gene<-list(geneset)

S_sub<-AddModuleScore_UCell(S_sub,features=gene)
colnames(S_sub@meta.data)
colnames(S_sub@meta.data)[22]<-'genset-9_UCell'
p <- VlnPlot(S_sub,features='genset-9_UCell',group.by = 'RNA_snn_res.0.3',pt.size = 0)
ggsave("./02.Mission/HCC/genset-9_VlnPlot.pdf",width = 16,height = 10)
ggsave("./02.Mission/HCC/genset-9_VlnPlot.png",width = 16,height = 10)
FeaturePlot(S_sub,features = c('genset-9_UCell'),reduction='umap',label = TRUE,raster=FALSE,min.cutoff = 0.1)
ggsave('./02.Mission/HCC/genset-9_umap.png',width =7,height =6)
ggsave('./02.Mission/HCC/genset-9_umap.pdf',width =7,height =6)

FeaturePlot(S_sub,features = c('ALB','CYP2C9',"AFP","KRT8",'APOA1',"APOA2"),reduction='umap',raster=FALSE)
ggsave('./02.Mission/HCC/Hepatocyte_marker.png',width =14,height =18)
ggsave('./02.Mission/HCC/Hepatocyte_marker.pdf',width =14,height =18)

geneset13<-read.gmt("02.Mission/genset-13.gmt")
geneset<-geneset13$gene
gene<-list(geneset)
S_sub<-AddModuleScore_UCell(S_sub,features=gene)
colnames(S_sub@meta.data)
colnames(S_sub@meta.data)[23]<-'genset-13_UCell'
p <- VlnPlot(S_sub,features='genset-13_UCell',group.by = 'RNA_snn_res.0.3',pt.size = 0)
ggsave("./02.Mission/HCC/genset-13_VlnPlot.pdf",width = 16,height = 10)
ggsave("./02.Mission/HCC/genset-13_VlnPlot.png",width = 16,height = 10)
FeaturePlot(S_sub,features = 'genset-13_UCell',reduction='umap',label = TRUE,raster=FALSE)
ggsave('./02.Mission/HCC/genset-13_umap.png',width =7,height =6)
ggsave('./02.Mission/HCC/genset-13_umap.pdf',width =7,height =6)


geneset9<-geneset9$gene
geneset13<-geneset13$gene
geneset<-intersect(geneset9, geneset13)
geneset9<-setdiff(geneset9,geneset)
geneset13<-setdiff(geneset13,geneset)

gene<-list(geneset9,geneset13)
S_sub<-AddModuleScore_UCell(S_sub,features=gene)
colnames(S_sub@meta.data)
colnames(S_sub@meta.data)[24]<-'genset-9-diff_UCell'
colnames(S_sub@meta.data)[25]<-'genset-13-diff_UCell'
p <- VlnPlot(S_sub,features=c('genset-9-diff_UCell','genset-13-diff_UCell'),group.by = 'RNA_snn_res.0.3',pt.size = 0)
ggsave("./02.Mission/HCC/genset-diff_VlnPlot.pdf",width = 24,height = 10)
ggsave("./02.Mission/HCC/genset-diff_VlnPlot.png",width = 24,height = 10)
FeaturePlot(S_sub,features = c('genset-9-diff_UCell','genset-13-diff_UCell'),reduction='umap',label = TRUE,raster=FALSE,min.cutoff = 0.1)
ggsave('./02.Mission/HCC/genset-diff_umap.png',width =14,height =6)
ggsave('./02.Mission/HCC/genset-diff_umap.pdf',width =14,height =6)
metadata <- S_sub@meta.data
write.table(metadata,"./02.Mission/HCC/all_cells_features.xls",sep="\t")
saveRDS(S_sub,"./02.Mission/HCC/HCC.rds")

#### 差异基因 ####
S_diffgene <- FindAllMarkers(object = S_sub, only.pos = TRUE)
write.table(S_diffgene,"./02.Mission/HCC/diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
S_diffgene <- read.table("./02.Mission/HCC/diffgene.xls",header = TRUE)

top10 <- S_diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
all.genes <- rownames(top10$gene)
S_sub <- ScaleData(S_sub, features = all.genes)
DoHeatmap(S_sub, features = top10$gene,raster = FALSE)+NoLegend()
ggsave(file="./02.Mission/HCC/top10_diffgene_heatmap.png", width=12, height=16, bg = 'white')
ggsave(file="./02.Mission/HCC/top10_diffgene_heatmap.pdf", width=12, height=16, bg = 'white')

mean_exp <- AverageExpression(S_sub,features = top10$gene,group.by = "RNA_snn_res.0.3",slot = "data")[[1]]
library(pheatmap)
p<-pheatmap(t(mean_exp),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "column",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         angle_col = 45,
         border_color = NA,
         color = colorRampPalette(c("#3D52A4", "white", "#D13830"))(100))
ggsave(file="./02.Mission/HCC/top10_diffgene_heatmap_cluster.png", p, width=28, height=8, bg = 'white')
ggsave(file="./02.Mission/HCC/top10_diffgene_heatmap_cluster.pdf", p, width=28, height=8, bg = 'white')

levels(Idents(S_sub)) <- c('HCC_C1','HCC_C2','HCC_C1','HCC_C2','HCC_C2','HCC_C1','HCC_C3','HCC_C3',
                           'HCC_C4','HCC_C5','HCC_C4','HCC_C2','HCC_C1','HCC_C1','HCC_C1')
S_sub$sub_cell_type<-Idents(S_sub)
DimPlot(S_sub,reduction='umap')
ggsave(file="./02.Mission/HCC/dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="./02.Mission/HCC/dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')

metadata <- S_sub@meta.data
write.table(metadata,"./02.Mission/HCC/all_cells_features.xls",sep="\t")
DimPlot(S_sub,reduction='umap',split.by = 'group')
ggsave(file="./02.Mission/HCC/dimplot_umap_group_celltype.png", width=28, height=6, bg = 'white')
ggsave(file="./02.Mission/HCC/dimplot_umap_group_celltype.pdf", width=28, height=6, bg = 'white')

count_df <- table(S_sub$group, S_sub$sub_cell_type)
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
ggsave(file="./02.Mission/HCC/celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="./02.Mission/HCC/celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

saveRDS(S_sub,'./02.Mission/HCC/recluster_HCC.rds')