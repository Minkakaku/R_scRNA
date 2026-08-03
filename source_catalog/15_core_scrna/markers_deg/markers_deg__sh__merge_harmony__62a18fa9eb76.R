library(Seurat)
library(Matrix)
library(ggplot2)
library(magrittr)
library(dplyr)
library(scales)
library(RColorBrewer)

S_sub<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/02.Mission/Main/Main.rds")
DimPlot(S_sub,reduction='umap',group.by = 'cell_type')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_proportion/dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_proportion/dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')

count_df <- table(S_sub$group, S_sub$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = hue_pal()(10)
# "#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"
# final_palette = colorRampPalette(my_colors)(length(unique(proportion_df$cell_type)))
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_proportion/celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_proportion/celltype_proportion_group.pdf", width=7, height=7,bg = 'white')



S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Fibroblast/Fibroblast.rds")
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Fibroblast/")
DimPlot(S,reduction='umap',group.by='seurat_clusters')
ggsave(file="dimplot_umap_celltype.png", width=7, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype.pdf", width=7, height=6, bg = 'white')
DimPlot(S,reduction='umap',group.by='seurat_clusters',split.by='group')
ggsave(file="dimplot_umap_celltype_group.png", width=14, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype_group.pdf", width=14, height=6, bg = 'white')

levels(Idents(S)) <- c('FABP5+Fibro','MTRNR2L8+Fibro','COL1A2+Fibro','LAMA2+Fibro','HBB+Fibro','S100A8+Fibro')
S$cell_type<-Idents(S)
DimPlot(S,reduction='umap')
ggsave(file="dimplot_umap_celltype.png", width=8, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype.pdf", width=8, height=6, bg = 'white')
DimPlot(S,reduction='umap',split.by='group')
ggsave(file="dimplot_umap_celltype_group.png", width=14, height=6, bg = 'white')
ggsave(file="dimplot_umap_celltype_group.pdf", width=14, height=6, bg = 'white')

rd_bu_11 <- colorRampPalette(brewer.pal(11, "RdBu"))(11)
DotPlot(S,features = c('COL1A1','COL3A1','ACTA2','FN1','SMAD3','CTNNB1','MYC','LRP6','TGFB1','TGFBR1','SPTBN'),dot.scale = 10) + scale_colour_gradientn(colours = rev(rd_bu_11)) + scale_fill_gradientn(colours   = rev(rd_bu_11)) + coord_flip()
ggsave(file="DotPlot_celltype.png", width=12, height=6, bg = 'white')
ggsave(file="DotPlot_celltype.pdf", width=12, height=6, bg = 'white')

count_df <- table(S$group, S$cell_type)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
my_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
final_palette = hue_pal()(6)
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

saveRDS(S,"Fibroblast.rds")


# 巨噬细胞细分，基本全是M1
setwd("/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Macrophage/")
getwd()
S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Monocyte/Macrophage.rds")
DimPlot(S,reduction='umap')
ggsave("umap_mapping.pdf",width = 8,height = 7)
ggsave("umap_mapping.png",width = 8,height = 7)

S <- FindClusters(S,resolution = 0.5)
S <- RunUMAP(S,reduction = "harmony",dims = 1:10)
# S <- RunTSNE(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction='umap',pt.size=3)
ggsave("sub_cluster_umap.pdf",width = 7,height = 7,bg='white')
ggsave("sub_cluster_umap.png",width = 7,height = 7,bg='white')
DimPlot(S,reduction='umap',split.by='group',pt.size=3)
ggsave("umap_split_groups.pdf",width = 14,height = 7,bg='white')
ggsave("umap_split_groups.png",width = 14,height = 7,bg='white')

select1<-c('ATG101','ATG12','ATG7','BECN1','MAP1LC3B','MTOR','PIK3C3','PRKAA1','RB1CC1','WIPI2')
DotPlot(S,features = select1,group.by='group')
ggsave(file="geneset1_dotplot.png", width=12, height=6,bg = 'white')
ggsave(file="geneset1_dotplot.pdf", width=12, height=6,bg = 'white')
select2<-c('AKT2','AKT3','CSF3R','CXCR4','IGF1R','IL6ST','JAK1','JAK2','JAK3','MAPK1','MAPK2K2','PIK3R1','PTPN11','SOCS2','SOCS4','STAT1','STAT2','STAT3','STAT4','STAT5A','STAT6')
DotPlot(S,features = select2,group.by='group')
ggsave(file="geneset2_dotplot.png", width=20, height=6,bg = 'white')
ggsave(file="geneset2_dotplot.pdf", width=20, height=6,bg = 'white')

S_diffgene <- FindAllMarkers(object = S)
write.table(S_diffgene,"diffgene.xls",sep="\t",quote=F,row.names=F,col.names=T)
top3 <- S_diffgene %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
all.genes <- unique(top3$gene)
DotPlot(S,features = all.genes)
ggsave(file="top3gene_dotplot.png", width=18, height=7,bg = 'white')
ggsave(file="top3gene_dotplot.pdf", width=18, height=7,bg = 'white')

count_df <- table(S$group, S$seurat_clusters)
proportion_df <- prop.table(count_df, margin = 1)
proportion_df <- as.data.frame(proportion_df)
colnames(proportion_df) <- c("group", "cell_type", "proportion")
final_palette = hue_pal()(5)
ggplot(proportion_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col(color = "black", size = 0.2) +  # 添加边框，颜色为黑色，线宽为 0.2
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 3) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = final_palette)+
  theme_minimal()
ggsave(file="celltype_proportion_group.png", width=7, height=7,bg = 'white')
ggsave(file="celltype_proportion_group.pdf", width=7, height=7,bg = 'white')

# 2. 按常用阈值过滤，并增加一列方向
markers.filt <- S_diffgene %>%
  filter(p_val_adj < 0.05) %>%               # 显著性
  mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))
# 3. 统计每类每方向的基因数
plot.dat <- markers.filt %>%
  count(cluster, direction) %>%
  mutate(direction = factor(direction, levels = c("Up", "Down")))
# 4. 堆积柱状图
ggplot(plot.dat, aes(x = cluster, y = n, fill = direction)) +
  geom_col(width = 0.7, color = "black", size = 0.25) +
  scale_fill_manual(values = c("Down" = "#377eb8", "Up" = "#e41a1c")) +
  labs(x = "Cell type", y = "Number of DEGs") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    # 去掉主要和次要网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加边框
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
ggsave(file="diffgene_barplot.png", width=10, height=7,bg = 'white')
ggsave(file="diffgene_barplot.pdf", width=10, height=7,bg = 'white')
saveRDS(S,'Macrophage.rds')


S<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/03.Mission/Main_mapping/Main_mapping.rds")

S_sub<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Fibroblast/Fibroblast.rds")
S_sub<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Macrophage/Macrophage.rds")
S_sub@meta.data$cell_type<-S_sub@meta.data$seurat_clusters
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
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_mapping/dimplot_umap_celltype.png", width=8, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_mapping/dimplot_umap_celltype.pdf", width=8, height=7,bg = 'white')
DimPlot(S,reduction='umap',group.by = 'sub_cell_type',raster=FALSE)
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_mapping/dimplot_umap_sub_celltype.png", width=9, height=7,bg = 'white')
ggsave(file="/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_mapping/dimplot_umap_sub_celltype.pdf", width=9, height=7,bg = 'white')

metadata <- S@meta.data
write.table(metadata,"/mnt/NexData_file/home/sh/Projects/FWAE20252866/04.Mission/Main_mapping/all_cells_features.xls",sep="\t")
levels(S@meta.data$sub_cell_type)
S$raw_cell_type<-S$cell_type
S$cell_type<-S$sub_cell_type
Idents(S)<-S$sub_cell_type
saveRDS(S,"Main_mapping.rds")