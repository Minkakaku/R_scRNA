#!/usr/bin/env Rscript
# conda activate diffgene
suppressMessages({
library(ggplot2)
library(reshape2)
library(Seurat)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
library(dplyr)
library(ggpubr)
library(ggbreak)
library(pheatmap)
# library(scCustomize)
})
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
# PRO$sample<- ifelse(PRO$sample == "YLY", "YLY0701", PRO$sample)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/"
marker <- c("THBS1")
PRO_sub <- subset(PRO,cell_type %in% "PTTG1 classical monocyte")
for (i in marker){
		p = VlnPlot(object = PRO_sub, features = i,pt.size = 0, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'PTTG1 classical monocyte-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'PTTG1 classical monocyte-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T,method = "t.test")
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'PTTG1 classical monocyte-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'PTTG1 classical monocyte-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}

marker <- c("CD47")
PRO_sub <- subset(PRO,cell_type %in% "Classical monocyte") #Dendritic cell,Non-classical monocyte,Classical monocyte
for (i in marker){
		p = VlnPlot(object = PRO_sub, features = i,pt.size = 0, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Classical monocyte-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'Classical monocyte-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Classical monocyte-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'Classical monocyte-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}



gene<-list(c("THBS1","CD47","CD36"))
# ann_colors=list(cell_type=c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"),group=c(Before='#1f77b4',Post='#ff7f0e'))
ann_colors=list(cell_type=c("Non-classical monocyte"='#f8766d',"PTTG1 classical monocyte"='#b79f00',"ZFP36L1 classical monocyte"='#00ba38',"S100A8 classical monocyte"="#00bfc4",'Classical monocyte'='#619cff','Dendritic cell'='#f564e3'),group=c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
# data_matrix <- FetchData(adata,vars = c(gene$x, "cell_type", "group"),slot = "data")
data_matrix <- FetchData(PRO,vars = c(gene[[1]], "cell_type", "group"),slot = "data")
data_table <- data_matrix %>%   group_by(cell_type,group) %>%  summarise(across(everything(), mean, na.rm = TRUE)) %>% as.data.frame()
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/Monocyte-cell_type-group.xls",sep="\t",row.names=FALSE)
data_all<-as.matrix(data_table[,3:5])
rownames(data_all)<-paste0("sample-",rownames(data_table))
# 绘制热图
anno_row <- data.frame(data_table[,1:2])
rownames(anno_row)<-rownames(data_all)
# 绘制热图
p2 <- pheatmap(data_all,  
               color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
               border_color = "black", 
							#  border=F, 
               scale = "none", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               legend_labels = c("low","","heigh"), 
               show_rownames = FALSE, 
               show_colnames = TRUE, 
               fontsize = 8, 
               annotation_row = anno_row,   #添加列注释信息
               annotation_legend = TRUE,    #是否显示注释信息图例
               annotation_names_row = FALSE,  #是否显示列注释的名称
               annotation_colors=ann_colors
							 )
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/Monocyte-cell_type-group-Heatmap.pdf",width = 6,height = 7.5)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/Monocyte-cell_type-group-Heatmap.png",width = 6*60,height = 7.5*60)
p2
dev.off()

data_all<-as.matrix(data_table[,3:4])
rownames(data_all)<-paste0("sample-",rownames(data_table))
p2 <- pheatmap(data_all,  
               color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
               border_color = "black", 
							#  border=F, 
               scale = "none", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               legend_labels = c("low","","heigh"), 
               show_rownames = FALSE, 
               show_colnames = TRUE, 
               fontsize = 8, 
               annotation_row = anno_row,   #添加列注释信息
               annotation_legend = TRUE,    #是否显示注释信息图例
               annotation_names_row = FALSE,  #是否显示列注释的名称
               annotation_colors=ann_colors
							 )
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/Monocyte-cell_type-group-Heatmap_2.pdf",width = 6,height = 7.5)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/marker_group/Monocyte-cell_type-group-Heatmap_2.png",width = 6*60,height = 7.5*60)
p2
dev.off()