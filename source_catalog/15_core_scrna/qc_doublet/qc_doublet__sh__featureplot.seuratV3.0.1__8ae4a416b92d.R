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
# library(scCustomize)
})
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

marker_file <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
marker <- marker_file$HALLMARK_GLYCOLYSIS
PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/all_data_sc_doublet.rds")
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/18.Mission/GLYCOLYSIS_gene/"

Idents(PRO)<-PRO$group
marker<-c("MIF","ZNF292","UGP2","TXN","HDLBP","G6PD","PGK1",'ENO1','HSPA5','PKM','TALDO1','LDHA','VCAN')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色

result <- data.frame(p$data) %>% filter(id=='LAHS') %>% arrange(pct.exp)
marker <- rownames(result)
p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.pdf'),width = 7,height = 8,bg = 'white')
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.png'),width = 7,height = 8,bg = 'white')



PRO <- subset(PRO,annotation %in% "Monocyte")
marker<-c("TXN","SLC16A3","PSMC4","PPIA","PKM","PGK1",'PGAM1','NASP','MIF','LDHA','IRS2','IDH1','HAX1','GUSB','ENO1')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色

result <- data.frame(p$data) %>% filter(id=='LAHS') %>% arrange(pct.exp)
marker <- rownames(result)

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.pdf'),width = 7,height = 8,bg = 'white')
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.png'),width = 7,height = 8,bg = 'white')



marker <- c("HIF1A","SLC2A3","PKM","LDHA","ATP6V0B",'ATP6V0C','ENO1')
PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_anno.rds")
PRO$sample<- ifelse(PRO$sample == "YLY", "YLY0701", PRO$sample)
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
PRO$sample<-factor(PRO$sample)
# PRO$sample<-factor(PRO$sample,levels=c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/18.Mission/marker_group/"
for (i in marker){
		p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}


PRO <- subset(PRO,annotation %in% "Monocyte")
for (i in marker){
		p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}
