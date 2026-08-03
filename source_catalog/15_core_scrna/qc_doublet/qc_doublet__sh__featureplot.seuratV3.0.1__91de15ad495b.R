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

# marker_file <- read.table("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/script/markergene.txt",header = FALSE,sep="\t",stringsAsFactors=FALSE)
marker <- c("HIF1A","SLC2A3","PKM","LDHA","ATP6V0B",'ATP6V0C','ENO1')
PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/all_data_sc_doublet.rds")
PRO$sample<- ifelse(PRO$sample == "YLY", "YLY0701", PRO$sample)
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
PRO$sample<-factor(PRO$sample)
# PRO$sample<-factor(PRO$sample,levels=c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))

cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/marker_group/"
for (i in marker){
		# p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'sample')+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)+ scale_fill_manual(values = c("#be348f","#be348f","#be348f","#be348f","#197ec0", "#197EC0", "#197ec0","#999999","#999999"))
		# pdf(paste0(outdir, 'mono-', i, '_VlnPlot.pdf'),width = 11)
    # print(p)
    # dev.off()
    # png(paste0(outdir, 'mono-', i, '_VlnPlot.png'),width = 11/7*480)
    # print(p)
    # dev.off()

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'Main-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE,scale.min = 0.5, scale.max = 0.5) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Main-group_DotPlot.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Main-group_DotPlot.png'),width = 7,height = 7,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.png'),width = 7,height = 7,bg = 'white')

PRO <- subset(PRO,annotation %in% "Monocyte")
for (i in marker){
		# p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'sample')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)+ scale_fill_manual(values = c("#be348f","#be348f","#be348f","#be348f","#197ec0", "#197EC0", "#197ec0","#999999","#999999"))
		# pdf(paste0(outdir, 'CD14mono-', i, '_VlnPlot.pdf'),width = 11)
    # print(p)
    # dev.off()
    # png(paste0(outdir, 'CD14mono-', i, '_VlnPlot.png'),width = 11/7*480)
    # print(p)
    # dev.off()

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, group.by = 'group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,'Monocyte-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE,scale.min = 0.5, scale.max = 0.5) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Monocyte-group_DotPlot.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Monocyte-group_DotPlot.png'),width = 7,height = 7,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.png'),width = 7,height = 7,bg = 'white')


marker_file <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
marker <- marker_file$HALLMARK_GLYCOLYSIS
PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/all_data_sc_doublet.rds")
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/GLYCOLYSIS_gene/"

marker<-intersect(marker, rownames(PRO))
average_expression <- AverageExpression(PRO, features = marker, assays = 'RNA3', group.by = "group")
write.table(average_expression$RNA3,paste0(outdir,"Main_average_expression.xls"),sep="\t",row.names=TRUE)
result <- data.frame(average_expression$RNA3) %>% filter(LAHS > NLAHS & NLAHS > NC)
marker<-rownames(result)
Idents(PRO)<-PRO$group
p<-VlnPlot(PRO, features = marker, 
             stack = TRUE, 
             sort = FALSE, 
             flip = FALSE)+theme(legend.position = "none")
ggsave(paste0(outdir,'Main-group_stackvlnPlot.pdf'),width = 25,height = 7,bg = 'white')
ggsave(paste0(outdir,'Main-group_stackvlnPlot.png'),width = 25,height = 7,bg = 'white')

marker<-c("MIF","ZNF292","UGP2","TXN","HDLBP","G6PD","PGK1",'ENO1','HSPA5','PKM','TALDO1','LDHA','VCAN')
p<-VlnPlot(PRO, features = marker, 
             stack = TRUE, 
             sort = FALSE, 
             flip = TRUE)+theme(legend.position = "none")+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,'Main-group_stackvlnPlot_select.pdf'),width = 7,height = 26,bg = 'white')
ggsave(paste0(outdir,'Main-group_stackvlnPlot_select.png'),width = 7,height = 26,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE,scale.min = 0.5, scale.max = 0.5) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Main-group_DotPlot.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Main-group_DotPlot.png'),width = 7,height = 7,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Main-group_DotPlot_raw.png'),width = 7,height = 7,bg = 'white')

PRO <- subset(PRO,annotation %in% "Monocyte")
marker_file <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
marker <- marker_file$HALLMARK_GLYCOLYSIS
marker<-intersect(marker, rownames(PRO))
average_expression <- AverageExpression(PRO, features = marker, assays = 'RNA3', group.by = "group")
write.table(average_expression$RNA3,paste0(outdir,"Monocyte_average_expression.xls"),sep="\t",row.names=TRUE)
result <- data.frame(average_expression$RNA3) %>% filter(LAHS > NLAHS & NLAHS > NC)
marker<-rownames(result)

p<-VlnPlot(PRO, features = marker, 
             stack = TRUE, 
             sort = FALSE, 
             flip = FALSE)+theme(legend.position = "none")
ggsave(paste0(outdir,'Monocyte-group_stackvlnPlot.pdf'),width = 25,height = 7,bg = 'white')
ggsave(paste0(outdir,'Monocyte-group_stackvlnPlot.png'),width = 25,height = 7,bg = 'white')

marker<-c("TXN","SLC16A3","PSMC4","PPIA","PKM","PGK1",'PGAM1','NASP','MIF','LDHA','IRS2','IDH1','HAX1','GUSB','ENO1')
p<-VlnPlot(PRO, features = marker, 
             stack = TRUE, 
             sort = FALSE, 
             flip = TRUE)+theme(legend.position = "none")+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,'Monocyte-group_stackvlnPlot_select.pdf'),width = 7,height = 26,bg = 'white')
ggsave(paste0(outdir,'Monocyte-group_stackvlnPlot_select.png'),width = 7,height = 26,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15,scale=TRUE,scale.min = 0.5, scale.max = 0.5) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Monocyte-group_DotPlot.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Monocyte-group_DotPlot.png'),width = 7,height = 7,bg = 'white')

p<-DotPlot(PRO, features = marker,group.by = 'group',dot.scale = 15) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.pdf'),width = 7,height = 7,bg = 'white')
		ggsave(paste0(outdir,'Monocyte-group_DotPlot_raw.png'),width = 7,height = 7,bg = 'white')
