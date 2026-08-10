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

PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/02.Mission/sc_cd8t.rds")
# PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/02.Mission/sc_cd4t.rds")
PRO[["RNA3"]] <- as(object = PRO[["RNA"]], Class = "Assay")
DefaultAssay(PRO) <- 'RNA3'
# PRO <- subset(PRO,subset=sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS'))
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)

# marker_file <- read.table("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/script/markergene.txt",header = FALSE,sep="\t",stringsAsFactors=FALSE)
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))

# marker <- unlist(strsplit(marker_file$V1, ","))
# marker <- unique(marker)
marker <- "IFNG"
# print("gene not find: ")
# print(marker[!(marker %in% rownames(PRO))])
# marker<-intersect(marker,rownames(PRO))
outdir<-"/data/NAS02/sh/Projects/FWSC2024-192/08.Mission/IFNG/"

for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('#072727','#195e6c','#eb2d4a','#f1681e','#ffeee5'),pt.size = pt_use)
    pdf(paste0(outdir,'cd8t-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,'cd8t-', i, '_Featureplot.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'cd8t-', i, '_VlnPlot.pdf'),width = 10,height = 10)
		ggsave(paste0(outdir,'cd8t-', i, '_VlnPlot.png'),width = 10,height = 10)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Main",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Main'),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'cd8t-', i, '_VlnPlot_sig_ana.pdf'),width = 10,height = 20)
		ggsave(paste0(outdir,'cd8t-', i, '_VlnPlot_sig_ana.png'),width = 10,height = 20)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'cd8t-group', i, '_VlnPlot.pdf'),width = 15,height = 10)
		ggsave(paste0(outdir,'cd8t-group', i, '_VlnPlot.png'),width = 15,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix="cd8t_group",map_sig_level=T,only_sig = T)
		get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix="cd8t_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'cd8t-group', i, '_VlnPlot_sig_ana.pdf'),width = 15,height = 20)
		ggsave(paste0(outdir,'cd8t-group', i, '_VlnPlot_sig_ana.png'),width = 15,height = 20)

}



marker <- c("HIF1A","ANXA2","PKM","LDHA")
PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds")
PRO[["RNA3"]] <- as(object = PRO[["RNA"]], Class = "Assay")
DefaultAssay(PRO) <- 'RNA3'
PRO <- subset(PRO,subset=sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
PRO$sample<-factor(PRO$sample)
PRO$sample<-factor(PRO$sample,levels=c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))

outdir<-"/data/NAS02/sh/Projects/FWSC2024-192/08.Mission/marker_group/"
for (i in marker){
		p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'sample')+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)+ scale_fill_manual(values = c("#be348f","#be348f","#be348f","#be348f","#197ec0", "#197EC0", "#197ec0","#999999","#999999"))
		pdf(paste0(outdir, 'mono-', i, '_VlnPlot.pdf'),width = 11)
    print(p)
    dev.off()
    png(paste0(outdir, 'mono-', i, '_VlnPlot.png'),width = 11/7*480)
    print(p)
    dev.off()

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, group.by = 'group')+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 12)
		ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 12)
}

PRO <- subset(PRO,ident="CD14 Monocyte")
for (i in marker){
		p = VlnPlot(object = PRO, features = i,pt.size = 0, group.by = 'sample')+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)+ scale_fill_manual(values = c("#be348f","#be348f","#be348f","#be348f","#197ec0", "#197EC0", "#197ec0","#999999","#999999"))
		pdf(paste0(outdir, 'CD14mono-', i, '_VlnPlot.pdf'),width = 11)
    print(p)
    dev.off()
    png(paste0(outdir, 'CD14mono-', i, '_VlnPlot.png'),width = 11/7*480)
    print(p)
    dev.off()

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, group.by = 'group')+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'CD14mono-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'CD14mono-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'CD14mono-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 12)
		ggsave(paste0(outdir,'CD14mono-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 12)
}