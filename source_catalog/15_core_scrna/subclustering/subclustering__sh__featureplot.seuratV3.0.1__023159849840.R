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

# PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/01.Mission/all_data_sc_doublet.rds")
# PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_cd8t.rds")
PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_cd4t.rds")
prefix <- "cd4t"
PRO[["RNA3"]] <- as(object = PRO[["RNA"]], Class = "Assay")
DefaultAssay(PRO) <- 'RNA3'
# PRO <- subset(PRO,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY")))
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
DefaultAssay(PRO) <-'RNA'
# PRO<-subset(PRO,idents=c("CD4T","CD8T","T proliferating"))
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
marker <- "HIF1A"
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/09.Mission/HIF1A/"

for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('lightgrey','red'),pt.size = pt_use)
    pdf(paste0(outdir,prefix,'-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,prefix,'-', i, '_Featureplot.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot.png'),width = 7,height = 7)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix=prefix,map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,'_',prefix),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot.pdf'),width = 15,height = 10)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot.png'),width = 15,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),map_sig_level=T,only_sig = T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot_sig_ana.pdf'),width = 15,height = 15)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot_sig_ana.png'),width = 15,height = 15)
}


PRO<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds")
PRO[["RNA3"]] <- as(object = PRO[["RNA"]], Class = "Assay")
DefaultAssay(PRO) <- 'RNA3'
PRO <- subset(PRO,subset=sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))
PRO$group<-factor(PRO$group)
PRO$group<-factor(PRO$group,levels=c("LAHS","NLAHS","NC"))
PRO$sample<-factor(PRO$sample)
PRO$sample<-factor(PRO$sample,levels=c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HC01','HC02'))

cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
marker <- "IFNGR1"
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/09.Mission/IFNGR/"

for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('lightgrey','red'),pt.size = pt_use)
    pdf(paste0(outdir,prefix,'-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,prefix,'-', i, '_Featureplot.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot.png'),width = 7,height = 7)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix=prefix,map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,'_',prefix),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir,prefix,'-', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot.pdf'),width = 15,height = 10)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot.png'),width = 15,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),map_sig_level=T,only_sig = T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot_sig_ana.pdf'),width = 15,height = 15)
		ggsave(paste0(outdir,prefix,'-group', i, '_VlnPlot_sig_ana.png'),width = 15,height = 15)
}