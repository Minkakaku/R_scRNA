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

PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/01.Mission/all_data_sc_doublet.rds")
PRO <- subset(PRO,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY")))
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
DefaultAssay(PRO) <-'RNA'
marker_file <- read.table("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/script/markergene.txt",header = FALSE,sep="\t",stringsAsFactors=FALSE)
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))

marker <- unlist(strsplit(marker_file$V1, ","))
marker <- unique(marker)

print("gene not find: ")
print(marker[!(marker %in% rownames(PRO))])
marker<-intersect(marker,rownames(PRO))
outumap<-"/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/UMAP/"
outvln<-"/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/Vln/"
for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('#072727','#195e6c','#eb2d4a','#f1681e','#ffeee5'),pt.size = pt_use)
    pdf(paste0(outumap, i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outumap, i, '_Featureplot.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outvln, i, '_VlnPlot.pdf'),width = 10,height = 10)
		ggsave(paste0(outvln, i, '_VlnPlot.png'),width = 10,height = 10)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Main",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Main'),outdir=outvln,map_sig_level=T,only_sig = T)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana.pdf'),width = 10,height = 20)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana.png'),width = 10,height = 20)
}

PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/02.Mission/sc_neut.rds")
PRO <- subset(PRO,subset=(sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS')))
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
DefaultAssay(PRO) <-'RNA'
for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('#072727','#195e6c','#eb2d4a','#f1681e','#ffeee5'),pt.size = pt_use)
    pdf(paste0(outumap, i, '_Featureplot_Neutrophil.pdf'))
    print(p)
    dev.off()
    png(paste0(outumap, i, '_Featureplot_Neutrophil.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outvln, i, '_VlnPlot_Neutrophil.pdf'),width = 7,height = 6)
		ggsave(paste0(outvln, i, '_VlnPlot_Neutrophil.png'),width = 7,height = 6)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Neutrophil",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Neutrophil'),outdir=outvln,map_sig_level=T,only_sig = T)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana_Neutrophil.pdf'),width = 7,height = 8)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana_Neutrophil.png'),width = 7,height = 8)
}


PRO<-readRDS("/data/NAS02/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds")
PRO <- subset(PRO,subset=(sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS')))
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
DefaultAssay(PRO) <-'RNA'
for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('#072727','#195e6c','#eb2d4a','#f1681e','#ffeee5'),pt.size = pt_use)
    pdf(paste0(outumap, i, '_Featureplot_Monocyte.pdf'))
    print(p)
    dev.off()
    png(paste0(outumap, i, '_Featureplot_Monocyte.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outvln, i, '_VlnPlot_Monocyte.pdf'),width = 7,height = 6)
		ggsave(paste0(outvln, i, '_VlnPlot_Monocyte.png'),width = 7,height = 6)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Monocyte",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Monocyte'),outdir=outvln,map_sig_level=T,only_sig = T)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana_Monocyte.pdf'),width = 7,height = 8)
		ggsave(paste0(outvln, i, '_VlnPlot_sig_ana_Monocyte.png'),width = 7,height = 8)
}