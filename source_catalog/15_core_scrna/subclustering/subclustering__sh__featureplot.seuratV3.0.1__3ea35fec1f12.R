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
})
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

PRO<-readRDS("/data/NAS02/sh/Projects/FWSC20250326/00.data/result.rds")

# PRO[["RNA3"]] <- as(object = PRO[["RNA"]], Class = "Assay")
# DefaultAssay(PRO) <- 'RNA3'

PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)

# marker_file <- read.table("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/script/markergene.txt",header = FALSE,sep="\t",stringsAsFactors=FALSE)
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))

# marker <- unlist(strsplit(marker_file$V1, ","))
# marker <- unique(marker)
marker <- "Men1"
# print("gene not find: ")
# print(marker[!(marker %in% rownames(PRO))])
# marker<-intersect(marker,rownames(PRO))
outdir<-"/data/NAS02/sh/Projects/FWSC20250326/02.Mission/Men1/"

for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('lightgrey','red'),pt.size = pt_use, raster=FALSE)
    pdf(paste0(outdir,'Main-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,'Main-', i, '_Featureplot.png'))
    print(p)
    dev.off()

# Idents(PRO)<-PRO$sample2
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Main-', i, '_VlnPlot.pdf'),width = 10,height = 10)
		ggsave(paste0(outdir,'Main-', i, '_VlnPlot.png'),width = 10,height = 10)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Main",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Main'),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Main-', i, '_VlnPlot_sig_ana.pdf'),width = 10,height = 20)
		ggsave(paste0(outdir,'Main-', i, '_VlnPlot_sig_ana.png'),width = 10,height = 20)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')+geom_boxplot(width=.1,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir,'Main-group-', i, '_VlnPlot.pdf'),width = 28,height = 10)
		ggsave(paste0(outdir,'Main-group-', i, '_VlnPlot.png'),width = 28,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = T)
		get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix="Main_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Main-group-', i, '_VlnPlot_sig_ana.pdf'),width = 28,height = 15)
		ggsave(paste0(outdir,'Main-group-', i, '_VlnPlot_sig_ana.png'),width = 28,height = 15)

}

PRO<-readRDS("/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Macrophages/Macrophages.rds")
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('lightgrey','red'),pt.size = pt_use, raster=FALSE)
    pdf(paste0(outdir,'Macrophages-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,'Macrophages-', i, '_Featureplot.png'))
    print(p)
    dev.off()
		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)+geom_boxplot(width=.1,col="black",fill="white",outlier.shape = NA)
		ggsave(paste0(outdir,'Macrophages-', i, '_VlnPlot.pdf'),width = 8,height = 10)
		ggsave(paste0(outdir,'Macrophages-', i, '_VlnPlot.png'),width = 8,height = 10)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="Macrophages",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' Macrophages'),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Macrophages-', i, '_VlnPlot_sig_ana.pdf'),width = 8,height = 10)
		ggsave(paste0(outdir,'Macrophages-', i, '_VlnPlot_sig_ana.png'),width = 8,height = 10)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')+geom_boxplot(width=.1,col="black",fill="white",outlier.shape = NA)
		ggsave(paste0(outdir,'Macrophages-group-', i, '_VlnPlot.pdf'),width = 12,height = 10)
		ggsave(paste0(outdir,'Macrophages-group-', i, '_VlnPlot.png'),width = 12,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix="Macrophages_group",map_sig_level=T,only_sig = T)
		get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix="Macrophages_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'Macrophages-group-', i, '_VlnPlot_sig_ana.pdf'),width = 12,height = 15)
		ggsave(paste0(outdir,'Macrophages-group-', i, '_VlnPlot_sig_ana.png'),width = 12,height = 15)

}


PRO<-readRDS("/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Pulmonary alveolar cells/AT.rds")
PRO$group<-factor(PRO$group)
PRO$sample<-factor(PRO$sample)
cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
for (i in marker){
    p <- FeaturePlot(PRO, features=i, cols = c('lightgrey','red'),pt.size = pt_use, raster=FALSE)
    pdf(paste0(outdir,'AT-', i, '_Featureplot.pdf'))
    print(p)
    dev.off()
    png(paste0(outdir,'AT-', i, '_Featureplot.png'))
    print(p)
    dev.off()

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1)+geom_boxplot(width=.1,col="black",fill="white",outlier.shape = NA)
		ggsave(paste0(outdir,'AT-', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'AT-', i, '_VlnPlot.png'),width = 7,height = 7)
		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="AT",map_sig_level=T,only_sig = T,reset_y=T)
		get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix=paste0(i,' AT'),outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'AT-', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir,'AT-', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)

		p = VlnPlot(object = PRO, features = i,pt.size = 0.1, split.by = 'group')+geom_boxplot(width=.1,col="black",fill="white",outlier.shape = NA)
		ggsave(paste0(outdir,'AT-group-', i, '_VlnPlot.pdf'),width = 10,height = 10)
		ggsave(paste0(outdir,'AT-group-', i, '_VlnPlot.png'),width = 10,height = 10)

		p <- get_sig_ana(p,as.formula(paste0(i,'~split')),group.by='ident',prefix="AT_group",map_sig_level=T,only_sig = T)
		get_sig_ana(p$data,as.formula(paste0(i,'~split')),group.by='ident',prefix="AT_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir,'AT-group-', i, '_VlnPlot_sig_ana.pdf'),width = 10,height = 10)
		ggsave(paste0(outdir,'AT-group-', i, '_VlnPlot_sig_ana.png'),width = 10,height = 10)

}