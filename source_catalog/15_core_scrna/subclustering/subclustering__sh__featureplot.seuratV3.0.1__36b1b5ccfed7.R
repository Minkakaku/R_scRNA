#!/usr/bin/env Rscript
suppressMessages({
library(ggplot2)
library(reshape2)
library(argparser)
library(Cairo)
library(Seurat)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
library(dplyr)
library(ggpubr)
library(ggbreak)
library(scCustomize)
})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input rds")
argv <- add_argument(argv,"--prefix", help="group name prefix")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--markergene", help="markergene name file")
#argv <- add_argument(argv,"--clustcol", help="color source", default='rdsmeta')

argv <- add_argument(argv,"--cluster", help="subcluster list ,split by ,", default='all')
argv <- add_argument(argv,"--sample", help="sample list ,split by ,", default='all')
argv <- add_argument(argv,'--group', help="group file, sample group. add group info to rds and plot")
argv <- add_argument(argv,'--ident', help='cluster or sample or group or colnames of meta data', default='cluster')
argv <- add_argument(argv,'--split', help='sample or group or colnames of meta data', default='F')
argv <- add_argument(argv,'--sampleorder', help="sampleorder for plot")
argv <- add_argument(argv,'--grouporder', help="grouporder for plot")
argv <- add_argument(argv, '--clusterorder', help='default origin cluster order')

argv <- add_argument(argv,'--groupby', help='sample or group or colnames of meta data', default='NULL')
argv <- add_argument(argv,'--splitby', help='sample or group or colnames of meta data', default='NULL')
argv <- add_argument(argv,"--slot", help="slot:data|counts|scale.data",default="data")
argv <- add_argument(argv,"--plotmethod", help="{tsne|umap}",default="umap")
argv <- add_argument(argv,"--mincutoff", help="Vector of minimum cutoff values for each feature",default='q1')
argv <- add_argument(argv,"--maxcutoff", help="Vector of maximum cutoff values for each feature",default='q99')
argv <- add_argument(argv,"--vlnpoint", help="Point size for geom_violin",default=0)

argv <- add_argument(argv,"--feature", help="only gene list featureplot", default='T')
argv <- add_argument(argv,"--dot", help="only gene list dotplot", default='F')
argv <- add_argument(argv,"--vln", help="only gene list vlnplot", default='F')
argv <- add_argument(argv,"--box", help="only gene list boxplot", default='F')
argv <- add_argument(argv,"--bar", help="only gene list barplot", default='F')
argv <- parse_args(argv)

argv$outprefix <- file.path(argv$outdir, argv$prefix)

source("F:/新格元/plotmarker/seuratSingleron.R")
source("F:/新格元/plotmarker/seurat_function.R")
source("F:/新格元/plotmarker/plot_Vln.R")
source("F:/新格元/plotmarker/boxplot_function.r")
source("F:/新格元/plotmarker/heatmap_function.R")
source("F:/新格元/plotmarker/Barplot_functionV3.R")
#color<-argv$clustcol
#if(color!='rdsmeta'){
    source('/Public/Script/shouhou/SCRIPT/lib/scRNA/config/colorConfig.R')
    clustcol <- colorNew
#}

outdir<-argv$outdir
prefix<-argv$prefix
createDir(argv$outdir)
PRO <- readRDS(argv$rds)
PRO <- modifyPRO(PRO, argv)
PRO$sample2<-factor(PRO$sample2,levels=c('HLH01-LAHS-B','HLH02-LAHS-T','HLH01-LAHS-B-Post','HLH03-LAHS-T','HLH03-LAHS-T-Post','HLH04-LAHS','HLH05-NLAHS-FHLH','HLH06-NLAHS-EBV-BM','HLH06-NLAHS-EBV-PB','HLH07-NLAHS','HLH07-NLAHS-Post','HC01','HC02'))
# pal <- viridis(n = 10, option = "C", direction = -1)
DefaultAssay(PRO) <-'RNA'
groupby<-argv$groupby
splitby<-argv$splitby
marker_file <- read.table(argv$markergene,header = TRUE,sep="\t",row.names = 1,stringsAsFactors=FALSE)
cell_name <- rownames(marker_file)

cell_number<-length(colnames(PRO))
pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
cluster_number <- length(levels(PRO))
pdfwidth <- ifelse(cluster_number<10, 6, (cluster_number-10)*0.15+6)
# pngwidth <- ifelse(cluster_number<10, 480, (cluster_number-10)*10+480)
maxnchar <- max(nchar(levels(PRO)))
pdfheight <- ifelse(maxnchar<6, 3, maxnchar/15+3)
# pngheight <- ifelse(maxnchar<6, 240, maxnchar*3+240)

allMarker <- c()
for (cell in cell_name){
    marker <- unlist(strsplit(marker_file[cell,1], ","))
    marker <- unique(marker)

    print("gene not find: ")
    print(marker[!(marker %in% rownames(PRO))])
    marker<-intersect(marker,rownames(PRO))
    allMarker <- c(allMarker, marker)
    for (i in marker){
        if(argv$feature!='F'){
            # filename<-paste(argv$outprefix, cell, i, 'Featureplot.pdf', sep='.')
						# filename<-paste(argv$outprefix, i, 'Featureplot.pdf', sep='_')
            p <- FeaturePlot_scCustom(PRO, features=i, colors_use = c("lightgrey","red"),split.by = "sample2",min.cutoff='q1',max.cutoff='q99',pt.size = pt_use,reduction=argv$plotmethod,slot = argv$slot,num_columns = 4)
            pdf(paste0(argv$outprefix,'/', i, '_Featureplot.pdf'),width = 26,height=26)
            print(p)
            dev.off()
            png(paste0(argv$outprefix,'/', i, '_Featureplot.png'),width = 26/7*480,height=26/7*480)
            print(p)
            dev.off()
        }
    if(argv$vln!='F'){
        # plot_Vln(PRO,features=i,outdir = outdir,prefix = prefix,groupby = "cluster",splitby =  "group")
        }
				argv$outprefix = "/data/NAS02/sh/Projects/FWSC2024-192/03.Mission/NLAHS_LAHS-NC/VlnPlot/"
				p = VlnPlot(object = PRO, features = i,pt.size = 0)+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)+ scale_fill_manual(values = c("#be348f","#be348f","#be348f","#be348f","#be348f","#be348f","#197ec0", "#197EC0", "#197ec0","#197ec0","#197ec0","#999999","#999999"))
				pdf(paste0(argv$outprefix,'/', i, '_VlnPlot.pdf'),width = 13)
        print(p)
        dev.off()
        png(paste0(argv$outprefix,'/', i, '_VlnPlot.png'),width = 13/7*480)
        print(p)
        dev.off()
    }
}

if(argv$dot!='F'){
    dotplotTop(PRO, rev(unique(allMarker)), paste0("F:/新格元/plotmarker/", '.genelist_Dotplot.pdf'), exp=T)
}

if(argv$box!='F'|argv$bar!='F'){
    avg <- t(as.data.frame(AverageExpression(PRO,features = marker,add.ident="group")$RNA))
    cg<-t(data.frame(matrix(unlist(strsplit(rownames(avg),split = "_")),2,length(rownames(avg)))))
    colnames(cg)<-c("cluster","group")
    avgExp<-as.data.frame(cbind(cg,avg))
    write.table(avgExp, file=paste0(argv$outprefix,"/","average.xls"),sep='\t',quote=F,row.names=F,col.names=T)
}

if(argv$box!='F'){
    num<-length(marker)+2
    avgExp[,c(3:num)]<-as.numeric(unlist(avgExp[,c(3:num)]))
for (i in marker){
    boxplot_function(data=avgExp,x="group",y=i,group="group",filltype="all", fig_save_name=paste("F:/新格元/plotmarker/boxplot_",i))
    }
}

if(argv$bar!='F'){
    
    for (i in marker){
    GroupBarPlot_function(data=avgExp, x="cluster", y=i, group="group", filltype="fill", errorbar=T, usecolor=PRO@meta.data$group_colors, fig_save_path="F:/新格元/plotmarker/",'barplot_',i)
    }
}
