rm(list = ls())
library(Seurat)
library(patchwork)
library(getopt)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
outdir = "/lvdata/wzb/scRNA/FW2023-504/2024.11.25/Epithelial_cells/Metastatic"
setwd(outdir)
forplot1 = read.csv("/lvdata/wzb/scRNA/FW2023-504/2024.11.25/Epithelial_cells/Metastatic/res.xls",sep="\t")

cols = c( "#5E4FA2" ,"#3288BD", "#66C2A5",
          "#ABDDA4" ,"#E6F598", "#FFFFBF",
          "#FEE08B", "#FDAE61" ,"#F46D43",
           "#D53E4F", "#9E0142"
          )
celltype = unique(forplot1$cell_type)
for (i in celltype){
  forplot = forplot1[forplot1$cell_type==i,]
  Group = unique(forplot$group)
  gg = ggplot(data=forplot,mapping=aes(x=group,y=Metastatic_Score))+
    geom_jitter(aes(colour  = group),
                width = 0.1,  size=2.5, stroke=0.5, shape=16,alpha=10) +
    geom_violin(alpha=0.5, linewidth=1, width=0.5,adjust = 0.9) +
    
    geom_boxplot(width=0.3, linewidth=1,outlier.color = NA, fill=NA) +
    scale_colour_manual(values=c("#7fc97f","#beaed4"))+
    stat_compare_means(method = "wilcox.test", comparisons = combn(Group, 2, simplify = FALSE), label = "p.format") +
    theme_bw()+theme(panel.grid=element_blank())
     # facet_wrap(~ celltype, nrow = 2)  # 按照 cell_type 分面
  gg
  ggsave(paste0(i,"_Metastatic_Score.pdf"),gg,width = 10,height = 8)
  ggsave(paste0(i,"_Metastatic_Score.png"),gg,width = 10,height = 8)
}