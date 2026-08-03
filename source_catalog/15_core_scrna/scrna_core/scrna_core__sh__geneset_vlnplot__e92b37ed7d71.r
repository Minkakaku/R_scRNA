library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(corrplot)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/AE2408071595-01/02.Mission/Monocytes/Monocyte.rds")
# adata$group<-factor(adata$group)
# adata$sample<-factor(adata$sample)

# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = ",",head=TRUE)
geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = ",",head=FALSE)
#转换成list
gene<-geneset$V1
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2408071595-01/02.Mission/Monocytes-group_cytokine_vlnplot/"
gene <- intersect(gene, rownames(adata))
for (g in gene){
# Idents(adata)<-adata$annotation
# p <- VlnPlot(adata,features=g,pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features=g,split.by = 'group',pt.size = 0.1)+scale_fill_manual(values = c(NC='#1f77b4',HLH='#ff7f0e','iMCD-TAFRO'='#279e68'))
ggsave(paste0(outdir,g,"_group.VlnPlot.pdf"),width = 16,height = 10)
ggsave(paste0(outdir,g,"_group.VlnPlot.png"),width = 16,height = 10)

p <- get_sig_ana(p,as.formula(paste0(g,' ~ split')),group.by='ident',prefix=paste0(g,'_group'),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(NC='#1f77b4',HLH='#ff7f0e','iMCD-TAFRO'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(g,' ~ split')),group.by='ident',prefix=paste0(g,'_group'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,g,"_group_sig_ana.VlnPlot.pdf"),width = 16,height = 10)
ggsave(paste0(outdir,g,"_group_sig_ana.VlnPlot.png"),width = 16,height = 10)
}