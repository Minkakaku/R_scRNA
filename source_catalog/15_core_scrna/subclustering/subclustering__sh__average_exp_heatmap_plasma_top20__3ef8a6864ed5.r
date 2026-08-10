# library(pheatmap)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.rds")

MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/script/pathway.csv",sep = "\t",head=FALSE)
genesets$V2 <- strsplit(genesets$V2, ",")
gene <- setNames(genesets$V2, genesets$V1)
#转换成list
gene<-c(as.list(MIF),gene)
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'
colnames(adata@meta.data)[18]<-'mTOR'
bdata<-subset(adata,group %in% c("P_after"))
# bdata$cell_type_group<-factor(paste0(bdata$raw_cell_type,"_",bdata$group),levels=c("Plasma cells_HD","Plasma cells_P","Plasma cells_P_after"))
# bdata$cell_type_group<-factor(paste0(bdata$raw_cell_type,"_",bdata$group),levels=rev(c("Plasma cells_HD","Plasma cells_P","Plasma cells_P_after")))
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")

outdir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/21.Mission/score_vlnplot/"
outdir_sub="/mnt/NexData_file/home/sh/Projects/FWAE20251216/21.Mission/score_vlnplot_gene/"
bdata_sub<-subset(bdata,cell_type %in% c("IgA Plasma","IgG Plasma C1","IgG Plasma C2","IgM Plasma"))

p <- VlnPlot(bdata_sub,features='rna_MIF',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Plasma-",'MIF',"-VlnPlot.pdf"),width = 9,height = 8)
ggsave(paste0(outdir_sub,"Plasma-",'MIF',"-VlnPlot.png"),width = 9,height = 8)

p <- get_sig_ana(p,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Plasma_group-",'MIF'),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Plasma_group-",'MIF'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Plasma_group_sig_ana-",'MIF',"-VlnPlot.pdf"),width = 9,height = 10)
ggsave(paste0(outdir_sub,"Plasma_group_sig_ana-",'MIF',"-VlnPlot.png"),width = 9,height = 10)

p <- VlnPlot(bdata_sub,features='mTOR',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Plasma-",'mTOR',"-VlnPlot.pdf"),width = 9,height = 8)
ggsave(paste0(outdir_sub,"Plasma-",'mTOR',"-VlnPlot.png"),width = 9,height = 8)

p <- get_sig_ana(p,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Plasma_group-",'mTOR'),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Plasma_group-",'mTOR'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Plasma_group_sig_ana-",'mTOR',"-VlnPlot.pdf"),width = 9,height = 10)
ggsave(paste0(outdir_sub,"Plasma_group_sig_ana-",'mTOR',"-VlnPlot.png"),width = 9,height = 10)

plot_data <- FetchData(bdata_sub, 
                       vars = c("cell_type", "MIF", "mTOR"))

# 重命名
colnames(plot_data) <- c("cell_type", "MIF", "mTOR")
# 转换成长格式
plot_long <- plot_data %>%
  pivot_longer(cols = c(MIF, mTOR),
               names_to = "geneset",
               values_to = "score")
# 绘制
p <- ggplot(plot_long, aes(x = cell_type, y = score, fill = geneset)) +
  geom_violin(trim = FALSE, position = position_dodge(0.7),width = 0.6, scale = "width") +
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(0.7), 
               size = 2, color = "black") +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(title = "Gene Set Expression by Cell Type",
       x = "Cell Type", y = "Expression Level") +
  scale_fill_manual(values = c("MIF" = "#E41A1C", 
                               "mTOR" = "#377EB8"),
                    labels = c("MIF", "mTOR"))
ggsave(paste0(outdir,"Plasma-VlnPlot.pdf"),width = 12,height = 10)
ggsave(paste0(outdir,"Plasma-VlnPlot.png"),width = 12,height = 10)

p1 <- p +  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     size = 7,                    # 字体大小
                     bracket.size = 1.2,          # 连线粗细
                     tip.length = 0.03,           # 连线端点长度
                     vjust = 0.5,
                     label.y = max(plot_long$score)+0.3)
ggsave(plot=p1, paste0(outdir,"Plasma_sig_ana-VlnPlot.pdf"),width = 12,height = 10)
ggsave(plot=p1, paste0(outdir,"Plasma_sig_ana-VlnPlot.png"),width = 12,height = 10)


bdata_sub<-subset(bdata,cell_type %in% c("Activated B","Memory B C1","Memory B C2","Naive B"))

p <- VlnPlot(bdata_sub,features='rna_MIF',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Bcells-",'MIF',"-VlnPlot.pdf"),width = 9,height = 8)
ggsave(paste0(outdir_sub,"Bcells-",'MIF',"-VlnPlot.png"),width = 9,height = 8)

p <- get_sig_ana(p,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Bcells_group-",'MIF'),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Bcells_group-",'MIF'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Bcells_group_sig_ana-",'MIF',"-VlnPlot.pdf"),width = 9,height = 10)
ggsave(paste0(outdir_sub,"Bcells_group_sig_ana-",'MIF',"-VlnPlot.png"),width = 9,height = 10)

p <- VlnPlot(bdata_sub,features='mTOR',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Bcells-",'mTOR',"-VlnPlot.pdf"),width = 9,height = 8)
ggsave(paste0(outdir_sub,"Bcells-",'mTOR',"-VlnPlot.png"),width = 9,height = 8)

p <- get_sig_ana(p,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Bcells_group-",'mTOR'),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Bcells_group-",'mTOR'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Bcells_group_sig_ana-",'mTOR',"-VlnPlot.pdf"),width = 9,height = 10)
ggsave(paste0(outdir_sub,"Bcells_group_sig_ana-",'mTOR',"-VlnPlot.png"),width = 9,height = 10)

plot_data <- FetchData(bdata_sub, 
                       vars = c("cell_type", "MIF", "mTOR"))
colnames(plot_data) <- c("cell_type", "MIF", "mTOR")
# 转换成长格式
plot_long <- plot_data %>%
  pivot_longer(cols = c(MIF, mTOR),
               names_to = "geneset",
               values_to = "score")
# 绘制
p <- ggplot(plot_long, aes(x = cell_type, y = score, fill = geneset)) +
  geom_violin(trim = FALSE, position = position_dodge(0.7),width = 0.6, scale = "width") +
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(0.7), 
               size = 2, color = "black") +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(title = "Gene Set Expression by Cell Type",
       x = "Cell Type", y = "Expression Level") +
  scale_fill_manual(values = c("MIF" = "#E41A1C", 
                               "mTOR" = "#377EB8"),
                    labels = c("MIF", "mTOR"))
ggsave(paste0(outdir,"Bcells-VlnPlot.pdf"),width = 12,height = 10)
ggsave(paste0(outdir,"Bcells-VlnPlot.png"),width = 12,height = 10)

p1 <- p +  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     size = 7,                    # 字体大小
                     bracket.size = 1.2,          # 连线粗细
                     tip.length = 0.03,           # 连线端点长度
                     vjust = 0.5,
                     label.y = max(plot_long$score)+0.3)
ggsave(plot=p1, paste0(outdir,"Bcells_sig_ana-VlnPlot.pdf"),width = 12,height = 10)
ggsave(plot=p1, paste0(outdir,"Bcells_sig_ana-VlnPlot.png"),width = 12,height = 10)


bdata_sub<-subset(bdata,cell_type %in% c("CCL2 Classical monocyte","cDCs","DUSP1 Classical monocyte","LRMDA Classical monocyte","Non-classical monocyte", "pDCS", "VCAN Classical monocyte"))

p <- VlnPlot(bdata_sub,features='rna_MIF',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Monocyte-",'MIF',"-VlnPlot.pdf"),width = 12,height = 8)
ggsave(paste0(outdir_sub,"Monocyte-",'MIF',"-VlnPlot.png"),width = 12,height = 8)

p <- get_sig_ana(p,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Monocyte_group-",'MIF'),map_sig_level=T,only_sig = T,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('rna_MIF','~ident')),group.by='ident',prefix=paste0("Monocyte_group-",'MIF'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Monocyte_group_sig_ana-",'MIF',"-VlnPlot.pdf"),width = 12,height = 10)
ggsave(paste0(outdir_sub,"Monocyte_group_sig_ana-",'MIF',"-VlnPlot.png"),width = 12,height = 10)

p <- VlnPlot(bdata_sub,features='mTOR',group.by = 'cell_type',pt.size = 0)
ggsave(paste0(outdir_sub,"Monocyte-",'mTOR',"-VlnPlot.pdf"),width = 12,height = 8)
ggsave(paste0(outdir_sub,"Monocyte-",'mTOR',"-VlnPlot.png"),width = 12,height = 8)

p <- get_sig_ana(p,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Monocyte_group-",'mTOR'),map_sig_level=T,only_sig = T,method = "wilcox.test",reset_y=T)
get_sig_ana(p$data,as.formula(paste0('mTOR','~ident')),group.by='ident',prefix=paste0("Monocyte_group-",'mTOR'),outdir=outdir_sub,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir_sub,"Monocyte_group_sig_ana-",'mTOR',"-VlnPlot.pdf"),width = 12,height = 10)
ggsave(paste0(outdir_sub,"Monocyte_group_sig_ana-",'mTOR',"-VlnPlot.png"),width = 12,height = 10)

plot_data <- FetchData(bdata_sub, 
                       vars = c("cell_type", "MIF", "mTOR"))
colnames(plot_data) <- c("cell_type", "MIF", "mTOR")
# 转换成长格式
plot_long <- plot_data %>%
  pivot_longer(cols = c(MIF, mTOR),
               names_to = "geneset",
               values_to = "score")
# 绘制
p <- ggplot(plot_long, aes(x = cell_type, y = score, fill = geneset)) +
  geom_violin(trim = FALSE, position = position_dodge(0.7),width = 0.6, scale = "width") +
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(0.7), 
               size = 2, color = "black") +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(title = "Gene Set Expression by Cell Type",
       x = "Cell Type", y = "Expression Level") +
  scale_fill_manual(values = c("MIF" = "#E41A1C", 
                               "mTOR" = "#377EB8"),
                    labels = c("MIF", "mTOR"))
ggsave(paste0(outdir,"Monocyte-VlnPlot.pdf"),width = 30,height = 10)
ggsave(paste0(outdir,"Monocyte-VlnPlot.png"),width = 30,height = 10)

p1 <- p +  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     size = 7,                    # 字体大小
                     bracket.size = 1.2,          # 连线粗细
                     tip.length = 0.03,           # 连线端点长度
                     vjust = 0.5,
                     label.y = max(plot_long$score)+0.3)
ggsave(plot=p1, paste0(outdir,"Monocyte_sig_ana-VlnPlot.pdf"),width = 30,height = 10)
ggsave(plot=p1, paste0(outdir,"Monocyte_sig_ana-VlnPlot.png"),width = 30,height = 10)