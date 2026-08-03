setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require6/')
rm(list = ls())
gc()
library(ggplot2)
library(dplyr)
outdir <- "path/to/your/directory"

# 获取所有 CSV 文件的绝对路径
# csv_files <- list.files(
#   path = '/mnt/NexData_file/home/zwl/Project/FWSC20241172/require4/group_DEG_enrich',
#   pattern = "\\.csv$",        # 匹配 .csv 后缀
#   full.names = TRUE,          # 返回完整路径
#   recursive = FALSE           # 是否递归搜索子目录
# )
for (ci in c('B cells','T cells','Endothelium','Myeloid','SMC')) {
  DEG = read.csv(file = paste0('../require4/group_DEG_enrich/1.',ci,'_EXP_DEG_plot.csv'),row.names = 1)
  DEG$threshold = factor(ifelse(DEG$pvals_adj>0.05|abs(DEG$logfoldchanges) <1,'Nosign',
                                ifelse(DEG$logfoldchanges>1,'Up','Down')),
                         levels = c('Up','Nosign','Down'))
  top_genes = DEG %>% dplyr::filter(threshold != 'Nosign') %>% 
    arrange(logfoldchanges) %>%
    slice(1:10, (n() - 9):n())
  DEG$y_value <- -log10(DEG$pvals_adj)
  max_y <- max(DEG$y_value[is.finite(DEG$y_value)]) * 1.1  # 设置 y 轴上限为最大值的 1.1 倍
  
  gg = ggplot(DEG, aes(x = logfoldchanges, -log10(pvals_adj), colour = threshold)) +
    geom_jitter(size=0.8,stroke=0.5, alpha = 0.6) +
    scale_color_manual(values = c('Up'="firebrick", 'Down'="navy",'Nosign'='grey'))+
    theme_classic()+labs(x='logFC EXP/CON',y='-log10 pvalue_adj')+
    coord_cartesian(ylim = c(0, max_y)) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = names), size = 3, max.overlaps = Inf,box.padding = 0.5 ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")
  
  ggsave(filename = paste0('./1.',ci,'_group_DEG_Valcno.png'),gg,w=7,h=5,dpi = 200)
  ggsave(filename = paste0('./1.',ci,'_group_DEG_Valcno.pdf'),gg,w=7,h=5,dpi = 200)
}

#################################
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require7')
library(sceasy)
library(Seurat)
library(reticulate)
library(scales)

sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241172/adata_sample5.h5ad', from="anndata", to="seurat",
                      outFile='adata_sample5.rds')
obj <-readRDS('./adata_sample5.rds')
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','group'))

#################################ALL group DEG########

allDEG = list()
for (i in cellfile){
  j = gsub(' ','_',gsub('_EXP_DEG_plot.csv','',gsub('1.','',strsplit(i,split = '/')[[1]][10])))
  data = read.csv(file = i,row.names = 1) %>% dplyr::filter(pvals_adj <0.05,abs(logfoldchanges)>1)
  allDEG[[j]] = data$names
}
library(UpSetR)
library(ggplot2)
pdf(file= '../require9/ALL_cluster_DEG.pdf') # or other device
UpSetR::upset(fromList(allDEG),nsets = 15)
dev.off()
png(file= '../require9/ALL_cluster_DEG.png',res = 200,width = 1500,height = 1500) # or other device
UpSetR::upset(fromList(allDEG),nsets = 15)
dev.off()
