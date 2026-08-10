# import scanpy as sc
# import os
# import pandas as pd
# from scanpy_fun_subcluster import find_group_marker
# adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/Myeloid.h5ad")
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Non-classical Monocytes']), :]
# output="/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/volcano/"
# sc.tl.rank_genes_groups(adata_sub,groupby='group',method='wilcoxon',pts=True)
# find_group_marker(adata=adata_sub,species='human',output=output, logfc=0.25, fdr = 0.05, pct = 0.1)
# all_diff = []
# groups = adata_sub.uns['rank_genes_groups']['names'].dtype.names
# for c in groups:
# 		tmp = sc.get.rank_genes_groups_df(adata_sub, group=c)
# 		tmp['cluster'] = c
# 		all_diff.append(tmp)
# all_diff = pd.concat(all_diff, axis=0)
# all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
# all_diff.to_csv(os.path.join(output, 'all_cluster_marker.xls'), sep='\t', index=None)


rm(list = ls())
# BiocManager::install("ggplot2")
library(ggplot2)
library(magrittr)
library(dplyr)
input = "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/volcano/"
output = "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/volcano/"
setwd(input)
set.seed(42)
data = read.csv("./all_cluster_marker.xls",sep = "\t")
###################################勿动#########################################
# 是否去除表达率较低的数据
data = data[data$pct_nz_group >0.1,]
################################################################################




# 依据数值对应颜色进行区分文字
data$threshold = factor(ifelse(data$pvals_adj < 0.05 & abs(data$logfoldchanges) >= 0.25,#p值小于0.05 log2FC绝对值大于0.25进行另一步判断
                               ifelse(data$logfoldchanges >= 0.25 ,'Up','Down'),    #p值大于0.25设置为上调，小于0.25定义为下调
                               'NoSignificant'),                                #其中的定义为不显著
                        levels=c('Up','Down','NoSignificant'))                  #定义向量排列顺序

# BiocManager::install("RImagePalette")
# 可以使用这个包查看颜色
# library(RImagePalette)
# scales::show_col(c("#DC143C","#00008B","#808080"))
# 绘图
###################################勿动#########################################
# 调整横轴顺序
data$cluster = factor(data$cluster,levels = c("control","case"))
################################################################################
jitter_position <- position_jitter(width = 0.2, height = 0, seed = 123)
set.seed(123)  
p=ggplot(data,aes(x=cluster,y=logfoldchanges,color=threshold))+                        #横轴为细胞类型，纵轴为log2FC数值，颜色以分类变量
  geom_point(aes(x=cluster
                 # ,alpha=-log10(pvals)
                 ),                                          #透明度使用-log10(p_value)
             # position="jitter",                                                 #添加微小扰动点
             position = jitter_position,
             pch=16, 
             #点的形状
             cex=1)+                                                            #点大小
  scale_color_manual(values=c("#DC143C",                                        #设置theshold第一种level up红色
                              "#00008B",                                        #设置theshold第二种level down蓝色
                              "#808080"))+                                      #设置theshold第二种level NoSignificant灰色
  theme_bw()+
  geom_hline(yintercept = c(-0.25,0.25),                                        #添加线的位置
             lty=3,                                                             #添加线的类型
             col="black",                                                       #添加线的颜色
             lwd=0.5)+                                                          #添加线的粗细
  theme(axis.text.x=element_text(hjust=1,                                       #水平对齐方式0.5居中，较小值向左，大值向右
                                 # vjust=0.5,                                       #垂直对齐方式0.5居中，较小值向上，大值向下
                                 angle=45,                                      #逆时针旋转角度
                                 size=10))                                      #调整横坐标轴标签位置以及角度大小
# 可以预览
p
# 输出
# width宽度 height高度
ggsave(paste0(output,"plot.png"),p,width=length(levels(data$type))*2.5,height=8)
ggsave(paste0(output,"plot.pdf"),p,width=length(levels(data$type))*2.5,height=8)

library(ggrepel)
Up_top_10 <- data %>%
  filter(threshold == 'Up') %>%
  arrange(desc(abs(logfoldchanges)), pvals_adj) %>%
  group_by(cluster) %>%
  slice_head(n = 10)

Down_top_10 <- data %>%
  filter(threshold == 'Down') %>%
  arrange(desc(abs(logfoldchanges)), pvals_adj) %>%
  group_by(cluster) %>%
  slice_head(n = 10)


p1 <- p + 
  geom_label_repel(data = Up_top_10,
                   aes(x = cluster, y = logfoldchanges, label = Symbol),
                   
                   position = jitter_position,  # 控制标签的扰动
                   size = 3, fill = "#FFCCCC", alpha = 0.65, color = "black", 
                   max.overlaps = 50,
                   segment.color = "transparent", # 去掉连线
                   segment.size = 0)+
  geom_label_repel(data = Down_top_10,
                   aes(x = cluster, y = logfoldchanges, label = Symbol),

                   position = jitter_position,  # 控制标签的扰动
                   size = 3, fill = "#CCFFFF", alpha = 0.65, color = "black",
                   max.overlaps = 50,
                   segment.color = "transparent", # 去掉连线
                   segment.size = 0)+
  labs(x="")
  
p1
ggsave(paste0(output,"plot_top10.png"),p1,width=length(levels(data$cluster))*4,height=8)
ggsave(paste0(output,"plot_top10.pdf"),p1,width=length(levels(data$cluster))*4,height=8)






