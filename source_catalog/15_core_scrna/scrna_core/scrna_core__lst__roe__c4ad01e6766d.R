.libPaths(c(.libPaths(),"/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")) 
library(dplyr)
library(ggplot2)
library(ggsci)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(gtable)
source("/lvdata/wzb/pipline/Roe/source.R")
setwd('/lvdata/wzb/scRNA/FW2024-225_05/ROE')
# adata = readRDS("/lvdata/wzb/scRNA/FW2023-461/2023.09.07/Monocytic/h5_trans.rds")
# df = adata@meta.data
df = read.csv("/lvdata/wzb/scRNA/FW2024-225_05/All/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('cell_type','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "Roe.xls",sep = "\t",quote=F)
# roe = read.csv('/lvdata/wzb/scRNA/FW2024-213/fig/fig1/fig1c/Roe.xls',row.names = 1)
roe$celltype = rownames(roe)



# bubble_data = melt(roe,variable.name = "condition",value.name = "Reo")
# # 依据Roe值进行定义
# bubble_data = bubble_data %>% mutate(class=ifelse(Reo>=1.05,"Enrichment",ifelse(Reo>=0.95,"No Change","Depletion")))
# bubble_data$class = factor(bubble_data$class,levels = c("Depletion","No Change","Enrichment"))     #转换因子型排序
# col.order = read.csv('/lvdata/wzb/scRNA/FW2024-213/fig/color_cluster.xls',sep="\t")
# if(!is.null(col.order)){
#   bubble_data$celltype = factor(bubble_data$celltype,levels = col.order[,1])
# }
# row.order = c('NC','iMCD')
# if(!is.null(row.order)){
#   bubble_data$condition = factor(bubble_data$condition,levels = row.order)
# }
# # 创建预设的样式
# ang = 60
# bubblelowcolor = '#4a6fe3'
# bubblehighcolor = '#8e063b'
#   bubblemidcolor =  'grey'
# bubblesizelength = 5
# bubble_theme = theme_bw()+
#   # 页边距设置margin，panel.grid.major主网格线颜色
#   theme(plot.margin = margin(.5,.5,.5,.5,'cm'),panel.grid.major=element_line(colour=NA))+
#   # axis.text.x 横轴标签调整 
#   # hjust水平调整，vjust垂直调整，angle角度调整，设置为element_blank()为不显示
#   # size 调整字体大小
#   # axis.text.y 纵轴标签调整 同axis.text.x
#   theme(axis.text.x=element_text(hjust=1,vjust=1,angle = ang,size=10))
# # 绘制bubble_plot
# p1 = ggplot(bubble_data,aes(x=condition,y=celltype))+
#   # color显示颜色数值大小
#   geom_point(aes(color=class,size=Reo))+
#   # 颜色设置针对上面的color参数low为低值对应颜色，high为高值对应颜色
#   scale_color_manual(values = c(bubblelowcolor,bubblemidcolor,bubblehighcolor))+
#   # 调整圈大小
#   scale_size(range = c(2,8),
#              limits = c(floor(min(bubble_data$Reo)),ceiling(max(bubble_data$Reo))),
#              breaks = seq(floor(min(bubble_data$Reo)),floor(max(bubble_data$Reo)),length.out = bubblesizelength))+
#   # # 调整圈大小
#   # scale_size(range = c(2,8),
#   #            limits = c(floor(min(bubble_data$Reo)),ceiling(max(bubble_data$Reo))),
#   #            breaks = seq(floor(min(bubble_data$Reo)),ceiling(max(bubble_data$Reo)),length.out = bubblesizelength))+
#   # 横纵轴转换
#   coord_flip()+
#   labs(x = "Group",y = "Celltype",color="",size="Ro/e")
# 
# 
# p1+bubble_theme
# ggsave(paste0("Roe_","bubbleplot.png"),p1+bubble_theme,width = 24,height = 14,units = "cm")
# ggsave(paste0("Roe_","bubbleplot.pdf"),p1+bubble_theme,width = 24,height = 14,units = "cm")

p = pheatmap::pheatmap(roe, 
             color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
             cluster_rows = F, 
             cluster_cols = F, 
             display_numbers = TRUE,
             number_color = "black",
             cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
             border_color = '#ffffff',
             angle_col='45')
dev.off()
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
p
dev.off()
png("heatmap.png",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
p
dev.off()
