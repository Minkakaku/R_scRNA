library(Seurat)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(gt)

diff1<-diffgroup[[m]][1]
diff2<-diffgroup[[m]][2]
print(paste0(diff1,"vs",diff2))
data = read.table(paste0(outdir,cluster_diff_sub,"\\",diff1,"-",diff2, "_diffgenes.csv"),sep = ",",head=TRUE)
colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
data$log10PValue = -log10(data$PValue)

data$log10PValue[is.infinite(data$log10PValue)] = 300
cut_off_PValue =0.05 #设置FDR的阈值
cut_off_log2FC =0.25 #设置log2FC的阈值
data$Sig = ifelse(data$PValue < cut_off_PValue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(data$log2FC > cut_off_log2FC ,'Up','Down'),'no')
data = data.frame(data)
table(data$Sig)
# data = data %>%
#   filter(pct_nz_group >= 0.2) %>%
#   filter(pct_nz_reference <= 0.6) %>%
#   arrange(desc(abs(log2FC)),FDR)
xlim_min<-min(data$log2FC)-0.3
xlim_max<-max(data$log2FC)+0.3
ylim_max <- min(ceiling(max(data$log10PValue) / 50) * 50,300)
###绘图——基础火山图###
p1 <- ggplot(data, aes(x =log2FC, y=log10PValue, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=3) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_PValue), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10PValue") +  #x、y轴标签
  ggtitle("Volcano Plot") + #标题
  scale_x_continuous(limits = c(xlim_min, xlim_max)) +
  scale_y_continuous(limits = c(0, ylim_max), breaks = seq(0, ylim_max, 50)) +
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 15),
        legend.position="right", 
        legend.title = element_blank()
  ) 

Up_top_10 =(     #筛选差异显著上调的前10个Gene
  data %>%
    filter(Sig == 'Up') %>%
    filter(pct_nz_group >= 0.1) %>%
    arrange(desc(abs(log2FC)),PValue) %>%
    head(10)
)
Up_top_10 %>% gt() #数据制成表
Down_top_10 = (       #筛选差异显著下调的前10个Gene
  data %>%
    filter(Sig == 'Down') %>%
    filter(pct_nz_reference >= 0.1) %>%
    arrange( desc(abs(log2FC)),PValue) %>%
    head(10)
)
Down_top_10 %>% gt() #数据制成表

p2 <- p1 +
  geom_label_repel(data = Up_top_10,
                   aes(log2FC, log10PValue, label = gene_name),
                   size = 3, fill="#CCFFFF",
                   alpha = 0.65, color = "black")+
  geom_label_repel(data = Down_top_10,
                   aes(log2FC, log10PValue, label = gene_name),
                   size = 3, fill="#FFCCCC",
                   alpha = 0.65, color = "black")

p2 #出图
ggsave(paste0(cluster_diff_sub,"\\",diff1,"-",diff2, ".Volcanoplot.pdf"),width = 10,height = 10)
ggsave(paste0(cluster_diff_sub,"\\",diff1,"-",diff2, ".Volcanoplot.png"),width = 10,height = 10)