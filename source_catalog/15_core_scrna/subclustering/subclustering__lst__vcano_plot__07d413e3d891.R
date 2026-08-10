rm(list = ls())
library(ggplot2)

library(ggrepel)
setwd('/lvdata/wzb/scRNA/FW2024-225_03/4.Recluster/NK_cells/celltype/DEG_GZMK')
data = read.csv("/lvdata/wzb/scRNA/FW2024-225_03/4.Recluster/NK_cells/celltype/DEG_GZMK/clusteriMCD/differentially_expressed_genes.xls",sep = "\t")
colnames(data)
colnames(data) <- c("gene_name", "scores", "log2FC", "PValue", "FDR", "pct_nz_group" ,    "pct_nz_reference")
data$log10FDR = -log10(data$FDR)

data$log10FDR[is.infinite(data$log10FDR)] = 300
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =1 #设置log2FC的阈值
data$Sig = ifelse(data$FDR < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(data$log2FC > cut_off_log2FC ,'Up','Down'),'no')
data = data.frame(data)
table(data$Sig)
# data = data %>%
#   filter(pct_nz_group >= 0.2) %>%
#   filter(pct_nz_reference <= 0.6) %>%
#   arrange(desc(abs(log2FC)),FDR)

###绘图——基础火山图###
p1 <- ggplot(data, aes(x =log2FC, y=log10FDR, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(-7,7) +  #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("Volcano Plot") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()
  ) 
p1
# install.packages("gt")
library(dplyr) # 用于数据处理
library(gt) # 制作表格
Up_top_10 =(     #筛选差异显著上调的前10个Gene
  data %>%
    filter(Sig == 'Up') %>%
    filter(pct_nz_group >= 0.2) %>%
    filter(pct_nz_reference <= 0.1) %>%
    arrange(desc(abs(log2FC)),FDR) %>%
    head(10)
)
Up_top_10 %>% gt() #数据制成表
Down_top_10 = (       #筛选差异显著下调的前10个Gene
  data %>%
    filter(Sig == 'Down') %>%
    filter(pct_nz_group >= 0.2) %>%
    arrange( desc(abs(log2FC)),FDR) %>%
    head(10)
)
Down_top_10 %>% gt() #数据制成表



p2 <- p1 +
  geom_label_repel(data = Up_top_10,
                            aes(log2FC, log10FDR, label = gene_name),
                            size = 3, fill="#CCFFFF",
                            alpha = 0.65, color = "black")+
  geom_label_repel(data = Down_top_10,
                   aes(log2FC, log10FDR, label = gene_name),
                   size = 3, fill="#FFCCCC",
                   alpha = 0.65, color = "black")


p2 #出图
ggsave("iMCDvsNC.Volcanoplot.pdf",width = 10,height = 10)
ggsave("iMCDvsNC.Volcanoplot.png",width = 10,height = 10)

