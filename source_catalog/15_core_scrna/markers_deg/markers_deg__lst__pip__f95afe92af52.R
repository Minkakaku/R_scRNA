rm(list = ls())
# BiocManager::install("ggplot2")
library(ggplot2)
input = "/lvdata/wzb/pipline/other/DEG_group_plot/Use/input"
output = "/lvdata/wzb/pipline/other/DEG_group_plot/Use/output"
setwd(input)
set.seed(42)
###################################勿动#########################################
filelist = list.files(pattern = ".xls$")                                        #input中包含txt格式文件
for(i in 1:length(filelist)){
  file <- gsub('.xls','',filelist[i])                                           #删除txt后缀
  soucefile = paste0(file,".xls")                                    
  input<-read.table(soucefile,header = T,check.names = F,sep='\t')
  input$type=file                                                               #添加文件名称到type列
  if(i==1){
    data=input
    }
  else{
    data=rbind(data,input)
    }
}
################################################################################

###################################勿动#########################################
# 是否去除表达率较低的数据
data = data[data$pct_nz_group >0.1,]
################################################################################




# 依据数值对应颜色进行区分文字
data$threshold = factor(ifelse(data$pvals < 0.05 & abs(data$logfoldchanges) >= 0.25,#p值小于0.05 log2FC绝对值大于0.25进行另一步判断
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
data$type = factor(data$type,levels = c("NC","iMCD-TAFRO","iMCD-TAFRO_after"))
################################################################################
jitter_position <- position_jitter(width = 0.2, height = 0, seed = 123)
set.seed(123)  
p=ggplot(data,aes(x=type,y=logfoldchanges,color=threshold))+                        #横轴为细胞类型，纵轴为log2FC数值，颜色以分类变量
  geom_point(aes(x=type
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
                                 vjust=0.5,                                       #垂直对齐方式0.5居中，较小值向上，大值向下
                                 angle=90,                                      #逆时针旋转角度
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
  arrange(desc(abs(logfoldchanges)), pvals) %>%
  group_by(type) %>%
  slice_head(n = 10)

# Down_top_10 <- data %>%
#   filter(threshold == 'Down') %>%
#   arrange(desc(abs(logfoldchanges)), pvals) %>%
#   group_by(type) %>%
#   slice_head(n = 10)


p1 <- p + 
  geom_label_repel(data = Up_top_10,
                   aes(x = type, y = logfoldchanges, label = names),
                   
                   position = jitter_position,  # 控制标签的扰动
                   size = 3, fill = "#FFCCCC", alpha = 0.65, color = "black", 
                   max.overlaps = 50,
                   segment.color = "transparent", # 去掉连线
                   segment.size = 0)
# +
#   geom_label_repel(data = Down_top_10,
#                    aes(x = type, y = logfoldchanges, label = names),
#                    
#                    position = jitter_position,  # 控制标签的扰动
#                    size = 3, fill = "#CCFFFF", alpha = 0.65, color = "black", 
#                    max.overlaps = 50,
#                    segment.color = "transparent", # 去掉连线
#                    segment.size = 0)
  
p1
ggsave(paste0(output,"plot_top10.png"),p1,width=length(levels(data$type))*3,height=8)
ggsave(paste0(output,"plot_top10.pdf"),p1,width=length(levels(data$type))*3,height=8)






