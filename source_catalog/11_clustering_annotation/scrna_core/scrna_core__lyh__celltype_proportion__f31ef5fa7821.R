library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsci)
library(cowplot)
args=commandArgs(T)

data <- read.table(args[1], sep = '\t', header = T, row.names = 1)
xaxis = args[2]
out_img = paste(args[3], 'celltype_proportion.png', sep = '/')

cellratio <- prop.table(table(data$cell_type, data[, xaxis]), margin = 2) %>%
  as.data.frame() %>%
  `colnames<-`(c('CellType', 'sample', 'proportion')) %>%
  filter(proportion!=0)

n = length(unique(cellratio$CellType))
if(n<10){
  mycolor = pal_npg()(10)
} else {
  mycolor = colorRampPalette(colors = pal_npg()(10))(n)
}

p <- ggplot(cellratio, aes(x=sample, y=proportion, fill=CellType))+
  geom_bar(stat = "identity",width = 0.95,colour = '#222222')+
  theme_classic()+
  labs(x=args[2],y = 'proportion',fill='CellType')+
  # coord_flip()+
  scale_fill_manual(values=rev(mycolor))+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=45, hjust=1))
ggsave(p, filename = out_img, width = as.numeric(args[4]), height = as.numeric(args[5]))
