library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)

S <- readRDS"counts.RDS"
Idents(S) <- S$celltype

Idents(S) <- S$celltype
ALL.markers <- FindAllMarkers(S, only.pos = TRUE, min.pct = 0.25)
ALL.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10


S <- ScaleData(S, features = rownames(S))
DoHeatmap(S, features = top10$gene,angle = 45,label = T)+ 
  scale_fill_gradientn(colors = c('#330066','#336699','#66CC66','#FFCC33'))


DotPlot(S, features = unique(top10$gene)) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色

