library(monocle)
setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou17/")
rm(list = ls())

################################################################################
### 需求1 monocle分支明确一下
load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

cds@reducedDimK %>% 
  t %>% 
  as.data.frame() %>%
  rownames_to_column("label") %>%
  
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(., label %in% c("Y_114", "Y_98")),
    aes(label = label),
    min.segment.length = 0,
    color = "red"
  )+
  theme_classic()

ggsave(plot = last_plot(),filename = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou8/label.png")

gc()


library(tidyverse)
library(ggrepel)

df <- cds@reducedDimK %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("label")

ggplot(df, aes(V1, V2)) +
  geom_point(
    colour = "grey60",
    size = 2,
    alpha = 0.45
  ) +
  geom_point(
    data = subset(df, label %in% c("Y_114","Y_98")),
    colour = "#D73027",
    size = 4
  ) +
  geom_label_repel(
    data = subset(df, label %in% c("Y_114","Y_98")),
    aes(label = label),
    fill = "white",
    colour = "#D73027",
    fontface = "bold",
    size = 4.5,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "grey40",
    label.size = 0.25
  ) +
  coord_equal() +
  # labs(
  #   x = "Component 1",
  #   y = "Component 2"
  # ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.border = element_blank()
  )


g <- minSpanningTree(cds)

plot(g)
igraph::degree(g)
cds@auxOrderingData$DDRTree$branch_points

df <- cds@reducedDimK %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("label")

neighbors(g, "Y_113")

g2 <- delete_vertices(g, "Y_113")

comp <- components(g2)

branch114_nodes <- names(comp$membership)[
  comp$membership ==
    comp$membership[which(names(comp$membership)=="Y_114")]
]

branch98_nodes <- names(comp$membership)[
  comp$membership ==
    comp$membership[which(names(comp$membership)=="Y_98")]
]

df$Branch <- "Others"

df$Branch[df$label %in% branch114_nodes] <- "Branch after Y_114"

df$Branch[df$label %in% branch98_nodes] <- "Branch after Y_98"

## 把真正branch point也单独标出来
df$Branch[df$label=="Y_113"] <- "Branch point"

ggplot(df,
       aes(V1,V2,color=Branch))+
  
  geom_point(
    size=1.5,
    alpha=0.9
  )+
  
  geom_label_repel(
    data=subset(df,label %in% c("Y_113","Y_114","Y_98")),
    aes(label=label),
    fill="white",
    fontface="bold",
    size=8 / .pt,
    label.size=.25,
    
    ## 标签与点的距离
    box.padding = 0.5,
    point.padding = 0.4,
    
    ## 引导线
    segment.color = "black",
    segment.size = 0.5
  )+
  
  scale_color_manual(
    values=c(
      
      "Others"="grey80",
      
      "Branch after Y_114"="#2C7BB6",
      
      "Branch after Y_98"="#D73027",
      
      "Branch point"="black"
      
    )
  )+
  
  coord_equal()+
  
  theme_classic(base_size = 8,base_family = "Arial")+
  
  theme(
    
    axis.text=element_blank(),
    
    axis.ticks=element_blank(),
    
    axis.title=element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    
    legend.position="right"
    
  )
ggsave(plot = last_plot(),
       filename = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou17/trajectory_nodelabel.pdf",
       width = 9,height = 6)

save.image("./trajectory_nodelabel.Rdata")
