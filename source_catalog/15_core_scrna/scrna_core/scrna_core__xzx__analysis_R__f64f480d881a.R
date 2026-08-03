library(monocle)
setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou8")
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

################################################################################
### 需求2 cellchat
rm(list = ls())
gc()
library(dplyr)    # 数据筛选与分组统计
library(ggplot2)  # 数据可视化
library(tidyr)    # 数据格式辅助处理

#load("../shouhou6/cellchat/pos_neg/tmp_object.RData")

pos <- read.table("../shouhou6/cellchat/pos_neg/pos-netPathways_lr.txt",sep = '\t',header = TRUE)
#pos_lr <- read.table("../shouhou6/cellchat/pos_neg/pos-net_lr.txt",sep = '\t',header = TRUE)
neg <- read.table("../shouhou6/cellchat/pos_neg/neg-netPathways_lr.txt",sep = '\t',header = TRUE)

pos$group <- "pos"
neg$group <- "neg"

pathway <- rbind(pos,neg)
  
# 筛选source或target列中包含"tumor_1"的行
# grepl()实现模糊匹配，|表示“或”逻辑，filter()筛选符合条件的行
filtered_df <- pathway %>%
  filter(grepl("tumor_1", source) | grepl("tumor_1", target))

# 注意：模拟数据中无"tumor_1"，替换为你的真实数据后会得到有效筛选结果

# 按pathway和group分组，统计数量和prob总和
# group_by()分组，summarise()聚合统计，n()计算每组行数（数量），sum()计算prob总和
stats_df <- filtered_df %>%
  group_by(pathway_name, group) %>%
  summarise(
    count = n(),  # 统计每组的行数（数量）
    prob_sum = sum(prob, na.rm = TRUE),  # 统计每组的prob总和，忽略缺失值
    .groups = "drop"  # 取消分组状态，返回普通数据框
  )

# 5. 绘制柱状图（分两张图：数量柱状图 + prob总和柱状图）
# 5.1 绘制“不同pathway&group的数量柱状图”
p1 <- ggplot(stats_df, aes(x = pathway_name, y = count, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +  # dodge表示分组并列显示，不堆叠
  scale_fill_brewer(palette = "Set2") +  # 设置配色方案
  labs(
    title = "",
    x = "",
    y = "Number of interactions"
  ) +
  theme_classic() +  # 白色背景主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 横坐标文字旋转，避免重叠
p1
ggsave(plot = p1,filename = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou8/cellchat.png",width = 12,height = 6,dpi = 300)
write.csv(filtered_df,"./tumor1_pathway.csv",row.names = FALSE)

# 5.2 绘制“不同pathway&group的prob总和柱状图”
p2 <- ggplot(stats_df, aes(x = pathway_name, y = prob_sum, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Sum of Prob per Pathway and Group",
    x = "Pathway Name",
    y = "Sum of Prob"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5.3 并列显示两张图（需提前安装patchwork包：install.packages("patchwork")）
library(patchwork)
p1 + p2
