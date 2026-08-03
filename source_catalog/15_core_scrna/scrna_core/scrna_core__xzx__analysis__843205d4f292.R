dir.create("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou12")
setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou12")

rm(list = ls())
library(igraph)
library(monocle)
library(Seurat)
library(ggplot2)
#library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

################################################################################
### 需求1 GJA1 各tumor亚群的组间表达箱型图

load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

sc <- readRDS("../shouhou4/Epithelial_all/epi_trans.rds")
DimPlot(sc,group.by = 'cell_type')

FeaturePlot(sc,features = c("GJA1"))
VlnPlot(sc,features = 'GJA1',group.by = 'cell_type',split.by = 'group')


# 提取表达矩阵 + 元数据
gene <- "GJA1"
expr_data <- FetchData(sc, vars = c(gene, "cell_type", "group"))

# ggplot2 绘制箱型图
ggplot(expr_data, aes(x = cell_type, y = .data[[gene]], fill = group)) +
  # 箱型图，分组并排显示
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.size = 0.5) +
  # 主题与标签
  labs(x = NULL, y = "GJA1 Expression", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = rel(1.2)),
        plot.title = element_text(hjust = 0.5))

ggsave(plot = last_plot(),"./GJA1_celltype_group_box.png",device = 'png',width = 9,height = 6,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"./GJA1_celltype_group_box.pdf",device = 'pdf',width = 9,height = 6,dpi = 300,bg = 'white')


### 需求2 Tumor1中各monocle state的比例图

# 1. 提取 Monocle 对象的 State 信息
# pData(cds) 提取的是 metadata 数据框
monocle_states <- pData(cds)$State

# 2. 为提取出的 State 命名（确保能和 Seurat 匹配）
names(monocle_states) <- colnames(cds)

# 3. 将 State 添加到 Seurat 对象的 metadata 中
# Seurat 会自动根据 names(monocle_states) 匹配细胞 ID
sc[["Monocle_State"]] <- monocle_states[colnames(sc)]

# 4. 验证是否成功
DimPlot(sc, group.by = "Monocle_State", label = TRUE)

# 1. 从Seurat对象提取元数据
meta_data <- sc@meta.data  # sc 是你的Seurat对象

# 2. 统计：每个细胞类型 → 各state的数量 & 百分比
plot_df <- meta_data %>%
  group_by(cell_type, Monocle_State) %>%  # 按细胞类型+state分组
  summarise(count = n(), .groups = "drop") %>%  # 计数
  group_by(cell_type) %>%
  mutate(percent = count / sum(count) * 100)  # 计算占比(%)

# 3. 绘制百分比堆叠柱状图
ggplot(plot_df, 
       aes(x = cell_type,           # X轴：细胞类型
           y = percent,             # Y轴：百分比
           fill = factor(Monocle_State))) +  # 填充色：monocle状态
  geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
  scale_y_continuous(labels = scales::percent) +  # Y轴显示%
  # 配色（可自定义）
  scale_fill_brewer(palette = "Paired") +
  # 标签与主题
  labs(x = "Cell Type", 
       y = "Proportion of Monocle State", 
       fill = "Monocle State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
ggsave(plot = last_plot(),"./celltype_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"./celltype_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')


### 需求3 单碳代谢的这些基因 在pos tumor中的 表达score UMAP

SGOC <- c(
  "DNMT3B", "SHMT1", "MTHFD2", "PHGDH", "PSAT1", "PSPH", "SFXN3",
  "SLC1A4", "SLC1A5", "SHMT2", "MTHFD1", "DNMT1", "TYMS", "MTHFD1L",
  "AHCY", "MTFMT", "MTHFS", "ATIC", "GART", "MAT2B", "MTR", "MTHFR",
  "GLRX", "TCN2", "ALDH1L1", "FTCD", "BHMT", "DHFR", "DNMT3A", "CHDH",
  "FOLH1"
)

SGOC <- list(SGOC)
sc <- AddModuleScore(
  object = sc,
  features = SGOC,
  ctrl = 100,
  nbin = 5,
  name = "SGOC"
)

# 对模块分数进行 Z-score 标准化
sc@meta.data$SGOC1_scaled <- scale(sc@meta.data$SGOC1)[,1]

# 查看标准化后的范围
summary(sc@meta.data$SGOC1_scaled)

FeaturePlot(
  object = sc,
  features = "SGOC1_scaled",        # 必须加 1！AddModuleScore 自动生成
  pt.size = 0.3,keep.scale = 'all'
) +
  # 美化
  scale_color_gradientn(
    colors = c("#0099B4","white","lightgrey", "red")  # 经典配色
  ) +
  ggtitle("SGOC Score") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# 2. 绘制漂亮的UMAP图
FeaturePlot(sc, features = "SGOC1_scaled", pt.size = 0.3) +
  scale_color_gradient2(
    low = "#0099B4",
    
    mid = "white",
    
    high = "#E41A1C",
    midpoint = 0,
    limits = c(-2, 3)  # 颜色范围锁定在-2~2
  ) +
  ggtitle("SGOC Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = last_plot(),"./UMAP_SGOCscore.png",device = 'png',width = 6,height = 6,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"./UMAP_SGOCscore.pdf",device = 'pdf',width = 6,height = 6,dpi = 300,bg = 'white')
