library(ggplot2)
library(Seurat)
library(ggpubr)
library(viridis)  # 高级配色
library(scales)   # 坐标轴格式化
library(tidyr, help, pos = 2, lib.loc = NULL)

setwd("/mnt//NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/")

# 读入seurat对象
sce <- readRDS("/mnt//NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/Anno//Epith/adata_anno.rds")
# 打分的基因集 ---
# 1. 耗竭，效应，活化
# 2. 3个代谢功能 糖酵解，氧化磷酸化，氨基酸代谢

table <- read.delim2("/mnt//NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/Anno//Epith_score.csv",
                     header = T, check.names = F, sep = ',')

gene_signatures <- list()
for(coll in colnames(table)){
    # 【关键修改】：剔除空字符串 "" 和 NA
    genes <- table[, coll]
    genes <- genes[genes != "" & !is.na(genes)] 
    gene_signatures[[coll]] <- genes
}
sce1 <- sce
# 确保在进入循环前，gene_signatures 里的空字符串已经被清理
for(i in 1:length(gene_signatures)){
    sig_name <- names(gene_signatures)[i]
    
    # 【关键点 1】：清理当前基因列表中的空字符串和 NA
    # 否则 AddModuleScore 会去搜索名为 "" 的基因，导致匹配失败
    current_genes <- gene_signatures[[i]]
    current_genes <- current_genes[current_genes != "" & !is.na(current_genes)]
    
    # 检查交集数量，用于调试
    common <- intersect(rownames(sce1), current_genes)
    print(paste0("Processing ", sig_name, ": Found ", length(common), " genes"))
    
    # 如果交集基因太少（比如少于 2 个），AddModuleScore 可能会报错
    if(length(common) < 2) {
        message(paste0("跳过 ", sig_name, "，因为匹配到的基因太少。"))
        next
    }
    
    # 【关键点 2】：features 必须传 list(current_genes)
    # 虽然 gene_signatures 本身是 list，但 gene_signatures[[i]] 提取出来后是 vector
    # AddModuleScore 要求的是 list of vectors
    sce1 <- AddModuleScore(
        object = sce1,
        features = list(current_genes), 
        name = sig_name
    )
    
    # 找到最近添加的那一列
    last_col <- ncol(sce1@meta.data)
    colnames(sce1@meta.data)[last_col] <- sig_name
}

meta_data <- sce1@meta.data
# 模块评分列名
score_cols <- names(gene_signatures)
# 转换为长格式
long_data <- pivot_longer(meta_data,
                          cols = all_of(score_cols),
                          names_to = "Module",
                          values_to = "Score")
# 清理模块名称
long_data$Module <- gsub("1", "", long_data$Module, perl = TRUE)
long_data$group <- factor(long_data$group,levels = c("CON","KO"))

# 1. 手动计算显著性 (Wilcoxon test)
# 我们创建一个小表格来存放每个 Module 的显著性标签
library(dplyr)

sig_labels <- long_data %>%
  group_by(Module) %>%
  summarise(
    p_val = wilcox.test(Score ~ group)$p.value,
    # 根据 p 值转换成星号
    label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    # 确定标签显示的高度（每个分面最大值之上一点）
    y_pos = max(Score) * 1.1 
  )

# 2. 绘图
p_publication <- ggplot(long_data, aes(x = group, y = Score, fill = group)) +
  # 保持原有的提琴图和箱线图逻辑
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, fill = "white") +
  
  # 【替代 stat_compare_means】：使用 geom_text 添加星号
  geom_text(data = sig_labels, 
            aes(x = 1.5, y = y_pos, label = label), 
            inherit.aes = FALSE, 
            size = 6, fontface = "bold") +
  
  # 分面
  facet_wrap(~ Module, ncol = 3, scales = "free_y") +
  
  # 颜色与主题 (保持你原有的出版级设置)
  scale_fill_manual(values = c("CON"='#2B739F', "KO" = "#E68525")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "top",
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "", y = "Module Score", title = "Function Scores by Group")

print(p_publication)

ggsave( "Epith_score/Epith_function_Group.pdf",
        p_publication, width = 10, height = 10, dpi = 600)

# 各个细胞类型
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. 提取并定义图片中的颜色
cell_type_colors <- c(
  'EMT'= '#a93a2b', 
  'Proliferative' = '#732d28', 
  'Stress' = '#75a8c4'
)

# 1. 计算多组间整体差异 (Kruskal-Wallis)
stats_data <- long_data %>%
  group_by(Module) %>%
  summarise(
    # 使用 kruskal.test 处理多组比较
    p_val = kruskal.test(Score ~ cell_type)$p.value, 
    .groups = 'drop'
  ) %>%
  mutate(
    label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

p_celltype_diff <- ggplot(long_data, aes(x = cell_type, y = Score, fill = cell_type)) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  
  # 在每个分面的顶部展示该 Module 的整体显著性 (P.value)
  geom_text(data = stats_data, 
            aes(x = 2, y = Inf, label = paste0("P-val: ", label)), 
            vjust = 2, inherit.aes = FALSE, size = 4.5, fontface = "italic") +
  
  facet_wrap(~ Module, scales = "free_y", ncol = 3) +
  
  # 使用提取自 UMAP 的颜色
  scale_fill_manual(values = cell_type_colors) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    panel.spacing = unit(1.2, "lines")
  ) +
  labs(x = "Cell Types", y = "Module Score", title = "Functional Scores Across Cell Types")

print(p_celltype_diff)

ggsave( "Epith_score/Epith_function_cellType.pdf",
        p_celltype_diff, width = 10, height = 6, dpi = 600)


saveRDS(sce1,"Epith_score/Epith_score.rds")
