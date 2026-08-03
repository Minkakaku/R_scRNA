setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou11")

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

rm(list = ls())

################################################################################
### 需求1 monocle分支明确一下
load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

sc <- readRDS("../shouhou4/Epithelial_all/epi_trans.rds")
DimPlot(sc,group.by = 'cell_type')

FeaturePlot(sc,features = c("SHMT1","SHMT2","MTHFD2","PHGDH","PSAT1","PSPH","SFXN3","SLC1A4","SLC1A5","FOLH1","GJA1","JUNB"))

gene_name <- c("SHMT1","SHMT2","MTHFD2","PHGDH","PSAT1","PSPH","SFXN3","SLC1A4","SLC1A5","FOLH1","GJA1","JUNB")
# 提取表达值并添加到 metadata
# FetchData 返回的是一个数据框，列名就是基因名
sc[[gene_name]] <- FetchData(sc, vars = gene_name)

for (i in gene_name) {
  sc[[i]] <- FetchData(sc, vars = i)
}

meta <- as.data.frame(sc@meta.data)
meta <- meta[,16:27]
write.csv(meta,"./gene_meta.csv")

load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

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

library(msigdbr)
# 提取人类 Hallmark 集合中的 EMT 基因
emt_genes <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

emt_genes <- emt_genes$human_gene_symbol

genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 100,
    nbin = 5,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    Monocle_State = sc@meta.data[[group_name]],  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = Monocle_State, y = score, fill = Monocle_State)) +
    geom_violin(scale = "area") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(gene_name," Module Score"), x = NULL, y = "Score")
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.png"),
         width = 6,height = 6,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.pdf"),
         width = 6,height = 6,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}


##
#sc <- readRDS('./B/celltype/adata_anno.rds')

EMT_Gene <- emt_genes
sc <- genesets_addmodulescore(gene_sets = EMT_Gene,gene_name = "EMT_Gene",
                              scRNA = sc,group_name = "Monocle_State",out_dir = "./")



cnv = read.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2507233889-01/895/celltype/CNV/metadata.csv',row.names = 1)
HSMM_myo = readRDS('./895/3.obj_reorder.rds')
HSMM_myo@phenoData@data$cnv = cnv[rownames(HSMM_myo@phenoData@data),'cnv_score_cluster']#
p = plot_cell_trajectory(HSMM_myo,x=1,y=2,color_by = 'cnv',cell_size = .5)+
  facet_wrap(~cell_type)+scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                                               values = c(0,0.5,0.8,1))
ggsave(filename = "0.895_monocle_pseudotime_celltype_infercnv.png",p,h=7,w=10,dpi = 150)




# 1. 准备库
library(Seurat)
library(ggplot2)
library(ggpubr) # 用于自动添加统计 P 值
library(dplyr)

# 假设你的 Seurat 对象名为 sce
# metadata 中应包含：'sample_id' (6个ID) 和 'group' (HIV+ / HIV-)

library(Seurat)
library(ggplot2)
library(ggpubr) # 用于添加 P 值
library(dplyr)
library(tidyr) # 用于 replace_na

# 假设你的 Seurat 对象名为 sc

# ==============================================================================
# 第一步：准备数据（从细胞水平到样本水平）
# ==============================================================================

# 1. 提取 meta.data
meta_data <- sc@meta.data

# 确保 tumor_1 在 cell_type 列中的拼写完全正确
# 检查：table(meta_data$cell_type)

## --- 1. 计算每个样本的总细胞数 ---
# 结果表格：sample | group | total_cells
sample_totals <- meta_data %>%
  group_by(sample, group) %>%
  summarise(total_cells = n(), .groups = "drop")

## --- 2. 计算每个样本中 'tumor_1' 的细胞数 ---
# 结果表格：sample | tumor_1_cells
tumor1_counts <- meta_data %>%
  # *** 核心筛选：只保留 tumor_1 亚群 ***
  filter(cell_type == "tumor_1") %>%
  group_by(sample) %>%
  summarise(tumor_1_cells = n(), .groups = "drop")

## --- 3. 合并数据并计算比例 ---
# 结果表格：sample | group | proportion
plot_data <- sample_totals %>%
  left_join(tumor1_counts, by = "sample") %>%
  # 计算比例
  mutate(proportion = tumor_1_cells / total_cells)

# 检查生成的数据框，确认是 6 行，包含了所有样本、分组和比例
# print(plot_data)

# ==============================================================================
# 第二步：按图需求进行可视化
# ==============================================================================

# 创建图形
ggplot(plot_data, aes(x = group, y = proportion, fill = group)) +
  # 1. 绘制箱线图展示分布（alpha设为0.5让它变透明）
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) + 
  # 2. 叠加热图/散点图，明确展示那 6 个病人样本的数据点
  geom_jitter(width = 0.1, size = 3, color = "black") + 
  # 3. 设置颜色 (假设你的分组名称是 pos / neg)
  scale_fill_manual(values = c("neg" = "#56B4E9", "pos" = "#E69F00"),
                    labels = c("neg" = "neg", "pos" = "pos")) +
  # 4. 添加统计学检验 P 值 (Mann-Whitney / Wilcoxon)
  stat_compare_means(method = "wilcox.test",paired = TRUE,
                     label = "p.format", 
                     label.x = 1.4) + # 调整P值位置
  # 5. 美化
  theme_classic() +
  labs(#title = "Proportion of tumor_1 per Patient",
       y = "Proportion (tumor_1 / Total)",
       x = "Group") +
  theme(legend.position = "right",
        # 如果需要 Panel 占很小位置，可以设置纵横比
        aspect.ratio = 1.2)



# 1. 优化 X 轴排序（可选但推荐）：
# 如果你希望 neg 的三个样本排在一起，pos 的三个排在一起
plot_data <- plot_data %>%
  arrange(group, sample) %>% # 先按组排，再按ID排
  mutate(sample = factor(sample, levels = sample)) # 将 sample 列设为 factor 以保持排序

# --- 1. 手动计算 Mann-Whitney P 值 ---
# 提取两组的比例向量
group_neg <- plot_data$proportion[plot_data$group == "neg"]
group_pos <- plot_data$proportion[plot_data$group == "pos"]

# 执行检验
stats_result <- wilcox.test(group_neg, group_pos,paired = TRUE)
p_val <- round(stats_result$p.value, 3)

# --- 2. 绘图并在上方添加 P 值大括号 ---
ggplot(plot_data, aes(x = sample, y = proportion, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("neg" = "#56B4E9", "pos" = "#E69F00")) +
  theme_classic() +
  # 添加跨组的大括号和 P 值
  # xmin 为第1个柱子，xmax 为第6个柱子
  geom_bracket(
    xmin = 1, xmax = 6, 
    y.position = max(plot_data$proportion) * 1.1, # 位于最高柱子上方10%
    label = paste0("p = ", p_val),
    tip.length = 0.03,
    inherit.aes = FALSE
  ) +
  labs(#title = "Tumor_1 Proportion per Patient",
       subtitle = "Group comparison: Mann-Whitney U test",
       y = "Proportion", x = "Patient ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "Tumor_1_Proportion.png",plot = last_plot(),h=7,w=7,dpi = 300)
