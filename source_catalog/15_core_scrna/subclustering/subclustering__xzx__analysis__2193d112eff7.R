setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/")

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

load("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/mono.Rdata")

#### 添加新的celltype分类
cnv = read.csv(file = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou15/meta.csv',row.names = 1)
cds@phenoData@data$cell_type = cnv[rownames(cds@phenoData@data),'cell_type']#
cds@phenoData@data$inferCNV = cnv[rownames(cds@phenoData@data),'inferCNV']
table(cds@phenoData@data$cell_type)

plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)

cell_type_colors <- c(
  "BGN+ malignant Epi" = "#1f77b4",
  "CCNL2+ malignant Epi" = "#ff7f0e",
  "KRT13+ malignant Epi" = "#2ca02c",
  "SPRR1B+ malignant Epi" = "#d62728",
  "Apoptosis Epithelium" = "#9467bd",
  "MUC5B+ Epithelium" = "#8c564b"
)


p <- plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)+
  facet_wrap(~cell_type)+
  scale_color_manual(values = cell_type_colors)
p
ggsave(filename = "./PPT_2/pseudotime_celltype_facet.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "./PPT_2/pseudotime_celltype_facet.pdf",p,h=7,w=10)



p <- plot_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type',cell_size = .5)+
  scale_color_manual(values = cell_type_colors)
p
ggsave(filename = "./PPT_2/pseudotime_celltype.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "./PPT_2/pseudotime_celltype.pdf",p,h=7,w=10)

################################################################################
# 读入 CytoTRACE 结果（默认是 tab/空格分隔）
cytotrace <- read.delim(
  file = "./CytoTRACE.value.txt", 
  header = TRUE, 
  sep = "", # 自动识别分隔符
  row.names = 1, # 把第一列（Cell）设为行名
  check.names = FALSE
)

# 检查读入是否正确
head(cytotrace)

cds@phenoData@data$CytoTRACE_value = cytotrace[rownames(cds@phenoData@data), 'CytoTRACE_value']

p <- plot_cell_trajectory(cds,x=1,y=2,color_by = 'CytoTRACE_value',cell_size = .5)+
  scale_color_viridis_c(option = "plasma")
p
ggsave(filename = "./PPT_2/pseudotime_CytoTRACE_value.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "./PPT_2/pseudotime_CytoTRACE_value.pdf",p,h=7,w=10)


################################################################################
sc <- readRDS("./adata_epi_anno.rds")
DimPlot(sc,group.by = 'cell_type')

################################################################################
# ===================== 1. 加载必需包 =====================
library(tidyverse)  # 包含数据处理(dplyr)、绘图(ggplot2)、文件读取(readr)

# ===================== 2. 设置文件路径（仅需修改这里！） =====================
# 核心路径：定位到 recluster 文件夹（cluster0~5的上级目录）
main_path <- "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/Epithelium/"
# 文件名（所有cluster下的文件名一致）
file_name <- "GO_Biological_Process_2023.Human.enrichr.reports.txt"

# ===================== 3. 批量读取6个cluster的GO结果 =====================
# 创建空列表存储数据
all_data <- list()

# 循环读取 cluster0 ~ cluster5
for (i in 0:5) {
  # 拼接每个文件的完整路径
  cluster_path <- file.path(main_path, paste0("cluster", i), file_name)
  
  # 读取txt文件（enrichr输出的制表符分隔文件）
  df <- read.delim(cluster_path, sep = "\t", stringsAsFactors = FALSE)
  
  # 添加分组列（标记cluster编号）
  df$cluster <- paste0("cluster", i)
  
  # 存入列表
  all_data[[i+1]] <- df
}

# 合并所有数据
combined_data <- bind_rows(all_data)

combined_data <- combined_data %>%
  mutate(cluster = recode(cluster,
                          "cluster0" = "BGN+ malignant Epi",
                          "cluster1" = "CCNL2+ malignant Epi",
                          "cluster2" = "KRT13+ malignant Epi",
                          "cluster3" = "SPRR1B+ malignant Epi",
                          "cluster4" = "Apoptosis Epithelium",
                          "cluster5" = "MUC5B+ Epithelium"))

# ===================== 4. 数据筛选：Adjusted P-value <0.05 + 每组TOP5 =====================
# 关键列名说明（enrichr标准输出）：
# Adjusted.P.value = 校正后P值 | Odds.Ratio = 优势比 | Term = GO功能条目
final_data <- combined_data %>%
  # 1. 筛选显著条目：Adjusted P-value < 0.05，每组取TOP5（按P值升序，最显著在前）
  dplyr::filter(Adjusted.P.value < 0.05) %>%
  group_by(cluster) %>%
  arrange(Adjusted.P.value, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  
  # 2. 【核心1】彻底删除Term中的(GO:xxxx)编号
  mutate(
    Term = str_remove(Term, " \\(GO:\\d+\\)"), # 正则匹配并删除" (GO:数字)"
    Term = str_trim(Term) # 清除多余空格
  )

# 3. 【核心2】按横轴细胞类型顺序排列纵轴Term
# 先按cluster顺序（normal→tumor1→tumor2→tumor3→tumor4→tumor5）+ P值升序排序
term_order <- final_data %>%
  mutate(cluster = factor(cluster, levels = c("MUC5B+ Epithelium", "BGN+ malignant Epi", 
                                              "CCNL2+ malignant Epi", "KRT13+ malignant Epi", 
                                              "SPRR1B+ malignant Epi", "Apoptosis Epithelium")
  )) %>%
  arrange(cluster, Adjusted.P.value) %>%
  pull(Term) %>%
  unique() # 去重，避免重复term

# 将Term设置为因子，固定顺序（反转后normal在最上方，符合阅读习惯）
final_data <- final_data %>%
  mutate(
    Term = factor(Term, levels = rev(term_order)),
    cluster = factor(cluster, levels = c("MUC5B+ Epithelium", "BGN+ malignant Epi", 
                                         "CCNL2+ malignant Epi", "KRT13+ malignant Epi", 
                                         "SPRR1B+ malignant Epi", "Apoptosis Epithelium")
    ) # 固定横轴顺序
  )

# ===================== 5. 科研风美化气泡图绘制 =====================
p <- ggplot(final_data, 
            aes(x = cluster, 
                y = Term, 
                size = Odds.Ratio, 
                fill = Adjusted.P.value)) + # 填充色映射P值
  # 气泡样式：带黑色边框+半透明，避免重叠，更清晰
  geom_point(
    shape = 21, # 带边框的圆形
    alpha = 0.9, 
    stroke = 0.8, # 边框粗细
    color = "black" # 边框统一黑色
  ) +
  
  # 颜色渐变：P值越小（越显著）→ 颜色越深（红→蓝，符合科研习惯）
  scale_fill_gradient(
    low = "#D73027",   # 低P值（显著）：深红色
    high = "#3274A1",  # 高P值（接近0.05）：蓝色
    limits = c(0, 0.05), # 固定颜色标尺，统一标准
    name = "Adjusted P-value",
    guide = guide_colorbar(reverse = TRUE) # 反转图例，0.00在最上方，更直观
  ) +
  
  # 气泡大小：映射Odds Ratio，调整范围让差异更明显
  scale_size_continuous(
    range = c(4, 12), # 最小4，最大12（可根据数据调整）
    name = "Odds Ratio",
    breaks = c(20, 40, 60, 80) # 对应图例刻度，和原图一致
  ) +
  
  # 纵轴标签自动换行，避免过长拥挤
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  
  # 图表标签与标题
  labs(
    x = "Cell Cluster",
    y = "GO Biological Process",
    title = "GO Biological Process Enrichment Analysis"
    #subtitle = "Top 5 Significant Terms (Adjusted P-value < 0.05)"
  ) +
  
  # 主题美化：极简科研风，符合期刊要求
  theme_minimal() +
  theme(
    # 字体统一为Arial（绝大多数期刊要求）
    text = element_text(size = 10,face = 'bold'),
    # 标题样式：居中、加粗、放大
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 10)),
    # 坐标轴标签样式
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, hjust = 1),
    # 网格线：仅保留纵向网格，横向删除，更简洁
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # 图例样式：右侧垂直排列，调整大小
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),
    # 边距调整，避免标签被截断
    plot.margin = margin(1, 1.5, 1, 1.5, "cm")
  ) +
  
  # 【可选】替换为Nature子刊配色（更符合高分期刊风格）
  # scale_fill_npg(guide = guide_colorbar(reverse = TRUE), name = "Adjusted P-value") +
  
  NULL

# 显示图表
print(p)

# ===================== 6. 保存高清图片 =====================
ggsave("./PPT_2/GOBP_dotplot.pdf", p, width = 16, height = 10, dpi = 300)
ggsave("./PPT_2/GOBP_dotplot.png", p, width = 16, height = 10, dpi = 300,bg = 'white')


#################################################################################
gene_list <- c(  "PHGDH", 
                 "SLC1A4", "MTHFD1", "DNMT1", "TYMS", 
                 "AHCY", "ATIC", "DNMT3A")

p <- VlnPlot(sc, 
             pt.size=0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = T, # 翻转
             group.by = 'cell_type',
             split.by = 'group',
             adjust=3, # 调整粗细
             fill.by = "ident", 
             cols = c('#3274A1','#D73027'), 
             features = gene_list) + 
  labs(x = "") + 
  theme(#legend.position = 'none', 
    axis.text.x = element_text(size = 12, color = "black",angle = 0,hjust = 0.5),
    axis.title.y = element_blank(),
    strip.text.y = element_text(hjust = 0))
p

ggsave(p, file = "./Fig4/Fig4_ocm_vlnplot.pdf", width = 16, height = 9)
ggsave(p, file = "./Fig4/Fig4_ocm_vlnplot.png", width = 16, height = 9, dpi = 300,bg = 'white')
#################################################################################

gene_list <- c(
  "BGN","LUM","SDC2",
  'CCNL2','SRSF5','PNISR',
  'KRT13','KRT15','SOX2',
  'SPRR1B','KRTDAP','CNFN',
  "NQO1","ALDH3A1","SPP1",  #"CAPN1", "CAPN2", "NFKBIA","PRKACA", "PRKAR1A", "TNFRSF1A"
  "MUC5B","PIGR",'FCGBP'
  )

p <- VlnPlot(sc, 
             pt.size=0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             #flip = T, # 翻转
             group.by = 'cell_type',
             #split.by = 'group',
             adjust=5, # 调整粗细
             fill.by = "ident", 
             #cols = c('#3274A1','#D73027'), 
             features = gene_list) + 
  labs(x = "") + 
  theme(#legend.position = 'none', 
    axis.text.x = element_text(size = 12, color = "black",angle = 45,hjust = 0.5),
    axis.title.y = element_blank(),
    strip.text.y = element_text(hjust = 0))+
  scale_fill_manual(values = cell_type_colors)
p
ggsave(p, file = "./PPT_2/Epimarker_vlnplot.pdf", width = 16, height = 9)
ggsave(p, file = "./PPT_2/Epimarker_vlnplot.png", width = 16, height = 9, dpi = 300,bg = 'white')


sc_scale <- ScaleData(sc, features = gene_list)

# 3. 绘制热图
p <- DoHeatmap(
  object = sc_scale,
  features = gene_list,
  group.by = "cell_type",
  group.colors = cell_type_colors,
  label = FALSE,
  size = 3,
  angle = 0
) +
  scale_fill_gradientn(colors = c("navy", "gray", "firebrick3")) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12,face = 'bold')
  )
p
ggsave(p, file = "./PPT_2/Epimarker_heatmap.pdf", width = 16, height = 10)
ggsave(p, file = "./PPT_2/Epimarker_heatmap.png", width = 16, height = 10, dpi = 300,bg = 'white')

p <- DoHeatmap(
  object = sc_scale,
  features = gene_list,
  group.by = "cell_type",
  group.colors = cell_type_colors,
  label = TRUE,
  size = 3,hjust = 0.5,
  angle = 0
) +
  scale_fill_gradientn(colors = c("navy", "gray", "firebrick3")) +
  theme(
    axis.line.x = element_blank(),  # 移除x轴黑线
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )
p
ggsave(p, file = "./PPT_2/Epimarker_heatmap.pdf", width = 18, height = 9,dpi = 300)
ggsave(p, file = "./PPT_2/Epimarker_heatmap.png", width = 18, height = 9, dpi = 300,bg = 'white')

################################################################################
p <- VlnPlot(sc,features = 'GJA1',pt.size = 0,group.by = 'cell_type')+
  scale_fill_manual(values = cell_type_colors)+
  theme(axis.title.x = element_blank())
p
ggsave(plot = p,filename = "./GJA1_exp_celltype.png",
       width = 12,height = 6,dpi = 300,bg = "white")
ggsave(plot = p,filename = "./GJA1_exp_celltype.pdf",
       width = 12,height = 6,dpi = 300,bg = "white")



################################################################################
sc_fib <- readRDS('/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Fibroblasts/celltype/h5_trans.rds')

VlnPlot(sc_fib,features = c('JAG1','GJA1','NOTCH1','NOTCH3'),group.by = 'group',pt.size = 0,log = TRUE) &
  stat_compare_means(
    comparisons = list(c("pos", "neg")),
    method = "wilcox.test",
    label = "p.signif",
    #hide.ns = TRUE,  # 隐藏不显著的标记
    size = 3,
    label.y.npc = 0.95  # 用“相对坐标”：0=底，1=顶，0.95=靠近顶部
  ) &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))  # 上方多留20%空白


gene_list <- c('JAG1','GJA1','NOTCH1','NOTCH3')
for (i in gene_list) {
  p <- VlnPlot(sc_fib,features = i,group.by = 'group',pt.size = 0) &
    stat_compare_means(
      comparisons = list(c("pos", "neg")),
      method = "wilcox.test",
      label = "p.signif",
      #hide.ns = TRUE,  # 隐藏不显著的标记
      size = 3,
      label.y.npc = 0.95  # 用“相对坐标”：0=底，1=顶，0.95=靠近顶部
    ) &
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
  
  ggsave(paste0("./",i,"_Fibroblast_group_vln.png"),p,width =5,height = 5,bg = 'white')
  ggsave(paste0("./",i,"_Fibroblast_group_vln.pdf"),p,width =5,height = 5,bg = 'white')
  
}
################################################################################
sc_cellchat <- readRDS('/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou6/cellchat.rds')
# 1. 定义你要改的映射
map <- c(
  "tumor_1" = "BGN+ malignant Epi",
  "tumor_2" = "CCNL2+ malignant Epi",
  "tumor_3" = "KRT13+ malignant Epi",
  "tumor_4" = "SPRR1B+ malignant Epi",
  "tumor_5" = "Apoptosis Epithelium",
  "normal"  = "MUC5B+ Epithelium"
)

# 2. 转字符，避免因子level报错、丢细胞
sc_cellchat$cell_type <- as.character(sc_cellchat$cell_type)

# 3. 只替换 map 里有的，其他不变
for (old in names(map)) {
  sc_cellchat$cell_type[sc_cellchat$cell_type == old] <- map[old]
}

# 4. 转回因子，并更新 Idents（绘图用）
sc_cellchat$cell_type <- factor(sc_cellchat$cell_type)
Idents(sc_cellchat) <- "cell_type"

# 5. 检查：所有细胞都在，只有6个被改名
table(sc_cellchat$cell_type)
saveRDS(sc_cellchat,"/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou16/cellchat.rds")
################################################################################
setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou16/")

gene_name <- c("BGN")
# 提取表达值并添加到 metadata
# FetchData 返回的是一个数据框，列名就是基因名
sc[[gene_name]] <- FetchData(sc, vars = gene_name)

for (i in gene_name) {
  sc[[i]] <- FetchData(sc, vars = i)
}

meta <- as.data.frame(sc@meta.data)
meta <- meta[,18:19]
write.csv(meta,"./gene_meta.csv")

colnames(meta)

cds@phenoData@data[['BGN']] <- meta[rownames(cds@phenoData@data),'BGN']

p = plot_cell_trajectory(cds,x=1,y=2,color_by = 'BGN',cell_size = .5)+
  scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                        values = c(0,0.5,1,1.5))
p
ggsave(filename = "monocle_pseudotime_BGN.png",p,h=7,w=10,dpi = 300)
ggsave(filename = "monocle_pseudotime_BGN.pdf",p,h=7,w=10,dpi = 300)
