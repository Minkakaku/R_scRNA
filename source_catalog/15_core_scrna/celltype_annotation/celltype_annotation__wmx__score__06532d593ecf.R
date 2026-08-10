# 对于genelist里的每个list进行单细胞addmodulescore
# 对于每个list的打分进行leiden对比小提琴图绘制，进行差异对比
# 对于打分绘制umap图，一张图展示

rm(list = ls())
rds_file <- "D:/Project/SC2508154388-01/0526/mye_final/score.rds"
geneset_file <- "D:/Project/SC2508154388-01/0526/mye_final/geneset.csv"
outdir <- "D:/Project/SC2508154388-01/0526/mye_final/Score/"
dir.create(outdir,showWarnings = F,recursive = T)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(grid)

set.seed(123)

# 读入文件 ------
obj <- readRDS(rds_file)

gs_raw <- read.csv(
  geneset_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 按列转成 list
pathway_list <- lapply(colnames(gs_raw), function(pw) {
  genes <- gs_raw[[pw]]
  genes <- genes[!is.na(genes)]
  genes <- trimws(as.character(genes))
  genes <- genes[genes != ""]
  unique(genes)
})
names(pathway_list) <- colnames(gs_raw)

## =========================
# AddmoduleScore 打分
## =========================
for(ss in names(pathway_list)){
obj <- AddModuleScore(
  object = obj,
  features = pathway_list[ss],
  name = ss
)
}

colnames(obj@meta.data) <- gsub(pattern = "1$",x = colnames(obj@meta.data),replacement = "")

## =========================
# 函数计算差异、颜色，主题定义
## =========================
p_to_label <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  return("ns")
}

theme_pub <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.position = "top",
      legend.text = element_text(size = 16,face = 'bold'),
      legend.title = element_text(size = 16,face = 'bold'),
      title = element_text(size = 18,face = 'bold'),
      panel.spacing = unit(1, "lines")
    )
}

## =========================
## Module score 小提琴图 + 差异统计
## =========================
group_col <- 'leiden'

score_long <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_col), all_of(names(pathway_list))) %>%
  pivot_longer(
    cols = all_of(names(pathway_list)),
    names_to = "Score_col",
    values_to = "Score"
  ) 

p_score_violin <- ggplot(score_long, aes(x = .data[[group_col]], y = Score, fill = .data[[group_col]])) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, fill = "white") +
  facet_wrap(~ Score_col, ncol = 3, scales = "free_y") +
  # scale_fill_manual(values = group_colors) +
  theme_pub() +
  labs(x = "", y = "Module Score", title = "Function Scores by leiden")

ggsave(
  file.path(outdir, "ModuleScore_violin_by_leiden.pdf"),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 3.5)
)

ggsave(
  file.path(outdir, "ModuleScore_violin_by_leiden.png"),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 3.5),
  dpi = 300
)

group_col <- 'cell_type'
score_long <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_col), all_of(names(pathway_list))) %>%
  pivot_longer(
    cols = all_of(names(pathway_list)),
    names_to = "Score_col",
    values_to = "Score"
  ) 

p_score_violin <- ggplot(score_long, aes(x = .data[[group_col]], y = Score, fill = .data[[group_col]])) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, fill = "white") +
  facet_wrap(~ Score_col, ncol = 3, scales = "free_y") +
  # scale_fill_manual(values = group_colors) +
  theme_pub() +
  theme(legend.position = 'right') +
  labs(x = "", y = "Module Score", title = "Function Scores by CellType")

ggsave(
  file.path(outdir, "ModuleScore_violin_by_cell_type.pdf"),
  p_score_violin,
  width = 19,
  height = max(8, ceiling(length(pathway_list) / 3) * 3.5)
)

ggsave(
  file.path(outdir, "ModuleScore_violin_by_cell_type.png"),
  p_score_violin,
  width = 19,
  height = max(8, ceiling(length(pathway_list) / 3) * 3.5),
  dpi = 300
)

## =========================
## Module score UMAP 图
## 所有score一张图
## =========================
theme_umap_pub <- theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(
      size = 16, face = "bold", hjust = 0.5, color = "black",
      margin = margin(b = 6)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.text = element_text(size = 10)
  )

score_feature_plots <- vector("list", length(pathway_list))

for (i in seq_along(pathway_list)) {
  sc_col <- names(pathway_list)[i]
  module_nm <- names(pathway_list)[i]
  
  p <- FeaturePlot(
    object = obj,
    features = sc_col,
    pt.size = 1,
    reduction = 'umap',
    raster = FALSE,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95"
  ) +
    scale_color_gradientn(
      colors = c("#08306B", "white", "#CB181D")
    ) +
    ggtitle(module_nm) +
    theme_umap_pub
  
  score_feature_plots[[i]] <- p
}

p_score_umap <- wrap_plots(score_feature_plots, ncol = 3) +
  plot_annotation(title = "Module Score UMAP")

ggsave(
  file.path(outdir, "ModuleScore_UMAP_all.pdf"),
  p_score_umap,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 4)
)

ggsave(
  file.path(outdir, "ModuleScore_UMAP_all.png"),
  p_score_umap,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 4),
  dpi = 300
)

