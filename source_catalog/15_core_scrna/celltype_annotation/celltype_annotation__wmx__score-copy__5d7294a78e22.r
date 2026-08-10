# 对于genelist里的每个list进行单细胞addmodulescore
# 对于每个list的打分进行对比小提琴图绘制，进行差异对比
# 对于打分绘制umap图，一张图展示

rm(list = ls())
rds_file <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou7/Microglia_fixed.rds"
geneset_file <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou7/geneset.csv"
outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou7/score/"
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
obj$Cmss1_exp = obj@assays$RNA@data['Cmss1',]
obj$Cmss1_group = ifelse(obj$Cmss1_exp > median(obj$Cmss1_exp), 'Cmss1_high', 'Cmss1_low')
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
safe_wilcox <- function(df, value_col = "value", group_col = "group") {
  df <- df %>% filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  g <- unique(as.character(df[[group_col]]))
  x <- df %>% filter(.data[[group_col]] == g[1]) %>% pull(.data[[value_col]])
  y <- df %>% filter(.data[[group_col]] == g[2]) %>% pull(.data[[value_col]])
  if (length(x) < 1 || length(y) < 1) {
    return(data.frame(p = NA_real_))
  }
  p <- tryCatch(
    wilcox.test(x, y)$p.value,
    error = function(e) NA_real_
  )
  data.frame(p = p)
}
## =========================
## Module score 小提琴图 + 差异统计
## =========================
group_col <- 'Cmss1_group'

score_long <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_col), all_of(names(pathway_list))) %>%
  pivot_longer(
    cols = all_of(names(pathway_list)),
    names_to = "Score_col",
    values_to = "Score"
  ) 
score_long$Cmss1_group <- factor(score_long$Cmss1_group,levels = c("Cmss1_low","Cmss1_high"))
group_colors <- c("Cmss1_low" = "#2B739F", "Cmss1_high" = "#E68525")
score_stats <- score_long %>%
  group_by(Score_col) %>%
  group_modify(~ safe_wilcox(.x, value_col = "Score", group_col = group_col)) %>%
  ungroup() %>%
  mutate(label = sapply(p, p_to_label))

score_ypos <- score_long %>%
  group_by(Score_col) %>%
  summarise(
    y_max = max(Score, na.rm = TRUE),
    y_min = min(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_pos = y_max + (y_max - y_min + 1e-6) * 0.08)

sig_labels_score <- left_join(score_stats, score_ypos, by = "Score_col")

p_score_violin <- ggplot(score_long, aes(x = .data[[group_col]], y = Score, fill = .data[[group_col]])) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, fill = "white") +
  geom_text(
    data = sig_labels_score,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  facet_wrap(~ Score_col, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_pub() +
  labs(x = "", y = "Module Score", title = "Function Scores")

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".pdf")),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5)
)

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".png")),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5),
  dpi = 300
)

# 分组绘图
group_col <- 'group'
score_long <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_col), all_of(names(pathway_list))) %>%
  pivot_longer(
    cols = all_of(names(pathway_list)),
    names_to = "Score_col",
    values_to = "Score"
  ) 
score_long$group<- factor(score_long$group,levels = c("M_3d","M_15d","M_30d"))
group_colors <- c(
  "M_3d" = "#66C2A5",
  "M_15d" = "#FC8D62",
  "M_30d"  = "#8DA0CB"
)
p_score_violin <- ggplot(score_long, aes(x = .data[[group_col]], y = Score, fill = .data[[group_col]])) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, fill = "white") +
  facet_wrap(~ Score_col, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_pub() +
  theme(legend.position = 'right',
        plot.margin = margin(t = 10, r = 10, b = 20, l = 40)) +
  labs(x = "", y = "Module Score", title = "Function Scores")

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".pdf")),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5)
)

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".png")),
  p_score_violin,
  width = 14,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5),
  dpi = 300
)

# celltype
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
  theme(legend.position = 'right',
        plot.margin = margin(t = 10, r = 10, b = 20, l = 40)) +
  labs(x = "", y = "Module Score", title = "Function Scores")

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".pdf")),
  p_score_violin,
  width = 18,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5)
)

ggsave(
  file.path(outdir, paste0("ModuleScore_violin_by_",group_col,".png")),
  p_score_violin,
  width = 18,
  height = max(8, ceiling(length(pathway_list) / 3) * 6.5),
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
  width = 18,
  height = max(8, ceiling(length(pathway_list) / 3) * 4)
)

ggsave(
  file.path(outdir, "ModuleScore_UMAP_all.png"),
  p_score_umap,
  width = 18,
  height = max(8, ceiling(length(pathway_list) / 3) * 4),
  dpi = 300
)

# 保存分组的UMAP
p <- FeaturePlot(
    object = obj,
    features = "Cmss1_exp",
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
    ggtitle('Cmss1 exp') +
    theme_umap_pub


ggsave(
  file.path(outdir, "Cmss1_exp_UMAP.pdf"),
  p,
  width =6,
  height = 6
)

ggsave(
  file.path(outdir, "Cmss1_exp_UMAP.png"),
  p,
  width = 6,
  height =6,
  dpi = 300
)
group_colors <- c("Cmss1_low" = "#2B739F", "Cmss1_high" = "#E68525")
p <- DimPlot(obj,group.by = 'Cmss1_group') +
  scale_color_manual(values = group_colors) +
  ggtitle('Cmss1 group') +
  theme_umap_pub


ggsave(
  file.path(outdir, "Cmss1_group_UMAP.pdf"),
  p,
  width = 7.5,
  height = 6
)

ggsave(
  file.path(outdir, "Cmss1_group_UMAP.png"),
  p,
  width = 7.5,
  height = 6,
  dpi = 300
)

saveRDS(object = obj,file = "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou7/score/score.rds")


## =========================
## 每个 list 的基因表达箱线图
## 一个 list 一张图，基因分面
## =========================
gene_boxplot_dir <- file.path(outdir, "Gene_boxplot_by_module")
dir.create(gene_boxplot_dir, showWarnings = FALSE, recursive = TRUE)

gene_stats_all <- list()

group_col = "Cmss1_group"
for (i in seq_along(names(pathway_list))) {
  module_nm <- names(pathway_list)[i]
  genes_use <- pathway_list[[module_nm]] %>% intersect(.,rownames(obj))
  
  if (length(genes_use) == 0) next
  
  expr_df <- FetchData(
    object = obj,
    vars = c(group_col, genes_use)
  ) %>%
    rownames_to_column("cell") %>%
    pivot_longer(
      cols = all_of(genes_use),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    dplyr::rename(group = !!group_col)
  
  # 统计
  gene_stats <- expr_df %>%
    group_by(Gene) %>%
    group_modify(~ safe_wilcox(.x, value_col = "Expression", group_col = "group")) %>%
    ungroup() %>%
    mutate(
      Module = module_nm,
      label = sapply(p, p_to_label)
    )
  
  ypos_df <- expr_df %>%
    group_by(Gene) %>%
    summarise(
      y_max = max(Expression, na.rm = TRUE),
      y_min = min(Expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_pos = y_max + (y_max - y_min + 1e-6) * 0.08)
  
  sig_labels_gene <- left_join(gene_stats, ypos_df, by = "Gene")
  
  gene_stats_all[[module_nm]] <- gene_stats
  
  expr_df$group <- factor(expr_df$group ,levels = c("Cmss1_low","Cmss1_high"))
  group_colors <- c("Cmss1_low" = "#2B739F", "Cmss1_high" = "#E68525")
  p_gene_box <- ggplot(expr_df, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
    # geom_jitter(width = 0.18, size = 0.2, alpha = 0.2, color = "grey30") +
    geom_text(
      data = sig_labels_gene,
      aes(x = 1.5, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    facet_wrap(~ Gene, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    theme_pub() +
    labs(
      x = "",
      y = "Expression",
      title = paste0(module_nm, " genes by group")
    )
  
  safe_nm <- make.names(module_nm)
  
  ggsave(
    file.path(gene_boxplot_dir, paste0(safe_nm, "_gene_boxplot.pdf")),
    p_gene_box,
    width = 16,
    height = max(6, ceiling(length(genes_use) / 4) * 3.2)
  )
  
  ggsave(
    file.path(gene_boxplot_dir, paste0(safe_nm, "_gene_boxplot.png")),
    p_gene_box,
    width = 16,
    height = max(6, ceiling(length(genes_use) / 4) * 3.2),
    dpi = 300
  )
}

gene_stats_all_df <- bind_rows(gene_stats_all)
write.csv(gene_stats_all_df, file.path(outdir, "gene_expression_wilcox_stats.csv"), row.names = FALSE)

## =========================
## 每个list基因的UMAP表达图
## 一个list一张图（优化主题版）
## =========================
gene_umap_dir <- file.path(outdir, "Gene_UMAP_by_module")
dir.create(gene_umap_dir, showWarnings = FALSE, recursive = TRUE)

for (module_nm in names(pathway_list) ){
  genes_use <- pathway_list[[module_nm]]
  if (length(genes_use) == 0) next
  
  gene_plots <- FeaturePlot(
    object = obj,
    pt.size = 1,
    features = genes_use,
    reduction = 'umap',
    raster = FALSE,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    combine = FALSE
  )
  
  gene_plots <- lapply(seq_along(gene_plots), function(i) {
    gene_plots[[i]] +
      scale_color_gradientn(
        colors = c("#08306B", "white", "#CB181D")
      )  +
      ggtitle(genes_use[i]) +
      theme_umap_pub
  })
  
  p_gene_umap <- wrap_plots(gene_plots, ncol = 4) +
    plot_annotation(
      title = paste0(module_nm, " gene expression UMAP"),
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )
  
  safe_nm <- make.names(module_nm)
  
  ggsave(
    file.path(gene_umap_dir, paste0(safe_nm, "_gene_UMAP.pdf")),
    p_gene_umap,
    width = 18,
    height = max(6, ceiling(length(genes_use) / 4) * 4)
  )
  
  ggsave(
    file.path(gene_umap_dir, paste0(safe_nm, "_gene_UMAP.png")),
    p_gene_umap,
    width = 18,
    height = max(6, ceiling(length(genes_use) / 4) * 4),
    dpi = 300
  )
}

obj@meta.data %>%  as.data.frame() %>% rownames_to_column("barcode") %>% write.csv("/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou7/score/Score_meta.csv",quote = F,row.names = F)

