library(Seurat)
# library(harmony)
library(ggplot2)
library(magrittr)
library(dplyr)
# library(clusterProfiler)
library(msigdbr)
library(UCell)
library(AUCell)
library(ggpubr)
library(rstatix)
library(ggsignif)

setwd("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/CD4_Tfh_scMetabolism/")
S_subf<-readRDS("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/Results/T/CD4_T/sub_CD4/CD4T_sub_harmony_umap_tsne_cell_type.rds")
# DimPlot(S_subf,reduction='umap',label = TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt",sep="\t",header = T)
lines <- readLines("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/script/KEGG_metabolism_nc.gmt")
geneset <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset[[set_name]] <- genes
  }
}
# geneset<-geneset$gene
# gene<-list(geneset)
S_subf<-subset(S_subf,idents='CD4_Tfh')
adata<-AddModuleScore_UCell(S_subf,features=geneset)
colnames(adata@meta.data)

# 假设基因集名称存储在向量中
# colnames(adata@meta.data) <- gsub(" / ", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub(" \\(", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("\\)", "", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("-", "_", colnames(adata@meta.data))
gene_set_names <- colnames(adata@meta.data)[19:103]  # 85个基因集的名称

# 创建输出目录
dir.create("./violin_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./umap_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./dot_plots", recursive = TRUE, showWarnings = FALSE)
# 循环处理每个基因集
adata$group <- factor(adata$group,levels = c("non_pcr_A","non_pcr_B","pcr_A","pcr_B"))
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")
for (i in seq_along(gene_set_names)) {
  score_col <- gene_set_names[i]
  if (!score_col %in% colnames(adata@meta.data)) {
    warning(paste("Column not found:", score_col, "- skipping"))
    next
  }
  cat("Processing", i, "/", length(gene_set_names), ":", score_col, "\n")
  # 绘制小提琴图
  p <- VlnPlot(adata, features = score_col, pt.size = 0, group.by = 'group') + NoLegend()+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.7),outlier.shape = NA)
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.pdf"), 
         width = 7, height =8, bg = 'white')
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.png"), 
         width = 7, height = 8, bg = 'white')
  p <- get_sig_ana(p,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
	get_sig_ana(p$data,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,outdir="./violin_plots/",map_sig_level=T,only_sig = F,method = "wilcox.test")
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.png"),width = 7,height = 10)
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.pdf"),width = 7,height = 10)
  # 绘制 UMAP 图
  p <- FeaturePlot(adata, features = score_col, reduction = 'umap', 
                   label = TRUE, max.cutoff = 0.99, cols = c("lightgrey","red"))
  ggsave(paste0("./umap_plots/", score_col, "_umap.png"), 
         width = 7, height = 6, bg = 'white')
  ggsave(paste0("./umap_plots/", score_col, "_umap.pdf"), 
         width = 7, height = 6, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'group', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.png"), 
         width = 7, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.pdf"), 
         width = 7, height = 4, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'cell_type', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.png"), 
         width = 12, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.pdf"), 
         width = 12, height = 4, bg = 'white')
}



setwd("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/CD8_Tm_C2_scMetabolism/")
S_subf<-readRDS("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/CD8T_sub_harmony_umap_tsne_cell_type.rds")
# DimPlot(S_subf,reduction='umap',label = TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt",sep="\t",header = T)
lines <- readLines("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/script/KEGG_metabolism_nc.gmt")
geneset <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset[[set_name]] <- genes
  }
}
# geneset<-geneset$gene
# gene<-list(geneset)
S_subf<-subset(S_subf,idents='CD8_Tm_C2')
adata<-AddModuleScore_UCell(S_subf,features=geneset)
colnames(adata@meta.data)

# 假设基因集名称存储在向量中
# colnames(adata@meta.data) <- gsub(" / ", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub(" \\(", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("\\)", "", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("-", "_", colnames(adata@meta.data))
gene_set_names <- colnames(adata@meta.data)[19:103]  # 85个基因集的名称

# 创建输出目录
dir.create("./violin_plots", recursive = TRUE, showWarnings = FALSE)
# dir.create("./umap_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./dot_plots", recursive = TRUE, showWarnings = FALSE)
# 循环处理每个基因集
adata$group <- factor(adata$group,levels = c("non_pcr_A","non_pcr_B","pcr_A","pcr_B"))
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")
for (i in seq_along(gene_set_names)) {
  score_col <- gene_set_names[i]
  if (!score_col %in% colnames(adata@meta.data)) {
    warning(paste("Column not found:", score_col, "- skipping"))
    next
  }
  cat("Processing", i, "/", length(gene_set_names), ":", score_col, "\n")
  # 绘制小提琴图
  p <- VlnPlot(adata, features = score_col, pt.size = 0, group.by = 'group') + NoLegend()+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.7),outlier.shape = NA)
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.pdf"), 
         width = 7, height =8, bg = 'white')
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.png"), 
         width = 7, height = 8, bg = 'white')
  p <- get_sig_ana(p,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
	get_sig_ana(p$data,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,outdir="./violin_plots/",map_sig_level=T,only_sig = F,method = "wilcox.test")
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.png"),width = 7,height = 10)
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.pdf"),width = 7,height = 10)
  # 绘制 UMAP 图
  # p <- FeaturePlot(adata, features = score_col, reduction = 'umap', 
  #                  label = TRUE, max.cutoff = 0.99, cols = c("lightgrey","red"))
  # ggsave(paste0("./umap_plots/", score_col, "_umap.png"), 
  #        width = 7, height = 6, bg = 'white')
  # ggsave(paste0("./umap_plots/", score_col, "_umap.pdf"), 
  #        width = 7, height = 6, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'group', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.png"), 
         width = 7, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.pdf"), 
         width = 7, height = 4, bg = 'white')
  # p <- DotPlot(adata, features = score_col, group.by = 'cell_type', cols = c("lightgrey","red")) + coord_flip() +
  #   theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
  #   ggtitle(score_col) 
  # ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.png"), 
  #        width = 12, height = 4, bg = 'white')
  # ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.pdf"), 
  #        width = 12, height = 4, bg = 'white')
}



setwd("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/KRT16_Tumor_scMetabolism/")
S_subf<-readRDS("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Epithelial/Epithelial.rds")
# DimPlot(S_subf,reduction='umap',label = TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2510235895/13.Mission/script/KEGG_metabolism_nc.gmt",sep="\t",header = T)
lines <- readLines("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/script/KEGG_metabolism_nc.gmt")
geneset <- list()
for (line in lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    geneset[[set_name]] <- genes
  }
}
# geneset<-geneset$gene
# gene<-list(geneset)
S_subf<-subset(S_subf,idents='KRT16+ Tumor')
adata<-AddModuleScore_UCell(S_subf,features=geneset)
colnames(adata@meta.data)

# 假设基因集名称存储在向量中
# colnames(adata@meta.data) <- gsub(" / ", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub(" \\(", "_", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("\\)", "", colnames(adata@meta.data))
# colnames(adata@meta.data) <- gsub("-", "_", colnames(adata@meta.data))
gene_set_names <- colnames(adata@meta.data)[15:99]  # 85个基因集的名称

# 创建输出目录
dir.create("./violin_plots", recursive = TRUE, showWarnings = FALSE)
# dir.create("./umap_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./dot_plots", recursive = TRUE, showWarnings = FALSE)
# 循环处理每个基因集
adata$group <- factor(adata$group,levels = c("non_pcr_A","non_pcr_B","pcr_A","pcr_B"))
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")
for (i in seq_along(gene_set_names)) {
  score_col <- gene_set_names[i]
  if (!score_col %in% colnames(adata@meta.data)) {
    warning(paste("Column not found:", score_col, "- skipping"))
    next
  }
  cat("Processing", i, "/", length(gene_set_names), ":", score_col, "\n")
  # 绘制小提琴图
  p <- VlnPlot(adata, features = score_col, pt.size = 0, group.by = 'group') + NoLegend()+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.7),outlier.shape = NA)
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.pdf"), 
         width = 7, height =8, bg = 'white')
  ggsave(paste0("./violin_plots/", score_col, ".VlnPlot.png"), 
         width = 7, height = 8, bg = 'white')
  p <- get_sig_ana(p,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)
	get_sig_ana(p$data,as.formula(paste0(score_col,'~ident')),group.by='ident',prefix=score_col,outdir="./violin_plots/",map_sig_level=T,only_sig = F,method = "wilcox.test")
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.png"),width = 7,height = 10)
	ggsave(paste0("./violin_plots/", score_col, ".VlnPlot_sig.pdf"),width = 7,height = 10)
  # 绘制 UMAP 图
  # p <- FeaturePlot(adata, features = score_col, reduction = 'umap', 
  #                  label = TRUE, max.cutoff = 0.99, cols = c("lightgrey","red"))
  # ggsave(paste0("./umap_plots/", score_col, "_umap.png"), 
  #        width = 7, height = 6, bg = 'white')
  # ggsave(paste0("./umap_plots/", score_col, "_umap.pdf"), 
  #        width = 7, height = 6, bg = 'white')
  p <- DotPlot(adata, features = score_col, group.by = 'group', cols = c("lightgrey","red")) + coord_flip() +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
    ggtitle(score_col) 
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.png"), 
         width = 7, height = 4, bg = 'white')
  ggsave(paste0("./dot_plots/", score_col, "_dotplot_group.pdf"), 
         width = 7, height = 4, bg = 'white')
  # p <- DotPlot(adata, features = score_col, group.by = 'cell_type', cols = c("lightgrey","red")) + coord_flip() +
  #   theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5,size = 10)) +
  #   ggtitle(score_col) 
  # ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.png"), 
  #        width = 15, height = 4, bg = 'white')
  # ggsave(paste0("./dot_plots/", score_col, "_dotplot_cell_type.pdf"), 
  #        width = 15, height = 4, bg = 'white')
}



setwd("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/")
S_subf<-readRDS("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/Results/T/CD4_T/sub_CD4/CD4T_sub_harmony_umap_tsne_cell_type.rds")
S_subf$group <- factor(S_subf$group,levels = c("non_pcr_A","non_pcr_B","pcr_A","pcr_B"))
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")

plot_data <- FetchData(S_subf, vars = c('CXCL13', 'group'))
plot_data <- plot_data %>%
  mutate(group = as.character(group))
head(plot_data)
summary_df <- plot_data %>%
  group_by(group) %>%
  summarise(
    mean = mean(CXCL13),
    sd = sd(CXCL13),
    n = n(),
    se = sd / sqrt(n)
  )
group_names <- unique(plot_data$group)
comparisons <- combn(group_names, 2, simplify = FALSE)
group_colors <- c("non_pcr_B" = "#1f77b4", "pcr_B" = "#ff7f0e", "non_pcr_A" = "#2ca02c", "pcr_A" = "#d62728")
stat_results <- plot_data %>%
  pairwise_t_test(CXCL13 ~ group, p.adjust.method = "none") %>%
  add_xy_position(x = "group")

p <- ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, size = 0.8, color = "black") +
    # 指定分组颜色
  scale_fill_manual(values = group_colors)+
  theme_classic(base_size = 14) +
  theme(
    # 外边框（加粗黑色边框）
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    legend.position = "none",
    # 坐标轴
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    # 标题
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
  ) +
  labs(
    x = "Group",
    y = "CXCL13 Expression",
    title = "CXCL13 Expression by Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave("CD4T_CXCL13_barplot.pdf", width = 7, height = 8, bg = 'white', dpi = 300)
ggsave("CD4T_CXCL13_barplot.png", width = 7, height = 8, bg = 'white', dpi = 300)


  # 使用预计算的统计结果
p <- ggplot() +
  geom_bar(data = summary_df, 
           aes(x = group, y = mean, fill = group), 
           stat = "identity", width = 0.7, color = "black", size = 0.5) +
  geom_errorbar(data = summary_df,
                aes(x = group, ymin = mean - se, ymax = mean + se), 
                width = 0.2, size = 0.8, color = "black") +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  # 使用 ggsignif 添加显著性
  geom_signif(
    data = plot_data,
    aes(x = group, y = CXCL13),
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.01,
		tip_length = 0.002,
    size = 0.8,
    textsize = 4,
    y_position = 0.2
  ) +
  labs(
    x = "Group",
    y = "CXCL13 Expression",
    title = "CXCL13 Expression by Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave("CD4T_CXCL13_barplot_sig.pdf", width = 7, height = 8, bg = 'white', dpi = 300)
ggsave("CD4T_CXCL13_barplot_sig.png", width = 7, height = 8, bg = 'white', dpi = 300)


S_subf<-readRDS("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/CD8T_sub_harmony_umap_tsne_cell_type.rds")
plot_data <- FetchData(S_subf, vars = c('CXCL13', 'group'))
plot_data <- plot_data %>%
  mutate(group = as.character(group))
head(plot_data)
summary_df <- plot_data %>%
  group_by(group) %>%
  summarise(
    mean = mean(CXCL13),
    sd = sd(CXCL13),
    n = n(),
    se = sd / sqrt(n)
  )
group_names <- unique(plot_data$group)
comparisons <- combn(group_names, 2, simplify = FALSE)
group_colors <- c("non_pcr_B" = "#1f77b4", "pcr_B" = "#ff7f0e", "non_pcr_A" = "#2ca02c", "pcr_A" = "#d62728")
stat_results <- plot_data %>%
  pairwise_t_test(CXCL13 ~ group, p.adjust.method = "none") %>%
  add_xy_position(x = "group")

p <- ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, size = 0.8, color = "black") +
    # 指定分组颜色
  scale_fill_manual(values = group_colors)+
  theme_classic(base_size = 14) +
  theme(
    # 外边框（加粗黑色边框）
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    legend.position = "none",
    # 坐标轴
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    # 标题
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
  ) +
  labs(
    x = "Group",
    y = "CXCL13 Expression",
    title = "CXCL13 Expression by Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave("CD8T_CXCL13_barplot.pdf", width = 7, height = 8, bg = 'white', dpi = 300)
ggsave("CD8T_CXCL13_barplot.png", width = 7, height = 8, bg = 'white', dpi = 300)


  # 使用预计算的统计结果
p <- ggplot() +
  geom_bar(data = summary_df, 
           aes(x = group, y = mean, fill = group), 
           stat = "identity", width = 0.7, color = "black", size = 0.5) +
  geom_errorbar(data = summary_df,
                aes(x = group, ymin = mean - se, ymax = mean + se), 
                width = 0.2, size = 0.8, color = "black") +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  # 使用 ggsignif 添加显著性
  geom_signif(
    data = plot_data,
    aes(x = group, y = CXCL13),
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.01,
		tip_length = 0.002,
    size = 0.8,
    textsize = 4,
    y_position = 0.7
  ) +
  labs(
    x = "Group",
    y = "CXCL13 Expression",
    title = "CXCL13 Expression by Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave("CD8T_CXCL13_barplot_sig.pdf", width = 7, height = 8, bg = 'white', dpi = 300)
ggsave("CD8T_CXCL13_barplot_sig.png", width = 7, height = 8, bg = 'white', dpi = 300)