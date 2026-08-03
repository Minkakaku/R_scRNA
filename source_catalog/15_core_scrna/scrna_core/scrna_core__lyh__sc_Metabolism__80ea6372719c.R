#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: run_scMetabolism.R
# Description: Single-cell metabolism pathway scoring (scMetabolism) and 
#              visualization (FeaturePlot, DotPlot, ViolinPlot)
# ==============================================================================

library(optparse)

# 1. Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to the input Seurat QS object (Required)", metavar = "character"),
  make_option(c("-w", "--work_dir"), type = "character", default = "./", 
              help = "Output working directory path [Default: %default]", metavar = "character"),
  make_option(c("-c", "--celltype_col"), type = "character", default = "Celltype", 
              help = "Column name in Metadata representing cell types [Default: %default]", metavar = "character"),
  make_option(c("-g", "--group_col"), type = "character", default = "Group", 
              help = "Column name in Metadata representing groups (e.g., Disease/Control) [Default: %default]", metavar = "character"),
  make_option(c("-s", "--score_method"), type = "character", default = "AUCell", 
              help = "Scoring method (AUCell, VISION, ssGSEA) [Default: %default]", metavar = "character"),
  make_option(c("-m", "--metabolism_type"), type = "character", default = "KEGG", 
              help = "Metabolism database (KEGG or REACTOME) [Default: %default]", metavar = "character"),
  make_option(c("-n", "--core_num"), type = "integer", default = 2, 
              help = "Number of threads used for sc.metabolism scoring [Default: %default]", metavar = "integer"),
  make_option(c("-t", "--plot_workers"), type = "integer", default = 8, 
              help = "Number of parallel workers for future.apply plotting [Default: %default]", metavar = "integer"),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Argument Validation
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Error: Input .qs file path (-i/--input) must be provided.\n", call. = FALSE)
}

if (!file.exists(opt$color_script)) {
  warning(sprintf("Warning: Color script %s not found. Please check this file if plotting fails later!", opt$color_script))
}

# 2. Load dependencies and environment
message(">>> Loading R packages...")
suppressMessages(require(pacman))
p_load(Seurat, tidyverse, qs, scMetabolism, rsvd, scRNAtoolVis, scop, future.apply)

# Load color script (Assuming it defines color3, color4, color_s2)
source("/home/lyh/Script/R_fuc/color.R")

# Set working directory
if (!dir.exists(opt$work_dir)) {
  dir.create(opt$work_dir, recursive = TRUE)
}
setwd(opt$work_dir)
message(">>> Working directory set to: ", getwd())

# 3. Read data and perform metabolism scoring
message(">>> Reading Seurat object: ", opt$input)
obj <- qread(opt$input)

message(sprintf(">>> Starting metabolism scoring (Method: %s, Database: %s, Threads: %s)...", 
                opt$score_method, opt$metabolism_type, opt$core_num))
exp_matrix <- LayerData(obj, layer = "counts")
df <- sc.metabolism(countexp = exp_matrix, 
                    method = opt$score_method, 
                    imputation = FALSE, 
                    ncores = opt$core_num, 
                    metabolism.type = opt$metabolism_type) %>% t()

# Save results and add to metadata
colnames(df) <- gsub("/", "_", colnames(df))
write.csv(df, "scMetabolism.csv", row.names = TRUE)

obj <- AddMetaData(obj, as.data.frame(df))
colnames(obj@meta.data) <- gsub("/", "_", colnames(obj@meta.data))
message(">>> Metabolism scoring completed, results saved to scMetabolism.csv")

# 4. Visualization - FeaturePlot
message(">>> Starting to plot FeaturePlots...")
plan(multisession, workers = opt$plot_workers)
dir.create("Metabolism_Featureplot", showWarnings = FALSE)

plot_func_feature <- function(i, obj) {
  p <- scop::FeatureDimPlot(obj, features = i, reduction = "umap", theme_use = "theme_blank")
  ggsave(paste0("Metabolism_Featureplot/", i, "_dim.pdf"), p, height = 6, width = 7)
  ggsave(paste0("Metabolism_Featureplot/", i, "_dim.png"), p, height = 6, width = 7, dpi = 300)
}
future_lapply(colnames(df), plot_func_feature, obj = obj, future.seed = TRUE)
plan(sequential)

# 5. Visualization - DotPlot
message(">>> Starting to plot DotPlots...")
obj$Celltype_group <- paste0(obj@meta.data[, opt$group_col], "_", obj@meta.data[, opt$celltype_col])

# Set Level order
groups <- unique(obj@meta.data[, opt$group_col]) %>% as.character()
celltypes <- unique(obj@meta.data[, opt$celltype_col]) %>% as.character()
ordered_levels <- as.vector(outer(groups, celltypes, function(x, y) paste0(x, "_", y)))
obj$Celltype_group <- factor(obj$Celltype_group, levels = ordered_levels)

# Set group colors
group_colors <- color3[1:length(groups)]
x_axis_colors <- rep(group_colors, times = length(celltypes))  

# Custom theme
my_theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, color = x_axis_colors, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"))

# Plot
width_num <- length(unique(obj$Celltype_group))
height_num <- length(colnames(df))
options(repr.plot.width = width_num * 1.5, repr.plot.height = height_num * 0.4) 

p_dot <- Seurat::DotPlot(obj, features = colnames(df), group.by = "Celltype_group") + 
  scale_color_gradientn(name = "Score", colours = color_s2) + 
  coord_flip() + labs(x = "", y = "") + my_theme 

ggsave("Metabolism_dotplot.pdf", p_dot, height = height_num * 0.4, width = width_num * 1.5)
ggsave("Metabolism_dotplot.png", p_dot, height = height_num * 0.4, width = width_num * 1.5, dpi = 300)

# 6. Visualization - ViolinPlot
message(">>> Starting to plot ViolinPlots...")
plan(multisession, workers = opt$plot_workers)
dir.create("Metabolism_vlnplot", showWarnings = FALSE)

plot_func_vln <- function(i, obj) {
  # Assuming plotthis package is installed
  p <- plotthis::ViolinPlot(obj@meta.data, x = opt$celltype_col, y = i, group_by = opt$group_col, 
                            comparisons = TRUE, sig_label = "p.signif", add_bg = TRUE, 
                            bg_palcolor = color3, palcolor = color4, 
                            legend.position = "top", legend.direction = "horizontal", 
                            add_box = TRUE, box_ptsize = 1)
  ggsave(paste0("Metabolism_vlnplot/", i, "_vln.pdf"), p, height = 5, width = 8)
  ggsave(paste0("Metabolism_vlnplot/", i, "_vln.png"), p, height = 5, width = 8, dpi = 300)
}
future_lapply(colnames(df), plot_func_vln, obj = obj, future.seed = TRUE)
plan(sequential)

message(">>> All tasks completed successfully!")