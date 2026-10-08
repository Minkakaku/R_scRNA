# 用 Seurat 自带小数据检查自定义 cluster.name、完整 QS 保存和读回兼容性。
suppressPackageStartupMessages(library(Seurat))
source("src/R/00_runtime.R")
source("src/R/steps/16_harmony_selected_pcs.R")
data("pbmc_small", package = "SeuratObject")
pbmc_small$sample <- rep(c("sample_1", "sample_2"), length.out = ncol(pbmc_small))
if ("seurat_clusters" %in% names(pbmc_small[[]])) pbmc_small$seurat_clusters <- NULL
PROJECT_DIR <- tempfile("harmony_compat_")
dir.create(PROJECT_DIR)
LOGS_DIR <- PROJECT_DIR
SEED <- 1001L
qs::qsave(pbmc_small, file.path(PROJECT_DIR, "input.qs"))
CONFIG <- list(harmony_reintegration = list(
  pca_checkpoint = "input.qs", integrated_output = "output.qs", integrated_results = "results",
  selected_pcs = 1:5, batch_col = "sample", theta = 2, lambda = 1,
  max_iter = 5, ncores = 1, resolution = .4, algorithm = 1,
  umap_neighbors = 15, umap_min_dist = .3, umap_metric = "cosine"
))
step_harmony_selected_pcs()
sc <- qs::qread(file.path(PROJECT_DIR, "output.qs"))
stopifnot(ncol(sc) == ncol(pbmc_small), identical(sc$harmony_cluster, sc$seurat_clusters),
          identical(as.character(Idents(sc)), as.character(sc$harmony_cluster)))
cat("HARMONY_SEURAT55_QS_COMPAT_OK\n")
