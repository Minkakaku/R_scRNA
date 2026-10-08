# 从人工选定 PC 的检查点继续整合，保存完整对象；绘图可单独重跑。
step_harmony_selected_pcs <- function() {
  suppressPackageStartupMessages(library(Seurat))
  future::plan(future::sequential)
  settings <- CONFIG$harmony_reintegration
  dims <- as.integer(unlist(settings$selected_pcs))
  stopifnot(length(dims) > 1L, !anyDuplicated(dims), all(dims > 0))
  destination <- file.path(PROJECT_DIR, settings$integrated_output)
  output <- file.path(PROJECT_DIR, settings$integrated_results)
  stopifnot(!file.exists(destination), !file.exists(paste0(destination, ".partial")))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  sc <- run_step("Read PCA checkpoint", qs::qread(file.path(PROJECT_DIR, settings$pca_checkpoint), nthreads = 4))
  stopifnot(all(dims <= ncol(Embeddings(sc, "pca"))),
            !"harmony" %in% Reductions(sc), !"harmony_cluster" %in% names(sc[[]]))
  original_metadata <- sc[[]]
  original_cells <- colnames(sc)
  original_genes <- rownames(sc)
  # 完整稀疏矩阵指纹用于核对 counts 和标准化表达，避免只比列和漏掉变化。
  expression_hash <- function(x) vapply(c("counts", "data"), function(layer)
    digest::digest(LayerData(x, assay = "RNA", layer = layer), algo = "xxhash64"), character(1))
  original_expression <- expression_hash(sc)
  original_reductions <- vapply(sc@reductions, digest::digest, character(1), algo = "xxhash64")
  batches <- as.character(sc[[settings$batch_col, drop = TRUE]])
  stopifnot(!anyNA(batches), all(nzchar(batches)), length(unique(batches)) > 1L)
  DefaultAssay(sc) <- "RNA"
  set.seed(SEED)
  # Harmony 2.x 的接口；仅校正低维坐标，不修改表达矩阵。
  sc <- run_step("Harmony on selected PCs", harmony::RunHarmony(
    sc, group.by.vars = settings$batch_col, reduction.use = "pca", dims.use = dims,
    reduction.save = "harmony", project.dim = FALSE, theta = settings$theta,
    lambda = settings$lambda, max_iter = settings$max_iter, ncores = settings$ncores,
    early_stop = TRUE, verbose = TRUE
  ))
  harmony_dims <- seq_along(dims)
  sc <- run_step("FindNeighbors on Harmony", FindNeighbors(
    sc, reduction = "harmony", dims = harmony_dims, k.param = 20,
    graph.name = c("harmony_nn", "harmony_snn"), verbose = TRUE
  ))
  sc <- run_step("FindClusters on Harmony", FindClusters(
    sc, graph.name = "harmony_snn", cluster.name = "harmony_cluster",
    resolution = settings$resolution, algorithm = settings$algorithm,
    random.seed = SEED, verbose = TRUE
  ))
  # Seurat 5.5 自定义 cluster.name 时不会自动写 seurat_clusters；显式同步。
  sc$seurat_clusters <- sc$harmony_cluster
  # 校正前后使用相同 PC 范围与 UMAP 参数，各自重新计算二维坐标。
  for (reduction in c("pca", "harmony")) {
    reduction_name <- if (reduction == "pca") "umap.pca" else "umap"
    sc <- run_step(paste("RunUMAP", reduction), RunUMAP(
      sc, reduction = reduction, dims = if (reduction == "pca") dims else harmony_dims,
      reduction.name = reduction_name,
      reduction.key = if (reduction == "pca") "UMAPPCA_" else "UMAP_",
      n.neighbors = settings$umap_neighbors, min.dist = settings$umap_min_dist,
      metric = settings$umap_metric, seed.use = SEED, verbose = TRUE
    ))
  }
  sc@misc$harmony_reintegration$stage <- "integrated"
  sc@misc$harmony_reintegration$selected_pcs <- dims
  sc@misc$harmony_reintegration$settings <- settings
  sc@misc$harmony_reintegration$cluster_column <- "harmony_cluster"
  sc@misc$harmony_reintegration$input_expression_hash <- original_expression
  stopifnot(identical(sc[[]][, names(original_metadata), drop = FALSE], original_metadata),
            identical(expression_hash(sc), original_expression))
  for (name in names(original_reductions)) {
    stopifnot(identical(digest::digest(sc[[name]], algo = "xxhash64"), original_reductions[[name]]))
  }
  for (name in c("harmony", "umap.pca", "umap")) {
    stopifnot(identical(rownames(Embeddings(sc, name)), original_cells),
              all(is.finite(Embeddings(sc, name))))
  }
  stopifnot(ncol(Embeddings(sc, "harmony")) == length(dims),
            !anyNA(sc$harmony_cluster), !anyNA(sc$seurat_clusters))
  expected_reductions <- vapply(sc@reductions, digest::digest, character(1), algo = "xxhash64")
  expected_clusters <- sc$harmony_cluster
  run_step("Save complete integrated QS", qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4))
  rm(sc)
  invisible(gc())
  check <- run_step("Read back complete integrated QS", qs::qread(paste0(destination, ".partial"), nthreads = 4))
  stopifnot(isTRUE(validObject(check, test = TRUE)), identical(colnames(check), original_cells),
            identical(rownames(check), original_genes), identical(expression_hash(check), original_expression),
            identical(check[[]][, names(original_metadata), drop = FALSE], original_metadata),
            identical(check$harmony_cluster, expected_clusters),
            identical(vapply(check@reductions, digest::digest, character(1), algo = "xxhash64"), expected_reductions),
            identical(as.character(Idents(check)), as.character(expected_clusters)),
            file.rename(paste0(destination, ".partial"), destination))
  parameters <- list(stage = "integrated", cells = ncol(check), genes = nrow(check),
    samples = length(unique(batches)), clusters = nlevels(droplevels(check$harmony_cluster)),
    selected_pcs = dims, seed = SEED, harmony_version = as.character(packageVersion("harmony")),
    settings = settings, expression_hash = as.list(original_expression))
  yaml::write_yaml(parameters, file.path(output, "integration_parameters.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "06_harmony_selected_pcs_sessionInfo.txt"))
  log_msg("HARMONY_INTEGRATION_OK: ", ncol(check), " cells; ", parameters$clusters,
          " clusters; selected PCs=", paste(dims, collapse = ","))
  invisible(destination)
}

step_plot_selected_harmony <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$harmony_reintegration
  output <- file.path(PROJECT_DIR, settings$integrated_results)
  sc <- run_step("Read integrated QS for plotting", qs::qread(file.path(PROJECT_DIR, settings$integrated_output), nthreads = 4))
  metadata <- sc[[]]
  umap <- Embeddings(sc, "umap")
  pca_umap <- Embeddings(sc, "umap.pca")
  frame <- data.frame(cell_id = colnames(sc), metadata, UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
                      PCA_UMAP_1 = pca_umap[, 1], PCA_UMAP_2 = pca_umap[, 2], check.names = FALSE)
  write.table(frame, gzfile(file.path(output, "cell_metadata_umap.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
  for (column in c("sample", "group", "cell_type")) {
    tab <- as.data.frame(table(cluster = metadata$harmony_cluster, label = metadata[[column]]))
    names(tab) <- c("harmony_cluster", column, "cells")
    write.table(tab, file.path(output, paste0("cluster_by_", column, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  counts <- as.data.frame(table(cluster = metadata$harmony_cluster))
  names(counts) <- c("harmony_cluster", "cells")
  write.table(counts, file.path(output, "cluster_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  counts <- as.data.frame(table(sample = metadata$sample))
  names(counts)[2] <- "cells"
  write.table(counts, file.path(output, "sample_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # 复用来源对象概览的细胞类型颜色规则，保证旧标签在本轮图中颜色一致。
  old_plot_data <- readRDS(file.path(PROJECT_DIR, CONFIG$object_overview$plot_data))
  common_cells <- Reduce(intersect, lapply(old_plot_data, function(x) rownames(x$metadata)))
  shared_labels <- sort(unique(as.character(old_plot_data[[1]]$metadata[common_cells, "cell_type"])))
  type_levels <- c(shared_labels, setdiff(sort(unique(as.character(metadata$cell_type))), shared_labels))
  rgb <- grDevices::col2rgb(HONGF_NPG) / 255
  extended <- c(HONGF_NPG, grDevices::rgb(t(rgb * .65 + .35)),
                grDevices::rgb(t(rgb * .70)), grDevices::rgb(t(rgb * .50 + .50)))
  type_colors <- setNames(extended[seq_along(type_levels)], type_levels)
  sc$cell_type <- factor(as.character(sc$cell_type), levels = type_levels)
  rm(old_plot_data)
  invisible(gc())
  make_plot <- function(column, title, reduction = "umap", ...) {
    p <- if (column == "cell_type") {
      DimPlot(sc, reduction = reduction, group.by = column, cols = type_colors,
              raster = FALSE, shuffle = TRUE, seed = SEED, ...)
    } else {
      npg_dimplot(sc, reduction = reduction, group.by = column,
                  raster = FALSE, shuffle = TRUE, seed = SEED, ...)
    }
    p + labs(title = title) + theme_classic(base_size = 11) +
      theme(plot.title = element_text(hjust = .5), legend.text = element_text(size = 9))
  }
  save_one <- function(p, stem, width, height) run_step(paste("Save figure", stem),
    save_umap_plot(p, file.path(output, stem), width = width, height = height))
  cluster <- make_plot("harmony_cluster", "Harmony clusters", label = TRUE, repel = TRUE) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))
  sample <- make_plot("sample", "Harmony by sample")
  group <- make_plot("group", "Harmony by group")
  celltype <- make_plot("cell_type", "Original cell types on Harmony UMAP") +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))
  save_one(cluster, "UMAP_harmony_clusters", 10, 7)
  save_one(sample, "UMAP_harmony_sample", 8, 7)
  save_one(group, "UMAP_harmony_group", 8, 7)
  save_one(celltype, "UMAP_harmony_original_cell_type", 12, 8)
  before <- make_plot("sample", "Before Harmony", reduction = "umap.pca")
  after <- make_plot("sample", "After Harmony")
  save_one(patchwork::wrap_plots(before, after, ncol = 2, guides = "collect"), "UMAP_sample_before_after", 14, 6)
  split_sample <- make_plot("harmony_cluster", "Harmony clusters by sample", split.by = "sample", ncol = 3) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))
  save_one(split_sample, "UMAP_harmony_sample_split", 16, 10)
  split_group <- make_plot("harmony_cluster", "Harmony clusters by group", split.by = "group", ncol = 2) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))
  save_one(split_group, "UMAP_harmony_group_split", 14, 6)
  log_msg("HARMONY_PLOTS_OK: 7 vector PDF / 600 dpi PNG pairs; ", nrow(frame), " cells")
}
