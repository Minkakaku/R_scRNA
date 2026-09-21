# PCA graph clustering and UMAP

# 由各 PC 的标准差计算方差、解释比例和累计解释比例。
pca_variance_table <- function(obj, reduction = "pca") {
  pc_sd <- Stdev(object = obj, reduction = reduction)
  variance <- pc_sd^2
  variance_percent <- variance / sum(variance) * 100
  data.frame(
    pc = seq_along(variance),
    variance = variance,
    variance_percent = variance_percent,
    cumulative_percent = cumsum(variance_percent)
  )
}

# 并列记录 PCAtools 与 nFactors 的 PC 建议，以 PCAtools 结果作为推荐值。
detect_pca_elbow <- function(obj, reduction = "pca") {
  variance <- pca_variance_table(obj, reduction)
  scree <- nFactors::nScree(x = variance$variance)
  selections <- data.frame(
    method = c("PCAtools", "nFactors_noc", "nFactors_naf"),
    pc = as.integer(c(
      PCAtools::findElbowPoint(variance$variance_percent),
      scree$Components$noc,
      scree$Components$naf
    ))
  )
  list(
    recommended_pc = selections$pc[selections$method == "PCAtools"],
    selections = selections,
    variance = variance
  )
}

# 保存 PCA 方差表、不同方法的拐点建议及带标线的配对图。
save_pca_elbow <- function(obj, path_no_ext, reduction = "pca") {
  result <- detect_pca_elbow(obj, reduction)
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  write.csv(result$variance, paste0(path_no_ext, "_variance.csv"), row.names = FALSE, quote = FALSE)
  write.csv(result$selections, paste0(path_no_ext, "_selection.csv"), row.names = FALSE, quote = FALSE)

  method_colors <- c(PCAtools = HONGF_NPG[1], nFactors_noc = HONGF_NPG[4], nFactors_naf = HONGF_NPG[3])
  p <- ggplot(result$variance, aes(pc, variance_percent)) +
    geom_line(linewidth = 0.55, color = "#4D4D4D") +
    geom_point(size = 1.6, color = "#4D4D4D") +
    geom_vline(
      data = result$selections,
      aes(xintercept = pc, color = method, linetype = method),
      linewidth = 0.7
    ) +
    scale_color_manual(values = method_colors) +
    labs(
      x = "Principal component",
      y = "Explained variance (%)",
      title = "PCA variance and automatic elbow",
      color = NULL,
      linetype = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))
  save_plot(p, path_no_ext)
  log_msg(
    "Automatic PCA elbow: ",
    paste0(result$selections$method, "=PC", result$selections$pc, collapse = ", ")
  )
  invisible(result)
}

# 在对象 misc 中记录本轮聚类实际使用的降维与 PC 维度。
store_clustering_metadata <- function(obj, context, reduction, dims) {
  if (is.null(obj@misc$pipeline)) {
    obj@misc$pipeline <- list()
  }
  obj@misc$pipeline[[paste0(context, "_reduction")]] <- reduction
  obj@misc$pipeline[[paste0(context, "_dims")]] <- dims
  obj
}

# 优先复用已记录的降维；缺失时按整合 CCA、Harmony、PCA 的顺序选择。
choose_clustering_reduction <- function(obj, context, fallback = c("integrated.cca", "harmony", "pca")) {
  stored <- obj@misc$pipeline[[paste0(context, "_reduction")]]
  if (!is.null(stored) && stored %in% names(obj@reductions)) {
    return(stored)
  }
  available <- fallback[fallback %in% names(obj@reductions)]
  if (length(available) > 0) {
    return(available[1])
  }
  stop("No usable dimensional reduction found for clustering.", call. = FALSE)
}

# 优先复用已记录的维度，否则使用传入维度，并去除非正值。
choose_clustering_dims <- function(obj, context, requested_dims) {
  stored <- obj@misc$pipeline[[paste0(context, "_dims")]]
  dims <- if (!is.null(stored))
    stored else requested_dims
  dims[dims > 0]
}

# 在指定降维和维度上构建邻接图、聚类及 UMAP，可分别关闭重算步骤。
run_graph_cluster_umap <- function(obj, reduction, dims, resolution, context_label, run_graph_cluster = TRUE,
  run_umap = TRUE) {
  if (run_graph_cluster) {
    obj <- run_step(paste0(context_label, " FindNeighbors"), FindNeighbors(obj, reduction = reduction,
      dims = dims, verbose = FALSE))
    obj <- run_step(paste0(context_label, " FindClusters"), FindClusters(obj, resolution = resolution,
      random.seed = SEED, verbose = FALSE))
  }
  if (run_umap) {
    obj <- run_step(paste0(context_label, " RunUMAP"), RunUMAP(obj, reduction = reduction, dims = dims))
  }
  obj
}

# 从未校正的 PCA 重新生成独立 UMAP，供批次校正前后比较。
run_before_batch_umap <- function(obj, dims) {
  if (BATCH_EFFECT_UMAP$before_reduction %in% Reductions(obj)) {
    obj[[BATCH_EFFECT_UMAP$before_reduction]] <- NULL
  }
  run_step("Global RunUMAP before batch correction", RunUMAP(
    obj,
    reduction = "pca",
    dims = dims,
    reduction.name = BATCH_EFFECT_UMAP$before_reduction,
    reduction.key = BATCH_EFFECT_UMAP$before_key,
    seed.use = SEED
  ))
}

# 执行 PCA、拐点评估、批次校正、图聚类与 UMAP，并记录实际参数。
run_global_clustering <- function(sc) {
  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  sc <- run_step("RunPCA", RunPCA(sc, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(sc, file.path(GLOBAL_DIR, "PCA_elbow"))
  sc <- run_before_batch_umap(sc, pca_settings$dims)
  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Global")
  sc <- corrected$object
  reduction <- corrected$reduction
  sc <- run_graph_cluster_umap(sc, reduction = reduction, dims = pca_settings$dims, resolution = CLUSTER$resolution,
    context_label = "Global")
  sc <- store_clustering_metadata(sc, "global", reduction, pca_settings$dims)
  sc
}
