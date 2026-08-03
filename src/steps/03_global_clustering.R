save_batch_effect_umap_outputs <- function(sc) {
  before <- BATCH_EFFECT_UMAP$before_reduction
  after <- BATCH_EFFECT_UMAP$after_reduction
  sample_col <- BATCH_EFFECT_UMAP$sample_col
  reductions <- list(before_batch_correction = before, after_batch_correction = after)

  for (stage in names(reductions)) {
    reduction <- reductions[[stage]]
    coordinates <- Embeddings(sc, reduction = reduction)
    xlim <- range(coordinates[, 1])
    ylim <- range(coordinates[, 2])

    p_sample <- DimPlot(
      sc,
      reduction = reduction,
      group.by = sample_col,
      pt.size = BATCH_EFFECT_UMAP$point_size,
      raster = TRUE
    ) + ggtitle(gsub("_", " ", stage))
    save_umap_plot(
      p_sample,
      file.path(GLOBAL_DIR, paste0("UMAP_", stage, "_by_sample")),
      xlim = xlim,
      ylim = ylim
    )

    p_split <- DimPlot(
      sc,
      reduction = reduction,
      group.by = "seurat_clusters",
      split.by = sample_col,
      ncol = BATCH_EFFECT_UMAP$split_ncol,
      pt.size = BATCH_EFFECT_UMAP$split_point_size,
      raster = TRUE
    )
    save_umap_plot(
      p_split,
      file.path(GLOBAL_DIR, paste0("UMAP_", stage, "_split_by_sample")),
      width = BATCH_EFFECT_UMAP$split_width,
      height = BATCH_EFFECT_UMAP$split_height,
      xlim = xlim,
      ylim = ylim
    )
  }
  sync_final_umap_aliases()
}

sync_final_umap_aliases <- function() {
  for (extension in c("png", "pdf")) {
    file.copy(
      file.path(GLOBAL_DIR, paste0("UMAP_after_batch_correction_by_sample.", extension)),
      file.path(GLOBAL_DIR, paste0("UMAP_by_sample.", extension)),
      overwrite = TRUE
    )
    file.copy(
      file.path(GLOBAL_DIR, paste0("UMAP_after_batch_correction_split_by_sample.", extension)),
      file.path(GLOBAL_DIR, paste0("UMAP_split_by_sample.", extension)),
      overwrite = TRUE
    )
  }
}

save_contamination_umap_output <- function(sc) {
  if (!"Contamination" %in% colnames(sc@meta.data)) return(invisible(NULL))

  contamination <- sc$Contamination
  threshold <- QC_DEFAULTS$contamination_max
  legend_max <- CONTAMINATION_UMAP$legend_max
  summary_table <- data.frame(
    cells = length(contamination),
    minimum = min(contamination),
    p05 = unname(quantile(contamination, 0.05)),
    median = median(contamination),
    mean = mean(contamination),
    p95 = unname(quantile(contamination, 0.95)),
    p99 = unname(quantile(contamination, 0.99)),
    maximum = max(contamination),
    filter_threshold = threshold,
    cells_ge_filter_threshold = sum(contamination >= threshold)
  )
  write.csv(summary_table, file.path(GLOBAL_DIR, "Contamination_summary.csv"), row.names = FALSE, quote = FALSE)

  p <- FeaturePlot(
    sc,
    features = "Contamination",
    pt.size = CONTAMINATION_UMAP$point_size,
    order = TRUE,
    raster = TRUE,
    min.cutoff = CONTAMINATION_UMAP$min_cutoff,
    max.cutoff = legend_max,
    cols = CONTAMINATION_UMAP$colors[c(1, length(CONTAMINATION_UMAP$colors))]
  ) + scale_color_gradientn(
    colours = CONTAMINATION_UMAP$colors,
    limits = c(CONTAMINATION_UMAP$min_cutoff, legend_max),
    breaks = seq(CONTAMINATION_UMAP$min_cutoff, legend_max, by = 0.25),
    labels = scales::label_number(accuracy = 0.01),
    oob = scales::squish
  ) + labs(
    title = sprintf("Contamination after QC (cutoff < %.2f)", threshold),
    subtitle = sprintf("Median %.3f | P95 %.3f | all %s cells below cutoff",
      summary_table$median, summary_table$p95, format(summary_table$cells, big.mark = ",")),
    color = "Contamination"
  )
  save_umap_plot(p, file.path(GLOBAL_DIR, "Feature_Contamination"))
  invisible(p)
}

save_global_clustering_outputs <- function(sc) {
  cluster_counts <- as.data.frame(table(cluster = sc$seurat_clusters))
  cluster_sample_counts <- as.data.frame(table(cluster = sc$seurat_clusters, sample = sc$samples))
  positive <- cluster_sample_counts[cluster_sample_counts$Freq > 0, ]
  coverage <- do.call(rbind, lapply(split(positive, positive$cluster), function(x) {
    top <- x[which.max(x$Freq), ]
    data.frame(
      cluster = x$cluster[1],
      samples_present = nrow(x),
      cells = sum(x$Freq),
      top_sample = top$sample,
      top_sample_cells = top$Freq
    )
  }))
  parameters <- data.frame(
    parameter = c("cells", "samples", "normalization", "sct_features", "regress_vars", "batch_method", "batch_group", "pca_dims", "clustering_dims", "resolution", "clusters"),
    value = c(ncol(sc), length(unique(sc$orig.ident)), NORMALIZE$method, NORMALIZE$nfeatures, paste(NORMALIZE$regress_vars, collapse = ","), CLUSTER$batch_method, CLUSTER$batch_group, CLUSTER$npcs, paste0(min(CLUSTER$dims), ":", max(CLUSTER$dims)), CLUSTER$resolution, nrow(cluster_counts))
  )
  write.csv(cluster_counts, file.path(GLOBAL_DIR, "cluster_cell_counts.csv"), row.names = FALSE, quote = FALSE)
  write.csv(cluster_sample_counts, file.path(GLOBAL_DIR, "cluster_by_sample_counts.csv"), row.names = FALSE, quote = FALSE)
  write.csv(coverage, file.path(GLOBAL_DIR, "cluster_sample_coverage.csv"), row.names = FALSE, quote = FALSE)
  write.csv(parameters, file.path(GLOBAL_DIR, "clustering_parameters.csv"), row.names = FALSE, quote = FALSE)

  save_umap_plot(DimPlot(sc, label = TRUE, repel = TRUE), file.path(GLOBAL_DIR, "UMAP_clusters"))
  save_batch_effect_umap_outputs(sc)
  save_contamination_umap_output(sc)

  save_qs2(sc, QS2_CLUSTER)
  invisible(sc)
}

step_global_clustering <- function() {
  log_msg("Step 03 normalization, batch correction and global clustering started")

  sc <- read_qs2(QS2_QC)
  sc <- assign_group_metadata(sc, required_sheet = TRUE)
  sc <- run_normalize_scale(sc)
  sc <- run_global_clustering(sc)
  save_global_clustering_outputs(sc)
  log_msg("Step 03 done. Clusters: ", length(unique(sc$seurat_clusters)))
}

step_global_recluster <- function() {
  log_msg("Step 03b reclustering started. PCs: ", paste0(min(CLUSTER$dims), ":", max(CLUSTER$dims)))

  sc <- read_qs2(QS2_CLUSTER)
  DefaultAssay(sc) <- NORMALIZE$assay
  sc[["pca"]] <- NULL
  sc[["harmony"]] <- NULL
  sc[["umap"]] <- NULL
  sc@graphs <- list()

  sc <- run_step("RunPCA", RunPCA(sc, npcs = CLUSTER$npcs, verbose = FALSE))
  save_pca_elbow(sc, file.path(GLOBAL_DIR, "PCA_elbow"))
  sc <- run_before_batch_umap(sc, CLUSTER$dims)
  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Global")
  sc <- run_graph_cluster_umap(
    corrected$object,
    reduction = corrected$reduction,
    dims = CLUSTER$dims,
    resolution = CLUSTER$resolution,
    context_label = "Global"
  )
  sc <- store_clustering_metadata(sc, "global", corrected$reduction, CLUSTER$dims)
  save_global_clustering_outputs(sc)
  log_msg("Step 03b done. Clusters: ", length(unique(sc$seurat_clusters)))
}
step_batch_effect_umap <- function() {
  log_msg("Step 03c batch-effect UMAP comparison started")
  sc <- read_qs2(QS2_CLUSTER)
  sc <- run_before_batch_umap(sc, CLUSTER$dims)
  save_batch_effect_umap_outputs(sc)
  save_qs2(sc, QS2_CLUSTER)
  log_msg("Step 03c done")
}

step_redraw_global_umaps <- function() {
  log_msg("Step 03d redraw batch-effect and contamination UMAPs started")
  sc <- read_qs2(QS2_CLUSTER)
  save_batch_effect_umap_outputs(sc)
  save_contamination_umap_output(sc)
  log_msg("Step 03d done")
}

step_redraw_contamination_umap <- function() {
  log_msg("Step 03e redraw continuous Contamination legend started")
  sc <- read_qs2(QS2_CLUSTER)
  save_contamination_umap_output(sc)
  log_msg("Step 03e done")
}
