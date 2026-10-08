# Prepare a reproducible PCA checkpoint; selection is a required user decision.
step_prepare_harmony_pca <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  future::plan(future::sequential)
  settings <- CONFIG$harmony_reintegration
  NORMALIZE <- utils::modifyList(NORMALIZE, settings$normalize, keep.null = TRUE)
  NORMALIZE$regress_vars <- unlist(NORMALIZE$regress_vars, use.names = FALSE)
  sys.source(file.path(WORKSPACE_DIR, "src/R/05_normalization.R"), envir = environment())
  sys.source(file.path(WORKSPACE_DIR, "src/R/07_clustering.R"), envir = environment())
  output <- file.path(PROJECT_DIR, settings$output_dir)
  checkpoint <- file.path(PROJECT_DIR, settings$pca_checkpoint)
  stopifnot(!file.exists(checkpoint), !file.exists(paste0(checkpoint, ".partial")),
            is.null(settings$selected_pcs))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  set.seed(SEED)
  sc <- run_step("Read merged counts", qs::qread(file.path(PROJECT_DIR, CONFIG$merged_output), nthreads = 4))
  original_metadata <- sc[[]]
  original_cells <- colnames(sc)
  original_genes <- rownames(sc)
  original_reductions <- sc@reductions
  original_counts <- Matrix::colSums(LayerData(sc, assay = "RNA", layer = "counts"))
  stopifnot(!anyDuplicated(original_cells), settings$batch_col %in% names(original_metadata),
            !anyNA(original_metadata[[settings$batch_col]]),
            all(nzchar(as.character(original_metadata[[settings$batch_col]]))))
  sc <- run_normalize_scale(sc)
  sc <- run_step("RunPCA candidate components", RunPCA(
    sc, assay = NORMALIZE$assay, features = VariableFeatures(sc),
    npcs = settings$candidate_npcs, seed.use = SEED, reduction.name = "pca", verbose = FALSE
  ))
  stopifnot(ncol(Embeddings(sc[["pca"]])) == settings$candidate_npcs,
            all(is.finite(Embeddings(sc[["pca"]]))), identical(colnames(sc), original_cells),
            identical(rownames(sc), original_genes), identical(sc[[]], original_metadata),
            identical(Matrix::colSums(LayerData(sc, assay = "RNA", layer = "counts")), original_counts))
  for (name in names(original_reductions)) stopifnot(identical(sc[[name]], original_reductions[[name]]))

  elbow <- detect_pca_elbow(sc)
  scaled <- LayerData(sc, assay = NORMALIZE$assay, layer = "scale.data")
  pca_features <- rownames(Loadings(sc[["pca"]]))
  total_variance <- sum(matrixStats::rowVars(scaled, rows = match(pca_features, rownames(scaled))))
  variance <- data.frame(
    pc = elbow$variance$pc, standard_deviation = Stdev(sc, reduction = "pca"),
    eigenvalue = elbow$variance$variance,
    explained_variance_percent = elbow$variance$variance / total_variance * 100,
    cumulative_variance_percent = cumsum(elbow$variance$variance) / total_variance * 100
  )
  stopifnot(all(diff(variance$standard_deviation) <= 1e-8),
            all(variance$cumulative_variance_percent > 0),
            max(variance$cumulative_variance_percent) <= 100 + 1e-6)
  write.table(variance, file.path(output, "PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(elbow$selections, file.path(output, "PCA_automatic_references.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(VariableFeatures(sc), file.path(output, "PCA_variable_features.txt"))
  write.table(Embeddings(sc[["pca"]]), gzfile(file.path(output, "PCA_embeddings.tsv.gz")),
              sep = "\t", quote = FALSE, col.names = NA)
  method_colors <- c(PCAtools = HONGF_NPG[1], nFactors_noc = HONGF_NPG[4], nFactors_naf = HONGF_NPG[3])
  elbow_plot <- ggplot(variance, aes(pc, standard_deviation)) +
    geom_line(color = "#555555", linewidth = .5) + geom_point(size = 1.5, color = "#555555") +
    geom_vline(data = elbow$selections, aes(xintercept = pc, color = method, linetype = method), linewidth = .6) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 5)) +
    labs(title = "PCA elbow", x = "Principal component", y = "Standard deviation", color = NULL, linetype = NULL) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom")
  cumulative_plot <- ggplot(variance, aes(pc, cumulative_variance_percent)) +
    geom_line(color = HONGF_NPG[4], linewidth = .6) + geom_point(size = 1.4, color = HONGF_NPG[4]) +
    scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 5)) +
    labs(title = "Cumulative explained variance", x = "Number of principal components", y = "Total variance explained (%)") +
    theme_classic(base_size = 11)
  save_plot(patchwork::wrap_plots(elbow_plot, cumulative_plot, ncol = 2),
            file.path(output, "PCA_selection_diagnostics"), width = 11, height = 4.5)
  sc@misc$harmony_reintegration <- list(
    stage = "awaiting_user_pc_selection", candidate_npcs = settings$candidate_npcs,
    selected_pcs = NULL, planned_batch_column = settings$batch_col,
    normalization = NORMALIZE, seed = SEED, automatic_references = elbow$selections,
    variance = variance, input = CONFIG$merged_output
  )
  stopifnot(!"harmony" %in% Reductions(sc), !"umap" %in% Reductions(sc))
  expected_pca <- Embeddings(sc[["pca"]])
  rm(scaled, original_reductions)
  invisible(gc())
  run_step("Save PCA checkpoint", qs::qsave(sc, paste0(checkpoint, ".partial"), preset = "high", nthreads = 4))
  rm(sc)
  invisible(gc())
  check <- run_step("Read back PCA checkpoint", qs::qread(paste0(checkpoint, ".partial"), nthreads = 4))
  stopifnot(identical(colnames(check), original_cells), identical(rownames(check), original_genes),
            identical(check[[]], original_metadata), identical(Embeddings(check[["pca"]]), expected_pca),
            identical(Matrix::colSums(LayerData(check, assay = "RNA", layer = "counts")), original_counts),
            is.null(check@misc$harmony_reintegration$selected_pcs),
            isTRUE(validObject(check, test = TRUE)), file.rename(paste0(checkpoint, ".partial"), checkpoint))
  yaml::write_yaml(list(stage = "awaiting_user_pc_selection", candidate_npcs = settings$candidate_npcs,
                       selected_pcs = NULL, batch_column = settings$batch_col,
                       checkpoint = settings$pca_checkpoint), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "05_prepare_harmony_pca_sessionInfo.txt"))
  print(elbow$selections, row.names = FALSE)
  print(variance[variance$pc %in% c(10, 15, 20, 25, 30, 40, 50), ], row.names = FALSE)
  message("PCA_CHECKPOINT_OK: ", ncol(check), " cells; ", settings$candidate_npcs,
          " candidate PCs. Harmony, clustering and UMAP await the user's PC selection.")
  invisible(checkpoint)
}
