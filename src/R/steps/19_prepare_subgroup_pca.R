# 各大类独立预处理，只运行至 PCA；任何自动拐点均不代替用户选 PC。
step_prepare_subgroup_pca <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(PCAtools))
  future::plan(future::sequential)
  options(Seurat.object.assay.version = "v3")
  settings <- CONFIG$subgroup_pca
  groups <- unlist(settings$groups, use.names = FALSE)
  stopifnot(setequal(names(settings$selected_pcs), groups),
            all(vapply(settings$selected_pcs, is.null, logical(1))))
  NORMALIZE <- utils::modifyList(NORMALIZE, settings$normalize, keep.null = TRUE)
  NORMALIZE$regress_vars <- unlist(NORMALIZE$regress_vars, use.names = FALSE)
  sys.source(file.path(WORKSPACE_DIR, "src/R/05_normalization.R"), envir = environment())
  sys.source(file.path(WORKSPACE_DIR, "src/R/07_clustering.R"), envir = environment())
  output <- file.path(PROJECT_DIR, settings$output_dir)
  checkpoint_dir <- file.path(PROJECT_DIR, settings$checkpoint_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- setNames(file.path(checkpoint_dir, paste0("12.", CONFIG$project_id, "_", groups,
                        "_pca", settings$candidate_npcs, "_pending_selection.qs")), groups)
  stopifnot(!any(file.exists(paths)), !any(file.exists(paste0(paths, ".partial"))))
  parent <- run_step("Read complete annotated parent", qs::qread(file.path(PROJECT_DIR, settings$input), nthreads = 4))
  metadata <- parent[[]]
  counts <- LayerData(parent, assay = "RNA", layer = "counts")
  parent_cells <- colnames(parent)
  parent_genes <- rownames(parent)
  stopifnot(!anyNA(metadata[[settings$parent_col]]), !anyDuplicated(parent_cells),
            setequal(unique(as.character(metadata[[settings$parent_col]])), groups))
  # 母对象留在原文件；子对象从原始 counts 新建，不混入全局 HVG、scale.data 或降维。
  rm(parent)
  invisible(gc())
  summary_rows <- variance_rows <- reference_rows <- completed_cells <- list()
  colors <- c(PCAtools = HONGF_NPG[1], nFactors_noc = HONGF_NPG[4], nFactors_naf = HONGF_NPG[3])
  layer_hashes <- function(x) vapply(c("counts", "data", "scale.data"), function(layer)
    digest::digest(LayerData(x, assay = "RNA", layer = layer), algo = "xxhash64"), character(1))

  for (group in groups) {
    log_msg("SUBGROUP_START: ", group)
    cell_ids <- parent_cells[as.character(metadata[[settings$parent_col]]) == group]
    original_metadata <- metadata[cell_ids, , drop = FALSE]
    original_counts <- counts[, cell_ids, drop = FALSE]
    original_counts_hash <- digest::digest(original_counts, algo = "xxhash64")
    sc <- CreateSeuratObject(counts = original_counts, project = group, min.cells = 0, min.features = 0,
                              meta.data = original_metadata)
    # 保留母对象元数据的原始值、类型及列顺序，去掉新建对象自动补的 orig.ident。
    sc@meta.data <- original_metadata
    Idents(sc) <- original_metadata[[settings$parent_col]]
    rm(original_counts)
    stopifnot(identical(rownames(sc), parent_genes), identical(colnames(sc), cell_ids),
              identical(digest::digest(LayerData(sc, assay = "RNA", layer = "counts"), algo = "xxhash64"), original_counts_hash))
    set.seed(SEED)
    sc <- run_normalize_scale(sc)
    sc <- run_step(paste(group, "RunPCA"), RunPCA(sc, assay = "RNA", features = VariableFeatures(sc),
      npcs = settings$candidate_npcs, seed.use = SEED, reduction.name = "pca", verbose = FALSE))
    stopifnot(identical(sc[[]], original_metadata), ncol(Embeddings(sc, "pca")) == settings$candidate_npcs,
              all(is.finite(Embeddings(sc, "pca"))), identical(Reductions(sc), "pca"), length(sc@graphs) == 0L)
    elbow <- detect_pca_elbow(sc)
    scaled <- LayerData(sc, assay = "RNA", layer = "scale.data")
    pca_features <- rownames(Loadings(sc[["pca"]]))
    total_variance <- sum(matrixStats::rowVars(scaled, rows = match(pca_features, rownames(scaled))))
    pc_sd <- Stdev(sc, "pca")
    variance <- data.frame(subgroup = group, pc = seq_along(pc_sd), standard_deviation = pc_sd,
      eigenvalue = pc_sd^2, explained_variance_percent = pc_sd^2 / total_variance * 100,
      cumulative_variance_percent = cumsum(pc_sd^2) / total_variance * 100)
    stopifnot(all(diff(pc_sd) <= 1e-8), all(variance$cumulative_variance_percent > 0),
              max(variance$cumulative_variance_percent) <= 100 + 1e-6)
    references <- data.frame(subgroup = group, elbow$selections)
    group_dir <- file.path(output, group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(variance, file.path(group_dir, "PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(references, file.path(group_dir, "PCA_automatic_references.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(VariableFeatures(sc), file.path(group_dir, "variable_features.txt"))
    write.table(as.data.frame(table(sample = original_metadata[[settings$sample_col]])),
      file.path(group_dir, "cell_counts_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    label <- gsub("_", " ", group)
    p1 <- ggplot(variance, aes(pc, standard_deviation)) + geom_line(color = "#555555", linewidth = .5) +
      geom_point(color = "#555555", size = 1.4) +
      geom_vline(data = references, aes(xintercept = pc, color = method, linetype = method), linewidth = .6) +
      scale_color_manual(values = colors) + scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 5)) +
      labs(title = paste(label, "PCA"), x = "Principal component", y = "Standard deviation", color = NULL, linetype = NULL) +
      theme_classic(base_size = 11) + theme(legend.position = "bottom")
    p2 <- ggplot(variance, aes(pc, cumulative_variance_percent)) + geom_line(color = HONGF_NPG[4], linewidth = .6) +
      geom_point(color = HONGF_NPG[4], size = 1.4) + scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 5)) +
      labs(title = "Cumulative explained variance", x = "Number of principal components", y = "Total variance explained (%)") +
      theme_classic(base_size = 11)
    save_plot(patchwork::wrap_plots(p1, p2, ncol = 2), file.path(group_dir, "PCA_selection_diagnostics"), width = 11, height = 4.5)
    sc@misc$subgroup_pca <- list(stage = "awaiting_user_pc_selection", subgroup = group,
      parent_input = settings$input, parent_column = settings$parent_col, candidate_npcs = settings$candidate_npcs,
      selected_pcs = NULL, planned_batch_column = settings$sample_col, normalization = NORMALIZE, seed = SEED,
      automatic_references = references, variance = variance, parent_counts_hash = original_counts_hash,
      inherited_cluster_columns_are_parent_labels = TRUE)
    expected_layers <- layer_hashes(sc)
    expected_pca <- digest::digest(sc[["pca"]], algo = "xxhash64")
    expected_features <- VariableFeatures(sc)
    stopifnot(identical(unname(expected_layers["counts"]), original_counts_hash))
    rm(scaled)
    invisible(gc())
    checkpoint <- paths[[group]]
    run_step(paste(group, "Save PCA QS"), qs::qsave(sc, paste0(checkpoint, ".partial"), preset = "high", nthreads = 4))
    rm(sc)
    invisible(gc())
    check <- run_step(paste(group, "Read back PCA QS"), qs::qread(paste0(checkpoint, ".partial"), nthreads = 4))
    stopifnot(isTRUE(validObject(check, test = TRUE)), identical(colnames(check), cell_ids),
      identical(rownames(check), parent_genes), identical(check[[]], original_metadata),
      identical(layer_hashes(check), expected_layers), identical(VariableFeatures(check), expected_features),
      identical(digest::digest(check[["pca"]], algo = "xxhash64"), expected_pca),
      identical(Reductions(check), "pca"), length(check@graphs) == 0L,
      is.null(check@misc$subgroup_pca$selected_pcs), file.rename(paste0(checkpoint, ".partial"), checkpoint))
    row <- data.frame(subgroup = group, cells = ncol(check), genes = nrow(check),
      samples = length(unique(original_metadata[[settings$sample_col]])), variable_genes = length(expected_features),
      pca_features = length(pca_features), candidate_pcs = ncol(Embeddings(check, "pca")), selected_pcs = NA_character_,
      PCAtools_reference = references$pc[references$method == "PCAtools"],
      nFactors_noc_reference = references$pc[references$method == "nFactors_noc"],
      nFactors_naf_reference = references$pc[references$method == "nFactors_naf"],
      checkpoint = file.path(settings$checkpoint_dir, basename(checkpoint)))
    for (n in c(10, 15, 20, 25, 30, 40, 50)) row[[paste0("variance_PC", n)]] <- variance$cumulative_variance_percent[variance$pc == n]
    summary_rows[[group]] <- row
    variance_rows[[group]] <- variance
    reference_rows[[group]] <- references
    completed_cells[[group]] <- colnames(check)
    yaml::write_yaml(list(stage = "awaiting_user_pc_selection", subgroup = group, cells = ncol(check),
      candidate_pcs = settings$candidate_npcs, selected_pcs = NULL, checkpoint = row$checkpoint),
      file.path(group_dir, "analysis_status.yaml"))
    log_msg("SUBGROUP_PCA_CHECKPOINT_OK: ", group, "; cells=", ncol(check), "; candidate PCs=", settings$candidate_npcs)
    rm(check)
    invisible(gc())
  }
  all_cells <- unlist(completed_cells, use.names = FALSE)
  stopifnot(!anyDuplicated(all_cells), setequal(all_cells, parent_cells), length(all_cells) == length(parent_cells))
  summary <- do.call(rbind, summary_rows)
  variance <- do.call(rbind, variance_rows)
  references <- do.call(rbind, reference_rows)
  write.table(summary, file.path(output, "subgroup_PCA_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(variance, file.path(output, "all_subgroups_PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(references, file.path(output, "all_subgroups_PCA_references.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  variance$subgroup <- factor(variance$subgroup, levels = groups)
  references$subgroup <- factor(references$subgroup, levels = groups)
  facet_labels <- setNames(gsub("_", " ", groups), groups)
  p <- ggplot(variance, aes(pc, standard_deviation)) + geom_line(color = "#555555", linewidth = .5) +
    geom_point(color = "#555555", size = 1.2) +
    geom_vline(data = references, aes(xintercept = pc, color = method, linetype = method), linewidth = .55) +
    scale_color_manual(values = colors) + facet_wrap(~subgroup, ncol = 3, scales = "free_y", labeller = labeller(subgroup = facet_labels)) +
    scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 10)) +
    labs(x = "Principal component", y = "Standard deviation", color = NULL, linetype = NULL) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom")
  save_plot(p, file.path(output, "All_subgroups_PCA_elbow"), width = 14, height = 8)
  p <- ggplot(variance, aes(pc, cumulative_variance_percent)) + geom_line(color = HONGF_NPG[4], linewidth = .6) +
    geom_point(color = HONGF_NPG[4], size = 1.2) +
    facet_wrap(~subgroup, ncol = 3, labeller = labeller(subgroup = facet_labels)) +
    scale_x_continuous(breaks = seq(0, settings$candidate_npcs, 10)) +
    labs(x = "Number of principal components", y = "Total variance explained (%)") + theme_classic(base_size = 11)
  save_plot(p, file.path(output, "All_subgroups_PCA_cumulative_variance"), width = 14, height = 8)
  yaml::write_yaml(list(stage = "awaiting_user_pc_selection", completed_groups = groups,
    total_cells = length(all_cells), candidate_pcs = settings$candidate_npcs,
    selected_pcs = settings$selected_pcs, harmony_run = FALSE, clustering_run = FALSE, umap_run = FALSE),
    file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "12_prepare_subgroup_pca_sessionInfo.txt"))
  print(summary, row.names = FALSE)
  log_msg("ALL_SUBGROUPS_PCA_OK: ", length(groups), " groups; ", length(all_cells), " cells; waiting for user PC choices")
  invisible(summary)
}
