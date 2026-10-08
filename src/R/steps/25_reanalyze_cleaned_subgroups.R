# 从清理后的counts重建各大类对象，完整重算预处理、PCA、Harmony、聚类和UMAP。
step_reanalyze_cleaned_subgroups <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  future::plan(future::sequential)
  options(Seurat.object.assay.version = "v3")
  settings <- utils::modifyList(CONFIG$subgroup_harmony, CONFIG$subgroup_reanalysis)
  selected <- setNames(rep(NA_integer_, length(CONFIG$subgroup_pca$groups)), unlist(CONFIG$subgroup_pca$groups))
  NORMALIZE <- utils::modifyList(NORMALIZE, CONFIG$subgroup_pca$normalize, keep.null = TRUE)
  NORMALIZE$regress_vars <- unlist(NORMALIZE$regress_vars, use.names = FALSE)
  sys.source(file.path(WORKSPACE_DIR, "src/R/05_normalization.R"), envir = environment())
  sys.source(file.path(WORKSPACE_DIR, "src/R/07_clustering.R"), envir = environment())
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest))
  stopifnot(setequal(manifest$subgroup, names(selected)), sum(manifest$retained_cells) == settings$expected_cells,
    identical(settings$pc_selection_method, "PCAtools"))
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROJECT_DIR, settings$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(settings$checkpoint_dir, paste0(settings$object_prefix, ".", CONFIG$project_id, "_",
    manifest$subgroup, "_reanalyzed_autoPC_res", settings$resolution, ".qs"))
  stopifnot(!any(file.exists(file.path(PROJECT_DIR, paths))), !any(file.exists(paste0(file.path(PROJECT_DIR, paths), ".partial"))))
  fingerprints <- function(x) c(RNA = digest::digest(x[["RNA"]], algo = "xxhash64"),
    metadata = digest::digest(x[[]], algo = "xxhash64"), reductions = digest::digest(x@reductions, algo = "xxhash64"),
    graphs = digest::digest(x@graphs, algo = "xxhash64"))
  summaries <- all_variance <- all_references <- cell_sets <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    folder <- file.path(output, group); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    input_hash <- digest::digest(file = input, algo = "sha256")
    old <- run_step(paste(group, "Read cleaned QS"), qs::qread(input, nthreads = 4))
    cells <- colnames(old); genes <- rownames(old); metadata <- old[[]]
    counts <- LayerData(old, assay = "RNA", layer = "counts")
    counts_hash <- digest::digest(counts, algo = "xxhash64")
    old_hvg <- VariableFeatures(old); old_umap <- Embeddings(old, "umap")
    stopifnot(ncol(old) == manifest$retained_cells[i], !anyNA(metadata$group),
      all(as.character(metadata$group) %in% CONFIG$group_levels))
    # 历史聚类/候选标签显式加前缀；不能将旧cluster注释冒充新cluster注释。
    historical <- names(metadata)[grepl("^candidate_", names(metadata)) |
      names(metadata) %in% c("subgroup_cluster", "seurat_clusters", "annotation_stage")]
    stopifnot(!any(paste0(settings$historical_prefix, historical) %in% names(metadata)))
    names(metadata)[match(historical, names(metadata))] <- paste0(settings$historical_prefix, historical)
    sc <- CreateSeuratObject(counts = counts, project = group, min.cells = 0, min.features = 0, meta.data = metadata)
    sc@meta.data <- metadata
    Idents(sc) <- "broad_cell_type"
    rm(old, counts); invisible(gc())
    set.seed(SEED)
    sc <- run_normalize_scale(sc)
    sc <- run_step(paste(group, "Fresh PCA"), RunPCA(sc, assay = "RNA", features = VariableFeatures(sc),
      npcs = CONFIG$subgroup_pca$candidate_npcs, seed.use = SEED, verbose = FALSE))
    stopifnot(identical(sc[[]], metadata), identical(Reductions(sc), "pca"), length(sc@graphs) == 0)
    elbow <- detect_pca_elbow(sc)
    selected[[group]] <- elbow$recommended_pc
    stopifnot(length(selected[[group]]) == 1, !is.na(selected[[group]]),
      selected[[group]] >= 2, selected[[group]] <= CONFIG$subgroup_pca$candidate_npcs)
    dims <- seq_len(selected[[group]])
    references <- data.frame(subgroup = group, elbow$selections,
      used_for_analysis = elbow$selections$method == settings$pc_selection_method)
    all_references[[group]] <- references
    write.table(references, file.path(folder, "PCA_automatic_references.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("CLEANED_AUTO_PC_SELECTED ", group, "; PCAtools=", length(dims),
      "; nFactors_noc=", references$pc[references$method == "nFactors_noc"])
    pc_sd <- Stdev(sc, "pca")
    scaled <- LayerData(sc, assay = "RNA", layer = "scale.data")
    total_var <- sum(matrixStats::rowVars(scaled, rows = match(rownames(Loadings(sc, "pca")), rownames(scaled))))
    variance <- data.frame(subgroup = group, pc = seq_along(pc_sd), standard_deviation = pc_sd,
      explained_variance_percent = pc_sd^2 / total_var * 100, cumulative_variance_percent = cumsum(pc_sd^2) / total_var * 100,
      selected_pcs = length(dims))
    write.table(variance, file.path(folder, "PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(VariableFeatures(sc), file.path(folder, "variable_features.txt"))
    rm(scaled)
    sc <- run_step(paste(group, "Fresh Harmony"), harmony::RunHarmony(sc,
      group.by.vars = settings$batch_col, reduction.use = "pca", dims.use = dims,
      reduction.save = "harmony", project.dim = FALSE, theta = settings$theta, lambda = settings$lambda,
      max_iter = settings$max_iter, ncores = settings$ncores, early_stop = TRUE, verbose = TRUE))
    sc <- run_step(paste(group, "Fresh neighbors"), FindNeighbors(sc, reduction = "harmony", dims = dims,
      k.param = settings$k_neighbors, graph.name = c("subgroup_nn", "subgroup_snn"), verbose = TRUE))
    sc <- run_step(paste(group, "Fresh clusters"), FindClusters(sc, graph.name = "subgroup_snn", cluster.name = "subgroup_cluster",
      resolution = settings$resolution, algorithm = settings$algorithm, random.seed = SEED, verbose = TRUE))
    sc$seurat_clusters <- sc$subgroup_cluster
    for (reduction in c("pca", "harmony")) sc <- run_step(paste(group, "Fresh UMAP", reduction), RunUMAP(sc,
      reduction = reduction, dims = dims, reduction.name = if (reduction == "pca") "umap.pca" else "umap",
      reduction.key = if (reduction == "pca") "UMAPPCA_" else "UMAP_", n.neighbors = settings$umap_neighbors,
      min.dist = settings$umap_min_dist, metric = settings$umap_metric, seed.use = SEED, verbose = TRUE))
    sc$annotation_stage <- settings$stage
    sc@misc$subgroup_reanalysis <- list(stage = settings$stage, input = manifest$checkpoint[i], input_sha256 = input_hash,
      selected_pcs = dims, selection_source = "user_authorized_fresh_PCAtools_elbow", automatic_references = references, settings = settings,
      normalization = NORMALIZE, candidate_npcs = CONFIG$subgroup_pca$candidate_npcs, seed = SEED,
      rebuilt_from_counts = TRUE, counts_hash = counts_hash, history_columns = paste0(settings$historical_prefix, historical))
    stopifnot(identical(sc[[]][, names(metadata), drop = FALSE], metadata), identical(colnames(sc), cells),
      identical(rownames(sc), genes), identical(digest::digest(LayerData(sc, assay = "RNA", layer = "counts"), algo = "xxhash64"), counts_hash),
      ncol(Embeddings(sc, "harmony")) == length(dims), !anyNA(sc$subgroup_cluster),
      !identical(as.numeric(Embeddings(sc, "umap")), as.numeric(old_umap)))
    for (reduction in Reductions(sc)) stopifnot(identical(rownames(Embeddings(sc, reduction)), cells), all(is.finite(Embeddings(sc, reduction))))
    for (column in c("pre_reanalysis_subgroup_cluster", "pre_reanalysis_candidate_annotation")) {
      d <- as.data.frame(table(subgroup_cluster = sc$subgroup_cluster, previous = sc[[]][[column]]))
      names(d)[2:3] <- c(column, "cells")
      write.table(d, file.path(folder, paste0("cluster_by_", column, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    before <- fingerprints(sc)
    destination <- file.path(PROJECT_DIR, paths[i])
    run_step(paste(group, "Save reanalyzed QS"), qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4))
    hvg_overlap <- length(intersect(VariableFeatures(sc), old_hvg)); rm(sc); invisible(gc())
    check <- run_step(paste(group, "Read back reanalyzed QS"), qs::qread(paste0(destination, ".partial"), nthreads = 4))
    stopifnot(isTRUE(validObject(check, test = TRUE)), identical(fingerprints(check), before),
      identical(colnames(check), cells), identical(rownames(check), genes),
      identical(as.character(Idents(check)), as.character(check$subgroup_cluster)),
      identical(digest::digest(file = input, algo = "sha256"), input_hash), file.rename(paste0(destination, ".partial"), destination))
    summaries[[group]] <- data.frame(subgroup = group, cells = ncol(check), genes = nrow(check),
      samples = length(unique(check$sample)), selected_pcs = length(dims), resolution = settings$resolution,
      clusters = length(unique(check$subgroup_cluster)), previous_clusters = manifest$retained_clusters[i],
      hvg_overlap_previous = hvg_overlap, checkpoint = paths[i], input_sha256 = input_hash,
      counts_preserved = TRUE, readback_verified = TRUE)
    all_variance[[group]] <- variance; cell_sets[[group]] <- cells
    write.table(do.call(rbind, summaries), file.path(output, "subgroup_integration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(do.call(rbind, all_references), file.path(output, "all_subgroups_PCA_references.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    yaml::write_yaml(list(stage = settings$stage, selected_pc_indices = dims, selection_method = settings$pc_selection_method,
      settings = settings, normalization = NORMALIZE, input = manifest$checkpoint[i], input_sha256 = input_hash),
      file.path(folder, "integration_parameters.yaml"))
    log_msg("CLEANED_REANALYSIS_OK ", group, "; cells=", ncol(check), "; clusters=", length(unique(check$subgroup_cluster)))
    rm(check); invisible(gc())
  }
  all_cells <- unlist(cell_sets, use.names = FALSE)
  stopifnot(!anyDuplicated(all_cells), length(all_cells) == settings$expected_cells)
  variance <- do.call(rbind, all_variance)
  write.table(variance, file.path(output, "all_subgroups_PCA_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  vlines <- unique(variance[, c("subgroup", "selected_pcs")])
  p <- ggplot(variance, aes(pc, explained_variance_percent)) + geom_line(color = "#555555", linewidth = .5) +
    geom_point(color = "#555555", size = 1.2) + geom_vline(data = vlines, aes(xintercept = selected_pcs), color = HONGF_NPG[1], linetype = 2) +
    facet_wrap(~subgroup, ncol = 3, scales = "free_y") + labs(x = "Principal component", y = "Explained variance (%)") + theme_classic()
  save_plot(p, file.path(output, "All_subgroups_recomputed_PCA_elbow"), width = 14, height = 8)
  yaml::write_yaml(list(stage = settings$stage, total_cells = length(all_cells), selected_pcs = as.list(selected),
    settings = settings, normalization = NORMALIZE, rebuilt_from_counts = TRUE, input_objects_unchanged = TRUE), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "24_reanalyze_cleaned_subgroups_sessionInfo.txt"))
  log_msg("ALL_CLEANED_SUBGROUPS_REANALYZED_OK ", length(all_cells), " cells")
}

# 在新cluster上输出同一文献marker面板；旧候选名称仅在组成表中对照，不按多数票继承。
step_plot_reanalyzed_markers <- function() {
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$subgroup_reanalysis
  output <- file.path(PROJECT_DIR, settings$output_dir)
  manifest <- read.delim(file.path(output, "subgroup_integration_summary.tsv"))
  panel <- read.delim(file.path(PROJECT_DIR, CONFIG$subgroup_marker_review$marker_panel))
  summaries <- plots <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; folder <- file.path(output, group)
    e <- readRDS(file.path(folder, "expression_evidence.rds"))
    counts <- table(e$metadata$subgroup_cluster)
    mapping <- data.frame(subgroup_cluster = as.integer(names(counts)), cells = as.integer(counts),
      candidate_annotation = "review")
    mapping <- mapping[order(mapping$subgroup_cluster), ]
    genes <- strsplit(panel$genes[panel$subgroup == group], ";", fixed = TRUE)[[1]]
    stopifnot(all(genes %in% rownames(e$average)), sum(mapping$cells) == manifest$cells[i])
    p <- candidate_dotplot(e, mapping, genes, gsub("_", " ", group))
    save_plot(p$plot, file.path(folder, "DotPlot_reanalyzed_cluster_markers"),
      width = max(12, length(genes) * .22 + 3.5), height = max(4, nrow(mapping) * .29 + 1.9))
    plots[[group]] <- p$plot
    write.table(p$data, file.path(folder, "DotPlot_reanalyzed_cluster_markers_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    composition <- as.data.frame(table(cluster = e$metadata$subgroup_cluster,
      previous_candidate = e$metadata$pre_reanalysis_candidate_annotation))
    names(composition)[3] <- "cells"; composition <- composition[composition$cells > 0, ]
    composition$fraction <- composition$cells / ave(composition$cells, composition$cluster, FUN = sum)
    summaries[[group]] <- data.frame(subgroup = group, composition)
    write.table(composition, file.path(folder, "previous_candidate_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("REANALYZED_MARKER_DOTPLOT_OK ", group, "; clusters=", nrow(mapping))
  }
  grDevices::pdf(file.path(output, "All_subgroups_reanalyzed_marker_DotPlots.pdf"), width = 18, height = 7)
  for (p in plots) print(p)
  grDevices::dev.off()
  write.table(do.call(rbind, summaries), file.path(output, "all_previous_candidate_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("ALL_REANALYZED_MARKER_DOTPLOTS_OK")
}
