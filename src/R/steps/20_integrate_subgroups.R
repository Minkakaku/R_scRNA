# 从各大类 PCA 检查点继续；PC 数完全来自用户逐类选择，不重新做 PCA 或过滤细胞。
step_integrate_subgroups <- function() {
  suppressPackageStartupMessages(library(Seurat))
  future::plan(future::sequential)
  settings <- CONFIG$subgroup_harmony
  selected <- unlist(CONFIG$subgroup_pca$selected_pcs)
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest), stringsAsFactors = FALSE)
  stopifnot(setequal(names(selected), manifest$subgroup), !anyNA(selected),
            all(selected == as.integer(selected)), all(selected >= 2),
            all(selected <= CONFIG$subgroup_pca$candidate_npcs))
  output <- file.path(PROJECT_DIR, settings$output_dir)
  checkpoint_dir <- file.path(PROJECT_DIR, settings$checkpoint_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  manifest$selected_pcs <- as.integer(selected[manifest$subgroup])
  manifest$integrated_checkpoint <- file.path(settings$checkpoint_dir,
    paste0(settings$object_prefix, ".", CONFIG$project_id, "_", manifest$subgroup, "_harmony_PC",
           manifest$selected_pcs, "_res", settings$resolution, ".qs"))
  paths <- file.path(PROJECT_DIR, manifest$integrated_checkpoint)
  stopifnot(all(file.exists(file.path(PROJECT_DIR, manifest$checkpoint))),
            !any(file.exists(paths)), !any(file.exists(paste0(paths, ".partial"))))
  group_mapping <- unlist(CONFIG$group_mapping)
  stopifnot(setequal(unname(group_mapping), CONFIG$group_levels))
  layer_hashes <- function(x) vapply(c("counts", "data", "scale.data"), function(layer)
    digest::digest(LayerData(x, assay = "RNA", layer = layer), algo = "xxhash64"), character(1))
  rows <- cell_sets <- list()

  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    dims <- seq_len(manifest$selected_pcs[i])
    destination <- paths[i]
    group_dir <- file.path(output, group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    log_msg("SUBGROUP_INTEGRATION_START: ", group, "; PCs=1-", length(dims))
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    input_sha256 <- digest::digest(file = input, algo = "sha256")
    sc <- run_step(paste(group, "Read PCA QS"), qs::qread(input, nthreads = 4))
    original_metadata <- sc[[]]
    original_cells <- colnames(sc)
    original_genes <- rownames(sc)
    expression_hash <- layer_hashes(sc)
    pca_hash <- digest::digest(sc[["pca"]], algo = "xxhash64")
    hvg <- VariableFeatures(sc)
    stopifnot(identical(Reductions(sc), "pca"), length(sc@graphs) == 0L,
      ncol(sc) == manifest$cells[i], nrow(sc) == manifest$genes[i],
      all(as.character(sc[[CONFIG$subgroup_pca$parent_col, drop = TRUE]]) == group),
      !any(c(settings$original_group_col, "parent_seurat_clusters", "subgroup_cluster") %in% names(original_metadata)))
    # 原 A/B 及全局 seurat_clusters 留作历史字段，当前 group 和 seurat_clusters 用于新分析。
    sc[[settings$original_group_col]] <- original_metadata[[settings$group_col]]
    sc[[settings$group_col]] <- factor(unname(group_mapping[as.character(original_metadata[[settings$group_col]])]),
                                        levels = CONFIG$group_levels)
    sc$parent_seurat_clusters <- original_metadata$seurat_clusters
    batches <- as.character(sc[[settings$batch_col, drop = TRUE]])
    stopifnot(!anyNA(sc[[settings$group_col, drop = TRUE]]), !anyNA(batches),
      all(nzchar(batches)), length(unique(batches)) == manifest$samples[i],
      all(vapply(split(as.character(sc[[settings$group_col, drop = TRUE]]), batches),
                 function(x) length(unique(x)) == 1L, logical(1))))
    DefaultAssay(sc) <- "RNA"
    set.seed(SEED)
    # 使用已验证的 Harmony 2.x 接口；只校正选定的低维坐标，不改表达矩阵。
    sc <- run_step(paste(group, "Harmony"), harmony::RunHarmony(sc,
      group.by.vars = settings$batch_col, reduction.use = "pca", dims.use = dims,
      reduction.save = "harmony", project.dim = FALSE, theta = settings$theta,
      lambda = settings$lambda, max_iter = settings$max_iter, ncores = settings$ncores,
      early_stop = TRUE, verbose = TRUE))
    sc <- run_step(paste(group, "FindNeighbors"), FindNeighbors(sc, reduction = "harmony", dims = dims,
      k.param = settings$k_neighbors, graph.name = c("subgroup_nn", "subgroup_snn"), verbose = TRUE))
    sc <- run_step(paste(group, "FindClusters"), FindClusters(sc, graph.name = "subgroup_snn",
      cluster.name = "subgroup_cluster", resolution = settings$resolution,
      algorithm = settings$algorithm, random.seed = SEED, verbose = TRUE))
    sc$seurat_clusters <- sc$subgroup_cluster
    for (reduction in c("pca", "harmony")) {
      sc <- run_step(paste(group, "RunUMAP", reduction), RunUMAP(sc, reduction = reduction, dims = dims,
        reduction.name = if (reduction == "pca") "umap.pca" else "umap",
        reduction.key = if (reduction == "pca") "UMAPPCA_" else "UMAP_",
        n.neighbors = settings$umap_neighbors, min.dist = settings$umap_min_dist,
        metric = settings$umap_metric, seed.use = SEED, verbose = TRUE))
    }
    sc@misc$subgroup_pca$selected_pcs <- dims
    sc@misc$subgroup_pca$stage <- "continued_to_harmony"
    sc@misc$subgroup_harmony <- list(stage = "clustered_not_annotated", subgroup = group,
      selected_pcs = dims, selection_source = "user", input = manifest$checkpoint[i], input_sha256 = input_sha256,
      settings = settings, group_mapping = group_mapping, seed = SEED, expression_hash = expression_hash,
      cluster_column = "subgroup_cluster", inherited_harmony_cluster_is_parent_label = TRUE)
    unchanged_columns <- setdiff(names(original_metadata), c(settings$group_col, "seurat_clusters"))
    stopifnot(identical(sc[[]][, unchanged_columns, drop = FALSE], original_metadata[, unchanged_columns, drop = FALSE]),
      identical(sc[[]][[settings$original_group_col]], original_metadata[[settings$group_col]]),
      identical(sc[[]]$parent_seurat_clusters, original_metadata$seurat_clusters),
      identical(layer_hashes(sc), expression_hash), identical(VariableFeatures(sc), hvg),
      identical(digest::digest(sc[["pca"]], algo = "xxhash64"), pca_hash),
      ncol(Embeddings(sc, "harmony")) == length(dims), !anyNA(sc$subgroup_cluster))
    for (reduction in c("harmony", "umap.pca", "umap")) stopifnot(
      identical(rownames(Embeddings(sc, reduction)), original_cells), all(is.finite(Embeddings(sc, reduction))))
    expected_metadata <- sc[[]]
    expected_reductions <- digest::digest(sc@reductions, algo = "xxhash64")
    expected_graphs <- digest::digest(sc@graphs, algo = "xxhash64")
    run_step(paste(group, "Save integrated QS"), qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4))
    rm(sc)
    invisible(gc())
    check <- run_step(paste(group, "Read back integrated QS"), qs::qread(paste0(destination, ".partial"), nthreads = 4))
    stopifnot(isTRUE(validObject(check, test = TRUE)), identical(colnames(check), original_cells),
      identical(rownames(check), original_genes), identical(check[[]], expected_metadata),
      identical(layer_hashes(check), expression_hash), identical(VariableFeatures(check), hvg),
      identical(digest::digest(check@reductions, algo = "xxhash64"), expected_reductions),
      identical(digest::digest(check@graphs, algo = "xxhash64"), expected_graphs),
      identical(check@misc$subgroup_harmony$selected_pcs, dims),
      identical(as.character(Idents(check)), as.character(check$subgroup_cluster)),
      file.rename(paste0(destination, ".partial"), destination))
    row <- data.frame(subgroup = group, cells = ncol(check), genes = nrow(check), samples = length(unique(batches)),
      selected_pcs = length(dims), resolution = settings$resolution,
      clusters = length(unique(check$subgroup_cluster)), checkpoint = manifest$integrated_checkpoint[i],
      input_sha256 = input_sha256, stringsAsFactors = FALSE)
    rows[[group]] <- row
    cell_sets[[group]] <- colnames(check)
    yaml::write_yaml(c(list(stage = "clustered_not_annotated", selected_pc_indices = dims,
      settings = settings, group_mapping = as.list(group_mapping)), as.list(row)),
      file.path(group_dir, "integration_parameters.yaml"))
    write.table(do.call(rbind, rows), file.path(output, "subgroup_integration_summary.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("SUBGROUP_INTEGRATION_OK: ", group, "; cells=", ncol(check), "; clusters=", row$clusters,
            "; PCs=1-", length(dims))
    rm(check)
    invisible(gc())
  }
  all_cells <- unlist(cell_sets, use.names = FALSE)
  stopifnot(!anyDuplicated(all_cells), length(all_cells) == sum(manifest$cells))
  yaml::write_yaml(list(stage = "clustered_not_annotated", total_cells = length(all_cells),
    selected_pcs = as.list(selected), group_mapping = as.list(group_mapping),
    settings = settings, input_manifest = settings$input_manifest), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "13_integrate_subgroups_sessionInfo.txt"))
  log_msg("ALL_SUBGROUPS_INTEGRATED_OK: ", nrow(manifest), " groups; ", length(all_cells), " cells")
  invisible(do.call(rbind, rows))
}

# 绘图独立于计算，复用已保存坐标；各 split 面板固定相同坐标范围，所有点保持矢量。
step_plot_subgroup_integrations <- function(settings = CONFIG$subgroup_harmony,
  session_log = "14_plot_subgroups_sessionInfo.txt") {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  output <- file.path(PROJECT_DIR, settings$output_dir)
  summary <- read.delim(file.path(output, "subgroup_integration_summary.tsv"), stringsAsFactors = FALSE)
  stopifnot(setequal(summary$subgroup, names(CONFIG$subgroup_pca$selected_pcs)))
  overview <- list()
  for (i in seq_len(nrow(summary))) {
    group <- summary$subgroup[i]
    group_dir <- file.path(output, group)
    sc <- run_step(paste(group, "Read for plots"), qs::qread(file.path(PROJECT_DIR, summary$checkpoint[i]), nthreads = 4))
    umap <- Embeddings(sc, "umap")
    before <- Embeddings(sc, "umap.pca")
    xlim <- range(umap[, 1]); ylim <- range(umap[, 2])
    frame <- data.frame(cell_id = colnames(sc), sc[[]], UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
      PCA_UMAP_1 = before[, 1], PCA_UMAP_2 = before[, 2], check.names = FALSE)
    connection <- gzfile(file.path(group_dir, "cell_metadata_umap.tsv.gz"), "wt")
    write.table(frame, connection, sep = "\t", quote = FALSE, row.names = FALSE)
    close(connection)
    counts <- as.data.frame(table(subgroup_cluster = sc$subgroup_cluster))
    names(counts)[2] <- "cells"
    write.table(counts, file.path(group_dir, "cluster_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    for (column in c(settings$batch_col, settings$group_col, "cell_type", "harmony_cluster")) {
      tab <- as.data.frame(table(subgroup_cluster = sc$subgroup_cluster, label = sc[[column, drop = TRUE]]))
      names(tab) <- c("subgroup_cluster", column, "cells")
      write.table(tab, file.path(group_dir, paste0("cluster_by_", column, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    title <- gsub("_", " ", group)
    make_plot <- function(column, reduction = "umap", ...) {
      p <- npg_dimplot(sc, reduction = reduction, group.by = column, raster = FALSE,
        shuffle = TRUE, seed = SEED, ...) + labs(title = title) + theme_classic(base_size = 11) +
        theme(plot.title = element_text(hjust = .5), legend.text = element_text(size = 9))
      p + guides(colour = guide_legend(override.aes = list(size = 3), ncol = if (column == "subgroup_cluster") 2 else 1))
    }
    cluster <- make_plot("subgroup_cluster", label = TRUE, repel = TRUE)
    run_step(paste(group, "Save cluster UMAP"), save_umap_plot(cluster,
      file.path(group_dir, "UMAP_clusters"), width = 9, height = 7, xlim = xlim, ylim = ylim))
    overview[[group]] <- standardize_umap_plot(cluster + theme(legend.position = "none"), xlim, ylim)
    for (column in c(settings$batch_col, settings$group_col)) {
      run_step(paste(group, "Save UMAP", column), save_umap_plot(make_plot(column),
        file.path(group_dir, paste0("UMAP_", column)), width = 8, height = 7, xlim = xlim, ylim = ylim))
      split_plot <- make_plot("subgroup_cluster", split.by = column, ncol = if (column == settings$batch_col) 3 else 2)
      run_step(paste(group, "Save split UMAP", column), save_umap_plot(split_plot,
        file.path(group_dir, paste0("UMAP_", column, "_split")), width = 15,
        height = if (column == settings$batch_col) 10 else 6, xlim = xlim, ylim = ylim))
    }
    comparison <- patchwork::wrap_plots(
      standardize_umap_plot(make_plot(settings$batch_col, "umap.pca") + labs(title = "Before Harmony")),
      standardize_umap_plot(make_plot(settings$batch_col) + labs(title = "After Harmony")), ncol = 2, guides = "collect")
    run_step(paste(group, "Save before-after"), save_plot(comparison,
      file.path(group_dir, "UMAP_sample_before_after"), width = 14, height = 6))
    log_msg("SUBGROUP_PLOTS_OK: ", group, "; 6 vector PDF and 600 dpi PNG pairs")
    rm(sc, frame)
    invisible(gc())
  }
  run_step("Save subgroup overview", save_plot(patchwork::wrap_plots(overview, ncol = 3),
    file.path(output, "All_subgroups_UMAP_clusters"), width = 16, height = 10))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, session_log))
  log_msg("ALL_SUBGROUP_PLOTS_OK: ", nrow(summary), " groups; ", nrow(summary) * 6 + 1, " figure pairs")
}
