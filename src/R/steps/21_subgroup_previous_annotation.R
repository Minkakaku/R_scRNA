# 根据新 cluster 内旧标签的实际组成生成暂定映射；没有过半标签时保留 Mixed。
summarize_previous_annotation <- function(metadata, group, settings) {
  cluster <- droplevels(factor(metadata[[settings$cluster_column]]))
  original <- as.character(metadata[[settings$source_column]])
  stopifnot(!anyNA(cluster), !anyNA(original), all(nzchar(original)))
  composition <- as.data.frame(table(subgroup_cluster = cluster, original_label = original), stringsAsFactors = FALSE)
  names(composition)[3] <- "cells"
  composition <- composition[composition$cells > 0, ]
  composition$fraction <- composition$cells / ave(composition$cells, composition$subgroup_cluster, FUN = sum)
  mapping <- lapply(levels(cluster), function(id) {
    d <- composition[composition$subgroup_cluster == id, ]
    d <- d[order(-d$cells, d$original_label), ]
    n <- sum(d$cells)
    top <- d$fraction[1]
    annotation <- if (top > settings$majority_fraction) d$original_label[1] else settings$mixed_label
    status <- if (group %in% settings$limited_source_groups) "limited_source_resolution" else
      if (top >= settings$dominant_fraction) "dominant_original_label" else
      if (top > settings$majority_fraction) "partial_majority_review" else "mixed_no_majority"
    flags <- c(if (top < settings$dominant_fraction) "mixed_labels_review",
               if (n < settings$small_cluster_below) "small_cluster_review",
               if (group %in% settings$limited_source_groups) "source_has_no_finer_labels")
    data.frame(subgroup = group, subgroup_cluster = id, cells = n,
      original_label_count = nrow(d), top_original_label = d$original_label[1],
      top_original_cells = d$cells[1], top_original_fraction = top,
      tied_top_labels = paste(d$original_label[d$cells == d$cells[1]], collapse = ";"),
      second_original_label = if (nrow(d) >= 2) d$original_label[2] else "",
      second_original_fraction = if (nrow(d) >= 2) d$fraction[2] else 0,
      annotation = annotation, cluster_annotation = paste0("C", id, ": ", annotation),
      status = status, review_flags = paste(flags, collapse = ";"), stringsAsFactors = FALSE)
  })
  list(mapping = do.call(rbind, mapping), composition = data.frame(subgroup = group, composition))
}

# 添加注释字段并另存完整 QS，保留逐细胞 cell_type、新聚类和全部表达及坐标。
step_annotate_subgroups_previous <- function() {
  suppressPackageStartupMessages(library(Seurat))
  future::plan(future::sequential)
  settings <- CONFIG$subgroup_previous_annotation
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest), stringsAsFactors = FALSE)
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROJECT_DIR, settings$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  destinations <- file.path(settings$checkpoint_dir, paste0(settings$object_prefix, ".", CONFIG$project_id,
    "_", manifest$subgroup, "_previous_annotation.qs"))
  stopifnot(!any(file.exists(file.path(PROJECT_DIR, destinations))),
            !any(file.exists(paste0(file.path(PROJECT_DIR, destinations), ".partial"))))
  fingerprints <- function(x) c(RNA = digest::digest(x[["RNA"]], algo = "xxhash64"),
    reductions = digest::digest(x@reductions, algo = "xxhash64"), graphs = digest::digest(x@graphs, algo = "xxhash64"))
  summaries <- mappings <- compositions <- cell_sets <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    group_dir <- file.path(output, group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    input_sha256 <- digest::digest(file = input, algo = "sha256")
    sc <- run_step(paste(group, "Read integrated QS"), qs::qread(input, nthreads = 4))
    old_metadata <- sc[[]]
    cells <- colnames(sc); genes <- rownames(sc)
    before <- fingerprints(sc)
    stopifnot(!"subgroup_annotation" %in% names(old_metadata), ncol(sc) == manifest$cells[i])
    result <- summarize_previous_annotation(old_metadata, group, settings)
    mapping <- result$mapping
    index <- match(as.character(old_metadata[[settings$cluster_column]]), mapping$subgroup_cluster)
    stopifnot(!anyNA(index), sum(mapping$cells) == ncol(sc), nrow(mapping) == manifest$clusters[i])
    sc$subgroup_annotation <- factor(mapping$annotation[index], levels = unique(mapping$annotation))
    sc$subgroup_annotation_label <- factor(mapping$cluster_annotation[index], levels = mapping$cluster_annotation)
    sc$subgroup_annotation_top_label <- mapping$top_original_label[index]
    sc$subgroup_annotation_purity <- mapping$top_original_fraction[index]
    sc$subgroup_annotation_status <- mapping$status[index]
    sc$subgroup_annotation_review <- mapping$review_flags[index]
    Idents(sc) <- "subgroup_annotation_label"
    sc@misc$subgroup_previous_annotation <- list(stage = "provisional_previous_annotation",
      method = "within_cluster_previous_label_majority", source_column = settings$source_column,
      settings = settings, mapping = mapping, input = manifest$checkpoint[i], input_sha256 = input_sha256,
      marker_validation_performed = FALSE, original_cell_labels_preserved = TRUE)
    expected_metadata <- sc[[]]
    stopifnot(identical(expected_metadata[, names(old_metadata), drop = FALSE], old_metadata),
              identical(fingerprints(sc), before))
    destination <- file.path(PROJECT_DIR, destinations[i])
    run_step(paste(group, "Save annotated QS"), qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4))
    rm(sc)
    invisible(gc())
    check <- run_step(paste(group, "Read back annotated QS"), qs::qread(paste0(destination, ".partial"), nthreads = 4))
    stopifnot(isTRUE(validObject(check, test = TRUE)), identical(colnames(check), cells),
      identical(rownames(check), genes), identical(check[[]], expected_metadata),
      identical(fingerprints(check), before),
      identical(as.character(Idents(check)), as.character(check$subgroup_annotation_label)),
      identical(check@misc$subgroup_previous_annotation$mapping, mapping),
      file.rename(paste0(destination, ".partial"), destination))
    write.table(mapping, file.path(group_dir, "cluster_annotation_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result$composition, file.path(group_dir, "cluster_previous_label_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    for (column in c("sample", "group")) {
      d <- as.data.frame(table(cluster_annotation = check$subgroup_annotation_label, label = check[[column, drop = TRUE]]))
      names(d) <- c("cluster_annotation", column, "cells")
      write.table(d, file.path(group_dir, paste0("annotation_by_", column, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    umap <- Embeddings(check, "umap")
    connection <- gzfile(file.path(group_dir, "cell_metadata_annotated_umap.tsv.gz"), "wt")
    write.table(data.frame(cell_id = colnames(check), check[[]], UMAP_1 = umap[, 1], UMAP_2 = umap[, 2]),
                connection, sep = "\t", quote = FALSE, row.names = FALSE)
    close(connection)
    summaries[[group]] <- data.frame(subgroup = group, cells = ncol(check), clusters = nrow(mapping),
      previous_labels = length(unique(check[[settings$source_column, drop = TRUE]])),
      assigned_clusters = sum(mapping$annotation != settings$mixed_label),
      mixed_clusters = sum(mapping$annotation == settings$mixed_label),
      partial_majority_clusters = sum(mapping$status == "partial_majority_review"),
      limited_source_clusters = sum(mapping$status == "limited_source_resolution"),
      checkpoint = destinations[i], input_sha256 = input_sha256, stringsAsFactors = FALSE)
    mappings[[group]] <- mapping; compositions[[group]] <- result$composition; cell_sets[[group]] <- cells
    write.table(do.call(rbind, summaries), file.path(output, "subgroup_annotation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("PREVIOUS_ANNOTATION_QS_OK: ", group, "; clusters=", nrow(mapping),
      "; mixed=", sum(mapping$annotation == settings$mixed_label), "; cells=", ncol(check))
    rm(check)
    invisible(gc())
  }
  all_cells <- unlist(cell_sets, use.names = FALSE)
  stopifnot(!anyDuplicated(all_cells), length(all_cells) == sum(manifest$cells))
  write.table(do.call(rbind, mappings), file.path(output, "all_cluster_annotation_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, compositions), file.path(output, "all_cluster_previous_label_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  yaml::write_yaml(list(stage = "provisional_previous_annotation", settings = settings,
    total_cells = length(all_cells), clusters = sum(manifest$clusters),
    original_labels_preserved = TRUE, RNA_reductions_graphs_preserved = TRUE,
    marker_validation_performed = FALSE), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "15_annotate_subgroups_previous_sessionInfo.txt"))
  log_msg("ALL_PREVIOUS_ANNOTATIONS_OK: ", length(all_cells), " cells; ", sum(manifest$clusters), " clusters")
}

# 将旧标签投射到现有 UMAP，并并列提供 cluster 暂定映射与原标签组成。
step_plot_subgroups_previous <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$subgroup_previous_annotation
  output <- file.path(PROJECT_DIR, settings$output_dir)
  manifest <- read.delim(file.path(output, "subgroup_annotation_summary.tsv"), stringsAsFactors = FALSE)
  # 复用早期全局 UMAP 的标签配色算法，使同一旧标签继续使用原颜色。
  old_plot_data <- readRDS(file.path(PROJECT_DIR, CONFIG$object_overview$plot_data))
  shared_cells <- Reduce(intersect, lapply(old_plot_data, function(x) rownames(x$metadata)))
  shared_labels <- sort(unique(as.character(old_plot_data[[1]]$metadata[shared_cells, settings$source_column])))
  all_labels <- sort(unique(unlist(lapply(old_plot_data, function(x) as.character(x$metadata[[settings$source_column]])))))
  label_levels <- c(shared_labels, setdiff(all_labels, shared_labels))
  rgb <- grDevices::col2rgb(HONGF_NPG) / 255
  extended <- c(HONGF_NPG, grDevices::rgb(t(rgb * .65 + .35)), grDevices::rgb(t(rgb * .70)), grDevices::rgb(t(rgb * .50 + .50)))
  palette <- c(setNames(extended[seq_along(label_levels)], label_levels), setNames("#999999", settings$mixed_label))
  write.table(data.frame(label = names(palette), color = unname(palette)),
    file.path(output, "previous_label_colors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rm(old_plot_data)
  overview <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    group_dir <- file.path(output, group)
    sc <- run_step(paste(group, "Read for annotation plots"), qs::qread(file.path(PROJECT_DIR, manifest$checkpoint[i]), nthreads = 4))
    mapping <- read.delim(file.path(group_dir, "cluster_annotation_mapping.tsv"), stringsAsFactors = FALSE)
    composition <- read.delim(file.path(group_dir, "cluster_previous_label_composition.tsv"), stringsAsFactors = FALSE)
    sc[[settings$source_column]] <- factor(as.character(sc[[settings$source_column, drop = TRUE]]), levels = label_levels)
    umap <- Embeddings(sc, "umap"); xlim <- range(umap[, 1]); ylim <- range(umap[, 2])
    title <- gsub("_", " ", group)
    plot_theme <- theme_classic(base_size = 11) + theme(plot.title = element_text(hjust = .5),
      legend.text = element_text(size = 9))
    original_plot <- function(...) DimPlot(sc, reduction = "umap", group.by = settings$source_column,
      cols = palette, raster = FALSE, shuffle = TRUE, seed = SEED, ...) + labs(title = title) + plot_theme +
      guides(colour = guide_legend(override.aes = list(size = 3)))
    p <- original_plot()
    save_umap_plot(p, file.path(group_dir, "UMAP_previous_labels"), width = 11, height = 7, xlim = xlim, ylim = ylim)
    overview[[group]] <- standardize_umap_plot(p, xlim, ylim)
    for (column in c("sample", "group")) save_umap_plot(
      original_plot(split.by = column, ncol = if (column == "sample") 3 else 2),
      file.path(group_dir, paste0("UMAP_previous_labels_", column, "_split")),
      width = 17, height = if (column == "sample") 10 else 6, xlim = xlim, ylim = ylim)
    cluster_palette <- setNames(unname(palette[mapping$annotation]), mapping$cluster_annotation)
    p <- DimPlot(sc, reduction = "umap", group.by = "subgroup_annotation_label", cols = cluster_palette,
      raster = FALSE, shuffle = TRUE, seed = SEED, combine = FALSE)[[1]]
    centers <- aggregate(data.frame(x = umap[, 1], y = umap[, 2]),
      list(label = paste0("C", sc[[settings$cluster_column, drop = TRUE]])), median)
    p <- p + ggrepel::geom_text_repel(data = centers, aes(x, y, label = label),
      inherit.aes = FALSE, size = 4, seed = SEED, max.overlaps = Inf) +
      labs(title = title) + plot_theme + guides(colour = guide_legend(override.aes = list(size = 3)))
    save_umap_plot(p, file.path(group_dir, "UMAP_cluster_annotation"), width = 12, height = 7, xlim = xlim, ylim = ylim)
    composition$subgroup_cluster <- factor(composition$subgroup_cluster, levels = rev(mapping$subgroup_cluster))
    composition$original_label <- factor(composition$original_label, levels = label_levels)
    ticks <- setNames(paste0("C", mapping$subgroup_cluster, " (n=", mapping$cells, ")"), mapping$subgroup_cluster)
    p <- ggplot(composition, aes(fraction, subgroup_cluster, fill = original_label)) + geom_col(width = .8) +
      scale_fill_manual(values = palette) + scale_x_continuous(labels = function(x) paste0(round(x * 100), "%")) +
      coord_cartesian(xlim = c(0, 1)) + scale_y_discrete(labels = ticks) +
      labs(title = title, x = "Previous label fraction within cluster", y = NULL, fill = NULL) + plot_theme
    save_plot(p, file.path(group_dir, "Cluster_previous_label_composition"), width = 11, height = max(4.5, nrow(mapping) * .32))
    log_msg("PREVIOUS_ANNOTATION_PLOTS_OK: ", group, "; 5 figure pairs")
    rm(sc)
    invisible(gc())
  }
  save_plot(patchwork::wrap_plots(overview, ncol = 3), file.path(output, "All_subgroups_previous_labels"), width = 24, height = 12)
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "16_plot_subgroups_previous_sessionInfo.txt"))
  log_msg("ALL_PREVIOUS_ANNOTATION_PLOTS_OK: ", nrow(manifest) * 5 + 1, " vector PDF and PNG pairs")
}
