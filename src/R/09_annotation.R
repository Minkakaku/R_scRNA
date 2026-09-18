# Manual annotation templates and annotation application

# 输出 cluster、标签和 top marker 模板，优先保留已有非空人工标签。
write_annotation_template <- function(cluster_ids, label_values, path, cluster_col = "cluster", label_col = "cell_type",
  markers = NULL, top_n = 10L) {
  cluster_ids <- sort_cluster_ids(cluster_ids)
  label_values <- unname(label_values[cluster_ids])
  label_values[is.na(label_values)] <- ""
  if (file.exists(path)) {
    existing <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c(cluster_col, label_col) %in% colnames(existing))) {
      existing_ids <- as.character(existing[[cluster_col]])
      existing_labels <- trimws(as.character(existing[[label_col]]))
      existing_map <- stats::setNames(existing_labels, existing_ids)
      reusable <- cluster_ids %in% names(existing_map) & has_text(existing_map[cluster_ids])
      label_values[reusable] <- unname(existing_map[cluster_ids[reusable]])
    }
  }
  marker_summary <- top_marker_summary(markers, top_n = top_n)
  top_markers <- unname(marker_summary[cluster_ids])
  top_markers[is.na(top_markers)] <- ""
  template <- data.frame(cluster = cluster_ids, cell_type = label_values, top_markers = top_markers,
    stringsAsFactors = FALSE, check.names = FALSE)
  colnames(template)[1:2] <- c(cluster_col, label_col)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(template, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

# 根据全局聚类生成 cell_type 注释模板及 marker 提示。
write_cluster_annotation_template <- function(sc, markers = NULL) {
  cluster_ids <- sort_cluster_ids(sc$seurat_clusters)
  write_annotation_template(cluster_ids = cluster_ids, label_values = empty_annotation_mapping(), path = ANNOTATION$template_path,
    cluster_col = "cluster", label_col = "cell_type", markers = markers, top_n = MARKER_PARAMS$top_n_for_template)
}

# 根据子聚类生成 sub_cell_type 注释模板及 marker 提示。
write_subcluster_annotation_template <- function(obj, markers = NULL) {
  cluster_ids <- sort_cluster_ids(obj$cluster_id)
  write_annotation_template(cluster_ids = cluster_ids, label_values = empty_annotation_mapping(), path = SUBCLUSTER$template_path,
    cluster_col = "cluster_id", label_col = "sub_cell_type", markers = markers, top_n = SUBCLUSTER$marker_params$top_n_for_template)
}

# 读取非空 cluster-to-label 映射；同一 cluster 重复出现时采用最后一条。
read_template_mapping <- function(path, cluster_col, label_col) {
  if (!file.exists(path)) {
    return(character(0))
  }
  template <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c(cluster_col, label_col)
  missing_cols <- setdiff(required_cols, colnames(template))
  if (length(missing_cols) > 0) {
    stop("Annotation template missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  labels <- trimws(as.character(template[[label_col]]))
  ids <- as.character(template[[cluster_col]])
  keep <- !is.na(ids) & nzchar(ids) & !is.na(labels) & nzchar(labels)
  mapping <- stats::setNames(labels[keep], ids[keep])
  mapping[!duplicated(names(mapping), fromLast = TRUE)]
}

# 优先采用注释表中的有效映射；表为空时回退到配置映射。
resolve_annotation_mapping <- function(config_mapping, template_path, cluster_col, label_col) {
  template_mapping <- read_template_mapping(template_path, cluster_col, label_col)
  if (length(template_mapping) > 0) {
    return(template_mapping)
  }
  config_mapping <- config_mapping[!is.na(config_mapping) & nzchar(config_mapping)]
  config_mapping
}

# 检查映射覆盖全部 cluster 后写入标签；缺项时输出模板并中止。
apply_annotation <- function(obj, cluster_values, config_mapping, template_path, output_template_fn,
  cluster_col, label_col) {
  clusters <- sort_cluster_ids(cluster_values)
  mapping <- resolve_annotation_mapping(config_mapping, template_path, cluster_col, label_col)
  missing_clusters <- setdiff(clusters, names(mapping))
  if (length(missing_clusters) > 0) {
    template <- output_template_fn()
    stop(label_col, " missing clusters: ", paste(missing_clusters, collapse = ", "), ". Fill the annotation template, then rerun. Template: ",
      template, call. = FALSE)
  }
  labels <- unname(mapping[as.character(cluster_values)])
  obj[[label_col]] <- labels
  obj
}

# 按配置移除指定 cluster、保存移除统计，并可重新计算邻接图或 UMAP。
cleanup_removed_clusters <- function(obj, cleanup, context, output_dir) {
  remove_clusters <- cleanup$remove_clusters
  if (is.null(remove_clusters) || length(remove_clusters) == 0) {
    return(obj)
  }
  Idents(obj) <- "seurat_clusters"
  cluster_values <- as.character(Idents(obj))
  present_clusters <- intersect(remove_clusters, unique(cluster_values))
  if (length(present_clusters) == 0) {
    warning("No configured clusters found for removal: ", paste(remove_clusters, collapse = ", "))
    return(obj)
  }
  backup_col <- cleanup$cluster_backup_col
  if (!backup_col %in% colnames(obj@meta.data)) {
    obj[[backup_col]] <- cluster_values
  }
  removed_cells <- colnames(obj)[cluster_values %in% present_clusters]
  removed_summary <- as.data.frame(table(cluster_values[cluster_values %in% present_clusters]))
  colnames(removed_summary) <- c("removed_cluster", "cell_count")
  write.csv(removed_summary, file.path(output_dir, paste0(context, "_removed_clusters.csv")), row.names = FALSE,
    quote = FALSE)
  keep_cells <- setdiff(colnames(obj), removed_cells)
  obj <- subset(obj, cells = keep_cells)
  obj$seurat_clusters <- droplevels(factor(as.character(obj$seurat_clusters)))
  Idents(obj) <- "seurat_clusters"
  log_msg("Removed clusters from ", context, ": ", paste(present_clusters, collapse = ", "), "; removed cells=",
    length(removed_cells))
  rerun_graph_cluster <- isTRUE(cleanup$rerun_neighbors_clusters)
  rerun_umap <- isTRUE(cleanup$rerun_umap)
  if (!rerun_graph_cluster && !rerun_umap) {
    return(obj)
  }
  reduction <- choose_clustering_reduction(obj, context)
  dims <- choose_clustering_dims(obj, context, cleanup$rerun_dims)
  reduction_dims <- ncol(Embeddings(obj, reduction = reduction))
  dims <- dims[dims <= reduction_dims]
  if (length(dims) == 0) {
    stop("No valid dims left for cleanup rerun using reduction: ", reduction, call. = FALSE)
  }
  obj <- run_graph_cluster_umap(obj, reduction = reduction, dims = dims, resolution = cleanup$rerun_resolution,
    context_label = paste0(context, " cleanup"), run_graph_cluster = rerun_graph_cluster, run_umap = rerun_umap)
  store_clustering_metadata(obj, context, reduction, dims)
}

# 输出细胞类型 UMAP 及按样本、实验组拆分的展示图。
save_annotation_plots <- function(sc) {
  save_umap_plot(DimPlot(sc, group.by = "cell_type", label = TRUE, repel = TRUE), file.path(ANNOTATION_DIR,
    "UMAP_by_celltype"))
  save_umap_plot(DimPlot(sc, group.by = "cell_type", split.by = "samples", ncol = 4), file.path(ANNOTATION_DIR,
    "UMAP_split_by_sample"))
  save_umap_plot(DimPlot(sc, group.by = "cell_type", split.by = "group", ncol = 3), file.path(ANNOTATION_DIR,
    "UMAP_split_by_group"))
}

# 输出子类型 UMAP 及按实验组、样本拆分的展示图。
save_subcluster_annotation_plots <- function(obj) {
  save_umap_plot(DimPlot(obj, group.by = "sub_cell_type", label = TRUE, repel = TRUE), file.path(SUBCLUSTER_DIR,
    "UMAP_by_sub_cell_type"))
  save_umap_plot(DimPlot(obj, group.by = "sub_cell_type", split.by = "group", ncol = 3), file.path(SUBCLUSTER_DIR,
    "UMAP_split_by_group"))
  save_umap_plot(DimPlot(obj, group.by = "sub_cell_type", split.by = "samples", ncol = 3), file.path(SUBCLUSTER_DIR,
    "UMAP_split_by_samples"))
}

# 提供空映射，使人工注释默认来自项目注释表。
empty_annotation_mapping <- function() {
  character(0)
}

# 读取项目注释表中的有效 cluster-to-label 映射。
annotation_mapping_from_template <- function(template_path, cluster_col, label_col) {
  mapping <- resolve_annotation_mapping(empty_annotation_mapping(), template_path, cluster_col, label_col)
  mapping
}

# 以项目注释表完成标签赋值，缺少 cluster 标签时要求先补齐模板。
ensure_annotation_mapping <- function(obj, cluster_values, template_path, output_template_fn, cluster_col,
  label_col) {
  apply_annotation(obj = obj, cluster_values = cluster_values, config_mapping = empty_annotation_mapping(),
    template_path = template_path, output_template_fn = output_template_fn, cluster_col = cluster_col,
    label_col = label_col)
}
