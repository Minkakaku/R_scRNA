# Marker detection and marker summaries

# 为 marker 检验准备表达 assay：SCT 尝试校正模型，其余 assay 尝试合并 layers。
prepare_marker_assay <- function(obj, assay) {
  DefaultAssay(obj) <- assay
  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "SCTAssay")) {
    obj <- tryCatch(run_step(paste0("PrepSCTFindMarkers assay=", assay), PrepSCTFindMarkers(obj,
      assay = assay, verbose = TRUE)), error = function(e) {
      warning("PrepSCTFindMarkers skipped: ", conditionMessage(e))
      obj
    })
    assay_obj <- obj[[assay]]
  }
  if (inherits(assay_obj, "SCTAssay")) {
    return(obj)
  }
  obj <- tryCatch(JoinLayers(obj, assay = assay), error = function(e) {
    warning("JoinLayers skipped for assay ", assay, ": ", conditionMessage(e))
    obj
  })
  obj
}

# 纯数字 cluster 按数值排序，其他标签按字符排序。
sort_cluster_ids <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(!is.na(numeric_x))) {
    return(x[order(numeric_x)])
  }
  sort(x)
}

# 在每个 cluster 内按 logFC 降序排列 marker，兼容两种 Seurat 列名。
sort_markers_by_cluster <- function(markers) {
  if (is.null(markers) || nrow(markers) == 0) {
    return(markers)
  }
  logfc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
  if (!is.na(logfc_col) && "cluster" %in% colnames(markers)) {
    markers <- arrange(group_by(markers, cluster), desc(.data[[logfc_col]]), .by_group = TRUE)
  }
  markers
}

# 为每个 cluster 提取去重后的前 N 个 marker，并拼成注释参考文本。
top_marker_summary <- function(markers, top_n = 10L) {
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers) || !"gene" %in% colnames(markers)) {
    return(character(0))
  }
  markers <- sort_markers_by_cluster(markers)
  genes_by_cluster <- split(as.character(markers$gene), as.character(markers$cluster))
  vapply(genes_by_cluster, function(genes) paste(head(unique(genes), top_n), collapse = ";"), character(1))
}

# 准备 assay 后按配置运行 FindAllMarkers，并整理 cluster 内的排序。
find_all_markers_wrapped <- function(obj, assay, params, label) {
  obj <- prepare_marker_assay(obj, assay)
  markers <- run_step(label, FindAllMarkers(obj, assay = assay, only.pos = params$only_pos, test.use = params$test_use,
    min.pct = params$min_pct, logfc.threshold = params$logfc_threshold, return.thresh = params$return_thresh))
  sort_markers_by_cluster(markers)
}

# 仅保留对象中存在的预设大类 marker，绘制用于人工注释的 DotPlot。
save_broad_marker_dotplot <- function(sc) {
  broad_markers <- lapply(BROAD_MARKER_GENES, function(x) intersect(x, rownames(sc)))
  broad_markers <- broad_markers[lengths(broad_markers) > 0]
  if (length(broad_markers) > 0) {
    p <- hongf_marker_dotplot(sc, features = broad_markers, assay = NORMALIZE$assay)
    save_plot(p, file.path(MARKER_DIR, "cell_markers"),
      max(11, length(unique(unlist(broad_markers))) * 0.34 + 3), 6)
  }
}
