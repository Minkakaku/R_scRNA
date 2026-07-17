# Marker detection and marker summaries

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

sort_cluster_ids <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(!is.na(numeric_x))) {
    return(x[order(numeric_x)])
  }
  sort(x)
}

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

top_marker_summary <- function(markers, top_n = 10L) {
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers) || !"gene" %in% colnames(markers)) {
    return(character(0))
  }
  markers <- sort_markers_by_cluster(markers)
  genes_by_cluster <- split(as.character(markers$gene), as.character(markers$cluster))
  vapply(genes_by_cluster, function(genes) paste(head(unique(genes), top_n), collapse = ";"), character(1))
}

find_all_markers_wrapped <- function(obj, assay, params, label) {
  obj <- prepare_marker_assay(obj, assay)
  markers <- run_step(label, FindAllMarkers(obj, assay = assay, only.pos = params$only_pos, test.use = params$test_use,
    min.pct = params$min_pct, logfc.threshold = params$logfc_threshold, return.thresh = params$return_thresh))
  sort_markers_by_cluster(markers)
}

save_broad_marker_dotplot <- function(sc) {
  broad_markers <- lapply(BROAD_MARKER_GENES, function(x) intersect(x, rownames(sc)))
  broad_markers <- broad_markers[lengths(broad_markers) > 0]
  if (length(broad_markers) > 0) {
    save_plot(DotPlot(sc, features = broad_markers, assay = NORMALIZE$assay) + RotatedAxis(), file.path(MARKER_DIR,
      "cell_markers"), 18, 6)
  }
}
