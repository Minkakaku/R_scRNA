step_subcluster_find_markers <- function() {
  log_msg("Step 08 subcluster markers started: ", SUBCLUSTER_NAME)
  ensure_subcluster_output_dirs()

  obj <- read_qs2(QS2_SUBCLUSTER_CLUSTERED)
  Idents(obj) <- "seurat_clusters"
  obj$cluster_id <- as.character(Idents(obj))

  markers <- find_all_markers_wrapped(
    obj,
    SUBCLUSTER$assay,
    SUBCLUSTER$marker_params,
    "Subcluster FindAllMarkers"
  )

  write.csv(markers, file.path(SUBCLUSTER_TABLE_DIR, "markers_all_clusters.csv"), row.names = FALSE)
  save_qs2(markers, QS2_SUBCLUSTER_MARKERS)
  template <- write_subcluster_annotation_template(obj, markers)

  log_msg("Step 08 done. Subcluster annotation template: ", template)
}
