step_find_markers <- function() {
  log_msg("Step 04 find markers started")

  sc <- read_qs2(QS2_CLUSTER)
  DefaultAssay(sc) <- MARKER_PARAMS$assay
  Idents(sc) <- "seurat_clusters"

  markers <- find_all_markers_wrapped(sc, MARKER_PARAMS$assay, MARKER_PARAMS, "FindAllMarkers")
  write.csv(markers, file.path(MARKER_TABLE_DIR, "AllClusters_markers.csv"), row.names = FALSE)
  save_qs2(markers, QS2_MARKERS)

  template <- write_cluster_annotation_template(sc, markers)
  save_broad_marker_dotplot(sc)

  log_msg("Step 04 done. Annotation template: ", template)
}
