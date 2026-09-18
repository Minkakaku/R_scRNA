# 应用人工子类型映射和可选 cluster 清理，保存映射、图及对象。
step_subcluster_annotation <- function() {
  log_msg("Step 09 subcluster annotation started: ", SUBCLUSTER_NAME)
  ensure_subcluster_output_dirs()

  obj <- read_qs2(QS2_SUBCLUSTER_CLUSTERED)
  obj <- prepare_subcluster_input(obj)
  Idents(obj) <- "seurat_clusters"
  obj$cluster_id <- as.character(Idents(obj))
  mapping_used <- annotation_mapping_from_template(
    SUBCLUSTER$template_path,
    "cluster_id",
    "sub_cell_type"
  )

  obj <- ensure_annotation_mapping(
    obj = obj,
    cluster_values = obj$cluster_id,
    template_path = SUBCLUSTER$template_path,
    output_template_fn = function() write_subcluster_annotation_template(obj),
    cluster_col = "cluster_id",
    label_col = "sub_cell_type"
  )
  obj <- cleanup_removed_clusters(obj, SUBCLUSTER, "subcluster", SUBCLUSTER_TABLE_DIR)

  write.csv(
    data.frame(cluster_id = names(mapping_used), sub_cell_type = unname(mapping_used)),
    file.path(SUBCLUSTER_TABLE_DIR, "subcluster_annotation_mapping.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  save_subcluster_annotation_plots(obj)
  save_qs2(obj, QS2_SUBCLUSTER_ANNOTATED)

  log_msg("Step 09 done.")
}
