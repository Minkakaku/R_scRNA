# 应用完整的人工注释和可选 cluster 清理，保存使用的映射、图和注释对象。
step_cell_annotation <- function() {
  log_msg("Step 05 cell annotation started")

  sc <- read_qs2(QS2_CLUSTER)
  sc <- assign_group_metadata(sc, required_sheet = TRUE)
  mapping_used <- annotation_mapping_from_template(
    ANNOTATION$template_path,
    "cluster",
    "cell_type"
  )

  sc <- ensure_annotation_mapping(
    obj = sc,
    cluster_values = as.character(sc$seurat_clusters),
    template_path = ANNOTATION$template_path,
    output_template_fn = function() write_cluster_annotation_template(sc),
    cluster_col = "cluster",
    label_col = "cell_type"
  )
  sc <- cleanup_removed_clusters(sc, ANNOTATION, "global", ANNOTATION_TABLE_DIR)

  mapping <- data.frame(cluster = names(mapping_used), cell_type = unname(mapping_used))
  write.table(mapping, ANNOTATION$mapping_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(
    as.data.frame(table(cluster = sc$seurat_clusters, cell_type = sc$cell_type)),
    file.path(ANNOTATION_TABLE_DIR, "cluster_celltype_count.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  save_annotation_plots(sc)
  save_qs2(sc, QS2_ANNOTATED)

  log_msg("Step 05 done.")
}
