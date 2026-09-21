# 创建当前亚群分析所需的图表输出目录。
ensure_subcluster_output_dirs <- function() {
  dir.create(SUBCLUSTER_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(SUBCLUSTER_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# 提取目标大类细胞并重新归一化和聚类，保存亚群 UMAP 与对象。
step_subcluster_clustering <- function() {
  log_msg("Step 07 subcluster clustering started: ", SUBCLUSTER_NAME)
  ensure_subcluster_output_dirs()

  sc <- read_qs2(SUBCLUSTER$input_qs2)
  sc <- prepare_subcluster_input(sc)
  obj <- subset_subcluster_target(sc)

  obj <- run_subcluster_normalize_scale(obj)
  obj <- run_subcluster_clustering(obj)

  save_umap_plot(
    npg_dimplot(obj, label = TRUE, repel = TRUE),
    file.path(SUBCLUSTER_DIR, "UMAP_clusters")
  )
  save_umap_plot(
    npg_dimplot(obj, group.by = "samples"),
    file.path(SUBCLUSTER_DIR, "UMAP_by_sample")
  )
  save_umap_plot(
    npg_dimplot(obj, group.by = "group"),
    file.path(SUBCLUSTER_DIR, "UMAP_by_group")
  )

  save_qs2(obj, QS2_SUBCLUSTER_CLUSTERED)
  log_msg("Step 07 done. Check ", file.path(SUBCLUSTER_DIR, "UMAP_clusters.pdf"))
}
