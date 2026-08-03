stepwise_dir <- function(kind) {
  base_dir <- if (identical(kind, "rds")) QS2_DIR else FIGURE_DIR
  dir.create(file.path(base_dir, "stepwise"), recursive = TRUE, showWarnings = FALSE)
  file.path(base_dir, "stepwise")
}

stepwise_rds <- function(name) {
  file.path(stepwise_dir("rds"), paste0(name, ".qs2"))
}

stepwise_fig <- function(name) {
  file.path(stepwise_dir("figures"), name)
}

normalize_global_01_preprocess_pca <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_global_route(route)

  sc <- read_qs2(QS2_QC)
  sc <- assign_group_metadata(sc, required_sheet = TRUE)
  sc <- run_normalize_scale(sc)

  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  sc <- run_step("Stepwise global NormalizeData RunPCA", RunPCA(sc, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(sc, stepwise_fig("normalize_global_PCA_elbow"))
  save_qs2(sc, stepwise_rds("normalize_global_01_preprocess_pca"))
}

normalize_global_02_cluster_umap <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_global_route(route)

  sc <- read_qs2(stepwise_rds("normalize_global_01_preprocess_pca"))
  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Stepwise global NormalizeData")
  sc <- corrected$object
  sc <- run_graph_cluster_umap(sc, corrected$reduction, pca_settings$dims, CLUSTER$resolution, "Stepwise global NormalizeData")
  sc <- store_clustering_metadata(sc, "global", corrected$reduction, pca_settings$dims)

  save_umap_plot(DimPlot(sc, label = TRUE, repel = TRUE), stepwise_fig("normalize_global_UMAP_clusters"))
  save_umap_plot(DimPlot(sc, group.by = "samples"), stepwise_fig("normalize_global_UMAP_by_sample"))
  save_qs2(sc, stepwise_rds("normalize_global_02_cluster_umap"))
  save_qs2(sc, QS2_CLUSTER)
}

sct_global_01_integrate_pca <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_global_route(route)

  sc <- read_qs2(QS2_QC)
  sc <- assign_group_metadata(sc, required_sheet = TRUE)
  sc <- run_normalize_scale(sc)

  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  sc <- run_step("Stepwise global SCT RunPCA", RunPCA(sc, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(sc, stepwise_fig("sct_global_PCA_elbow"))
  save_qs2(sc, stepwise_rds("sct_global_01_integrate_pca"))
}

sct_global_02_cluster_umap <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_global_route(route)

  sc <- read_qs2(stepwise_rds("sct_global_01_integrate_pca"))
  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Stepwise global SCT")
  sc <- corrected$object
  sc <- run_graph_cluster_umap(sc, corrected$reduction, pca_settings$dims, CLUSTER$resolution, "Stepwise global SCT")
  sc <- store_clustering_metadata(sc, "global", corrected$reduction, pca_settings$dims)

  save_umap_plot(DimPlot(sc, label = TRUE, repel = TRUE), stepwise_fig("sct_global_UMAP_clusters"))
  save_umap_plot(DimPlot(sc, group.by = "samples"), stepwise_fig("sct_global_UMAP_by_sample"))
  save_qs2(sc, stepwise_rds("sct_global_02_cluster_umap"))
  save_qs2(sc, QS2_CLUSTER)
}

normalize_subcluster_01_preprocess_pca <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_subcluster_route(route)

  sc <- read_qs2(SUBCLUSTER$input_qs2)
  sc <- prepare_subcluster_input(sc)
  obj <- subset_subcluster_target(sc)
  obj <- run_subcluster_normalize_scale(obj)

  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)
  obj <- run_step("Stepwise subcluster NormalizeData RunPCA", RunPCA(obj, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(obj, stepwise_fig("normalize_subcluster_PCA_elbow"))
  save_qs2(obj, stepwise_rds("normalize_subcluster_01_preprocess_pca"))
}

normalize_subcluster_02_cluster_umap <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_subcluster_route(route)

  obj <- read_qs2(stepwise_rds("normalize_subcluster_01_preprocess_pca"))
  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)
  corrected <- run_batch_correction(obj, SUBCLUSTER, SUBCLUSTER, "Stepwise subcluster NormalizeData")
  obj <- corrected$object
  obj <- run_graph_cluster_umap(obj, corrected$reduction, pca_settings$dims, SUBCLUSTER$resolution, "Stepwise subcluster NormalizeData")
  obj <- store_clustering_metadata(obj, "subcluster", corrected$reduction, pca_settings$dims)
  obj$subcluster_parent <- paste(SUBCLUSTER$target_celltype, collapse = ",")
  obj$cluster_id <- as.character(Idents(obj))

  save_umap_plot(DimPlot(obj, label = TRUE, repel = TRUE), stepwise_fig("normalize_subcluster_UMAP_clusters"))
  save_umap_plot(DimPlot(obj, group.by = "samples"), stepwise_fig("normalize_subcluster_UMAP_by_sample"))
  save_umap_plot(DimPlot(obj, group.by = "group"), stepwise_fig("normalize_subcluster_UMAP_by_group"))
  save_qs2(obj, stepwise_rds("normalize_subcluster_02_cluster_umap"))
  save_qs2(obj, QS2_SUBCLUSTER_CLUSTERED)
}

sct_subcluster_01_integrate_pca <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_subcluster_route(route)

  sc <- read_qs2(SUBCLUSTER$input_qs2)
  sc <- prepare_subcluster_input(sc)
  obj <- subset_subcluster_target(sc)
  obj <- run_subcluster_normalize_scale(obj)

  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)
  obj <- run_step("Stepwise subcluster SCT RunPCA", RunPCA(obj, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(obj, stepwise_fig("sct_subcluster_PCA_elbow"))
  save_qs2(obj, stepwise_rds("sct_subcluster_01_integrate_pca"))
}

sct_subcluster_02_cluster_umap <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_subcluster_route(route)

  obj <- read_qs2(stepwise_rds("sct_subcluster_01_integrate_pca"))
  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)
  corrected <- run_batch_correction(obj, SUBCLUSTER, SUBCLUSTER, "Stepwise subcluster SCT")
  obj <- corrected$object
  obj <- run_graph_cluster_umap(obj, corrected$reduction, pca_settings$dims, SUBCLUSTER$resolution, "Stepwise subcluster SCT")
  obj <- store_clustering_metadata(obj, "subcluster", corrected$reduction, pca_settings$dims)
  obj$subcluster_parent <- paste(SUBCLUSTER$target_celltype, collapse = ",")
  obj$cluster_id <- as.character(Idents(obj))

  save_umap_plot(DimPlot(obj, label = TRUE, repel = TRUE), stepwise_fig("sct_subcluster_UMAP_clusters"))
  save_umap_plot(DimPlot(obj, group.by = "samples"), stepwise_fig("sct_subcluster_UMAP_by_sample"))
  save_umap_plot(DimPlot(obj, group.by = "group"), stepwise_fig("sct_subcluster_UMAP_by_group"))
  save_qs2(obj, stepwise_rds("sct_subcluster_02_cluster_umap"))
  save_qs2(obj, QS2_SUBCLUSTER_CLUSTERED)
}
