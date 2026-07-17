# PCA graph clustering and UMAP

store_clustering_metadata <- function(obj, context, reduction, dims) {
  if (is.null(obj@misc$pipeline)) {
    obj@misc$pipeline <- list()
  }
  obj@misc$pipeline[[paste0(context, "_reduction")]] <- reduction
  obj@misc$pipeline[[paste0(context, "_dims")]] <- dims
  obj
}

choose_clustering_reduction <- function(obj, context, fallback = c("integrated.cca", "harmony", "pca")) {
  stored <- obj@misc$pipeline[[paste0(context, "_reduction")]]
  if (!is.null(stored) && stored %in% names(obj@reductions)) {
    return(stored)
  }
  available <- fallback[fallback %in% names(obj@reductions)]
  if (length(available) > 0) {
    return(available[1])
  }
  stop("No usable dimensional reduction found for clustering.", call. = FALSE)
}

choose_clustering_dims <- function(obj, context, requested_dims) {
  stored <- obj@misc$pipeline[[paste0(context, "_dims")]]
  dims <- if (!is.null(stored))
    stored else requested_dims
  dims[dims > 0]
}

run_graph_cluster_umap <- function(obj, reduction, dims, resolution, context_label, run_graph_cluster = TRUE,
  run_umap = TRUE) {
  if (run_graph_cluster) {
    obj <- run_step(paste0(context_label, " FindNeighbors"), FindNeighbors(obj, reduction = reduction,
      dims = dims, verbose = FALSE))
    obj <- run_step(paste0(context_label, " FindClusters"), FindClusters(obj, resolution = resolution,
      random.seed = SEED, verbose = FALSE))
  }
  if (run_umap) {
    obj <- run_step(paste0(context_label, " RunUMAP"), RunUMAP(obj, reduction = reduction, dims = dims))
  }
  obj
}

run_before_batch_umap <- function(obj, dims) {
  if (BATCH_EFFECT_UMAP$before_reduction %in% Reductions(obj)) {
    obj[[BATCH_EFFECT_UMAP$before_reduction]] <- NULL
  }
  run_step("Global RunUMAP before batch correction", RunUMAP(
    obj,
    reduction = "pca",
    dims = dims,
    reduction.name = BATCH_EFFECT_UMAP$before_reduction,
    reduction.key = BATCH_EFFECT_UMAP$before_key,
    seed.use = SEED
  ))
}

run_global_clustering <- function(sc) {
  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)
  sc <- run_step("RunPCA", RunPCA(sc, npcs = pca_settings$npcs, verbose = FALSE))
  save_plot(ElbowPlot(sc, ndims = pca_settings$npcs), file.path(GLOBAL_DIR, "PCA_elbow"))
  sc <- run_before_batch_umap(sc, pca_settings$dims)
  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Global")
  sc <- corrected$object
  reduction <- corrected$reduction
  sc <- run_graph_cluster_umap(sc, reduction = reduction, dims = pca_settings$dims, resolution = CLUSTER$resolution,
    context_label = "Global")
  sc <- store_clustering_metadata(sc, "global", reduction, pca_settings$dims)
  sc
}
