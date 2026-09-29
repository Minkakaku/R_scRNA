apply_global_route <- function(route) {
  if (is.null(route)) route <- list()
  method <- route$method %||% NORMALIZE$method

  NORMALIZE$method <<- method
  NORMALIZE$assay <<- route$assay %||% if (identical(method, "SCT")) "SCT" else "RNA"
  NORMALIZE$integrated_assay <<- route$integrated_assay %||% NORMALIZE$integrated_assay
  NORMALIZE$per_sample_sct <<- route$per_sample_sct %||% NORMALIZE$per_sample_sct

  CLUSTER$batch_method <<- route$batch_method %||% CLUSTER$batch_method
  CLUSTER$batch_group <<- route$batch_group %||% CLUSTER$batch_group
  invisible(route)
}

apply_subcluster_route <- function(route) {
  if (is.null(route)) route <- list()
  method <- route$method %||% SUBCLUSTER$method

  SUBCLUSTER$method <<- method
  SUBCLUSTER$assay <<- route$assay %||% if (identical(method, "SCT")) "SCT" else "RNA"
  SUBCLUSTER$batch_method <<- route$batch_method %||% SUBCLUSTER$batch_method
  SUBCLUSTER$batch_group <<- route$batch_group %||% SUBCLUSTER$batch_group
  invisible(route)
}

run_global_normalize_route <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_global_route(route)
  step_global_clustering()
}

run_global_sct_route <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_global_route(route)
  step_global_clustering()
}

run_subcluster_normalize_route <- function(project = PROJECT, route = NORMALIZE_ROUTE) {
  apply_subcluster_route(route)
  step_subcluster_clustering()
}

run_subcluster_sct_route <- function(project = PROJECT, route = SCT_ROUTE) {
  apply_subcluster_route(route)
  step_subcluster_clustering()
}
