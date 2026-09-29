# Batch correction helpers

valid_cluster_dims <- function(obj, requested_npcs, requested_dims) {
  variable_feature_count <- length(VariableFeatures(obj))
  max_npcs <- max(1L, min(requested_npcs, variable_feature_count, ncol(obj) - 1L))
  dims <- requested_dims[requested_dims <= max_npcs]
  if (length(dims) == 0) {
    stop("No valid clustering dimensions. Check CLUSTER$npcs and CLUSTER$dims in config/config.R.",
      call. = FALSE)
  }
  list(npcs = max_npcs, dims = dims)
}

ensure_batch_column <- function(obj, batch_col, context_label) {
  if (!batch_col %in% colnames(obj@meta.data)) {
    stop(context_label, " batch column not found: ", batch_col, call. = FALSE)
  }
  batches <- unique(as.character(obj[[batch_col, drop = TRUE]]))
  batches <- batches[!is.na(batches) & nzchar(batches)]
  if (length(batches) < 2) {
    warning(context_label, " has fewer than two batches in ", batch_col, "; using PCA without correction.")
    return(FALSE)
  }
  TRUE
}

run_harmony_correction <- function(obj, settings, assay, context_label) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop(context_label, " batch_method='harmony' requires the harmony package.", call. = FALSE)
  }
  if (!ensure_batch_column(obj, settings$batch_group, context_label)) {
    return(list(object = obj, reduction = "pca"))
  }
  obj <- run_step(paste0(context_label, " RunHarmony group=", settings$batch_group), tryCatch(harmony::RunHarmony(obj,
    group.by.vars = settings$batch_group, reduction.use = "pca", dims.use = settings$dims,
    assay.use = assay, verbose = FALSE),
    error = function(e) {
      if (!grepl("assay.use", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }
      harmony::RunHarmony(obj, group.by.vars = settings$batch_group, reduction.use = "pca",
        dims.use = settings$dims, verbose = FALSE)
    }))
  list(object = obj, reduction = "harmony")
}

run_cca_correction <- function(obj, settings, normalize_settings, context_label) {
  if (!exists("IntegrateLayers", mode = "function") || !exists("CCAIntegration", mode = "function")) {
    stop(context_label, " batch_method='cca' requires Seurat v5 IntegrateLayers/CCAIntegration.",
      call. = FALSE)
  }
  if (!ensure_batch_column(obj, settings$batch_group, context_label)) {
    return(list(object = obj, reduction = "pca"))
  }
  assay <- normalize_settings$assay
  DefaultAssay(obj) <- assay
  batch_values <- as.factor(as.character(obj[[settings$batch_group, drop = TRUE]]))
  obj <- run_step(paste0(context_label, " JoinLayers before CCA"), JoinLayers(obj, assay = assay))
  obj[[assay]] <- split(obj[[assay]], f = batch_values)
  normalization <- if (identical(normalize_settings$method, "SCT"))
    "SCT" else "LogNormalize"
  obj <- run_step(paste0(context_label, " IntegrateLayers CCA group=", settings$batch_group), IntegrateLayers(object = obj,
    method = CCAIntegration, orig.reduction = "pca", new.reduction = settings$cca_reduction, normalization.method = normalization,
    verbose = FALSE))
  obj <- run_step(paste0(context_label, " JoinLayers after CCA"), JoinLayers(obj, assay = assay))
  list(object = obj, reduction = settings$cca_reduction)
}

run_batch_correction <- function(obj, settings, normalize_settings, context_label) {
  method <- settings$batch_method
  if (identical(method, "none")) {
    return(list(object = obj, reduction = "pca"))
  }
  if (identical(method, "harmony")) {
    return(run_harmony_correction(obj, settings, normalize_settings$assay, context_label))
  }
  if (identical(method, "cca")) {
    return(run_cca_correction(obj, settings, normalize_settings, context_label))
  }
  if (identical(method, "sct_integration")) {
    return(list(object = obj, reduction = "pca"))
  }
  stop("Unknown batch correction method: ", method, call. = FALSE)
}
