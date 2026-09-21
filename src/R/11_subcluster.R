# Cell-type subclustering

# 统一亚群对象中的样本列、分组映射和显示顺序，检查分组是否缺失。
prepare_subcluster_input <- function(sc) {
  if (!SUBCLUSTER$sample_col %in% colnames(sc@meta.data)) {
    sc$samples <- as.character(sc$orig.ident)
  } else {
    sc$samples <- as.character(sc[[SUBCLUSTER$sample_col, drop = TRUE]])
  }
  group_info <- resolve_group_mapping(sc$samples, required_sheet = TRUE)
  sc$group <- unname(group_info$mapping[sc$samples])
  if (!is.null(group_info$sample_levels) && length(group_info$sample_levels) > 0) {
    sc$samples <- factor(sc$samples, levels = group_info$sample_levels)
  }
  if (!is.null(group_info$group_levels) && length(group_info$group_levels) > 0) {
    sc$group <- factor(sc$group, levels = group_info$group_levels)
  }
  if (anyNA(sc$group)) {
    stop("Some subcluster samples do not have group values in ", SAMPLE_SHEET_PATH, call. = FALSE)
  }
  sc
}

# 按父级注释提取配置指定的大类细胞，找不到目标时中止。
subset_subcluster_target <- function(sc) {
  if (!SUBCLUSTER$parent_col %in% colnames(sc@meta.data)) {
    stop("Missing parent annotation column: ", SUBCLUSTER$parent_col, call. = FALSE)
  }
  values <- as.character(sc[[SUBCLUSTER$parent_col, drop = TRUE]])
  cells <- colnames(sc)[!is.na(values) & values %in% SUBCLUSTER$target_celltype]
  if (length(cells) == 0) {
    stop("No cells found for SUBCLUSTER_TARGET=", paste(SUBCLUSTER$target_celltype, collapse = ","),
      call. = FALSE)
  }
  subset(sc, cells = cells)
}

# 在目标亚群内重新归一化、选择高变基因并缩放，支持 SCT 或 LogNormalize。
run_subcluster_normalize_scale <- function(obj) {
  if (identical(SUBCLUSTER$method, "SCT")) {
    return(run_sctransform(obj, SUBCLUSTER, "Subcluster"))
  }
  if (!identical(SUBCLUSTER$method, "NormalizeData")) {
    stop("SUBCLUSTER$method must be NormalizeData or SCT.", call. = FALSE)
  }
  DefaultAssay(obj) <- SUBCLUSTER$assay
  obj <- run_step("Subcluster NormalizeData", NormalizeData(obj, verbose = TRUE))
  obj <- run_step(paste0("Subcluster FindVariableFeatures nfeatures=", SUBCLUSTER$nfeatures), FindVariableFeatures(obj,
    nfeatures = SUBCLUSTER$nfeatures, verbose = TRUE))
  regress_vars <- intersect(SUBCLUSTER$regress_vars, colnames(obj@meta.data))
  scale_mode <- tolower(SUBCLUSTER$scale_features)
  if (scale_mode %in% c("variable", "hvg", "variable_features")) {
    scale_features <- VariableFeatures(obj)
  } else if (scale_mode %in% c("all", "all_genes", "all_features")) {
    scale_features <- rownames(obj[[SUBCLUSTER$assay]])
  } else {
    stop("SUBCLUSTER$scale_features must be 'variable' or 'all'.", call. = FALSE)
  }
  if (length(scale_features) == 0) {
    stop("No features found before subcluster ScaleData.", call. = FALSE)
  }
  if (length(regress_vars) > 0) {
    obj <- run_step(paste0("Subcluster ScaleData method=", SUBCLUSTER$method, " assay=", SUBCLUSTER$assay,
      " scale_features=", scale_mode, " nfeatures=", length(scale_features), " regress=", format_list_value(regress_vars)),
      ScaleData(object = obj, features = scale_features, vars.to.regress = regress_vars, verbose = TRUE))
  } else {
    obj <- run_step(paste0("Subcluster ScaleData method=", SUBCLUSTER$method, " assay=", SUBCLUSTER$assay,
      " scale_features=", scale_mode, " nfeatures=", length(scale_features), " regress=<none>"),
      ScaleData(object = obj, features = scale_features, verbose = TRUE))
  }
  obj
}

# 对亚群重新执行 PCA、批次校正、聚类和 UMAP，并保留父类及子 cluster 编号。
run_subcluster_clustering <- function(obj) {
  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)
  obj <- run_step("Subcluster RunPCA", RunPCA(obj, npcs = pca_settings$npcs, verbose = FALSE))
  save_pca_elbow(obj, file.path(SUBCLUSTER_DIR, "PCA_elbow"))
  save_umap_plot(npg_dimplot(run_step("Subcluster raw RunUMAP", RunUMAP(obj, reduction = "pca", dims = pca_settings$dims)),
    group.by = "samples", pt.size = 0.2), file.path(SUBCLUSTER_DIR, "UMAP_raw_by_sample"))
  corrected <- run_batch_correction(obj, SUBCLUSTER, SUBCLUSTER, "Subcluster")
  obj <- corrected$object
  reduction <- corrected$reduction
  obj <- run_graph_cluster_umap(obj, reduction = reduction, dims = pca_settings$dims, resolution = SUBCLUSTER$resolution,
    context_label = "Subcluster")
  obj <- store_clustering_metadata(obj, "subcluster", reduction, pca_settings$dims)
  obj$subcluster_parent <- paste(SUBCLUSTER$target_celltype, collapse = ",")
  obj$cluster_id <- as.character(Idents(obj))
  obj
}
