# 返回 Seurat 对象已有的 assay 名称，供后续方法检查输入。
seurat_assay_names <- function(obj) {
  names(obj@assays)
}

# 按配置的源 assay、回归变量和拟合方法执行 SCTransform。
run_sctransform <- function(obj, settings, context_label) {
  if (!settings$source_assay %in% seurat_assay_names(obj)) {
    stop(context_label, " SCT requires source assay: ", settings$source_assay, call. = FALSE)
  }

  DefaultAssay(obj) <- settings$source_assay
  regress_vars <- intersect(settings$regress_vars, colnames(obj@meta.data))
  has_sct_method <- nzchar(settings$sct_method)

  if (has_sct_method &&
      identical(settings$sct_method, "glmGamPoi") &&
      !requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop(context_label, " SCT method glmGamPoi requires the glmGamPoi package.", call. = FALSE)
  }

  sct_label <- paste0(
    context_label,
    " SCTransform source_assay=",
    settings$source_assay,
    " assay=",
    settings$assay,
    " nfeatures=",
    settings$nfeatures,
    " regress=",
    format_list_value(regress_vars),
    " flavor=",
    settings$sct_vst_flavor,
    " method=",
    ifelse(has_sct_method, settings$sct_method, "<default>")
  )

  # 只在配置有效时传入回归变量和 SCT 拟合方法，保留 Seurat 的其余默认行为。
  run_call <- function(vst_flavor = settings$sct_vst_flavor, conserve_memory = settings$sct_conserve_memory) {
    if (length(regress_vars) > 0 && has_sct_method) {
      return(SCTransform(
        object = obj,
        assay = settings$source_assay,
        new.assay.name = settings$assay,
        variable.features.n = settings$nfeatures,
        vst.flavor = vst_flavor,
        return.only.var.genes = settings$sct_return_only_var_genes,
        conserve.memory = conserve_memory,
        vars.to.regress = regress_vars,
        method = settings$sct_method,
        verbose = settings$verbose
      ))
    }

    if (length(regress_vars) > 0) {
      return(SCTransform(
        object = obj,
        assay = settings$source_assay,
        new.assay.name = settings$assay,
        variable.features.n = settings$nfeatures,
        vst.flavor = vst_flavor,
        return.only.var.genes = settings$sct_return_only_var_genes,
        conserve.memory = conserve_memory,
        vars.to.regress = regress_vars,
        verbose = settings$verbose
      ))
    }

    if (has_sct_method) {
      return(SCTransform(
        object = obj,
        assay = settings$source_assay,
        new.assay.name = settings$assay,
        variable.features.n = settings$nfeatures,
        vst.flavor = vst_flavor,
        return.only.var.genes = settings$sct_return_only_var_genes,
        conserve.memory = conserve_memory,
        method = settings$sct_method,
        verbose = settings$verbose
      ))
    }

    SCTransform(
      object = obj,
      assay = settings$source_assay,
      new.assay.name = settings$assay,
      variable.features.n = settings$nfeatures,
      vst.flavor = vst_flavor,
      return.only.var.genes = settings$sct_return_only_var_genes,
      conserve.memory = conserve_memory,
      verbose = settings$verbose
    )
  }

  obj <- run_step(sct_label, run_call())
  DefaultAssay(obj) <- settings$assay
  obj
}

# 尝试合并指定 assay 的分层数据；接口不支持时保留原对象并提示。
join_assay_layers_if_possible <- function(obj, assay, context_label) {
  if (!exists("JoinLayers", mode = "function") || !assay %in% seurat_assay_names(obj)) {
    return(obj)
  }

  tryCatch(
    JoinLayers(obj, assay = assay),
    error = function(e) {
      warning(context_label, " JoinLayers skipped for assay ", assay, ": ", conditionMessage(e))
      obj
    }
  )
}

# 逐样本 SCT 后选择共同特征、寻找锚点并生成整合 assay。
run_per_sample_sct_integration <- function(obj, settings, cluster_settings, context_label) {
  required_fns <- c(
    "SplitObject",
    "SelectIntegrationFeatures",
    "PrepSCTIntegration",
    "FindIntegrationAnchors",
    "IntegrateData"
  )
  missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function")]
  if (length(missing_fns) > 0) {
    stop(
      context_label,
      " per-sample SCT integration requires missing Seurat functions: ",
      paste(missing_fns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!cluster_settings$batch_group %in% colnames(obj@meta.data)) {
    stop(context_label, " batch column not found: ", cluster_settings$batch_group, call. = FALSE)
  }

  obj <- join_assay_layers_if_possible(obj, settings$source_assay, context_label)
  sample_values <- as.character(obj[[cluster_settings$batch_group, drop = TRUE]])
  sample_values <- sample_values[!is.na(sample_values) & nzchar(sample_values)]
  if (length(unique(sample_values)) < 2) {
    warning(context_label, " has fewer than two samples in ", cluster_settings$batch_group, "; running regular SCT.")
    return(run_sctransform(obj, settings, context_label))
  }

  obj_list <- run_step(
    paste0(context_label, " SplitObject by ", cluster_settings$batch_group),
    SplitObject(obj, split.by = cluster_settings$batch_group)
  )
  rm(obj)
  invisible(gc())

  sample_names <- names(obj_list)
  obj_list <- setNames(lapply(sample_names, function(sample_name) {
    run_sctransform(obj_list[[sample_name]], settings, paste0(context_label, " sample=", sample_name))
  }), sample_names)

  features <- run_step(
    paste0(context_label, " SelectIntegrationFeatures nfeatures=", settings$nfeatures),
    SelectIntegrationFeatures(object.list = obj_list, nfeatures = settings$nfeatures)
  )

  obj_list <- run_step(
    paste0(context_label, " PrepSCTIntegration"),
    PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = settings$verbose)
  )

  integration_reduction <- tolower(cluster_settings$sct_integration_reduction)
  if (integration_reduction %in% c("rpca", "pca")) {
    sample_names <- names(obj_list)
    obj_list <- setNames(lapply(sample_names, function(sample_name) {
      run_step(
        paste0(context_label, " sample=", sample_name, " RunPCA for SCT anchors"),
        RunPCA(
          obj_list[[sample_name]],
          features = features,
          npcs = cluster_settings$npcs,
          verbose = FALSE
        )
      )
    }), sample_names)
  }

  anchors <- run_step(
    paste0(
      context_label,
      " FindIntegrationAnchors SCT reduction=",
      integration_reduction,
      " dims=",
      format_list_value(cluster_settings$dims)
    ),
    FindIntegrationAnchors(
      object.list = obj_list,
      normalization.method = "SCT",
      anchor.features = features,
      reduction = integration_reduction,
      dims = cluster_settings$dims,
      verbose = settings$verbose
    )
  )

  integrated <- run_step(
    paste0(context_label, " IntegrateData SCT assay=", settings$integrated_assay),
    IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = cluster_settings$dims,
      new.assay.name = settings$integrated_assay,
      verbose = settings$verbose
    )
  )

  DefaultAssay(integrated) <- settings$integrated_assay
  VariableFeatures(integrated) <- features
  integrated@misc$pipeline$per_sample_sct_integration <- list(
    batch_group = cluster_settings$batch_group,
    source_assay = settings$source_assay,
    sct_assay = settings$assay,
    integrated_assay = settings$integrated_assay,
    anchor_features = length(features),
    anchor_reduction = integration_reduction,
    dims = cluster_settings$dims
  )
  integrated
}

# 逐样本 SCT 后合并对象并保留共同特征，供后续降维与批次校正使用。
run_per_sample_sct_merge <- function(obj, settings, cluster_settings, context_label) {
  obj <- join_assay_layers_if_possible(obj, settings$source_assay, context_label)
  obj_list <- run_step(
    paste0(context_label, " SplitObject by ", cluster_settings$batch_group),
    SplitObject(obj, split.by = cluster_settings$batch_group)
  )
  rm(obj)
  invisible(gc())

  sample_names <- names(obj_list)
  obj_list <- setNames(lapply(sample_names, function(sample_name) {
    run_sctransform(obj_list[[sample_name]], settings, paste0(context_label, " sample=", sample_name))
  }), sample_names)

  features <- run_step(
    paste0(context_label, " SelectIntegrationFeatures nfeatures=", settings$nfeatures),
    SelectIntegrationFeatures(object.list = obj_list, nfeatures = settings$nfeatures)
  )
  obj_list <- run_step(
    paste0(context_label, " PrepSCTIntegration for merge"),
    PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = settings$verbose)
  )
  merged <- run_step(
    paste0(context_label, " Merge per-sample SCT objects"),
    merge(obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
  )
  DefaultAssay(merged) <- settings$assay
  VariableFeatures(merged) <- features
  merged@misc$pipeline$per_sample_sct_harmony <- list(
    batch_group = cluster_settings$batch_group,
    source_assay = settings$source_assay,
    sct_assay = settings$assay,
    features = length(features)
  )
  merged
}

# 按配置选择 SCT 或 LogNormalize，并明确缩放基因范围和回归变量。
run_normalize_scale <- function(sc) {
  if (identical(NORMALIZE$method, "SCT")) {
    if (isTRUE(NORMALIZE$per_sample_sct)) {
      if (identical(CLUSTER$batch_method, "sct_integration")) {
        return(run_per_sample_sct_integration(sc, NORMALIZE, CLUSTER, "Global"))
      }
      return(run_per_sample_sct_merge(sc, NORMALIZE, CLUSTER, "Global"))
    }
    return(run_sctransform(sc, NORMALIZE, "Global"))
  }

  if (!identical(NORMALIZE$method, "NormalizeData")) {
    stop("NORMALIZE$method must be NormalizeData or SCT.", call. = FALSE)
  }

  DefaultAssay(sc) <- NORMALIZE$assay
  sc <- run_step("NormalizeData", NormalizeData(sc, verbose = NORMALIZE$verbose))
  sc <- run_step(
    paste0("FindVariableFeatures nfeatures=", NORMALIZE$nfeatures),
    FindVariableFeatures(sc, nfeatures = NORMALIZE$nfeatures, verbose = NORMALIZE$verbose)
  )
  regress_vars <- intersect(NORMALIZE$regress_vars, colnames(sc@meta.data))
  scale_mode <- tolower(NORMALIZE$scale_features)
  if (scale_mode %in% c("variable", "hvg", "variable_features")) {
    scale_features <- VariableFeatures(sc)
  } else if (scale_mode %in% c("all", "all_genes", "all_features")) {
    scale_features <- rownames(sc[[NORMALIZE$assay]])
  } else {
    stop("NORMALIZE$scale_features must be 'variable' or 'all'.", call. = FALSE)
  }

  if (length(scale_features) == 0) {
    stop("No features found before ScaleData.", call. = FALSE)
  }

  if (length(regress_vars) > 0) {
    sc <- run_step(
      paste0(
        "ScaleData method=",
        NORMALIZE$method,
        " assay=",
        NORMALIZE$assay,
        " scale_features=",
        scale_mode,
        " nfeatures=",
        length(scale_features),
        " regress=",
        format_list_value(regress_vars)
      ),
      ScaleData(
        object = sc,
        features = scale_features,
        vars.to.regress = regress_vars,
        verbose = NORMALIZE$verbose
      )
    )
  } else {
    sc <- run_step(
      paste0(
        "ScaleData method=",
      NORMALIZE$method,
      " assay=",
      NORMALIZE$assay,
      " scale_features=",
      scale_mode,
      " nfeatures=",
      length(scale_features),
        " regress=<none>"
      ),
      ScaleData(
        object = sc,
        features = scale_features,
        verbose = NORMALIZE$verbose
      )
    )
  }
  sc
}
