# 加载公共 YAML，并用项目 YAML 递归覆盖同名参数。
# 配置为 NULL 或空值时使用回退值。
.config_or <- function(x, fallback) if (is.null(x) || length(x) == 0) fallback else x

# 按 UTF-8 读取 YAML，空文件视为空配置列表。
.read_yaml <- function(path) .config_or(yaml::read_yaml(path, fileEncoding = "UTF-8"), list())

# 功能：加载一个项目的配置、路径和输出约定。
load_project_config <- function(project_dir) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  workspace_dir <- normalizePath(file.path(project_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  config_files <- c(
    file.path(workspace_dir, "config", "default.yaml"),
    file.path(workspace_dir, "config", "reference.yaml"),
    file.path(project_dir, "config", "config.yaml")
  )
  # 按公共默认、参考、项目配置的顺序递归覆盖；项目层具有最高优先级。
  config <- Reduce(
    function(current, path) utils::modifyList(current, .read_yaml(path), keep.null = TRUE),
    config_files,
    init = list()
  )

  config_dir <- file.path(project_dir, "config")
  data_dir <- file.path(project_dir, "data")
  data_raw_dir <- file.path(data_dir, "raw")
  data_processed_dir <- file.path(data_dir, "processed")
  data_reference_dir <- file.path(data_dir, "reference")
  results_dir <- file.path(project_dir, "results")
  report_dir <- file.path(project_dir, "docs", "reports")
  logs_dir <- file.path(project_dir, "logs")

  qc_dir <- file.path(results_dir, "01_qc")
  global_dir <- file.path(results_dir, "02_global_clustering")
  marker_dir <- file.path(results_dir, "03_markers")
  annotation_dir <- file.path(results_dir, "04_annotation")
  proportion_dir <- file.path(results_dir, "05_proportion")
  subcluster_dir <- file.path(results_dir, "06_subcluster")

  qs2_dir <- data_processed_dir
  qs2_pre_qc <- file.path(qs2_dir, "00.sc_before_QC_merged.qs2")
  qs2_qc <- file.path(qs2_dir, "01.sc_after_QC_merged.qs2")
  qs2_cluster <- file.path(qs2_dir, "02.sc_global_clustered.qs2")
  qs2_markers <- file.path(qs2_dir, "03.all_markers.qs2")
  qs2_annotated <- file.path(qs2_dir, "04.sc_annotated.qs2")
  qs2_final <- file.path(qs2_dir, "05.sc_final.qs2")

  normalize <- config$normalize
  normalize$regress_vars <- unlist(normalize$regress_vars, use.names = FALSE)
  cluster <- config$cluster
  cluster$dims <- as.integer(unlist(cluster$dims, use.names = FALSE))
  marker_params <- config$marker_params
  marker_params$assay <- normalize$assay
  group_levels <- unlist(config$group_levels, use.names = FALSE)
  proportion_comparisons <- lapply(config$proportion_comparisons, unlist, use.names = FALSE)

  annotation <- config$annotation
  annotation$template_path <- file.path(config_dir, "cluster_annotation.tsv")
  annotation$mapping_path <- file.path(annotation_dir, "cluster_celltype_mapping.tsv")
  annotation$rerun_dims <- as.integer(unlist(annotation$rerun_dims, use.names = FALSE))

  proportion <- config$proportion
  proportion$comparisons <- proportion_comparisons

  contamination_umap <- config$contamination_umap
  contamination_umap$colors <- unlist(contamination_umap$colors, use.names = FALSE)

  subcluster <- config$subcluster
  subcluster$input_qs2 <- qs2_annotated
  subcluster$target_celltype <- unlist(subcluster$target_celltype, use.names = FALSE)
  subcluster$method <- .config_or(subcluster$method, normalize$method)
  subcluster$assay <- .config_or(subcluster$assay, if (identical(subcluster$method, "SCT")) "SCT" else "RNA")
  subcluster$sct_method <- .config_or(subcluster$sct_method, normalize$sct_method)
  subcluster$sct_vst_flavor <- .config_or(subcluster$sct_vst_flavor, normalize$sct_vst_flavor)
  subcluster$sct_return_only_var_genes <- .config_or(
    subcluster$sct_return_only_var_genes,
    normalize$sct_return_only_var_genes
  )
  subcluster$sct_conserve_memory <- .config_or(subcluster$sct_conserve_memory, normalize$sct_conserve_memory)
  subcluster$regress_vars <- unlist(subcluster$regress_vars, use.names = FALSE)
  subcluster$dims <- as.integer(unlist(subcluster$dims, use.names = FALSE))
  subcluster$rerun_dims <- as.integer(unlist(subcluster$rerun_dims, use.names = FALSE))
  subcluster$template_path <- file.path(config_dir, "subcluster_annotation.tsv")
  subcluster$comparisons <- proportion_comparisons

  subcluster_name <- gsub("[^A-Za-z0-9_]+", "_", paste(subcluster$target_celltype, collapse = "_"))
  qs2_subcluster_clustered <- file.path(qs2_dir, paste0("06.", subcluster_name, ".subcluster_clustered.qs2"))
  qs2_subcluster_markers <- file.path(qs2_dir, paste0("07.", subcluster_name, ".subcluster_markers.qs2"))
  qs2_subcluster_annotated <- file.path(qs2_dir, paste0("08.", subcluster_name, ".subcluster_annotated.qs2"))

  project <- list(
    workspace_dir = workspace_dir,
    project_dir = project_dir,
    config_dir = config_dir,
    data_dir = data_dir,
    data_raw_dir = data_raw_dir,
    data_processed_dir = data_processed_dir,
    data_reference_dir = data_reference_dir,
    results_dir = results_dir,
    report_dir = report_dir,
    logs_dir = logs_dir,
    qc_dir = qc_dir,
    global_dir = global_dir,
    marker_dir = marker_dir,
    annotation_dir = annotation_dir,
    proportion_dir = proportion_dir,
    subcluster_dir = subcluster_dir,
    qs2_dir = qs2_dir
  )

  values <- list(
    CONFIG = config,
    WORKSPACE_DIR = workspace_dir,
    PROJECT_DIR = project_dir,
    CONFIG_DIR = config_dir,
    DATA_DIR = data_dir,
    DATA_RAW_DIR = data_raw_dir,
    DATA_PROCESSED_DIR = data_processed_dir,
    DATA_REFERENCE_DIR = data_reference_dir,
    RESULTS_DIR = results_dir,
    TABLE_DIR = results_dir,
    FIGURE_DIR = results_dir,
    REPORT_DIR = report_dir,
    LOGS_DIR = logs_dir,
    QC_DIR = qc_dir,
    QC_TABLE_DIR = qc_dir,
    GLOBAL_DIR = global_dir,
    MARKER_DIR = marker_dir,
    MARKER_TABLE_DIR = marker_dir,
    ANNOTATION_DIR = annotation_dir,
    ANNOTATION_TABLE_DIR = annotation_dir,
    PROPORTION_DIR = proportion_dir,
    PROPORTION_TABLE_DIR = proportion_dir,
    SUBCLUSTER_DIR = subcluster_dir,
    SUBCLUSTER_TABLE_DIR = subcluster_dir,
    QS2_DIR = qs2_dir,
    SAMPLE_SHEET_PATH = file.path(config_dir, "samples.tsv"),
    QC_FILTER_SHEET_PATH = file.path(config_dir, "qc_thresholds.tsv"),
    QS2_PRE_QC = qs2_pre_qc,
    QS2_QC = qs2_qc,
    QS2_CLUSTER = qs2_cluster,
    QS2_MARKERS = qs2_markers,
    QS2_ANNOTATED = qs2_annotated,
    QS2_FINAL = qs2_final,
    PACKAGES = unlist(config$packages, use.names = FALSE),
    OPTIONAL_PACKAGES = unlist(config$optional_packages, use.names = FALSE),
    SEED = as.integer(config$seed),
    QC_READ = config$qc_read,
    QC_DEFAULTS = config$qc_defaults,
    NORMALIZE = normalize,
    CLUSTER = cluster,
    UMAP_PLOT = config$umap_plot,
    BATCH_EFFECT_UMAP = config$batch_effect_umap,
    CONTAMINATION_UMAP = contamination_umap,
    GROUP_LEVELS = group_levels,
    MARKER_PARAMS = marker_params,
    BROAD_MARKER_GENES = lapply(config$broad_marker_genes, unlist, use.names = FALSE),
    ANNOTATION = annotation,
    PROPORTION_COMPARISONS = proportion_comparisons,
    PROPORTION = proportion,
    SUBCLUSTER = subcluster,
    SUBCLUSTER_NAME = subcluster_name,
    QS2_SUBCLUSTER_CLUSTERED = qs2_subcluster_clustered,
    QS2_SUBCLUSTER_MARKERS = qs2_subcluster_markers,
    QS2_SUBCLUSTER_ANNOTATED = qs2_subcluster_annotated,
    NORMALIZE_ROUTE = config$normalize_route,
    SCT_ROUTE = config$sct_route,
    PROJECT = project
  )
  # 将解析后的配置注册到全局环境，供现有 R 步骤函数共同使用。
  list2env(values, envir = .GlobalEnv)
  options(future.globals.maxSize = config$future_globals_max_size_gb * 1024^3)
  directories <- unique(c(unlist(project, use.names = FALSE), report_dir))
  # 已整理的项目按实际分析创建结果模块，加载配置时不重新生成空结果目录。
  if (identical(config$initialize_result_subdirs, FALSE)) {
    directories <- setdiff(directories, c(qc_dir, global_dir, marker_dir, annotation_dir, proportion_dir, subcluster_dir))
  }
  invisible(lapply(directories, dir.create,
    recursive = TRUE, showWarnings = FALSE))
  invisible(project)
}

# 功能：加载项目配置，再加载全部 R 复用函数。
load_project <- function(project_dir) {
  load_project_config(project_dir)
  sys.source(file.path(WORKSPACE_DIR, "src", "R", "load_all.R"), envir = .GlobalEnv)
  invisible(PROJECT)
}
