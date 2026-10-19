# 项目配置：只放路径和需要调整的分析参数。
# workflow 脚本应先把工作目录切换到项目根目录。

PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
CONFIG_DIR <- "config"
DATA_DIR <- "data"
DATA_RAW_DIR <- DATA_DIR
DATA_PROCESSED_DIR <- file.path(DATA_DIR, "processed")
DATA_REFERENCE_DIR <- file.path(DATA_DIR, "reference")
RESULTS_DIR <- "results"
REPORT_DIR <- file.path("docs", "reports")

QC_DIR <- file.path(RESULTS_DIR, "01_qc")
GLOBAL_DIR <- file.path(RESULTS_DIR, "02_global_clustering")
MARKER_DIR <- file.path(RESULTS_DIR, "03_markers")
ANNOTATION_DIR <- file.path(RESULTS_DIR, "04_annotation")
PROPORTION_DIR <- file.path(RESULTS_DIR, "05_proportion")
SUBCLUSTER_DIR <- file.path(RESULTS_DIR, "06_subcluster")

TABLE_DIR <- RESULTS_DIR
FIGURE_DIR <- RESULTS_DIR
QC_TABLE_DIR <- QC_DIR
MARKER_TABLE_DIR <- MARKER_DIR
ANNOTATION_TABLE_DIR <- ANNOTATION_DIR
PROPORTION_TABLE_DIR <- PROPORTION_DIR
SUBCLUSTER_TABLE_DIR <- SUBCLUSTER_DIR

QS2_DIR <- file.path(DATA_PROCESSED_DIR, "rds")
SAMPLE_SHEET_PATH <- file.path(CONFIG_DIR, "sample_groups.csv")
QC_FILTER_SHEET_PATH <- file.path(CONFIG_DIR, "qc_thresholds.csv")

QS2_PRE_QC <- file.path(QS2_DIR, "00.sc_before_QC_merged.qs2")
QS2_QC <- file.path(QS2_DIR, "01.sc_after_QC_merged.qs2")
QS2_CLUSTER <- file.path(QS2_DIR, "02.sc_global_clustered.qs2")
QS2_MARKERS <- file.path(QS2_DIR, "03.all_markers.qs2")
QS2_ANNOTATED <- file.path(QS2_DIR, "04.sc_annotated.qs2")
QS2_FINAL <- file.path(QS2_DIR, "05.sc_final.qs2")

PACKAGES <- c(
  "Seurat",
  "SingleCellExperiment",
  "scDblFinder",
  "celda",
  "dplyr",
  "ggplot2",
  "scales",
  "RColorBrewer",
  "qs2"
)

OPTIONAL_PACKAGES <- c("harmony", "ggalluvial", "glmGamPoi")
SEED <- 1001L
options(future.globals.maxSize = 24 * 1024^3)

# 样本读取与 QC -----------------------------------------------------------

QC_READ <- list(
  create_min_cells = 3,
  create_min_features = 200,
  skip_sample_cells_below = 1000,
  run_doublet = TRUE,
  run_decontx = TRUE
)

# Step 01 先用默认值生成 config/qc_thresholds.csv，再逐样本人工调整。
QC_DEFAULTS <- list(
  nFeature_RNA_min = 300,
  nFeature_RNA_max = 10000,
  percent_mt_max = 10,
  percent_hb_max = 1,
  nCount_RNA_max = 60000,
  log10GenesPerUMI_min = 0.8,
  keep_singlet_only = TRUE,
  contamination_max = 0.2
)

# 归一化、批次校正与聚类 -------------------------------------------------

PIPELINE_NORMALIZATION <- "SCT" # "SCT" or "NormalizeData"
PIPELINE_ASSAY <- if (PIPELINE_NORMALIZATION == "SCT") "SCT" else "RNA"

NORMALIZE <- list(
  method = PIPELINE_NORMALIZATION,
  assay = PIPELINE_ASSAY,
  integrated_assay = "integrated",
  per_sample_sct = TRUE,
  source_assay = "RNA",
  nfeatures = 3000L,
  regress_vars = "percent.mt",
  scale_features = "variable",
  sct_method = "",
  sct_vst_flavor = "v2",
  sct_return_only_var_genes = TRUE,
  sct_conserve_memory = TRUE,
  verbose = TRUE
)

CLUSTER <- list(
  npcs = 20L,
  dims = seq_len(20L),
  resolution = 0.4,
  batch_method = "harmony", # "harmony", "cca", "sct_integration", "none"
  batch_group = "orig.ident",
  cca_reduction = "integrated.cca",
  sct_integration_reduction = "rpca"
)

UMAP_PLOT <- list(width = 12, height = 9)

BATCH_EFFECT_UMAP <- list(
  before_reduction = "umap.pca",
  before_key = "UMAPPCA_",
  after_reduction = "umap",
  sample_col = "samples",
  point_size = 1,
  split_point_size = 1.5,
  split_ncol = 5L,
  split_width = 15,
  split_height = 9
)

CONTAMINATION_UMAP <- list(
  min_cutoff = 0,
  legend_max = 1,
  point_size = 0.2,
  colors = c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019")
)

# 分组、marker 与注释 ----------------------------------------------------

GROUP_LEVELS <- c("CON", "EXP")

MARKER_PARAMS <- list(
  assay = NORMALIZE$assay,
  only_pos = TRUE,
  test_use = "wilcox",
  min_pct = 0.25,
  logfc_threshold = 0.25,
  return_thresh = 0.05,
  top_n_for_template = 10L
)

BROAD_MARKER_GENES <- list(
  Epithelial = c("EPCAM", "KRT18", "KRT19"),
  Endothelial = c("PECAM1", "CLDN5", "FLT1"),
  Fibroblast = c("COL1A1", "DCN", "PDGFRA"),
  Myeloid = c("TYROBP", "AIF1", "LST1"),
  T_cell = c("CD3D", "CD3E", "TRAC"),
  NK = c("NKG7", "KLRD1", "GNLY"),
  B_cell = c("CD79A", "MS4A1", "CD74"),
  Plasma = c("JCHAIN", "MZB1", "SDC1"),
  Mast = c("MS4A2", "KIT", "CPA3"),
  Pericyte = c("RGS5", "PDGFRB", "CSPG4"),
  SMC = c("ACTA2", "TAGLN", "MYH11")
)

ANNOTATION <- list(
  template_path = file.path(CONFIG_DIR, "cell_annotation.csv"),
  mapping_path = file.path(ANNOTATION_TABLE_DIR, "cluster_celltype_mapping.csv"),
  remove_clusters = character(0),
  rerun_neighbors_clusters = FALSE,
  rerun_umap = FALSE,
  rerun_resolution = CLUSTER$resolution,
  rerun_dims = CLUSTER$dims,
  cluster_backup_col = "global_cluster_before_annotation_cleanup"
)

PROPORTION_COMPARISONS <- list(c("CON", "EXP"))

PROPORTION <- list(
  comparisons = PROPORTION_COMPARISONS,
  make_sankey = TRUE,
  sankey_min_count = 1L
)

# 子聚类 -----------------------------------------------------------------

SUBCLUSTER_NORMALIZATION <- NORMALIZE$method
SUBCLUSTER_ASSAY <- if (SUBCLUSTER_NORMALIZATION == "SCT") "SCT" else "RNA"

SUBCLUSTER <- list(
  input_qs2 = QS2_ANNOTATED,
  parent_col = "cell_type",
  target_celltype = c("Monocyte", "Macrophage"),
  sample_col = "samples",
  group_col = "group",
  method = SUBCLUSTER_NORMALIZATION,
  assay = SUBCLUSTER_ASSAY,
  source_assay = "RNA",
  nfeatures = 3000L,
  regress_vars = "percent.mt",
  scale_features = "variable",
  sct_method = NORMALIZE$sct_method,
  sct_vst_flavor = NORMALIZE$sct_vst_flavor,
  sct_return_only_var_genes = NORMALIZE$sct_return_only_var_genes,
  sct_conserve_memory = NORMALIZE$sct_conserve_memory,
  verbose = TRUE,
  npcs = 30L,
  dims = seq_len(30L),
  resolution = 0.5,
  batch_method = "harmony",
  batch_group = "orig.ident",
  cca_reduction = "subcluster.integrated.cca",
  remove_clusters = character(0),
  rerun_neighbors_clusters = FALSE,
  rerun_umap = FALSE,
  rerun_resolution = 0.5,
  rerun_dims = seq_len(30L),
  cluster_backup_col = "subcluster_before_annotation_cleanup",
  marker_params = list(
    only_pos = TRUE,
    test_use = "wilcox",
    min_pct = 0.25,
    logfc_threshold = 0.25,
    return_thresh = 0.05,
    top_n_for_template = 10L
  ),
  template_path = file.path(CONFIG_DIR, "subcluster_annotation.csv"),
  comparisons = PROPORTION_COMPARISONS,
  make_sankey = TRUE
)

SUBCLUSTER_NAME <- gsub("[^A-Za-z0-9_]+", "_", paste(SUBCLUSTER$target_celltype, collapse = "_"))
QS2_SUBCLUSTER_CLUSTERED <- file.path(QS2_DIR, paste0("06.", SUBCLUSTER_NAME, ".subcluster_clustered.qs2"))
QS2_SUBCLUSTER_MARKERS <- file.path(QS2_DIR, paste0("07.", SUBCLUSTER_NAME, ".subcluster_markers.qs2"))
QS2_SUBCLUSTER_ANNOTATED <- file.path(QS2_DIR, paste0("08.", SUBCLUSTER_NAME, ".subcluster_annotated.qs2"))

# 可选的路线封装，供 src/steps/00_routes.R 和 stepwise 方法使用。
NORMALIZE_ROUTE <- list(
  method = "NormalizeData",
  assay = "RNA",
  batch_method = CLUSTER$batch_method,
  batch_group = CLUSTER$batch_group
)

SCT_ROUTE <- list(
  method = "SCT",
  assay = "SCT",
  integrated_assay = NORMALIZE$integrated_assay,
  per_sample_sct = TRUE,
  batch_method = "harmony",
  batch_group = CLUSTER$batch_group
)

PROJECT <- list(
  project_dir = PROJECT_DIR,
  config_dir = CONFIG_DIR,
  data_dir = DATA_DIR,
  data_raw_dir = DATA_RAW_DIR,
  data_processed_dir = DATA_PROCESSED_DIR,
  data_reference_dir = DATA_REFERENCE_DIR,
  results_dir = RESULTS_DIR,
  table_dir = TABLE_DIR,
  figure_dir = FIGURE_DIR,
  report_dir = REPORT_DIR,
  qc_dir = QC_DIR,
  qc_table_dir = QC_TABLE_DIR,
  global_dir = GLOBAL_DIR,
  marker_dir = MARKER_DIR,
  marker_table_dir = MARKER_TABLE_DIR,
  annotation_dir = ANNOTATION_DIR,
  annotation_table_dir = ANNOTATION_TABLE_DIR,
  proportion_dir = PROPORTION_DIR,
  proportion_table_dir = PROPORTION_TABLE_DIR,
  subcluster_dir = SUBCLUSTER_DIR,
  subcluster_table_dir = SUBCLUSTER_TABLE_DIR,
  qs2_dir = QS2_DIR
)

PROJECT_DIRS <- unique(c(
  CONFIG_DIR,
  DATA_RAW_DIR,
  DATA_PROCESSED_DIR,
  DATA_REFERENCE_DIR,
  REPORT_DIR,
  QC_DIR,
  GLOBAL_DIR,
  MARKER_DIR,
  ANNOTATION_DIR,
  PROPORTION_DIR,
  SUBCLUSTER_DIR,
  QS2_DIR,
  "logs"
))
invisible(lapply(PROJECT_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
