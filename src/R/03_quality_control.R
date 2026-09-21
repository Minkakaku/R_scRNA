# 根据实际基因名选择人或小鼠常见的线粒体基因前缀。
detect_mito_pattern <- function(features) {
  if (sum(grepl("^MT-", features)) > 0) return("^MT-")
  if (sum(grepl("^mt-", features)) > 0) return("^mt-")
  "^[Mm][Tt]-"
}

# 构造血红蛋白基因的完整名称匹配式，兼容人和小鼠命名。
detect_hb_pattern <- function() {
  hb_genes <- c(
    "HBA[0-9]*",
    "HBB",
    "HBD",
    "HBE[0-9]*",
    "HBG[0-9]*",
    "HBM",
    "HBQ[0-9]*",
    "HBZ",
    "Hba-[[:alnum:]]+",
    "Hbb-[[:alnum:]]+",
    "Hba[0-9]*",
    "Hbb[0-9]*",
    "Hbd",
    "Hbe[0-9]*",
    "Hbg[0-9]*",
    "Hbm",
    "Hbq[0-9]*",
    "Hbz"
  )
  paste0("^(", paste(hb_genes, collapse = "|"), ")$")
}

# 基于 RNA assay 计算线粒体、血红蛋白比例及 log10GenesPerUMI 复杂度。
add_qc_metrics <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = detect_mito_pattern(rownames(obj)))
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = detect_hb_pattern())
  obj$log10GenesPerUMI <- ifelse(
    obj$nCount_RNA > 1,
    log10(obj$nFeature_RNA) / log10(obj$nCount_RNA),
    NA_real_
  )
  obj
}

# 在发现的矩阵清单中定位一个样本；没有匹配时返回 NULL。
resolve_sample_matrix <- function(sample) {
  manifest <- discover_samples(DATA_DIR)
  hits <- manifest[manifest$sample == sample, , drop = FALSE]
  if (nrow(hits) == 0) {
    return(NULL)
  }
  hits[1, , drop = FALSE]
}

# 读取单样本计数并创建 Seurat 对象，按配置附加双细胞和污染评估。
read_one_sample <- function(sample) {
  matrix_record <- resolve_sample_matrix(sample)
  if (is.null(matrix_record)) {
    warning("Skip sample without matrix: ", sample)
    return(NULL)
  }

  matrix_path <- matrix_record$matrix_dir[[1]]
  counts <- run_step(paste("Read10X", sample), {
    if (identical(matrix_record$matrix_type[[1]], "h5")) {
      Read10X_h5(matrix_path)
    } else {
      Read10X(matrix_path)
    }
  })
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = QC_READ$create_min_cells,
    min.features = QC_READ$create_min_features
  )
  obj <- add_qc_metrics(obj)

  if (ncol(obj) < QC_READ$skip_sample_cells_below) {
    warning("Skip small sample: ", sample)
    return(NULL)
  }

  if (QC_READ$run_doublet) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = NORMALIZE$nfeatures, verbose = FALSE)
    sce <- scDblFinder(as.SingleCellExperiment(obj))
    obj$doublet_class <- sce$scDblFinder.class
    obj$doublet_score <- sce$scDblFinder.score
  }

  if (QC_READ$run_decontx) {
    sce <- as.SingleCellExperiment(obj)
    decont <- decontX(counts(sce))
    obj$Contamination <- decont$contamination
  }

  obj
}

# 为每个样本生成独立的 QC 对象缓存路径。
qc_sample_cache_path <- function(sample) {
  file.path(QS2_DIR, "qc_samples", paste0(sample, ".qs2"))
}

# 优先复用已有 QC 缓存；否则读取样本、计算指标并保存缓存。
read_qc_sample <- function(sample) {
  path <- qc_sample_cache_path(sample)
  if (file.exists(path)) return(read_qs2(path))
  obj <- read_one_sample(sample)
  save_qs2(obj, path)
  obj
}
# 排除无效样本对象，并为多样本细胞编号添加样本前缀后合并。
merge_samples <- function(obj_list) {
  obj_list <- Filter(Negate(is.null), obj_list)
  if (length(obj_list) == 0) stop("No valid sample objects.", call. = FALSE)
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
}

# 双细胞类别采用固定语义配色；过滤后只剩 singlet 时仍为蓝色。
plot_qc_doublet_ratio <- function(metadata, sample_col = "orig.ident", class_col = "doublet_class") {
  stopifnot(all(as.character(metadata[[class_col]]) %in% c("singlet", "doublet")))
  ggplot(metadata, aes(.data[[sample_col]], fill = .data[[class_col]])) +
    geom_bar(position = "fill", show.legend = TRUE) +
    scale_fill_manual(values = c(singlet = "#4DBBD5", doublet = "#E64B35"),
      limits = c("singlet", "doublet"), drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Cell proportion", fill = NULL) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 构造基因数、UMI、线粒体等 QC 图，并按已有列附加双细胞与污染图。
make_qc_plots <- function(obj) {
  sample_levels <- if (is.factor(obj$orig.ident)) levels(droplevels(obj$orig.ident)) else sort(unique(obj$orig.ident))
  obj$orig.ident <- factor(obj$orig.ident, levels = sample_levels)
  palette <- npg_palette(sample_levels)
  plots <- list(
    violin = VlnPlot(
      obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "log10GenesPerUMI"),
      ncol = 5,
      pt.size = 0, group.by = "orig.ident", cols = unname(palette), raster = FALSE
    ),
    scatter_gene_umi = FeatureScatter(obj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", cols = unname(palette), raster = FALSE),
    scatter_mt = FeatureScatter(obj, "nCount_RNA", "percent.mt", group.by = "orig.ident", cols = unname(palette), raster = FALSE),
    cells_per_sample = ggplot(obj@meta.data, aes(orig.ident, fill = orig.ident)) +
      geom_bar(width = 0.75) + scale_fill_manual(values = palette, guide = "none") +
      labs(x = NULL, y = "Cell count") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  if ("doublet_class" %in% colnames(obj@meta.data)) {
    plots$doublet_ratio <- plot_qc_doublet_ratio(obj@meta.data)
  }

  if ("Contamination" %in% colnames(obj@meta.data)) {
    plots$contamination <- VlnPlot(obj, features = "Contamination", group.by = "orig.ident", pt.size = 0,
      cols = unname(palette), raster = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  plots
}

# 按过滤前后阶段标签逐一保存 QC 图，保持 PNG/PDF 同名配对。
save_qc_plots <- function(obj, stage) {
  plots <- make_qc_plots(obj)
  for (nm in names(plots)) {
    save_plot(
      plots[[nm]],
      file.path(QC_DIR, paste0("QC_", nm, "_", stage)),
      width = ifelse(nm == "violin", 16, 10),
      height = 5
    )
  }
}
