detect_mito_pattern <- function(features) {
  if (sum(grepl("^MT-", features)) > 0) return("^MT-")
  if (sum(grepl("^mt-", features)) > 0) return("^mt-")
  "^[Mm][Tt]-"
}

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

resolve_sample_matrix <- function(sample) {
  manifest <- discover_samples(DATA_DIR)
  hits <- manifest[manifest$sample == sample, , drop = FALSE]
  if (nrow(hits) == 0) {
    return(NULL)
  }
  hits[1, , drop = FALSE]
}

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

qc_sample_cache_path <- function(sample) {
  file.path(QS2_DIR, "qc_samples", paste0(sample, ".qs2"))
}

read_qc_sample <- function(sample) {
  path <- qc_sample_cache_path(sample)
  if (file.exists(path)) return(read_qs2(path))
  obj <- read_one_sample(sample)
  save_qs2(obj, path)
  obj
}
merge_samples <- function(obj_list) {
  obj_list <- Filter(Negate(is.null), obj_list)
  if (length(obj_list) == 0) stop("No valid sample objects.", call. = FALSE)
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
}

make_qc_plots <- function(obj) {
  plots <- list(
    violin = VlnPlot(
      obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "log10GenesPerUMI"),
      ncol = 5,
      pt.size = 0.1
    ),
    scatter_gene_umi = FeatureScatter(obj, "nCount_RNA", "nFeature_RNA"),
    scatter_mt = FeatureScatter(obj, "nCount_RNA", "percent.mt"),
    cells_per_sample = ggplot(obj@meta.data, aes(orig.ident)) +
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  if ("doublet_class" %in% colnames(obj@meta.data)) {
    plots$doublet_ratio <- ggplot(obj@meta.data, aes(orig.ident, fill = doublet_class)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = percent) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  if ("Contamination" %in% colnames(obj@meta.data)) {
    plots$contamination <- VlnPlot(obj, features = "Contamination", group.by = "orig.ident", pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  plots
}

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
