# Per-sample QC threshold table and filtering.

# 返回项目配置中用于初始化逐样本阈值表的 QC 默认值。
default_qc_filter_thresholds <- function() QC_DEFAULTS

# 返回阈值表必须包含的参数列名。
qc_filter_columns <- function() {
  names(default_qc_filter_thresholds())
}

# 转为数值并排除缺失值后计算分位数；无有效值时返回 NA。
metric_quantile <- function(x, prob) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  unname(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

# 按样本汇总细胞数及 QC 分位数，为人工选择过滤阈值提供参考。
sample_qc_summary <- function(sc, samples) {
  sample_values <- as.character(sc$orig.ident)
  rows <- lapply(samples, function(sample) {
    cells <- sample_values == sample
    summary_row <- data.frame(
      sample = sample,
      cells_before_qc = sum(cells),
      obs_nFeature_RNA_p01 = round(metric_quantile(sc$nFeature_RNA[cells], 0.01), 2),
      obs_nFeature_RNA_p99 = round(metric_quantile(sc$nFeature_RNA[cells], 0.99), 2),
      obs_nCount_RNA_p99 = round(metric_quantile(sc$nCount_RNA[cells], 0.99), 2),
      obs_percent_mt_p95 = round(metric_quantile(sc$percent.mt[cells], 0.95), 2),
      obs_percent_hb_p95 = round(metric_quantile(sc$percent.hb[cells], 0.95), 2),
      obs_log10GenesPerUMI_p01 = round(metric_quantile(sc$log10GenesPerUMI[cells], 0.01), 3),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if ("Contamination" %in% colnames(sc@meta.data)) {
      summary_row$obs_contamination_p50 <- round(metric_quantile(sc$Contamination[cells], 0.50), 4)
      summary_row$obs_contamination_p75 <- round(metric_quantile(sc$Contamination[cells], 0.75), 4)
      summary_row$obs_contamination_p95 <- round(metric_quantile(sc$Contamination[cells], 0.95), 4)
      summary_row$obs_contamination_p99 <- round(metric_quantile(sc$Contamination[cells], 0.99), 4)
    }

    summary_row
  })
  Reduce(rbind, rows)
}

# 解析 TSV 中常见的真假值写法；遇到非空且无法识别的值时报错。
parse_bool_value <- function(x, column_name) {
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n")] <- FALSE
  bad <- is.na(out) & !is.na(value) & nzchar(value)
  if (any(bad)) {
    stop("Invalid boolean values in ", column_name, ": ", paste(unique(value[bad]), collapse = ", "), call. = FALSE)
  }
  out
}

# 读取并校验逐样本 QC 阈值表的列、样本唯一性及参数类型。
read_qc_filter_sheet <- function(path = QC_FILTER_SHEET_PATH, required = TRUE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Sample QC threshold table not found: ", path, ". Run qc1 first, edit thresholds, then rerun qc2.", call. = FALSE)
    }
    return(NULL)
  }

  sheet <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("sample", qc_filter_columns())
  missing_cols <- setdiff(required_cols, colnames(sheet))
  if (length(missing_cols) > 0) {
    stop("Sample QC threshold table missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  sheet$sample <- trimws(as_character_column(sheet, "sample", "Sample QC threshold table"))
  if (any(!has_text(sheet$sample))) {
    stop("Sample QC threshold table contains blank sample names: ", path, call. = FALSE)
  }

  duplicate_samples <- unique(sheet$sample[duplicated(sheet$sample)])
  if (length(duplicate_samples) > 0) {
    stop("Sample QC threshold table has duplicated samples: ", paste(duplicate_samples, collapse = ", "), call. = FALSE)
  }

  numeric_cols <- setdiff(qc_filter_columns(), "keep_singlet_only")
  for (col in numeric_cols) {
    sheet[[col]] <- suppressWarnings(as.numeric(sheet[[col]]))
    if (any(is.na(sheet[[col]]))) {
      stop("Sample QC threshold table has missing or non-numeric values in ", col, ".", call. = FALSE)
    }
  }
  sheet$keep_singlet_only <- parse_bool_value(sheet$keep_singlet_only, "keep_singlet_only")
  if (any(is.na(sheet$keep_singlet_only))) {
    stop("Sample QC threshold table has missing values in keep_singlet_only.", call. = FALSE)
  }

  sheet
}

# 生成 QC 指标与阈值表，刷新统计量时保留已填写的样本阈值。
write_qc_filter_sheet <- function(sc, path = QC_FILTER_SHEET_PATH) {
  samples <- unique(as.character(sc$orig.ident))
  thresholds <- default_qc_filter_thresholds()
  summary_df <- sample_qc_summary(sc, samples)
  threshold_df <- data.frame(sample = samples, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in names(thresholds)) {
    threshold_df[[col]] <- thresholds[[col]]
  }

  existing <- read_qc_filter_sheet(path, required = FALSE)
  if (!is.null(existing)) {
    existing_sample <- existing$sample
    existing_index <- stats::setNames(seq_len(nrow(existing)), existing_sample)
    for (col in qc_filter_columns()) {
      reusable <- samples %in% existing_sample
      threshold_df[[col]][reusable] <- existing[[col]][unname(existing_index[samples[reusable]])]
    }
  }

  out <- merge(summary_df, threshold_df, by = "sample", all.x = TRUE, sort = FALSE)
  out <- out[match(samples, out$sample), , drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("Sample QC threshold table written: ", path, ". Edit threshold columns before running qc2 if needed.")
  invisible(path)
}

# 对每个样本应用其阈值，并按配置过滤双细胞和高污染细胞。
apply_qc_filter <- function(sc) {
  thresholds <- read_qc_filter_sheet(required = TRUE)
  sample_values <- as.character(sc$orig.ident)
  missing_samples <- setdiff(unique(sample_values), thresholds$sample)
  if (length(missing_samples) > 0) {
    stop(
      "Sample QC threshold table is missing samples: ",
      paste(missing_samples, collapse = ", "),
      ". Rerun qc1 to refresh ",
      QC_FILTER_SHEET_PATH,
      ".",
      call. = FALSE
    )
  }

  keep <- rep(FALSE, ncol(sc))
  for (i in seq_len(nrow(thresholds))) {
    th <- thresholds[i, , drop = FALSE]
    cells <- sample_values == th$sample[[1]]
    if (!any(cells)) next

    sample_keep <- sc$nFeature_RNA[cells] > th$nFeature_RNA_min &
      sc$nFeature_RNA[cells] < th$nFeature_RNA_max &
      sc$percent.mt[cells] < th$percent_mt_max &
      sc$percent.hb[cells] < th$percent_hb_max &
      sc$nCount_RNA[cells] < th$nCount_RNA_max &
      sc$log10GenesPerUMI[cells] > th$log10GenesPerUMI_min

    if (isTRUE(th$keep_singlet_only[[1]]) && "doublet_class" %in% colnames(sc@meta.data)) {
      sample_keep <- sample_keep & sc$doublet_class[cells] == "singlet"
    }

    if ("Contamination" %in% colnames(sc@meta.data) && is.finite(th$contamination_max[[1]])) {
      sample_keep <- sample_keep & !is.na(sc$Contamination[cells]) & sc$Contamination[cells] < th$contamination_max
    }

    sample_keep[is.na(sample_keep)] <- FALSE
    keep[cells] <- sample_keep
  }

  sc[, keep]
}
# 功能：计算并汇总 QC 过滤 指标（函数：qc_retention_table）。
qc_retention_table <- function(before, after) {
  samples <- names(before)
  kept <- after[samples]
  kept[is.na(kept)] <- 0
  data.frame(
    sample = samples,
    cells_before = as.integer(before),
    cells_after = as.integer(kept),
    retention_pct = round(100 * as.integer(kept) / as.integer(before), 1),
    row.names = NULL
  )
}
# 功能：计算并汇总 QC 过滤 指标（函数：qc_filter_pass_table）。
qc_filter_pass_table <- function(sc) {
  thresholds <- read_qc_filter_sheet()
  th <- thresholds[match(as.character(sc$orig.ident), thresholds$sample), , drop = FALSE]
  checks <- data.frame(
    nFeature_min = sc$nFeature_RNA > th$nFeature_RNA_min,
    nFeature_max = sc$nFeature_RNA < th$nFeature_RNA_max,
    percent_mt = sc$percent.mt < th$percent_mt_max,
    percent_hb = sc$percent.hb < th$percent_hb_max,
    nCount_max = sc$nCount_RNA < th$nCount_RNA_max,
    complexity = sc$log10GenesPerUMI > th$log10GenesPerUMI_min
  )
  if ("doublet_class" %in% colnames(sc@meta.data)) {
    checks$singlet <- !th$keep_singlet_only | sc$doublet_class == "singlet"
  }
  if ("Contamination" %in% colnames(sc@meta.data) && any(is.finite(th$contamination_max))) {
    checks$contamination <- !is.finite(th$contamination_max) |
      (!is.na(sc$Contamination) & sc$Contamination < th$contamination_max)
  }
  pass <- colSums(checks, na.rm = TRUE)
  data.frame(
    criterion = names(pass),
    pass_cells = as.integer(pass),
    fail_cells = ncol(sc) - as.integer(pass),
    pass_pct = round(100 * pass / ncol(sc), 1),
    row.names = NULL
  )
}
