# Sample groups and Seurat metadata.
# 判断元数据值是否同时满足非缺失和非空字符串。
has_text <- function(x) !is.na(x) & nzchar(x)

# 将指定元数据列转成字符向量，供样本和分组匹配使用。
as_character_column <- function(df, column_name, context_label = "table") {
  as.character(df[[column_name]])
}

# 读取 TSV/CSV 样本表，校验重复样本；必需表还须包含完整分组。
read_sample_sheet <- function(path = SAMPLE_SHEET_PATH, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Sample sheet not found: ", path, call. = FALSE)
    return(NULL)
  }
  reader <- if (tolower(tools::file_ext(path)) == "csv") read.csv else read.delim
  sheet <- reader(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_columns <- setdiff(c("sample", "group"), colnames(sheet))
  if (length(missing_columns)) stop("Sample sheet is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  if (any(!has_text(sheet$sample)) || anyDuplicated(sheet$sample)) {
    stop("Sample sheet contains missing or duplicated sample IDs.", call. = FALSE)
  }
  if (required && any(!has_text(sheet$group))) stop("Sample sheet contains missing groups.", call. = FALSE)
  sheet
}

# 刷新当前样本清单，并保留已有样本表中人工填写的分组。
write_sample_group_sheet <- function(sc, path = SAMPLE_SHEET_PATH, metadata_columns = NULL) {
  samples <- unique(as.character(sc$orig.ident))
  metadata <- sc[[]]
  # 现成对象可直接提供分组和可选样本信息，无需重新访问原始矩阵目录。
  if (is.null(metadata_columns) && "group" %in% names(metadata)) metadata_columns <- c(group = "group")
  if (length(metadata_columns)) {
    if (is.null(names(metadata_columns)) || !"group" %in% names(metadata_columns) ||
        !all(unname(metadata_columns) %in% names(metadata))) stop("Invalid sample metadata columns.", call. = FALSE)
    rows <- lapply(samples, function(sample) {
      hit <- metadata[as.character(metadata$orig.ident) == sample, , drop = FALSE]
      values <- lapply(unname(metadata_columns), function(column) {
        observed <- unique(as.character(hit[[column]]))
        observed <- observed[has_text(observed)]
        if (length(observed) > 1L) stop("Inconsistent ", column, " metadata for sample ", sample, call. = FALSE)
        if (length(observed)) observed[[1L]] else NA_character_
      })
      names(values) <- names(metadata_columns)
      as.data.frame(c(list(sample = sample), values), stringsAsFactors = FALSE)
    })
    manifest <- do.call(rbind, rows)
    if (any(!has_text(manifest$group))) stop("Cannot write missing sample groups.", call. = FALSE)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(manifest, path, sep = if (tolower(tools::file_ext(path)) == "csv") "," else "\t",
                row.names = FALSE, quote = tolower(tools::file_ext(path)) == "csv")
    return(invisible(path))
  }
  manifest <- discover_samples(DATA_DIR)
  manifest <- manifest[match(samples, manifest$sample), , drop = FALSE]

  existing <- read_sample_sheet(path)
  if (!is.null(existing)) {
    old_group <- setNames(existing$group, existing$sample)
    keep <- manifest$sample %in% names(old_group)
    manifest$group[keep] <- old_group[manifest$sample[keep]]
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(manifest, path, sep = if (tolower(tools::file_ext(path)) == "csv") "," else "\t",
              row.names = FALSE, quote = tolower(tools::file_ext(path)) == "csv")
  invisible(path)
}

# 从样本表生成分组映射和显示顺序；无表时以样本名作为分组。
resolve_group_mapping <- function(samples, required_sheet = FALSE) {
  sheet <- read_sample_sheet(required = required_sheet)
  if (is.null(sheet)) {
    sample_levels <- unique(as.character(samples))
    return(list(
      mapping = setNames(sample_levels, sample_levels),
      sample_levels = sample_levels,
      group_levels = sample_levels
    ))
  }

  group_levels <- c(GROUP_LEVELS, setdiff(unique(sheet$group), GROUP_LEVELS))
  list(
    mapping = setNames(sheet$group, sheet$sample),
    sample_levels = sheet$sample,
    group_levels = group_levels
  )
}

# 按 orig.ident 匹配样本表，写入有固定顺序的 samples 和 group 列。
assign_group_metadata <- function(sc, required_sheet = FALSE) {
  group_info <- resolve_group_mapping(sc$orig.ident, required_sheet)
  sc$samples <- factor(as.character(sc$orig.ident), levels = group_info$sample_levels)
  sc$group <- factor(
    unname(group_info$mapping[as.character(sc$orig.ident)]),
    levels = group_info$group_levels
  )
  if (required_sheet && (anyNA(sc$samples) || anyNA(sc$group))) {
    missing_samples <- unique(as.character(sc$orig.ident)[is.na(sc$samples) | is.na(sc$group)])
    stop("Sample sheet does not cover object samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  sc
}
