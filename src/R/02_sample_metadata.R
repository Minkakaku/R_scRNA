# Sample groups and Seurat metadata.
# 判断元数据值是否同时满足非缺失和非空字符串。
has_text <- function(x) !is.na(x) & nzchar(x)

# 将指定元数据列转成字符向量，供样本和分组匹配使用。
as_character_column <- function(df, column_name, context_label = "table") {
  as.character(df[[column_name]])
}

# 读取样本 TSV；必需表缺失时报错，可选表缺失时返回 NULL。
read_sample_sheet <- function(path = SAMPLE_SHEET_PATH, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Sample sheet not found: ", path, call. = FALSE)
    return(NULL)
  }
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# 刷新当前样本清单，并保留已有样本表中人工填写的分组。
write_sample_group_sheet <- function(sc, path = SAMPLE_SHEET_PATH) {
  samples <- unique(as.character(sc$orig.ident))
  manifest <- discover_samples(DATA_DIR)
  manifest <- manifest[match(samples, manifest$sample), , drop = FALSE]

  existing <- read_sample_sheet(path)
  if (!is.null(existing)) {
    old_group <- setNames(existing$group, existing$sample)
    keep <- manifest$sample %in% names(old_group)
    manifest$group[keep] <- old_group[manifest$sample[keep]]
  }

  write.table(manifest, path, sep = "\t", row.names = FALSE, quote = FALSE)
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
  sc
}
