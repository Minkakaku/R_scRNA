# Sample groups and Seurat metadata.
has_text <- function(x) !is.na(x) & nzchar(x)

as_character_column <- function(df, column_name, context_label = "table") {
  as.character(df[[column_name]])
}

read_sample_sheet <- function(path = SAMPLE_SHEET_PATH, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Sample group CSV not found: ", path, call. = FALSE)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

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

  write.csv(manifest, path, row.names = FALSE, quote = FALSE)
  invisible(path)
}

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

assign_group_metadata <- function(sc, required_sheet = FALSE) {
  group_info <- resolve_group_mapping(sc$orig.ident, required_sheet)
  sc$samples <- factor(as.character(sc$orig.ident), levels = group_info$sample_levels)
  sc$group <- factor(
    unname(group_info$mapping[as.character(sc$orig.ident)]),
    levels = group_info$group_levels
  )
  sc
}
