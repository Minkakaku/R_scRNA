# 独立重计数核验每个类别、每个样本和分组；检查ROE守恒关系和输入缓存哈希。
validate_composition_and_roe <- function(s) {
  e <- readRDS(file.path(PROJECT_DIR, s$evidence_dir, 'composition_evidence.rds'))
  meta <- e$metadata
  stopifnot(nrow(meta) == s$expected_cells, !anyDuplicated(meta$cell_id))
  for (i in seq_len(nrow(e$input_hashes))) stopifnot(
    digest::digest(file = e$input_hashes$path[i], algo = 'sha256') == e$input_hashes$sha256[i])
  same <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10))
  rows_checked <- 0L
  for (scope in names(e$tables)) {
    main <- scope == 'Main_lineages'
    m <- if (main) meta else meta[meta$analysis_container == scope, ]
    category <- if (main) s$major_column else s$subtype_column
    out <- file.path(PROJECT_DIR, if (main) s$main_results else file.path(s$subgroup_results, scope))
    formal <- read.delim(file.path(out, 'Composition_ROE.tsv'))
    for (i in seq_len(nrow(formal))) {
      d <- formal[i, ]; split <- s[[paste0(d$split_by, '_column')]]
      in_type <- as.character(m[[category]]) == d$cell_type
      in_stratum <- as.character(m[[split]]) == d$stratum
      observed <- sum(in_type & in_stratum)
      expected <- (sum(in_type) / nrow(m)) * sum(in_stratum)
      stopifnot(d$observed == observed, d$type_total == sum(in_type), d$stratum_total == sum(in_stratum),
        same(d$expected, expected), same(d$fraction, observed / sum(in_stratum)), same(d$ROE, observed / expected))
    }
    rows_checked <- rows_checked + nrow(formal)
    # 每个柱子的比例和为1；按类别总体占比加权的ROE均为1。
    for (kind in unique(formal$split_by)) {
      d <- formal[formal$split_by == kind, ]
      stopifnot(sum(d$observed) == nrow(m))
      for (a in unique(d$stratum)) {
        z <- d[d$stratum == a, ]
        stopifnot(same(sum(z$fraction), 1), same(sum(z$ROE * z$type_total / nrow(m)), 1))
      }
    }
    group_table <- formal[formal$split_by == 'group', ]; sample_table <- formal[formal$split_by == 'sample', ]
    for (i in seq_len(nrow(group_table))) {
      d <- group_table[i, ]
      samples <- e$sample_groups[[s$sample_column]][e$sample_groups[[s$group_column]] == d$stratum]
      stopifnot(d$observed == sum(sample_table$observed[sample_table$cell_type == d$cell_type & sample_table$stratum %in% samples]))
    }
  }
  result <- list(cells = nrow(meta), main_lineages = length(unique(meta[[s$major_column]])),
    analysis_containers = length(unique(meta$analysis_container)), subtype_labels = length(unique(meta[[s$subtype_column]])),
    sample_group_mapping = e$sample_groups, rows_independently_verified = rows_checked,
    input_cache_hashes_verified = nrow(e$input_hashes), all_counts_fractions_ROE_verified = TRUE,
    group_counts_equal_sample_sums = TRUE, no_duplicate_or_dropped_cells = TRUE)
  jsonlite::write_json(result, file.path(PROJECT_DIR, s$evidence_dir, 'independent_validation.json'), pretty = TRUE, auto_unbox = TRUE)
  print(result); log_msg('COMPOSITION_INDEPENDENT_VALIDATION_OK')
}
