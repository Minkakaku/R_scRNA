# 与STARTRAC的chisq分支定义一致：只比较observed/expected，不将细胞级P值当生物重复检验。
validate_roe_dotplots <- function(s) {
  folder <- file.path(PROJECT_DIR, s$evidence_dir)
  e <- readRDS(file.path(folder, 'composition_evidence.rds'))
  dots <- readRDS(file.path(folder, 'roe_dotplot_data.rds'))
  errors <- numeric(); n <- 0L
  for (key in names(dots)) {
    keys <- strsplit(key, '|', fixed = TRUE)[[1]]
    d <- dots[[key]]
    src <- e$tables[[keys[1]]]; src <- src[src$split_by == keys[2], ]
    stopifnot(identical(d$observed, src$observed), identical(d$ROE, src$ROE))
    expected_grade <- cut(d$ROE, breaks = c(-Inf, 0, .2, .8, 1, Inf),
      labels = c('0', '(0, 0.2]', '(0.2, 0.8]', '(0.8, 1]', '>1'), right = TRUE)
    stopifnot(identical(as.character(d$ROE_grade), as.character(expected_grade)))
    p <- roe_dotplot(d, levels(d$cell_type), levels(d$stratum), '')
    points <- ggplot2::ggplot_build(p)$data[[1]]
    stopifnot(nrow(points) == sum(d$ROE > 0), all(points$size > 0))
    m <- xtabs(observed ~ cell_type + stratum, d)
    checked <- suppressWarnings(chisq.test(m))
    ratio <- checked$observed / checked$expected
    ix <- cbind(match(d$cell_type, rownames(ratio)), match(d$stratum, colnames(ratio)))
    errors <- c(errors, max(abs(ratio[ix] - d$ROE)))
    stopifnot(tail(errors, 1) < 1e-12)
    n <- n + nrow(d)
  }
  result <- list(figures = length(dots), rows = n, ROE_values_unchanged = TRUE,
    dot_sizes_match_STARTRAC_classes = TRUE, STARTRAC_chisq_formula_verified = TRUE,
    max_formula_error = max(errors), p_values_used = FALSE, zero_ROE_has_no_point = TRUE)
  jsonlite::write_json(result, file.path(folder, 'roe_dotplot_validation.json'), pretty = TRUE, auto_unbox = TRUE)
  print(result); log_msg('ROE_DOTPLOT_VALIDATION_OK')
}
