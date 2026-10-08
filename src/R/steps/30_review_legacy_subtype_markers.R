# 以旧名称对应的命名基因和谱系组合核对当前cluster；不按旧标签多数票赋值。
step_review_legacy_subtype_markers <- function() {
  s <- CONFIG$legacy_marker_annotation
  panel <- read.delim(file.path(PROJECT_DIR, s$legacy_markers), encoding = 'UTF-8')
  output <- file.path(PROJECT_DIR, s$evidence_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  all_values <- list()
  for (group in unique(panel$subgroup)) {
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    genes <- unique(unlist(strsplit(panel$genes[panel$subgroup == group], ';', fixed = TRUE)))
    stopifnot(all(genes %in% rownames(e$average)))
    d <- expand.grid(gene = genes, cluster = colnames(e$average), stringsAsFactors = FALSE)
    d$subgroup <- group
    d$average <- as.vector(e$average[genes, , drop = FALSE])
    d$percent <- as.vector(e$percent[genes, , drop = FALSE])
    all_values[[group]] <- d
    cat('\n', group, 'DETECTION_PERCENT\n')
    print(round(e$percent[genes, , drop = FALSE], 1))
    cat('\n', group, 'LOG1P_MEAN_EXPRESSION\n')
    print(round(log1p(e$average[genes, , drop = FALSE]), 2))
  }
  write.table(do.call(rbind, all_values), file.path(output, 'legacy_marker_expression.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  cat('LEGACY_MARKER_EVIDENCE_OK\n')
}
