# 独立会话从QS重取RNA，以base R公式复算全部分数及样本汇总。
validate_gene_module_scores <- function(s) {
  suppressPackageStartupMessages(library(Seurat))
  folder <- file.path(PROJECT_DIR, s$evidence_dir)
  e <- readRDS(file.path(folder, 'score_evidence.rds'))
  sc <- qs::qread(file.path(PROJECT_DIR, s$input), nthreads = 4)
  scored <- qs::qread(file.path(PROJECT_DIR, s$output), nthreads = 4)
  sha <- function(p) digest::digest(file = p, algo = 'sha256')
  same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-10, check.attributes = FALSE))
  stopifnot(sha(file.path(PROJECT_DIR, s$input)) == s$input_sha256,
    identical(sc[[]], e$metadata), identical(annotation_fingerprints(sc), annotation_fingerprints(scored)),
    identical(colnames(sc), colnames(scored)), identical(rownames(sc), rownames(scored)), identical(Idents(sc), Idents(scored)))
  old_names <- setdiff(names(sc[[]]), 'annotation_stage')
  stopifnot(identical(sc[[]][, old_names], scored[[]][, old_names]), identical(sc$annotation_stage, scored$pre_score_annotation_stage))
  x <- t(as.matrix(LayerData(sc, assay = 'RNA', layer = 'data')[rownames(e$expression), , drop = FALSE]))
  z <- scale(x)
  recalculated <- vapply(e$sets, function(g) rowMeans(z[, g, drop = FALSE]), numeric(nrow(x)))
  filtered <- vapply(e$filtered_sets, function(g) rowMeans(z[, g, drop = FALSE]), numeric(nrow(x)))
  stopifnot(same(recalculated, e$scores), same(filtered, e$filtered),
    same(colMeans(recalculated), rep(0, ncol(recalculated))),
    same(as.matrix(scored[[]][, paste0('figure_score_', names(e$sets))]), recalculated))
  ids <- as.character(sc[[s$cluster_column, drop = TRUE]])
  samples <- as.character(sc[[s$sample_column, drop = TRUE]])
  groups <- as.character(sc[[s$group_column, drop = TRUE]])
  cs <- read.delim(file.path(folder, 'scores_by_cluster_sample.tsv'))
  tab <- read.delim(file.path(folder, 'scores_by_cluster.tsv'))
  stopifnot(nrow(cs) == s$clusters * length(unique(samples)) * length(e$sets),
    !anyDuplicated(paste(cs$cluster, cs$sample, cs$module)))
  for (i in seq_len(nrow(cs))) {
    ix <- which(ids == as.character(cs$cluster[i]) & samples == cs$sample[i]); k <- cs$module[i]
    stopifnot(cs$cells[i] == length(ix), cs$small_sample[i] == (length(ix) < s$small_sample_cells))
    if (!length(ix)) {
      stopifnot(all(is.na(cs[i, c('mean_z', 'rare_filtered_mean_z', 'mean_log_expression', 'mean_gene_detection_fraction')]))); next
    }
    stopifnot(same(cs$mean_z[i], mean(recalculated[ix, k])), same(cs$rare_filtered_mean_z[i], mean(filtered[ix, k])),
      same(cs$mean_log_expression[i], mean(x[ix, e$sets[[k]], drop = FALSE])),
      same(cs$mean_gene_detection_fraction[i], mean(x[ix, e$sets[[k]], drop = FALSE] > 0)),
      identical(as.character(cs$group[i]), unique(groups[ix])))
  }
  for (i in seq_len(nrow(tab))) {
    a <- cs[cs$cluster == tab$cluster[i] & cs$module == tab$module[i], ]; keep <- a$cells > 0
    good <- a$cells >= s$small_sample_cells
    stopifnot(tab$cells[i] == sum(a$cells),
      same(tab$cell_mean_z[i], weighted.mean(a$mean_z[keep], a$cells[keep])),
      same(tab$sample_balanced_mean_z[i], mean(a$mean_z[keep])),
      same(tab$sample_balanced_min_cells_z[i], mean(a$mean_z[good])),
      same(tab$sample_balanced_rare_filtered_z[i], mean(a$rare_filtered_mean_z[keep])))
  }
  # DotPlot的均值、检出率及缩放值与Seurat原生逐值一致。
  native <- Seurat::DotPlot(sc, features = colnames(x), group.by = s$cluster_column, assay = 'RNA')$data
  plotted <- readRDS(file.path(folder, 'figure_plot_data.rds'))$dotplot
  ix <- match(paste(plotted$gene, plotted$cluster), paste(native$features.plot, native$id))
  stopifnot(!anyNA(ix), same(plotted$percent_expressed, native$pct.exp[ix]),
    same(plotted$average_expression, native$avg.exp[ix]), same(plotted$color_value, native$avg.exp.scaled[ix]))
  result <- list(cells = ncol(sc), clusters = length(unique(ids)), labels = length(unique(sc[[s$annotation_column, drop = TRUE]])),
    samples = length(unique(samples)), modules = length(e$sets), scores_recomputed = length(recalculated),
    maximum_score_error = max(abs(recalculated - e$scores)),
    missing_cluster_sample_pairs = length(unique(paste(cs$cluster[cs$cells == 0], cs$sample[cs$cells == 0]))),
    small_present_cluster_sample_pairs = length(unique(paste(cs$cluster[cs$cells > 0 & cs$small_sample], cs$sample[cs$cells > 0 & cs$small_sample]))),
    all_annotations_RNA_reductions_graphs_preserved = TRUE, independent_sample_summary_verified = TRUE,
    dotplot_matches_native_Seurat = TRUE, input_sha256_unchanged = TRUE)
  jsonlite::write_json(result, file.path(folder, 'independent_validation.json'), pretty = TRUE, auto_unbox = TRUE)
  print(result); log_msg('GENE_MODULE_INDEPENDENT_VALIDATION_OK')
}
