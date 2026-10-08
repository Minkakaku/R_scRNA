# 基于现有RNA归一化层计算描述性模块分数，并保留逐细胞和逐样本证据。
step_score_gene_modules <- function(s) {
  suppressPackageStartupMessages(library(Seurat))
  stopifnot(isTRUE(s$run_scoring), !is.null(s$confirmed_on))
  folder <- file.path(PROJECT_DIR, s$evidence_dir)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  input <- file.path(PROJECT_DIR, s$input)
  sha <- function(p) digest::digest(file = p, algo = 'sha256')
  stopifnot(sha(input) == s$input_sha256,
    sha(file.path(PROJECT_DIR, s$mapping)) == s$mapping_sha256,
    sha(file.path(PROJECT_DIR, s$user_reference)) == s$reference_sha256)
  sc <- qs::qread(input, nthreads = 4)
  old <- sc[[]]; before <- annotation_fingerprints(sc)
  ids <- as.character(old[[s$cluster_column]])
  samples <- as.character(old[[s$sample_column]])
  groups <- as.character(old[[s$group_column]])
  clusters <- sort(unique(ids)); sample_levels <- sort(unique(samples))
  map <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8', check.names = FALSE)
  map <- map[map$subgroup == s$scope, ]
  map <- map[match(clusters, as.character(map$subgroup_cluster)), ]
  map$cells <- as.integer(table(factor(ids, levels = clusters)))
  map$display <- paste0(map$annotation, ' (n=', format(map$cells, big.mark = ',', trim = TRUE), ')')
  stopifnot(ncol(sc) == s$expected_cells, length(clusters) == s$clusters,
    identical(as.character(old[[s$annotation_column]]), map$annotation[match(ids, clusters)]))
  programs <- read.delim(file.path(PROJECT_DIR, s$modules), encoding = 'UTF-8', check.names = FALSE)
  sets <- setNames(strsplit(programs$genes, ';', fixed = TRUE), programs$module)
  expr <- LayerData(sc, assay = 'RNA', layer = 'data')
  genes <- unique(unlist(sets)); stopifnot(all(genes %in% rownames(expr)))
  x <- as.matrix(expr[genes, , drop = FALSE])
  mu <- rowMeans(x); sig <- matrixStats::rowSds(x); names(sig) <- genes
  detection <- rowMeans(x > 0) * 100
  usable <- genes[sig > 0]
  z <- sweep(sweep(x[usable, , drop = FALSE], 1, mu[usable], '-'), 1, sig[usable], '/')
  filtered_sets <- lapply(sets, function(g) intersect(g, usable[detection[usable] >= s$rare_detection_percent]))
  stopifnot(all(lengths(filtered_sets) > 0), all(lengths(lapply(sets, intersect, usable)) > 0))
  scores <- sapply(sets, function(g) colMeans(z[intersect(g, usable), , drop = FALSE]))
  filtered <- sapply(filtered_sets, function(g) colMeans(z[g, , drop = FALSE]))
  raw <- sapply(sets, function(g) colMeans(x[g, , drop = FALSE]))
  detected <- sapply(sets, function(g) colMeans(x[g, , drop = FALSE] > 0))
  stopifnot(all(is.finite(scores)), all(is.finite(filtered)))
  coverage <- do.call(rbind, lapply(names(sets), function(k) data.frame(
    module = k, specified_genes = length(sets[[k]]), variable_genes = sum(sets[[k]] %in% usable),
    retained_after_rare_filter = length(filtered_sets[[k]]),
    excluded_rare = paste(setdiff(sets[[k]], filtered_sets[[k]]), collapse = ';'))))
  gene_stats <- data.frame(gene = genes, percent_expressed = detection, mean_log_expression = mu,
    sd_log_expression = sig, retained_after_rare_filter = detection >= s$rare_detection_percent)
  # 缺失的cluster×sample保留零细胞和NA分数，不填充成零分。
  sample_summary <- do.call(rbind, lapply(clusters, function(k) do.call(rbind, lapply(sample_levels, function(a) {
    ix <- which(ids == k & samples == a)
    means <- function(m) if (length(ix)) colMeans(m[ix, , drop = FALSE]) else rep(NA_real_, ncol(m))
    data.frame(cluster = k, annotation = map$annotation[match(k, clusters)], sample = a,
      group = unique(groups[samples == a]), cells = length(ix), small_sample = length(ix) < s$small_sample_cells,
      module = names(sets), mean_z = means(scores), rare_filtered_mean_z = means(filtered),
      mean_log_expression = means(raw), mean_gene_detection_fraction = means(detected))
  }))))
  summary <- do.call(rbind, lapply(clusters, function(k) {
    ix <- which(ids == k)
    do.call(rbind, lapply(names(sets), function(m) {
      d <- sample_summary[sample_summary$cluster == k & sample_summary$module == m, ]
      good <- d$cells >= s$small_sample_cells
      data.frame(cluster = k, annotation = map$annotation[match(k, clusters)], cells = length(ix), module = m,
        genes = paste(sets[[m]], collapse = ';'), samples_present = sum(d$cells > 0),
        samples_at_least_min_cells = sum(good),
        cell_mean_z = mean(scores[ix, m]), sample_balanced_mean_z = mean(d$mean_z, na.rm = TRUE),
        sample_balanced_min_cells_z = if (any(good)) mean(d$mean_z[good]) else NA_real_,
        cell_rare_filtered_z = mean(filtered[ix, m]), sample_balanced_rare_filtered_z = mean(d$rare_filtered_mean_z, na.rm = TRUE),
        mean_log_expression = mean(raw[ix, m]), mean_gene_detection_fraction = mean(detected[ix, m]))
    }))
  }))
  corr <- data.frame(module = names(sets),
    rho_nFeature = apply(scores, 2, function(v) cor(v, old$nFeature_RNA, method = 'spearman')),
    rho_nCount = apply(scores, 2, function(v) cor(v, old$nCount_RNA, method = 'spearman')),
    rho_rare_filtered = vapply(names(sets), function(k) cor(scores[, k], filtered[, k], method = 'spearman'), numeric(1)))
  # 逐基因点图使用Seurat定义的线性表达均值，而不是直接平均log值。
  average <- sapply(clusters, function(k) rowMeans(expm1(x[, ids == k, drop = FALSE])))
  percent <- sapply(clusters, function(k) rowMeans(x[, ids == k, drop = FALSE] > 0) * 100)
  cell <- data.frame(cell_id = colnames(sc), cluster = ids, annotation = as.character(old[[s$annotation_column]]),
    sample = samples, group = groups)
  for (pair in list(c('score_', 'scores'), c('rare_filtered_', 'filtered'), c('mean_log_', 'raw'), c('detected_fraction_', 'detected'))) {
    values <- as.data.frame(get(pair[2])); names(values) <- paste0(pair[1], names(values)); cell <- cbind(cell, values)
  }
  write_tsv <- function(d, name) write.table(d, file.path(folder, name), sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  write_tsv(sample_summary, 'scores_by_cluster_sample.tsv')
  write_tsv(summary, 'scores_by_cluster.tsv')
  write_tsv(gene_stats, 'gene_detection.tsv'); write_tsv(coverage, 'module_coverage.tsv'); write_tsv(corr, 'score_correlations.tsv')
  con <- gzfile(file.path(folder, 'scores_by_cell.tsv.gz'), 'wt')
  write.table(cell, con, sep = '\t', quote = FALSE, row.names = FALSE); close(con)
  new_meta <- as.data.frame(scores); names(new_meta) <- paste0('figure_score_', names(new_meta))
  new_filtered <- as.data.frame(filtered); names(new_filtered) <- paste0('figure_score_rare_filtered_', names(new_filtered))
  sc <- AddMetaData(sc, cbind(new_meta, new_filtered))
  sc$pre_score_annotation_stage <- sc$annotation_stage
  sc$annotation_stage <- paste0('user_confirmed_', s$confirmed_on)
  sc@misc$figure_scoring <- list(settings = s, modules = programs, coverage = coverage,
    annotation_confirmed = TRUE, reference_source_identified = FALSE,
    module_sha256 = sha(file.path(PROJECT_DIR, s$modules)), input_sha256 = sha(input),
    method = 'Per-gene sample SD z of existing log-normalized RNA across input cells, then equal gene mean; no control genes',
    interpretation = 'RNA expression programs; not N1/N2 probabilities, ROS/NET activity, or proven tumor direction')
  stopifnot(identical(sc[[]][, setdiff(names(old), 'annotation_stage')], old[, setdiff(names(old), 'annotation_stage')]))
  save_annotation_verified(sc, file.path(PROJECT_DIR, s$output), before)
  plot_object <- CreateSeuratObject(LayerData(sc, assay = 'RNA', layer = 'counts')[1, , drop = FALSE], meta.data = sc[[]])
  plot_object[['umap']] <- sc[['umap']]
  saveRDS(list(scores = scores, filtered = filtered, raw = raw, detection = detected, expression = x,
    gene_mean = mu, gene_sd = sig, programs = programs, sets = sets, filtered_sets = filtered_sets,
    mapping = map, metadata = old, plot_object = plot_object, average = average, percent = percent,
    summary = summary, sample_summary = sample_summary, settings = s, fingerprints = before),
    file.path(folder, 'score_evidence.rds'))
  stopifnot(sha(input) == s$input_sha256)
  write_tsv(data.frame(input_sha256 = sha(input), output_sha256 = sha(file.path(PROJECT_DIR, s$output)),
    cells = ncol(sc), genes = nrow(sc), clusters = length(clusters), labels = length(unique(map$annotation)),
    samples = length(sample_levels), preserved_RNA_reductions_graphs = TRUE, full_readback_verified = TRUE), 'object_validation.tsv')
  write.table(summary, file.path(PROJECT_DIR, s$results_dir, 'Figure_scores_summary.tsv'),
    sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  capture.output(sessionInfo(), file = file.path(folder, 'sessionInfo.txt'))
  print(summary[, c('annotation', 'cells', 'module', 'cell_mean_z', 'sample_balanced_mean_z', 'sample_balanced_rare_filtered_z')], row.names = FALSE)
  log_msg('GENE_MODULE_SCORING_OK cells=', ncol(sc), '; clusters=', length(clusters))
}
