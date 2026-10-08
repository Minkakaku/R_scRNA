# 文献基因集的描述性比较；不把程序得分当成抗瘤/促瘤功能实验或处理组检验。
step_neutrophil_literature_evidence <- function() {
  suppressPackageStartupMessages(library(Seurat))
  settings <- CONFIG$literature_annotation
  folder <- file.path(PROJECT_DIR, settings$evidence_dir)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest))
  input <- file.path(PROJECT_DIR, manifest$checkpoint[manifest$subgroup == settings$focus_subgroup])
  sc <- qs::qread(input, nthreads = 4)
  programs <- read.delim(file.path(PROJECT_DIR, settings$programs), check.names = FALSE)
  sets <- setNames(strsplit(programs$genes, ';', fixed = TRUE), programs$program)
  expr <- LayerData(sc, assay = 'RNA', layer = 'data')
  genes <- intersect(unique(unlist(sets)), rownames(expr))
  x <- as.matrix(expr[genes, , drop = FALSE])
  # 每个基因在全部中性粒内对log1p归一化表达作z标准化，然后等权平均。
  # 与原文Scanpy score_genes的对照基因扣除法不同；并存raw mean和逐样本结果。
  sd_gene <- matrixStats::rowSds(x)
  names(sd_gene) <- rownames(x)
  usable <- rownames(x)[sd_gene > 0]
  z <- sweep(sweep(x[usable, , drop = FALSE], 1, rowMeans(x[usable, , drop = FALSE]), '-'), 1, sd_gene[usable], '/')
  scores <- sapply(sets, function(g) colMeans(z[intersect(g, usable), , drop = FALSE]))
  raw <- sapply(sets, function(g) colMeans(x[intersect(g, genes), , drop = FALSE]))
  coverage <- data.frame(program = names(sets), specified = lengths(sets),
    present = vapply(sets, function(g) sum(g %in% genes), integer(1)),
    variable = vapply(sets, function(g) sum(g %in% usable), integer(1)),
    missing = vapply(sets, function(g) paste(setdiff(g, genes), collapse = ';'), character(1)))
  stopifnot(all(coverage$variable >= 2), all(is.finite(scores)))
  meta <- sc[[]]; ids <- as.character(meta$subgroup_cluster)
  keys <- paste(ids, meta$sample, sep = '|')
  cluster_sample <- do.call(rbind, lapply(unique(keys), function(k) {
    ix <- which(keys == k)
    data.frame(cluster = ids[ix[1]], sample = meta$sample[ix[1]], group = as.character(meta$group[ix[1]]), cells = length(ix),
      program = colnames(scores), mean_z = colMeans(scores[ix, , drop = FALSE]),
      mean_log_expression = colMeans(raw[ix, , drop = FALSE]))
  }))
  cluster <- do.call(rbind, lapply(sort(unique(ids)), function(k) {
    ix <- which(ids == k)
    data.frame(cluster = k, cells = length(ix), program = colnames(scores),
      mean_z = colMeans(scores[ix, , drop = FALSE]), mean_log_expression = colMeans(raw[ix, , drop = FALSE]))
  }))
  # 每个样本内比较cluster顺序；不把数千细胞误当作独立生物重复。
  cluster_sample$within_sample_rank <- ave(-cluster_sample$mean_z, paste(cluster_sample$sample, cluster_sample$program), FUN = rank)
  panel <- unique(c(genes, 'Sell', 'Siglecf', 'Ccl3', 'Csf1', 'Ptgs2', 'Fgl2', 'Tgm2', 'Irf1', 'Csf3r', 'Ly6g', 'Retnlg', 'Mmp9', 'Lcn2', 'Arg1', 'Nos2', 'Cd80', 'Cd86'))
  panel <- intersect(panel, rownames(expr))
  pct <- do.call(rbind, lapply(unique(keys), function(k) {
    ix <- which(keys == k); xx <- expr[panel, ix, drop = FALSE]
    data.frame(cluster = ids[ix[1]], sample = meta$sample[ix[1]], cells = length(ix), gene = panel,
      percent = Matrix::rowMeans(xx > 0) * 100, mean_log_expression = Matrix::rowMeans(xx))
  }))
  ifn <- Matrix::colSums(expr[intersect(sets[['Interferon signaling']], rownames(expr)), , drop = FALSE] > 0) >= 3
  positive <- function(g) as.numeric(expr[g, ]) > 0
  co <- cbind(Sell_and_IFN3 = positive('Sell') & ifn,
    Ccl3_and_Csf1 = positive('Ccl3') & positive('Csf1'),
    Sell_IFN3_and_Ptgs2 = positive('Sell') & ifn & positive('Ptgs2'),
    Cd74_and_MHCII2 = positive('Cd74') & Matrix::colSums(expr[c('H2-Aa','H2-Ab1','H2-Eb1'), , drop = FALSE] > 0) >= 2)
  co_cluster <- aggregate(co * 100, list(cluster = ids), mean)
  corr <- data.frame(program = colnames(scores), rho_nFeature = apply(scores, 2, function(s) cor(s, meta$nFeature_RNA, method = 'spearman')))
  for (key in c('coverage', 'cluster', 'cluster_sample', 'pct', 'co_cluster', 'corr'))
    write.table(get(key), file.path(folder, paste0('neutrophil_', key, '.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
  saveRDS(list(scores = scores, raw_scores = raw, programs = programs, metadata = meta,
    umap = Embeddings(sc, 'umap'), input = input, source_sha256 = digest::digest(file = input, algo = 'sha256'),
    method = 'mean per-gene z of log-normalized RNA; no control genes; not exact Scanpy score_genes replication'),
    file.path(folder, 'neutrophil_program_evidence.rds'))
  print(reshape(cluster[, c('cluster', 'program', 'mean_z')], idvar = 'program', timevar = 'cluster', direction = 'wide'), row.names = FALSE)
  print(co_cluster, row.names = FALSE)
  print(coverage, row.names = FALSE)
  log_msg('NEUTROPHIL_LITERATURE_EVIDENCE_OK cells=', ncol(sc))
}
