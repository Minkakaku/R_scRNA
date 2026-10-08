# 复核现有注释：表达阈值仅用于描述共表达，不是过滤、doublet判定或功能验证。
step_review_annotation_evidence <- function(settings = CONFIG$annotation_review) {
  suppressPackageStartupMessages(library(Seurat))
  s <- settings
  manifest <- read.delim(file.path(PROJECT_DIR, s$input_manifest))
  panels <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  rules <- read.delim(file.path(PROJECT_DIR, s$rules), check.names = FALSE)
  output <- file.path(PROJECT_DIR, s$evidence_dir)
  stopifnot(!dir.exists(output))
  dir.create(output, recursive = TRUE)
  marker <- sample_marker <- coexpression <- qc <- missing_markers <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    stopifnot(identical(digest::digest(file = input, algo = 'sha256'), manifest$output_sha256[i]))
    sc <- qs::qread(input, nthreads = 4)
    genes <- strsplit(panels$genes[panels$subgroup == group], ';', fixed = TRUE)[[1]]
    absent <- setdiff(genes, rownames(sc))
    if (!isTRUE(s$allow_missing_markers)) stopifnot(!length(absent))
    if (length(absent)) missing_markers[[group]] <- data.frame(subgroup = group, gene = absent,
      interpretation = 'not_in_input_gene_table_not_a_measured_negative')
    genes <- intersect(genes, rownames(sc))
    x <- LayerData(sc, assay = 'RNA', layer = 'data')[genes, , drop = FALSE]
    linear <- x; linear@x <- expm1(linear@x)
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    ids <- as.character(sc$subgroup_cluster)
    # 规则是项目配置内人工审核的R逻辑表达式，基因变量均表示RNA data > 0。
    rr <- rules[rules$subgroup == group, ]
    stopifnot(all(unique(unlist(lapply(rr$expression, function(z) all.vars(str2lang(z))))) %in% genes))
    detected <- setNames(lapply(genes, function(g) as.numeric(x[g, ]) > 0), genes)
    flags <- setNames(lapply(rr$expression, function(z) eval(str2lang(z), envir = detected)), rr$rule)
    stopifnot(all(lengths(flags) == ncol(sc)))
    for (cluster in sort(unique(ids))) {
      ix <- which(ids == cluster); key <- paste(group, cluster, sep = '|')
      m <- data.frame(subgroup = group, cluster = cluster, cells = length(ix), gene = genes,
        percent = as.numeric(Matrix::rowMeans(x[, ix, drop = FALSE] > 0) * 100),
        average = as.numeric(Matrix::rowMeans(linear[, ix, drop = FALSE])))
      stopifnot(max(abs(m$percent - e$percent[genes, cluster])) < 1e-8,
        max(abs(m$average - e$average[genes, cluster])) < 1e-8)
      marker[[key]] <- m
      md <- sc[[]][ix, , drop = FALSE]
      qc[[key]] <- data.frame(subgroup = group, cluster = cluster, cells = length(ix),
        annotation = as.character(md$subtype_annotation[1]), samples = length(unique(md$sample)),
        max_sample_fraction = max(table(md$sample)) / nrow(md),
        median_nCount_RNA = median(md$nCount_RNA), median_nFeature_RNA = median(md$nFeature_RNA),
        median_input_pct_counts_mt = median(md[[s$mitochondrial_column]]))
      for (sample in c('ALL', sort(unique(as.character(md$sample))))) {
        jj <- if (sample == 'ALL') ix else ix[as.character(sc$sample[ix]) == sample]
        k <- paste(key, sample, sep = '|')
        coexpression[[k]] <- data.frame(subgroup = group, cluster = cluster, sample = sample,
          cells = length(jj), rule = rr$rule, expression = rr$expression,
          percent = vapply(flags, function(f) mean(f[jj]) * 100, numeric(1)))
        if (sample != 'ALL') sample_marker[[k]] <- data.frame(subgroup = group, cluster = cluster,
          sample = sample, cells = length(jj), gene = genes,
          percent = as.numeric(Matrix::rowMeans(x[, jj, drop = FALSE] > 0) * 100),
          average = as.numeric(Matrix::rowMeans(linear[, jj, drop = FALSE])))
      }
    }
    stopifnot(identical(Embeddings(sc, 'umap'), e$umap))
    log_msg('ANNOTATION_REVIEW_OK ', group, '; cells=', ncol(sc), '; genes=', length(genes))
    rm(sc, x, linear, e, detected, flags); invisible(gc())
  }
  for (pair in list(c('marker', 'review_marker_expression.tsv'), c('sample_marker', 'review_sample_markers.tsv'),
    c('coexpression', 'review_coexpression.tsv'), c('qc', 'review_cluster_QC.tsv')))
    write.table(do.call(rbind, get(pair[1])), file.path(output, pair[2]), sep = '\t', quote = FALSE, row.names = FALSE)
  if (length(missing_markers)) write.table(do.call(rbind, missing_markers), file.path(output, 'missing_identity_markers.tsv'),
    sep = '\t', quote = FALSE, row.names = FALSE)
  log_msg('ALL_ANNOTATION_REVIEW_OK')
}

# 先保存marker证据及未改变的中性粒程序，再复用已有QS写入和完整读回校验。
step_apply_reviewed_annotations <- function() {
  step_review_functional_markers(CONFIG$reviewed_annotation)
  step_apply_functional_annotations(CONFIG$reviewed_annotation)
}
