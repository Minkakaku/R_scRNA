# 在既有亚群聚类上计算 RNA marker；使用全部细胞，不改变表达、降维或聚类。
step_subgroup_marker_evidence <- function(settings = CONFIG$subgroup_marker_review,
  session_log = "17_subgroup_marker_evidence_sessionInfo.txt") {
  suppressPackageStartupMessages(library(Seurat))
  future::plan(future::sequential)
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest))
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  summaries <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    folder <- file.path(output, group)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    stopifnot(!file.exists(file.path(folder, "expression_evidence.rds")))
    sc <- run_step(paste(group, "Read"), qs::qread(file.path(PROJECT_DIR, manifest$checkpoint[i]), nthreads = 4))
    Idents(sc) <- factor(as.character(sc[[]][[settings$cluster_column]]),
      levels = as.character(sort(unique(as.integer(as.character(sc[[]][[settings$cluster_column]]))))))
    stopifnot(ncol(sc) == manifest$cells[i], nlevels(Idents(sc)) == manifest$clusters[i])
    markers <- run_step(paste(group, "RNA Wilcoxon all cells"), FindAllMarkers(sc, assay = "RNA",
      slot = "data", only.pos = TRUE, test.use = "wilcox", min.pct = settings$min_pct,
      logfc.threshold = settings$logfc_threshold, max.cells.per.ident = Inf,
      return.thresh = 1, random.seed = SEED, verbose = TRUE))
    markers$pct_difference <- markers$pct.1 - markers$pct.2
    con <- gzfile(file.path(folder, "all_positive_markers.tsv.gz"), "wt")
    write.table(markers, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
    candidates <- markers[markers$p_val_adj < settings$adjusted_p_cutoff & markers$pct_difference > 0 &
      !grepl("^(mt-|Rpl|Rps|Gm[0-9]|[0-9].*Rik$)", markers$gene), ]
    candidates <- candidates[order(-candidates$pct_difference, -candidates$avg_log2FC), ]
    top <- do.call(rbind, lapply(levels(Idents(sc)), function(id) head(candidates[as.character(candidates$cluster) == id, ], 50)))
    write.table(top, file.path(folder, "top50_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # 全基因的线性归一化均值和非零比例支持后续独立选择 marker，避免只看差异显著基因。
    expr <- LayerData(sc, assay = "RNA", layer = "data")
    ids <- Idents(sc); sizes <- as.integer(table(ids))
    membership <- Matrix::sparse.model.matrix(~ 0 + ids)
    linear <- expr; linear@x <- expm1(linear@x)
    avg <- sweep(as.matrix(linear %*% membership), 2, sizes, "/")
    pct <- sweep(as.matrix((expr > 0) %*% membership), 2, sizes, "/") * 100
    colnames(avg) <- colnames(pct) <- levels(ids)
    sample_counts <- as.data.frame(table(cluster = ids, sample = sc$sample))
    write.table(sample_counts, file.path(folder, "cluster_sample_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    saveRDS(list(average = avg, percent = pct, metadata = sc[[]], umap = Embeddings(sc, "umap"),
      input = manifest$checkpoint[i], input_sha256 = digest::digest(file = file.path(PROJECT_DIR, manifest$checkpoint[i]), algo = "sha256"),
      settings = settings, cells = ncol(sc), clusters = levels(ids)), file.path(folder, "expression_evidence.rds"))
    summaries[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nlevels(ids),
      positive_markers = nrow(markers), significant_markers = sum(markers$p_val_adj < settings$adjusted_p_cutoff))
    write.table(do.call(rbind, summaries), file.path(output, "marker_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("MARKER_EVIDENCE_OK ", group, " cells=", ncol(sc), " clusters=", nlevels(ids))
    rm(sc, expr, linear, markers, candidates, top, avg, pct); invisible(gc())
  }
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, session_log))
}
