# 谱系整体重新拟合小鼠pseudobulk；亚群复用原有全基因检验，筛选FerrDb不会重算BH。
run_ferrdb_pseudobulk_comparison <- function(cfg) {
  feature <- qs2::qs_read(cfg$feature_cache)
  previous <- qs2::qs_read(cfg$subtype_cache)
  stopifnot(identical(previous$source_md5, tools::md5sum(names(previous$source_md5))),
    identical(feature$source_md5, tools::md5sum(names(feature$source_md5))))
  total <- Reduce(`+`, lapply(previous$results, `[[`, "counts"))
  info <- previous$results[[1]]$all_samples
  info$n_cells <- Reduce(`+`, lapply(previous$results, function(a)
    if (is.null(a$all_samples)) a$samples$n_cells else a$all_samples$n_cells))
  group <- factor(info$group, levels = cfg$group_levels)
  stopifnot(identical(info$sample, colnames(total)), all(table(group) == 3L))
  design <- model.matrix(~ group)
  y <- edgeR::DGEList(total, group = group)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  de <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = "groupEXP"), n = Inf, sort.by = "none")$table
  de$gene <- rownames(de)
  stats <- c("logFC", "logCPM", "F", "PValue", "FDR")
  results <- c(list(`All EMT Epi` = list(de = de, counts = total, samples = info,
    keep = keep, design = design, library = y$samples, logcpm = edgeR::cpm(y, log = TRUE, prior.count = 2))), previous$results)
  tables <- lapply(names(results), function(scope) {
    a <- results[[scope]]
    m <- if (is.null(a$all_samples)) a$samples else a$all_samples
    n <- table(factor(m$group[m$used], levels = cfg$group_levels))
    if (scope == "All EMT Epi") n <- table(group)
    out <- data.frame(scope, gene = feature$catalogue$gene, roles = feature$catalogue$roles,
      n_CON = unname(n[1]), n_EXP = unname(n[2]), test_status = "tested")
    if (is.null(a$de)) {
      out$test_status <- "insufficient_mice"; out[, stats] <- NA_real_
    } else {
      index <- match(out$gene, a$de$gene)
      out[, stats] <- a$de[index, stats]
      out$test_status[is.na(index)] <- "not_tested_expression_filter"
    }
    out$direction <- ifelse(out$test_status != "tested", "Not tested",
      ifelse(out$FDR < cfg$fdr & out$logFC >= cfg$logfc, "EXP higher",
        ifelse(out$FDR < cfg$fdr & out$logFC <= -cfg$logfc, "CON higher", "Not significant")))
    out
  })
  table <- do.call(rbind, tables)
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(table, file.path(cfg$output, "FerrDb_EXP_vs_CON_pseudobulk.csv"), bom = TRUE)
  result <- list(results = results, table = table, catalogue = feature$catalogue,
    config = cfg, source_md5 = tools::md5sum(c(cfg$feature_cache, cfg$subtype_cache)), session = sessionInfo())
  qs2::qs_save(result, cfg$cache)
  plot_ferrdb_comparison(result, cfg)
  for (scope in names(results)) {
    d <- table[table$scope == scope, ]
    message(scope, " tested=", sum(d$test_status == "tested"), " FDR<0.05=", sum(d$FDR < cfg$fdr, na.rm = TRUE),
      " EXP_higher=", sum(d$direction == "EXP higher"), " CON_higher=", sum(d$direction == "CON higher"))
  }
  invisible(result)
}

plot_ferrdb_comparison <- function(x, cfg) {
  d <- x$table[x$table$scope == "All EMT Epi" & x$table$test_status == "tested", ]
  colors <- c(`CON higher` = cfg$colors[1], `Not significant` = "#B5B5B5", `EXP higher` = cfg$colors[2])
  names(colors) <- c("CON higher", "Not significant", "EXP higher")
  labels <- head(d[order(d$FDR, d$PValue), ], 8)
  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, -log10(pmax(FDR, .Machine$double.xmin)), color = direction)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-cfg$logfc, cfg$logfc), linetype = 2, color = "#999999") +
    ggplot2::geom_hline(yintercept = -log10(cfg$fdr), linetype = 2, color = "#999999") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = gene), color = "#333333", show.legend = FALSE,
      size = 3, max.overlaps = Inf, seed = 106) +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::labs(title = "EMT Epi | FerrDb genes", x = "log2 fold change (EXP - CON)", y = "-log10(BH FDR)") + cnmf_theme()
  cnmf_save_plot(p, file.path(cfg$output, "Fig1_FerrDb_pseudobulk_volcano"), 8, 5.5)
  plot_ferrdb_heatmaps(x, cfg)
}

# 只改变基因行序，不改变列序、数值或NA；全NA列不参与距离计算。
ferrdb_row_cluster <- function(mat, method = "complete") {
  observed <- mat[, colSums(is.finite(mat)) > 0L, drop = FALSE]
  stopifnot(nrow(observed) > 1L, ncol(observed) > 0L)
  # dist按每对基因共有的非缺失列计算，并按总列数/共有列数缩放平方距离。
  distance <- stats::dist(observed, method = "euclidean")
  if (any(!is.finite(distance))) stop("Some gene pairs have no shared finite values for row clustering")
  stats::hclust(distance, method = method)
}

plot_ferrdb_heatmaps <- function(x, cfg) {
  output <- list()
  cluster_rows <- isTRUE(cfg$cluster_rows)
  method <- if (is.null(cfg$cluster_method)) "complete" else cfg$cluster_method
  d <- x$table[x$table$scope == "All EMT Epi" & x$table$test_status == "tested", ]
  genes <- head(d$gene[order(d$FDR, d$PValue)], cfg$heatmap_top)
  a <- x$results[["All EMT Epi"]]
  z <- t(scale(t(a$logcpm[genes, , drop = FALSE]))); z[!is.finite(z)] <- 0
  tree <- if (cluster_rows) ferrdb_row_cluster(z, method) else NULL
  if (cluster_rows) z <- z[tree$order, , drop = FALSE]
  p <- cnmf_heat(z, setNames(genes, genes), paste0(a$samples$sample, "\n", a$samples$group),
    legend = "logCPM\nrow Z-score", divergent = TRUE, limits = c(-2, 2)) +
    ggplot2::labs(title = "Top-ranked FerrDb genes")
  cnmf_save_plot(p, file.path(cfg$output, "Fig2_FerrDb_top_genes_by_mouse"), 8.5, 10)
  output$Fig2 <- list(matrix = z, cluster = tree, plot = p)
  tested <- x$table[x$table$test_status == "tested", ]
  best <- tapply(tested$FDR, tested$gene, min)
  genes <- head(names(sort(best)), cfg$heatmap_top)
  scopes <- names(x$results)
  mat <- matrix(NA_real_, length(genes), length(scopes), dimnames = list(genes, scopes))
  stars <- list()
  for (scope in scopes) {
    z <- x$table[x$table$scope == scope, ]; z <- z[match(genes, z$gene), ]
    mat[, scope] <- z$logFC
    hit <- which(!is.na(z$FDR) & z$FDR < cfg$fdr)
    if (length(hit)) stars[[scope]] <- data.frame(row = genes[hit], column = scope)
  }
  column_labels <- ifelse(vapply(x$results, function(r) is.null(r$de), logical(1)),
    paste0(scopes, " (not tested)"), scopes)
  tree <- if (cluster_rows) ferrdb_row_cluster(mat, method) else NULL
  if (cluster_rows) mat <- mat[tree$order, , drop = FALSE]
  p <- cnmf_heat(mat, setNames(genes, genes), column_labels,
    legend = "EXP - CON\nlog2FC\nGray: not tested", divergent = TRUE, limits = c(-2, 2)) +
    ggplot2::labs(title = "FerrDb genes | overall and within subtypes")
  if (length(stars)) {
    # cnmf_heat的坐标字段是row/column，星号仅表示原全基因BH FDR<0.05。
    s <- do.call(rbind, stars)
    p <- p + ggplot2::geom_text(data = s, ggplot2::aes(column, row), label = "*", inherit.aes = FALSE, size = 4)
  }
  cnmf_save_plot(p, file.path(cfg$output, "Fig3_FerrDb_effect_by_subtype"), 12, 10)
  output$Fig3 <- list(matrix = mat, cluster = tree, plot = p)
  invisible(output)
}
