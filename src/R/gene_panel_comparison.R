# 从既有表达及全基因pseudobulk结果制作候选基因专题，不改变原统计检验范围。
prepare_gene_panel <- function(feature, de, cfg) {
  genes <- cfg$genes$gene
  stopifnot(all(genes %in% rownames(feature$expression)), !anyDuplicated(genes))
  expression <- as.matrix(feature$expression[genes, , drop = FALSE])
  meta <- feature$metadata[, cfg$metadata_columns, drop = FALSE]
  names(meta) <- names(cfg$metadata_columns)
  meta[] <- lapply(meta, as.character)
  meta$cell <- rownames(meta)
  stopifnot(identical(colnames(expression), meta$cell), identical(rownames(feature$umap), meta$cell))
  annotation <- read.csv(cfg$annotation, stringsAsFactors = FALSE)
  scopes <- c(cfg$overall_label, annotation$annotation)
  stopifnot(identical(scopes, names(de$results)))
  rows <- lapply(scopes, function(s) {
    d <- de$table[de$table$scope == s, ]; d[match(genes, d$gene), ]
  })
  table <- do.call(rbind, rows); rownames(table) <- NULL
  stopifnot(!anyNA(table$gene), nrow(table) == length(scopes) * length(genes))
  table$protein <- cfg$genes$protein[match(table$gene, genes)]
  table$fold_change <- 2^table$logFC
  samples <- unique(meta[c("sample", "group")])
  samples <- samples[order(match(samples$group, cfg$group_levels), samples$sample), ]
  stopifnot(!anyDuplicated(samples$sample))
  by_sample <- list()
  for (scope in scopes) for (sample in samples$sample) {
    use <- meta$sample == sample & (scope == cfg$overall_label | meta$subtype == scope)
    m <- expression[, use, drop = FALSE]
    stopifnot(ncol(m) > 0L)
    by_sample[[paste(scope, sample)]] <- data.frame(scope, sample,
      group = samples$group[match(sample, samples$sample)], gene = genes, n_cells = ncol(m),
      mean_normalized_counts = rowMeans(expm1(m)), pct_expressed = rowMeans(m > 0) * 100)
  }
  by_sample <- do.call(rbind, by_sample); rownames(by_sample) <- NULL
  # 样本等权，先在原线性归一化尺度求均值，再log1p用于点图颜色；不按细胞数加权。
  by_group <- stats::aggregate(cbind(mean_normalized_counts, pct_expressed) ~ scope + group + gene,
    data = by_sample, FUN = mean)
  by_group$log1p_mean <- log1p(by_group$mean_normalized_counts)
  for (group in cfg$group_levels) {
    d <- by_group[by_group$group == group, ]
    i <- match(paste(table$scope, table$gene), paste(d$scope, d$gene))
    table[[paste0("RNA_mean_CP10K_", group, "_all_mice")]] <- d$mean_normalized_counts[i]
    table[[paste0("RNA_pct_", group, "_all_mice")]] <- d$pct_expressed[i]
  }
  overall <- de$results[[cfg$overall_label]]
  pseudobulk <- do.call(rbind, lapply(genes, function(g)
    data.frame(gene = g, overall$samples[c("sample", "group")],
      logcpm = as.numeric(overall$logcpm[g, overall$samples$sample]))))
  list(expression = expression, metadata = meta, umap = feature$umap, point_size = feature$point_size,
    table = table, sample_expression = by_sample, group_expression = by_group,
    pseudobulk_expression = pseudobulk, annotation = annotation, config = cfg)
}

plot_gene_panel <- function(r, cfg) {
  save <- function(p, stem, width, height) cnmf_save_plot(p, file.path(cfg$output, stem), width, height)
  genes <- cfg$genes$gene
  labels <- setNames(cfg$genes$label, genes)
  contrast <- paste(rev(cfg$group_levels), collapse = " - ")
  d <- r$pseudobulk_expression
  d$group <- factor(d$group, cfg$group_levels); d$gene <- factor(d$gene, genes)
  means <- stats::aggregate(logcpm ~ gene + group, d, mean)
  a <- r$table[r$table$scope == cfg$overall_label, ]
  a$gene <- factor(a$gene, genes)
  a$text <- sprintf("edgeR P = %.3f\nBH FDR = %.3f", a$PValue, a$FDR)
  p <- ggplot2::ggplot(d, ggplot2::aes(group, logcpm, color = group)) +
    ggplot2::geom_point(size = 2.8, position = ggplot2::position_jitter(width = 0.065, height = 0, seed = cfg$step)) +
    ggplot2::geom_segment(data = means, ggplot2::aes(x = as.numeric(group) - 0.2,
      xend = as.numeric(group) + 0.2, y = logcpm, yend = logcpm), inherit.aes = FALSE, color = "black") +
    ggplot2::geom_text(data = a, ggplot2::aes(x = 1.5, y = Inf, label = text), inherit.aes = FALSE,
      vjust = 1.15, size = 3, color = "#222222") +
    ggplot2::facet_wrap(~ gene, ncol = 4, scales = "free_y", labeller = ggplot2::labeller(gene = labels)) +
    ggplot2::scale_color_manual(values = cfg$group_colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.38))) +
    ggplot2::labs(x = NULL, y = "RNA pseudobulk log2 CPM") + cnmf_theme() +
    ggplot2::theme(legend.position = "none")
  save(p, "Fig1_Expression_by_mouse", 12, 7.5)
  d <- r$group_expression[r$group_expression$scope != cfg$overall_label, ]
  d$gene <- factor(d$gene, genes); d$group <- factor(d$group, cfg$group_levels)
  d$scope <- factor(d$scope, rev(r$annotation$annotation))
  scope_labels <- setNames(r$annotation$plot_label, r$annotation$annotation)
  p <- hongf_dotplot_style(ggplot2::ggplot(d, ggplot2::aes(gene, scope, color = log1p_mean, size = pct_expressed)) +
    ggplot2::geom_point()) +
    ggplot2::facet_wrap(~ group, nrow = 1) +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::scale_y_discrete(labels = scope_labels) +
    ggplot2::scale_size_continuous(range = c(0, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
      name = "Mean % expressing\n(equal mouse weight)") +
    ggplot2::labs(x = NULL, y = NULL, color = "log1p(mean normalized RNA)\n(equal mouse weight)") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 12, 8, 8))
  save(p, "Fig2_Expression_DotPlot_by_subtype_and_group", 10, 6.5)
  scopes <- c(cfg$overall_label, r$annotation$annotation)
  mat <- matrix(NA_real_, length(genes), length(scopes), dimnames = list(genes, scopes))
  for (scope in scopes) {
    d <- r$table[r$table$scope == scope, ]
    mat[, scope] <- d$logFC[match(genes, d$gene)]
  }
  columns <- c(cfg$overall_label, r$annotation$plot_label)
  columns[colSums(is.finite(mat)) == 0] <- paste0(columns[colSums(is.finite(mat)) == 0], " (not tested)")
  limit <- max(abs(mat), na.rm = TRUE)
  p <- cnmf_heat(mat, labels, columns, paste0(contrast, "\nlog2FC\nGray: not tested"), TRUE, c(-limit, limit))
  p$data$text <- ifelse(is.finite(p$data$value), sprintf("%.2f", p$data$value), "")
  p <- p + ggplot2::geom_text(ggplot2::aes(label = text), size = 3) +
    ggplot2::labs(title = "Gene-level pseudobulk differences")
  save(p, "Fig3_Pseudobulk_log2FC_by_subtype", 13, 7.5)
  plot_gene_panel_features(r, cfg)
  invisible(NULL)
}

plot_gene_panel_features <- function(r, cfg) {
  genes <- cfg$genes$gene
  labels <- setNames(cfg$genes$label, genes)
  # 原UMAP、全细胞、默认自动点径；每个基因在两组间统一色阶，不截断表达值。
  pages <- split(genes, ceiling(seq_along(genes) / cfg$feature_genes_per_page))
  for (page in seq_along(pages)) {
    panels <- list()
    for (gene in pages[[page]]) for (group in cfg$group_levels) {
      use <- r$metadata$group == group
      d <- data.frame(r$umap[use, , drop = FALSE], expression = r$expression[gene, use])
      names(d)[1:2] <- c("UMAP_1", "UMAP_2")
      d <- d[order(d$expression), ]
      feature <- list(umap = r$umap, point_size = r$point_size,
        color_limits = c(0, max(1, ceiling(max(r$expression[gene, ])))) )
      panels[[paste(gene, group)]] <- ferrdb_feature_panel(d, paste(labels[gene], group, sep = " | "), feature, cfg)
    }
    cnmf_save_plot(patchwork::wrap_plots(panels, ncol = cfg$feature_ncol),
      file.path(cfg$output, sprintf("Fig4_FeaturePlot_EXP_vs_CON_%02d", page)),
      cfg$feature_width, cfg$feature_height)
  }
  invisible(NULL)
}

run_gene_panel_comparison <- function(cfg) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  protected <- tools::md5sum(c(cfg$feature_cache, cfg$de_cache))
  feature <- qs2::qs_read(cfg$feature_cache)
  de <- qs2::qs_read(cfg$de_cache)
  stopifnot(identical(feature$source_md5, tools::md5sum(names(feature$source_md5))),
    identical(de$source_md5, tools::md5sum(names(de$source_md5))))
  result <- prepare_gene_panel(feature, de, cfg)
  result$source_md5 <- protected
  result$session <- sessionInfo()
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(result, cfg$cache)
  data.table::fwrite(result$table, file.path(cfg$output, "Selected_genes_EXP_vs_CON.csv"), bom = TRUE)
  plot_gene_panel(result, cfg)
  stopifnot(identical(protected, tools::md5sum(names(protected))))
  print(result$table[result$table$scope == cfg$overall_label, c("protein", "gene", "logFC", "fold_change", "PValue", "FDR")])
  message("GENE_PANEL_OK cells=", ncol(result$expression), " genes=", nrow(result$expression),
    " comparisons=", nrow(result$table), " tested=", sum(result$table$test_status == "tested"))
  invisible(result)
}
