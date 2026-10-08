# Summarize the existing population; recompute count-based metrics without filtering.
step_sample_qc <- function() {
  suppressPackageStartupMessages(library(SeuratObject))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$sample_qc
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  object <- qs::qread(file.path(PROJECT_DIR, CONFIG$merged_output), nthreads = 4)
  metadata <- object[[]]
  counts <- LayerData(object, assay = "RNA", layer = "counts")
  stopifnot(!anyDuplicated(colnames(object)), identical(colnames(counts), rownames(metadata)),
            all(is.finite(counts@x)), all(counts@x >= 0))
  integer_counts <- all(abs(counts@x - round(counts@x)) < 1e-8)
  stopifnot(integer_counts)
  mito <- grep(detect_mito_pattern(rownames(counts)), rownames(counts), value = TRUE)
  ribo <- grep(settings$ribosomal_pattern, rownames(counts), value = TRUE)
  stopifnot(length(mito) > 0L, length(ribo) > 0L)
  qc <- data.frame(
    cell = colnames(object), sample = as.character(metadata[[settings$sample_col]]),
    group = as.character(metadata[[settings$group_col]]),
    genes = Matrix::colSums(counts > 0), umi = Matrix::colSums(counts),
    mito_counts = Matrix::colSums(counts[mito, , drop = FALSE]),
    ribo_counts = Matrix::colSums(counts[ribo, , drop = FALSE]),
    doublet_score = metadata$doublet_score, predicted_doublet = metadata$predicted_doublet,
    cell_type = as.character(metadata$cell_type), stringsAsFactors = FALSE
  )
  stopifnot(!anyNA(qc), all(nzchar(qc$sample)), all(nzchar(qc$group)), all(qc$umi > 0))
  qc$mito_percent <- qc$mito_counts / qc$umi * 100
  qc$ribo_percent <- qc$ribo_counts / qc$umi * 100
  samples <- sort(unique(qc$sample))
  qc$sample <- factor(qc$sample, levels = samples)
  comparisons <- list(genes = metadata$nFeature_RNA, umi = metadata$nCount_RNA,
                      mito_percent = metadata$pct_counts_mt, ribo_percent = metadata$pct_counts_ribo)
  validation <- do.call(rbind, lapply(names(comparisons), function(metric) {
    delta <- abs(qc[[metric]] - comparisons[[metric]])
    data.frame(metric = metric, max_absolute_difference = max(delta),
               cells_difference_above_0_001 = sum(delta > .001))
  }))
  legacy_delta <- abs(qc$umi - metadata$total_counts)
  validation <- rbind(validation, data.frame(metric = "umi_vs_legacy_total_counts",
                      max_absolute_difference = max(legacy_delta),
                      cells_difference_above_0_001 = sum(legacy_delta > .001)))
  stopifnot(all(abs(qc$mito_counts - metadata$total_counts_mt) < 1e-8),
            all(abs(qc$ribo_counts - metadata$total_counts_ribo) < 1e-8))
  summary <- do.call(rbind, lapply(samples, function(sample) {
    x <- qc[qc$sample == sample, ]
    stopifnot(length(unique(x$group)) == 1L)
    row <- data.frame(sample = sample, group = unique(x$group), cells = nrow(x),
                      percent_of_all_cells = nrow(x) / nrow(qc) * 100,
                      predicted_doublets = sum(x$predicted_doublet),
                      genes_below_300 = sum(x$genes < 300), umi_below_1000 = sum(x$umi < 1000))
    for (metric in c("genes", "umi", "mito_percent", "ribo_percent", "doublet_score")) {
      for (quantile in c(.05, .25, .5, .75, .95)) {
        row[[paste0(metric, "_p", quantile * 100)]] <- unname(stats::quantile(x[[metric]], quantile))
      }
    }
    for (cutoff in unlist(settings$reference_mito_percent)) {
      row[[paste0("mito_above_", cutoff, "_cells")]] <- sum(x$mito_percent > cutoff)
      row[[paste0("mito_above_", cutoff, "_percent")]] <- mean(x$mito_percent > cutoff) * 100
    }
    row
  }))
  stopifnot(sum(summary$cells) == ncol(object), !anyDuplicated(qc$cell))
  npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
           "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
  stopifnot(length(samples) <= length(npg))
  colors <- setNames(npg[seq_along(samples)], samples)
  style <- theme_classic(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11),
          axis.title.x = element_blank(), plot.margin = margin(8, 8, 8, 8))
  count_plot <- ggplot(summary, aes(factor(sample, levels = samples), cells, fill = sample)) +
    geom_col(width = .72) + geom_text(aes(label = format(cells, big.mark = ",")), vjust = -.5, size = 3) +
    scale_fill_manual(values = colors) + scale_y_continuous(expand = expansion(mult = c(0, .13))) +
    labs(title = "Cells per sample", y = "Cells") + style
  labels <- c(genes = "Detected genes", umi = "UMI counts", mito_percent = "Mitochondrial counts (%)",
              ribo_percent = "Ribosomal counts (%)", doublet_score = "Existing doublet score")
  violins <- lapply(names(labels), function(metric) {
    ggplot(qc, aes(sample, .data[[metric]], fill = sample)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = .25) +
      geom_boxplot(width = .12, outlier.shape = NA, fill = "white", linewidth = .3) +
      scale_fill_manual(values = colors) + labs(title = labels[[metric]], y = labels[[metric]]) + style
  })
  overview <- patchwork::wrap_plots(c(list(count_plot), violins), ncol = 3)
  scatter <- ggplot(qc, aes(umi, genes, color = mito_percent)) +
    geom_point(size = .16, alpha = .55, stroke = 0) +
    facet_wrap(~sample, ncol = 3) +
    scale_color_gradientn(colors = c("#DAEEF3", "#4DBBD5", "#3C5488", "#E64B35"),
                          name = "Mitochondrial\ncounts (%)", limits = c(0, max(qc$mito_percent))) +
    guides(color = guide_colorbar(display = "rectangles")) +
    labs(x = "UMI counts", y = "Detected genes") + theme_classic(base_size = 10) +
    theme(strip.background = element_rect(fill = "white", linewidth = .3),
          strip.text = element_text(face = "bold"))
  plots <- list(QC_sample_overview = overview, QC_genes_vs_UMI_by_sample = scatter)
  for (name in names(plots)) {
    width <- if (name == "QC_sample_overview") 12 else 11
    height <- if (name == "QC_sample_overview") 7.5 else 7
    ggsave(file.path(output, paste0(name, ".png")), plots[[name]], width = width, height = height,
           dpi = 600, bg = "white", device = ragg::agg_png)
    ggsave(file.path(output, paste0(name, ".pdf")), plots[[name]], width = width, height = height,
           bg = "white", device = grDevices::cairo_pdf)
  }
  write.table(summary, file.path(output, "sample_qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(validation, file.path(output, "metric_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  connection <- gzfile(file.path(output, "cell_qc_metrics.tsv.gz"), "wt")
  tryCatch(write.table(qc, connection, sep = "\t", quote = FALSE, row.names = FALSE),
           finally = close(connection))
  writeLines(c("Mitochondrial genes:", mito, "", "Ribosomal genes:", ribo),
             file.path(output, "qc_gene_sets.txt"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "04_sample_qc_sessionInfo.txt"))
  print(summary[, c("sample", "group", "cells", "genes_p50", "umi_p50", "mito_percent_p50",
                      "mito_above_20_percent", "predicted_doublets")], row.names = FALSE)
  print(validation, row.names = FALSE)
  message("SAMPLE_QC_OK: ", length(samples), " samples; ", nrow(qc),
          " unique cells; ", length(mito), " mitochondrial genes; no cells filtered.")
  invisible(summary)
}
