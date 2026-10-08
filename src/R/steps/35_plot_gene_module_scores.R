# 复用注释、细胞数和原UMAP，展示模块分数及逐样本不均衡，不做组间显著性检验。
step_plot_gene_module_scores <- function(s) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  folder <- file.path(PROJECT_DIR, s$evidence_dir)
  out <- file.path(PROJECT_DIR, s$results_dir)
  e <- readRDS(file.path(folder, 'score_evidence.rds'))
  clusters <- as.character(e$mapping$subgroup_cluster)
  labels <- setNames(e$mapping$display, clusters)
  modules <- e$programs$module
  module_labels <- setNames(e$programs$label, modules)
  d <- do.call(rbind, lapply(c('cell_mean_z', 'sample_balanced_mean_z', 'sample_balanced_rare_filtered_z'), function(k) {
    data.frame(e$summary[, c('cluster', 'module')], value = e$summary[[k]], weighting = k)
  }))
  d$weighting <- factor(d$weighting,
    levels = c('cell_mean_z', 'sample_balanced_mean_z', 'sample_balanced_rare_filtered_z'),
    labels = c('Cell-weighted', 'Sample-balanced', paste0('Sample-balanced: genes >=', s$rare_detection_percent, '%')))
  lim <- max(abs(d$value))
  p <- ggplot(d, aes(module, as.character(cluster), fill = value)) + geom_tile(colour = 'white', linewidth = .3) +
    geom_text(aes(label = sprintf('%.2f', value), colour = abs(value) > .75 * lim), size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = c('FALSE' = 'black', 'TRUE' = 'white')) +
    scale_fill_gradient2(low = HONGF_NPG_DIVERGING[1], mid = HONGF_NPG_DIVERGING[2], high = HONGF_NPG_DIVERGING[3],
      limits = c(-lim, lim), guide = guide_colorbar(display = 'rectangles')) +
    scale_x_discrete(limits = modules, labels = module_labels) +
    scale_y_discrete(limits = rev(clusters), labels = labels) + facet_grid(~weighting) +
    labs(x = NULL, y = NULL, fill = 'Mean gene z',
      caption = 'Compare clusters within each module. RNA programs do not establish N1/N2 identity or tumor function.') +
    theme_bw(base_size = 10) + theme(panel.grid = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1), axis.text.y = element_text(colour = 'black'),
      strip.background = element_blank(), legend.position = 'bottom', plot.caption = element_text(hjust = 0, size = 9))
  save_plot(p, file.path(out, 'Heatmap_figure_scores'), width = 14, height = 5.8)
  # 每个点为一个实际样本；稀少的cluster×sample用空心标记，缺失不补点。
  a <- e$sample_summary[e$sample_summary$cells > 0, ]
  sample_levels <- sort(unique(a$sample))
  a$y <- match(as.character(a$cluster), rev(clusters)) + (match(a$sample, sample_levels) - (length(sample_levels) + 1) / 2) * .095
  a$module <- factor(a$module, levels = modules)
  a$sample <- factor(a$sample, levels = sample_levels)
  avg <- e$summary; avg$y <- match(as.character(avg$cluster), rev(clusters)); avg$module <- factor(avg$module, levels = modules)
  p_sample <- ggplot(a, aes(mean_z, y)) + geom_vline(xintercept = 0, colour = '#DDDDDD', linewidth = .4) +
    geom_point(aes(colour = sample, shape = small_sample), size = 1.6, stroke = .65) +
    geom_point(data = avg, aes(x = sample_balanced_mean_z, y = y), shape = 124, size = 3.5, colour = 'black') +
    scale_colour_manual(values = npg_palette(sample_levels)) +
    scale_shape_manual(values = c('FALSE' = 16, 'TRUE' = 1), labels = c('FALSE' = paste0('>=', s$small_sample_cells, ' cells'), 'TRUE' = paste0('<', s$small_sample_cells, ' cells'))) +
    scale_y_continuous(breaks = seq_along(clusters), labels = unname(labels[rev(clusters)]), limits = c(.5, length(clusters) + .5)) +
    facet_wrap(~module, nrow = 1, scales = 'free_x', labeller = as_labeller(setNames(gsub(' ', '\n', module_labels), modules))) +
    labs(x = 'Sample mean gene z', y = NULL, colour = 'Sample', shape = 'Cluster within sample',
      caption = 'Black ticks: equal-sample means. A1-A3 = CON; B1-B3 = TREAT. Missing combinations have no point.') +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = '#EEEEEE'),
      panel.grid.major.x = element_blank(), axis.text.y = element_text(colour = 'black'), strip.background = element_blank(),
      legend.position = 'bottom', plot.caption = element_text(hjust = 0, size = 9))
  save_plot(p_sample, file.path(out, 'Scores_figure_by_sample'), width = 15, height = 5.6)
  m <- e$mapping; m$candidate_annotation <- m$annotation
  dot <- candidate_dotplot(e, m, rownames(e$expression), '', scaled = TRUE)
  save_plot(dot$plot, file.path(out, 'DotPlot_figure_score_genes'), width = 12, height = 4.6)
  saveRDS(list(heatmap = d, sample_points = a, dotplot = dot$data), file.path(folder, 'figure_plot_data.rds'))
  # 沿用Seurat默认点径，关闭点和色条的栅格化；不裁剪极端分数。
  selected <- modules[seq_len(min(4, length(modules)))]
  umaps <- lapply(selected, function(k) {
    p <- FeaturePlot(e$plot_object, features = paste0('figure_score_', k), reduction = 'umap',
      raster = FALSE, order = TRUE, combine = FALSE)[[1]]
    p + scale_colour_gradient2(low = HONGF_NPG_DIVERGING[1], mid = HONGF_NPG_DIVERGING[2], high = HONGF_NPG_DIVERGING[3],
      midpoint = 0, guide = guide_colorbar(display = 'rectangles')) +
      labs(title = module_labels[[k]], colour = 'Mean gene z') + theme(plot.title = element_text(size = 11))
  })
  save_umap_plot(wrap_plots(umaps, ncol = 2), file.path(out, 'UMAP_figure_scores'), width = 10, height = 8)
  log_msg('GENE_MODULE_FIGURES_OK figures=4; vector UMAP; no cluster-number labels')
}
