# 参照STARTRAC-dist的Ro/e分级设置圆点大小；颜色保留连续值，0值不加伪计数。
roe_dotplot <- function(d, category_levels, split_levels, title) {
  d$cell_type <- factor(d$cell_type, levels = category_levels)
  d$stratum <- factor(d$stratum, levels = split_levels)
  limit <- max(2, d$ROE, na.rm = TRUE)
  grade_levels <- c('0', '(0, 0.2]', '(0.2, 0.8]', '(0.8, 1]', '>1')
  d$ROE_grade <- factor(ifelse(d$ROE == 0, '0', ifelse(d$ROE <= .2, '(0, 0.2]',
    ifelse(d$ROE <= .8, '(0.2, 0.8]', ifelse(d$ROE <= 1, '(0.8, 1]', '>1')))), levels = grade_levels)
  breaks <- c(0, .5, 1, 2, 4, 8); breaks <- breaks[breaks <= limit]
  ggplot(d, aes(stratum, cell_type)) +
    geom_point(data = d[!is.na(d$ROE) & d$ROE > 0, ], aes(size = ROE_grade, colour = ROE)) +
    scale_size_manual(values = setNames(c(0, 2, 4, 6, 8), grade_levels), breaks = grade_levels[-1], drop = FALSE) +
    scale_colour_gradient2(low = HONGF_NPG_DIVERGING[1], mid = '#D6D6D6', high = HONGF_NPG_DIVERGING[3],
      midpoint = 1, limits = c(0, limit), breaks = breaks, guide = guide_colorbar(display = 'rectangles')) +
    scale_x_discrete(limits = split_levels, drop = FALSE) +
    scale_y_discrete(limits = rev(category_levels), drop = FALSE) +
    labs(x = NULL, y = NULL, colour = 'Ro/e', size = 'Ro/e class', title = title,
      caption = 'Size: STARTRAC-dist classes; colour: Ro/e; no dot = 0.') +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank(),
      panel.border = element_rect(colour = 'black', linewidth = .5), axis.ticks = element_blank(),
      axis.text = element_text(colour = 'black'), plot.title = element_text(size = 12),
      plot.caption = element_text(hjust = 0, size = 8), legend.box = 'vertical')
}

step_replot_roe_dots <- function(s) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  e <- readRDS(file.path(PROJECT_DIR, s$evidence_dir, 'composition_evidence.rds'))
  pdf(file.path(PROJECT_DIR, s$main_results, 'All_composition_ROE.pdf'), width = 18, height = 6.2, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  dot_data <- list()
  for (scope in names(e$tables)) {
    main <- scope == 'Main_lineages'
    levels <- if (main) unlist(s$major_levels) else e$palette$label[e$palette$subgroup == scope]
    palette <- if (main) npg_palette(levels) else setNames(e$palette$color[e$palette$subgroup == scope], levels)
    out <- file.path(PROJECT_DIR, if (main) s$main_results else file.path(s$subgroup_results, scope))
    title <- if (main) 'Main cell lineages' else gsub('_', ' ', scope)
    for (kind in c('sample', 'group')) {
      d <- e$tables[[scope]]; d <- d[d$split_by == kind, ]
      split_levels <- unlist(s[[paste0(kind, '_levels')]])
      p <- composition_roe_plots(d, levels, split_levels, palette, title, kind)
      p$ROE <- roe_dotplot(d, levels, split_levels, title)
      save_plot(p$ROE, file.path(out, paste0('ROE_by_', kind)),
        width = if (kind == 'sample') 8.4 else 6.8, height = max(4.3, length(levels) * .38 + 1.6))
      # 保存实际绘图映射以核对分级大小和连续颜色；不改变原始ROE。
      dot_data[[paste(scope, kind, sep = '|')]] <- p$ROE$data
      overview <- wrap_plots(p$Counts + theme(legend.position = 'none'), p$Proportions, p$ROE,
        nrow = 1, widths = c(1, 1.5, 1.3)) + plot_annotation(title = paste0(title, ' - by ', kind))
      print(overview)
    }
    log_msg('ROE_DOTPLOTS_OK scope=', scope)
  }
  saveRDS(dot_data, file.path(PROJECT_DIR, s$evidence_dir, 'roe_dotplot_data.rds'))
  log_msg('ALL_ROE_DOTPLOTS_OK plots=', length(dot_data), '; counts and proportions preserved')
}
