# 同名注释合并计数，按样本/合并分组展示绝对细胞数、构成比例和观察/期望比。
composition_roe_table <- function(meta, category, split, category_levels, split_levels) {
  observed <- table(factor(meta[[category]], levels = category_levels), factor(meta[[split]], levels = split_levels))
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  d <- expand.grid(cell_type = category_levels, stratum = split_levels, stringsAsFactors = FALSE)
  d$observed <- as.vector(observed); d$expected <- as.vector(expected)
  d$stratum_total <- rep(colSums(observed), each = nrow(observed))
  d$type_total <- rep(rowSums(observed), times = ncol(observed))
  d$fraction <- ifelse(d$stratum_total > 0, d$observed / d$stratum_total, NA_real_)
  d$ROE <- ifelse(d$expected > 0, d$observed / d$expected, NA_real_)
  d
}

composition_roe_plots <- function(d, category_levels, split_levels, palette, title, split_kind) {
  d$cell_type <- factor(d$cell_type, levels = category_levels)
  d$stratum <- factor(d$stratum, levels = split_levels)
  totals <- unique(d[, c('stratum', 'stratum_total')])
  common <- ggplot(d, aes(stratum, fill = cell_type)) +
    scale_fill_manual(values = palette, drop = FALSE) + labs(x = NULL, fill = NULL, title = title) +
    theme_classic(base_size = 11) + theme(axis.text = element_text(colour = 'black'),
      plot.title = element_text(size = 12), legend.key.height = grid::unit(4.5, 'mm'),
      legend.text = element_text(size = 9), plot.margin = margin(8, 12, 8, 8))
  bar_count <- common + geom_col(aes(y = observed), width = .7, colour = 'white', linewidth = .2) +
    geom_text(data = totals, aes(x = stratum, y = stratum_total, label = format(stratum_total, big.mark = ',', trim = TRUE)),
      inherit.aes = FALSE, vjust = -.5, size = 3.2) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, .10))) +
    labs(y = 'Observed cells', caption = if (split_kind == 'group') 'Counts pooled across samples within each group.' else NULL)
  bar_prop <- common + geom_col(aes(y = fraction), width = .7, colour = 'white', linewidth = .2) +
    scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, .25), limits = c(0, 1), expand = c(0, 0)) +
    labs(y = 'Cell proportion', caption = if (split_kind == 'group') 'Pooled cell proportions; not the mean of sample proportions.' else NULL)
  limit <- max(2, d$ROE, na.rm = TRUE)
  heatmap <- ggplot(d, aes(stratum, cell_type, fill = ROE)) + geom_tile(colour = 'white', linewidth = .5) +
    geom_text(aes(label = ifelse(is.na(ROE), 'NA', sprintf('%.2f', ROE)), colour = ROE > 1 + .65 * (limit - 1)),
      size = 3.3, show.legend = FALSE) +
    scale_colour_manual(values = c('FALSE' = 'black', 'TRUE' = 'white'), na.value = 'black') +
    scale_y_discrete(limits = rev(category_levels), drop = FALSE) + scale_x_discrete(drop = FALSE) +
    scale_fill_gradient2(low = HONGF_NPG_DIVERGING[1], mid = HONGF_NPG_DIVERGING[2], high = HONGF_NPG_DIVERGING[3],
      midpoint = 1, limits = c(0, limit), na.value = '#DDDDDD', guide = guide_colorbar(display = 'rectangles')) +
    labs(x = NULL, y = NULL, title = title, fill = 'ROE', caption = 'ROE = observed / expected; 1 = expected composition.') +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), axis.text = element_text(colour = 'black'),
      plot.title = element_text(size = 12), plot.caption = element_text(size = 9, hjust = 0))
  list(Counts = bar_count, Proportions = bar_prop, ROE = heatmap)
}

step_composition_and_roe <- function(s) {
  suppressPackageStartupMessages(library(SeuratObject))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  folder <- file.path(PROJECT_DIR, s$evidence_dir)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROJECT_DIR, s$main_results), recursive = TRUE, showWarnings = FALSE)
  manifest <- read.delim(file.path(PROJECT_DIR, s$input_manifest))
  palette_table <- read.delim(file.path(PROJECT_DIR, s$palette))
  annotated <- readRDS(file.path(PROJECT_DIR, s$annotation_evidence, 'current_cell_annotations.rds'))
  all_meta <- hashes <- list()
  for (scope in unlist(s$scopes)) {
    input <- file.path(PROJECT_DIR, s$annotation_evidence, scope, 'annotation_plot_data.rds')
    cache <- readRDS(input); m <- cache$plot_object[[]]
    m$cell_id <- rownames(m); m$analysis_container <- scope
    ix <- match(m$cell_id, annotated$cell_id)
    stopifnot(!anyNA(ix), nrow(m) == manifest$cells[manifest$subgroup == scope],
      identical(as.character(m[[s$subtype_column]]), as.character(annotated$subtype_annotation[ix])),
      identical(as.character(m[[s$major_column]]), as.character(annotated$biological_lineage[ix])))
    all_meta[[scope]] <- m[, c('cell_id', 'analysis_container', s$subtype_column, s$major_column, s$sample_column, s$group_column)]
    hashes[[scope]] <- data.frame(path = input, sha256 = digest::digest(file = input, algo = 'sha256'))
  }
  meta <- do.call(rbind, all_meta); rownames(meta) <- meta$cell_id
  stopifnot(!anyDuplicated(meta$cell_id), nrow(meta) == s$expected_cells, setequal(meta$cell_id, annotated$cell_id))
  sample_groups <- unique(meta[, c(s$sample_column, s$group_column)])
  stopifnot(nrow(sample_groups) == length(s$sample_levels), !anyDuplicated(sample_groups[[s$sample_column]]))
  tables <- figures <- list()
  pdf(file.path(PROJECT_DIR, s$main_results, 'All_composition_ROE.pdf'), width = 18, height = 6.2, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  for (scope in c('Main_lineages', unlist(s$scopes))) {
    main <- scope == 'Main_lineages'
    part <- if (main) meta else all_meta[[scope]]
    category <- if (main) s$major_column else s$subtype_column
    levels <- if (main) unlist(s$major_levels) else palette_table$label[palette_table$subgroup == scope]
    palette <- if (main) npg_palette(levels) else setNames(palette_table$color[palette_table$subgroup == scope], levels)
    stopifnot(setequal(as.character(part[[category]]), levels), !anyDuplicated(unname(palette)))
    out <- file.path(PROJECT_DIR, if (main) s$main_results else file.path(s$subgroup_results, scope))
    title <- if (main) 'Main cell lineages' else gsub('_', ' ', scope)
    scope_tables <- list()
    for (kind in c('sample', 'group')) {
      split <- s[[paste0(kind, '_column')]]; split_levels <- unlist(s[[paste0(kind, '_levels')]])
      d <- composition_roe_table(part, category, split, levels, split_levels)
      d$scope <- scope; d$split_by <- kind
      d$denominator <- if (main) 'all_retained_cells' else 'cells_in_analysis_container'
      scope_tables[[kind]] <- d
      p <- composition_roe_plots(d, levels, split_levels, palette, title, kind)
      for (name in names(p)) {
        stem <- paste0(name, '_by_', kind)
        path <- file.path(out, stem)
        width <- if (name == 'ROE') if (kind == 'sample') 8.4 else 6.8 else if (kind == 'sample') 9.8 else 8.4
        height <- if (name == 'ROE') max(4.3, length(levels) * .38 + 1.6) else 5.6
        stopifnot(!file.exists(paste0(path, '.pdf')), !file.exists(paste0(path, '.png')))
        save_plot(p[[name]], path, width = width, height = height)
        figures[[length(figures) + 1L]] <- data.frame(scope = scope, kind = kind, plot = name, path = path)
      }
      overview <- wrap_plots(p$Counts + theme(legend.position = 'none'), p$Proportions, p$ROE,
        nrow = 1, widths = c(1, 1.5, 1.3)) + plot_annotation(title = paste0(title, ' - by ', kind))
      print(overview)
    }
    combined <- do.call(rbind, scope_tables); rownames(combined) <- NULL
    write.table(combined, file.path(out, 'Composition_ROE.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    tables[[scope]] <- combined
    log_msg('COMPOSITION_ROE_OK scope=', scope, '; cells=', nrow(part), '; labels=', length(levels))
  }
  saveRDS(list(metadata = meta, tables = tables, sample_groups = sample_groups, settings = s,
    input_hashes = do.call(rbind, hashes), palette = palette_table), file.path(folder, 'composition_evidence.rds'))
  write.table(do.call(rbind, hashes), file.path(folder, 'input_hashes.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, figures), file.path(folder, 'figure_manifest.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  capture.output(table(meta[[s$sample_column]], meta[[s$major_column]]), file = file.path(folder, 'major_counts.txt'))
  log_msg('ALL_COMPOSITION_ROE_OK cells=', nrow(meta), '; plots=', length(figures), '; overview_pages=14')
}
