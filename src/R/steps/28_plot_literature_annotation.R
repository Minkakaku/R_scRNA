# 按细胞名称放置一个标签；复用同一注释的cluster共享名称和颜色。
cellname_annotation_umap <- function(sc, title) {
  umap <- Embeddings(sc, 'umap')
  labels <- levels(droplevels(sc$subtype_annotation))
  centers <- aggregate(data.frame(x = umap[,1], y = umap[,2]),
    list(annotation = as.character(sc$subtype_annotation)), median)
  p <- DimPlot(sc, reduction = 'umap', group.by = 'subtype_annotation',
    cols = npg_palette(labels), raster = FALSE, shuffle = TRUE, seed = SEED, combine = FALSE)[[1]] +
    theme_classic(base_size = 11) + theme(legend.text = element_text(size = 9), plot.subtitle = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    ggrepel::geom_text_repel(data = centers, aes(x, y, label = annotation),
      inherit.aes = FALSE, size = 3.4, max.overlaps = Inf, seed = SEED,
      box.padding = .5, min.segment.length = 0, segment.color = '#777777') +
    labs(title = title)
  stopifnot(setequal(centers$annotation, labels), !anyDuplicated(centers$annotation))
  p
}

# 仅更新UMAP_annotated及其合辑；使用已核验的轻量绘图缓存，不读写分析QS。
step_update_annotated_umap_labels <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  s <- CONFIG$literature_annotation
  evidence <- file.path(PROJECT_DIR, s$evidence_dir)
  output <- file.path(PROJECT_DIR, CONFIG$current_results$subgroups)
  manifest <- read.delim(file.path(evidence, 'subgroup_annotation_summary.tsv'))
  plots <- list()
  for (group in manifest$subgroup) {
    input <- file.path(evidence, group, 'annotation_plot_data.rds')
    before <- digest::digest(file = input, algo = 'sha256')
    cached <- readRDS(input); sc <- cached$plot_object; umap <- Embeddings(sc, 'umap')
    stopifnot(identical(levels(sc$subtype_annotation), unique(cached$mapping$annotation)))
    p <- cellname_annotation_umap(sc, gsub('_', ' ', group))
    save_umap_plot(p, file.path(output, group, 'UMAP_annotated'), width = 10, height = 6,
      xlim = range(umap[,1]), ylim = range(umap[,2]))
    plots[[group]] <- standardize_umap_plot(p, range(umap[,1]), range(umap[,2]))
    stopifnot(identical(before, digest::digest(file = input, algo = 'sha256')))
    log_msg('CELLNAME_UMAP_OK ', group, '; labels=', length(unique(sc$subtype_annotation)), '; input_unchanged=TRUE')
  }
  grDevices::pdf(file.path(output, 'All_subgroups_UMAP_annotated.pdf'), width = 10, height = 6)
  for (p in plots) print(p)
  grDevices::dev.off()
  log_msg('ALL_CELLNAME_UMAPS_OK')
}

# 注释图复用第24步的坐标及全细胞表达缓存；按当前cluster展示可复用的生物标签。
step_plot_literature_annotation <- function(settings = CONFIG$literature_annotation) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- settings; evidence <- file.path(PROJECT_DIR, s$evidence_dir)
  output <- file.path(PROJECT_DIR, CONFIG$current_results$subgroups)
  manifest <- read.delim(file.path(evidence, 'subgroup_annotation_summary.tsv'))
  panels <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  overviews <- dots <- palette_records <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; folder <- file.path(output, group)
    cached <- readRDS(file.path(evidence, group, 'annotation_plot_data.rds'))
    sc <- cached$plot_object; m <- cached$mapping
    umap <- Embeddings(sc, 'umap'); xlim <- range(umap[,1]); ylim <- range(umap[,2])
    labels <- unique(m$annotation); palette <- npg_palette(labels)
    sc$subtype_annotation <- factor(as.character(sc$subtype_annotation), levels = labels)
    palette_records[[group]] <- data.frame(subgroup = group, label = labels, color = unname(palette))
    make_umap <- function(...) DimPlot(sc, reduction = 'umap', group.by = 'subtype_annotation',
      cols = palette, raster = FALSE, shuffle = TRUE, seed = SEED, ...) &
      theme_classic(base_size = 11) & theme(legend.text = element_text(size = 9), plot.subtitle = element_blank()) &
      guides(colour = guide_legend(override.aes = list(size = 3)))
    p <- cellname_annotation_umap(sc, gsub('_', ' ', group))
    save_umap_plot(p, file.path(folder, 'UMAP_annotated'), width = 10, height = 6, xlim = xlim, ylim = ylim)
    overviews[[group]] <- standardize_umap_plot(p, xlim, ylim)
    for (column in c('sample', 'group')) {
      p <- make_umap(split.by = column, ncol = if (column == 'sample') 3 else 2)
      save_umap_plot(p, file.path(folder, paste0('UMAP_', column, '_split')),
        width = if (column == 'sample') 17 else 13, height = if (column == 'sample') 10 else 6,
        xlim = xlim, ylim = ylim)
    }
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    m$candidate_annotation <- m$annotation
    genes <- strsplit(panels$genes[panels$subgroup == group], ';', fixed = TRUE)[[1]]
    p <- candidate_dotplot(e, m, genes, gsub('_', ' ', group))
    reference <- file.path(evidence, group, 'Seurat_DotPlot_reference.tsv')
    if (file.exists(reference)) {
      direct <- read.delim(reference)
      ix <- match(paste(p$data$cluster, p$data$gene), paste(direct$id, direct$features.plot))
      stopifnot(!anyNA(ix), max(abs(p$data$percent_expressed - direct$pct.exp[ix])) < 1e-8,
        max(abs(p$data$average_expression - direct$avg.exp[ix])) < 1e-8,
        max(abs(p$data$color_value - direct$avg.exp.scaled[ix])) < 1e-8)
    }
    save_plot(p$plot, file.path(folder, 'DotPlot_annotation_markers'),
      width = max(12, length(genes) * .23 + 5), height = max(4.5, nrow(m) * .3 + 2.1))
    dots[[group]] <- p$plot
    write.table(p$data, file.path(evidence, group, 'DotPlot_annotation_data.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    if (group == s$focus_subgroup) {
      genes <- strsplit(panels$genes[panels$subgroup == 'Neutrophil_function'], ';', fixed = TRUE)[[1]]
      f <- candidate_dotplot(e, m, genes, 'Neutrophil functional markers')
      save_plot(f$plot, file.path(folder, 'DotPlot_functional_markers'), width = 14, height = 4.8)
      write.table(f$data, file.path(evidence, 'neutrophil_functional_DotPlot_data.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
      levels_function <- unique(m$neutrophil_function)
      colors <- npg_palette(levels_function)
      sc$neutrophil_function <- factor(sc$neutrophil_function, levels = levels_function)
      centers <- aggregate(data.frame(x = umap[,1], y = umap[,2]),
        list(annotation = as.character(sc$subtype_annotation)), median)
      f <- DimPlot(sc, group.by = 'neutrophil_function', cols = colors, raster = FALSE, shuffle = TRUE, seed = SEED) +
        theme_classic(base_size = 11) + labs(title = 'Literature-inferred tendencies') +
        ggrepel::geom_text_repel(data = centers, aes(x, y, label = annotation), inherit.aes = FALSE,
          size = 3.5, max.overlaps = Inf, seed = SEED)
      save_umap_plot(f, file.path(folder, 'UMAP_functional_candidates'), width = 11, height = 6, xlim = xlim, ylim = ylim)
      programs <- read.delim(file.path(PROJECT_DIR, s$programs))
      d <- read.delim(file.path(evidence, 'neutrophil_cluster.tsv'))
      d$program <- factor(d$program, levels = programs$program)
      row_labels <- setNames(m$annotation, m$subgroup_cluster)
      lim <- max(abs(d$mean_z))
      f <- ggplot(d, aes(program, factor(cluster, levels = rev(m$subgroup_cluster)), fill = mean_z)) +
        geom_tile(color = 'white', linewidth = .25) +
        scale_fill_gradient2(low = '#3C5488', mid = '#F7F7F7', high = '#E64B35', midpoint = 0, limits = c(-lim, lim),
          guide = guide_colorbar(display = 'rectangles')) +
        scale_y_discrete(labels = row_labels) + labs(x = NULL, y = NULL, fill = 'Mean gene z') + theme_bw(base_size = 10) +
        theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1, color = 'black'),
          axis.text.y = element_text(color = 'black'), legend.position = 'right')
      save_plot(f, file.path(folder, 'Heatmap_literature_programs'), width = 14, height = 5.6)
      d <- read.delim(file.path(evidence, 'neutrophil_cluster_sample.tsv'))
      d <- d[d$program %in% c('CCL3 aged TAN core', 'Interferon signaling', 'NeuMap IS-I core'), ]
      d$annotation <- factor(m$annotation[match(d$cluster, m$subgroup_cluster)], levels = labels)
      d$cluster <- factor(d$cluster, levels = m$subgroup_cluster)
      d$sample_coverage <- factor(ifelse(d$cells >= 20, '>=20 cells', '<20 cells'), levels = c('>=20 cells', '<20 cells'))
      f <- ggplot(d, aes(cluster, mean_z, color = annotation)) + geom_hline(yintercept = 0, color = '#CCCCCC', linewidth = .3) +
        geom_point(aes(shape = sample_coverage), position = position_jitter(width = .15, height = 0, seed = SEED), size = 2.8) +
        scale_color_manual(values = palette) + scale_shape_manual(values = c(16, 1)) + facet_wrap(~program, nrow = 1) +
        scale_x_discrete(labels = row_labels) +
        labs(x = NULL, y = 'Sample mean gene z', color = NULL, shape = NULL) + theme_classic(base_size = 10) +
        theme(legend.position = 'bottom', legend.box = 'vertical',
          axis.text.x = element_text(angle = 75, hjust = 1, size = 8))
      save_plot(f, file.path(folder, 'Sample_literature_programs'), width = 17, height = 7)
    }
    log_msg('ANNOTATION_FIGURES_OK ', group)
    rm(sc, e, cached, p); invisible(gc())
  }
  # 合辑为每大类独立页面，避免把49个标签塞入一个缩小面板。
  grDevices::pdf(file.path(output, 'All_subgroups_UMAP_annotated.pdf'), width = 10, height = 6)
  for (p in overviews) print(p)
  grDevices::dev.off()
  grDevices::pdf(file.path(output, 'All_subgroups_annotated_marker_DotPlots.pdf'), width = 16, height = 6.5)
  for (p in dots) print(p)
  grDevices::dev.off()
  write.table(do.call(rbind, palette_records), file.path(evidence, 'annotation_colors.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  mapping <- read.delim(file.path(evidence, 'cluster_annotations.tsv'), encoding = 'UTF-8')
  write.table(mapping, file.path(output, 'Cluster_annotations.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  log_msg('ALL_LITERATURE_ANNOTATION_FIGURES_OK')
}
