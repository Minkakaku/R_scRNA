# 从已确认对象重绘；不调用标准化、聚类、marker 检验或评分，不重存输入。
refresh_levels <- function(x) {
  if (is.factor(x)) levels(droplevels(x)) else sort_cluster_ids(x)
}

# 仅保留绘图所需的 SCT assay，释放 integrated 的大矩阵；只改内存副本。
read_figure_object <- function(path) {
  obj <- read_qs2(path)
  if ("SCT" %in% names(obj@assays)) {
    DefaultAssay(obj) <- "SCT"
    obj@assays <- obj@assays["SCT"]
    obj@assays$SCT@scale.data <- matrix(numeric(), 0, 0)
  }
  gc(verbose = FALSE)
  obj
}

# 图名由项目配置传入；同一对象的所有 split 使用完整坐标范围和同一标签顺序。
refresh_umaps <- function(obj, specs) {
  xy <- Embeddings(obj, "umap")
  for (spec in specs) {
    obj[[spec$group]] <- factor(obj[[spec$group, drop = TRUE]], levels = refresh_levels(obj[[spec$group, drop = TRUE]]))
    args <- list(object = obj, group.by = spec$group, label = isTRUE(spec$label), repel = TRUE)
    if (!is.null(spec$split)) args <- c(args, list(split.by = spec$split, ncol = spec$ncol))
    if (!is.null(spec$point_size)) args$pt.size <- spec$point_size
    p <- do.call(npg_dimplot, args)
    if (!is.null(spec$legend_labels)) {
      p <- p & scale_color_manual(values = npg_palette(levels(obj[[spec$group, drop = TRUE]])),
        labels = spec$legend_labels, drop = FALSE)
    }
    p <- p & labs(title = NULL, color = NULL, x = "UMAP 1", y = "UMAP 2") & theme(legend.position = "right") &
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    save_umap_plot(p, spec$stem, spec$width, spec$height, range(xy[, 1]), range(xy[, 2]))
    log_msg("REFRESHED ", spec$stem)
  }
}

# marker 集合保留配置顺序，同一个基因只出现一次，表达量沿用 SCT data。
refresh_marker_dotplot <- function(obj, markers, group_col, stem, width = 16, height = 7) {
  genes <- unlist(markers, use.names = FALSE)
  groups <- rep(names(markers), lengths(markers))
  keep <- !duplicated(genes) & genes %in% rownames(obj[["SCT"]])
  markers <- split(genes[keep], factor(groups[keep], levels = names(markers)))
  markers <- markers[lengths(markers) > 0]
  stopifnot(length(markers) > 0)
  p <- hongf_marker_dotplot(obj, features = markers, assay = "SCT", group.by = group_col)
  save_plot(p, stem, max(width, length(unique(unlist(markers))) * 0.34 + 3),
    max(height, length(unique(obj[[group_col, drop = TRUE]])) * 0.32 + 2))
  log_msg("REFRESHED ", stem)
}

# 组间比例先逐样本计算，再等权平均；ROE 使用 observed / expected，无推断性检验。
step_refresh_reference_dotplots <- function(cfg, target = "all") {
  if (target %in% c("all", "global")) {
    obj <- read_figure_object(cfg$global$path)
    original <- digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))
    for (spec in cfg$global$dotplots) refresh_marker_dotplot(obj, spec$markers, spec$group,
      spec$stem, spec$width, spec$height)
    stopifnot(identical(original, digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))))
    rm(obj); gc(verbose = FALSE)
  }
  for (key in names(cfg$subclusters)) if (target %in% c("all", key)) {
    spec <- cfg$subclusters[[key]]
    obj <- read_figure_object(spec$annotated_qs2)
    original <- digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))
    mapping <- read.csv(spec$mapping_path, stringsAsFactors = FALSE, check.names = FALSE)
    markers <- setNames(strsplit(mapping$key_markers, ";", fixed = TRUE), mapping$plot_label)
    refresh_marker_dotplot(obj, markers, paste0(spec$metadata_prefix, "_annotation"),
      file.path(spec$output_dir, paste0(spec$plot_prefix, "_marker_DotPlot_by_annotation")),
      width = 11, height = max(3.6, nrow(mapping) * 0.32 + 2))
    stopifnot(identical(original, digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))))
    rm(obj); gc(verbose = FALSE)
  }
  if (target %in% c("all", "AUCell")) refresh_aucell_figures(cfg$aucell, dotplots_only = TRUE)
  log_msg("REFERENCE_DOTPLOTS_OK: ", target, "; metadata and UMAP unchanged")
}

refresh_composition <- function(metadata, annotation_col, output_dir, prefix, group_levels,
                                sample_col = "samples", group_col = "group", roe = TRUE, stacked = TRUE) {
  d <- data.frame(Sample = metadata[[sample_col]], Group = metadata[[group_col]],
                  CellType = metadata[[annotation_col]])
  stopifnot(!anyNA(d), all(d$Group %in% group_levels),
    all(tapply(as.character(d$Group), d$Sample, function(x) length(unique(x))) == 1L))
  d$Sample <- factor(d$Sample, levels = refresh_levels(d$Sample))
  d$CellType <- factor(d$CellType, levels = refresh_levels(d$CellType))
  palette <- npg_palette(levels(d$CellType))
  counts <- as.data.frame(table(Sample = d$Sample, CellType = d$CellType))
  names(counts)[3] <- "Number"
  counts$Proportion <- ave(counts$Number, counts$Sample, FUN = function(x) x / sum(x))
  counts$Group <- factor(d$Group[match(counts$Sample, d$Sample)], levels = group_levels)
  means <- aggregate(Proportion ~ Group + CellType, counts, mean)
  totals <- aggregate(Number ~ Group + CellType, counts, sum)
  by_group <- merge(means, totals, by = c("Group", "CellType"), sort = FALSE)
  stopifnot(sum(counts$Number) == nrow(metadata),
    all(abs(tapply(counts$Proportion, counts$Sample, sum) - 1) < 1e-10),
    all(abs(tapply(by_group$Proportion, by_group$Group, sum) - 1) < 1e-10))
  if (stacked) for (scope in c("sample", "group")) for (kind in c("proportion", "count")) {
    z <- if (scope == "sample") counts else by_group
    z$X <- z[[if (scope == "sample") "Sample" else "Group"]]
    z$Value <- z[[if (kind == "proportion") "Proportion" else "Number"]]
    p <- ggplot(z, aes(X, Value, fill = CellType)) + geom_col(width = 0.75) +
      scale_fill_manual(values = palette, drop = FALSE) +
      scale_y_continuous(labels = if (kind == "proportion") scales::percent else scales::comma) +
      labs(x = NULL, y = if (kind == "count") "Cell count" else if (scope == "group") "Mean sample proportion" else "Cell proportion", fill = NULL) +
      theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.position = "right",
        axis.text.x = element_text(angle = if (scope == "sample") 45 else 0, hjust = if (scope == "sample") 1 else 0.5))
    save_plot(p, file.path(output_dir, paste0(prefix, "_", kind, "_by_", scope)),
      if (scope == "sample") 10 else 7, 6.5)
  }
  if (roe) {
    observed <- table(factor(d$Group, levels = group_levels), d$CellType)
    expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
    stopifnot(all(expected > 0), sum(observed) == nrow(metadata),
      max(abs(rowSums(expected) - rowSums(observed))) < 1e-8,
      max(abs(colSums(expected) - colSums(observed))) < 1e-8)
    z <- as.data.frame(as.table(observed / expected))
    names(z) <- c("Group", "CellType", "ROE")
    p <- ggplot(z, aes(Group, CellType, fill = ROE)) + geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", ROE)), size = 3.5) +
      scale_fill_gradient2(low = HONGF_NPG_DIVERGING[1], mid = HONGF_NPG_DIVERGING[2],
        high = HONGF_NPG_DIVERGING[3], midpoint = 1, guide = guide_colorbar(display = "rectangles")) +
      labs(x = NULL, y = NULL, fill = "ROE") + coord_fixed() + theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank())
    save_plot(p, file.path(output_dir, paste0(prefix, "_ROE_by_group")), 7, 6.5)
    log_msg("ROE_OK: ", ncol(observed), " cell types; ", nrow(observed), " groups; O/E from pooled cell counts")
    print(observed)
    print(round(observed / expected, 3))
  }
  log_msg("COMPOSITION_OK ", prefix, ": ", nrow(metadata), " cells; ", nlevels(d$Sample), " samples")
}

# 仅补充组间 ROE；不重绘堆积图，不写入对象或额外支撑表。
step_refresh_group_roe <- function(cfg) {
  obj <- read_figure_object(cfg$path)
  stopifnot(ncol(obj) == cfg$expected_cells)
  refresh_composition(obj[[]], cfg$annotation_col, cfg$proportion_dir, cfg$proportion_prefix,
    cfg$group_levels, roe = TRUE, stacked = FALSE)
}

# QC 保留原始过滤前、常规过滤后和去污染后三个阶段，避免混淆对象。
refresh_qc_stage <- function(cfg) {
  obj <- read_qs2(cfg$path)
  stopifnot(ncol(obj) == cfg$expected_cells)
  if (isTRUE(cfg$contamination_only)) {
    p <- make_qc_plots(obj)$contamination +
      scale_y_continuous(limits = c(0, cfg$cutoff * 1.05)) +
      geom_hline(yintercept = cfg$cutoff, linetype = 2, color = HONGF_NPG[1])
    save_plot(p, cfg$stem, 10, 5)
  } else save_qc_plots(obj, cfg$stage)
  log_msg("QC_OK ", cfg$stage, ": ", ncol(obj), " cells")
}

refresh_global_figures <- function(cfg) {
  obj <- read_figure_object(cfg$path)
  original <- digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))
  stopifnot(ncol(obj) == cfg$expected_cells)
  save_pca_elbow(obj, cfg$pca_stem, save_tables = FALSE, selected_pc = cfg$pcs)
  refresh_umaps(obj, cfg$umaps)
  for (spec in cfg$dotplots) refresh_marker_dotplot(obj, spec$markers, spec$group, spec$stem, spec$width, spec$height)
  refresh_composition(obj[[]], cfg$annotation_col, cfg$proportion_dir, cfg$proportion_prefix,
    cfg$group_levels, roe = FALSE)
  stopifnot(identical(original, digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))))
  log_msg("GLOBAL_OK: metadata and coordinates unchanged")
}

refresh_subcluster_figures <- function(cfg, group_levels) {
  obj <- read_figure_object(cfg$annotated_qs2)
  original <- digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))
  mapping <- read.csv(cfg$mapping_path, stringsAsFactors = FALSE, check.names = FALSE)
  annotation_col <- paste0(cfg$metadata_prefix, "_annotation")
  stopifnot(setequal(as.character(obj[[annotation_col, drop = TRUE]]), mapping$annotation))
  save_pca_elbow(obj, file.path(cfg$output_dir, "PCA_elbow"), save_tables = FALSE, selected_pc = cfg$pcs)
  legend_labels <- setNames(mapping$annotation, mapping$plot_label)
  plot_obj <- obj
  plot_obj$refresh_label <- factor(mapping$plot_label[match(as.character(obj[[annotation_col, drop = TRUE]]), mapping$annotation)],
    levels = mapping$plot_label)
  stem <- function(suffix) file.path(cfg$output_dir, paste0("UMAP_", cfg$plot_prefix, "_", suffix))
  specs <- list(
    list(group = cfg$cluster_col, label = TRUE, stem = stem("clusters"), width = 10, height = 7.5),
    list(group = "samples", stem = stem("by_sample"), width = 10, height = 7.5),
    list(group = "refresh_label", label = TRUE, legend_labels = legend_labels, stem = stem("annotation"), width = 12, height = 8),
    list(group = "refresh_label", split = "samples", ncol = 3, legend_labels = legend_labels,
      stem = stem("annotation_split_by_sample"), width = 15, height = 10),
    list(group = "refresh_label", split = "group", ncol = 2, legend_labels = legend_labels,
      stem = stem("annotation_split_by_group"), width = 12, height = 8))
  refresh_umaps(plot_obj, specs)
  rm(plot_obj)
  markers <- setNames(strsplit(mapping$key_markers, ";", fixed = TRUE), mapping$plot_label)
  refresh_marker_dotplot(obj, markers, annotation_col,
    file.path(cfg$output_dir, paste0(cfg$plot_prefix, "_marker_DotPlot_by_annotation")),
    width = max(11, min(19, length(unique(unlist(markers))) * 0.25 + 3)), height = max(3.4, nrow(mapping) * 0.28 + 2))
  refresh_composition(obj[[]], annotation_col, cfg$proportion_dir, paste0(cfg$plot_prefix, "_subtype"), group_levels)
  stopifnot(identical(original, digest::digest(list(obj@meta.data, Embeddings(obj, "umap")))))
  log_msg("SUBCLUSTER_OK ", cfg$label, ": ", ncol(obj), " cells; PC1-", cfg$pcs, "; metadata and coordinates unchanged")
}

# AUCell 仅从缓存读取样本均值；组差固定为 treatment - reference。
refresh_aucell_figures <- function(cfg, dotplots_only = FALSE) {
  res <- read_qs2(cfg$output_qs2)
  format_data <- function(d) {
    d$pathway <- factor(d$pathway, levels = cfg$pathway_labels)
    d$subtype <- factor(d$subtype, levels = rev(res$subtype_levels))
    d
  }
  mean_data <- format_data(res$by_subtype)
  group_data <- format_data(res$by_group)
  group_data$group <- factor(group_data$group, levels = cfg$group_levels)
  upper <- max(group_data$score, mean_data$score)
  dot <- function(d) hongf_dotplot_style(ggplot(d, aes(pathway, subtype, size = score, color = score)) + geom_point()) +
    scale_size_continuous(name = "Mean AUCell", range = c(0.5, 8), limits = c(0, upper)) +
    scale_color_gradientn(name = "Mean AUCell", colours = HONGF_DOTPLOT_COLORS, limits = c(0, upper),
      guide = guide_colorbar(display = "rectangles")) + labs(x = NULL, y = NULL) + theme(legend.position = "bottom")
  save_plot(dot(mean_data), file.path(cfg$output_dir, paste0(cfg$prefix, "_DotPlot_by_subtype")), 12.5, 7.5)
  save_plot(dot(group_data) + facet_grid(group ~ .),
    file.path(cfg$output_dir, paste0(cfg$prefix, "_DotPlot_by_subtype_and_group")), 12.5, 10.5)
  if (dotplots_only) return(invisible(NULL))
  mean_data$z_score <- ave(mean_data$score, mean_data$pathway, FUN = function(x) if (sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x)))
  delta <- merge(subset(res$by_group, group == cfg$group_levels[1]), subset(res$by_group, group == cfg$group_levels[2]),
    by = c("subtype", "gene_set", "pathway"), suffixes = c("_reference", "_treatment"))
  delta$delta <- delta$score_treatment - delta$score_reference
  delta <- format_data(delta)
  heatmap <- function(d, value, label) ggplot(d, aes(pathway, subtype, fill = .data[[value]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = HONGF_NPG_DIVERGING[1], mid = HONGF_NPG_DIVERGING[2], high = HONGF_NPG_DIVERGING[3],
      midpoint = 0, name = label, guide = guide_colorbar(display = "rectangles")) +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_plot(heatmap(mean_data, "z_score", "AUCell\npathway Z-score"),
    file.path(cfg$output_dir, paste0(cfg$prefix, "_Heatmap_by_subtype")), 13.5, 6.5)
  save_plot(heatmap(delta, "delta", paste0("Mean AUCell\n", cfg$group_levels[2], " - ", cfg$group_levels[1])),
    file.path(cfg$output_dir, paste0(cfg$prefix, "_Heatmap_group_delta")), 13.5, 6.5)
  log_msg("AUCELL_OK: reused cached scores and sample means")
}

# PCD 小提琴沿用缓存 UCell 分数与原始样本级 p/q；细胞仅决定分布形状。
refresh_pcd_violins <- function(cfg) {
  res <- readRDS(file.path(cfg$model_dir, "PCD_results.rds"))
  mode_ids <- res$reference$modes$mode
  pages <- split(mode_ids, ceiling(seq_along(mode_ids) / cfg$violin_modes_per_page))
  for (key in names(cfg$sources)) {
    scored <- read_qs2(file.path(cfg$model_dir, paste0("scores_", key, ".qs2")))
    for (page in seq_along(pages)) {
      ids <- pages[[page]]
      d <- do.call(rbind, lapply(ids, function(id) transform(scored$metadata, mode = id,
        score = as.numeric(scored$scores$UCell[id, ]))))
      sample <- subset(res$samples, lineage == key & method == "UCell" & mode %in% ids)
      stats <- subset(res$stats, lineage == key & method == "UCell" & mode %in% ids)
      stats$annotation <- ifelse(stats$testable, sprintf("p=%.2f\nq=%.2f", stats$p_value, stats$fdr), "low n")
      high <- tapply(d$score, d$mode, max); low <- tapply(d$score, d$mode, min)
      stats$y <- high[stats$mode] + 0.09 * pmax(high[stats$mode] - low[stats$mode], 0.01)
      format_data <- function(z) {
        z$subtype <- factor(z$subtype, levels = scored$subtype_levels)
        z$mode <- factor(z$mode, levels = ids, labels = res$reference$modes$label[match(ids, mode_ids)])
        if ("group" %in% names(z)) z$group <- factor(z$group, levels = cfg$group_levels)
        z
      }
      d <- format_data(d); sample <- format_data(sample); stats <- format_data(stats)
      sample$qualified <- sample$n_cells >= cfg$min_sample_cells
      p <- ggplot(d, aes(subtype, score, fill = group)) +
        geom_violin(position = position_dodge(0.82), width = 0.8, scale = "width", trim = TRUE,
          color = "#555555", linewidth = 0.25, alpha = 0.8) +
        geom_point(data = sample, aes(shape = qualified, group = group), size = 1.7, stroke = 0.5, color = "black",
          position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.82, seed = cfg$seed)) +
        geom_text(data = stats, aes(subtype, y, label = annotation), inherit.aes = FALSE, size = 2.4, lineheight = 0.95) +
        scale_fill_manual(values = cfg$group_colors, name = NULL) +
        scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 4), name = "Mouse mean", labels = c(`FALSE` = "<10 cells", `TRUE` = ">=10 cells")) +
        facet_grid(mode ~ ., scales = "free_y") + scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
        labs(x = NULL, y = "UCell association score", title = paste0(cfg$sources[[key]]$label, " | PCD programs ", page, "/", length(pages)),
          caption = sprintf("Violins: cells; points: mouse means. p: exact two-sided sample-label permutation; q: BH across subtypes and programs.\nLow n: fewer than %d mice with >=%d cells in either group. PCD-associated expression does not establish cell death.", cfg$min_samples_per_group, cfg$min_sample_cells)) +
        theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
          strip.background = element_rect(fill = "#F2F2F2"), strip.text.y = element_text(angle = 0, size = 9),
          plot.title = element_text(face = "bold"), legend.position = "top", plot.caption = element_text(hjust = 0, size = 8.5))
      save_plot(p, file.path(cfg$sources[[key]]$output_dir, "PCD", paste0("Fig5_", key, "_PCD_violin_", sprintf("%02d", page))), 15, 12.5)
    }
    log_msg("PCD_OK ", key, ": cached scores; unchanged sample-level p/q")
  }
}

# 单独重绘双细胞比例，复用统一的类别配色与导出方法。
refresh_qc_doublet_figures <- function(stages) {
  for (cfg in stages) {
    obj <- read_qs2(cfg$path)
    metadata <- obj[[]]
    rm(obj); gc(verbose = FALSE)
    p <- plot_qc_doublet_ratio(metadata)
    built <- ggplot_build(p)
    stopifnot(identical(built$plot$scales$get_scales("fill")$map(c("singlet", "doublet")),
      c("#4DBBD5", "#E64B35")))
    save_plot(p, file.path(QC_DIR, paste0("QC_doublet_ratio_", cfg$stage)), 10, 5)
    log_msg("QC_DOUBLET_COLORS_OK ", cfg$stage, ": singlet=#4DBBD5; doublet=#E64B35")
  }
}

# 每个任务可单独重绘，以便大对象逐一释放内存。
step_refresh_existing_figures <- function(cfg, task = "all") {
  choices <- c(names(cfg$qc), "global", names(cfg$subclusters), "AUCell", "PCD")
  stopifnot(task %in% c("all", choices))
  tasks <- if (task == "all") choices else task
  for (key in tasks) {
    log_msg("REFRESH_START ", key)
    if (key %in% names(cfg$qc)) refresh_qc_stage(cfg$qc[[key]])
    if (key == "global") refresh_global_figures(cfg$global)
    if (key %in% names(cfg$subclusters)) refresh_subcluster_figures(cfg$subclusters[[key]], cfg$group_levels)
    if (key == "AUCell") refresh_aucell_figures(cfg$aucell)
    if (key == "PCD") refresh_pcd_violins(cfg$pcd)
    gc(verbose = FALSE)
    log_msg("REFRESH_DONE ", key)
  }
  invisible(NULL)
}
