# 文献程序 AUCell 与分块热图；项目的基因集选择、路径和标签全部由 cfg 传入。
prepare_program_reference <- function(cfg) {
  mouse <- readRDS(cfg$mouse_msig)
  if (!file.exists(cfg$human_msig)) {
    human <- readRDS(cfg$human_msig_cache)
    human <- human[human$gs_name %in% cfg$human_sets, ]
    dir.create(dirname(cfg$human_msig), recursive = TRUE, showWarnings = FALSE)
    saveRDS(human, cfg$human_msig)
  }
  human <- readRDS(cfg$human_msig)
  selected_mouse <- setdiff(cfg$rows$id, c(cfg$human_sets, cfg$factors$id))
  stopifnot(all(selected_mouse %in% mouse$gs_name), all(cfg$human_sets %in% human$gs_name))
  pieces <- list()
  for (id in selected_mouse) {
    d <- mouse[mouse$gs_name == id, ]
    pieces[[id]] <- data.frame(id, source_gene = unique(d$gene_symbol), human_symbol = NA_character_,
      mouse_symbol = unique(d$gene_symbol), mapping = "MSigDB_mouse", source = unique(d$db_version),
      source_detail = unique(d$gs_url), weight = NA_real_)
  }
  # HGNC 精确现名优先，旧名/别名只采用唯一匹配；MGI 多小鼠同源基因不任意择一。
  hgnc <- read.delim(cfg$hgnc, check.names = FALSE, stringsAsFactors = FALSE)
  hgnc <- hgnc[hgnc$status == "Approved", ]
  aliases <- lapply(seq_len(nrow(hgnc)), function(i)
    unique(strsplit(paste(hgnc$prev_symbol[i], hgnc$alias_symbol[i], sep = "|"), "|", fixed = TRUE)[[1]]))
  alias_map <- split(rep(hgnc$symbol, lengths(aliases)), unlist(aliases))
  mgi <- read.delim(cfg$mgi, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  mm <- split(mgi$Symbol[mgi$`NCBI Taxon ID` == 10090], mgi$`DB Class Key`[mgi$`NCBI Taxon ID` == 10090])
  hs <- mgi[mgi$`NCBI Taxon ID` == 9606, ]
  map_human <- function(original) {
    symbol <- if (original %in% hgnc$symbol) original else {
      candidates <- unique(alias_map[[original]])
      if (length(candidates) == 1L) candidates else NA_character_
    }
    candidates <- if (is.na(symbol)) character() else
      unique(unlist(mm[as.character(hs$`DB Class Key`[hs$Symbol == symbol])], use.names = FALSE))
    data.frame(source_gene = original, human_symbol = symbol,
      mouse_symbol = if (length(candidates) == 1L) candidates else NA_character_,
      mapping = if (is.na(symbol)) "unresolved_HGNC" else if (length(candidates) == 0L)
        "no_MGI_mouse" else if (length(candidates) > 1L) "ambiguous_MGI_mouse" else "MGI_single_mouse")
  }
  for (id in cfg$human_sets) {
    d <- human[human$gs_name == id, ]
    out <- do.call(rbind, lapply(unique(d$db_gene_symbol), map_human))
    pieces[[id]] <- data.frame(id, out, source = unique(d$db_version),
      source_detail = paste(unique(d$gs_exact_source), collapse = "; "), weight = NA_real_)
  }
  weights <- as.data.frame(readxl::read_xlsx(cfg$weights, sheet = cfg$weights_sheet))
  for (i in seq_len(nrow(cfg$factors))) {
    f <- cfg$factors[i, ]
    d <- weights[weights$factor == f$factor & weights[[f$weight]] > 0, ]
    d <- d[order(-d[[f$weight]], d$symbol), ]
    d <- head(d[!duplicated(d$symbol), ], cfg$weights_n)
    stopifnot(nrow(d) == cfg$weights_n)
    out <- do.call(rbind, lapply(d$symbol, map_human))
    pieces[[f$id]] <- data.frame(id = f$id, out, source = cfg$weights_doi,
      source_detail = paste0(cfg$weights_sheet, ": factor ", f$factor, ", top ", cfg$weights_n,
        " positive weights; derived program, not original classifier"), weight = d[[f$weight]])
  }
  members <- do.call(rbind, pieces[cfg$rows$id]); rownames(members) <- NULL
  reference_files <- c(cfg$mouse_msig, cfg$human_msig, cfg$weights, cfg$hgnc, cfg$mgi)
  list(rows = cfg$rows, members = members, source_md5 = tools::md5sum(reference_files),
    weights_url = cfg$weights_url, config = cfg)
}

program_aucell_block <- function(expr, sets, cutoff) {
  # 沿用已验证的 AUCell 稀疏矩阵零值屏蔽，避免未检出基因进入排名得分。
  ranking <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, splitByBlocks = TRUE,
    keepZeroesAsNA = FALSE, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  ranks <- AUCell::getRanking(ranking)
  detected <- Matrix::colSums(expr > 0)
  ranks[sweep(ranks, 2, detected, ">")] <- NA_integer_
  stopifnot(all(colSums(!is.na(ranks)) == detected))
  SummarizedExperiment::assay(ranking, "ranking") <- ranks
  result <- AUCell::AUCell_calcAUC(sets, ranking, aucMaxRank = cutoff,
    normAUC = TRUE, nCores = 1, verbose = FALSE)
  AUCell::getAUC(result)
}

plot_program_heatmap <- function(result, cfg, values = result$z_score,
  legend_title = "Row Z-score", color_breaks = c(-2, -1, 0, 1, 2),
  legend_at = c(-2, 0, 2), legend_labels = c("<= -2", "0", ">= 2")) {
  z <- values
  colnames(z) <- result$annotation$plot_label[match(colnames(z), result$annotation$annotation)]
  split <- factor(cfg$rows$block, levels = names(cfg$block_colors))
  colors <- circlize::colorRamp2(color_breaks, cfg$heatmap_colors)
  block <- ComplexHeatmap::rowAnnotation(block = ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = unname(cfg$block_colors), col = NA), labels = levels(split),
    labels_rot = 90, labels_gp = grid::gpar(col = cfg$block_text_colors, fontsize = 10, fontface = "bold")),
    width = grid::unit(8, "mm"), show_annotation_name = FALSE)
  ht <- ComplexHeatmap::Heatmap(z, name = legend_title, col = colors,
    cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
    row_split = split, row_gap = grid::unit(1.7, "mm"), row_title = NULL,
    row_labels = cfg$rows$label, row_names_side = "right", row_names_gp = grid::gpar(fontsize = 10),
    row_names_max_width = grid::unit(55, "mm"), column_names_rot = 55,
    column_names_gp = grid::gpar(fontsize = 10), column_names_max_height = grid::unit(40, "mm"),
    left_annotation = block, rect_gp = grid::gpar(col = NA), border = FALSE, use_raster = FALSE,
    width = grid::unit(125, "mm"), height = grid::unit(nrow(z) * 5.5 + (nlevels(split) - 1) * 1.7, "mm"),
    heatmap_legend_param = list(at = legend_at, labels = legend_labels,
      title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(45, "mm"), color_bar = "continuous"))
  dir.create(dirname(cfg$figure), recursive = TRUE, showWarnings = FALSE)
  draw <- function() ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
    align_heatmap_legend = "heatmap_top", padding = grid::unit(c(8, 5, 6, 5), "mm"))
  grDevices::pdf(paste0(cfg$figure, ".pdf"), width = cfg$width, height = cfg$height,
    family = "Helvetica", useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  ragg::agg_png(paste0(cfg$figure, ".png"), width = cfg$width, height = cfg$height,
    units = "in", res = 600, background = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
}

run_program_heatmap <- function(cfg) {
  needed <- c("Seurat", "qs2", "AUCell", "ComplexHeatmap", "circlize", "ragg", "readxl")
  stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly = TRUE)))
  before <- file.info(cfg$input)[, c("size", "mtime")]
  reference <- prepare_program_reference(cfg)
  message("Reading annotated object; no changes will be written to it")
  obj <- qs2::qs_read(cfg$input)
  md <- obj[[]]
  meta <- data.frame(cell = rownames(md), sample = as.character(md[[cfg$sample_col]]),
    group = as.character(md[[cfg$group_col]]), subtype = as.character(md[[cfg$subtype_col]]))
  annotation <- read.csv(cfg$annotation, stringsAsFactors = FALSE)
  levels <- annotation$annotation
  stopifnot(!anyNA(meta), !anyDuplicated(meta$cell), setequal(meta$subtype, levels),
    all(vapply(split(meta$group, meta$sample), function(x) length(unique(x)), integer(1)) == 1L))
  layers <- SeuratObject::Layers(obj[[cfg$assay]], search = "^counts")
  counts <- setNames(lapply(layers, function(x) SeuratObject::LayerData(obj, assay = cfg$assay, layer = x)), layers)
  rm(obj, md); gc()
  genes <- rownames(counts[[1]])
  stopifnot(!anyDuplicated(genes), length(counts) > 0L,
    all(vapply(counts, function(x) identical(rownames(x), genes), logical(1))),
    setequal(unlist(lapply(counts, colnames)), meta$cell), !anyDuplicated(unlist(lapply(counts, colnames))))
  reference$members$present <- !is.na(reference$members$mouse_symbol) & reference$members$mouse_symbol %in% genes
  sets <- setNames(lapply(cfg$rows$id, function(id) unique(reference$members$mouse_symbol[
    reference$members$id == id & reference$members$present])), cfg$rows$id)
  coverage <- data.frame(id = cfg$rows$id,
    source_genes = as.integer(table(factor(reference$members$id, levels = cfg$rows$id))),
    mapped_mouse = vapply(cfg$rows$id, function(id) length(unique(na.omit(reference$members$mouse_symbol[reference$members$id == id]))), integer(1)),
    present = lengths(sets), row.names = NULL)
  print(coverage)
  stopifnot(all(lengths(sets) >= cfg$minimum_genes))
  cutoff <- ceiling(length(genes) * cfg$auc_fraction)
  set.seed(cfg$seed)
  scores <- matrix(NA_real_, nrow = length(sets), ncol = nrow(meta), dimnames = list(names(sets), meta$cell))
  detected <- setNames(rep(NA_integer_, nrow(meta)), meta$cell)
  for (i in seq_along(counts)) {
    expr <- counts[[i]]
    detected[colnames(expr)] <- Matrix::colSums(expr > 0)
    blocks <- split(seq_len(ncol(expr)), ceiling(seq_len(ncol(expr)) / cfg$chunk_cells))
    for (j in seq_along(blocks)) {
      message(format(Sys.time(), "%H:%M:%S"), " AUCell ", names(counts)[i], " block ", j, "/", length(blocks))
      b <- expr[, blocks[[j]], drop = FALSE]
      scores[, colnames(b)] <- program_aucell_block(b, sets, cutoff)[names(sets), colnames(b), drop = FALSE]
      rm(b); gc(verbose = FALSE)
    }
    counts[i] <- list(NULL); rm(expr); gc(verbose = FALSE)
  }
  meta$n_genes_detected <- unname(detected[meta$cell])
  stopifnot(all(is.finite(scores)), all(scores >= 0 & scores <= 1))
  sample_levels <- sort(unique(meta$sample))
  n_cells <- table(factor(meta$subtype, levels), factor(meta$sample, sample_levels))
  stopifnot(all(n_cells > 0))
  # 先对每只样本的每个亚群取细胞均值，再跨样本等权平均；不把细胞作为生物重复。
  sample_means <- array(NA_real_, c(nrow(scores), length(levels), length(sample_levels)),
    dimnames = list(rownames(scores), levels, sample_levels))
  for (s in sample_levels) for (type in levels) {
    cells <- meta$cell[meta$sample == s & meta$subtype == type]
    sample_means[, type, s] <- rowMeans(scores[, cells, drop = FALSE])
  }
  balanced <- apply(sample_means, c(1, 2), mean)
  stopifnot(all(apply(balanced, 1, sd) > 0))
  z <- t(scale(t(balanced)))
  comparison <- NULL
  if (!is.null(cfg$previous_scores)) {
    old <- qs2::qs_read(cfg$previous_scores)
    common <- intersect(names(sets), names(old$gene_sets))
    stopifnot(identical(genes, old$universe), identical(meta$cell, old$metadata$cell),
      identical(meta$subtype, old$metadata$subtype), identical(meta$sample, old$metadata$sample),
      identical(meta$group, old$metadata$group), cutoff == old$auc_max_rank,
      all(vapply(common, function(id) setequal(sets[[id]], old$gene_sets[[id]]), logical(1))))
    comparison <- data.frame(id = common,
      max_abs_difference = vapply(common, function(id) max(abs(scores[id, ] - old$scores[id, meta$cell])), numeric(1)))
    print(comparison)
    # 相同随机种子、分块、基因宇宙必须复现既有评分，非静默混用两种尺度。
    stopifnot(all(comparison$max_abs_difference < 1e-10))
  }
  result <- list(scores = scores, metadata = meta, annotation = annotation,
    sample_means = sample_means, n_cells = n_cells, balanced_mean = balanced, z_score = z,
    gene_sets = sets, reference = reference, coverage = coverage, universe = genes,
    auc_max_rank = cutoff, previous_comparison = comparison, config = cfg, session = sessionInfo())
  if (!is.null(cfg$contrast)) {
    group_names <- c(cfg$contrast$reference, cfg$contrast$test)
    group_samples <- setNames(lapply(group_names, function(g) unique(meta$sample[meta$group == g])), group_names)
    stopifnot(all(lengths(group_samples) >= 2L), !anyDuplicated(unlist(group_samples)))
    group_means <- lapply(group_samples, function(s) apply(sample_means[, , s, drop = FALSE], c(1, 2), mean))
    result$contrast <- list(reference = cfg$contrast$reference, test = cfg$contrast$test,
      samples = group_samples, group_means = group_means,
      delta = group_means[[cfg$contrast$test]] - group_means[[cfg$contrast$reference]])
    stopifnot(all(is.finite(result$contrast$delta)))
  }
  dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(result, cfg$output)
  saveRDS(reference, cfg$reference_output)
  plot_program_heatmap(result, cfg)
  if (!is.null(cfg$contrast)) {
    plot_cfg <- cfg
    plot_cfg$figure <- cfg$contrast$figure
    delta <- result$contrast$delta
    limit <- max(0.001, ceiling(max(abs(delta)) * 1000) / 1000)
    at <- c(-limit, 0, limit)
    plot_program_heatmap(result, plot_cfg, values = delta,
      legend_title = paste0(cfg$contrast$test, " - ", cfg$contrast$reference, "\nMean AUCell"),
      color_breaks = seq(-limit, limit, length.out = 5), legend_at = at,
      legend_labels = format(at, trim = TRUE, scientific = FALSE))
    cat("GROUP_CONTRAST_OK:", cfg$contrast$test, "-", cfg$contrast$reference,
      "; sample counts", lengths(result$contrast$samples), "; delta range", range(delta), "\n")
  }
  stopifnot(identical(before, file.info(cfg$input)[, c("size", "mtime")]))
  print(n_cells)
  cat("PROGRAM_HEATMAP_OK:", nrow(scores), "programs;", ncol(scores), "cells;",
      length(levels), "subtypes;", length(sample_levels), "samples\n")
  cat("Cells with fewer detected genes than AUC cutoff:", sum(detected < cutoff), "/", length(detected), "\n")
  invisible(result)
}

# 从既有每样本AUC均值计算独立组的Hedges' g；不重算单细胞评分或改写旧图。
run_program_effect_heatmap <- function(cfg) {
  r <- qs2::qs_read(cfg$input)
  protected <- c(cfg$input, paste0(r$config$figure, c(".png", ".pdf")),
    paste0(r$config$contrast$figure, c(".png", ".pdf")))
  before <- tools::md5sum(protected)
  sample_meta <- unique(r$metadata[c("sample", "group")])
  stopifnot(!anyDuplicated(sample_meta$sample), !anyNA(sample_meta),
    all(is.finite(r$sample_means)), cfg$reference != cfg$test)
  groups <- setNames(lapply(c(cfg$reference, cfg$test), function(g)
    sample_meta$sample[sample_meta$group == g]), c(cfg$reference, cfg$test))
  n <- lengths(groups)
  stopifnot(all(n >= 2L), all(unlist(groups) %in% dimnames(r$sample_means)[[3]]))
  means <- lapply(groups, function(s) apply(r$sample_means[, , s, drop = FALSE], c(1, 2), mean))
  variances <- lapply(groups, function(s) apply(r$sample_means[, , s, drop = FALSE], c(1, 2), var))
  df <- sum(n) - 2L
  pooled_sd <- sqrt(((n[[cfg$reference]] - 1) * variances[[cfg$reference]] +
    (n[[cfg$test]] - 1) * variances[[cfg$test]]) / df)
  # 方差为零时效应量不可估计；停止而非用任意epsilon制造极大数值。
  if (any(!is.finite(pooled_sd) | pooled_sd <= 0)) stop("Nonpositive pooled SD; effect cannot be estimated")
  delta <- means[[cfg$test]] - means[[cfg$reference]]
  # 正态、独立样本模型下的精确小样本修正J(df)，而非大样本近似。
  correction <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  hedges_g <- correction * delta / pooled_sd
  stopifnot(all(is.finite(hedges_g)), all(sign(hedges_g) == sign(delta)),
    identical(cfg$reference, r$contrast$reference), identical(cfg$test, r$contrast$test),
    max(abs(delta - r$contrast$delta)) < 1e-12)
  plot_cfg <- r$config
  plot_cfg$figure <- cfg$figure
  result <- list(hedges_g = hedges_g, delta = delta, pooled_sd = pooled_sd,
    group_means = means, group_variances = variances, correction = correction, df = df,
    n = n, sample_meta = sample_meta, sample_means = r$sample_means, rows = r$config$rows,
    config = cfg, source_md5 = before, session = sessionInfo())
  dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(result, cfg$output)
  limit <- cfg$display_limit
  sample_label <- if (length(unique(n)) == 1L) paste0(n[1], " mice/group") else paste(names(n), n, collapse = "; ")
  plot_program_heatmap(r, plot_cfg, values = hedges_g,
    legend_title = paste0("Hedges' g\n", cfg$test, " - ", cfg$reference, "\n", sample_label),
    color_breaks = seq(-limit, limit, length.out = 5), legend_at = c(-limit, 0, limit),
    legend_labels = c(paste0("<= -", limit), "0", paste0(">= ", limit)))
  stopifnot(identical(before, tools::md5sum(protected)))
  cat("HEDGES_G_OK: dimensions", dim(hedges_g), "; sample counts", n, "; df", df,
    "; J", correction, "; untruncated range", range(hedges_g),
    "; display-saturated cells", sum(abs(hedges_g) > limit), "\n")
  cat("INPUT_AND_PREVIOUS_FIGURE_HASHES_UNCHANGED\n")
  invisible(result)
}
