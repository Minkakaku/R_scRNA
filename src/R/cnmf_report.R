# 已完成 cNMF 模型的功能富集、样本级比较和可视化；不重新拟合或修改源对象。
cnmf_read_matrix <- function(path) {
  d <- data.table::fread(path, header = TRUE, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- d[[1]]
  storage.mode(m) <- "double"
  m
}

cnmf_exact_p <- function(reference, test) {
  x <- c(reference, test)
  differences <- apply(combn(seq_along(x), length(test)), 2, function(i) mean(x[i]) - mean(x[-i]))
  mean(abs(differences) >= abs(mean(test) - mean(reference)) - 1e-12)
}

cnmf_reference <- function(cfg) {
  h <- readRDS(cfg$hallmark)
  g <- readRDS(cfg$go); g <- g[g$gs_subcollection == "GO:BP", ]
  list(Hallmark = split(h$db_gene_symbol, h$gs_name), GO_BP = split(g$db_gene_symbol, g$gs_name),
    version = c(Hallmark = unique(h$db_version), GO_BP = unique(g$db_version)),
    md5 = tools::md5sum(c(cfg$hallmark, cfg$go)))
}

cnmf_ora <- function(scores, universe, reference, cfg) {
  rows <- list()
  for (collection in c("Hallmark", "GO_BP")) {
    sets <- lapply(reference[[collection]], intersect, y = universe)
    sets <- sets[lengths(sets) >= cfg$set_size[1] & lengths(sets) <= cfg$set_size[2]]
    for (program in rownames(scores)) {
      ranked <- sort(scores[program, universe], decreasing = TRUE)
      selected <- head(names(ranked)[ranked > 0], cfg$enrichment_top_n)
      overlap <- lapply(sets, intersect, y = selected)
      hits <- lengths(overlap); sizes <- lengths(sets)
      p <- phyper(hits - 1L, sizes, length(universe) - sizes, length(selected), lower.tail = FALSE)
      rows[[paste(collection, program)]] <- data.frame(program, collection, pathway = names(sets),
        overlap = unname(hits), set_size = unname(sizes), selected_n = length(selected),
        universe_n = length(universe), fold_enrichment = unname(hits / length(selected) / (sizes / length(universe))),
        p_value = unname(p), fdr = unname(p.adjust(p, "BH")),
        genes = vapply(overlap, paste, character(1), collapse = ";"), row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

prepare_cnmf_lineage <- function(cfg, key, reference) {
  item <- cfg$lineages[[key]]
  model <- file.path(cfg$model_dir, key)
  suffix <- sprintf("k_%d.dt_%s", item$k, cfg$density_key)
  files <- c(usage = file.path(model, paste0(key, ".usages.", suffix, ".consensus.txt")),
    scores = file.path(model, paste0(key, ".gene_spectra_score.", suffix, ".txt")),
    tpm = file.path(model, paste0(key, ".gene_spectra_tpm.", suffix, ".txt")),
    metadata = file.path(cfg$metadata_dir, paste0(key, "_metadata.csv")), object = item$object)
  usage <- cnmf_read_matrix(files[["usage"]]); scores <- cnmf_read_matrix(files[["scores"]])
  tpm <- cnmf_read_matrix(files[["tpm"]])
  stopifnot(identical(colnames(usage), rownames(scores)), identical(dimnames(scores), dimnames(tpm)),
    ncol(usage) == item$k, nrow(usage) == item$cells, !anyDuplicated(rownames(usage)),
    all(is.finite(usage)), all(usage >= 0), all(rowSums(usage) > 0))
  programs <- sprintf("GEP%02d", seq_len(item$k))
  original_program_ids <- colnames(usage)
  colnames(usage) <- rownames(scores) <- rownames(tpm) <- programs
  usage <- usage / rowSums(usage)
  meta <- read.csv(files[["metadata"]], stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(meta$cell), setequal(meta$cell, rownames(usage)))
  meta <- meta[match(rownames(usage), meta$cell), ]; rownames(meta) <- meta$cell
  object <- qs2::qs_read(item$object)
  stopifnot(setequal(colnames(object), meta$cell))
  om <- object[[]][meta$cell, ]
  stopifnot(identical(as.character(om$samples), meta$sample), identical(as.character(om$group), meta$group),
    identical(as.character(om[[item$subtype_col]]), meta$subtype))
  umap <- SeuratObject::Embeddings(object, "umap")[meta$cell, , drop = FALSE]
  rm(object); gc(FALSE)
  subtypes <- read.csv(item$annotation, stringsAsFactors = FALSE)
  stopifnot(setequal(meta$subtype, subtypes$annotation))
  sample_info <- unique(meta[, c("sample", "group")])
  sample_info <- sample_info[order(match(sample_info$group, cfg$group_levels), sample_info$sample), ]
  stopifnot(!anyDuplicated(sample_info$sample), identical(as.integer(table(factor(sample_info$group,
    levels = cfg$group_levels))), cfg$mice_per_group))
  sample_usage <- t(vapply(sample_info$sample, function(s) colMeans(usage[meta$sample == s, , drop = FALSE]),
    numeric(item$k))); dimnames(sample_usage) <- list(sample_info$sample, programs)
  subtype_usage <- t(vapply(subtypes$annotation, function(s) colMeans(usage[meta$subtype == s, , drop = FALSE]),
    numeric(item$k))); dimnames(subtype_usage) <- list(subtypes$annotation, programs)
  comparison <- do.call(rbind, lapply(programs, function(p) {
    a <- sample_usage[sample_info$group == cfg$group_levels[1], p]
    b <- sample_usage[sample_info$group == cfg$group_levels[2], p]
    data.frame(program = p, reference_mean = mean(a), test_mean = mean(b), difference = mean(b) - mean(a),
      p_value = cnmf_exact_p(a, b))
  }))
  comparison$fdr <- p.adjust(comparison$p_value, "BH")
  universe <- colnames(tpm)[colSums(tpm) > 0 & colSums(!is.finite(scores)) == 0]
  top <- do.call(rbind, lapply(programs, function(p) {
    genes <- head(names(sort(scores[p, ], decreasing = TRUE)), cfg$top_n)
    data.frame(program = p, rank = seq_along(genes), gene = genes, score = scores[p, genes], row.names = NULL)
  }))
  result <- list(key = key, programs = programs, original_program_ids = original_program_ids,
    usage = usage, scores = scores, metadata = meta, umap = umap, top = top,
    sample_info = sample_info, sample_usage = sample_usage, subtype_usage = subtype_usage,
    subtypes = subtypes, comparison = comparison, universe = universe,
    enrichment = cnmf_ora(scores, universe, reference, cfg),
    correlation = cor(sample_usage, method = "spearman"),
    reference = reference, source_md5 = tools::md5sum(files), config = cfg, session = sessionInfo())
  qs2::qs_save(result, file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_report.qs2")))
  message(key, ": ", nrow(usage), " cells, K=", item$k, "; universe=", length(universe),
    "; significant enrichment tests=", sum(result$enrichment$fdr < 0.05))
  invisible(result)
}

cnmf_save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in",
    dpi = 600, bg = "white", device = ragg::agg_png, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
}

cnmf_theme <- function() {
  ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", color = NA),
      strip.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(10, 12, 10, 10))
}

cnmf_heat <- function(matrix, row_labels, column_labels = colnames(matrix), legend, divergent = FALSE,
  limits = NULL, numbers = FALSE) {
  d <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
  names(d) <- c("row", "column", "value")
  d$row <- factor(d$row, levels = rev(rownames(matrix)))
  d$column <- factor(d$column, levels = colnames(matrix))
  p <- ggplot2::ggplot(d, ggplot2::aes(column, row, fill = value)) + ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(labels = setNames(column_labels, colnames(matrix)), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = row_labels, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, fill = legend) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(display = "rectangles")) + cnmf_theme() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.4))
  if (divergent) {
    p <- p + ggplot2::scale_fill_gradient2(low = "#3C5488", mid = "#FFFFFF", high = "#E64B35",
      limits = limits, oob = scales::squish)
  } else {
    p <- p + ggplot2::scale_fill_gradientn(colors = c("#F5F5F5", "#91D1C2", "#4DBBD5", "#3C5488"),
      limits = limits, oob = scales::squish)
  }
  if (numbers) p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 2.7)
  p
}

plot_cnmf_lineage <- function(result, cfg) {
  r <- result; key <- r$key; out <- file.path(cfg$output_dir, key)
  a <- read.csv(cfg$annotations, stringsAsFactors = FALSE, encoding = "UTF-8")
  a <- a[a$lineage == key, ]; a <- a[match(r$programs, a$program), ]
  stopifnot(identical(a$program, r$programs), !anyNA(a$label))
  labels <- setNames(paste(a$program, a$label, sep = " | "), a$program)
  short_labels <- setNames(paste(a$program, a$label, sep = "\n"), a$program)
  save <- function(p, stem, width, height) cnmf_save_plot(p, file.path(out, stem), width, height)
  # 保持既有K选择；不重新选择模型，也不使用双轴造成数值误读。
  k <- read.csv(file.path(out, "cNMF_K_selection.csv"))
  kp <- function(y, label, color) ggplot2::ggplot(k, ggplot2::aes(k, .data[[y]])) +
    ggplot2::geom_line(color = color) + ggplot2::geom_point(color = color, size = 2) +
    ggplot2::geom_vline(xintercept = cfg$lineages[[key]]$k, linetype = 2, color = "#777777") +
    ggplot2::scale_x_continuous(breaks = k$k) + ggplot2::labs(x = "Number of programs (K)", y = label) + cnmf_theme()
  p <- patchwork::wrap_plots(kp("silhouette", "Stability", "#3C5488"),
    kp("prediction_error", "Reconstruction error", "#E64B35"), ncol = 2)
  save(p, "Fig1_cNMF_K_selection", 8.8, 3.5)
  genes <- unique(r$top$gene[r$top$rank <= 5])
  heat <- scale(r$scores[, genes, drop = FALSE])
  p <- cnmf_heat(heat, labels, legend = "Gene-score\nZ-score", divergent = TRUE, limits = c(-2.5, 2.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
  save(p, "Fig2_cNMF_program_top_genes", max(10, length(genes) * 0.21 + 3), 4.8)
  heat <- t(scale(r$subtype_usage)); heat[!is.finite(heat)] <- 0
  p <- cnmf_heat(heat, labels, r$subtypes$plot_label, "Mean usage\nrow Z-score", TRUE, c(-2.5, 2.5))
  save(p, "Fig3_cNMF_program_usage_by_subtype", 9.2, 5.2)
  sample_labels <- paste0(r$sample_info$sample, "\n", r$sample_info$group)
  p <- cnmf_heat(t(r$sample_usage), labels, sample_labels, "Mean normalized\nusage", limits = c(0, 1))
  save(p, "Fig4_cNMF_program_usage_by_sample", 8.2, 5.2)
  d <- r$comparison; d$program <- factor(d$program, levels = rev(r$programs))
  d$direction <- ifelse(d$difference >= 0, cfg$group_levels[2], cfg$group_levels[1])
  comparison_name <- paste(cfg$group_levels[2], "-", cfg$group_levels[1])
  p <- ggplot2::ggplot(d, ggplot2::aes(difference, program, fill = direction)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::scale_y_discrete(labels = labels) + ggplot2::scale_fill_manual(values = cfg$group_colors) +
    ggplot2::labs(x = paste0("Mean normalized usage difference (", comparison_name, ")"), y = NULL) +
    cnmf_theme() + ggplot2::theme(legend.position = "none")
  save(p, "Fig5_cNMF_program_group_difference", 8.2, 4.8)
  # 固定已有亚群UMAP坐标，全细胞绘制；与Seurat自动点径一致，不栅格化或下采样。
  um <- as.data.frame(r$umap); names(um) <- c("UMAP_1", "UMAP_2")
  point_size <- Seurat::AutoPointSize(um, raster = FALSE)
  plots <- lapply(r$programs, function(program) {
    z <- um; z$usage <- r$usage[, program]; z <- z[order(z$usage), ]
    ggplot2::ggplot(z, ggplot2::aes(UMAP_1, UMAP_2, color = usage)) +
      ggplot2::geom_point(size = point_size) + ggplot2::coord_equal() +
      ggplot2::scale_color_gradientn(colors = c("#E9E9E9", "#F39B7F", "#E64B35", "#7E1D26"),
        values = c(0, 0.15, 0.45, 1), limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Normalized\nusage") +
      ggplot2::labs(title = short_labels[[program]], x = "UMAP 1", y = "UMAP 2") +
      ggplot2::guides(color = ggplot2::guide_colorbar(display = "rectangles")) + cnmf_theme() +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5))
  })
  p <- patchwork::wrap_plots(plots, ncol = 3, guides = "collect")
  save(p, "Fig6_cNMF_program_usage_UMAP", 11, ceiling(length(plots) / 3) * 3.2)
  # 每点是一只小鼠；短横线表示3只小鼠均值，不连接无配对的样本。
  d <- as.data.frame(as.table(r$sample_usage), stringsAsFactors = FALSE)
  names(d) <- c("sample", "program", "usage")
  d$group <- factor(r$sample_info$group[match(d$sample, r$sample_info$sample)], levels = cfg$group_levels)
  d$program <- factor(d$program, levels = r$programs)
  stat <- r$comparison; stat$program <- factor(stat$program, levels = r$programs)
  stat$label <- sprintf("P = %.2f; FDR = %.2f", stat$p_value, stat$fdr)
  p <- ggplot2::ggplot(d, ggplot2::aes(group, usage, color = group)) +
    ggplot2::geom_point(size = 2.4, position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1001)) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 9, color = "#333333") +
    ggplot2::geom_text(data = stat, ggplot2::aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
      vjust = 1.3, size = 2.8, color = "#333333") +
    ggplot2::facet_wrap(~ program, ncol = 3, scales = "free_y", labeller = ggplot2::as_labeller(short_labels)) +
    ggplot2::scale_color_manual(values = cfg$group_colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.35))) +
    ggplot2::labs(x = NULL, y = "Mean normalized usage per mouse") + cnmf_theme() +
    ggplot2::theme(legend.position = "none")
  save(p, "Fig7_cNMF_program_usage_by_group", 10, ceiling(length(r$programs) / 3) * 3)
  # 每GEP/每集合取FDR最小的2个显著条目；完整显著结果另存一张生物学结果表。
  e <- r$enrichment[r$enrichment$fdr < 0.05, ]
  chosen <- unique(unlist(lapply(split(e, paste(e$program, e$collection)), function(d)
    head(d$pathway[order(d$fdr, -d$overlap, d$pathway)], 2))))
  ed <- e[e$pathway %in% chosen, ]
  row_order <- unique(ed$pathway[order(match(ed$collection, c("Hallmark", "GO_BP")),
    match(ed$program, r$programs), ed$fdr)])
  ed$pathway <- factor(ed$pathway, levels = rev(row_order)); ed$program <- factor(ed$program, levels = r$programs)
  ed$collection <- factor(ed$collection, levels = c("Hallmark", "GO_BP"))
  ed$strength <- pmin(-log10(ed$fdr), 10)
  path_labels <- vapply(row_order, function(s) paste(strwrap(tolower(gsub("_", " ",
    sub("^(HALLMARK|GOBP)_", "", s))), width = 54), collapse = "\n"), character(1))
  p <- ggplot2::ggplot(ed, ggplot2::aes(program, pathway, size = overlap, color = strength)) +
    ggplot2::geom_point() + ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_y_discrete(labels = setNames(path_labels, row_order)) +
    ggplot2::scale_color_gradientn(colors = c("#F7C4B6", "#E64B35", "#7E1D26"),
      limits = c(0, 10), breaks = c(0, 5, 10), labels = c("0", "5", ">=10")) +
    ggplot2::scale_size_continuous(range = c(1.5, 5.5)) +
    ggplot2::labs(x = NULL, y = NULL, color = "-log10(FDR)", size = "Overlap genes") +
    ggplot2::guides(color = ggplot2::guide_colorbar(display = "rectangles")) + cnmf_theme() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#ECECEC", linewidth = 0.3),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.4), strip.text.y = ggplot2::element_text(angle = 0))
  save(p, "Fig8_cNMF_Hallmark_GO_enrichment", 10, max(6.5, length(row_order) * 0.29 + 1.2))
  p <- cnmf_heat(r$correlation, labels, short_labels, "Spearman rho\n6 mice", TRUE, c(-1, 1), TRUE)
  save(p, "Fig9_cNMF_GEP_correlation", 9.5, 6.8)
  a <- merge(a, r$comparison, by = "program", sort = FALSE); a <- a[match(r$programs, a$program), ]
  a$top_Hallmark <- vapply(a$program, function(p) {
    d <- e[e$program == p & e$collection == "Hallmark", ]; paste(head(d$pathway[order(d$fdr)], 3), collapse = ";")
  }, character(1))
  a$top_GO_BP <- vapply(a$program, function(p) {
    d <- e[e$program == p & e$collection == "GO_BP", ]; paste(head(d$pathway[order(d$fdr)], 3), collapse = ";")
  }, character(1))
  a$contrast <- comparison_name
  data.table::fwrite(a, file.path(out, "cNMF_program_annotation.csv"), bom = TRUE)
  data.table::fwrite(e[order(e$program, e$collection, e$fdr), ],
    file.path(out, "cNMF_pathway_enrichment.csv"), bom = TRUE)
  message(key, ": 9 PNG/PDF pairs and functional annotations written")
  invisible(r)
}

run_cnmf_report <- function(cfg) {
  reference <- cnmf_reference(cfg)
  for (key in names(cfg$lineages)) {
    path <- file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_report.qs2"))
    if (file.exists(path)) {
      r <- qs2::qs_read(path)
      stopifnot(identical(r$source_md5, tools::md5sum(names(r$source_md5))),
        identical(r$reference$md5, reference$md5), identical(r$config, cfg))
    } else r <- prepare_cnmf_lineage(cfg, key, reference)
    plot_cnmf_lineage(r, cfg)
  }
}
