# 按原始RNA UMI汇总小鼠×亚群，edgeR稳健QL检验；富集背景仅为该比较实际检验的基因。
subtype_ora <- function(de, reference, cfg) {
  universe <- de$gene
  out <- list()
  for (collection in c("Hallmark", "GO_BP")) {
    sets <- lapply(reference[[collection]], intersect, y = universe)
    sets <- sets[lengths(sets) >= cfg$set_size[1] & lengths(sets) <= cfg$set_size[2]]
    for (direction in c("Up", "Down")) {
      selected <- de$gene[de$direction == direction]
      overlap <- lapply(sets, intersect, y = selected)
      hits <- lengths(overlap); sizes <- lengths(sets)
      p <- phyper(hits - 1L, sizes, length(universe) - sizes, length(selected), lower.tail = FALSE)
      out[[paste(collection, direction)]] <- data.frame(collection, direction, pathway = names(sets),
        overlap = unname(hits), selected_n = length(selected), set_size = unname(sizes), universe_n = length(universe),
        p_value = unname(p), fdr = unname(p.adjust(p, "BH")),
        genes = vapply(overlap, paste, character(1), collapse = ";"), row.names = NULL)
    }
  }
  do.call(rbind, out)
}

# 排序两端非常不平衡时，增加初始置换估计精度；不调整FDR阈值或改变检验方向。
subtype_fgsea <- function(sets, ranks, cfg, seed) {
  for (permutations in c(1000L, 10000L, 100000L)) {
    set.seed(seed)
    g <- as.data.frame(fgsea::fgseaMultilevel(sets, ranks,
      minSize = cfg$set_size[1], maxSize = cfg$set_size[2], eps = 0,
      nPermSimple = permutations, sampleSize = 101L, nproc = 1L))
    missing <- sum(is.na(g$pval))
    if (!missing) {
      if (permutations > 1000L) message("GSEA_PRECISION_OK nPermSimple=", permutations)
      break
    }
    message("GSEA_PRECISION nPermSimple=", permutations, " missing_p=", missing)
  }
  if (missing) warning("GSEA retains ", missing, " unestimable pathways at nPermSimple=100000; do not interpret them as nonsignificant")
  g
}

subtype_gsea <- function(de, reference, cfg, seed) {
  # 一自由度QL F值开根号后按logFC定向，用全部通过表达过滤的基因排序；不只选DEG。
  # QL F在接近零时可有浮点负误差；保留基因并将该排序统计量归零，不静默丢弃。
  stopifnot(all(is.finite(de$F)), min(de$F) >= -1e-6)
  stat <- setNames(sign(de$logFC) * sqrt(pmax(de$F, 0)), de$gene)
  stat <- sort(stat, decreasing = TRUE)
  out <- list()
  for (collection in c("Hallmark", "GO_BP")) {
    g <- subtype_fgsea(reference[[collection]], stat, cfg, seed)
    g$leadingEdge <- vapply(g$leadingEdge, paste, character(1), collapse = ";")
    g$collection <- collection
    out[[collection]] <- g
  }
  do.call(rbind, out)
}

plot_subtype_de <- function(result, cfg, directory, label) {
  d <- result$de; d$direction <- factor(d$direction, levels = c("Down", "Not significant", "Up"))
  colors <- c(Down = cfg$group_colors[[cfg$group_levels[1]]], `Not significant` = "#C4C4C4",
    Up = cfg$group_colors[[cfg$group_levels[2]]])
  chosen <- head(d$gene[d$direction != "Not significant"][order(d$FDR[d$direction != "Not significant"])], 8)
  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, -log10(pmax(FDR, .Machine$double.xmin)), color = direction)) +
    ggplot2::geom_point(size = 0.65, stroke = 0, shape = 16) +
    ggplot2::geom_vline(xintercept = c(-cfg$logfc, cfg$logfc), linetype = 2, color = "#777777", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(cfg$fdr), linetype = 2, color = "#777777", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = d[d$gene %in% chosen, ], ggplot2::aes(label = gene),
      size = 2.8, seed = cfg$seed, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::labs(title = label, x = paste0("log2 fold change (", cfg$group_levels[2], " - ", cfg$group_levels[1], ")"),
      y = "-log10(BH FDR)", color = NULL) + cnmf_theme() + ggplot2::theme(legend.position = "bottom")
  # 同图显示最显著30个检验基因的样本表达；无DEG时不把Top-ranked命名为显著。
  genes <- head(d$gene[order(d$FDR, d$PValue)], 30)
  heat <- t(scale(t(result$logcpm[genes, , drop = FALSE]))); heat[!is.finite(heat)] <- 0
  info <- result$samples
  labels <- paste0(info$sample, "\n", info$group, "; n=", info$n_cells)
  p2 <- cnmf_heat(heat, setNames(genes, genes), labels, "logCPM\nrow Z-score", TRUE, c(-2, 2)) +
    ggplot2::labs(title = "Top-ranked genes") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), axis.text.x = ggplot2::element_text(size = 8, angle = 55, hjust = 1))
  cnmf_save_plot(patchwork::wrap_plots(p, p2, ncol = 2, widths = c(1.1, 1)),
    file.path(directory, "Fig1_pseudobulk_DE"), 12, 7.2)
  g <- result$gsea
  g <- g[!is.na(g$padj) & g$padj < cfg$fdr, ]
  if (nrow(g)) {
    chosen <- do.call(rbind, lapply(split(g, paste(g$collection, sign(g$NES))), function(x)
      head(x[order(x$padj, -abs(x$NES), x$pathway), ], 5)))
    chosen$pathway <- factor(chosen$pathway, levels = unique(chosen$pathway[order(chosen$collection, chosen$NES)]))
    path_labels <- vapply(levels(chosen$pathway), function(s) paste(strwrap(tolower(gsub("_", " ",
      sub("^(HALLMARK|GOBP)_", "", s))), width = 52), collapse = "\n"), character(1))
    chosen$direction <- ifelse(chosen$NES > 0, cfg$group_levels[2], cfg$group_levels[1])
    chosen$strength <- pmin(-log10(chosen$padj), 10)
    p <- ggplot2::ggplot(chosen, ggplot2::aes(NES, pathway, color = direction, size = strength)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#AAAAAA") + ggplot2::geom_point() +
      ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_y_discrete(labels = setNames(path_labels, levels(chosen$pathway))) +
      ggplot2::scale_color_manual(values = cfg$group_colors, name = "Enriched in") +
      ggplot2::scale_size_continuous(range = c(2, 5), limits = c(0, 10), breaks = c(2, 5, 10), labels = c("2", "5", ">=10")) +
      ggplot2::labs(title = label, x = "Normalized enrichment score (EXP positive)", y = NULL, size = "-log10(FDR)") +
      cnmf_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#ECECEC", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.4), axis.line = ggplot2::element_blank())
    height <- max(5, nrow(chosen) * 0.33 + 1.4)
  } else {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No Hallmark / GO BP pathway at BH FDR < 0.05", size = 4) +
      ggplot2::labs(title = label) + ggplot2::theme_void(base_family = "Arial")
    height <- 4
  }
  cnmf_save_plot(p, file.path(directory, "Fig2_pseudobulk_GSEA"), 10, height)
}

run_subtype_pseudobulk <- function(cfg, keys = names(cfg$lineages)) {
  stopifnot(all(keys %in% names(cfg$lineages)))
  dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
  reference <- cnmf_reference(cfg)
  for (key in keys) {
    item <- cfg$lineages[[key]]
    input <- qs2::qs_read(item$input); meta <- input$metadata
    stopifnot(identical(colnames(input$counts), meta$cell), !anyDuplicated(meta$cell), !anyNA(meta),
      !anyDuplicated(rownames(input$counts)), all(input$counts@x >= 0),
      all(input$counts@x == floor(input$counts@x)))
    annotation <- read.csv(item$annotation, stringsAsFactors = FALSE)
    stopifnot(setequal(annotation$annotation, meta$subtype))
    # 明确核对当前正式对象中的细胞、分组与亚群，防止缓存注释过期。
    obj <- qs2::qs_read(item$object); om <- obj[[]][meta$cell, ]
    stopifnot(setequal(colnames(obj), meta$cell), identical(as.character(om$samples), meta$sample),
      identical(as.character(om$group), meta$group), identical(as.character(om[[item$subtype_col]]), meta$subtype))
    rm(obj, om); gc(FALSE)
    samples <- unique(meta[, c("sample", "group")])
    samples <- samples[order(match(samples$group, cfg$group_levels), samples$sample), ]
    stopifnot(!anyDuplicated(samples$sample), all(samples$group %in% cfg$group_levels))
    design_cells <- Matrix::sparseMatrix(i = seq_len(nrow(meta)), j = match(meta$sample, samples$sample),
      x = 1, dims = c(nrow(meta), nrow(samples)))
    cache_path <- file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_EXP_vs_CON.qs2"))
    hashes <- tools::md5sum(c(item$input, item$object, item$annotation))
    if (file.exists(cache_path)) {
      cache <- qs2::qs_read(cache_path)
      stopifnot(identical(cache$source_md5, hashes), identical(cache$config, cfg), identical(cache$reference$md5, reference$md5))
    } else cache <- list(results = list(), status = list(), config = cfg, source_md5 = hashes,
      reference = reference, session = sessionInfo())
    dir.create(item$output, recursive = TRUE, showWarnings = FALSE)
    for (j in seq_len(nrow(annotation))) {
      subtype <- annotation$annotation[j]
      if (subtype %in% names(cache$status)) next
      cells <- meta$subtype == subtype
      aggregated <- as.matrix(input$counts[, cells, drop = FALSE] %*% design_cells[cells, , drop = FALSE])
      colnames(aggregated) <- samples$sample
      nc <- as.integer(table(factor(meta$sample[cells], levels = samples$sample)))
      info <- samples; info$n_cells <- nc; info$used <- nc >= cfg$min_cells
      n <- table(factor(info$group[info$used], levels = cfg$group_levels))
      status <- data.frame(lineage = key, cluster = annotation$cluster[j], subtype,
        n_reference = unname(n[1]), n_test = unname(n[2]), status = "tested",
        tested_genes = 0L, up = 0L, down = 0L, ora_significant = 0L, gsea_significant = 0L,
        stringsAsFactors = FALSE)
      status <- cbind(status, as.data.frame(as.list(setNames(nc, paste0("cells_", samples$sample)))))
      message(key, " / ", subtype, ": n=", n[1], " vs ", n[2])
      if (any(n < cfg$min_mice)) {
        status$status <- "insufficient_mice_after_cell_filter"
        cache$results[[subtype]] <- list(counts = aggregated, samples = info)
      } else {
        used <- info[info$used, ]
        group <- factor(used$group, levels = cfg$group_levels)
        design <- model.matrix(~ group)
        y <- edgeR::DGEList(aggregated[, info$used, drop = FALSE], group = group)
        keep <- edgeR::filterByExpr(y, design = design)
        stopifnot(sum(keep) >= cfg$min_genes, qr(design)$rank == ncol(design))
        y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
        y <- edgeR::estimateDisp(y, design, robust = TRUE)
        fit <- edgeR::glmQLFit(y, design, robust = TRUE)
        test <- edgeR::glmQLFTest(fit, coef = 2)
        de <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
        de$gene <- rownames(de)
        de$direction <- ifelse(de$FDR < cfg$fdr & de$logFC >= cfg$logfc, "Up",
          ifelse(de$FDR < cfg$fdr & de$logFC <= -cfg$logfc, "Down", "Not significant"))
        de <- de[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "direction")]
        ora <- subtype_ora(de, reference, cfg)
        gsea <- subtype_gsea(de, reference, cfg, cfg$seed + match(key, names(cfg$lineages)) * 100L + j)
        result <- list(counts = aggregated, samples = used, all_samples = info, de = de, ora = ora, gsea = gsea,
          logcpm = edgeR::cpm(y, log = TRUE, prior.count = 2), keep = keep,
          library = y$samples, dispersion = y$tagwise.dispersion, design = design)
        status$tested_genes <- nrow(de); status$up <- sum(de$direction == "Up"); status$down <- sum(de$direction == "Down")
        status$ora_significant <- sum(ora$fdr < cfg$fdr)
        status$gsea_significant <- sum(gsea$padj < cfg$fdr, na.rm = TRUE)
        slug <- gsub("[^A-Za-z0-9]+", "_", subtype)
        directory <- file.path(item$output, sprintf("C%02d_%s", annotation$cluster[j], slug))
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(de[order(de$FDR, de$PValue), ], file.path(directory, "DEG_EXP_vs_CON.csv"), bom = TRUE)
        data.table::fwrite(ora[ora$fdr < cfg$fdr, ], file.path(directory, "ORA_Hallmark_GO_BP.csv"), bom = TRUE)
        data.table::fwrite(gsea[order(gsea$collection, gsea$padj), ], file.path(directory, "GSEA_Hallmark_GO_BP.csv"), bom = TRUE)
        plot_subtype_de(result, cfg, directory, subtype)
        cache$results[[subtype]] <- result
        message("  genes=", nrow(de), "; Up=", status$up, "; Down=", status$down, "; GSEA=", status$gsea_significant)
      }
      cache$status[[subtype]] <- status
      qs2::qs_save(cache, cache_path)
      data.table::fwrite(do.call(rbind, cache$status), file.path(item$output, "pseudobulk_summary.csv"), bom = TRUE)
    }
    message(key, " PSEUDOBULK_COMPLETE: ", length(cache$status), " subtypes")
  }
  invisible(TRUE)
}
