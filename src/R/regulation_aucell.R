# 对预先选择的原生物种MSigDB调控集进行AUCell；基因集/路径/标签均由配置指定。
prepare_regulation_reference <- function(cfg, universe) {
  database <- readRDS(cfg$reference_cache)
  members <- database[database$gs_name %in% cfg$sets$id, ]
  stopifnot(setequal(members$gs_name, cfg$sets$id),
    all(members$db_version == cfg$reference_version), all(members$gs_subcollection == "GO:BP"))
  online <- setNames(lapply(cfg$sets$id, function(id) {
    url <- paste0(cfg$download_base, "?geneSetName=", id, "&fileType=json")
    x <- jsonlite::fromJSON(url)[[id]]
    stopifnot(setequal(x$geneSymbols, members$db_gene_symbol[members$gs_name == id]),
      identical(x$exactSource, cfg$sets$go_id[match(id, cfg$sets$id)]))
    list(url = url, content = x)
  }), cfg$sets$id)
  mapping <- unique(members[c("db_gene_symbol", "db_ncbi_gene")])
  names(mapping) <- c("source_gene", "entrez")
  mapping$input_gene <- ifelse(mapping$source_gene %in% universe, mapping$source_gene, NA_character_)
  mapping$method <- ifelse(is.na(mapping$input_gene), "absent", "exact_symbol")
  db <- getExportedValue(cfg$orgdb, cfg$orgdb)
  for (i in which(is.na(mapping$input_gene))) {
    a <- AnnotationDbi::select(db, keys = as.character(mapping$entrez[i]),
      keytype = "ENTREZID", columns = "ALIAS")
    hit <- intersect(a$ALIAS, universe)
    if (length(hit) == 1L) {
      reverse <- AnnotationDbi::select(db, keys = hit, keytype = "ALIAS", columns = "ENTREZID")
      if (setequal(reverse$ENTREZID, as.character(mapping$entrez[i]))) {
        mapping$input_gene[i] <- hit
        mapping$method[i] <- "unique_Entrez_alias"
      }
    }
  }
  members$input_gene <- mapping$input_gene[match(members$db_gene_symbol, mapping$source_gene)]
  sets <- setNames(lapply(cfg$sets$id, function(id)
    unique(na.omit(members$input_gene[members$gs_name == id]))), cfg$sets$id)
  stopifnot(all(lengths(sets) >= cfg$minimum_genes))
  coverage <- data.frame(id = names(sets), source_n = vapply(names(sets), function(id)
    length(unique(members$db_gene_symbol[members$gs_name == id])), integer(1)),
    used_n = lengths(sets), row.names = NULL)
  result <- list(members = members, mapping = mapping, sets = sets, coverage = coverage,
    online = online, retrieved_at = as.character(Sys.time()),
    source_md5 = tools::md5sum(cfg$reference_cache), orgdb_metadata = AnnotationDbi::metadata(db))
  dir.create(dirname(cfg$reference_output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, cfg$reference_output)
  print(coverage)
  print(mapping[mapping$method != "exact_symbol", ])
  result
}

summarize_regulation_aucell <- function(scores, meta, annotation, cfg) {
  scopes <- c(cfg$overall_label, annotation$annotation)
  sample_info <- unique(meta[c("sample", "group")])
  sample_info <- sample_info[order(match(sample_info$group, cfg$group_levels), sample_info$sample), ]
  stopifnot(!anyDuplicated(sample_info$sample), setequal(sample_info$group, cfg$group_levels))
  values <- list()
  for (scope in scopes) for (s in sample_info$sample) {
    cells <- meta$cell[meta$sample == s & (scope == cfg$overall_label | meta$subtype == scope)]
    values[[paste(scope, s)]] <- data.frame(scope, sample = s,
      group = sample_info$group[match(s, sample_info$sample)], id = rownames(scores),
      n_cells = length(cells), auc = if (length(cells)) rowMeans(scores[, cells, drop = FALSE]) else NA_real_)
  }
  sample_scores <- do.call(rbind, values); rownames(sample_scores) <- NULL
  sample_scores$used <- sample_scores$n_cells >= cfg$minimum_cells
  comparisons <- list()
  for (scope in scopes) for (id in rownames(scores)) {
    d <- sample_scores[sample_scores$scope == scope & sample_scores$id == id, ]
    all_means <- vapply(cfg$group_levels, function(g) mean(d$auc[d$group == g], na.rm = TRUE), numeric(1))
    a <- d$auc[d$group == cfg$group_levels[1] & d$used]
    b <- d$auc[d$group == cfg$group_levels[2] & d$used]
    enough <- length(a) >= cfg$minimum_mice & length(b) >= cfg$minimum_mice
    delta <- g <- p <- sd_pool <- NA_real_
    if (enough) {
      delta <- mean(b) - mean(a)
      df <- length(a) + length(b) - 2L
      sd_pool <- sqrt((sum((a - mean(a))^2) + sum((b - mean(b))^2)) / df)
      correction <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
      if (sd_pool > 0) g <- correction * delta / sd_pool
      p <- cnmf_exact_p(a, b)
    }
    comparisons[[paste(scope, id)]] <- data.frame(scope, id,
      label = cfg$sets$label[match(id, cfg$sets$id)],
      mean_reference_all_mice = unname(all_means[1]), mean_test_all_mice = unname(all_means[2]),
      descriptive_delta = unname(all_means[2] - all_means[1]),
      n_reference = length(a), n_test = length(b),
      mean_reference = if (length(a)) mean(a) else NA_real_,
      mean_test = if (length(b)) mean(b) else NA_real_, delta, pooled_sd = sd_pool, hedges_g = g,
      p_value = p, status = if (enough) "tested" else "insufficient_mice")
  }
  comparison <- do.call(rbind, comparisons); rownames(comparison) <- NULL
  comparison$fdr <- NA_real_
  tested <- comparison$status == "tested"
  comparison$fdr[tested] <- p.adjust(comparison$p_value[tested], "BH")
  list(sample_scores = sample_scores, comparison = comparison)
}

plot_regulation_aucell <- function(r, cfg) {
  save <- function(p, name, w, h) cnmf_save_plot(p, file.path(cfg$output, paste0(name, "_", cfg$label_prefix)), w, h)
  plots <- list()
  point_size <- Seurat::AutoPointSize(as.data.frame(r$umap))
  for (id in cfg$sets$id) for (group in cfg$group_levels) {
    cells <- r$metadata$cell[r$metadata$group == group]
    d <- data.frame(r$umap[cells, , drop = FALSE], auc = as.numeric(r$scores[id, cells]))
    names(d)[1:2] <- c("UMAP_1", "UMAP_2")
    d <- d[order(d$auc), ]
    plots[[paste(id, group)]] <- ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, color = auc)) +
      ggplot2::geom_point(size = point_size, shape = 16, alpha = 1) +
      ggplot2::coord_equal(xlim = range(r$umap[, 1]), ylim = range(r$umap[, 2])) +
      ggplot2::scale_color_gradientn(colors = cfg$feature_colors, limits = c(0, max(r$scores[id, ])), name = "AUCell") +
      ggplot2::guides(color = ggplot2::guide_colorbar(display = "rectangles")) +
      ggplot2::labs(title = paste(cfg$sets$label[match(id, cfg$sets$id)], group, sep = " | "), x = "UMAP 1", y = "UMAP 2") +
      cnmf_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 11, hjust = 0.5))
  }
  save(patchwork::wrap_plots(plots, ncol = 2), "Fig1_UMAP_EXP_vs_CON", 12, 9)
  # 细胞分布用于可视化；黑点为每只小鼠均值，P值来自小鼠级检验。
  d <- do.call(rbind, lapply(cfg$sets$id, function(id)
    data.frame(r$metadata, id, auc = as.numeric(r$scores[id, r$metadata$cell]))))
  d$group <- factor(d$group, cfg$group_levels)
  label_by_scope <- setNames(r$annotation$plot_label, r$annotation$annotation)
  label_by_id <- setNames(cfg$sets$label, cfg$sets$id)
  d$subtype <- factor(d$subtype, r$annotation$annotation)
  d$id <- factor(d$id, cfg$sets$id)
  s <- r$sample_scores[r$sample_scores$scope != cfg$overall_label, ]
  s$subtype <- factor(s$scope, r$annotation$annotation); s$group <- factor(s$group, cfg$group_levels)
  s$id <- factor(s$id, cfg$sets$id)
  a <- r$comparison[r$comparison$scope != cfg$overall_label, ]
  a$subtype <- factor(a$scope, r$annotation$annotation); a$id <- factor(a$id, cfg$sets$id)
  a$text <- ifelse(is.na(a$p_value), "Not tested", sprintf("P = %.2f", a$p_value))
  p <- ggplot2::ggplot(d, ggplot2::aes(group, auc, fill = group)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.85, linewidth = 0.3) +
    ggplot2::geom_point(data = s, ggplot2::aes(shape = used), color = "black", fill = "white", size = 1.6,
      position = ggplot2::position_jitter(width = 0.065, height = 0, seed = cfg$step)) +
    ggplot2::geom_text(data = a, ggplot2::aes(x = 1.5, y = Inf, label = text), vjust = 1.3,
      size = 2.5, inherit.aes = FALSE) +
    ggplot2::facet_grid(id ~ subtype, scales = "free_y",
      labeller = ggplot2::labeller(id = label_by_id, subtype = label_by_scope)) +
    ggplot2::scale_fill_manual(values = cfg$group_colors, name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "Mouse mean: >=20 cells", `FALSE` = "Mouse mean: <20 cells"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(x = NULL, y = "Single-cell AUCell") + cnmf_theme() +
    ggplot2::theme(legend.position = "bottom", strip.text.x = ggplot2::element_text(size = 8),
      strip.text.y = ggplot2::element_text(size = 9), axis.text.x = ggplot2::element_text(size = 8))
  save(p, "Fig2_Violin_EXP_vs_CON", 16, 7)
  scopes <- c(cfg$overall_label, r$annotation$annotation)
  column_labels <- c(cfg$overall_label, r$annotation$plot_label)
  tested_scope <- vapply(scopes, function(s) any(r$comparison$status[r$comparison$scope == s] == "tested"), logical(1))
  column_labels[!tested_scope] <- paste0(column_labels[!tested_scope], " (not tested)")
  heat <- function(field, title, legend, digits) {
    mat <- matrix(NA_real_, nrow(cfg$sets), length(scopes), dimnames = list(cfg$sets$id, scopes))
    for (i in seq_len(nrow(r$comparison))) {
      row <- r$comparison[i, ]; mat[row$id, row$scope] <- row[[field]]
    }
    lim <- max(abs(mat), na.rm = TRUE); if (lim == 0) lim <- 0.001
    p <- cnmf_heat(mat, label_by_id, column_labels, paste0(legend, "\nGray: not tested"), TRUE, c(-lim, lim))
    p$data$text <- ifelse(is.finite(p$data$value), sprintf(paste0("%.", digits, "f"), p$data$value), "")
    p + ggplot2::geom_text(ggplot2::aes(label = text), size = 2.7) + ggplot2::labs(title = title)
  }
  h1 <- heat("delta", "Mean AUCell difference (EXP - CON)", "Delta AUCell", 4)
  h2 <- heat("hedges_g", "Standardized effect (EXP - CON)", "Hedges' g", 2)
  save(patchwork::wrap_plots(h1, h2, ncol = 1), "Fig3_Heatmap_EXP_vs_CON", 13.5, 8)
  invisible(NULL)
}

run_regulation_aucell <- function(cfg) {
  stopifnot(requireNamespace("Seurat", quietly = TRUE), requireNamespace("Matrix", quietly = TRUE))
  protected <- tools::md5sum(c(cfg$input, cfg$control_cache))
  old <- qs2::qs_read(cfg$control_cache)
  reference <- prepare_regulation_reference(cfg, old$universe)
  message("Reading original annotated RNA object; no write-back")
  obj <- qs2::qs_read(cfg$input)
  md <- obj[[]]
  meta <- data.frame(cell = rownames(md), sample = as.character(md[[cfg$sample_col]]),
    group = as.character(md[[cfg$group_col]]), subtype = as.character(md[[cfg$subtype_col]]))
  stopifnot(identical(meta, old$metadata[names(meta)]))
  umap <- SeuratObject::Embeddings(obj, cfg$reduction)[meta$cell, , drop = FALSE]
  annotation <- read.csv(cfg$annotation, stringsAsFactors = FALSE)
  layers <- SeuratObject::Layers(obj[[cfg$assay]], search = "^counts")
  counts <- lapply(layers, function(x) SeuratObject::LayerData(obj, assay = cfg$assay, layer = x))
  rm(obj, md); gc(FALSE)
  stopifnot(all(vapply(counts, function(x) identical(rownames(x), old$universe), logical(1))),
    setequal(unlist(lapply(counts, colnames)), meta$cell), !anyDuplicated(unlist(lapply(counts, colnames))))
  sets <- c(reference$sets, old$gene_sets[cfg$control_set])
  cutoff <- ceiling(length(old$universe) * cfg$auc_fraction)
  stopifnot(cutoff == old$auc_max_rank)
  scores <- matrix(NA_real_, length(sets), nrow(meta), dimnames = list(names(sets), meta$cell))
  set.seed(cfg$seed)
  for (i in seq_along(counts)) {
    expr <- counts[[i]]
    blocks <- split(seq_len(ncol(expr)), ceiling(seq_len(ncol(expr)) / cfg$chunk_cells))
    for (j in seq_along(blocks)) {
      message(format(Sys.time(), "%H:%M:%S"), " AUCell ", layers[i], " block ", j, "/", length(blocks))
      b <- expr[, blocks[[j]], drop = FALSE]
      scores[, colnames(b)] <- program_aucell_block(b, sets, cutoff)[names(sets), colnames(b), drop = FALSE]
      rm(b); gc(FALSE)
    }
    counts[i] <- list(NULL); rm(expr); gc(FALSE)
  }
  control_error <- max(abs(scores[cfg$control_set, ] - old$scores[cfg$control_set, meta$cell]))
  stopifnot(control_error < 1e-10, all(is.finite(scores)), all(scores >= 0 & scores <= 1))
  selected <- scores[cfg$sets$id, , drop = FALSE]
  summary <- summarize_regulation_aucell(selected, meta, annotation, cfg)
  result <- c(list(scores = selected, metadata = meta, umap = umap, annotation = annotation,
    reference = reference, universe = old$universe, auc_max_rank = cutoff,
    control_max_error = control_error, control_scores = scores[cfg$control_set, ],
    source_md5 = protected, config = cfg, session = sessionInfo()), summary)
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(result, cfg$cache)
  data.table::fwrite(result$comparison, file.path(cfg$output, paste0(cfg$label_prefix, "_EXP_vs_CON.csv")), bom = TRUE)
  plot_regulation_aucell(result, cfg)
  stopifnot(identical(protected, tools::md5sum(names(protected))))
  print(result$comparison[, c("scope", "label", "delta", "hedges_g", "p_value", "fdr", "status")])
  message("REGULATION_AUCELL_OK cells=", ncol(selected), " sets=", nrow(selected),
    " cutoff=", cutoff, " Hallmark_control_max_error=", control_error)
  invisible(result)
}
