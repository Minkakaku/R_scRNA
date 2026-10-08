# 下载FerrDb官方三类记录，只提取指定死亡模式；以HGNC及MGI映射到小鼠基因。
prepare_ferrdb_reference <- function(cfg) {
  dir.create(cfg$reference_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(cfg$reference_dir, paste0("FerrDb_v3_", cfg$categories, ".csv"))
  records <- lapply(seq_along(paths), function(i) {
    if (!file.exists(paths[i])) {
      h <- curl::new_handle(postfields = paste0("selectDataset=", cfg$categories[i], "&selectFormat=csv"))
      curl::handle_setheaders(h, "Content-Type" = "application/x-www-form-urlencoded")
      curl::curl_fetch_disk(cfg$download_url, paths[i], handle = h)
    }
    d <- data.table::fread(paths[i], data.table = FALSE, na.strings = c("_NA_", ""))
    stopifnot(all(d$Category == cfg$categories[i]))
    d[d$RCD == cfg$death_mode, ]
  })
  records <- do.call(rbind, records)
  hgnc <- data.table::fread(cfg$hgnc, data.table = FALSE)
  hgnc <- hgnc[hgnc$status == "Approved", ]
  mgi <- data.table::fread(cfg$mgi, data.table = FALSE, quote = "")
  mouse <- split(mgi$Symbol[mgi$`NCBI Taxon ID` == 10090], mgi$`DB Class Key`[mgi$`NCBI Taxon ID` == 10090])
  human <- mgi[mgi$`NCBI Taxon ID` == 9606, ]
  aliases <- lapply(seq_len(nrow(hgnc)), function(i)
    unique(strsplit(paste(hgnc$prev_symbol[i], hgnc$alias_symbol[i], sep = "|"), "|", fixed = TRUE)[[1]]))
  alias_map <- split(rep(hgnc$symbol, lengths(aliases)), unlist(aliases))
  genes <- unique(records[, c("Symbol", "HGNC_ID", "Gene_type")])
  mapping <- lapply(seq_len(nrow(genes)), function(i) {
    g <- genes[i, ]; hs <- NA_character_; candidates <- character()
    # 3'基因计数不能测量成熟miRNA或circRNA；不把它们替换成宿主基因。
    incompatible <- grepl("circ|micro|mirc|mirco|transfer", g$Gene_type, ignore.case = TRUE)
    if (is.na(incompatible)) incompatible <- FALSE
    if (!incompatible) {
      id <- match(g$HGNC_ID, hgnc$hgnc_id)
      if (!is.na(id)) hs <- hgnc$symbol[id]
      else if (!is.na(g$Symbol) && g$Symbol %in% hgnc$symbol) hs <- g$Symbol
      else if (!is.na(g$Symbol)) {
        a <- unique(alias_map[[g$Symbol]])
        if (length(a) == 1L) hs <- a
      }
      if (!is.na(hs)) {
        hgnc_id <- hgnc$hgnc_id[match(hs, hgnc$symbol)]
        classes <- human$`DB Class Key`[human$`HGNC ID` == hgnc_id | human$Symbol == hs]
        candidates <- sort(unique(unlist(mouse[as.character(classes)], use.names = FALSE)))
      }
    }
    roles <- cfg$categories[cfg$categories %in% records$Category[records$Symbol == g$Symbol & !is.na(records$Symbol)]]
    data.frame(source_gene = g$Symbol, human_symbol = hs, gene_type = g$Gene_type,
      mouse_symbol = if (length(candidates)) candidates else NA_character_,
      roles = paste(roles, collapse = ";"), n_mouse_homologs = length(candidates),
      mapping_status = if (incompatible) "not_measured_by_gene_count_assay" else if (is.na(hs)) "unresolved_HGNC"
        else if (!length(candidates)) "no_MGI_mouse_homolog" else if (length(candidates) > 1L) "multiple_MGI_homologs_retained" else "single_MGI_mouse_homolog")
  })
  list(records = records, mapping = do.call(rbind, mapping), version = cfg$version,
    source_page = cfg$source_page, download_url = cfg$download_url, retrieved = Sys.Date(),
    source_md5 = tools::md5sum(c(paths, cfg$hgnc, cfg$mgi)))
}

run_ferrdb_featureplots <- function(cfg, pages = NULL) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  if (file.exists(cfg$cache)) {
    result <- qs2::qs_read(cfg$cache)
    stopifnot(identical(result$config, cfg), identical(result$source_md5, tools::md5sum(names(result$source_md5))),
      identical(result$reference$source_md5, tools::md5sum(names(result$reference$source_md5))))
  } else {
    reference <- prepare_ferrdb_reference(cfg)
    obj <- qs2::qs_read(cfg$object)
    umap <- Seurat::Embeddings(obj, cfg$reduction)
    stopifnot(identical(rownames(umap), colnames(obj)))
    map <- reference$mapping
    all_genes <- sort(unique(na.omit(map$mouse_symbol)))
    genes <- intersect(all_genes, rownames(obj[[cfg$assay]]))
    genes <- c(intersect(cfg$priority_genes, genes), setdiff(genes, cfg$priority_genes))
    # 仅提取既有RNA data层，不使用SCT/integrated，也不重新标准化或改变原对象。
    layers <- SeuratObject::Layers(obj[[cfg$assay]], search = "^data")
    stopifnot(length(layers) > 0)
    expression <- do.call(cbind, lapply(layers, function(layer)
      SeuratObject::LayerData(obj, assay = cfg$assay, layer = layer)[genes, , drop = FALSE]))
    stopifnot(!anyDuplicated(colnames(expression)), setequal(colnames(expression), colnames(obj)))
    expression <- expression[, colnames(obj), drop = FALSE]
    counts <- do.call(cbind, lapply(SeuratObject::Layers(obj[[cfg$assay]], search = "^counts"), function(layer)
      SeuratObject::LayerData(obj, assay = cfg$assay, layer = layer)[genes, , drop = FALSE]))
    counts <- counts[, colnames(obj), drop = FALSE]
    stopifnot(all(expression >= 0), all(is.finite(expression@x)),
      all((expression > 0) == (counts > 0)))
    labels <- vapply(genes, function(g) {
      roles <- unique(unlist(strsplit(map$roles[!is.na(map$mouse_symbol) & map$mouse_symbol == g], ";", fixed = TRUE)))
      paste(cfg$categories[cfg$categories %in% roles], collapse = ";")
    }, character(1))
    catalogue <- data.frame(gene = genes, roles = unname(labels),
      detected_cells = Matrix::rowSums(counts > 0), mean_log_expression = Matrix::rowMeans(expression),
      max_log_expression = apply(as.matrix(expression), 1, max),
      page = ceiling(seq_along(genes) / cfg$per_page))
    reference$mapping$in_expression_matrix <- !is.na(map$mouse_symbol) & map$mouse_symbol %in% genes
    result <- list(expression = expression, counts = counts, umap = umap,
      metadata = obj[[]][, c(cfg$sample_col, "group", cfg$subtype_col), drop = FALSE], catalogue = catalogue,
      reference = reference, color_limits = c(0, ceiling(max(expression))),
      point_size = Seurat::AutoPointSize(as.data.frame(umap), raster = FALSE),
      source_md5 = tools::md5sum(cfg$object), config = cfg, session = sessionInfo())
    qs2::qs_save(result, cfg$cache)
    rm(obj, expression, counts); gc(FALSE)
  }
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  all_pages <- sort(unique(result$catalogue$page))
  if (is.null(pages)) pages <- all_pages
  stopifnot(all(pages %in% all_pages))
  message("FERRDB_READY source_symbols=", length(unique(result$reference$records$Symbol)),
    " plotted_mouse_genes=", nrow(result$catalogue), " all_zero=", sum(result$catalogue$detected_cells == 0),
    " cells=", nrow(result$umap), " pages=", length(all_pages), " color_max=", result$color_limits[2])
  um <- as.data.frame(result$umap); names(um) <- c("UMAP_1", "UMAP_2")
  for (page in pages) {
    stem <- file.path(cfg$output, sprintf("FerrDb_FeaturePlot_%03d", page))
    if (all(file.exists(paste0(stem, c(".png", ".pdf"))))) next
    catalog <- result$catalogue[result$catalogue$page == page, ]
    plots <- lapply(seq_len(nrow(catalog)), function(i) {
      gene <- catalog$gene[i]; d <- um; d$expression <- as.numeric(result$expression[gene, ])
      d <- d[order(d$expression), ]
      tags <- paste(toupper(substr(strsplit(catalog$roles[i], ";", fixed = TRUE)[[1]], 1, 1)), collapse = "/")
      label <- paste0(gene, " [", tags, "]", if (catalog$detected_cells[i] == 0) " (not detected)" else "")
      ferrdb_feature_panel(d, label, result, cfg)
    })
    p <- patchwork::wrap_plots(plots, ncol = cfg$ncol, guides = "collect")
    cnmf_save_plot(p, stem, cfg$width, cfg$height)
    message("FERRDB_PAGE_COMPLETE ", page, "/", length(all_pages), " ", paste(catalog$gene, collapse = ";"))
  }
  invisible(result$catalogue)
}

ferrdb_feature_panel <- function(d, label, result, cfg) {
  ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, color = expression)) +
    ggplot2::geom_point(size = result$point_size, shape = 16, alpha = 1) +
    ggplot2::coord_equal(xlim = range(result$umap[, 1]), ylim = range(result$umap[, 2])) +
    ggplot2::scale_color_gradientn(colors = cfg$colors, values = cfg$color_values,
      limits = result$color_limits, breaks = pretty(result$color_limits, n = 4), name = "RNA\nlog-normalized") +
    ggplot2::guides(color = ggplot2::guide_colorbar(display = "rectangles")) +
    ggplot2::labs(title = label, x = "UMAP 1", y = "UMAP 2") + cnmf_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5, face = "bold"))
}

plot_ferrdb_group_split <- function(cfg, pages = NULL) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  result <- qs2::qs_read(cfg$source_cache)
  stopifnot(identical(result$source_md5, tools::md5sum(names(result$source_md5))))
  um <- as.data.frame(result$umap); names(um) <- c("UMAP_1", "UMAP_2")
  um$group <- as.character(result$metadata$group)
  stopifnot(setequal(um$group, cfg$group_levels), identical(rownames(um), rownames(result$metadata)))
  catalogue <- result$catalogue
  catalogue$page <- ceiling(seq_len(nrow(catalogue)) / cfg$per_page)
  if (is.null(pages)) pages <- sort(unique(catalogue$page))
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
  message("SPLIT_READY genes=", nrow(catalogue), " pages=", max(catalogue$page), " CON=", sum(um$group == cfg$group_levels[1]),
    " EXP=", sum(um$group == cfg$group_levels[2]))
  for (page in pages) {
    stem <- file.path(cfg$output, sprintf("FerrDb_FeaturePlot_EXP_vs_CON_%03d", page))
    if (all(file.exists(paste0(stem, c(".png", ".pdf"))))) next
    cat <- catalogue[catalogue$page == page, ]; plots <- list()
    for (i in seq_len(nrow(cat))) {
      gene <- cat$gene[i]
      tags <- paste(toupper(substr(strsplit(cat$roles[i], ";", fixed = TRUE)[[1]], 1, 1)), collapse = "/")
      for (group in cfg$group_levels) {
        d <- um; d$expression <- as.numeric(result$expression[gene, ])
        d <- d[d$group == group, ]; d <- d[order(d$expression), ]
        label <- paste0(gene, " [", tags, "] | ", group,
          if (!any(d$expression > 0)) " (not detected)" else "")
        plots[[paste(gene, group)]] <- ferrdb_feature_panel(d, label, result, cfg)
      }
    }
    p <- patchwork::wrap_plots(plots, ncol = cfg$ncol, guides = "collect")
    cnmf_save_plot(p, stem, cfg$width, cfg$height)
    message("SPLIT_PAGE_COMPLETE ", page, "/", max(catalogue$page), " ", paste(cat$gene, collapse = ";"))
  }
  invisible(catalogue)
}
