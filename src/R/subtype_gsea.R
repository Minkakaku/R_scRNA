# 官方KEGG小鼠通路与MSigDB GO全集的GSEA扩展；不重复拟合差异模型。
prepare_subtype_gsea_reference <- function(cfg) {
  if (file.exists(cfg$reference_file)) return(readRDS(cfg$reference_file))
  go <- readRDS(cfg$go)
  sets <- lapply(c("GO:CC", "GO:MF"), function(s) {
    d <- go[go$gs_subcollection == s, ]; lapply(split(d$db_gene_symbol, d$gs_name), unique)
  }); names(sets) <- c("GO_CC", "GO_MF")
  # KEGGREST直接读取物种通路及NCBI Gene转换；不以其他数据库代称KEGG。
  paths <- KEGGREST::keggList("pathway", cfg$organism)
  links <- KEGGREST::keggLink("pathway", cfg$organism)
  conv <- KEGGREST::keggConv("ncbi-geneid", cfg$organism)
  ids <- sub("^ncbi-geneid:", "", unname(conv[match(names(links), names(conv))]))
  path_ids <- sub("^path:", "", unname(links))
  sets$KEGG <- lapply(split(ids[!is.na(ids)], path_ids[!is.na(ids)]), unique)
  names(paths) <- sub("^path:", "", names(paths))
  db <- getExportedValue(cfg$orgdb, cfg$orgdb)
  mapping <- AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, "SYMBOL"),
    columns = "ENTREZID", keytype = "SYMBOL")
  mapping <- unique(mapping[!is.na(mapping$ENTREZID), ])
  # 仅保留双向一对一的精确SYMBOL映射，避免同一个基因被重复加权。
  mapping <- mapping[!duplicated(mapping$SYMBOL) & !duplicated(mapping$SYMBOL, fromLast = TRUE) &
    !duplicated(mapping$ENTREZID) & !duplicated(mapping$ENTREZID, fromLast = TRUE), ]
  descriptions <- c(paths, setNames(go$gs_description, go$gs_name))
  descriptions <- descriptions[!duplicated(names(descriptions))]
  out <- list(sets = sets, descriptions = descriptions, gene_map = mapping,
    kegg_raw = list(pathways = paths, links = links, ncbi_conversion = conv),
    go_version = unique(go$db_version), go_md5 = tools::md5sum(cfg$go),
    orgdb_package = cfg$orgdb, orgdb_version = as.character(utils::packageVersion(cfg$orgdb)),
    organism = cfg$organism, retrieved = Sys.time(),
    urls = paste0("https://rest.kegg.jp/", c(paste0("list/pathway/", cfg$organism),
      paste0("link/pathway/", cfg$organism), paste0("conv/ncbi-geneid/", cfg$organism))))
  stopifnot(length(sets$KEGG) > 100L, length(sets$GO_CC) > 100L, length(sets$GO_MF) > 100L,
    !anyNA(names(sets$KEGG)), !anyDuplicated(mapping$SYMBOL), !anyDuplicated(mapping$ENTREZID))
  saveRDS(out, cfg$reference_file)
  message("REFERENCE_READY: ", paste(names(sets), lengths(sets), collapse = "; "))
  out
}

extend_subtype_gsea <- function(de, old_gsea, reference, cfg, seed) {
  stopifnot(all(is.finite(de$F)), min(de$F) >= -1e-6)
  stat <- setNames(sign(de$logFC) * sqrt(pmax(de$F, 0)), de$gene)
  out <- list(old_gsea)
  mapped <- reference$gene_map[reference$gene_map$SYMBOL %in% names(stat), ]
  collections <- if (is.null(cfg$additional_collections)) c("KEGG", "GO_CC", "GO_MF") else cfg$additional_collections
  for (collection in collections) {
    if (collection == "KEGG") {
      ranks <- setNames(stat[mapped$SYMBOL], mapped$ENTREZID)
    } else ranks <- stat
    ranks <- sort(ranks, decreasing = TRUE)
    stopifnot(!anyDuplicated(names(ranks)), all(is.finite(ranks)))
    # 此扩展文件也被独立入口加载，不依赖上游差异分析模块中的函数。
    for (permutations in c(1000L, 10000L, 100000L)) {
      set.seed(seed)
      g <- as.data.frame(fgsea::fgseaMultilevel(reference$sets[[collection]], ranks,
        minSize = cfg$set_size[1], maxSize = cfg$set_size[2], eps = 0,
        nPermSimple = permutations, sampleSize = 101L, nproc = 1L))
      missing <- sum(is.na(g$pval))
      if (!missing) {
        if (permutations > 1000L) message("GSEA_PRECISION_OK nPermSimple=", permutations)
        break
      }
      message("GSEA_PRECISION nPermSimple=", permutations, " missing_p=", missing)
    }
    if (missing) warning("GSEA retains ", missing, " unestimable pathways at nPermSimple=100000")
    if (collection == "KEGG") g$leadingEdge <- lapply(g$leadingEdge, function(z)
      mapped$SYMBOL[match(z, mapped$ENTREZID)])
    stopifnot(!anyNA(unlist(g$leadingEdge)))
    g$leadingEdge <- vapply(g$leadingEdge, paste, character(1), collapse = ";")
    g$collection <- collection
    out[[collection]] <- g
  }
  g <- do.call(rbind, out)
  # 展示标准通路名称；保留pathway字段的GO/MSigDB/KEGG标识用于追溯。
  g$pathway_name <- gsub("_", " ", sub("^(HALLMARK|GOBP|GOCC|GOMF)_", "", g$pathway))
  which_kegg <- g$collection == "KEGG"
  g$pathway_name[which_kegg] <- sub(" - Mus musculus.*$", "", reference$descriptions[g$pathway[which_kegg]])
  stopifnot(!anyNA(g$pathway_name))
  list(gsea = g, kegg_mapped_genes = mapped, ranked_genes = length(stat))
}

select_subtype_gsea_terms <- function(gsea, cfg) {
  g <- gsea[gsea$collection %in% cfg$collection_order & !is.na(gsea$padj) &
    gsea$padj < cfg$fdr & gsea$NES != 0, ]
  if (!nrow(g)) return(g)
  selected <- lapply(cfg$collection_order, function(collection) {
    d <- g[g$collection == collection, ]
    d <- d[order(d$padj, -abs(d$NES), d$pathway), ]
    initial <- c(head(which(d$NES > 0), cfg$display_terms_per_direction),
      head(which(d$NES < 0), cfg$display_terms_per_direction))
    # 新配置允许完成两个方向初选后，按同一FDR顺序补足总数；旧调用不补位。
    if (!is.null(cfg$display_max_total)) {
      stopifnot(cfg$display_max_total >= 2L * cfg$display_terms_per_direction)
      initial <- c(initial, head(setdiff(seq_len(nrow(d)), initial), cfg$display_max_total - length(initial)))
    }
    d[initial, , drop = FALSE]
  })
  do.call(rbind, selected)
}

plot_subtype_gsea_collections <- function(gsea, cfg, stem, label) {
  chosen <- select_subtype_gsea_terms(gsea, cfg)
  g <- chosen
  if (nrow(g)) {
    chosen$collection <- factor(chosen$collection, levels = cfg$collection_order)
    chosen$pathway <- factor(chosen$pathway, levels = unique(chosen$pathway[order(chosen$collection, chosen$NES)]))
    labels <- setNames(vapply(tolower(chosen$pathway_name), function(s)
      paste(strwrap(s, width = 54), collapse = "\n"), character(1)), chosen$pathway)
    chosen$direction <- factor(ifelse(chosen$NES > 0, cfg$group_levels[2], cfg$group_levels[1]), levels = cfg$group_levels)
    chosen$strength <- pmin(-log10(chosen$padj), 10)
    p <- ggplot2::ggplot(chosen, ggplot2::aes(NES, pathway, color = direction, size = strength)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#AAAAAA") + ggplot2::geom_point() +
      ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_y_discrete(labels = labels) +
      ggplot2::scale_color_manual(values = cfg$group_colors, drop = FALSE, name = "Enriched in") +
      ggplot2::scale_size_continuous(range = c(2, 5), limits = c(0, 10), breaks = c(2, 5, 10), labels = c("2", "5", ">=10")) +
      ggplot2::labs(title = label, x = paste0("GSEA NES (", cfg$group_levels[2], " positive)"), y = NULL, size = "-log10(FDR)") +
      cnmf_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#ECECEC", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.4), axis.line = ggplot2::element_blank())
    height <- max(5, sum(lengths(strsplit(labels, "\n", fixed = TRUE))) * 0.25 + 1.3)
  } else {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No Hallmark / KEGG / GO pathway at BH FDR < 0.05", size = 4) +
      ggplot2::labs(title = label) + ggplot2::theme_void(base_family = "Arial")
    height <- 4
  }
  cnmf_save_plot(p, stem, 10, height)
}

run_subtype_gsea_extension <- function(cfg, keys = names(cfg$lineages)) {
  reference <- prepare_subtype_gsea_reference(cfg)
  for (key in keys) {
    item <- cfg$lineages[[key]]
    de_path <- file.path(cfg$cache_dir, paste0(cfg$de_step, ".", key, "_EXP_vs_CON.qs2"))
    previous <- qs2::qs_read(de_path)
    cache_path <- file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_GSEA.qs2"))
    hashes <- tools::md5sum(c(de_path, cfg$reference_file))
    if (file.exists(cache_path)) {
      cache <- qs2::qs_read(cache_path)
      stopifnot(identical(cache$source_md5, hashes), identical(cache$config, cfg))
    } else cache <- list(results = list(), source_md5 = hashes, config = cfg,
      reference_file = cfg$reference_file, session = sessionInfo())
    annotation <- read.csv(item$annotation)
    for (j in seq_len(nrow(annotation))) {
      subtype <- annotation$annotation[j]
      if (is.null(previous$results[[subtype]]$de) || subtype %in% names(cache$results)) next
      p <- previous$results[[subtype]]
      message(key, " / ", subtype, ": extending GSEA")
      result <- extend_subtype_gsea(p$de, p$gsea, reference, cfg,
        cfg$seed + match(key, names(cfg$lineages)) * 100L + j)
      directory <- file.path(item$output, sprintf("C%02d_%s", annotation$cluster[j], gsub("[^A-Za-z0-9]+", "_", subtype)))
      g <- result$gsea
      data.table::fwrite(g[order(match(g$collection, cfg$collection_order), g$padj), ],
        file.path(directory, "GSEA_pathways.csv"), bom = TRUE)
      plot_subtype_gsea_collections(g, cfg, file.path(directory, "Fig2_pseudobulk_GSEA"), subtype)
      cache$results[[subtype]] <- result
      qs2::qs_save(cache, cache_path)
      message("  KEGG mapped=", nrow(result$kegg_mapped_genes), "/", result$ranked_genes,
        "; significant=", sum(g$padj < cfg$fdr, na.rm = TRUE))
    }
    summary <- do.call(rbind, previous$status)
    for (subtype in names(cache$results)) summary$gsea_significant[summary$subtype == subtype] <-
      sum(cache$results[[subtype]]$gsea$padj < cfg$fdr, na.rm = TRUE)
    data.table::fwrite(summary, file.path(item$output, "pseudobulk_summary.csv"), bom = TRUE)
    message(key, " GSEA_EXTENSION_COMPLETE: ", length(cache$results), " subtypes")
  }
  invisible(TRUE)
}

refresh_subtype_gsea_display <- function(cfg, keys = names(cfg$lineages)) {
  for (key in keys) {
    item <- cfg$lineages[[key]]
    input <- file.path(cfg$cache_dir, paste0(cfg$input_gsea_step, ".", key, "_GSEA.qs2"))
    r <- qs2::qs_read(input)
    stopifnot(identical(r$source_md5, tools::md5sum(names(r$source_md5))))
    annotation <- read.csv(item$annotation)
    summary <- data.table::fread(file.path(item$output, "pseudobulk_summary.csv"), data.table = FALSE)
    selected <- list()
    for (j in seq_len(nrow(annotation))) {
      subtype <- annotation$annotation[j]
      if (is.null(r$results[[subtype]])) next
      g <- r$results[[subtype]]$gsea
      g <- g[g$collection %in% cfg$collection_order, ]
      directory <- file.path(item$output, sprintf("C%02d_%s", annotation$cluster[j], gsub("[^A-Za-z0-9]+", "_", subtype)))
      data.table::fwrite(g[order(match(g$collection, cfg$collection_order), g$padj), ],
        file.path(directory, "GSEA_pathways.csv"), bom = TRUE)
      plot_subtype_gsea_collections(g, cfg, file.path(directory, "Fig2_pseudobulk_GSEA"), subtype)
      selected[[subtype]] <- select_subtype_gsea_terms(g, cfg)
      summary$gsea_significant[summary$subtype == subtype] <- sum(g$padj < cfg$fdr, na.rm = TRUE)
    }
    qs2::qs_save(list(selected = selected, config = cfg, source_md5 = tools::md5sum(input)),
      file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_GSEA_display.qs2")))
    data.table::fwrite(summary, file.path(item$output, "pseudobulk_summary.csv"), bom = TRUE)
    message(key, " GSEA_DISPLAY_COMPLETE: ", length(selected), " subtypes")
  }
}
