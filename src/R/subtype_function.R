# 基于已验证原始UMI求和矩阵，按真实小鼠阻断的亚群-vs-rest全基因排名与GSEA。
run_subtype_function <- function(cfg, keys = names(cfg$lineages)) {
  reference <- cnmf_reference(cfg)
  extended <- prepare_subtype_gsea_reference(cfg)
  for (key in keys) {
    item <- cfg$lineages[[key]]
    input <- file.path(cfg$cache_dir, paste0(cfg$input_step, ".", key, "_EXP_vs_CON.qs2"))
    x <- qs2::qs_read(input)
    stopifnot(identical(x$source_md5, tools::md5sum(names(x$source_md5))))
    ann <- read.csv(item$annotation, stringsAsFactors = FALSE)
    stopifnot(setequal(names(x$results), ann$annotation))
    all_counts <- lapply(x$results, `[[`, "counts")
    stopifnot(all(vapply(all_counts, function(m) identical(dimnames(m), dimnames(all_counts[[1]])), logical(1))))
    total <- Reduce(`+`, all_counts)
    metadata <- lapply(x$results, function(r) if (is.null(r$all_samples)) r$samples else r$all_samples)
    samples <- metadata[[1]]$sample
    stopifnot(all(vapply(metadata, function(m) identical(m$sample, samples), logical(1))), identical(colnames(total), samples))
    total_cells <- Reduce(`+`, lapply(metadata, `[[`, "n_cells"))
    cache_path <- file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_function.qs2"))
    hashes <- tools::md5sum(c(input, cfg$reference_file, item$annotation))
    if (file.exists(cache_path)) {
      cache <- qs2::qs_read(cache_path)
      stopifnot(identical(cache$config, cfg), identical(cache$source_md5, hashes), identical(cache$reference$md5, reference$md5))
    } else cache <- list(results = list(), status = list(), config = cfg, source_md5 = hashes,
      reference = reference, session = sessionInfo())
    dir.create(item$output, recursive = TRUE, showWarnings = FALSE)
    for (j in seq_len(nrow(ann))) {
      subtype <- ann$annotation[j]
      if (subtype %in% names(cache$status)) {
        prior <- cache$results[[subtype]]
        if (is.null(prior$de)) next
        rank_ok <- !any(prior$de$F < 0) || isTRUE(prior$negative_F_clamped)
        precision_ok <- !anyNA(prior$gsea$pval) || isTRUE(prior$gsea_precision_checked)
        if (rank_ok && precision_ok) next
        message("Recomputing GSEA to resolve numerical precision/ranking: ", subtype)
      }
      target <- all_counts[[subtype]]; rest <- total - target
      m <- metadata[[subtype]]; m$rest_cells <- total_cells - m$n_cells
      m$used <- m$n_cells >= cfg$min_cells & m$rest_cells >= cfg$min_cells
      pairs <- sum(m$used)
      status <- data.frame(lineage = key, cluster = ann$cluster[j], subtype,
        status = if (pairs >= cfg$min_pairs) "tested" else "insufficient_paired_mice",
        n_mice = pairs, used_samples = paste(m$sample[m$used], collapse = ";"),
        source_conditions = paste(unique(m$group[m$used]), collapse = ";"),
        tested_genes = 0L, positive_markers = 0L, enriched_pathways = 0L, depleted_pathways = 0L)
      status <- cbind(status, as.data.frame(as.list(setNames(m$n_cells, paste0("cells_", m$sample)))))
      message(key, " / ", subtype, ": ", pairs, " paired mice")
      result <- list(target_counts = target, rest_counts = rest, samples = m)
      if (pairs >= cfg$min_pairs) {
        # 两种细胞集合来自同一只小鼠，是真实配对；CON/EXP不是本次比较因子。
        counts <- cbind(rest[, m$used, drop = FALSE], target[, m$used, drop = FALSE])
        info <- data.frame(sample = factor(rep(samples[m$used], 2), levels = samples[m$used]),
          population = factor(rep(c("Rest", "Target"), each = pairs), levels = c("Rest", "Target")))
        rownames(info) <- colnames(counts) <- paste(info$sample, info$population, sep = "__")
        design <- model.matrix(~ sample + population, info)
        stopifnot(qr(design)$rank == ncol(design), nrow(design) - ncol(design) == pairs - 1L)
        y <- edgeR::DGEList(counts, samples = info, group = info$population)
        # 按两种population确定表达过滤的重复数，而不是按每个小鼠固定效应的杠杆值过滤。
        keep <- edgeR::filterByExpr(y, group = info$population)
        stopifnot(sum(keep) >= cfg$min_genes)
        y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
        y <- edgeR::estimateDisp(y, design, robust = TRUE)
        fit <- edgeR::glmQLFit(y, design, robust = TRUE)
        test <- edgeR::glmQLFTest(fit, coef = "populationTarget")
        de <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
        de$gene <- rownames(de); de <- de[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
        seed <- cfg$seed + match(key, names(cfg$lineages)) * 100L + j
        base_gsea <- subtype_gsea(de, reference, cfg, seed)
        gsea_result <- extend_subtype_gsea(de, base_gsea, extended, cfg, seed)
        g <- gsea_result$gsea
        stopifnot(setequal(g$collection, cfg$collection_order))
        result <- c(result, list(de = de, gsea = g, library = y$samples, design = design,
          keep = keep, kegg_mapped_genes = gsea_result$kegg_mapped_genes,
          negative_F_clamped = TRUE, gsea_precision_checked = TRUE,
          logcpm = edgeR::cpm(y, log = TRUE, prior.count = 2),
          selected = select_subtype_gsea_terms(g, cfg)))
        status$tested_genes <- nrow(de)
        status$positive_markers <- sum(de$FDR < cfg$fdr & de$logFC >= cfg$logfc)
        status$enriched_pathways <- sum(g$padj < cfg$fdr & g$NES > 0, na.rm = TRUE)
        status$depleted_pathways <- sum(g$padj < cfg$fdr & g$NES < 0, na.rm = TRUE)
        directory <- file.path(item$output, sprintf("C%02d_%s", ann$cluster[j], gsub("[^A-Za-z0-9]+", "_", subtype)))
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(de[order(de$FDR, -de$logFC), ], file.path(directory, "Markers_subtype_vs_rest.csv"), bom = TRUE)
        data.table::fwrite(g[order(match(g$collection, cfg$collection_order), g$padj), ], file.path(directory, "GSEA_pathways.csv"), bom = TRUE)
        plot_subtype_gsea_collections(g, cfg, file.path(directory, "GSEA_subtype_function"), subtype)
        message("  genes=", nrow(de), "; target-enriched pathways=", status$enriched_pathways)
      }
      cache$results[[subtype]] <- result; cache$status[[subtype]] <- status
      qs2::qs_save(cache, cache_path)
      data.table::fwrite(do.call(rbind, cache$status), file.path(item$output, "subtype_function_summary.csv"), bom = TRUE)
    }
    message(key, " SUBTYPE_FUNCTION_COMPLETE: ", length(cache$status))
  }
  invisible(TRUE)
}

# 将人工审阅后的功能描述合并到已有核心汇总表；所有证据必须是目标亚群正向显著GSEA。
annotate_subtype_function <- function(cfg, annotation_path) {
  annotations <- data.table::fread(annotation_path, data.table = FALSE, encoding = "UTF-8", na.strings = NULL)
  stopifnot(!anyDuplicated(paste(annotations$lineage, annotations$cluster)))
  for (key in names(cfg$lineages)) {
    x <- qs2::qs_read(file.path(cfg$cache_dir, paste0(cfg$step, ".", key, "_function.qs2")))
    out <- do.call(rbind, x$status)
    out$unestimable_pathways <- vapply(x$results, function(r)
      if (is.null(r$gsea)) NA_integer_ else sum(is.na(r$gsea$pval)), integer(1))
    ann <- annotations[annotations$lineage == key, ]
    stopifnot(setequal(out$cluster, ann$cluster))
    ann <- ann[match(out$cluster, ann$cluster), ]
    out$main_function <- ann$main_function
    out$interpretation_limit <- ann$interpretation_limit
    out$evidence_pathways <- out$leading_edge_genes <- ""
    for (i in seq_len(nrow(out))) {
      a <- x$results[[out$subtype[i]]]
      if (out$status[i] != "tested") {
        stopifnot(!nzchar(ann$evidence_pathways[i]))
        next
      }
      ids <- strsplit(ann$evidence_pathways[i], ";", fixed = TRUE)[[1]]
      index <- match(ids, a$gsea$pathway)
      stopifnot(length(ids) > 0, !anyNA(index))
      g <- a$gsea[index, ]
      stopifnot(all(g$NES > 0), all(g$padj < cfg$fdr))
      out$evidence_pathways[i] <- paste(sprintf("%s [%s; NES=%.2f; FDR=%.3g]",
        g$pathway_name, g$pathway, g$NES, g$padj), collapse = " | ")
      genes <- unique(unlist(strsplit(g$leadingEdge, ";", fixed = TRUE)))
      d <- a$de[a$de$gene %in% genes & a$de$logFC > 0, ]
      out$leading_edge_genes[i] <- paste(head(d$gene[order(d$FDR, -d$logFC)], 10), collapse = ";")
    }
    data.table::fwrite(out, file.path(cfg$lineages[[key]]$output, "subtype_function_summary.csv"), bom = TRUE)
    message(key, " FUNCTION_ANNOTATION_COMPLETE: ", nrow(out))
  }
}
