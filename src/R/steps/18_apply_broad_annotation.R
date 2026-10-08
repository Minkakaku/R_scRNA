# 将经 marker/文献复核的原标签映射写入完整对象，保留原细分、cluster 和坐标。
step_apply_broad_annotation <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$broad_annotation
  output <- file.path(PROJECT_DIR, settings$output_dir)
  destination <- file.path(PROJECT_DIR, settings$annotated_object)
  stopifnot(!file.exists(destination), !file.exists(paste0(destination, ".partial")))
  mapping <- read.delim(file.path(PROJECT_DIR, settings$mapping), stringsAsFactors = FALSE)
  review <- read.delim(file.path(PROJECT_DIR, settings$review_table), stringsAsFactors = FALSE)
  references <- jsonlite::fromJSON(file.path(PROJECT_DIR, settings$literature_records))
  stopifnot(!anyDuplicated(mapping$original_cell_type), !anyDuplicated(review$broad_cell_type),
            all(unique(unlist(strsplit(review$reference_pmids, ";", fixed = TRUE))) %in% references$pmid))
  sc <- run_step("Read complete Harmony object", qs::qread(file.path(PROJECT_DIR, settings$input), nthreads = 4))
  original_metadata <- sc[[]]
  original_cells <- colnames(sc)
  original_genes <- rownames(sc)
  reduction_hash <- vapply(sc@reductions, digest::digest, character(1), algo = "xxhash64")
  graph_hash <- vapply(sc@graphs, digest::digest, character(1), algo = "xxhash64")
  layer_hashes <- function(x) vapply(c("counts", "data", "scale.data"), function(layer)
    digest::digest(LayerData(x, assay = "RNA", layer = layer), algo = "xxhash64"), character(1))
  expression_hash <- layer_hashes(sc)
  major <- mapping$broad_cell_type[match(as.character(sc$cell_type), mapping$original_cell_type)]
  stopifnot(!anyNA(major), setequal(unique(major), review$broad_cell_type))
  sc$broad_cell_type <- factor(major, levels = unlist(settings$broad_levels))
  sc$broad_annotation_status <- review$review_status[match(major, review$broad_cell_type)]
  sc$broad_annotation_source <- settings$annotation_method
  Idents(sc) <- "broad_cell_type"
  sc@misc$broad_annotation <- list(species = settings$species, method = settings$annotation_method,
    source_column = "cell_type", mapping = mapping, review = review, literature = references,
    input = settings$input, expression_hash = expression_hash, new_clustering = FALSE)
  expected_metadata <- sc[[]]
  run_step("Save complete annotated QS", qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4))
  rm(sc)
  invisible(gc())
  check <- run_step("Read back annotated QS", qs::qread(paste0(destination, ".partial"), nthreads = 4))
  stopifnot(isTRUE(validObject(check, test = TRUE)), identical(colnames(check), original_cells),
            identical(rownames(check), original_genes), identical(check[[]], expected_metadata),
            identical(check[[]][, names(original_metadata), drop = FALSE], original_metadata),
            identical(layer_hashes(check), expression_hash),
            identical(vapply(check@reductions, digest::digest, character(1), algo = "xxhash64"), reduction_hash),
            identical(vapply(check@graphs, digest::digest, character(1), algo = "xxhash64"), graph_hash),
            identical(as.character(Idents(check)), as.character(check$broad_cell_type)),
            file.rename(paste0(destination, ".partial"), destination))
  log_msg("BROAD_ANNOTATED_QS_OK: all cells, RNA layers, original metadata, reductions and graphs preserved")
  write.table(mapping, file.path(output, "original_to_broad_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(review, file.path(output, "broad_annotation_review.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  connection <- gzfile(file.path(output, "cell_metadata_broad_annotation.tsv.gz"), "wt")
  write.table(data.frame(cell_id = colnames(check), check[[]]), connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(connection)
  for (column in c("sample", "group", settings$cluster_col, "cell_type")) {
    tab <- as.data.frame(table(broad_cell_type = check$broad_cell_type, label = check[[column, drop = TRUE]]))
    names(tab) <- c("broad_cell_type", column, "cells")
    write.table(tab, file.path(output, paste0("broad_by_", column, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # 仅在绘图内存对象中添加显示标签；磁盘完整对象采用稳定的大类字段。
  check$broad_plot_label <- factor(gsub("_", " ", as.character(check$broad_cell_type)),
                                   levels = gsub("_", " ", unlist(settings$broad_levels)))
  p <- npg_dimplot(check, reduction = "umap", group.by = "broad_plot_label", label = TRUE,
                   repel = TRUE, raster = FALSE, shuffle = TRUE, seed = SEED) +
    labs(title = "Major cell types") + theme_classic(base_size = 11)
  save_umap_plot(p, file.path(output, "UMAP_broad_cell_types"), width = 11, height = 8)
  for (column in c("sample", "group")) {
    p <- npg_dimplot(check, reduction = "umap", group.by = "broad_plot_label", split.by = column,
                    ncol = if (column == "sample") 3 else 2, raster = FALSE,
                    shuffle = TRUE, seed = SEED) + labs(title = "Major cell types") + theme_classic(base_size = 11)
    save_umap_plot(p, file.path(output, paste0("UMAP_broad_", column, "_split")),
                  width = 16, height = if (column == "sample") 10 else 6)
  }
  stats <- read.delim(file.path(output, "canonical_marker_by_broad_cell_type.tsv"), stringsAsFactors = FALSE)
  selected <- unlist(settings$overview_markers)
  stats <- stats[stats$gene %in% selected, ]
  stats$scaled_average <- ave(log1p(stats$average_expression), stats$gene, FUN = function(x)
    if (sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x)))
  stats$gene <- factor(stats$gene, levels = selected)
  stats$label <- factor(gsub("_", " ", stats$label), levels = rev(gsub("_", " ", unlist(settings$broad_levels))))
  p <- ggplot(stats, aes(gene, label)) +
    geom_point(aes(size = pct_expressed, color = pmax(-2.5, pmin(2.5, scaled_average)))) +
    scale_size_area(max_size = 4.5, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    labs(x = NULL, y = NULL, size = "Expressed (%)", color = "Scaled mean")
  save_plot(hongf_dotplot_style(p), file.path(output, "DotPlot_broad_marker_evidence"), width = 12, height = 4)
  yaml::write_yaml(list(stage = "complete", cells = ncol(check), genes = nrow(check), broad_classes = nlevels(check$broad_cell_type),
    original_labels = length(unique(check$cell_type)), clusters = length(unique(check$harmony_cluster)),
    method = settings$annotation_method, annotated_object = settings$annotated_object,
    preserved_RNA_layers = TRUE, preserved_reductions_graphs = TRUE), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "11_apply_broad_annotation_sessionInfo.txt"))
  log_msg("BROAD_ANNOTATION_COMPLETE: six classes; original labels and all coordinates preserved")
  invisible(destination)
}
