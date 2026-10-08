# 复用原细分标签归并大类；以经典 marker 复核，不重新聚类或全基因找 marker。
step_broad_annotation_evidence <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$broad_annotation
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  mapping <- read.delim(file.path(PROJECT_DIR, settings$mapping), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(mapping$original_cell_type))
  sc <- run_step("Read integrated object", qs::qread(file.path(PROJECT_DIR, settings$input), nthreads = 4))
  metadata <- sc[[]]
  metadata$broad_cell_type <- factor(mapping$broad_cell_type[match(as.character(metadata$cell_type), mapping$original_cell_type)],
                                     levels = unlist(settings$broad_levels))
  stopifnot(!anyNA(metadata$broad_cell_type))
  panel <- read.delim(file.path(PROJECT_DIR, settings$marker_panel), stringsAsFactors = FALSE)
  features <- unique(unlist(strsplit(panel$genes, ";", fixed = TRUE)))
  missing <- setdiff(features, rownames(sc))
  writeLines(missing, file.path(output, "marker_genes_not_found.txt"))
  features <- intersect(features, rownames(sc))
  expression <- LayerData(sc, assay = "RNA", layer = "data")[features, ]
  umap <- Embeddings(sc, "umap")
  rm(sc)
  invisible(gc())
  linear <- expression
  linear@x <- expm1(linear@x)
  all_stats <- list()
  for (column in c("broad_cell_type", settings$cluster_col)) {
    group <- metadata[[column]]
    if (column == settings$cluster_col) group <- factor(as.character(group), levels = sort(unique(as.integer(as.character(group)))))
    sizes <- as.integer(table(group))
    membership <- Matrix::sparse.model.matrix(~ 0 + group)
    average <- sweep(as.matrix(linear %*% membership), 2, sizes, "/")
    pct <- sweep(as.matrix((expression > 0) %*% membership), 2, sizes, "/") * 100
    stats <- expand.grid(gene = features, label = levels(group), stringsAsFactors = FALSE)
    stats$average_expression <- as.vector(average)
    stats$pct_expressed <- as.vector(pct)
    write.table(stats, file.path(output, paste0("canonical_marker_by_", column, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    all_stats[[column]] <- stats
    stats$scaled_average <- ave(log1p(stats$average_expression), stats$gene, FUN = function(x)
      if (sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x)))
    stats$gene <- factor(stats$gene, levels = features)
    stats$label <- factor(stats$label, levels = rev(levels(group)))
    p <- ggplot(stats, aes(gene, label)) +
      geom_point(aes(size = pct_expressed, color = pmax(-2.5, pmin(2.5, scaled_average)))) +
      scale_size_area(max_size = 3.6, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      labs(x = NULL, y = NULL, size = "Expressed (%)", color = "Scaled mean")
    save_plot(hongf_dotplot_style(p), file.path(output, paste0("DotPlot_canonical_", column)),
              width = 22, height = if (column == "broad_cell_type") 4 else 8)
  }
  saveRDS(list(metadata = metadata, umap = umap, marker_expression = expression,
               marker_stats = all_stats, input = settings$input), file.path(PROJECT_DIR, settings$evidence_object))
  counts <- as.data.frame(table(broad_cell_type = metadata$broad_cell_type))
  names(counts)[2] <- "cells"
  counts$percent <- counts$cells / nrow(metadata) * 100
  write.table(counts, file.path(output, "broad_cell_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  print(counts, row.names = FALSE)
  log_msg("BROAD_EVIDENCE_OK: ", nrow(metadata), " cells; ", nlevels(metadata$broad_cell_type), " broad classes")
  invisible(output)
}
