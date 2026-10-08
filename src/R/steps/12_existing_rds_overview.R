# Inspect supplied Seurat objects without changing or re-saving them.
step_inspect_existing_rds <- function() {
  stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
  output <- file.path(PROJECT_DIR, CONFIG$object_overview$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  objects <- list()
  for (id in names(CONFIG$input_objects)) {
    path <- file.path(PROJECT_DIR, CONFIG$input_objects[[id]])
    message("Reading: ", id)
    object <- readRDS(path)
    stopifnot(inherits(object, "Seurat"))
    metadata <- object[[]]
    cells <- colnames(object)
    stopifnot(identical(rownames(metadata), cells), !anyDuplicated(cells))
    reductions <- SeuratObject::Reductions(object)
    umaps <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
    embeddings <- setNames(lapply(umaps, function(reduction) {
      x <- SeuratObject::Embeddings(object[[reduction]])
      stopifnot(ncol(x) >= 2L, all(is.finite(x)), all(cells %in% rownames(x)))
      x[cells, 1:2, drop = FALSE]
    }), umaps)
    assays <- vapply(names(object@assays), function(assay) {
      paste0(assay, " [", nrow(object[[assay]]), " features; ",
             paste(SeuratObject::Layers(object[[assay]]), collapse = ","), "]")
    }, character(1))
    metadata$active_ident <- as.character(SeuratObject::Idents(object))
    summary <- data.frame(
      object = id, cells = ncol(object), features = nrow(object),
      default_assay = SeuratObject::DefaultAssay(object),
      assays = paste(assays, collapse = "; "),
      reductions = paste(reductions, collapse = "; "),
      umaps = paste(umaps, collapse = "; "),
      commands = paste(names(object@commands), collapse = "; "),
      stringsAsFactors = FALSE
    )
    objects[[id]] <- list(summary = summary, metadata = metadata, embeddings = embeddings)
    print(summary[, 1:4], row.names = FALSE)
    cat("Assays:", summary$assays, "\nReductions:", summary$reductions, "\n")
    cat("Metadata:", paste(names(metadata), collapse = ", "), "\n")
    for (column in names(metadata)) {
      values <- metadata[[column]]
      if (length(unique(values)) <= 30L) {
        cat("  ", column, ": ", sep = "")
        counts <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
        cat(paste(paste(names(counts), counts, sep = "="), collapse = "; "), "\n")
      }
    }
    rm(object)
    invisible(gc())
  }
  summary <- do.call(rbind, lapply(objects, `[[`, "summary"))
  write.table(summary, file.path(output, "object_summary.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(objects, file.path(PROJECT_DIR, CONFIG$object_overview$plot_data))
  message("Inspection complete: ", output)
  invisible(objects)
}

# Plot the stored coordinates, retaining every cell and consistent label colors.
step_plot_existing_umaps <- function() {
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$object_overview
  objects <- readRDS(file.path(PROJECT_DIR, settings$plot_data))
  output <- file.path(PROJECT_DIR, settings$output_dir)
  columns <- c(cell_type = settings$cell_type_col, sample = settings$sample_col,
               group = settings$group_col)
  npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
           "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
  rgb <- grDevices::col2rgb(npg) / 255
  extended <- c(npg, grDevices::rgb(t(rgb * .65 + .35)),
                grDevices::rgb(t(rgb * .70)), grDevices::rgb(t(rgb * .50 + .50)))
  shared <- Reduce(intersect, lapply(objects, function(x) rownames(x$metadata)))
  shared_labels <- sort(unique(as.character(objects[[1]]$metadata[shared, columns[["cell_type"]]])))
  palettes <- lapply(names(columns), function(kind) {
    labels <- sort(unique(unlist(lapply(objects, function(x) as.character(x$metadata[[columns[[kind]]]])))))
    if (kind == "cell_type") labels <- c(shared_labels, setdiff(labels, shared_labels))
    stopifnot(length(labels) <= length(extended), !anyNA(labels))
    setNames(extended[seq_along(labels)], labels)
  })
  names(palettes) <- names(columns)
  frames <- lapply(objects, function(x) {
    embedding <- x$embeddings[[settings$reduction]]
    stopifnot(!is.null(embedding), identical(rownames(embedding), rownames(x$metadata)))
    data.frame(UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2],
               x$metadata[, columns, drop = FALSE], check.names = FALSE)
  })
  xlim <- range(unlist(lapply(frames, `[[`, "UMAP_1")))
  ylim <- range(unlist(lapply(frames, `[[`, "UMAP_2")))
  make_plot <- function(id, kind) {
    frame <- frames[[id]]
    set.seed(SEED)
    frame <- frame[sample.int(nrow(frame)), , drop = FALSE]
    frame$label <- factor(as.character(frame[[columns[[kind]]]]), levels = names(palettes[[kind]]))
    stopifnot(!anyNA(frame$label), nrow(frame) == objects[[id]]$summary$cells)
    ggplot(frame, aes(UMAP_1, UMAP_2, color = label)) +
      geom_point(size = .20, stroke = 0, alpha = .9) +
      scale_color_manual(values = palettes[[kind]], drop = TRUE) +
      coord_fixed(xlim = xlim, ylim = ylim) +
      labs(title = paste0(id, " (n = ", format(nrow(frame), big.mark = ","), ")"),
           x = "UMAP 1", y = "UMAP 2", color = NULL) +
      guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 8),
            legend.key.height = grid::unit(.34, "cm"),
            legend.key.width = grid::unit(.35, "cm"),
            plot.margin = margin(8, 8, 8, 8))
  }
  plots <- lapply(names(columns), function(kind) {
    setNames(lapply(names(objects), make_plot, kind = kind), names(objects))
  })
  names(plots) <- names(columns)
  exports <- lapply(plots, patchwork::wrap_plots, ncol = 2)
  names(exports) <- paste0("UMAP_", names(exports), "_overview")
  exports <- c(exports, setNames(plots$cell_type, paste0("UMAP_", names(objects), "_cell_type")))
  for (name in names(exports)) {
    combined <- grepl("_overview$", name)
    width <- if (combined) 16 else 8
    height <- if (combined) 11 else 5.5
    stem <- file.path(output, name)
    ggsave(paste0(stem, ".png"), exports[[name]], width = width, height = height,
           dpi = 600, bg = "white", device = ragg::agg_png)
    ggsave(paste0(stem, ".pdf"), exports[[name]], width = width, height = height,
           bg = "white", device = grDevices::cairo_pdf)
    message("Saved: ", name)
  }
  counts <- do.call(rbind, lapply(names(objects), function(id) {
    do.call(rbind, lapply(names(columns), function(kind) {
      values <- table(as.character(objects[[id]]$metadata[[columns[[kind]]]]))
      stopifnot(sum(values) == objects[[id]]$summary$cells)
      data.frame(object = id, field = kind, label = names(values), cells = as.integer(values))
    }))
  }))
  summary <- do.call(rbind, lapply(names(objects), function(id) {
    x <- objects[[id]]
    m <- x$metadata
    data.frame(x$summary, samples = length(unique(m[[columns[["sample"]]]])),
               cell_types = length(unique(m[[columns[["cell_type"]]]])),
               shared_cells = length(shared),
               median_genes = median(m$nFeature_RNA),
               median_counts = median(m$nCount_RNA),
               median_mt_percent = median(m$pct_counts_mt),
               predicted_doublets = sum(m$predicted_doublet))
  }))
  overlaps <- do.call(rbind, lapply(combn(names(objects), 2, simplify = FALSE), function(pair) {
    a <- objects[[pair[1]]]; b <- objects[[pair[2]]]
    cells <- intersect(rownames(a$metadata), rownames(b$metadata))
    same_metadata <- all(vapply(columns, function(column) {
      identical(as.character(a$metadata[cells, column]), as.character(b$metadata[cells, column]))
    }, logical(1)))
    difference <- max(abs(a$embeddings[[settings$reduction]][cells, ] -
                           b$embeddings[[settings$reduction]][cells, ]))
    data.frame(object_1 = pair[1], object_2 = pair[2], shared_cells = length(cells),
               same_labels_sample_group = same_metadata, max_umap_difference = difference)
  }))
  for (name in c("counts", "summary", "overlaps")) {
    write.table(get(name), file.path(output, paste0("object_", name, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  }
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "02_umap_sessionInfo.txt"))
  message("UMAP_CHECK_OK: ", length(exports), " PNG/PDF pairs; cells retained; original coordinates used.")
  invisible(summary)
}
