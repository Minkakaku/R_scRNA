# Deduplicate matching cells, retain original expression and provenance, write qs.
step_merge_existing_rds <- function() {
  suppressPackageStartupMessages(library(Seurat))
  stopifnot(requireNamespace("qs", quietly = TRUE))
  output <- file.path(PROJECT_DIR, CONFIG$merged_output)
  stopifnot(!file.exists(output), !file.exists(paste0(output, ".partial")))
  cache <- readRDS(file.path(PROJECT_DIR, CONFIG$object_overview$plot_data))
  source_cells <- lapply(cache, function(x) rownames(x$metadata))
  shared <- Reduce(intersect, source_cells)
  all_cells <- unique(unlist(source_cells, use.names = FALSE))
  objects <- list()
  seen <- character()
  checks <- list()
  pca_compatible <- TRUE
  for (id in names(CONFIG$input_objects)) {
    message("Reading for merge: ", id)
    object <- readRDS(file.path(PROJECT_DIR, CONFIG$input_objects[[id]]))
    stopifnot(inherits(object, "Seurat"), identical(colnames(object), source_cells[[id]]),
              identical(names(object@assays), "RNA"),
              setequal(Layers(object[["RNA"]]), c("counts", "data")))
    current <- lapply(c("counts", "data"), function(layer) {
      LayerData(object, assay = "RNA", layer = layer)[, shared, drop = FALSE]
    })
    umap <- Embeddings(object[[CONFIG$object_overview$reduction]])
    pca <- Embeddings(object[["pca"]])
    if (!length(objects)) {
      features <- rownames(object)
      reference <- current
      reference_umap <- umap[shared, , drop = FALSE]
      reference_pca <- pca[shared, , drop = FALSE]
    }
    stopifnot(identical(features, rownames(object)))
    layer_equal <- vapply(seq_along(current), function(i) {
      isTRUE(all.equal(reference[[i]], current[[i]], tolerance = 0))
    }, logical(1))
    stopifnot(all(layer_equal), identical(reference_umap, umap[shared, , drop = FALSE]))
    pca_equal <- identical(reference_pca, pca[shared, , drop = FALSE])
    pca_compatible <- pca_compatible && pca_equal
    duplicate <- intersect(seen, colnames(object))
    stopifnot(!length(duplicate) || setequal(duplicate, shared))
    keep <- setdiff(colnames(object), seen)
    checks[[id]] <- data.frame(object = id, original_cells = ncol(object),
                               retained_cells = length(keep), duplicate_cells_removed = length(duplicate),
                               shared_counts_equal = layer_equal[1], shared_data_equal = layer_equal[2],
                               shared_umap_equal = TRUE, shared_pca_equal = pca_equal)
    objects[[id]] <- subset(object, cells = keep)
    seen <- c(seen, keep)
    message(id, ": retained ", length(keep), "; removed duplicates ", length(duplicate))
    rm(object, current)
    invisible(gc())
  }
  stopifnot(identical(seen, all_cells), !anyDuplicated(seen))
  reduction_names <- c(CONFIG$object_overview$reduction, if (pca_compatible) "pca")
  source_embeddings <- lapply(reduction_names, function(reduction) {
    do.call(rbind, lapply(objects, function(x) Embeddings(x[[reduction]])))[all_cells, , drop = FALSE]
  })
  names(source_embeddings) <- c("umap.source", if (pca_compatible) "pca.source")
  # These are source coordinates, not evidence of a newly integrated embedding.
  merged <- merge(objects[[1]], y = objects[-1], merge.data = TRUE,
                  merge.dr = FALSE, project = CONFIG$project_id)
  if (inherits(merged[["RNA"]], "Assay5")) merged <- JoinLayers(merged, assay = "RNA")
  stopifnot(setequal(Layers(merged[["RNA"]]), c("counts", "data")))
  merged <- merged[, all_cells]
  membership_matrix <- vapply(source_cells, function(cells) all_cells %in% cells,
                               logical(length(all_cells)))
  memberships <- apply(membership_matrix, 1, function(member) paste(names(source_cells)[member], collapse = ";"))
  merged$source_objects <- unname(memberships)
  merged$n_source_objects <- lengths(strsplit(memberships, ";", fixed = TRUE))
  Idents(merged) <- "cell_type"
  for (name in names(source_embeddings)) {
    embedding <- source_embeddings[[name]]
    key <- if (name == "umap.source") "UMAPSRC_" else "PCASRC_"
    colnames(embedding) <- paste0(key, seq_len(ncol(embedding)))
    merged[[name]] <- CreateDimReducObject(embeddings = embedding, key = key, assay = "RNA")
  }
  merged@misc$merge_provenance <- list(
    source_paths = CONFIG$input_objects, checks = do.call(rbind, checks),
    shared_cells = length(shared), original_active_identities =
      lapply(cache, function(x) unique(x$metadata$active_ident)),
    note = "Original expression and metadata retained; duplicate cells verified and kept once. Source reductions are not a newly computed joint embedding. Leiden labels retain source meaning."
  )
  expected_metadata <- do.call(rbind, lapply(objects, function(x) x[[]]))
  rownames(expected_metadata) <- all_cells
  stopifnot(ncol(merged) == length(all_cells), nrow(merged) == length(features),
            identical(colnames(merged), all_cells), !anyDuplicated(colnames(merged)),
            all(vapply(names(expected_metadata), function(column) {
              identical(as.character(merged[[]][[column]]), as.character(expected_metadata[[column]]))
            }, logical(1))))
  expected_sums <- lapply(c("counts", "data"), function(layer) {
    Matrix::colSums(LayerData(merged, assay = "RNA", layer = layer))
  })
  expected_meta <- merged[[]]
  expected_umap <- Embeddings(merged[["umap.source"]])
  rm(objects, reference, cache, source_embeddings, expected_metadata)
  invisible(gc())
  message("Writing qs: ", nrow(merged), " genes x ", ncol(merged), " cells")
  temporary <- paste0(output, ".partial")
  qs::qsave(merged, temporary, preset = "high", nthreads = 4)
  rm(merged)
  invisible(gc())
  checked <- qs::qread(temporary, nthreads = 4)
  stopifnot(inherits(checked, "Seurat"), ncol(checked) == length(all_cells),
            identical(rownames(checked), features), identical(colnames(checked), all_cells),
            identical(checked[[]], expected_meta),
            identical(Embeddings(checked[["umap.source"]]), expected_umap))
  for (i in seq_along(expected_sums)) {
    stopifnot(identical(Matrix::colSums(LayerData(checked, assay = "RNA", layer = c("counts", "data")[i])),
                        expected_sums[[i]]))
  }
  stopifnot(file.rename(temporary, output))
  write.table(do.call(rbind, checks),
              file.path(PROJECT_DIR, CONFIG$object_overview$output_dir, "merge_checks.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "03_merge_sessionInfo.txt"))
  message("MERGE_QS_CHECK_OK: ", nrow(checked), " genes; ", ncol(checked),
          " unique cells; shared=", length(shared), "; file=", output)
  invisible(output)
}
