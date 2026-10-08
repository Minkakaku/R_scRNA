# 独立会话以Seurat原生DotPlot重算正式点图，并核对新标签、中文和输入完整性。
step_validate_literature_annotation <- function() {
  suppressPackageStartupMessages(library(Seurat))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- CONFIG$literature_annotation; folder <- file.path(PROJECT_DIR, s$evidence_dir)
  manifest <- read.delim(file.path(folder, 'subgroup_annotation_summary.tsv'))
  mapping <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8')
  merged <- readRDS(file.path(folder, 'current_cell_annotations.rds'))
  validations <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; sc <- qs::qread(file.path(PROJECT_DIR, manifest$checkpoint[i]), nthreads = 4)
    m <- mapping[mapping$subgroup == group, ]; idx <- match(as.character(sc$subgroup_cluster), m$subgroup_cluster)
    stopifnot(identical(as.character(sc$subtype_annotation), m$annotation[idx]),
      identical(sc[[]]$annotation_limitations, m$limitations_CN[idx]),
      !any(grepl('<U\\+', sc$annotation_limitations)),
      identical(as.character(Idents(sc)), as.character(sc$subtype_annotation)))
    ix <- match(colnames(sc), merged$cell_id)
    stopifnot(!anyNA(ix), identical(as.character(sc$subgroup_cluster), merged$current_subgroup_cluster[ix]),
      identical(as.character(sc$subtype_annotation), as.character(merged$subtype_annotation[ix])))
    source_path <- file.path(PROJECT_DIR, sc@misc$literature_annotation$input)
    stopifnot(identical(digest::digest(file = source_path, algo = 'sha256'), manifest$input_sha256[i]))
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    unchanged <- setdiff(names(e$metadata), 'annotation_stage')
    stopifnot(identical(sc[[]][, unchanged], e$metadata[, unchanged]), identical(Embeddings(sc, 'umap'), e$umap))
    files <- file.path(folder, group, 'DotPlot_annotation_data.tsv')
    if (group == s$focus_subgroup) files <- c(files, file.path(folder, 'neutrophil_functional_DotPlot_data.tsv'))
    Idents(sc) <- 'subgroup_cluster'
    for (file in files) {
      d <- read.delim(file)
      actual <- suppressWarnings(DotPlot(sc, assay = 'RNA', features = unique(d$gene)))$data
      ii <- match(paste(d$cluster, d$gene), paste(actual$id, actual$features.plot))
      stopifnot(!anyNA(ii), max(abs(d$percent_expressed - actual$pct.exp[ii])) < 1e-8,
        max(abs(d$average_expression - actual$avg.exp[ii])) < 1e-8,
        max(abs(d$color_value - actual$avg.exp.scaled[ii])) < 1e-8)
    }
    validations[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nrow(m),
      annotation_matches_config = TRUE, Chinese_preserved = TRUE, DotPlot_matches_Seurat = TRUE,
      all_prior_metadata_retained = TRUE, source_sha256_unchanged = TRUE, UMAP_unchanged = TRUE)
    log_msg('INDEPENDENT_ANNOTATION_VALIDATION_OK ', group)
    rm(sc, e, actual); invisible(gc())
  }
  stopifnot(nrow(merged) == sum(manifest$cells), !anyDuplicated(merged$cell_id))
  write.table(do.call(rbind, validations), file.path(folder, 'independent_annotation_validation.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  log_msg('ALL_INDEPENDENT_ANNOTATION_VALIDATIONS_OK')
}
