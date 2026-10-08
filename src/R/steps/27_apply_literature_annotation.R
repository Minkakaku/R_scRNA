# 在既有cluster上写入有来源的人工marker注释；RNA、图结构、降维和历史标签均保留。
annotation_fingerprints <- function(x) c(RNA = digest::digest(x[['RNA']], algo = 'xxhash64'),
  reductions = digest::digest(x@reductions, algo = 'xxhash64'), graphs = digest::digest(x@graphs, algo = 'xxhash64'))

save_annotation_verified <- function(sc, destination, expected_fingerprints) {
  stopifnot(!file.exists(destination), !file.exists(paste0(destination, '.partial')),
    identical(annotation_fingerprints(sc), expected_fingerprints))
  meta <- sc[[]]; cells <- colnames(sc); genes <- rownames(sc); ids <- Idents(sc)
  qs::qsave(sc, paste0(destination, '.partial'), preset = 'high', nthreads = 4)
  check <- qs::qread(paste0(destination, '.partial'), nthreads = 4)
  stopifnot(isTRUE(validObject(check, test = TRUE)), identical(check[[]], meta), identical(Idents(check), ids),
    identical(colnames(check), cells), identical(rownames(check), genes),
    identical(annotation_fingerprints(check), expected_fingerprints),
    file.rename(paste0(destination, '.partial'), destination))
  invisible(TRUE)
}

step_apply_literature_annotation <- function() {
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  suppressPackageStartupMessages(library(Seurat))
  s <- CONFIG$literature_annotation
  output <- file.path(PROJECT_DIR, s$evidence_dir)
  dir.create(file.path(PROJECT_DIR, s$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  manifest <- read.delim(file.path(PROJECT_DIR, s$input_manifest))
  mapping <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8', check.names = FALSE)
  panel <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  stopifnot(nrow(mapping) == sum(manifest$clusters), !anyDuplicated(paste(mapping$subgroup, mapping$subgroup_cluster)))
  all_meta <- summaries <- all_mappings <- all_support <- all_previous <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; log_msg('ANNOTATE ', group)
    folder <- file.path(output, group); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i]); input_hash <- digest::digest(file = input, algo = 'sha256')
    sc <- qs::qread(input, nthreads = 4); before <- annotation_fingerprints(sc); old <- sc[[]]
    m <- mapping[mapping$subgroup == group, ]; m <- m[order(m$subgroup_cluster), ]
    ids <- as.character(sc$subgroup_cluster); idx <- match(ids, m$subgroup_cluster)
    stopifnot(!anyNA(idx), setequal(ids, as.character(m$subgroup_cluster)), ncol(sc) == manifest$cells[i])
    m$cells <- as.integer(table(factor(ids, levels = m$subgroup_cluster)))
    counts <- as.data.frame(table(cluster = ids, sample = sc$sample))
    m$samples_present <- vapply(as.character(m$subgroup_cluster), function(k) sum(counts$Freq[counts$cluster == k] > 0), integer(1))
    m$max_sample_fraction <- vapply(as.character(m$subgroup_cluster), function(k) max(counts$Freq[counts$cluster == k]) / sum(counts$Freq[counts$cluster == k]), numeric(1))
    new <- data.frame(row.names = colnames(sc), subtype_annotation = m$annotation[idx],
      biological_lineage = m$biological_lineage[idx], annotation_confidence = m$confidence[idx],
      neutrophil_function = m$neutrophil_function[idx], annotation_PMIDs = m$PMIDs[idx],
      annotation_limitations = m$limitations_CN[idx],
      literature_cluster_label = paste0('C', ids, ': ', m$annotation[idx]),
      literature_annotation_stage = s$stage)
    stopifnot(!any(names(new) %in% names(old)))
    sc <- AddMetaData(sc, new)
    # annotation_stage的旧值另存；其他旧列逐列保持一致。
    sc$pre_literature_annotation_stage <- sc$annotation_stage
    sc$annotation_stage <- s$stage
    sc$subtype_annotation <- factor(sc$subtype_annotation, levels = unique(m$annotation))
    sc$literature_cluster_label <- factor(sc$literature_cluster_label, levels = paste0('C', m$subgroup_cluster, ': ', m$annotation))
    Idents(sc) <- 'subtype_annotation'
    sc@misc$literature_annotation <- list(mapping = m, settings = s, input = manifest$checkpoint[i],
      input_sha256 = input_hash, mapping_sha256 = digest::digest(file = file.path(PROJECT_DIR, s$mapping), algo = 'sha256'),
      original_cluster_and_RNA_preserved = TRUE, functional_labels_are_inferences = TRUE)
    stopifnot(identical(sc[[]][, setdiff(names(old), 'annotation_stage')], old[, setdiff(names(old), 'annotation_stage')]),
      identical(sc[[]]$pre_literature_annotation_stage, old$annotation_stage))
    checkpoint <- file.path(s$checkpoint_dir, paste0(s$object_prefix, '.', CONFIG$project_id, '_', group, '_literature_annotated.qs'))
    save_annotation_verified(sc, file.path(PROJECT_DIR, checkpoint), before)
    stopifnot(identical(digest::digest(file = input, algo = 'sha256'), input_hash))
    # 图用轻量对象只保存一个真实RNA基因、原元数据及原UMAP，不产生替代数据坐标。
    plot_object <- CreateSeuratObject(LayerData(sc, assay = 'RNA', layer = 'counts')[1, , drop = FALSE], meta.data = sc[[]])
    plot_object[['umap']] <- sc[['umap']]
    saveRDS(list(plot_object = plot_object, mapping = m), file.path(folder, 'annotation_plot_data.rds'))
    genes <- unique(c(unlist(strsplit(m$support_markers, ';', fixed = TRUE)),
      strsplit(panel$genes[panel$subgroup == group], ';', fixed = TRUE)[[1]]))
    stopifnot(all(genes %in% rownames(sc)))
    expr <- LayerData(sc, assay = 'RNA', layer = 'data')[genes, , drop = FALSE]
    keys <- paste(ids, sc$sample, sep = '|')
    sample_evidence <- do.call(rbind, lapply(unique(keys), function(k) {
      ix <- which(keys == k); x <- expr[, ix, drop = FALSE]; linear <- x; linear@x <- expm1(linear@x)
      data.frame(cluster = ids[ix[1]], sample = sc$sample[ix[1]], cells = length(ix), gene = genes,
        percent = Matrix::rowMeans(x > 0) * 100, average = Matrix::rowMeans(linear))
    }))
    write.table(sample_evidence, file.path(folder, 'marker_by_cluster_sample.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    stopifnot(identical(e$umap, Embeddings(sc, 'umap')), identical(e$metadata, old))
    support <- do.call(rbind, lapply(seq_len(nrow(m)), function(j) {
      g <- strsplit(m$support_markers[j], ';', fixed = TRUE)[[1]]; k <- as.character(m$subgroup_cluster[j])
      data.frame(subgroup = group, cluster = k, annotation = m$annotation[j], gene = g,
        percent = e$percent[g, k], average = e$average[g, k])
    }))
    previous <- do.call(rbind, lapply(c('cell_type', 'pre_reanalysis_candidate_annotation'), function(column) {
      d <- as.data.frame(table(cluster = ids, old_label = old[[column]])); names(d)[3] <- 'cells'
      d <- d[d$cells > 0, ]; d$source_column <- column; d$subgroup <- group
      d$fraction <- d$cells / ave(d$cells, d$cluster, FUN = sum); d
    }))
    cell_meta <- data.frame(cell_id = colnames(sc), analysis_container = group,
      current_subgroup_cluster = ids, sc[[]][, names(new)],
      subgroup_UMAP_1 = Embeddings(sc, 'umap')[,1], subgroup_UMAP_2 = Embeddings(sc, 'umap')[,2])
    all_meta[[group]] <- cell_meta; all_mappings[[group]] <- m; all_support[[group]] <- support; all_previous[[group]] <- previous
    summaries[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nrow(m),
      labels = length(unique(m$annotation)), checkpoint = checkpoint, input_sha256 = input_hash,
      output_sha256 = digest::digest(file = file.path(PROJECT_DIR, checkpoint), algo = 'sha256'),
      RNA_reductions_graphs_preserved = TRUE, readback_verified = TRUE)
    write.table(do.call(rbind, summaries), file.path(output, 'subgroup_annotation_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    log_msg('LITERATURE_ANNOTATION_QS_OK ', group, '; cells=', ncol(sc), '; clusters=', nrow(m))
    rm(sc, plot_object, e, expr, old, sample_evidence); invisible(gc())
  }
  metadata <- do.call(rbind, all_meta); rownames(metadata) <- metadata$cell_id
  stopifnot(!anyDuplicated(metadata$cell_id), nrow(metadata) == sum(manifest$cells))
  saveRDS(metadata, file.path(output, 'current_cell_annotations.rds'))
  for (pair in list(c('all_mappings', 'cluster_annotations.tsv'), c('all_support', 'support_marker_expression.tsv'), c('all_previous', 'previous_label_composition.tsv')))
    write.table(do.call(rbind, get(pair[1])), file.path(output, pair[2]), sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, '31_literature_annotation_sessionInfo.txt'))
  log_msg('ALL_LITERATURE_ANNOTATION_QS_OK cells=', nrow(metadata))
}

step_merge_literature_annotations <- function() {
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  suppressPackageStartupMessages(library(Seurat))
  s <- CONFIG$literature_annotation; output <- file.path(PROJECT_DIR, s$evidence_dir)
  metadata <- readRDS(file.path(output, 'current_cell_annotations.rds'))
  input <- file.path(PROJECT_DIR, s$full_input); source_hash <- digest::digest(file = input, algo = 'sha256')
  sc <- qs::qread(input, nthreads = 4); before <- annotation_fingerprints(sc); old <- sc[[]]
  ix <- match(colnames(sc), metadata$cell_id)
  stopifnot(!anyNA(ix), ncol(sc) == nrow(metadata), !anyDuplicated(ix))
  new <- metadata[ix, setdiff(names(metadata), 'cell_id'), drop = FALSE]
  new[] <- lapply(new, function(x) if (is.factor(x)) as.character(x) else x)
  stopifnot(!any(names(new) %in% names(old)))
  sc <- AddMetaData(sc, new)
  sc$pre_literature_annotation_stage <- sc$annotation_stage
  sc$annotation_stage <- s$stage
  Idents(sc) <- 'subtype_annotation'
  sc@misc$literature_annotation <- list(settings = s, full_source_sha256 = source_hash,
    subgroup_manifest = read.delim(file.path(output, 'subgroup_annotation_summary.tsv')),
    current_cluster_column = 'current_subgroup_cluster', global_reductions_are_historical = TRUE,
    local_coordinates = 'subgroup_UMAP_1/2 valid only within analysis_container; not a common UMAP')
  stopifnot(identical(sc[[]][, setdiff(names(old), 'annotation_stage')], old[, setdiff(names(old), 'annotation_stage')]))
  save_annotation_verified(sc, file.path(PROJECT_DIR, s$full_output), before)
  stopifnot(identical(digest::digest(file = input, algo = 'sha256'), source_hash))
  write.table(data.frame(cells = ncol(sc), genes = nrow(sc), samples = length(unique(sc$sample)),
    checkpoint = s$full_output, input_sha256 = source_hash,
    output_sha256 = digest::digest(file = file.path(PROJECT_DIR, s$full_output), algo = 'sha256'),
    annotations_matched_by_cell_ID = TRUE, RNA_reductions_graphs_preserved = TRUE, readback_verified = TRUE),
    file.path(output, 'complete_object_validation.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  log_msg('COMPLETE_LITERATURE_ANNOTATED_QS_OK cells=', ncol(sc))
}
