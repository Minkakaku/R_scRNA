# 在当前cluster上按文献类型/功能状态命名；未匹配的状态保留谱系_基因名。
step_review_functional_markers <- function(settings = CONFIG$functional_annotation) {
  s <- settings
  mapping <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8')
  panel <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  output <- file.path(PROJECT_DIR, s$evidence_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  values <- list()
  for (group in unique(mapping$subgroup)) {
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    m <- mapping[mapping$subgroup == group, ]
    genes <- unique(c(unlist(strsplit(m$support_markers, ';', fixed = TRUE)),
      strsplit(panel$genes[panel$subgroup == group], ';', fixed = TRUE)[[1]]))
    stopifnot(all(genes %in% rownames(e$average)))
    d <- expand.grid(gene = genes, cluster = colnames(e$average), stringsAsFactors = FALSE)
    d$subgroup <- group
    d$average <- as.vector(e$average[genes, , drop = FALSE])
    d$percent <- as.vector(e$percent[genes, , drop = FALSE])
    values[[group]] <- d
    log_msg('FUNCTIONAL_MARKER_REVIEW_OK ', group, '; genes=', length(genes))
  }
  write.table(do.call(rbind, values), file.path(output, 'marker_expression.tsv'),
    sep = '\t', quote = FALSE, row.names = FALSE)
  # 中性粒文献程序沿用已验证的数值；仅重命名图内文字。
  for (name in c('neutrophil_cluster.tsv', 'neutrophil_cluster_sample.tsv')) {
    input <- file.path(PROJECT_DIR, s$input_annotation_evidence, name)
    destination <- file.path(output, name)
    stopifnot(!file.exists(destination), file.copy(input, destination),
      identical(digest::digest(file = input, algo = 'sha256'), digest::digest(file = destination, algo = 'sha256')))
  }
}

step_apply_functional_annotations <- function(settings = CONFIG$functional_annotation) {
  suppressPackageStartupMessages(library(Seurat))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- settings; output <- file.path(PROJECT_DIR, s$evidence_dir)
  mapping <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8')
  prior <- read.delim(file.path(PROJECT_DIR, s$input_annotation_evidence, 'cluster_annotations.tsv'), encoding = 'UTF-8')
  manifest <- read.delim(file.path(PROJECT_DIR, s$input_manifest))
  panel <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  metadata <- readRDS(file.path(PROJECT_DIR, s$input_annotation_evidence, 'current_cell_annotations.rds'))
  stopifnot(nrow(mapping) == sum(manifest$clusters), !anyDuplicated(paste(mapping$subgroup, mapping$subgroup_cluster)))
  dir.create(file.path(PROJECT_DIR, s$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  replaced <- c('subtype_annotation', 'biological_lineage', 'annotation_confidence', 'annotation_PMIDs',
    'annotation_limitations', 'literature_cluster_label', 'annotation_stage',
    'functional_annotation_stage', 'annotation_name_type', 'functional_interpretation', 'neutrophil_function')
  summaries <- mappings <- checks <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; log_msg('FUNCTIONAL_ANNOTATE ', group)
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    source_hash <- digest::digest(file = input, algo = 'sha256')
    stopifnot(identical(source_hash, manifest$output_sha256[i]))
    sc <- qs::qread(input, nthreads = 4); old <- sc[[]]; fingerprints <- annotation_fingerprints(sc)
    m <- mapping[mapping$subgroup == group, ]; m <- m[order(m$subgroup_cluster), ]
    pm <- prior[prior$subgroup == group, ]; j <- match(m$subgroup_cluster, pm$subgroup_cluster)
    for (column in c('cells', 'samples_present', 'max_sample_fraction')) m[[column]] <- pm[[column]][j]
    stopifnot(identical(m$previous_annotation, pm$annotation[j]))
    idx <- match(as.character(sc$subgroup_cluster), m$subgroup_cluster)
    stopifnot(!anyNA(idx), identical(as.character(sc$subtype_annotation), m$previous_annotation[idx]),
      identical(as.character(sc$neutrophil_function), pm$neutrophil_function[j][idx]))
    for (column in intersect(replaced, names(old))) sc[[paste0(s$historical_prefix, column)]] <- old[[column]]
    sc$subtype_annotation <- factor(m$annotation[idx], levels = unique(m$annotation))
    sc$biological_lineage <- m$biological_lineage[idx]
    sc$annotation_confidence <- m$confidence[idx]
    sc$annotation_PMIDs <- m$PMIDs[idx]
    sc$annotation_limitations <- m$limitations_CN[idx]
    sc$literature_cluster_label <- m$annotation[idx]
    sc$annotation_stage <- s$stage
    sc$functional_annotation_stage <- s$stage
    sc$annotation_name_type <- m$name_type[idx]
    sc$functional_interpretation <- m$functional_interpretation_CN[idx]
    # 功能方向随经审核的新映射更新，旧判断按本轮历史前缀完整保留。
    sc$neutrophil_function <- m$neutrophil_function[idx]
    Idents(sc) <- 'subtype_annotation'
    if (!is.null(sc@misc$functional_annotation)) sc@misc[[paste0(s$historical_prefix, 'functional_annotation')]] <- sc@misc$functional_annotation
    sc@misc$functional_annotation <- list(mapping = m, settings = s, input = manifest$checkpoint[i],
      input_sha256 = source_hash, original_clusters_and_RNA_preserved = TRUE,
      positive_labels_indicate_RNA_enrichment_not_protein_gating = TRUE)
    stopifnot(identical(sc[[]][, setdiff(names(old), replaced)], old[, setdiff(names(old), replaced)]))
    for (column in intersect(replaced, names(old))) stopifnot(identical(sc[[]][[paste0(s$historical_prefix, column)]], old[[column]]))
    destination <- file.path(s$checkpoint_dir, paste0(s$object_prefix, '.', CONFIG$project_id, '_', group, '_functional_annotated.qs'))
    save_annotation_verified(sc, file.path(PROJECT_DIR, destination), fingerprints)
    group_dir <- file.path(output, group); dir.create(group_dir, showWarnings = FALSE)
    plot_object <- CreateSeuratObject(LayerData(sc, assay = 'RNA', layer = 'counts')[1, , drop = FALSE], meta.data = sc[[]])
    plot_object[['umap']] <- sc[['umap']]
    saveRDS(list(plot_object = plot_object, mapping = m), file.path(group_dir, 'annotation_plot_data.rds'))
    genes <- strsplit(panel$genes[panel$subgroup == group], ';', fixed = TRUE)[[1]]
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    direct <- suppressWarnings(DotPlot(sc, assay = 'RNA', features = genes, group.by = 'subgroup_cluster'))$data
    ii <- cbind(match(as.character(direct$features.plot), rownames(e$average)), match(as.character(direct$id), colnames(e$average)))
    stopifnot(max(abs(direct$pct.exp - e$percent[ii])) < 1e-8,
      max(abs(direct$avg.exp - e$average[ii])) < 1e-8, identical(Embeddings(sc, 'umap'), e$umap))
    write.table(direct, file.path(group_dir, 'Seurat_DotPlot_reference.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    # 样本分别记录实际表达；不把细胞当作处理组统计重复。
    expr <- LayerData(sc, assay = 'RNA', layer = 'data')[genes, , drop = FALSE]
    keys <- paste(sc$subgroup_cluster, sc$sample, sep = '|')
    sample_values <- do.call(rbind, lapply(unique(keys), function(key) {
      ix <- which(keys == key); x <- expr[, ix, drop = FALSE]; linear <- x; linear@x <- expm1(linear@x)
      data.frame(cluster = as.character(sc$subgroup_cluster[ix[1]]), sample = unname(sc$sample[ix[1]]), cells = length(ix),
        gene = genes, percent = as.numeric(Matrix::rowMeans(x > 0) * 100), average = as.numeric(Matrix::rowMeans(linear)))
    }))
    write.table(sample_values, file.path(group_dir, 'marker_by_cluster_sample.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    ix <- match(colnames(sc), metadata$cell_id); stopifnot(!anyNA(ix))
    for (column in replaced) {
      if (!column %in% names(metadata)) metadata[[column]] <- NA_character_
      metadata[[column]] <- as.character(metadata[[column]])
      metadata[ix, column] <- as.character(sc[[]][[column]])
    }
    stopifnot(identical(source_hash, digest::digest(file = input, algo = 'sha256')))
    summaries[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nrow(m), labels = length(unique(m$annotation)),
      checkpoint = destination, input_sha256 = source_hash, output_sha256 = digest::digest(file = file.path(PROJECT_DIR, destination), algo = 'sha256'),
      RNA_reductions_graphs_preserved = TRUE, readback_verified = TRUE)
    mappings[[group]] <- m
    checks[[group]] <- data.frame(subgroup = group, prior_metadata_preserved = TRUE,
      source_QS_unchanged = TRUE, marker_values_match_Seurat = TRUE, UMAP_unchanged = TRUE)
    log_msg('FUNCTIONAL_QS_OK ', group, '; cells=', ncol(sc))
    rm(sc, old, e, expr, plot_object, direct, sample_values); invisible(gc())
  }
  saveRDS(metadata, file.path(output, 'current_cell_annotations.rds'))
  for (pair in list(c('summaries', 'subgroup_annotation_summary.tsv'), c('mappings', 'cluster_annotations.tsv'), c('checks', 'annotation_validation.tsv')))
    write.table(do.call(rbind, get(pair[1])), file.path(output, pair[2]), sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  log_msg('ALL_FUNCTIONAL_QS_OK')
}

step_merge_functional_annotations <- function(settings = CONFIG$functional_annotation) {
  suppressPackageStartupMessages(library(Seurat))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- settings; output <- file.path(PROJECT_DIR, s$evidence_dir)
  new <- readRDS(file.path(output, 'current_cell_annotations.rds'))
  input <- file.path(PROJECT_DIR, s$full_input); source_hash <- digest::digest(file = input, algo = 'sha256')
  sc <- qs::qread(input, nthreads = 4); old <- sc[[]]; fingerprints <- annotation_fingerprints(sc)
  ix <- match(colnames(sc), new$cell_id); stopifnot(!anyNA(ix), nrow(new) == ncol(sc), !anyDuplicated(new$cell_id))
  columns <- c('subtype_annotation', 'biological_lineage', 'annotation_confidence', 'annotation_PMIDs',
    'annotation_limitations', 'literature_cluster_label', 'annotation_stage',
    'functional_annotation_stage', 'annotation_name_type', 'functional_interpretation', 'neutrophil_function')
  for (column in columns) {
    if (column %in% names(old)) sc[[paste0(s$historical_prefix, column)]] <- old[[column]]
    sc[[column]] <- as.character(new[[column]][ix])
  }
  Idents(sc) <- 'subtype_annotation'
  if (!is.null(sc@misc$functional_annotation)) sc@misc[[paste0(s$historical_prefix, 'functional_annotation')]] <- sc@misc$functional_annotation
  sc@misc$functional_annotation <- list(settings = s, input_sha256 = source_hash,
    current_cluster_column = 'current_subgroup_cluster', global_reductions_are_historical = TRUE, metadata_join = 'cell_ID')
  stopifnot(identical(sc[[]][, setdiff(names(old), columns)], old[, setdiff(names(old), columns)]))
  for (column in intersect(columns, names(old))) stopifnot(identical(sc[[]][[paste0(s$historical_prefix, column)]], old[[column]]))
  save_annotation_verified(sc, file.path(PROJECT_DIR, s$full_output), fingerprints)
  stopifnot(identical(digest::digest(file = input, algo = 'sha256'), source_hash))
  write.table(data.frame(cells = ncol(sc), genes = nrow(sc), samples = length(unique(sc$sample)), checkpoint = s$full_output,
    input_sha256 = source_hash, output_sha256 = digest::digest(file = file.path(PROJECT_DIR, s$full_output), algo = 'sha256'),
    annotations_matched_by_cell_ID = TRUE, RNA_reductions_graphs_preserved = TRUE, readback_verified = TRUE),
    file.path(output, 'complete_object_validation.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  log_msg('COMPLETE_FUNCTIONAL_QS_OK')
}
