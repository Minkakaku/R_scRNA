# 旧亚群名由当前marker支持后回填；明确保留上一版注释，排除中性粒更新。
step_update_legacy_marker_annotations <- function() {
  suppressPackageStartupMessages(library(Seurat))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- CONFIG$legacy_marker_annotation; folder <- file.path(PROJECT_DIR, s$evidence_dir)
  manifest <- read.delim(file.path(PROJECT_DIR, s$input_manifest))
  mapping <- read.delim(file.path(PROJECT_DIR, s$mapping), encoding = 'UTF-8')
  previous <- read.delim(file.path(PROJECT_DIR, s$input_annotation_evidence, 'cluster_annotations.tsv'), encoding = 'UTF-8')
  catalog <- read.delim(file.path(PROJECT_DIR, s$legacy_markers), encoding = 'UTF-8')
  panels <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  dir.create(file.path(PROJECT_DIR, s$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  excluded <- unlist(s$excluded_subgroups)
  keep <- mapping$subgroup %in% excluded
  stopifnot(identical(mapping[keep, ], previous[previous$subgroup %in% excluded, names(mapping)]))
  replaced <- c('subtype_annotation', 'biological_lineage', 'annotation_confidence', 'annotation_PMIDs',
    'annotation_limitations', 'literature_cluster_label', 'annotation_stage')
  metadata <- readRDS(file.path(PROJECT_DIR, s$input_annotation_evidence, 'current_cell_annotations.rds'))
  metadata$annotation_name_source <- 'unchanged_literature'
  metadata$annotation_stage <- metadata$literature_annotation_stage
  metadata$legacy_marker_annotation_stage <- NA_character_
  updated_manifest <- all_mappings <- checks <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    m <- mapping[mapping$subgroup == group, ]; m <- m[order(m$subgroup_cluster), ]
    pm <- previous[previous$subgroup == group, ]; j <- match(m$subgroup_cluster, pm$subgroup_cluster)
    m$cells <- pm$cells[j]; m$samples_present <- pm$samples_present[j]; m$max_sample_fraction <- pm$max_sample_fraction[j]
    m$previous_annotation <- pm$annotation[j]
    m$name_source <- if (group %in% excluded) 'unchanged_literature' else
      catalog$source[match(paste(group, m$annotation), paste(catalog$subgroup, catalog$old_annotation))]
    m$name_source[is.na(m$name_source)] <- 'marker_supported_exception'
    all_mappings[[group]] <- m
    if (group %in% excluded) {
      updated_manifest[[group]] <- data.frame(subgroup = group, cells = manifest$cells[i], clusters = manifest$clusters[i],
        labels = manifest$labels[i], checkpoint = manifest$checkpoint[i], input_sha256 = manifest$output_sha256[i],
        output_sha256 = manifest$output_sha256[i], changed_cluster_names = 0, status = 'retained_unchanged')
      next
    }
    log_msg('LEGACY_MARKER_ANNOTATE ', group)
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i]); input_hash <- digest::digest(file = input, algo = 'sha256')
    stopifnot(identical(input_hash, manifest$output_sha256[i]))
    sc <- qs::qread(input, nthreads = 4); before <- annotation_fingerprints(sc); old <- sc[[]]
    idx <- match(as.character(sc$subgroup_cluster), m$subgroup_cluster)
    stopifnot(!anyNA(idx))
    for (column in replaced) {
      stopifnot(!paste0('pre_legacy_marker_', column) %in% names(old))
      sc[[paste0('pre_legacy_marker_', column)]] <- old[[column]]
    }
    sc$subtype_annotation <- factor(m$annotation[idx], levels = unique(m$annotation))
    sc$biological_lineage <- m$biological_lineage[idx]
    sc$annotation_confidence <- m$confidence[idx]
    sc$annotation_PMIDs <- m$PMIDs[idx]
    sc$annotation_limitations <- m$limitations_CN[idx]
    sc$literature_cluster_label <- paste0('C', sc$subgroup_cluster, ': ', m$annotation[idx])
    sc$annotation_stage <- s$stage
    sc$legacy_marker_annotation_stage <- s$stage
    sc$annotation_name_source <- m$name_source[idx]
    Idents(sc) <- 'subtype_annotation'
    sc@misc$legacy_marker_annotation <- list(mapping = m, settings = s, input = manifest$checkpoint[i], input_sha256 = input_hash,
      method = 'manual matching of legacy name marker programs to current cluster expression; not old-label majority',
      prior_annotation_prefix = 'pre_legacy_marker_')
    stopifnot(identical(sc[[]][, setdiff(names(old), replaced)], old[, setdiff(names(old), replaced)]))
    for (column in replaced) stopifnot(identical(sc[[]][[paste0('pre_legacy_marker_', column)]], old[[column]]))
    outpath <- file.path(s$checkpoint_dir, paste0(s$object_prefix, '.', CONFIG$project_id, '_', group, '_legacy_marker_annotated.qs'))
    save_annotation_verified(sc, file.path(PROJECT_DIR, outpath), before)
    group_dir <- file.path(folder, group); dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    plot_object <- CreateSeuratObject(LayerData(sc, assay = 'RNA', layer = 'counts')[1, , drop = FALSE], meta.data = sc[[]])
    plot_object[['umap']] <- sc[['umap']]
    saveRDS(list(plot_object = plot_object, mapping = m), file.path(group_dir, 'annotation_plot_data.rds'))
    genes <- strsplit(panels$genes[panels$subgroup == group], ';', fixed = TRUE)[[1]]
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    direct <- suppressWarnings(DotPlot(sc, assay = 'RNA', features = genes, group.by = 'subgroup_cluster'))$data
    stopifnot(max(abs(direct$pct.exp - e$percent[cbind(match(as.character(direct$features.plot), rownames(e$percent)),
      match(as.character(direct$id), colnames(e$percent)))])) < 1e-8,
      identical(Embeddings(sc, 'umap'), e$umap))
    write.table(direct, file.path(group_dir, 'Seurat_DotPlot_reference.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    # 各样本分别检查命名marker；保存为计算证据，不做逐细胞组间显著性推断。
    expr <- LayerData(sc, assay = 'RNA', layer = 'data')[genes, , drop = FALSE]
    ids <- as.character(sc$subgroup_cluster); keys <- paste(ids, sc$sample, sep = '|')
    samples <- do.call(rbind, lapply(unique(keys), function(key) {
      ix <- which(keys == key); x <- expr[, ix, drop = FALSE]; linear <- x; linear@x <- expm1(linear@x)
      data.frame(cluster = ids[ix[1]], sample = sc$sample[ix[1]], cells = length(ix), gene = genes,
        percent = as.numeric(Matrix::rowMeans(x > 0) * 100), average = as.numeric(Matrix::rowMeans(linear)))
    }))
    write.table(samples, file.path(group_dir, 'marker_by_cluster_sample.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    ix <- match(colnames(sc), metadata$cell_id)
    stopifnot(!anyNA(ix))
    for (column in c(replaced, 'legacy_marker_annotation_stage', 'annotation_name_source')) {
      metadata[[column]] <- as.character(metadata[[column]])
      metadata[ix, column] <- as.character(sc[[]][[column]])
    }
    stopifnot(identical(input_hash, digest::digest(file = input, algo = 'sha256')))
    updated_manifest[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nrow(m), labels = length(unique(m$annotation)),
      checkpoint = outpath, input_sha256 = input_hash, output_sha256 = digest::digest(file = file.path(PROJECT_DIR, outpath), algo = 'sha256'),
      changed_cluster_names = sum(m$annotation != m$previous_annotation), status = 'updated_readback_verified')
    checks[[group]] <- data.frame(subgroup = group, RNA_reductions_graphs_preserved = TRUE, prior_metadata_preserved = TRUE,
      current_marker_values_match_Seurat = TRUE, source_QS_unchanged = TRUE)
    log_msg('LEGACY_MARKER_QS_OK ', group, '; renamed_clusters=', sum(m$annotation != m$previous_annotation))
    rm(sc, e, expr, plot_object, old, direct, samples); invisible(gc())
  }
  saveRDS(metadata, file.path(folder, 'current_cell_annotations.rds'))
  for (pair in list(c('updated_manifest', 'subgroup_annotation_summary.tsv'), c('all_mappings', 'cluster_annotations.tsv'), c('checks', 'annotation_validation.tsv')))
    write.table(do.call(rbind, get(pair[1])), file.path(folder, pair[2]), sep = '\t', quote = FALSE, row.names = FALSE, fileEncoding = 'UTF-8')
  log_msg('ALL_LEGACY_MARKER_QS_OK')
}

step_merge_legacy_marker_annotations <- function() {
  suppressPackageStartupMessages(library(Seurat))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- CONFIG$legacy_marker_annotation; output <- file.path(PROJECT_DIR, s$evidence_dir)
  new <- readRDS(file.path(output, 'current_cell_annotations.rds'))
  input <- file.path(PROJECT_DIR, s$full_input); source_hash <- digest::digest(file = input, algo = 'sha256')
  sc <- qs::qread(input, nthreads = 4); old <- sc[[]]; before <- annotation_fingerprints(sc)
  ix <- match(colnames(sc), new$cell_id); stopifnot(!anyNA(ix), nrow(new) == ncol(sc))
  update <- !old$analysis_container %in% unlist(s$excluded_subgroups)
  cols <- c('subtype_annotation', 'biological_lineage', 'annotation_confidence', 'annotation_PMIDs',
    'annotation_limitations', 'literature_cluster_label', 'annotation_stage')
  for (column in cols) {
    sc[[paste0('pre_legacy_marker_', column)]] <- old[[column]]
    values <- as.character(old[[column]]); values[update] <- as.character(new[[column]][ix[update]])
    sc[[column]] <- values
  }
  sc$annotation_name_source <- new$annotation_name_source[ix]
  sc$legacy_marker_annotation_stage <- new$legacy_marker_annotation_stage[ix]
  Idents(sc) <- 'subtype_annotation'
  sc@misc$legacy_marker_annotation <- list(settings = s, input_sha256 = source_hash,
    current_cluster_column = 'current_subgroup_cluster', excluded_subgroups = s$excluded_subgroups,
    global_reductions_are_historical = TRUE, metadata_join = 'cell_ID')
  stopifnot(identical(sc[[]][, setdiff(names(old), cols)], old[, setdiff(names(old), cols)]),
    identical(sc[[]][!update, names(old)], old[!update, ]))
  save_annotation_verified(sc, file.path(PROJECT_DIR, s$full_output), before)
  stopifnot(identical(digest::digest(file = input, algo = 'sha256'), source_hash))
  write.table(data.frame(cells = ncol(sc), genes = nrow(sc), checkpoint = s$full_output,
    input_sha256 = source_hash, output_sha256 = digest::digest(file = file.path(PROJECT_DIR, s$full_output), algo = 'sha256'),
    updated_cells = sum(update), excluded_cells = sum(!update), neutrophil_prior_metadata_unchanged = TRUE,
    RNA_reductions_graphs_preserved = TRUE, readback_verified = TRUE),
    file.path(output, 'complete_object_validation.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  log_msg('COMPLETE_LEGACY_MARKER_QS_OK')
}

step_plot_legacy_marker_annotations <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
  s <- CONFIG$legacy_marker_annotation; evidence <- file.path(PROJECT_DIR, s$evidence_dir)
  output <- file.path(PROJECT_DIR, CONFIG$current_results$subgroups)
  manifest <- read.delim(file.path(evidence, 'subgroup_annotation_summary.tsv'))
  panels <- read.delim(file.path(PROJECT_DIR, s$marker_panel))
  for (group in setdiff(manifest$subgroup, unlist(s$excluded_subgroups))) {
    cached <- readRDS(file.path(evidence, group, 'annotation_plot_data.rds')); sc <- cached$plot_object; m <- cached$mapping
    umap <- Embeddings(sc, 'umap'); xlim <- range(umap[,1]); ylim <- range(umap[,2])
    folder <- file.path(output, group); palette <- npg_palette(levels(sc$subtype_annotation))
    p <- cellname_annotation_umap(sc, gsub('_', ' ', group))
    save_umap_plot(p, file.path(folder, 'UMAP_annotated'), width = 10, height = 6, xlim = xlim, ylim = ylim)
    for (column in c('sample', 'group')) {
      p <- DimPlot(sc, group.by = 'subtype_annotation', split.by = column, cols = palette, raster = FALSE,
        shuffle = TRUE, seed = SEED, ncol = if (column == 'sample') 3 else 2) &
        theme_classic(base_size = 11) & theme(legend.text = element_text(size = 9)) &
        labs(title = gsub('_', ' ', group), subtitle = NULL) & guides(colour = guide_legend(override.aes = list(size = 3)))
      save_umap_plot(p, file.path(folder, paste0('UMAP_', column, '_split')),
        width = if (column == 'sample') 17 else 13, height = if (column == 'sample') 10 else 6, xlim = xlim, ylim = ylim)
    }
    e <- readRDS(file.path(PROJECT_DIR, s$input_evidence, group, 'expression_evidence.rds'))
    m$candidate_annotation <- m$annotation
    genes <- strsplit(panels$genes[panels$subgroup == group], ';', fixed = TRUE)[[1]]
    p <- candidate_dotplot(e, m, genes, gsub('_', ' ', group))
    direct <- read.delim(file.path(evidence, group, 'Seurat_DotPlot_reference.tsv'))
    ix <- match(paste(p$data$cluster, p$data$gene), paste(direct$id, direct$features.plot))
    stopifnot(!anyNA(ix), max(abs(p$data$average_expression - direct$avg.exp[ix])) < 1e-8,
      max(abs(p$data$color_value - direct$avg.exp.scaled[ix])) < 1e-8)
    save_plot(p$plot, file.path(folder, 'DotPlot_annotation_markers'),
      width = max(12, length(genes) * .23 + 5), height = max(4.5, nrow(m) * .3 + 2.1))
    write.table(p$data, file.path(evidence, group, 'DotPlot_annotation_data.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    log_msg('LEGACY_MARKER_PLOTS_OK ', group)
  }
  file.copy(file.path(evidence, 'cluster_annotations.tsv'), file.path(output, 'Cluster_annotations.tsv'), overwrite = TRUE)
}
