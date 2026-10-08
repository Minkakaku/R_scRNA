# 用户确认删除指定混杂/低RNA群；仅按已有cluster筛选，不重新降维或重新聚类。
step_clean_candidate_clusters <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$subgroup_cleaning
  review <- CONFIG$subgroup_marker_review
  output <- file.path(PROJECT_DIR, settings$output_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROJECT_DIR, settings$checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
  exclusions <- read.delim(file.path(PROJECT_DIR, settings$exclusions))
  manifest <- read.delim(file.path(PROJECT_DIR, review$input_manifest))
  panel <- read.delim(file.path(PROJECT_DIR, review$marker_panel))
  stopifnot(nrow(exclusions) == settings$expected_removed_clusters)
  destinations <- file.path(PROJECT_DIR, settings$checkpoint_dir,
    paste0(settings$object_prefix, ".", CONFIG$project_id, "_", manifest$subgroup, "_cleaned_candidate.qs"))
  stopifnot(!any(file.exists(c(destinations, file.path(PROJECT_DIR, settings$full_output)))))
  # 每个完整QS先写临时文件并读回校验，验证成功后再发布正式文件名。
  save_checked <- function(sc, destination) {
    expected <- c(RNA = digest::digest(sc[["RNA"]], algo = "xxhash64"),
      reductions = digest::digest(sc@reductions, algo = "xxhash64"),
      metadata = digest::digest(sc[[]], algo = "xxhash64"),
      graphs = digest::digest(sc@graphs, algo = "xxhash64"))
    stopifnot(!file.exists(destination), !file.exists(paste0(destination, ".partial")))
    qs::qsave(sc, paste0(destination, ".partial"), preset = "high", nthreads = 4)
    check <- qs::qread(paste0(destination, ".partial"), nthreads = 4)
    got <- c(RNA = digest::digest(check[["RNA"]], algo = "xxhash64"),
      reductions = digest::digest(check@reductions, algo = "xxhash64"),
      metadata = digest::digest(check[[]], algo = "xxhash64"),
      graphs = digest::digest(check@graphs, algo = "xxhash64"))
    stopifnot(identical(got, expected), isTRUE(validObject(check, test = TRUE)),
      identical(colnames(check), colnames(sc)), file.rename(paste0(destination, ".partial"), destination))
    invisible(TRUE)
  }
  all_metadata <- all_removed <- all_maps <- summaries <- all_sample_counts <- plots <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]
    source_folder <- file.path(PROJECT_DIR, settings$evidence_dir, group)
    folder <- file.path(output, group); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    candidate <- readRDS(file.path(source_folder, "candidate_annotation.rds"))
    sc <- run_step(paste(group, "Read to remove excluded clusters"), qs::qread(file.path(PROJECT_DIR, candidate$source), nthreads = 4))
    stopifnot(identical(colnames(sc), candidate$metadata$cell_id))
    excluded <- exclusions[exclusions$subgroup == group, ]
    original_cluster <- as.character(sc[[]][[review$cluster_column]])
    removed <- original_cluster %in% as.character(excluded$subgroup_cluster)
    stopifnot(setequal(unique(original_cluster[removed]), as.character(excluded$subgroup_cluster)))
    removed_metadata <- data.frame(subgroup = group, candidate$metadata[removed, ])
    removed_metadata$exclusion_reason <- excluded$reason[match(original_cluster[removed], excluded$subgroup_cluster)]
    metadata <- candidate$metadata[!removed, ]
    old_cells <- colnames(sc); old_genes <- rownames(sc); source_meta <- sc[[]][!removed, , drop = FALSE]
    clean <- subset(sc, cells = old_cells[!removed])
    stopifnot(identical(colnames(clean), metadata$cell_id), identical(rownames(clean), old_genes),
      identical(clean[[]], source_meta))
    for (layer in Layers(sc[["RNA"]])) stopifnot(isTRUE(all.equal(
      LayerData(clean, assay = "RNA", layer = layer),
      LayerData(sc, assay = "RNA", layer = layer)[, !removed, drop = FALSE], check.attributes = TRUE)))
    for (reduction in Reductions(sc)) stopifnot(identical(as.numeric(Embeddings(clean, reduction)),
      as.numeric(Embeddings(sc, reduction)[!removed, , drop = FALSE])))
    rm(sc); invisible(gc())
    m <- candidate$mapping[!candidate$mapping$subgroup_cluster %in% excluded$subgroup_cluster, ]
    m$stage <- settings$stage
    index <- match(as.character(clean[[]][[review$cluster_column]]), m$subgroup_cluster)
    clean$candidate_annotation <- m$candidate_annotation[index]
    clean$candidate_annotation_label <- factor(paste0("C", m$subgroup_cluster[index], ": ", m$candidate_annotation[index]),
      levels = paste0("C", m$subgroup_cluster, ": ", m$candidate_annotation))
    clean$candidate_confidence <- m$confidence[index]
    clean$candidate_review_flags <- m$review_flags[index]
    clean$annotation_stage <- settings$stage
    Idents(clean) <- "candidate_annotation_label"
    clean@misc$subgroup_candidate_annotation <- list(mapping = m, annotation_stage = settings$stage,
      exclusions = excluded, original_clusters_preserved = TRUE, reclustering_performed = FALSE,
      source = candidate$source, source_sha256 = candidate$input_sha256)
    run_step(paste(group, "Save and verify cleaned QS"), save_checked(clean, destinations[i]))
    # 同一cluster内未移除任何单细胞，故均值和比例与原证据完全相同；仅对保留群重新缩放颜色。
    e <- readRDS(file.path(source_folder, "expression_evidence.rds"))
    genes <- strsplit(panel$genes[panel$subgroup == group], ";", fixed = TRUE)[[1]]
    p <- candidate_dotplot(e, m, genes, gsub("_", " ", group))
    direct <- suppressWarnings(DotPlot(clean, features = genes, assay = "RNA", group.by = review$cluster_column))$data
    match_index <- match(paste(p$data$cluster, p$data$gene), paste(direct$id, direct$features.plot))
    stopifnot(!anyNA(match_index), max(abs(p$data$average_expression - direct$avg.exp[match_index])) < 1e-8,
      max(abs(p$data$percent_expressed - direct$pct.exp[match_index])) < 1e-8,
      max(abs(p$data$color_value - direct$avg.exp.scaled[match_index])) < 1e-8)
    save_plot(p$plot, file.path(folder, "DotPlot_cleaned_candidate_markers"),
      width = max(12, length(genes) * .22 + 5.2), height = max(4, nrow(m) * .29 + 1.9))
    plots[[group]] <- p$plot
    write.table(p$data, file.path(folder, "DotPlot_cleaned_candidate_markers_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(m, file.path(folder, "retained_cluster_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    metadata$annotation_stage <- settings$stage
    connection <- gzfile(file.path(folder, "cell_candidate_annotation.tsv.gz"), "wt")
    write.table(metadata, connection, sep = "\t", quote = FALSE, row.names = FALSE); close(connection)
    sample_counts <- as.data.frame(table(sample = candidate$metadata$sample, retained = !removed))
    sample_counts$subgroup <- group
    summaries[[group]] <- data.frame(subgroup = group, before_cells = nrow(candidate$metadata),
      removed_cells = sum(removed), retained_cells = ncol(clean), removed_clusters = nrow(excluded),
      retained_clusters = nrow(m), checkpoint = substring(destinations[i], nchar(PROJECT_DIR) + 2),
      readback_verified = TRUE, RNA_reductions_preserved = TRUE, native_DotPlot_verified = TRUE)
    all_metadata[[group]] <- metadata; all_removed[[group]] <- removed_metadata
    all_maps[[group]] <- m; all_sample_counts[[group]] <- sample_counts
    write.table(do.call(rbind, summaries), file.path(output, "cleaning_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("CLEAN_SUBGROUP_OK ", group, "; removed=", sum(removed), "; retained=", ncol(clean))
    rm(clean, e, p, candidate, direct); invisible(gc())
  }
  kept <- do.call(rbind, all_metadata); dropped <- do.call(rbind, all_removed)
  stopifnot(nrow(kept) == settings$expected_retained_cells, nrow(dropped) == settings$expected_removed_cells,
    !anyDuplicated(kept$cell_id), !anyDuplicated(dropped$cell_id), !any(kept$cell_id %in% dropped$cell_id))
  connection <- gzfile(file.path(output, "removed_cells.tsv.gz"), "wt")
  write.table(dropped, connection, sep = "\t", quote = FALSE, row.names = FALSE); close(connection)
  removed_counts <- as.data.frame(table(subgroup = dropped$subgroup, subgroup_cluster = dropped$subgroup_cluster))
  removed_counts <- removed_counts[removed_counts$Freq > 0, ]; names(removed_counts)[3] <- "removed_cells"
  removed_counts <- merge(exclusions, removed_counts, by = c("subgroup", "subgroup_cluster"))
  write.table(removed_counts, file.path(output, "removed_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, all_maps), file.path(output, "retained_cluster_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write.table(do.call(rbind, all_sample_counts), file.path(output, "sample_cell_retention.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  grDevices::pdf(file.path(output, "All_subgroups_cleaned_DotPlots.pdf"), width = 18, height = 7)
  for (p in plots) print(p)
  grDevices::dev.off()
  # 完整合并对象保留母对象的全局降维；各大类局部坐标另存metadata，避免误认为统一空间。
  full <- run_step("Read full parent for identical cell filtering", qs::qread(file.path(PROJECT_DIR, settings$full_input), nthreads = 4))
  stopifnot(setequal(colnames(full), c(kept$cell_id, dropped$cell_id)))
  retained_order <- colnames(full)[colnames(full) %in% kept$cell_id]
  before_umap <- Embeddings(full, "umap")[retained_order, , drop = FALSE]
  clean_full <- subset(full, cells = retained_order)
  stopifnot(identical(colnames(clean_full), retained_order), identical(rownames(clean_full), rownames(full)))
  for (layer in Layers(full[["RNA"]])) stopifnot(isTRUE(all.equal(LayerData(clean_full, assay = "RNA", layer = layer),
    LayerData(full, assay = "RNA", layer = layer)[, retained_order, drop = FALSE])))
  stopifnot(identical(as.numeric(Embeddings(clean_full, "umap")), as.numeric(before_umap)))
  rm(full); invisible(gc())
  k <- kept[match(colnames(clean_full), kept$cell_id), ]
  stopifnot(identical(as.character(clean_full$sample), as.character(k$sample)))
  clean_full$group_original <- as.character(clean_full$group)
  clean_full$group <- factor(k$group, levels = unlist(CONFIG$group_levels))
  clean_full$subgroup_cluster <- as.character(k$subgroup_cluster)
  clean_full$subgroup_cluster_id <- paste(k$broad_cell_type, paste0("C", k$subgroup_cluster), sep = ":")
  for (column in c("candidate_annotation", "candidate_confidence", "candidate_review_flags", "annotation_stage")) clean_full[[column]] <- k[[column]]
  clean_full$subgroup_umap_1 <- k$UMAP_1; clean_full$subgroup_umap_2 <- k$UMAP_2
  Idents(clean_full) <- "candidate_annotation"
  clean_full@misc$subgroup_candidate_annotation <- list(stage = settings$stage, settings = settings,
    selected_pcs = CONFIG$subgroup_pca$selected_pcs, exclusions = removed_counts,
    global_reductions_source = settings$full_input, subgroup_reductions_in_separate_checkpoints = TRUE,
    reclustering_performed = FALSE, subgroup_markers_recomputed_after_filtering = FALSE,
    annotation_confirmed_by_user = FALSE, removal_confirmed_by_user = TRUE)
  run_step("Save and verify complete cleaned QS", save_checked(clean_full, file.path(PROJECT_DIR, settings$full_output)))
  write.table(as.data.frame(table(sample = clean_full$sample, group = clean_full$group)),
    file.path(output, "complete_cleaned_sample_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  yaml::write_yaml(list(stage = settings$stage, input_cells = nrow(kept) + nrow(dropped), retained_cells = nrow(kept),
    removed_cells = nrow(dropped), removed_clusters = nrow(exclusions), retained_clusters = sum(vapply(all_maps, nrow, integer(1))),
    full_object = settings$full_output, readback_verified = TRUE, original_inputs_preserved = TRUE,
    reclustering_performed = FALSE, marker_confirmation_pending = TRUE), file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "23_clean_candidate_clusters_sessionInfo.txt"))
  log_msg("ALL_CLEANED_OBJECTS_VERIFIED; retained=", nrow(kept), "; removed=", nrow(dropped))
}
