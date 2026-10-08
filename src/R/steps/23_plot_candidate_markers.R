# 以全细胞表达均值/比例绘制候选注释点图；保持cluster逐行展示，不合并复用标签的群。
candidate_dotplot <- function(evidence, mapping, genes, title, scaled = TRUE) {
  stopifnot(all(genes %in% rownames(evidence$average)))
  d <- expand.grid(gene = genes, cluster = as.character(mapping$subgroup_cluster), stringsAsFactors = FALSE)
  d$average_expression <- as.vector(evidence$average[genes, as.character(mapping$subgroup_cluster), drop = FALSE])
  d$percent_expressed <- as.vector(evidence$percent[genes, as.character(mapping$subgroup_cluster), drop = FALSE])
  d$color_value <- log1p(d$average_expression)
  if (scaled) d$color_value <- ave(d$color_value, d$gene, FUN = function(x)
    if (sd(x) == 0) rep(0, length(x)) else pmax(-2.5, pmin(2.5, as.numeric(scale(x)))))
  # cluster仍作为独立行键；同名群不合并，图上只显示生物学名称。
  row_labels <- setNames(mapping$candidate_annotation, mapping$subgroup_cluster)
  p <- ggplot(d, aes(gene, cluster, size = percent_expressed, color = color_value)) +
    geom_point() + scale_size_area(max_size = 4.4, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_x_discrete(limits = genes, drop = FALSE) +
    scale_y_discrete(limits = rev(as.character(mapping$subgroup_cluster)), labels = row_labels, drop = FALSE) +
    labs(x = NULL, y = NULL, title = title, size = "Expressed (%)",
      color = if (scaled) "Scaled mean" else "log1p mean")
  p <- hongf_dotplot_style(p) + theme(panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = .25),
    axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9), plot.title = element_text(size = 12))
  if (scaled) p <- p + scale_color_gradientn(colours = HONGF_DOTPLOT_COLORS, limits = c(-2.5, 2.5),
    guide = guide_colorbar(display = "rectangles"))
  bars <- ggplot(mapping, aes(cells, as.character(subgroup_cluster))) +
    geom_col(width = .62, fill = "#F39B7F", colour = "#A13B2F", linewidth = .25) +
    geom_text(aes(label = format(cells, big.mark = ",", trim = TRUE)), hjust = -.1, size = 3) +
    scale_y_discrete(limits = rev(as.character(mapping$subgroup_cluster)), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, .65))) + coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(5.5, 12, 5.5, 2))
  p <- patchwork::wrap_plots(p, bars, widths = c(9, 1.4), guides = "collect") & theme(legend.position = "bottom")
  list(plot = p, data = d)
}

step_plot_candidate_markers <- function() {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  settings <- CONFIG$subgroup_marker_review
  output <- file.path(PROJECT_DIR, settings$output_dir)
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest))
  mapping <- read.delim(file.path(PROJECT_DIR, settings$mapping), encoding = "UTF-8", stringsAsFactors = FALSE)
  panel <- read.delim(file.path(PROJECT_DIR, settings$marker_panel), stringsAsFactors = FALSE)
  programs <- read.delim(file.path(PROJECT_DIR, settings$lineage_programs), stringsAsFactors = FALSE)
  stopifnot(nrow(mapping) == sum(manifest$clusters), !anyDuplicated(paste(mapping$subgroup, mapping$subgroup_cluster)))
  all_maps <- all_support <- all_missing <- all_cells <- list()
  pdf_path <- file.path(output, "All_subgroups_candidate_DotPlots.pdf")
  grDevices::pdf(pdf_path, width = 18, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; folder <- file.path(output, group)
    e <- readRDS(file.path(folder, "expression_evidence.rds"))
    m <- mapping[mapping$subgroup == group, ]; m <- m[order(m$subgroup_cluster), ]
    ids <- as.character(e$metadata[[settings$cluster_column]])
    m$cells <- as.integer(table(factor(ids, levels = m$subgroup_cluster)))
    m$stage <- settings$stage
    stopifnot(all(m$cells > 0), sum(m$cells) == manifest$cells[i], setequal(m$subgroup_cluster, e$clusters))
    genes <- unique(strsplit(panel$genes[panel$subgroup == group], ";", fixed = TRUE)[[1]])
    missing <- setdiff(genes, rownames(e$average))
    all_missing[[group]] <- data.frame(subgroup = rep(group, length(missing)), gene = missing)
    genes <- intersect(genes, rownames(e$average))
    p <- candidate_dotplot(e, m, genes, gsub("_", " ", group))
    save_plot(p$plot, file.path(folder, "DotPlot_candidate_markers"),
      width = max(12, length(genes) * .22 + 5.2), height = max(4, nrow(m) * .29 + 1.9))
    print(p$plot)
    write.table(p$data, file.path(folder, "DotPlot_candidate_markers_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # 候选名称的每个支持基因均保留实际检测比例与平均表达，低表达不被隐藏。
    support <- do.call(rbind, lapply(seq_len(nrow(m)), function(j) {
      features <- strsplit(m$support_markers[j], ";", fixed = TRUE)[[1]]
      stopifnot(all(features %in% rownames(e$average)))
      data.frame(subgroup = group, cluster = m$subgroup_cluster[j], candidate = m$candidate_annotation[j], gene = features,
        percent_expressed = e$percent[features, as.character(m$subgroup_cluster[j])],
        average_expression = e$average[features, as.character(m$subgroup_cluster[j])])
    }))
    all_support[[group]] <- support
    # 逐样本核对panel表达，并检验混合谱系是否出现在同一个条码；不由此直接判定双细胞。
    input <- file.path(PROJECT_DIR, manifest$checkpoint[i])
    sc <- run_step(paste(group, "Read for sample and lineage validation"), qs::qread(input, nthreads = 4))
    stopifnot(identical(sc[[]], e$metadata), identical(Embeddings(sc, "umap"), e$umap))
    tissue_genes <- if (group == "Epithelial") strsplit(panel$genes[panel$subgroup == "Tissue_origin"], ";", fixed = TRUE)[[1]] else character()
    program_genes <- unique(unlist(strsplit(programs$genes, ";", fixed = TRUE)))
    required <- intersect(unique(c(genes, tissue_genes, program_genes)), rownames(sc))
    expr <- LayerData(sc, assay = "RNA", layer = "data")[required, , drop = FALSE]
    keys <- paste(ids, sc$sample, sep = "|")
    sample_stats <- do.call(rbind, lapply(unique(keys), function(key) {
      ix <- which(keys == key); x <- expr[intersect(c(genes, tissue_genes), rownames(expr)), ix, drop = FALSE]
      linear <- x; linear@x <- expm1(linear@x)
      data.frame(cluster = ids[ix[1]], sample = sc$sample[ix[1]], cells = length(ix), gene = rownames(x),
        percent_expressed = Matrix::rowMeans(x > 0) * 100, average_expression = Matrix::rowMeans(linear))
    }))
    write.table(sample_stats, file.path(folder, "marker_expression_by_cluster_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    positive <- sapply(seq_len(nrow(programs)), function(j) {
      g <- strsplit(programs$genes[j], ";", fixed = TRUE)[[1]]
      stopifnot(all(g %in% rownames(expr)))
      Matrix::colSums(expr[g, , drop = FALSE] > 0) >= programs$minimum_detected[j]
    })
    colnames(positive) <- programs$program
    pairs <- combn(programs$program, 2)
    paired <- sapply(seq_len(ncol(pairs)), function(j) positive[, pairs[1,j]] & positive[, pairs[2,j]])
    colnames(paired) <- apply(pairs, 2, paste, collapse = "_and_")
    lineage <- aggregate(cbind(positive, paired) * 100, list(cluster = ids), mean)
    write.table(lineage, file.path(folder, "lineage_coexpression_percent.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    qc_columns <- intersect(c("nCount_RNA", "nFeature_RNA", "pct_counts_mt", "doublet_score"), names(e$metadata))
    qc <- aggregate(e$metadata[, qc_columns, drop = FALSE], list(subgroup_cluster = ids), median)
    m <- merge(m, qc, by = "subgroup_cluster", sort = FALSE); m <- m[order(as.integer(m$subgroup_cluster)), ]
    counts <- as.data.frame(table(cluster = ids, sample = sc$sample)); counts <- counts[counts$Freq > 0, ]
    m$samples_present <- as.integer(table(factor(counts$cluster, levels = m$subgroup_cluster)))
    m$max_sample_fraction <- vapply(as.character(m$subgroup_cluster), function(id) max(counts$Freq[counts$cluster == id]) / sum(counts$Freq[counts$cluster == id]), numeric(1))
    index <- match(ids, m$subgroup_cluster)
    metadata <- data.frame(cell_id = colnames(sc), e$metadata,
      candidate_annotation = m$candidate_annotation[index], candidate_confidence = m$confidence[index],
      candidate_review_flags = m$review_flags[index], annotation_stage = settings$stage,
      UMAP_1 = e$umap[,1], UMAP_2 = e$umap[,2])
    stopifnot(!anyNA(metadata$candidate_annotation), nrow(metadata) == ncol(sc))
    connection <- gzfile(file.path(folder, "cell_candidate_annotation.tsv.gz"), "wt")
    write.table(metadata, connection, sep = "\t", quote = FALSE, row.names = FALSE); close(connection)
    # 候选层仅保存紧凑元数据及映射；人工确认前不覆盖旧QS或提升为最终注释。
    saveRDS(list(metadata = metadata, mapping = m, source = manifest$checkpoint[i],
      input_sha256 = digest::digest(file = input, algo = "sha256"),
      stage = settings$stage, unchanged_RNA_and_reductions = TRUE), file.path(folder, "candidate_annotation.rds"))
    if (group == "Epithelial") {
      missing_tissue <- setdiff(tissue_genes, rownames(e$average))
      all_missing$Tissue_origin <- data.frame(subgroup = rep("Tissue_origin", length(missing_tissue)), gene = missing_tissue)
      p_tissue <- candidate_dotplot(e, m, intersect(tissue_genes, rownames(e$average)), "Tissue origin markers", scaled = FALSE)
      save_plot(p_tissue$plot, file.path(output, "DotPlot_tissue_origin"), width = 13, height = 5)
      write.table(p_tissue$data, file.path(output, "tissue_origin_marker_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      weight <- m$cells / sum(m$cells)
      d <- data.frame(gene = tissue_genes, present_in_object = tissue_genes %in% rownames(e$average), percent_expressed = NA_real_, average_expression = NA_real_)
      ix <- which(d$present_in_object)
      d$percent_expressed[ix] <- as.vector(e$percent[d$gene[ix], as.character(m$subgroup_cluster)] %*% weight)
      d$average_expression[ix] <- as.vector(e$average[d$gene[ix], as.character(m$subgroup_cluster)] %*% weight)
      write.table(d, file.path(output, "tissue_origin_weighted_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    all_maps[[group]] <- m; all_cells[[group]] <- colnames(sc)
    write.table(m, file.path(folder, "candidate_cluster_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    log_msg("CANDIDATE_MARKER_PLOT_OK ", group, "; cells=", ncol(sc), "; markers=", length(genes))
    rm(sc, expr, metadata, sample_stats, e, p); invisible(gc())
  }
  stopifnot(!anyDuplicated(unlist(all_cells)), length(unlist(all_cells)) == sum(manifest$cells))
  write.table(do.call(rbind, all_maps), file.path(output, "all_candidate_cluster_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write.table(do.call(rbind, all_support), file.path(output, "candidate_support_marker_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, all_missing), file.path(output, "markers_not_in_object.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  yaml::write_yaml(list(stage = settings$stage, total_cells = sum(manifest$cells), clusters = nrow(mapping),
    input_objects_modified = FALSE, previous_labels_used_for_annotation = FALSE,
    treatment_used_for_annotation = FALSE, marker_test = "Seurat RNA Wilcoxon, all cells, Bonferroni adjusted",
    settings = settings, tissue_hypothesis = "intestinal/colorectal-like epithelial origin; model unresolved"),
    file.path(output, "analysis_status.yaml"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "18_plot_candidate_markers_sessionInfo.txt"))
  log_msg("ALL_CANDIDATE_DOTPLOTS_OK; awaiting user marker confirmation")
}

# 用Seurat原生DotPlot独立复算图中数值，并核对逐细胞注释、来源哈希及所有坐标。
validate_candidate_marker_review <- function() {
  suppressPackageStartupMessages(library(Seurat))
  settings <- CONFIG$subgroup_marker_review
  output <- file.path(PROJECT_DIR, settings$output_dir)
  manifest <- read.delim(file.path(PROJECT_DIR, settings$input_manifest))
  previous <- read.delim(file.path(PROJECT_DIR, CONFIG$subgroup_previous_annotation$output_dir, "subgroup_annotation_summary.tsv"))
  config_mapping <- read.delim(file.path(PROJECT_DIR, settings$mapping), encoding = "UTF-8")
  completed <- list()
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$subgroup[i]; folder <- file.path(output, group)
    candidate <- readRDS(file.path(folder, "candidate_annotation.rds"))
    stopifnot(identical(candidate$input_sha256, previous$input_sha256[match(group, previous$subgroup)]))
    sc <- qs::qread(file.path(PROJECT_DIR, candidate$source), nthreads = 4)
    stopifnot(identical(colnames(sc), candidate$metadata$cell_id),
      identical(as.numeric(as.matrix(candidate$metadata[, c("UMAP_1", "UMAP_2")])), as.numeric(Embeddings(sc, "umap"))),
      isTRUE(all.equal(sc[[]], candidate$metadata[, names(sc[[]]), drop = FALSE], check.attributes = FALSE)))
    m <- read.delim(file.path(folder, "candidate_cluster_annotation.tsv"), encoding = "UTF-8")
    reference <- config_mapping[config_mapping$subgroup == group, ]
    stopifnot(identical(m$note_CN, reference$note_CN[match(m$subgroup_cluster, reference$subgroup_cluster)]))
    observed <- read.delim(file.path(folder, "DotPlot_candidate_markers_data.tsv"))
    Idents(sc) <- settings$cluster_column
    direct <- suppressWarnings(DotPlot(sc, features = unique(observed$gene), assay = "RNA"))$data
    key <- paste(observed$cluster, observed$gene)
    index <- match(key, paste(direct$id, direct$features.plot))
    stopifnot(!anyNA(index), !anyDuplicated(key),
      max(abs(observed$percent_expressed - direct$pct.exp[index])) < 1e-8,
      max(abs(observed$average_expression - direct$avg.exp[index])) < 1e-8,
      max(abs(observed$color_value - direct$avg.exp.scaled[index])) < 1e-8,
      sum(m$cells) == ncol(sc), nrow(m) == manifest$clusters[i],
      all(candidate$metadata$annotation_stage == settings$stage))
    completed[[group]] <- data.frame(subgroup = group, cells = ncol(sc), clusters = nrow(m),
      plot_values = nrow(observed), Seurat_DotPlot_match = TRUE, source_hash_unchanged = TRUE,
      metadata_UMAP_unchanged = TRUE, Chinese_notes_match = TRUE)
    log_msg("CANDIDATE_REVIEW_VERIFIED ", group, " rows=", nrow(observed))
    rm(sc, direct, candidate); invisible(gc())
  }
  write.table(do.call(rbind, completed), file.path(output, "verification_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
