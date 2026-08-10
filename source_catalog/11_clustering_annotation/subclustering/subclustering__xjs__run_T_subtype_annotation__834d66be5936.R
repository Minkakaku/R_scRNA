suppressPackageStartupMessages({
  required_pkgs <- c(
    "Seurat", "tidyverse", "harmony", "future", "openxlsx",
    "tibble", "data.table", "patchwork"
  )
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nNo packages were installed. Please install/activate them first if needed.",
      call. = FALSE
    )
  }
  
  library(Seurat)
  library(tidyverse)
  library(harmony)
  library(future)
})

set.seed(1234)
options(future.globals.maxSize = 10 * 1024^3)

input_rds <- "all_data2_TLS.rds"
out_dir <- "T_subtype_codex_result"
qc_dir <- file.path(out_dir, "01_qc")
cluster_dir <- file.path(out_dir, "02_cluster")
marker_dir <- file.path(out_dir, "03_marker")
anno_dir <- file.path(out_dir, "04_annotation")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(anno_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading: ", input_rds)
all_data <- readRDS(input_rds)

if (!inherits(all_data, "Seurat")) {
  stop("The input RDS is not a Seurat object: ", input_rds, call. = FALSE)
}

meta_cols <- colnames(all_data[[]])
if (!"annotation" %in% meta_cols) {
  stop("Metadata column `annotation` was not found in all_data2_TLS.rds.", call. = FALSE)
}
if (!"sample" %in% meta_cols) {
  stop("Metadata column `sample` was not found; it is required for Harmony batch correction.", call. = FALSE)
}

message("annotation table:")
print(table(all_data$annotation, useNA = "ifany"))

theme_90angle <- theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

join_layers_if_available <- function(object) {
  if (!exists("JoinLayers", mode = "function")) {
    return(object)
  }
  tryCatch(
    JoinLayers(object),
    error = function(e) {
      message("JoinLayers skipped: ", conditionMessage(e))
      object
    }
  )
}

select_pca_dims <- function(object, reduction = "pca", max_dims = 50, min_dims = 10, max_selected = 30) {
  stdev <- object[[reduction]]@stdev
  max_dims <- min(max_dims, length(stdev))
  stdev <- stdev[seq_len(max_dims)]
  variance <- stdev^2
  pct_variance <- variance / sum(variance) * 100
  cumulative_variance <- cumsum(pct_variance)
  
  x <- seq_along(stdev)
  y <- stdev
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  line_y <- y_scaled[1] + (y_scaled[length(y_scaled)] - y_scaled[1]) * x_scaled
  elbow_dim <- which.max(abs(y_scaled - line_y))
  
  variance_dim <- which(cumulative_variance >= 90 & pct_variance < 5)[1]
  if (is.na(variance_dim)) {
    variance_dim <- elbow_dim
  }
  
  selected_dim <- max(min_dims, min(max_selected, max(elbow_dim, variance_dim)))
  
  list(
    selected_dim = selected_dim,
    elbow_dim = elbow_dim,
    variance_dim = variance_dim,
    stats = tibble(
      PC = x,
      stdev = stdev,
      pct_variance = pct_variance,
      cumulative_variance = cumulative_variance,
      selected = x == selected_dim,
      elbow = x == elbow_dim,
      variance_rule = x == variance_dim
    )
  )
}

filter_features <- function(object, marker_vec) {
  marker_vec <- marker_vec[marker_vec %in% rownames(object)]
  duplicated_genes <- unique(unname(marker_vec)[duplicated(unname(marker_vec))])
  if (length(duplicated_genes) > 0) {
    message("Duplicated DotPlot markers removed: ", paste(duplicated_genes, collapse = ", "))
    marker_vec <- marker_vec[!duplicated(unname(marker_vec))]
  }
  marker_vec
}

safe_dotplot <- function(object, group.by, features, title = NULL, dot.scale = 6) {
  features_use <- filter_features(object, features)
  if (length(features_use) == 0) {
    warning("No requested features found for group.by=", group.by)
    return(ggplot() + theme_void() + ggtitle(paste("No markers found:", title %||% group.by)))
  }
  DotPlot(object, group.by = group.by, features = features_use, dot.scale = dot.scale) +
    RotatedAxis() +
    ggtitle(title %||% group.by)
}

save_pdf <- function(path, plot_expr, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(plot_expr)
}

major_marker_ls <- c(
  "Epithelial" = "EPCAM", "Epithelial" = "KRT5", "Epithelial" = "KRT14",
  "Endothelial" = "PECAM1", "Endothelial" = "VWF",
  "Fibroblast" = "PDGFRA", "Fibroblast" = "CCDC80", "Fibroblast" = "SFRP4",
  "Fibroblast" = "COL1A1", "Fibroblast" = "COL1A2", "Fibroblast" = "DCN",
  "B cell" = "CD19", "B cell" = "CD79A", "B cell" = "MS4A1",
  "Plasma" = "TNFRSF17", "Plasma" = "JCHAIN", "Plasma" = "MZB1",
  "T cell" = "CD3D", "T cell" = "CD3E", "T cell" = "CD3G", "T cell" = "TRAC",
  "NK" = "GNLY", "NK" = "NKG7", "NK" = "NCAM1",
  "Proliferating" = "MKI67", "Proliferating" = "TOP2A",
  "Monocyte" = "CD14", "Monocyte" = "FCGR3A",
  "Macrophage" = "CD68", "Macrophage" = "CD163", "Macrophage" = "C1QA",
  "Dendritic cell" = "CD83", "Dendritic cell" = "CLEC9A", "Dendritic cell" = "CD1C",
  "Mast cell" = "CPA3", "Mast cell" = "TPSAB1", "Mast cell" = "KIT"
)

T_major_marker_ls <- c(
  "CD8" = "CD8A", "CD8" = "CD8B", "CD8" = "GZMB", "CD8" = "GZMA",
  "CD8" = "GZMK", "CD8" = "GZMH", "CD8" = "CCL5",
  "Naive" = "LEF1", "Naive" = "TCF7", "Naive" = "SELL", "Naive" = "CCR7",
  "Memory" = "GPR183", "Memory" = "IL7R", "Memory" = "FOS", "Memory" = "JUN",
  "TRM" = "CXCR6", "TRM" = "ITGAE", "TRM" = "CD69", "TRM" = "ZNF683",
  "Effector" = "FGFBP2", "Effector" = "TYROBP", "Effector" = "CX3CR1",
  "Effector" = "FCGR3A", "Effector" = "GNLY",
  "Cytotoxic" = "XCL1", "Cytotoxic" = "IFNG", "Cytotoxic" = "CCL4",
  "Exhausted" = "CXCL13", "Exhausted" = "LAYN", "Exhausted" = "LAG3",
  "Exhausted" = "HAVCR2", "Exhausted" = "TNFRSF9", "Exhausted" = "CTLA4",
  "Exhausted" = "TIGIT", "Exhausted" = "TOX", "Exhausted" = "PDCD1",
  "CD4" = "FOXP3", "CD4" = "CD4", "CD4" = "IL2RA",
  "Tfh" = "CXCR5", "Tfh" = "ICOS", "Tfh" = "BATF", "Tfh" = "PDCD1",
  "gdT" = "TRDC",
  "Proliferation" = "MKI67", "Proliferation" = "TOP2A", "Proliferation" = "STMN1"
)

CD4_major_marker_ls <- c(
  "T cell" = "CD3D", "T cell" = "CD3E", "T cell" = "CD3G", "T cell" = "TRAC",
  "CD4" = "CD4",
  "CD4 Treg" = "TNFRSF4", "CD4 Treg" = "TNFRSF9", "CD4 Treg" = "FOXP3",
  "CD4 Treg" = "IL2RA", "CD4 Treg" = "LAYN", "CD4 Treg" = "CTLA4",
  "CD4 Treg" = "TIGIT",
  "CD4 naive" = "LEF1", "CD4 naive" = "CCR7", "CD4 naive" = "TCF7",
  "CD4 Tem" = "GZMK", "CD4 Tem" = "GZMA", "CD4 Tem" = "NKG7",
  "CD4 Th1" = "CXCL13", "CD4 Th1" = "BHLHE40", "CD4 Th1" = "IFNG",
  "CD4 Tfh" = "TOX2", "CD4 Tfh" = "PDCD1", "CD4 Tfh" = "CXCR5",
  "CD4 Th17" = "IL17A", "CD4 Th17" = "IL17F", "CD4 Th17" = "IL23R",
  "CD4 Tcm" = "TCF7", "CD4 Tcm" = "CCR7", "CD4 Tcm" = "S1PR1"
)

CD8_major_marker_ls <- c(
  "CD8" = "CD8A", "CD8" = "CD8B",
  "CD8 Tex" = "HAVCR2", "CD8 Tex" = "CTLA4", "CD8 Tex" = "LAYN",
  "CD8 Tex" = "LAG3", "CD8 Tex" = "TIGIT",
  "CD8 Tem" = "GZMK", "CD8 Tem" = "GZMA", "CD8 Tem" = "IFNG",
  "CD8 Tem" = "KLRG1", "CD8 Tem" = "GZMB", "CD8 Tem" = "PRF1",
  "CD8 Tem" = "NKG7",
  "CD8 Tnaive" = "SELL", "CD8 Tnaive" = "LEF1", "CD8 Tnaive" = "TCF7",
  "CD8 Tnaive" = "CCR7",
  "CD8 MAIT" = "NCR3", "CD8 MAIT" = "SLC4A10", "CD8 MAIT" = "KLRB1",
  "CD8 Trm" = "MYADM", "CD8 Trm" = "CD69", "CD8 Trm" = "ITGAE",
  "CD8 Trm" = "CXCR6", "CD8 Trm" = "ITGA1",
  "CD8 Tcm" = "IL7R", "CD8 Tcm" = "CD28", "CD8 Tcm" = "GPR183"
)

marker_vector_to_sets <- function(marker_vec) {
  marker_vec <- marker_vec[!is.na(names(marker_vec)) & names(marker_vec) != ""]
  split(unname(marker_vec), names(marker_vec))
}

combined_annotation_marker_sets <- c(
  marker_vector_to_sets(major_marker_ls),
  marker_vector_to_sets(T_major_marker_ls),
  marker_vector_to_sets(CD4_major_marker_ls),
  marker_vector_to_sets(CD8_major_marker_ls)
)

annotation_marker_sets <- split(
  seq_along(combined_annotation_marker_sets),
  names(combined_annotation_marker_sets)
) %>%
  imap(function(idx, nm) {
    unique(unlist(combined_annotation_marker_sets[idx], use.names = FALSE))
  })

celltype_marker_sets <- list(
  "CD8 Tem" = c("CD8A", "CD8B", "GZMK", "GZMA", "KLRG1", "NKG7", "PRF1"),
  "CD4 Trm" = c("CD4", "CXCR6", "ITGAE", "CD69", "ZNF683"),
  "CD8 T" = c("CD8A", "CD8B", "CD3D", "CD3E"),
  "CD8 Tnaive" = c("CD8A", "SELL", "LEF1", "TCF7", "CCR7"),
  "gdT" = c("TRDC", "TRGC1", "TRGC2"),
  "CD4 Treg" = c("CD4", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cycling T" = c("MKI67", "TOP2A", "STMN1"),
  "CD4 Tfh" = c("CD4", "CXCR5", "ICOS", "PDCD1", "BATF"),
  "CD8 Tex" = c("CD8A", "HAVCR2", "LAG3", "TIGIT", "TOX", "PDCD1")
)

celltype_marker_ls <- unlist(
  imap(celltype_marker_sets, function(genes, nm) {
    stats::setNames(genes, rep(nm, length(genes)))
  })
)

t_values <- c("T cell", "T cells")
All_T <- subset(all_data, subset = annotation %in% t_values)
message("T cells retained: ", ncol(All_T))
if (ncol(All_T) == 0) {
  stop("No cells found with annotation in: ", paste(t_values, collapse = ", "), call. = FALSE)
}

DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
vars_to_regress <- intersect("nCount_RNA", colnames(All_T[[]]))
All_T <- ScaleData(All_T, vars.to.regress = vars_to_regress)
All_T <- RunPCA(All_T, npcs = 50)

pca_dim_selection <- select_pca_dims(All_T, max_dims = 50, min_dims = 10, max_selected = 30)
dims <- pca_dim_selection$elbow_dim  # selected_dim
write.csv(pca_dim_selection$stats, file.path(qc_dir, "pca_dim_stats.csv"), row.names = FALSE)
message(
  "Selected PCA dims: ", dims,
  " (elbow=", pca_dim_selection$elbow_dim,
  ", variance_rule=", pca_dim_selection$variance_dim, ")"
)

save_pdf(
  file.path(qc_dir, "pca_elbow.pdf"),
  print(
    ElbowPlot(All_T, ndims = 50) +
      geom_vline(xintercept = dims, linetype = "dashed", color = "red") +
      ggtitle(str_glue("PCA elbow; selected dims = {dims}"))
  ),
  width = 8,
  height = 6
)

if (.Platform$OS.type == "windows") {
  future::plan(future::sequential)
} else {
  future::plan(future::multisession, workers = min(10, future::availableCores()))
}

All_T <- RunUMAP(All_T, dims = 1:dims, return.model = TRUE)
save_pdf(
  file.path(cluster_dir, sprintf("Dimplot_dim%d_sample.pdf", dims)),
  print(DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()),
  width = 8,
  height = 6
)

All_T <- harmony::RunHarmony(
  object = All_T,
  dims.use = 1:dims,
  group.by.vars = "sample",
  reduction.save = "harmony",
  verbose = TRUE
)

All_T <- RunUMAP(
  All_T,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key = "harmonyUMAP_",
  return.model = TRUE
)

save_pdf(
  file.path(cluster_dir, sprintf("Dimplot_dim%d_sample_harmony.pdf", dims)),
  print(DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()),
  width = 8,
  height = 6
)

All_T <- FindNeighbors(All_T, reduction = "harmony", dims = 1:dims, k.param = 20)
All_T <- FindClusters(All_T, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0))

cluster_cols <- c("RNA_snn_res.0.1", "RNA_snn_res.0.2", "RNA_snn_res.0.4",
                  "RNA_snn_res.0.6", "RNA_snn_res.0.8", "RNA_snn_res.1")

pdf(file.path(cluster_dir, sprintf("Dimplot_dim%d_clusters.pdf", dims)), width = 8, height = 6)
for (cluster_col in cluster_cols) {
  print(DimPlot(All_T, group.by = cluster_col, label = TRUE, reduction = "harmonyumap") + NoLegend())
}
dev.off()

pdf(file.path(marker_dir, sprintf("dotplot_dim%d_multi_res.pdf", dims)), width = 35, height = 8)
for (cluster_col in cluster_cols) {
  print(safe_dotplot(All_T, cluster_col, major_marker_ls, paste(cluster_col, "major markers")))
  print(safe_dotplot(All_T, cluster_col, T_major_marker_ls, paste(cluster_col, "T markers")))
}
dev.off()

anno_resolution <- "RNA_snn_res.0.6"
Idents(All_T) <- anno_resolution

message("Finding T cell cluster markers at ", anno_resolution)
T_markers <- FindAllMarkers(
  join_layers_if_available(All_T),
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE
)

openxlsx::write.xlsx(T_markers, file.path(marker_dir, "T_cluster_diff_gene.xlsx"))
write.csv(T_markers, file.path(marker_dir, "T_cluster_diff_gene.csv"), row.names = FALSE)

T_top5_markers <- T_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  deframe() %>%
  imap(function(genes, cluster) stats::setNames(genes, rep(paste0("cluster_", cluster), length(genes)))) %>%
  unlist() %>%
  .[!duplicated(.)]

pdf(file.path(marker_dir, sprintf("dotplot_dim%d_major_marker_res0.6.pdf", dims)), width = 35, height = 8)
print(safe_dotplot(All_T, anno_resolution, major_marker_ls, "major markers"))
print(safe_dotplot(All_T, anno_resolution, T_major_marker_ls, "T markers"))
print(safe_dotplot(All_T, anno_resolution, CD4_major_marker_ls, "CD4 markers"))
print(safe_dotplot(All_T, anno_resolution, CD8_major_marker_ls, "CD8 markers"))
if (length(T_top5_markers) > 0) {
  print(safe_dotplot(All_T, anno_resolution, T_top5_markers, "top cluster markers"))
}
dev.off()

avg_exp <- AverageExpression(
  All_T,
  assays = DefaultAssay(All_T),
  group.by = anno_resolution,
  slot = "data",
  verbose = FALSE
)[[DefaultAssay(All_T)]]

score_cluster <- function(genes) {
  genes_use <- intersect(genes, rownames(avg_exp))
  if (length(genes_use) == 0) {
    return(rep(0, ncol(avg_exp)))
  }
  colMeans(as.matrix(avg_exp[genes_use, , drop = FALSE]))
}

major_marker_sets <- marker_vector_to_sets(major_marker_ls)
major_score_mat <- vapply(major_marker_sets, score_cluster, numeric(ncol(avg_exp)))
rownames(major_score_mat) <- colnames(avg_exp)
major_score_mat[is.na(major_score_mat)] <- 0

score_mat <- vapply(annotation_marker_sets, score_cluster, numeric(ncol(avg_exp)))
rownames(score_mat) <- colnames(avg_exp)
score_mat[is.na(score_mat)] <- 0

T_subtype_candidates <- intersect(
  c(
    "CD4 naive", "CD4 Tcm", "CD4 Tem", "CD4 Tfh", "CD4 Th1", "CD4 Th17", "CD4 Treg",
    "CD8 MAIT", "CD8 Tcm", "CD8 Tem", "CD8 Tex", "CD8 Tnaive", "CD8 Trm",
    "Cytotoxic", "Effector", "Exhausted", "gdT", "Proliferating", "Proliferation",
    "Tfh", "TRM"
  ),
  colnames(score_mat)
)

if (length(T_subtype_candidates) == 0) {
  stop("No T subtype marker sets were available for annotation.", call. = FALSE)
}

subtype_score_mat <- score_mat[, T_subtype_candidates, drop = FALSE]
write.csv(
  tibble(T_subtype_candidate = T_subtype_candidates),
  file.path(anno_dir, "T_subtype_annotation_candidates.csv"),
  row.names = FALSE
)

major_annotation <- apply(major_score_mat, 1, function(x) {
  if (all(x == 0)) {
    return("Unresolved")
  }
  names(which.max(x))
})

subtype_annotation <- apply(subtype_score_mat, 1, function(x) {
  if (all(x == 0)) {
    return("Unresolved T")
  }
  names(which.max(x))
})

cluster_annotation <- stats::setNames(
  ifelse(
    major_annotation == "T cell",
    subtype_annotation[names(major_annotation)],
    major_annotation
  ),
  names(major_annotation)
)

cluster_annotation_df <- tibble(
  cluster = names(cluster_annotation),
  major_celltype = unname(major_annotation[names(cluster_annotation)]),
  refined_T_subtype = unname(subtype_annotation[names(cluster_annotation)]),
  T_subtype = unname(cluster_annotation)
) %>%
  
  bind_cols(
    as_tibble(major_score_mat, rownames = NULL, .name_repair = function(x) paste0("major_score_", x)),
    as_tibble(score_mat, rownames = NULL, .name_repair = function(x) paste0("score_", x))
  )

write.csv(cluster_annotation_df, file.path(anno_dir, "T_cluster_annotation_scores.csv"), row.names = FALSE)


cell_cluster <- stats::setNames(cluster_annotation_df$T_subtype, cluster_annotation_df$cluster)
All_T <- RenameIdents(All_T, cell_cluster)
All_T$T_subtype <- Idents(All_T)

pdf(file.path(anno_dir, sprintf("dotplot_dim%d_T_subtype_marker.pdf", dims)), width = 20, height = 8)
print(safe_dotplot(All_T, "T_subtype", celltype_marker_ls, "T subtype markers", dot.scale = 10) + theme_90angle)
dev.off()

save_pdf(
  file.path(anno_dir, sprintf("Dimplot_dim%d_T_subtype.pdf", dims)),
  print(DimPlot(All_T, group.by = "T_subtype", label = TRUE, reduction = "harmonyumap", repel = TRUE) + NoLegend()),
  width = 8,
  height = 6
)

if ("group" %in% colnames(All_T[[]])) {
  save_pdf(
    file.path(anno_dir, sprintf("Dimplot_dim%d_T_subtype_split_group.pdf", dims)),
    print(DimPlot(All_T, group.by = "T_subtype", split.by = "group", label = TRUE,
                  reduction = "harmonyumap", repel = TRUE) + NoLegend()),
    width = 20,
    height = 7
  )
}

subtype_meta <- All_T[["T_subtype"]] %>%
  rownames_to_column("cells")

all_data_meta <- all_data[["annotation"]] %>%
  rownames_to_column("cells") %>%
  left_join(subtype_meta, by = "cells") %>%
  mutate(
    T_annotation = ifelse(
      as.character(annotation) %in% t_values & !is.na(T_subtype),
      as.character(T_subtype),
      as.character(annotation)
    )
  ) %>%
  select(cells, T_annotation) %>%
  column_to_rownames("cells")

all_data <- AddMetaData(all_data, all_data_meta["T_annotation"])

save_pdf(
  file.path(anno_dir, "Dimplot_all_data_T_annotation.pdf"),
  print(DimPlot(all_data, group.by = "T_annotation", label = TRUE, repel = TRUE)),
  width = 10,
  height = 6
)

write.csv(
  All_T[[]] %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, sample, any_of("group"), all_of(anno_resolution), T_subtype),
  file.path(anno_dir, "T_cell_subtype_annotation.csv"),
  row.names = FALSE
)

write.csv(
  all_data[[]] %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, annotation, T_annotation),
  file.path(anno_dir, "all_data_T_annotation.csv"),
  row.names = FALSE
)

saveRDS(All_T, file.path(out_dir, "All_T_subtype.rds"))
saveRDS(all_data, file.path(out_dir, "all_data2_TLS_T_annotation.rds"))

message("Done. Results written to: ", normalizePath(out_dir, winslash = "/"))
suppressPackageStartupMessages({
  required_pkgs <- c(
    "Seurat", "tidyverse", "harmony", "future", "openxlsx",
    "tibble", "data.table", "patchwork"
  )
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nNo packages were installed. Please install/activate them first if needed.",
      call. = FALSE
    )
  }
  
  library(Seurat)
  library(tidyverse)
  library(harmony)
  library(future)
})

set.seed(1234)
options(future.globals.maxSize = 10 * 1024^3)

input_rds <- "all_data2_TLS.rds"
out_dir <- "T_subtype_codex_result"
qc_dir <- file.path(out_dir, "01_qc")
cluster_dir <- file.path(out_dir, "02_cluster")
marker_dir <- file.path(out_dir, "03_marker")
anno_dir <- file.path(out_dir, "04_annotation")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(anno_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading: ", input_rds)
all_data <- readRDS(input_rds)

if (!inherits(all_data, "Seurat")) {
  stop("The input RDS is not a Seurat object: ", input_rds, call. = FALSE)
}

meta_cols <- colnames(all_data[[]])
if (!"annotation" %in% meta_cols) {
  stop("Metadata column `annotation` was not found in all_data2_TLS.rds.", call. = FALSE)
}
if (!"sample" %in% meta_cols) {
  stop("Metadata column `sample` was not found; it is required for Harmony batch correction.", call. = FALSE)
}

message("annotation table:")
print(table(all_data$annotation, useNA = "ifany"))

theme_90angle <- theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

join_layers_if_available <- function(object) {
  if (!exists("JoinLayers", mode = "function")) {
    return(object)
  }
  tryCatch(
    JoinLayers(object),
    error = function(e) {
      message("JoinLayers skipped: ", conditionMessage(e))
      object
    }
  )
}

select_pca_dims <- function(object, reduction = "pca", max_dims = 50, min_dims = 10, max_selected = 30) {
  stdev <- object[[reduction]]@stdev
  max_dims <- min(max_dims, length(stdev))
  stdev <- stdev[seq_len(max_dims)]
  variance <- stdev^2
  pct_variance <- variance / sum(variance) * 100
  cumulative_variance <- cumsum(pct_variance)
  
  x <- seq_along(stdev)
  y <- stdev
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  line_y <- y_scaled[1] + (y_scaled[length(y_scaled)] - y_scaled[1]) * x_scaled
  elbow_dim <- which.max(abs(y_scaled - line_y))
  
  variance_dim <- which(cumulative_variance >= 90 & pct_variance < 5)[1]
  if (is.na(variance_dim)) {
    variance_dim <- elbow_dim
  }
  
  selected_dim <- max(min_dims, min(max_selected, max(elbow_dim, variance_dim)))
  
  list(
    selected_dim = selected_dim,
    elbow_dim = elbow_dim,
    variance_dim = variance_dim,
    stats = tibble(
      PC = x,
      stdev = stdev,
      pct_variance = pct_variance,
      cumulative_variance = cumulative_variance,
      selected = x == selected_dim,
      elbow = x == elbow_dim,
      variance_rule = x == variance_dim
    )
  )
}

filter_features <- function(object, marker_vec) {
  marker_vec[marker_vec %in% rownames(object)]
}

safe_dotplot <- function(object, group.by, features, title = NULL, dot.scale = 6) {
  features_use <- filter_features(object, features)
  if (length(features_use) == 0) {
    warning("No requested features found for group.by=", group.by)
    return(ggplot() + theme_void() + ggtitle(paste("No markers found:", title %||% group.by)))
  }
  DotPlot(object, group.by = group.by, features = features_use, dot.scale = dot.scale) +
    RotatedAxis() +
    ggtitle(title %||% group.by)
}

save_pdf <- function(path, plot_expr, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(plot_expr)
}

major_marker_ls <- c(
  "Epithelial" = "EPCAM", "Epithelial" = "KRT5", "Epithelial" = "KRT14",
  "Endothelial" = "PECAM1", "Endothelial" = "VWF",
  "Fibroblast" = "PDGFRA", "Fibroblast" = "CCDC80", "Fibroblast" = "SFRP4",
  "Fibroblast" = "COL1A1", "Fibroblast" = "COL1A2", "Fibroblast" = "DCN",
  "B cell" = "CD19", "B cell" = "CD79A", "B cell" = "MS4A1",
  "Plasma" = "TNFRSF17", "Plasma" = "JCHAIN", "Plasma" = "MZB1",
  "T cell" = "CD3D", "T cell" = "CD3E", "T cell" = "CD3G", "T cell" = "TRAC",
  "NK" = "GNLY", "NK" = "NKG7", "NK" = "NCAM1",
  "Proliferating" = "MKI67", "Proliferating" = "TOP2A",
  "Monocyte" = "CD14", "Monocyte" = "FCGR3A",
  "Macrophage" = "CD68", "Macrophage" = "CD163", "Macrophage" = "C1QA",
  "Dendritic cell" = "CD83", "Dendritic cell" = "CLEC9A", "Dendritic cell" = "CD1C",
  "Mast cell" = "CPA3", "Mast cell" = "TPSAB1", "Mast cell" = "KIT"
)

T_major_marker_ls <- c(
  "CD8" = "CD8A", "CD8" = "CD8B", "CD8" = "GZMB", "CD8" = "GZMA",
  "CD8" = "GZMK", "CD8" = "GZMH", "CD8" = "CCL5",
  "Naive" = "LEF1", "Naive" = "TCF7", "Naive" = "SELL", "Naive" = "CCR7",
  "Memory" = "GPR183", "Memory" = "IL7R", "Memory" = "FOS", "Memory" = "JUN",
  "TRM" = "CXCR6", "TRM" = "ITGAE", "TRM" = "CD69", "TRM" = "ZNF683",
  "Effector" = "FGFBP2", "Effector" = "TYROBP", "Effector" = "CX3CR1",
  "Effector" = "FCGR3A", "Effector" = "GNLY",
  "Cytotoxic" = "XCL1", "Cytotoxic" = "IFNG", "Cytotoxic" = "CCL4",
  "Exhausted" = "CXCL13", "Exhausted" = "LAYN", "Exhausted" = "LAG3",
  "Exhausted" = "HAVCR2", "Exhausted" = "TNFRSF9", "Exhausted" = "CTLA4",
  "Exhausted" = "TIGIT", "Exhausted" = "TOX", "Exhausted" = "PDCD1",
  "CD4" = "FOXP3", "CD4" = "CD4", "CD4" = "IL2RA",
  "Tfh" = "CXCR5", "Tfh" = "ICOS", "Tfh" = "BATF", "Tfh" = "PDCD1",
  "gdT" = "TRDC",
  "Proliferation" = "MKI67", "Proliferation" = "TOP2A", "Proliferation" = "STMN1"
)

CD4_major_marker_ls <- c(
  "T cell" = "CD3D", "T cell" = "CD3E", "T cell" = "CD3G", "T cell" = "TRAC",
  "CD4" = "CD4",
  "CD4 Treg" = "TNFRSF4", "CD4 Treg" = "TNFRSF9", "CD4 Treg" = "FOXP3",
  "CD4 Treg" = "IL2RA", "CD4 Treg" = "LAYN", "CD4 Treg" = "CTLA4",
  "CD4 Treg" = "TIGIT",
  "CD4 naive" = "LEF1", "CD4 naive" = "CCR7", "CD4 naive" = "TCF7",
  "CD4 Tem" = "GZMK", "CD4 Tem" = "GZMA", "CD4 Tem" = "NKG7",
  "CD4 Th1" = "CXCL13", "CD4 Th1" = "BHLHE40", "CD4 Th1" = "IFNG",
  "CD4 Tfh" = "TOX2", "CD4 Tfh" = "PDCD1", "CD4 Tfh" = "CXCR5",
  "CD4 Th17" = "IL17A", "CD4 Th17" = "IL17F", "CD4 Th17" = "IL23R",
  "CD4 Tcm" = "TCF7", "CD4 Tcm" = "CCR7", "CD4 Tcm" = "S1PR1"
)

CD8_major_marker_ls <- c(
  "CD8" = "CD8A", "CD8" = "CD8B",
  "CD8 Tex" = "HAVCR2", "CD8 Tex" = "CTLA4", "CD8 Tex" = "LAYN",
  "CD8 Tex" = "LAG3", "CD8 Tex" = "TIGIT",
  "CD8 Tem" = "GZMK", "CD8 Tem" = "GZMA", "CD8 Tem" = "IFNG",
  "CD8 Tem" = "KLRG1", "CD8 Tem" = "GZMB", "CD8 Tem" = "PRF1",
  "CD8 Tem" = "NKG7",
  "CD8 Tnaive" = "SELL", "CD8 Tnaive" = "LEF1", "CD8 Tnaive" = "TCF7",
  "CD8 Tnaive" = "CCR7",
  "CD8 MAIT" = "NCR3", "CD8 MAIT" = "SLC4A10", "CD8 MAIT" = "KLRB1",
  "CD8 Trm" = "MYADM", "CD8 Trm" = "CD69", "CD8 Trm" = "ITGAE",
  "CD8 Trm" = "CXCR6", "CD8 Trm" = "ITGA1",
  "CD8 Tcm" = "IL7R", "CD8 Tcm" = "CD28", "CD8 Tcm" = "GPR183"
)

marker_vector_to_sets <- function(marker_vec) {
  marker_vec <- marker_vec[!is.na(names(marker_vec)) & names(marker_vec) != ""]
  split(unname(marker_vec), names(marker_vec))
}

combined_annotation_marker_sets <- c(
  marker_vector_to_sets(major_marker_ls),
  marker_vector_to_sets(T_major_marker_ls),
  marker_vector_to_sets(CD4_major_marker_ls),
  marker_vector_to_sets(CD8_major_marker_ls)
)

annotation_marker_sets <- split(
  seq_along(combined_annotation_marker_sets),
  names(combined_annotation_marker_sets)
) %>%
  imap(function(idx, nm) {
    unique(unlist(combined_annotation_marker_sets[idx], use.names = FALSE))
  })

celltype_marker_sets <- list(
  "CD8 Tem" = c("CD8A", "CD8B", "GZMK", "GZMA", "KLRG1", "NKG7", "PRF1"),
  "CD4 Trm" = c("CD4", "CXCR6", "ITGAE", "CD69", "ZNF683"),
  "CD8 T" = c("CD8A", "CD8B", "CD3D", "CD3E"),
  "CD8 Tnaive" = c("CD8A", "SELL", "LEF1", "TCF7", "CCR7"),
  "gdT" = c("TRDC", "TRGC1", "TRGC2"),
  "CD4 Treg" = c("CD4", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cycling T" = c("MKI67", "TOP2A", "STMN1"),
  "CD4 Tfh" = c("CD4", "CXCR5", "ICOS", "PDCD1", "BATF"),
  "CD8 Tex" = c("CD8A", "HAVCR2", "LAG3", "TIGIT", "TOX", "PDCD1")
)

celltype_marker_ls <- unlist(
  imap(celltype_marker_sets, function(genes, nm) {
    stats::setNames(genes, rep(nm, length(genes)))
  })
)

t_values <- c("T cell", "T cells")
All_T <- subset(all_data, subset = annotation %in% t_values)
message("T cells retained: ", ncol(All_T))
if (ncol(All_T) == 0) {
  stop("No cells found with annotation in: ", paste(t_values, collapse = ", "), call. = FALSE)
}

DefaultAssay(All_T) <- "RNA"
All_T <- FindVariableFeatures(All_T)
vars_to_regress <- intersect("nCount_RNA", colnames(All_T[[]]))
All_T <- ScaleData(All_T, vars.to.regress = vars_to_regress)
All_T <- RunPCA(All_T, npcs = 50)

pca_dim_selection <- select_pca_dims(All_T, max_dims = 50, min_dims = 10, max_selected = 30)
dims <- pca_dim_selection$selected_dim
write.csv(pca_dim_selection$stats, file.path(qc_dir, "pca_dim_stats.csv"), row.names = FALSE)
message(
  "Selected PCA dims: ", dims,
  " (elbow=", pca_dim_selection$elbow_dim,
  ", variance_rule=", pca_dim_selection$variance_dim, ")"
)

save_pdf(
  file.path(qc_dir, "pca_elbow.pdf"),
  print(
    ElbowPlot(All_T, ndims = 50) +
      geom_vline(xintercept = dims, linetype = "dashed", color = "red") +
      ggtitle(str_glue("PCA elbow; selected dims = {dims}"))
  ),
  width = 8,
  height = 6
)

if (.Platform$OS.type == "windows") {
  future::plan(future::sequential)
} else {
  future::plan(future::multisession, workers = min(10, future::availableCores()))
}

All_T <- RunUMAP(All_T, dims = 1:dims, return.model = TRUE)
save_pdf(
  file.path(cluster_dir, sprintf("Dimplot_dim%d_sample.pdf", dims)),
  print(DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()),
  width = 8,
  height = 6
)

All_T <- harmony::RunHarmony(
  object = All_T,
  dims.use = 1:dims,
  group.by.vars = "sample",
  reduction.save = "harmony",
  verbose = TRUE
)

All_T <- RunUMAP(
  All_T,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key = "harmonyUMAP_",
  return.model = TRUE
)

save_pdf(
  file.path(cluster_dir, sprintf("Dimplot_dim%d_sample_harmony.pdf", dims)),
  print(DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()),
  width = 8,
  height = 6
)

All_T <- FindNeighbors(All_T, reduction = "harmony", dims = 1:dims, k.param = 20)
All_T <- FindClusters(All_T, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0))

cluster_cols <- c("RNA_snn_res.0.1", "RNA_snn_res.0.2", "RNA_snn_res.0.4",
                  "RNA_snn_res.0.6", "RNA_snn_res.0.8", "RNA_snn_res.1")

pdf(file.path(cluster_dir, sprintf("Dimplot_dim%d_clusters.pdf", dims)), width = 8, height = 6)
for (cluster_col in cluster_cols) {
  print(DimPlot(All_T, group.by = cluster_col, label = TRUE, reduction = "harmonyumap") + NoLegend())
}
dev.off()

pdf(file.path(marker_dir, sprintf("dotplot_dim%d_multi_res.pdf", dims)), width = 35, height = 8)
for (cluster_col in cluster_cols) {
  print(safe_dotplot(All_T, cluster_col, major_marker_ls, paste(cluster_col, "major markers")))
  print(safe_dotplot(All_T, cluster_col, T_major_marker_ls, paste(cluster_col, "T markers")))
}
dev.off()

anno_resolution <- "RNA_snn_res.0.6"
Idents(All_T) <- anno_resolution

message("Finding T cell cluster markers at ", anno_resolution)
T_markers <- FindAllMarkers(
  join_layers_if_available(All_T),
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE
)

openxlsx::write.xlsx(T_markers, file.path(marker_dir, "T_cluster_diff_gene.xlsx"))
write.csv(T_markers, file.path(marker_dir, "T_cluster_diff_gene.csv"), row.names = FALSE)

T_top5_markers <- T_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  deframe() %>%
  imap(function(genes, cluster) stats::setNames(genes, rep(paste0("cluster_", cluster), length(genes)))) %>%
  unlist() %>%
  .[!duplicated(.)]

pdf(file.path(marker_dir, sprintf("dotplot_dim%d_major_marker_res0.6.pdf", dims)), width = 35, height = 8)
print(safe_dotplot(All_T, anno_resolution, major_marker_ls, "major markers"))
print(safe_dotplot(All_T, anno_resolution, T_major_marker_ls, "T markers"))
print(safe_dotplot(All_T, anno_resolution, CD4_major_marker_ls, "CD4 markers"))
print(safe_dotplot(All_T, anno_resolution, CD8_major_marker_ls, "CD8 markers"))
if (length(T_top5_markers) > 0) {
  print(safe_dotplot(All_T, anno_resolution, T_top5_markers, "top cluster markers"))
}
dev.off()

avg_exp <- AverageExpression(
  All_T,
  assays = DefaultAssay(All_T),
  group.by = anno_resolution,
  slot = "data",
  verbose = FALSE
)[[DefaultAssay(All_T)]]

score_cluster <- function(genes) {
  genes_use <- intersect(genes, rownames(avg_exp))
  if (length(genes_use) == 0) {
    return(rep(NA_real_, ncol(avg_exp)))
  }
  colMeans(as.matrix(avg_exp[genes_use, , drop = FALSE]))
}

score_mat <- vapply(annotation_marker_sets, score_cluster, numeric(ncol(avg_exp)))
rownames(score_mat) <- colnames(avg_exp)

cluster_annotation <- apply(score_mat, 1, function(x) {
  if (all(is.na(x))) {
    return("T cell")
  }
  names(which.max(x))
})

cluster_annotation_df <- tibble(
  cluster = names(cluster_annotation),
  T_subtype = unname(cluster_annotation)
) %>%
  bind_cols(as_tibble(score_mat, rownames = NULL))

write.csv(cluster_annotation_df, file.path(anno_dir, "T_cluster_annotation_scores.csv"), row.names = FALSE)

cell_cluster <- stats::setNames(cluster_annotation_df$T_subtype, cluster_annotation_df$cluster)
All_T <- RenameIdents(All_T, cell_cluster)
All_T$T_subtype <- Idents(All_T)

pdf(file.path(anno_dir, sprintf("dotplot_dim%d_T_subtype_marker.pdf", dims)), width = 20, height = 8)
print(safe_dotplot(All_T, "T_subtype", celltype_marker_ls, "T subtype markers", dot.scale = 10) + theme_90angle)
dev.off()

save_pdf(
  file.path(anno_dir, sprintf("Dimplot_dim%d_T_subtype.pdf", dims)),
  print(DimPlot(All_T, group.by = "T_subtype", label = TRUE, reduction = "harmonyumap", repel = TRUE) + NoLegend()),
  width = 8,
  height = 6
)

if ("group" %in% colnames(All_T[[]])) {
  save_pdf(
    file.path(anno_dir, sprintf("Dimplot_dim%d_T_subtype_split_group.pdf", dims)),
    print(DimPlot(All_T, group.by = "T_subtype", split.by = "group", label = TRUE,
                  reduction = "harmonyumap", repel = TRUE) + NoLegend()),
    width = 20,
    height = 7
  )
}

subtype_meta <- All_T[["T_subtype"]] %>%
  rownames_to_column("cells")

all_data_meta <- all_data[["annotation"]] %>%
  rownames_to_column("cells") %>%
  left_join(subtype_meta, by = "cells") %>%
  mutate(
    T_annotation = ifelse(
      as.character(annotation) %in% t_values & !is.na(T_subtype),
      as.character(T_subtype),
      as.character(annotation)
    )
  ) %>%
  select(cells, T_annotation) %>%
  column_to_rownames("cells")

all_data <- AddMetaData(all_data, all_data_meta["T_annotation"])

save_pdf(
  file.path(anno_dir, "Dimplot_all_data_T_annotation.pdf"),
  print(DimPlot(all_data, group.by = "T_annotation", label = TRUE, repel = TRUE)),
  width = 10,
  height = 6
)

write.csv(
  All_T[[]] %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, sample, any_of("group"), all_of(anno_resolution), T_subtype),
  file.path(anno_dir, "T_cell_subtype_annotation.csv"),
  row.names = FALSE
)

write.csv(
  all_data[[]] %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, annotation, T_annotation),
  file.path(anno_dir, "all_data_T_annotation.csv"),
  row.names = FALSE
)

saveRDS(All_T, file.path(out_dir, "All_T_subtype.rds"))
saveRDS(all_data, file.path(out_dir, "all_data2_TLS_T_annotation.rds"))

message("Done. Results written to: ", normalizePath(out_dir, winslash = "/"))
