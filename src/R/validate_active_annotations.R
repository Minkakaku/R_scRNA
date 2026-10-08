# 独立会话核对当前配置、完整QS、绘图缓存与逐细胞注释的连接；仅输出校验日志。
args <- commandArgs(trailingOnly = TRUE)
source('src/R/project_config.R')
load_project_config(args[1])
suppressPackageStartupMessages(library(Seurat))
if (.Platform$OS.type == 'windows') Sys.setlocale('LC_CTYPE', 'English_United States.UTF-8')
s <- CONFIG$current_results
folder <- file.path(PROJECT_DIR, s$annotation_evidence)
mapping <- read.delim(file.path(PROJECT_DIR, s$annotation_mapping), encoding = 'UTF-8')
manifest <- read.delim(file.path(PROJECT_DIR, s$subgroup_manifest))
cells <- readRDS(file.path(folder, 'current_cell_annotations.rds'))
stopifnot(!anyDuplicated(cells$cell_id), sum(manifest$cells) == nrow(cells))
sc <- qs::qread(file.path(PROJECT_DIR, s$complete_object), nthreads = 4)
md <- sc[[]]
ix <- match(colnames(sc), cells$cell_id)
key <- paste(md$analysis_container, md$current_subgroup_cluster)
mi <- match(key, paste(mapping$subgroup, mapping$subgroup_cluster))
stopifnot(!anyNA(ix), !anyNA(mi), ncol(sc) == nrow(cells),
  identical(as.character(Idents(sc)), as.character(md$subtype_annotation)),
  identical(as.character(md$subtype_annotation), mapping$annotation[mi]),
  identical(md$annotation_PMIDs, mapping$PMIDs[mi]),
  identical(as.character(md$neutrophil_function), mapping$neutrophil_function[mi]),
  identical(md$annotation_limitations, mapping$limitations_CN[mi]),
  all(md$annotation_stage == s$annotation_stage))
columns <- c('subtype_annotation', 'biological_lineage', 'annotation_confidence',
  'annotation_PMIDs', 'annotation_limitations', 'literature_cluster_label',
  'annotation_stage', 'functional_annotation_stage', 'annotation_name_type', 'functional_interpretation', 'neutrophil_function')
for (column in columns) stopifnot(identical(as.character(md[[column]]), as.character(cells[[column]][ix])))
for (i in seq_len(nrow(manifest))) {
  group <- manifest$subgroup[i]
  cached <- readRDS(file.path(folder, group, 'annotation_plot_data.rds'))
  pm <- cached$plot_object[[]]
  ci <- match(rownames(pm), rownames(md))
  stopifnot(!anyNA(ci), nrow(pm) == manifest$cells[i],
    startsWith(manifest$checkpoint[i], paste0(s$subgroup_checkpoints, '/')))
  for (column in columns) stopifnot(identical(as.character(pm[[column]]), as.character(md[[column]][ci])))
  umap <- Embeddings(cached$plot_object, 'umap')
  stopifnot(max(abs(umap[, 1] - md$subgroup_UMAP_1[ci])) < 1e-12,
    max(abs(umap[, 2] - md$subgroup_UMAP_2[ci])) < 1e-12)
}
report <- list(cells = ncol(sc), genes = nrow(sc), samples = length(unique(md$sample)),
  clusters = nrow(mapping), active_object = s$complete_object, stage = s$annotation_stage,
  full_metadata_and_plot_caches_match_by_cell_ID = TRUE,
  mapping_PMIDs_limitations_match = TRUE, active_config_consistent = TRUE,
  current_annotations = sort(unique(as.character(md$subtype_annotation))))
jsonlite::write_json(report, file.path(PROJECT_DIR, args[2]), pretty = TRUE, auto_unbox = TRUE)
cat('ACTIVE_ANNOTATION_STATE_VALIDATED\n')
