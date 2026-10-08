# 用合成对象验证项目提取的改进，不运行真实 QC、整合或聚类。
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
source("src/R/00_runtime.R")
source("src/R/02_sample_metadata.R")
source("src/R/05_normalization.R")
source("src/R/07_clustering.R")
source("src/R/gene_id_harmonization.R")
source("src/R/annotation_validation.R")

check_project_updates <- function() {
  scratch <- tempfile("project-function-check-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  pdf(file.path(scratch, "check.pdf"))
  on.exit(dev.off(), add = TRUE)
  counts <- Matrix::Matrix(matrix(c(3, 2, 1, 2, 4, 1, 1, 3, 5, 2, 1, 3), nrow = 3), sparse = TRUE)
  rownames(counts) <- c("ENSG000001.1", "ENSG000002.1", "GeneC")
  colnames(counts) <- paste0("cell", 1:4)
  obj <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  obj$orig.ident <- rep(c("s1", "s2"), each = 2)
  obj$group <- rep(c("CON", "EXP"), each = 2)
  obj$tissue <- rep(c("site, one", "site two"), each = 2)
  sheet <- file.path(scratch, "samples.csv")
  write_sample_group_sheet(obj, sheet, c(group = "group", tissue = "tissue"))
  observed <- read_sample_sheet(sheet, required = TRUE)
  stopifnot(identical(observed$sample, c("s1", "s2")), observed$tissue[1] == "site, one")
  write.csv(rbind(observed, observed[1, ]), sheet, row.names = FALSE)
  stopifnot(inherits(try(read_sample_sheet(sheet), silent = TRUE), "try-error"))
  write.csv(observed[1, ], sheet, row.names = FALSE)
  SAMPLE_SHEET_PATH <<- sheet
  GROUP_LEVELS <<- c("CON", "EXP")
  stopifnot(inherits(try(assign_group_metadata(obj, required_sheet = TRUE), silent = TRUE), "try-error"))

  first <- obj[, 1:2]
  second <- obj[, 3:4]
  VariableFeatures(first) <- rownames(first)[c(1, 2, 3)]
  VariableFeatures(second) <- rownames(second)[c(2, 1)]
  selected <- select_integration_features_robust(list(first, second), 2, "RNA")
  stopifnot(identical(selected, sort(rownames(obj)[1:2])))
  # 强制原生接口失败，确认两个正式调用链可复用相同回退函数。
  fallback <- new.env(parent = globalenv())
  fallback$SelectIntegrationFeatures <- function(...) stop("synthetic interface failure")
  safe <- select_integration_features_safe
  environment(safe) <- fallback
  stopifnot(identical(safe(list(first, second), 2, "RNA", "check"), selected))

  variance <- data.frame(pc = 1:4, variance = c(5, 3, 1, 1), variance_percent = c(50, 30, 10, 10), cumulative_percent = c(50, 80, 90, 100))
  result <- detect_pca_elbow_from_variance(variance)
  stopifnot(result$recommended_pc %in% 1:4,
            result$selections$pc[result$selections$method == "cumulative_80"] == 2)
  singleton <- detect_pca_elbow_from_variance(data.frame(pc = 1, variance = 1, variance_percent = 100, cumulative_percent = 100))
  stopifnot(singleton$recommended_pc == 1)
  stem <- file.path(scratch, "pca")
  save_pca_elbow_result(result, stem, selected_pc = 2)
  repair_pca_elbow_outputs_from_csv(stem, selected_pc = 2)
  stopifnot(file.exists(paste0(stem, ".png")), file.exists(paste0(stem, ".pdf")))

  mapping <- file.path(scratch, "hgnc.tsv")
  write.table(data.frame(symbol = c("GeneA", "GeneA"), status = "Approved", ensembl_gene_id = c("ENSG000001", "ENSG000002")), mapping, sep = "\t", row.names = FALSE, quote = FALSE)
  mapped <- harmonize_seurat_gene_ids(obj, list(enabled = TRUE, mapping_path = mapping,
    keep_unmapped = TRUE, reference_name = "synthetic", reference_sha256 = "synthetic"), "check")$object
  actual <- LayerData(mapped, assay = "RNA", layer = "counts")
  stopifnot(identical(colnames(actual), colnames(counts)),
            all(actual["GeneA", ] == Matrix::colSums(counts[1:2, ])), all(actual["GeneC", ] == counts[3, ]))
  labels <- c("a", "a", "b", "b")
  stopifnot(adjusted_rand_index(labels, labels) == 1, normalized_mutual_information(labels, labels) == 1,
    cluster_purity(labels, labels) == 1, adjusted_rand_index(rep("a", 4), rep("b", 4)) == 1,
    is.na(cluster_purity(character(), character())))
  cat("PROJECT_FUNCTION_UPDATES_OK\n")
}
check_project_updates()
