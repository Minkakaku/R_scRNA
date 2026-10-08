# 人类 HGNC 映射的可选工具；显式传入参考与设置，不自动运行或改写源对象。

.gene_id_mapping_cache <- new.env(parent = emptyenv())

# 读取用户指定的 HGNC Approved Ensembl 到 symbol 映射。
read_hgnc_ensembl_symbol_map <- function(settings) {
  path <- settings$mapping_path
  if (!file.exists(path)) stop("HGNC mapping file not found: ", path, call. = FALSE)
  cache_key <- paste(normalizePath(path, winslash = "/", mustWork = TRUE), tools::md5sum(path), sep = ":")
  if (exists(cache_key, envir = .gene_id_mapping_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_id_mapping_cache, inherits = FALSE))
  }

  hgnc <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  required <- c("symbol", "status", "ensembl_gene_id")
  missing <- setdiff(required, colnames(hgnc))
  if (length(missing) > 0) {
    stop("HGNC mapping is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  hgnc <- hgnc[
    hgnc$status == "Approved" & nzchar(hgnc$symbol) & nzchar(hgnc$ensembl_gene_id),
    c("ensembl_gene_id", "symbol"),
    drop = FALSE
  ]
  hgnc$ensembl_gene_id <- sub("\\..*$", "", hgnc$ensembl_gene_id)
  hgnc <- hgnc[!duplicated(hgnc$ensembl_gene_id), , drop = FALSE]
  mapping <- stats::setNames(hgnc$symbol, hgnc$ensembl_gene_id)
  assign(cache_key, mapping, envir = .gene_id_mapping_cache)
  mapping
}

# 按指定映射合并同名基因的原始 counts，保留细胞顺序并返回审计信息。
harmonize_seurat_gene_ids <- function(obj, settings, sample_name, assay = "RNA") {
  if (is.null(settings) || !isTRUE(settings$enabled)) {
    return(list(object = obj, audit = NULL))
  }
  if (!assay %in% names(obj@assays)) {
    stop("Assay not found for gene ID harmonization: ", assay, call. = FALSE)
  }

  counts <- LayerData(obj, assay = assay, layer = "counts")
  expressed <- Matrix::rowSums(counts) > 0
  counts <- counts[expressed, , drop = FALSE]
  input_ids <- rownames(counts)
  stable_ids <- sub("\\..*$", "", input_ids)
  is_ensembl <- grepl("^ENSG[0-9]+$", stable_ids)
  mapping <- read_hgnc_ensembl_symbol_map(settings)
  mapped_symbols <- unname(mapping[stable_ids])
  mapped <- is_ensembl & !is.na(mapped_symbols) & nzchar(mapped_symbols)

  output_ids <- input_ids
  output_ids[mapped] <- mapped_symbols[mapped]
  if (!isTRUE(settings$keep_unmapped)) {
    keep <- !is_ensembl | mapped
    counts <- counts[keep, , drop = FALSE]
    input_ids <- input_ids[keep]
    stable_ids <- stable_ids[keep]
    is_ensembl <- is_ensembl[keep]
    mapped <- mapped[keep]
    output_ids <- output_ids[keep]
  }

  unique_output_ids <- unique(output_ids)
  collapse <- Matrix::sparseMatrix(
    i = match(output_ids, unique_output_ids),
    j = seq_along(output_ids),
    x = 1,
    dims = c(length(unique_output_ids), length(output_ids))
  )
  collapsed_counts <- collapse %*% counts
  rownames(collapsed_counts) <- unique_output_ids
  colnames(collapsed_counts) <- colnames(counts)
  obj[[assay]] <- SeuratObject::CreateAssay5Object(counts = collapsed_counts)
  DefaultAssay(obj) <- assay

  audit <- data.frame(
    sample = sample_name,
    input_union_features = length(expressed),
    expressed_features = length(output_ids),
    ensembl_features = sum(is_ensembl),
    mapped_ensembl = sum(mapped),
    ensembl_mapping_rate_pct = if (sum(is_ensembl) > 0) round(100 * sum(mapped) / sum(is_ensembl), 2) else NA_real_,
    duplicate_symbols_collapsed = length(output_ids) - length(unique_output_ids),
    output_features = length(unique_output_ids),
    reference_name = settings$reference_name,
    reference_sha256 = settings$reference_sha256,
    stringsAsFactors = FALSE
  )
  obj@misc$pipeline$gene_id_harmonization <- as.list(audit[1, , drop = FALSE])
  list(object = obj, audit = audit)
}

# 把各样本的映射统计写到调用者指定位置。
write_gene_id_harmonization_audit <- function(audit, path) {
  audit <- Filter(Negate(is.null), audit)
  if (length(audit) == 0) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, audit), path, row.names = FALSE, quote = FALSE)
  invisible(path)
}
