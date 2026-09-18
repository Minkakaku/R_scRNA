# 按样本与细胞类型统计细胞数量，输出后续比例分析所需的计数表。
save_proportion_tables <- function(obj, sample_col, celltype_col, output_dir, prefix) {
  counts <- as.data.frame.matrix(
    table(obj[[sample_col, drop = TRUE]], obj[[celltype_col, drop = TRUE]])
  )
  counts <- data.frame(sample = rownames(counts), counts, check.names = FALSE)
  colnames(counts)[1] <- sample_col
  write.csv(
    counts,
    file.path(output_dir, paste0(prefix, "_count.csv")),
    quote = FALSE,
    row.names = FALSE
  )
}

# 输出样本及组别的组成图、计数表和可选桑基图，并逐个处理指定比较。
save_proportion_outputs <- function(
    obj,
    comparisons,
    sample_col,
    group_col,
    celltype_col,
    output_dir,
    prefix,
    table_dir = output_dir,
    figure_dir = output_dir,
    make_sankey = TRUE,
    sankey_min_count = 1L) {
  plot_celltype_distribution(
    obj,
    sample_col = sample_col,
    celltype_col = celltype_col,
    output_file = paste0(prefix, "_by_sample.pdf"),
    output_dir = figure_dir
  )
  plot_celltype_distribution(
    obj,
    sample_col = group_col,
    celltype_col = celltype_col,
    output_file = paste0(prefix, "_by_group.pdf"),
    output_dir = figure_dir
  )
  save_proportion_tables(obj, sample_col, celltype_col, table_dir, paste0(prefix, "_by_sample"))

  if (make_sankey) {
    plot_celltype_sankey(
      obj,
      source_col = group_col,
      target_col = celltype_col,
      path_no_ext = file.path(figure_dir, paste0(prefix, "_sankey_group_to_celltype")),
      min_count = sankey_min_count,
      width = 12,
      height = 8
    )
  }

  for (groups in comparisons) {
    group_values <- as.character(obj[[group_col, drop = TRUE]])
    missing_groups <- setdiff(groups, unique(group_values))
    if (length(missing_groups) > 0) {
      warning("Skip comparison with missing groups: ", paste(groups, collapse = " vs "))
      next
    }

    obj_sub <- subset(obj, cells = colnames(obj)[group_values %in% groups])
    obj_sub[[group_col]] <- factor(obj_sub[[group_col, drop = TRUE]], levels = groups)
    tag <- make_tag(groups)

    plot_celltype_distribution(
      obj_sub,
      sample_col = sample_col,
      celltype_col = celltype_col,
      output_file = paste0(prefix, "_by_sample_", tag, ".pdf"),
      output_dir = figure_dir
    )
    plot_celltype_distribution(
      obj_sub,
      sample_col = group_col,
      celltype_col = celltype_col,
      output_file = paste0(prefix, "_by_group_", tag, ".pdf"),
      output_dir = figure_dir
    )
    save_proportion_tables(
      obj_sub,
      sample_col,
      celltype_col,
      table_dir,
      paste0(prefix, "_by_sample_", tag)
    )

    if (make_sankey) {
      plot_celltype_sankey(
        obj_sub,
        source_col = group_col,
        target_col = celltype_col,
        path_no_ext = file.path(figure_dir, paste0(prefix, "_sankey_group_to_celltype_", tag)),
        min_count = sankey_min_count,
        width = 12,
        height = 8
      )
    }
  }
}
