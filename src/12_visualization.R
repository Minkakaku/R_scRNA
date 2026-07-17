make_tag <- function(x) {
  tag <- gsub("[^A-Za-z0-9_]+", "_", paste(x, collapse = "_vs_"))
  gsub("^_+|_+$", "", tag)
}

plot_celltype_distribution <- function(
    sc,
    sample_col,
    celltype_col,
    output_file,
    output_dir = PROPORTION_DIR,
    width = 12,
    height = 8) {
  meta <- sc@meta.data[, c(sample_col, celltype_col), drop = FALSE]
  colnames(meta) <- c("Sample", "CellType")
  sample_levels <- if (is.factor(meta$Sample)) levels(droplevels(meta$Sample)) else unique(meta$Sample)
  meta$Sample <- factor(meta$Sample, levels = sample_levels)
  meta$CellType <- factor(meta$CellType)

  palette <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(meta$CellType))
  names(palette) <- levels(meta$CellType)
  df <- as.data.frame(table(meta$Sample, meta$CellType))
  colnames(df) <- c("Sample", "CellType", "Number")

  p <- ggplot(df, aes(Sample, Number, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, output_file), p, width = width, height = height)
}

plot_celltype_sankey <- function(
    sc,
    source_col,
    target_col,
    path_no_ext,
    min_count = 1L,
    width = 10,
    height = 7) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    warning("Skip Sankey plot because the optional ggalluvial package is not installed.")
    return(invisible(NULL))
  }

  if (!all(c(source_col, target_col) %in% colnames(sc@meta.data))) {
    warning("Skip Sankey plot because columns are missing: ", source_col, ", ", target_col)
    return(invisible(NULL))
  }

  meta <- sc@meta.data[, c(source_col, target_col), drop = FALSE]
  colnames(meta) <- c("Source", "Target")
  meta <- meta[!is.na(meta$Source) & !is.na(meta$Target), , drop = FALSE]
  df <- as.data.frame(table(meta$Source, meta$Target), stringsAsFactors = FALSE)
  colnames(df) <- c("Source", "Target", "Number")
  df <- df[df$Number >= min_count, , drop = FALSE]
  if (nrow(df) == 0) {
    warning("Skip Sankey plot because no flow reaches min_count=", min_count)
    return(invisible(NULL))
  }

  target_levels <- sort(unique(as.character(df$Target)))
  palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(target_levels))
  names(palette) <- target_levels

  p <- ggplot(
    df,
    aes(axis1 = Source, axis2 = Target, y = Number)
  ) +
    ggalluvial::geom_alluvium(aes(fill = Target), width = 1 / 12, alpha = 0.85) +
    ggalluvial::geom_stratum(width = 1 / 12, fill = "grey95", color = "grey45") +
    ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c(source_col, target_col), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  save_plot(p, path_no_ext, width = width, height = height)
  invisible(path_no_ext)
}

save_pheatmap_plot <- function(p, path_no_ext, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)

  png(paste0(path_no_ext, ".png"), width = width * dpi, height = height * dpi, res = dpi)
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = dev.off())

  pdf(paste0(path_no_ext, ".pdf"), width = width, height = height)
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = dev.off())
}
