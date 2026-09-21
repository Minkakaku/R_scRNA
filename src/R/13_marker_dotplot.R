# 参考新文章思路 step79 的点图+细胞数侧栏；表达统计仍由 Seurat 计算。
# 不修改对象、assay、marker 顺序或分析结果；重复 marker 仅绘制一次。
hongf_marker_dotplot <- function(object, features, group.by = NULL,
                                assay = Seurat::DefaultAssay(object), dot.range = c(1, 7.2)) {
  genes <- unique(unlist(features, use.names = FALSE))
  genes <- genes[genes %in% rownames(object[[assay]])]
  if (!length(genes)) stop("No marker genes present in assay: ", assay)
  ids <- if (is.null(group.by)) Seurat::Idents(object) else object[[group.by, drop = TRUE]]
  if (anyNA(ids)) stop("Missing marker group labels")
  labels <- if (is.factor(ids)) levels(droplevels(ids)) else unique(as.character(ids))
  d <- Seurat::DotPlot(object, features = genes, group.by = group.by, assay = assay)$data
  # 单一组或常量基因的组间 Z 分数未定义，显示为中性 0，不更改原始表达值。
  d$avg.exp.scaled[!is.finite(d$avg.exp.scaled)] <- 0
  p_dot <- hongf_dotplot_style(ggplot2::ggplot(d,
    ggplot2::aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::geom_point()) +
    ggplot2::scale_size_continuous(name = "% expressed", range = dot.range,
      limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::scale_x_discrete(limits = genes, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = rev(labels), drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, color = "Average expression")
  counts <- data.frame(id = labels, n_cells = as.integer(table(factor(ids, levels = labels))))
  stopifnot(sum(counts$n_cells) == ncol(object))
  counts$label <- format(counts$n_cells, big.mark = ",", scientific = FALSE, trim = TRUE)
  p_count <- ggplot2::ggplot(counts, ggplot2::aes(x = n_cells, y = id)) +
    ggplot2::geom_col(width = 0.62, fill = "#F28E75", color = "#A13B2F", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.08, size = 3, color = "#333333") +
    ggplot2::scale_y_discrete(limits = rev(labels), drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.5))) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 12, 5.5, 2))
  # 收集图例至整图底部，侧栏与矩阵采用完全相同的 y 轴范围。
  (p_dot + p_count + patchwork::plot_layout(widths = c(8.7, 1.3), guides = "collect")) &
    ggplot2::theme(legend.position = "bottom")
}
