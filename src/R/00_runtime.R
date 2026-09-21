# 当左侧为 NULL 或空向量时使用默认值；其余值原样返回。
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

HONGF_NPG <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
HONGF_NPG_SEQUENTIAL <- c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019")
HONGF_NPG_DIVERGING <- c("#3C5488", "#F7F7F7", "#E64B35")

# 按固定 NPG 顺序生成分类配色；超过 10 类时从 NPG 锚点扩展。
npg_palette <- function(labels) {
  labels <- as.character(labels)
  colors <- HONGF_NPG
  if (length(labels) > length(colors)) {
    extra_n <- length(labels) - length(colors)
    extra <- grDevices::colorRampPalette(HONGF_NPG)(2L * extra_n + 1L)[seq.int(2L, 2L * extra_n, by = 2L)]
    colors <- c(colors, extra)
  }
  stats::setNames(colors[seq_along(labels)], labels)
}

options(
  ggplot2.discrete.colour = function(...) ggplot2::discrete_scale("colour",
    palette = function(n) unname(npg_palette(seq_len(n))), ...),
  ggplot2.discrete.fill = function(...) ggplot2::discrete_scale("fill",
    palette = function(n) unname(npg_palette(seq_len(n))), ...)
)

# 输出带时间戳的进度信息，便于定位各分析步骤。
log_msg <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

# 执行传入表达式并记录起止时间，返回该表达式的计算结果。
run_step <- function(label, expr) {
  start <- Sys.time()
  log_msg("START ", label)
  out <- force(expr)
  log_msg("DONE ", label, " (", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), " min)")
  out
}

# 将参数向量拼接为日志文本；缺失参数显示为 <none>。
format_list_value <- function(x) {
  if (is.null(x) || length(x) == 0) return("<none>")
  paste(as.character(x), collapse = ",")
}

# 创建父目录并保存 QS2 中间对象，同时记录输出路径和耗时。
save_qs2 <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  run_step(paste("Save", path), qs_save(obj, path))
  cat("结果文件已输出: ", path, "\n", sep = "")
  invisible(path)
}

# 确认文件存在后读取 QS2 中间对象，并记录读取耗时。
read_qs2 <- function(path) {
  if (!file.exists(path)) stop("Missing QS2: ", path, call. = FALSE)
  run_step(paste("Read", path), qs_read(path))
}

# 按同一文件名前缀保存白底 PNG 和 PDF，尺寸单位为英寸。
save_plot <- function(p, path_no_ext, width = 8, height = 6) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  if (inherits(p, "patchwork")) {
    p <- p & ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  } else if (inherits(p, "ggplot")) {
    p <- p + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  }
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 600, limitsize = FALSE, bg = "white")
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, limitsize = FALSE, bg = "white")
}

# 默认输出非栅格 UMAP，并以标签的稳定顺序应用 NPG 分类色。
npg_dimplot <- function(object, group.by = NULL, ..., raster = FALSE) {
  values <- if (is.null(group.by)) Seurat::Idents(object) else object[[group.by, drop = TRUE]]
  labels <- if (is.factor(values)) levels(droplevels(values)) else sort(unique(as.character(values)))
  Seurat::DimPlot(object, group.by = group.by, ..., raster = raster, cols = npg_palette(labels))
}

# 参考点图：白底细边框、无内部网格、纵向标签、底部双图例。
HONGF_DOTPLOT_COLORS <- c("#FDEBE5", "#F6B29D", "#E85C47", "#A50F2D", "#67000D")
hongf_dotplot_style <- function(p) {
  p +
    ggplot2::scale_color_gradientn(colours = HONGF_DOTPLOT_COLORS,
      guide = ggplot2::guide_colorbar(display = "rectangles")) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, title.position = "top", nrow = 1,
        override.aes = list(colour = "black")),
      colour = ggplot2::guide_colorbar(order = 2, title.position = "top", display = "rectangles",
        barwidth = grid::unit(30, "mm"), barheight = grid::unit(3, "mm"))) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.55),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5)
    )
}

# 为普通图或 patchwork 各面板设置等比例坐标及可选统一范围。
standardize_umap_plot <- function(p, xlim = NULL, ylim = NULL) {
  fixed_coords <- coord_fixed(xlim = xlim, ylim = ylim)
  if (inherits(p, "patchwork")) {
    return(p & fixed_coords)
  }
  p + fixed_coords
}

# 统一 UMAP 坐标比例后，使用项目配置的尺寸输出 PNG/PDF。
save_umap_plot <- function(p, path_no_ext, width = UMAP_PLOT$width, height = UMAP_PLOT$height,
  xlim = NULL, ylim = NULL) {
  save_plot(standardize_umap_plot(p, xlim, ylim), path_no_ext, width = width, height = height)
}
