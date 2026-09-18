# 当左侧为 NULL 或空向量时使用默认值；其余值原样返回。
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

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
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, limitsize = FALSE, bg = "white")
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
