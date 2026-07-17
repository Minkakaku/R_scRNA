`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

log_msg <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

run_step <- function(label, expr) {
  start <- Sys.time()
  log_msg("START ", label)
  out <- force(expr)
  log_msg("DONE ", label, " (", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), " min)")
  out
}

format_list_value <- function(x) {
  if (is.null(x) || length(x) == 0) return("<none>")
  paste(as.character(x), collapse = ",")
}

save_qs2 <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  run_step(paste("Save", path), qs_save(obj, path))
  cat("结果文件已输出: ", path, "\n", sep = "")
  invisible(path)
}

read_qs2 <- function(path) {
  if (!file.exists(path)) stop("Missing QS2: ", path, call. = FALSE)
  run_step(paste("Read", path), qs_read(path))
}

save_plot <- function(p, path_no_ext, width = 8, height = 6) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, limitsize = FALSE, bg = "white")
}

standardize_umap_plot <- function(p, xlim = NULL, ylim = NULL) {
  fixed_coords <- coord_fixed(xlim = xlim, ylim = ylim)
  if (inherits(p, "patchwork")) {
    return(p & fixed_coords)
  }
  p + fixed_coords
}

save_umap_plot <- function(p, path_no_ext, width = UMAP_PLOT$width, height = UMAP_PLOT$height,
  xlim = NULL, ylim = NULL) {
  save_plot(standardize_umap_plot(p, xlim, ylim), path_no_ext, width = width, height = height)
}
