# Discover Cell Ranger matrices under data/.
discover_samples <- function(data_dir = DATA_DIR) {
  matrix_dirs <- list.dirs(data_dir, recursive = TRUE, full.names = TRUE)
  matrix_dirs <- matrix_dirs[
    basename(matrix_dirs) == "filtered_feature_bc_matrix" &
      basename(dirname(matrix_dirs)) == "outs"
  ]

  matrix_h5 <- list.files(
    data_dir,
    pattern = "^filtered_feature_bc_matrix[.]h5$",
    recursive = TRUE,
    full.names = TRUE
  )
  matrix_h5 <- matrix_h5[basename(dirname(matrix_h5)) == "outs"]

  matrix_path <- c(matrix_dirs, matrix_h5)
  sample_dir <- dirname(dirname(matrix_path))
  sample <- basename(sample_dir)

  manifest <- data.frame(
    sample = sample,
    group = sample,
    sample_dir = sample_dir,
    outs_dir = dirname(matrix_path),
    matrix_dir = matrix_path,
    matrix_type = ifelse(grepl("[.]h5$", matrix_path), "h5", "mtx_dir"),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$sample, match(manifest$matrix_type, c("h5", "mtx_dir"))), ]
  manifest[!duplicated(manifest$sample), , drop = FALSE]
}
