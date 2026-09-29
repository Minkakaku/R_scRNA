step_subcluster_proportion <- function() {
  log_msg("Step 10 subcluster proportion started: ", SUBCLUSTER_NAME)
  ensure_subcluster_output_dirs()

  obj <- read_qs2(QS2_SUBCLUSTER_ANNOTATED)
  obj <- prepare_subcluster_input(obj)

  if (!"sub_cell_type" %in% colnames(obj@meta.data)) {
    stop("Missing sub_cell_type. Run the subcluster annotation step first.", call. = FALSE)
  }

  save_proportion_outputs(
    obj,
    comparisons = SUBCLUSTER$comparisons,
    sample_col = "samples",
    group_col = "group",
    celltype_col = "sub_cell_type",
    output_dir = SUBCLUSTER_DIR,
    table_dir = SUBCLUSTER_TABLE_DIR,
    figure_dir = SUBCLUSTER_DIR,
    prefix = paste0(SUBCLUSTER_NAME, "_subtype"),
    make_sankey = SUBCLUSTER$make_sankey,
    sankey_min_count = PROPORTION$sankey_min_count
  )

  log_msg("Step 10 done.")
}
