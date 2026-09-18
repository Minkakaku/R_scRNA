# 读取已注释对象，输出样本及组别组成结果并保存最终对象。
step_cell_proportion <- function() {
  log_msg("Step 06 cell proportion started")

  sc <- read_qs2(QS2_ANNOTATED)
  sc <- assign_group_metadata(sc, required_sheet = TRUE)

  if (!"cell_type" %in% colnames(sc@meta.data)) {
    stop("Missing cell_type. Run the cell annotation step first.", call. = FALSE)
  }

  save_proportion_outputs(
    sc,
    comparisons = PROPORTION$comparisons,
    sample_col = "samples",
    group_col = "group",
    celltype_col = "cell_type",
    output_dir = PROPORTION_DIR,
    table_dir = PROPORTION_TABLE_DIR,
    figure_dir = PROPORTION_DIR,
    prefix = "sc_proportion",
    make_sankey = PROPORTION$make_sankey,
    sankey_min_count = PROPORTION$sankey_min_count
  )

  save_qs2(sc, QS2_FINAL)
  log_msg("Step 06 done.")
}
