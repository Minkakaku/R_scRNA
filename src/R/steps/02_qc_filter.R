# 使用逐样本阈值过滤细胞，保存保留率、过滤后 QC 图和中间对象。
step_qc_filter <- function() {
  log_msg("Step 02 QC filtering started")

  sc <- read_qs2(QS2_PRE_QC)
  cells_before <- table(sc$orig.ident)
  criteria <- qc_filter_pass_table(sc)
  sc <- run_step("Apply sample QC thresholds", apply_qc_filter(sc))
  sc <- assign_group_metadata(sc)

  if (ncol(sc) == 0) {
    stop("No cells left after sample QC thresholds. Relax parameters in ", QC_FILTER_SHEET_PATH, ".", call. = FALSE)
  }

  retention <- qc_retention_table(cells_before, table(sc$orig.ident))
  write.csv(criteria, file.path(QC_DIR, "QC_filter_criteria_pass_rate.csv"), row.names = FALSE, quote = FALSE)
  write.csv(retention, file.path(QC_DIR, "QC_filter_retention_by_sample.csv"), row.names = FALSE, quote = FALSE)
  write.csv(retention, file.path(QC_DIR, "QC_predicted_retention_current_thresholds.csv"), row.names = FALSE, quote = FALSE)
  write_sample_group_sheet(sc)
  save_qc_plots(sc, "after")
  save_qs2(sc, QS2_QC)

  log_msg("Step 02 QC filtering done. Kept cells: ", ncol(sc), ". Sample groups: ", SAMPLE_SHEET_PATH)
}
