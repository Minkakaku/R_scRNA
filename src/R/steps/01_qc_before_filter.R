# 读取和合并样本，保存过滤前 QC 图、对象及待人工确认的阈值表。
step_qc_before_filter <- function() {
  log_msg("Step 01 QC before filtering started")

  if (!dir.exists(DATA_DIR)) stop("Data dir not found: ", DATA_DIR, call. = FALSE)
  sample_manifest <- discover_samples(DATA_DIR)
  samples <- sample_manifest$sample
  if (length(samples) == 0) stop("No sample matrices under ", DATA_DIR, call. = FALSE)

  obj_list <- setNames(lapply(samples, read_qc_sample), samples)
  sc <- run_step("Merge samples", merge_samples(obj_list))
  sc <- assign_group_metadata(sc)

  save_qc_plots(sc, "before")
  write_qc_filter_sheet(sc)
  save_qs2(sc, QS2_PRE_QC)

  log_msg("Step 01 QC before filtering done. Edit ", QC_FILTER_SHEET_PATH, ", then run the QC filtering step.")
}
