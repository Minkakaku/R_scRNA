# 最小配置测试：确认公共配置、项目覆盖和路径路由均可用。
source("src/R/project_config.R")
project <- load_project_config("projects/project1")

stopifnot(
  basename(project$project_dir) == "project1",
  identical(CLUSTER$dims, 1:20),
  identical(GROUP_LEVELS, c("CON", "EXP")),
  identical(SUBCLUSTER$target_celltype, c("Monocyte", "Macrophage")),
  dirname(SAMPLE_SHEET_PATH) == CONFIG_DIR,
  QS2_DIR == DATA_PROCESSED_DIR,
  startsWith(RESULTS_DIR, PROJECT_DIR)
)

cat("project config: OK\n")
