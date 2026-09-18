# 使用人工确认后的 qc_thresholds.tsv 过滤细胞。
source("src/R/project_config.R")
load_project("projects/project1")
step_qc_filter()
