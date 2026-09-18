# 读取各样本矩阵，计算 QC 指标并生成逐样本阈值表。
source("src/R/project_config.R")
load_project("projects/project1")
step_qc_before_filter()
