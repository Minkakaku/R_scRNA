# 提取 config.yaml 指定的大类细胞并重新聚类。
source("src/R/project_config.R")
load_project("projects/project1")
step_subcluster_clustering()
