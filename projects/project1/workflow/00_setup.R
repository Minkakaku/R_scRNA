# 初始化项目目录并检查 YAML 能否读取；不运行分析。
source("src/R/project_config.R")
load_project_config("projects/project1")
message("Project ready: ", PROJECT_DIR)
