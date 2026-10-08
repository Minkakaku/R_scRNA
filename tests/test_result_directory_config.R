# 使用临时项目核对目录初始化开关，避免依赖真实数据或修改现有项目配置。
source("src/R/project_config.R")
test_root <- tempfile("project_config_")
dir.create(file.path(test_root, "config"), recursive = TRUE)
file.copy(c("config/default.yaml", "config/reference.yaml"), file.path(test_root, "config"))
for (initialize in c(FALSE, TRUE)) {
  test_project <- file.path(test_root, "projects", as.character(initialize))
  dir.create(file.path(test_project, "config"), recursive = TRUE)
  yaml::write_yaml(list(initialize_result_subdirs = initialize), file.path(test_project, "config/config.yaml"))
  paths <- load_project_config(test_project)
  stopifnot(dir.exists(paths$results_dir), dir.exists(paths$data_processed_dir),
            all(vapply(paths[c("qc_dir", "global_dir", "marker_dir", "annotation_dir", "proportion_dir", "subcluster_dir")],
                       dir.exists, logical(1)) == initialize))
}
cat("RESULT_DIRECTORY_CONFIG_OK\n")
