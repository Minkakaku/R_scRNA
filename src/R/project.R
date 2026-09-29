# 加载共享 R 默认值及 YAML 覆盖；保留原有分析函数接口。
load_project <- function(project_dir, load_methods = TRUE) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  workspace <- normalizePath(file.path(project_dir, "../.."), mustWork = TRUE)
  previous <- setwd(project_dir)
  on.exit(setwd(previous), add = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(workspace, "config/default.R"), envir = env)
  for (path in c(file.path(workspace, "config/default.yaml"),
                 file.path(workspace, "config/reference.yaml"),
                 file.path(project_dir, "config/config.yaml"))) {
    overrides <- yaml::read_yaml(path)
    for (key in names(overrides)) {
      value <- overrides[[key]]
      old <- if (exists(key, env, inherits = FALSE)) env[[key]] else NULL
      env[[key]] <- if (is.list(old) && is.list(value)) utils::modifyList(old, value, keep.null = TRUE) else value
    }
  }
  # 路径统一绝对化，调用者无须 setwd；配置中的相对路径基于项目目录。
  resolve <- function(x) {
    if (!is.character(x)) return(x)
    ifelse(grepl("^(/|[A-Za-z]:)", x), x, file.path(project_dir, x))
  }
  for (key in ls(env)) {
    if (grepl("(_DIR|_PATH)$|^QS2_", key) && is.character(env[[key]])) env[[key]] <- resolve(env[[key]])
  }
  for (key in c("ANNOTATION", "SUBCLUSTER")) {
    for (field in intersect(names(env[[key]]), c("template_path", "mapping_path", "input_qs2"))) {
      env[[key]][[field]] <- resolve(env[[key]][[field]])
    }
  }
  env$PROJECT <- lapply(env$PROJECT, resolve)
  env$PROJECT_DIRS <- resolve(env$PROJECT_DIRS)
  env$WORKSPACE_DIR <- workspace
  env$SAMPLE_SHEET_PATH <- file.path(project_dir, "config/samples.tsv")
  list2env(as.list(env), envir = .GlobalEnv)
  if (load_methods) sys.source(file.path(workspace, "src/R/load_all.R"), envir = .GlobalEnv)
  invisible(PROJECT)
}
