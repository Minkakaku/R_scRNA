# 加载单细胞复用函数；先通过 project_config.R 设置 WORKSPACE_DIR 和项目参数。
suppressPackageStartupMessages(invisible(lapply(PACKAGES, library, character.only = TRUE)))
set.seed(SEED)

src_dir <- file.path(WORKSPACE_DIR, "src", "R")
core_files <- sort(list.files(src_dir, pattern = "^[0-9]{2}_.+[.]R$", full.names = TRUE))
step_files <- sort(list.files(file.path(src_dir, "steps"), pattern = "^[0-9]{2}_.+[.]R$", full.names = TRUE))
invisible(lapply(c(core_files, step_files), sys.source, envir = .GlobalEnv))
