# Run from repository root: Rscript tests/test_project.R
root <- normalizePath(".")
tmp <- tempfile("project-layout-")
dir.create(tmp)
workspace <- file.path(tmp, "workspace")
dir.create(workspace)
file.copy(file.path(root, "config"), workspace, recursive = TRUE)
for (name in c("a", "b")) {
  project <- file.path(workspace, "projects", name)
  dir.create(project, recursive = TRUE)
  file.copy(file.path(root, "projects/project1/config"), project, recursive = TRUE)
}
source(file.path(root, "src/R/project.R"))
project <- file.path(workspace, "projects/a")
load_project(project, load_methods = FALSE)
a <- RESULTS_DIR
stopifnot(PROJECT_DIR == normalizePath(project), startsWith(a, PROJECT_DIR))
project <- file.path(workspace, "projects/b")
load_project(project, load_methods = FALSE)
stopifnot(RESULTS_DIR != a, startsWith(RESULTS_DIR, PROJECT_DIR))
stopifnot(identical(utils::modifyList(list(a = list(x = 1, y = 2)), list(a = list(y = 3))), list(a = list(x = 1, y = 3))))
unlink(tmp, recursive = TRUE)
cat("Project isolation checks passed\n")
