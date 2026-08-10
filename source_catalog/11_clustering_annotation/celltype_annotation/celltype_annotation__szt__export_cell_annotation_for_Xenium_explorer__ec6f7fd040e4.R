library(tidyverse)
library(Seurat)

all_data <- readRDS('output/data/adata_filtered_add_celltype.rds')

out_dir_export <- "export annotation"
dir.create(out_dir_export)


all_data[[]] %>% 
  dplyr::select(
    fov, 
    cell_id = cell_id2, 
    group = annotation
  ) %>%
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      file = file.path(out_dir_export, str_glue("{name}.csv")),
      row.names = FALSE
    )
  })
