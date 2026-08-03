library(Seurat)
library(SeuratDisk)
library(tidyverse)

# 读取 Seurat 对象
sce <- readRDS("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/T_cells_R/CD8_T/sub_CD8/CD8T_sub_harmony_umap_tsne_cell_type.rds")

# 如果是 Seurat v5，建议确认默认 assay
DefaultAssay(sce) <- "RNA"

# 保存为 h5Seurat
SaveH5Seurat(
  sce,
  filename = "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/CD8T.h5Seurat",
  overwrite = TRUE
)

# 转换为 h5ad
Convert(
  "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/CD8T.h5Seurat",
  dest = "h5ad",
  overwrite = TRUE
)

sce@meta.data %>% as.data.frame() %>% rownames_to_column('cell_barcode') %>% write.csv("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/CD8T_meta.csv")
