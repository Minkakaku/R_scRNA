cd /mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium_tnbc/sub/celltype_subsets/combine
conda activate tools

R

library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

adata_celltype_anno <-
  sceasy::convertFormat(
    "merged_single_cell_data_modified.h5ad", 
    from = "anndata", 
    to = "seurat",
    main_layer = "scale.data" # <----- 
  )

saveRDS(adata_celltype_anno, 'adata_celltype_anno.rds')
