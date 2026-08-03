library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

adata_celltype_anno <-
  sceasy::convertFormat(
    "merged_single_cell_data_with_group_new.h5ad", 
    from = "anndata", 
    to = "seurat",
    main_layer = "scale.data" # <----- 
  )

saveRDS(adata_celltype_anno, 'adata_celltype_anno.rds')
