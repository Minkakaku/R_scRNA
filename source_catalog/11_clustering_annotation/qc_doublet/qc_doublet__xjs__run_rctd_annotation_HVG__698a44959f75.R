library(tidyverse)
library(Seurat)


# only used fibroblast
all_data_rctd <- readRDS('RCTD/SC2510275964.rds') #  sc_data # sc
all_data_rctd$annotation <- all_data_rctd$cell_type  
all_data_rctd[['RNA']] <- as(all_data_rctd[['RNA']], "Assay5")


all_data_sp_rctd <- readRDS(file.path("R_OUTPUT", "all_data.rds"))  # all_data   # st
cells_position_ls <-
  all_data_sp_rctd[[]] %>%
  dplyr::select(
    project,
    cell_id = cell_id_sample,
    x_centroid = x,
    y_centroid = y
  ) %>%
  group_by(project) %>%
  nest() %>%
  deframe()


stopifnot(class(all_data_rctd[['RNA']]) == "Assay5")
stopifnot(class(all_data_sp_rctd[['RNA']]) == "Assay5")



# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD_HVG_5K' 
dir.create(rctd_out_dir)

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

#2. single data
message(1)
Idents(all_data_rctd) <- "annotation"
counts_sc  <- all_data_rctd %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_rctd$annotation)
nUMI_sc    <- all_data_rctd$nCount_RNA

sc_cell_filter_idx <- ! cluster_sc == "Unknown"

message(2)
# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
message(3)
all_data_sp_rctd_ls <- all_data_sp_rctd %>% JoinLayers() %>% SplitObject('project')

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    # counts_hd <- xenium %>% LayerData(layer = "counts")
    
    xenium <- NormalizeData(xenium, scale.factor = median(xenium$nCount_RNA))
    xenium <- FindVariableFeatures(xenium, nfeatures = 5000)
    counts_hd <- xenium %>% LayerData(layer = "counts", features = VariableFeatures(xenium))
    
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

message(4)
#4. run RCTD
rctd_result_ls <- 
  rctd_query_ls %>% # .[. %in% c('LFX')] %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 100, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    saveRDS(RCTD, str_glue("{rctd_out_dir}/RCTD_{name}.rds"))
    
    RCTD
  })

saveRDS(rctd_result_ls, str_glue('{rctd_out_dir}/rctd_result_ls.rds'))
        
        
        