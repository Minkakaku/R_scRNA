library(tidyverse)
library(Seurat)

load("Filter_Basicinfo.rData")

# only used fibroblast
all_data_rctd <- readRDS("sc_stma.rds") %>%
  .[, .$annotation %in% c(
    'S3 fibroblasts',
    'Myofibroblasts',
    # Pericytes                                  
    # Endothelial cells
    'S1 fibroblasts',
    'Immediate early-response (IER) fibroblasts',
    'Peri-cryptal S2 fibroblasts'
  )]

all_data_sp_rctd <- readRDS("sp_subcluster/sp_fibro.rds")

cells_position_ls <- readRDS('data_for_analysis/cells_position_ls.rds')


# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD' 
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
# all_data_rctd <- readRDS("all_data_sc_dim15_anno.rds")
# all_data_rctd <- all_data


Idents(all_data_rctd) <- "annotation"
counts_sc  <- all_data_rctd %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_rctd$annotation)
nUMI_sc    <- all_data_rctd$nCount_RNA

sc_cell_filter_idx <- ! cluster_sc == "other"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
all_data_sp_rctd_ls <- all_data_sp_rctd %>% JoinLayers() %>% SplitObject('sample')

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
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

#4. run RCTD
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 10, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

# saveRDS(rctd_result_ls, 'RCTD/rctd_result_ls.rds')

#5. addMetaData
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })




#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data_sp_rctd %>% 
  AddMetaData(meta_rctd_result)

saveRDS(all_data_sp_rctd_anno, file.path(rctd_out_dir, "xenium_rctd.rds"))

# spatial plot --------------

all_data_meta <-
  all_data_sp_rctd_anno[[]] %>%
  rownames_to_column("cell_id")
cells_boundary_ls<-readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls<-readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))
image_wh_r_ls<-readRDS(file.path("data_for_analysis", "image_wh_r_ls.rds"))

sp_p_ls1 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      dplyr::filter(!is.na(first_type))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = first_type, group = cell_id))  +
      scale_y_reverse() +
      coord_equal()+
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

# export image
sp_dir <- "spatial plot" # dir.create(sp_dir)

sp_p_ls1 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })



