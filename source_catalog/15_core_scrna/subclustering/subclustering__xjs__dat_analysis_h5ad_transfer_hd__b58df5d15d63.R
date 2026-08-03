library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
data_dir  <- "/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20251583-HKPG-HG/"
list.files(data_dir)
anno_file <- "data/Annotation.csv"

dir.create("R_OUTPUT")
dir.create("R_OUTPUT/data_for_analysis")


# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
#library(RBioFormats)
# source("xenium.R")

# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# 
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()


# read data -------
sample_df <- 
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('FWSP20251583-HKPG-HG'),
    group  = c('FWSP20251583-HKPG-HG')
  )

file_list <- 
  list.files(
    data_dir,
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]




visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, image.name = "tissue_lowres_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir, "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(
          coordinates[, c("imagerow", "imagecol")], 
          type = "centroids", 
          radius = scale.factors[["spot"]], 
          assay = assay, 
          key = Key(slice, quiet = TRUE)
        )
        visium.fov <- new(
          Class = "VisiumV2", 
          boundaries = fov@boundaries, 
          molecules = fov@molecules, 
          assay = fov@assay,
          key = fov@key, 
          image = image, 
          scale.factors = scale.factors
        )
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(file.path(path, "binned_outputs/square_008um/spatial"), assay = 'Spatial.008um')
    )
  })

# visim_ls <-
#   file_list %>% 
#   imap(function(path, name){
#     path <- file.path(path, "outs")

#     Read10X_Image <-
#       function (image.dir, spatial.dir = image.dir, image.name = "tissue_lowres_image.png", 
#                 assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {

#         image <- png::readPNG(source = file.path(image.dir, image.name))
#         scale.factors <- Read10X_ScaleFactors(filename = file.path(spatial.dir, 
#                                                                    "scalefactors_json.json"))
#         coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(spatial.dir, 
#                                                                          "*tissue_positions*parquet")), filter.matrix)
#         fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
#                          type = "centroids", radius = scale.factors[["spot"]], 
#                          assay = assay, key = Key(slice, quiet = TRUE))
#         visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
#                           molecules = fov@molecules, assay = fov@assay, key = fov@key, 
#                           image = image, scale.factors = scale.factors)
#         return(visium.fov)
#       }

#     Load10X_Spatial(
#       path, 
#       bin.size = 8, 
#       slice = name, 
#       image = Read10X_Image(
#         file.path(path, "spatial"),
#         spatial.dir = file.path(path, "binned_outputs/square_008um/spatial"),
#         assay = 'Spatial.008um'
#       )
#     )
#   })

# add sample_id or sample_name
visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 

all_data <- merge(visim_ls[[1]], visim_ls[-1])


## read coords -------------------
# cell coordinate and cell boundary
cells_position_ls <- 
  all_data@images %>%
  map(GetTissueCoordinates) 

# export data
saveRDS(cells_position_ls, file.path("R_OUTPUT", "data_for_analysis", "cells_position_ls.rds"))

# FORMAT CONVERSION ---------
# config python
library(reticulate)
use_condaenv(condaenv = "/home/szt/miniforge3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


slide_name <- 'FWSP20251583-HKPG-HG'
stopifnot(length(slide_name) == 1)

sc <- import("scanpy")
adata <- sc$read_h5ad("/home/xjs/case/FWSP20251583-02/output/data/adata_celltype.h5ad")

adata_meta <- adata$obs
adata_umap <- adata$obsm['X_umap']

# cell_id match
adata_meta$cell_id_raw <- adata_meta %>%  rownames() %>% str_sub(end = 21)
adata_meta$cell_id <- str_c(slide_name, adata_meta$cell_id_raw, sep = "_")

stopifnot(all(adata_meta$cell_id %in% colnames(all_data)))


# add METADATA
all_data2 <- all_data %>% subset(cells = adata_meta$cell_id)
all_data2@meta.data <-
  all_data2@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(adata_meta, by = "cell_id") %>%
  column_to_rownames("cell_id")


# add UMAP
adata_umap2 <- adata_umap
colnames(adata_umap2) <- str_c("UMAP", 1:ncol(adata_umap2), sep = "_")
rownames(adata_umap2) <- adata_meta$cell_id
adata_umap2 <- adata_umap2[colnames(all_data2), ]

all_data2[['umap']] <- 
  CreateDimReducObject(
    adata_umap2,
    key    = "UMAP_",
    assay  = DefaultAssay(all_data2),
    global = TRUE
  )

# add useful columns
cluster_res_cols <- all_data2[[]] %>% colnames() %>% str_subset("^res_")
for(col_name in cluster_res_cols){
  col_name_new <- col_name %>% str_replace("res_", str_glue("{DefaultAssay(all_data2)}_snn_res.")) 
  all_data2[[col_name_new]] <- all_data2[[col_name]][[1]]
}
all_data2$annotation <- all_data2$celltype


# plot
all_data2 %>% DimPlot(reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "annotation", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_1", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "Spatial.008um_snn_res.0.8", label = TRUE) + NoLegend()


# export 
saveRDS(all_data2, file = file.path("R_OUTPUT", "all_data.rds"))

