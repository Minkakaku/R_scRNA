library(tidyverse)
library(Seurat)


# read data -------
# sample_df <- 
#   # read.csv("data/Annoation.csv")
#   data.frame(
#     sample = c('1W2W3W'),
#     group  = c('1W2W3W')
#   )
  
file_list <- 
  list.files(
    "/mnt/NexData_file/outputdata/outputdata/spaceranger_out/FWSP20250659/",
    pattern = paste(sample_df$sample, collapse = "|"),
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]

names(file_list) <- basename(file_list)


visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, image.name = "tissue_lowres_image.png", 
        assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir, 
            "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, 
            "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
            type = "centroids", radius = scale.factors[["spot"]], 
            assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
            molecules = fov@molecules, assay = fov@assay, key = fov@key, 
            image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(file.path(path, "binned_outputs/square_008um/spatial"), assay = 'Spatial.008um')
    )
  })

visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })

# all_data
all_data <- readRDS('/mnt/NexData_file/home/wyt/S203_project/project/FWSP20250659_W123_PBS_mouse_heart/rds/lower_W1234CTR_celltype_rctd_anno.rds')

all_data$fov <- ifelse(all_data$sample %in% c("KD1W", "KD2W", "KD3W"), '1W2W3W', 'Others')
all_data$cell_id_ori <- all_data %>% colnames() %>% str_extract("s.*$")
all_data$cell_id2 <- str_c(all_data$fov, all_data$cell_id_ori, sep = "_")


table( colnames(visim_ls$`1W2W3W`) %in% all_data$cel_id2 )


# combine data
visium <- visim_ls$`1W2W3W`

visium@meta.data <-
  visium@meta.data %>%
  rownames_to_column('cell_id') %>%
  left_join(
    all_data[[]],
    by = c(cell_id = 'cell_id2')
  ) %>%
  column_to_rownames("cell_id")

visium <- visium[, !is.na(visium$fov)]


saveRDS(visium, 'visium.rds')


# visim_ls$`1W2W3W`@images$X1W2W3W.008um@image %>% dim
# visim_ls$`1W2W3W`@images$X1W2W3W.008um %>% GetTissueCoordinates(scale = "lowres") %>% {c(range(.$x), range(.$y))}
# 
# 
# visim_ls$`1W2W3W`@images$X1W2W3W.008um@image %>% grid::grid.raster()




