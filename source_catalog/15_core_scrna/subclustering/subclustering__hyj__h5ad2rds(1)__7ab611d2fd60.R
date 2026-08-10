Read10X_Image_spatial <- function(
    object,
    image.dir,
    label,
    cov = NULL,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice,
    filter.matrix = TRUE,
    image.type = "VisiumV2"
) {
  # Validate the `image.type` parameter.
  
  image.type <- match.arg(image.type, choices = c("VisiumV1", "VisiumV2"))
  
  # Read in the H&E stain image.
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  # Read in the scale factors.
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )
  
  # Read in the tissue coordinates as a data.frame.
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )
  object$sample = as.character(object$sample)
  object = subset(object,sample %in% label)
  rownames(coordinates) <- paste0(label,"_",rownames(coordinates))
  row_indices <- match(colnames(object), rownames(coordinates))
  indices <- coordinates[row_indices, ]
  key <- Key(slice, quiet = TRUE)
  
  image_X_min <- ceiling(min(indices$imagerow)*scale.factors$lowres)
  image_Y_min <- ceiling(min(indices$imagecol)*scale.factors$lowres) 
  image_X_max <- ceiling(max(indices$imagerow)*scale.factors$lowres) 
  image_Y_max <- ceiling(max(indices$imagecol)*scale.factors$lowres)  
  
  sub_image <- image[image_X_min:image_X_max, image_Y_min:image_Y_max,]
  
  # Use the `slice` value to populate a Seurat-style identifier for the image.
  key <- Key(slice, quiet = TRUE)
  
  # Return the specified `image.type`.
  if (image.type == "VisiumV1") {
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = indices[colnames(spatial),],
      scale.factors = scale.factors,
      image = sub_image
    )
    
    # As of v5.1.0 `Radius.VisiumV1` no longer returns the value of the 
    # `spot.radius` slot and instead calculates the value on the fly, but we 
    # can populate the static slot in case it's depended on.
    visium.v1@spot.radius <- Radius(visium.v1)
    
    return(visium.v1)
  }
  
  # If `image.type` is not "VisiumV1" then it must be "VisiumV2".
  stopifnot(image.type == "VisiumV2")
  # Create an `sp` compatible `FOV` instanc
  fov <- CreateFOV(
    indices[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = key
  )
  
  new_coord <- fov@boundaries$centroids@coords
  new_coord[,1] <- new_coord[,1] - image_X_min/scale.factors$lowres
  new_coord[,2] <- new_coord[,2] - image_Y_min/scale.factors$lowres
  fov@boundaries$centroids@coords <- new_coord
  
  # Build the final `VisiumV2` instance, essentially just adding `image` and
  # `scale.factors` to the `fov`.
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = sub_image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}



Read10X_Image_spatial_hires <- function(
    object,
    image.dir,
    label,
    cov = NULL,
    image.name = "tissue_hires_image.png",
    assay = "Spatial",
    slice,
    filter.matrix = TRUE,
    image.type = "VisiumV2"
) {
  # Validate the `image.type` parameter.
  
  image.type <- match.arg(image.type, choices = c("VisiumV1", "VisiumV2"))
  
  # Read in the H&E stain image.
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  # Read in the scale factors.
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )
  
  # Read in the tissue coordinates as a data.frame.
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )
  object$sample = as.character(object$sample)
  object = subset(object,sample %in% label)
  rownames(coordinates) <- paste0(label,"_",rownames(coordinates))
  row_indices <- match(colnames(object), rownames(coordinates))
  indices <- coordinates[row_indices, ]
  key <- Key(slice, quiet = TRUE)
  
  image_X_min <- ceiling(min(indices$imagerow)*scale.factors$hires)
  image_Y_min <- ceiling(min(indices$imagecol)*scale.factors$hires) 
  image_X_max <- ceiling(max(indices$imagerow)*scale.factors$hires) 
  image_Y_max <- ceiling(max(indices$imagecol)*scale.factors$hires)  
  
  sub_image <- image[image_X_min:image_X_max, image_Y_min:image_Y_max,]
  
  # Use the `slice` value to populate a Seurat-style identifier for the image.
  key <- Key(slice, quiet = TRUE)
  
  # Return the specified `image.type`.
  if (image.type == "VisiumV1") {
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = indices[colnames(spatial),],
      scale.factors = scale.factors,
      image = sub_image
    )
    
    # As of v5.1.0 `Radius.VisiumV1` no longer returns the value of the 
    # `spot.radius` slot and instead calculates the value on the fly, but we 
    # can populate the static slot in case it's depended on.
    visium.v1@spot.radius <- Radius(visium.v1)
    
    return(visium.v1)
  }
  
  # If `image.type` is not "VisiumV1" then it must be "VisiumV2".
  stopifnot(image.type == "VisiumV2")
  # Create an `sp` compatible `FOV` instanc
  fov <- CreateFOV(
    indices[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = key
  )
  
  new_coord <- fov@boundaries$centroids@coords
  new_coord[,1] <- (new_coord[,1] - image_X_min/scale.factors$hires) * (scale.factors$hires/scale.factors$lowres)
  new_coord[,2] <- (new_coord[,2] - image_Y_min/scale.factors$hires) * (scale.factors$hires/scale.factors$lowres)
  fov@boundaries$centroids@coords <- new_coord
  
  # Build the final `VisiumV2` instance, essentially just adding `image` and
  # `scale.factors` to the `fov`.
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = sub_image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}

Add_meta_umapData_visiumHD <- function(
    object,
    umapfile,
    metedata
){
  ###############Add_metadata##################
  row_cell <- match(rownames(metadata), rownames(object@meta.data))
  meta <- metadata[row_cell, ]
  object@meta.data <- meta
  ################Add_umap#####################
  new_umap <- umap[row_cell, ] %>% as.matrix()
  colnames(new_umap) <- c("UMAP_1", "UMAP_2")
  
  new_umap_reduc <- CreateDimReducObject(
    embeddings = new_umap, 
    key = "UMAP_",
    assay = "Spatial"
  )
  object@reductions$umap <- new_umap_reduc
  colnames(object@assays$Spatial@counts) = colnames(object)
  colnames(object@assays$Spatial) = colnames(object)
  return(object)
}

Read_fov_spatial <- function(data.dir,assay,mols.qv.threshold = 20,label,object,cov = NULL) {
  if(is.null(cov)){
    cells_data <- data.frame(cell_id = paste0(paste(label,collapse = ":"),"_",object@meta.data$cell_id),sample = object@meta.data$sample)
    cells_data_sub <- subset(cells_data,sample %in% label)
    cells <- as.character(cells_data_sub$cell_id)
  } else {
    cells_data <- data.frame(cell_id = paste0(paste(label,collapse = ":"),"_",object@meta.data$cell_id),sample = object@meta.data$sample)
    cells_total <- paste0(paste(label,collapse = ":"),"_",cov$`Cell ID`)
    cells_data_sub <- subset(cells_data,sample %in% label & cell_id %in% cells_total)
    cells <- as.character(cells_data_sub$cell_id)
  }
  features <- rownames(object)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  
  centroids = {
    if (has_dt) {
      centroids <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                             "cells.csv.gz")))
    } else {
      centroids <- read.csv(file.path(data.dir, "cells.csv.gz"))
    }
    centroids <- data.frame(x = centroids$x_centroid, 
                            y = centroids$y_centroid, cell = paste0(paste(label,collapse = ":"),"_",centroids$cell_id), 
                            stringsAsFactors = FALSE)}
  
  segmentations = {
    if (has_dt) {
      segmentations <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                 "cell_boundaries.csv.gz")))
    } else {
      segmentations <- read.csv(file.path(data.dir, 
                                          "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
    }
    segmentations <- data.frame(x = segmentations$vertex_x,y = segmentations$vertex_y,
                                cell = paste0(paste(label,collapse = ":"),"_",segmentations$cell_id))}
  
  # microns = {
  #   if (has_dt) {
  #     microns <- data.table::fread(file.path(data.dir, "transcripts.csv"))
  #     microns <- subset(microns, qv >= mols.qv.threshold)
  #   } else {
  #     microns <- read.csv(file.path(data.dir, "transcripts.csv"))
  #     microns <- subset(microns, qv >= mols.qv.threshold)
  #   }
  #   microns <- data.frame(x = microns$x_location, y = microns$y_location,
  #                         cell = paste0(paste(label,collapse = ":"),"_",microns$cell_id),
  #                         gene = microns$feature_name, stringsAsFactors = FALSE)
  # }
  # 
  microns = {
  microns <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  microns <- subset(microns, qv >= mols.qv.threshold)
  microns <- data.frame(x = microns$x_location, y = microns$y_location,
  cell = paste0(paste(label,collapse = ":"),"_",microns$cell_id),
  gene = microns$feature_name, stringsAsFactors = FALSE)
  }

  segmentations <- subset(segmentations,cell %in% cells)
  centroids <- subset(centroids,cell %in% cells)
  microns <- subset(microns,gene %in% features & cell %in% cells)
  
  segmentations.data <- list(centroids = CreateCentroids(centroids), 
                             segmentation = CreateSegmentation(segmentations))
  
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), molecules = microns, assay = assay)
  return(coords)
}

Add_meta_umapData_Xenium <- function(
    object,
    umapfile,
    metedata
){
  ###############Add_metadata##################
  row_cell <- match(rownames(metadata), rownames(object@meta.data))
  meta <- metadata[row_cell, ]
  object@meta.data <- meta
  ################Add_umap#####################
  new_umap <- umap[row_cell, ] %>% as.matrix()
  colnames(new_umap) <- c("UMAP_1", "UMAP_2")
  
  new_umap_reduc <- CreateDimReducObject(
    embeddings = new_umap, 
    key = "UMAP_",
    assay = "Xenium"
  )
  object@reductions$umap <- new_umap_reduc
  colnames(object@assays$Xenium) = colnames(object)
  return(object)
}





read_visium_cellsegemnt <- function(object,gdf,label,assay = "Xenium",microns = NULL){
  boundries <- gdf %>%  rename(cell_id = id) %>% 
    select(cell_id, geometry) %>% 
    mutate(
      geometry = geometry %>% 
        str_replace('.*\\((.*)\\){2}', '\\1') 
    ) %>% 
    separate_rows(geometry, sep = ',') %>%
    map(str_trim) %>%
    as.data.frame() %>% 
    separate(geometry, c('vertex_x', 'vertex_y'), sep = ' ') %>%
    mutate(
      vertex_x = vertex_x %>% as.numeric(), 
      vertex_y = vertex_y %>% as.numeric(),
      cell_id = paste0(label, "_", cell_id)
    )
  
  cell_csv <- gdf %>% select(-geometry,-X) %>% rename(cell_id = id) %>%
    mutate(centroid = sub("POINT \\((.*)\\)", "\\1", centroid)) %>%
    separate(centroid, c('centroid_x', 'centroid_y'), sep = ' ') %>%
    mutate(
      centroid_x = as.numeric(centroid_x), 
      centroid_y = as.numeric(centroid_y),
      cell_id = paste0(label, "_", cell_id)
    )
  
  centroids <- data.frame(x = cell_csv$centroid_x, 
                          y = cell_csv$centroid_y, cell = cell_csv$cell_id, 
                          stringsAsFactors = FALSE)
  object$sample = as.character(object$sample)
  object = subset(object,sample %in% label)
  row_indices <- match(colnames(object),centroids$cell)
  centroids <- centroids[row_indices, ]
  
  segmentations <- data.frame(x = boundries$vertex_x,y = boundries$vertex_y,
                              cell = boundries$cell_id)
  row_indices_seg <- match(colnames(object),segmentations$cell)
  segmentations <- segmentations[row_indices_seg, ]
  
  segmentations.data <- list(centroids = CreateCentroids(centroids), 
                             segmentation = CreateSegmentation(segmentations))
  
  
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), molecules = microns, assay = assay)
   
  return(coords)
}

stardist_Image_spatial <- function(
    object,
    image.dir,
    label,
    gdf,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice,
    filter.matrix = TRUE,
    image.type = "VisiumV2"
) {
  ##################################################
  #object = cellsegment
  #image.dir = "/home/wyt/R_project/FWSP20250659_1W2W3W/image"
  #label = "1W"
  #gdf = gdf_1W
  #image.name = "tissue_lowres_image.png"
  #assay = "Spatial"
  #slice = "1W"
  #filter.matrix = TRUE
  #image.type = "VisiumV2"
  
  # Validate the `image.type` parameter.
  
  image.type <- match.arg(image.type, choices = c("VisiumV1", "VisiumV2"))
  
  # Read in the H&E stain image.
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  # Read in the scale factors.
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )
  
  # Read in the tissue coordinates as a data.frame.
  cell_csv <- gdf %>% select(-geometry,-X) %>% rename(cell_id = id) %>%
    mutate(centroid = sub("POINT \\((.*)\\)", "\\1", centroid)) %>%
    separate(centroid, c('centroid_x', 'centroid_y'), sep = ' ') %>%
    mutate(
      centroid_x = as.numeric(centroid_x), 
      centroid_y = as.numeric(centroid_y),
      cell_id = paste0(label, "_", cell_id)
    )
  
  boundries <- gdf %>%  rename(cell_id = id) %>% 
    select(cell_id, geometry) %>% 
    mutate(
      geometry = geometry %>% 
        str_replace('.*\\((.*)\\){2}', '\\1') 
    ) %>% 
    separate_rows(geometry, sep = ',') %>%
    map(str_trim) %>%
    as.data.frame() %>% 
    separate(geometry, c('vertex_x', 'vertex_y'), sep = ' ') %>%
    mutate(
      vertex_x = vertex_x %>% as.numeric(), 
      vertex_y = vertex_y %>% as.numeric(),
      cell_id = paste0(label, "_", cell_id)
    )
  
  coordinates <- data.frame(x = cell_csv$centroid_x, 
                            y = cell_csv$centroid_y,row.names = cell_csv$cell_id,
                            cell = cell_csv$cell_id,
                            stringsAsFactors = FALSE)
  
  segmentations <- data.frame(x = boundries$vertex_x,y = boundries$vertex_y,
                              cell = boundries$cell_id)
  
  object$sample = as.character(object$sample)
  object = subset(object,sample %in% label)
  cells = colnames(object)
  indices  <- subset(coordinates,cell %in% cells)
  segmentations <- subset(segmentations,cell %in% cells)
  
  key <- Key(slice, quiet = TRUE)
  
  image_Y_min <- ceiling(min(indices$x)*scale.factors$lowres) 
  image_X_min <- ceiling(min(indices$y)*scale.factors$lowres) 
  image_Y_max <- ceiling(max(indices$x)*scale.factors$lowres) 
  image_X_max <- ceiling(max(indices$y)*scale.factors$lowres) 
  
  if (image_X_min < 0 ) {image_X_min = 0}
  if (image_Y_min < 0 ) {image_X_min = 0}
  if (image_X_max > dim(image)[1] ) {image_X_max = dim(image)[1]}
  if (image_Y_max > dim(image)[2] ) {image_Y_max = dim(image)[2]}
  sub_image <- image[image_X_min:image_X_max,image_Y_min:image_Y_max,]
  
  indices['imagecol'] <- indices['x'] - image_Y_min/scale.factors$lowres
  indices['imagerow'] <- indices['y'] - image_X_min/scale.factors$lowres
  segmentations['imagecol'] <- segmentations['x'] - image_Y_min/scale.factors$lowres
  segmentations['imagerow'] <- segmentations['y'] - image_X_min/scale.factors$lowres
  segmentations.data <- list(centroids = CreateCentroids(indices[, c("imagerow", "imagecol")]), 
                             segmentation = CreateSegmentation(segmentations))
  
  # Return the specified `image.type`.
  if (image.type == "VisiumV1") {
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = indices[colnames(spatial),],
      scale.factors = scale.factors,
      image = sub_image
    )
    
    # As of v5.1.0 `Radius.VisiumV1` no longer returns the value of the 
    # `spot.radius` slot and instead calculates the value on the fly, but we 
    # can populate the static slot in case it's depended on.
    visium.v1@spot.radius <- Radius(visium.v1)
    
    return(visium.v1)
  }
  
  # If `image.type` is not "VisiumV1" then it must be "VisiumV2".
  stopifnot(image.type == "VisiumV2")
  # Create an `sp` compatible `FOV` instanc

  fov <- CreateFOV(coords = segmentations.data, 
                      type = c("segmentation",  "centroids"), 
                      assay = assay,
                      radius = scale.factors[["spot"]],
                      key = key)
  
  # Build the final `VisiumV2` instance, essentially just adding `image` and
  # `scale.factors` to the `fov`.
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = sub_image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}

