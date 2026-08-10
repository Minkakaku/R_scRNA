library(Seurat)
library(dplyr)
library(tibble)
library(tidyverse)
library(openxlsx)


Anno_Seuratobject <- function(object=pbmc_filt, Anno_mod="cluster", Anno_path, Group_path=NULL) {
  if (Anno_mod == "cluster") {
    new.idends <- read.xlsx(Anno_path) %>%
    dplyr::select(1,2) %>%
    dplyr::rename(seurat_clusters=1, Celltype=2) 
    new.idends$seurat_clusters <- as.character(new.idends$seurat_clusters)
    object$seurat_clusters <- as.character(object$seurat_clusters)
    object@meta.data <- object@meta.data %>%  tibble::rownames_to_column("barcodes")  %>% 
      left_join(., new.idends, by="seurat_clusters")   %>% column_to_rownames("barcodes")
    Idents(object) <- "Celltype"
  } else if (Anno_mod == "cells") {
    new.idends <- read.xlsx({Anno_path}) %>%
      dplyr::select(1,2) %>%
      dplyr::rename(barcodes=1, Celltype=2)
    
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      left_join(., new.idends, by="barcodes") %>%
      column_to_rownames("barcodes")
    Idents(object) <- "Celltype"
    cells_id <- rownames(object@meta.data)[which(!is.na(object@meta.data$Celltype))] 
    if (length(cells_id) > 0) {
      object <- subset(object, cells = cells_id)
    } else {
      message("No cells with non-NA celltype found.")
    }
  }
  
  # 处理 Group_path 为空的情况
  if (!is.null(Group_path)) {
    Group <- read.xlsx(Group_path) %>%
      dplyr::select(1,2) %>%
      dplyr::rename(orig.ident=1, Group=2)
    
    object@meta.data <- object@meta.data %>%  tibble::rownames_to_column("barcodes")  %>% 
      left_join(., Group, by="orig.ident")   %>% column_to_rownames("barcodes")
  }
  
  return(object)
}
