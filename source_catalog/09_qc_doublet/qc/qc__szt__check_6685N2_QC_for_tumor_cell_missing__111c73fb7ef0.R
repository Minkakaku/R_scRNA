library(Seurat)
library(tidyverse)

# qc metric for 06685N2
# sample,  nFeature_RNA, nCount_RNA, log10GenesPerUMI
# 06685N2, <100,         <200,       >0.8:<0.985

# some tumor of 6685N2 were kicked out

sc_6685n2 <- CreateSeuratObject(
  counts       = Read10X_h5("data/06685N2/cell_feature_matrix.h5")[[1]],
  min.cells    = 10,
  min.features = 0
)


sc_6685n2$log10GenesPerUMI <- log10(sc_6685n2$nFeature_RNA) / log10(sc_6685n2$nCount_RNA)


coords_6685N2 <- read.csv("data/06685N2/cells.csv")

# QC plot
VlnPlot(
  sc_6685n2, 
  features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"), 
  ncol = 4, 
  log = FALSE, 
  pt.size = 0,
  raster = TRUE
) 



# spatial plot for QC plot
# one of qc metric kicked out the tumor cells
# nFeature_RNA
sc_6685n2[[]] %>% 
  rownames_to_column("cell_id") %>%
  right_join(coords_6685N2[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  dplyr::filter(
    x_centroid < 2000,
    y_centroid < 3000
  ) %>%
  
  dplyr::filter(nFeature_RNA < 100) %>%
  
  ggplot(aes(x = x_centroid, y = y_centroid, color = nFeature_RNA)) +
  # ggplot(aes(x = x_centroid, y = y_centroid, color = nCount_RNA)) +
  # ggplot(aes(x = x_centroid, y = y_centroid, color = log10GenesPerUMI)) +
  scattermore::geom_scattermore() +
  scale_color_gradient(low = "gray", high = "red")

# nCount_RNA
sc_6685n2[[]] %>% 
  rownames_to_column("cell_id") %>%
  right_join(coords_6685N2[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  dplyr::filter(
    x_centroid < 2000,
    y_centroid < 3000
  ) %>%
  
  dplyr::filter(nCount_RNA < 200) %>%
  
  # ggplot(aes(x = x_centroid, y = y_centroid, color = nFeature_RNA)) +
  ggplot(aes(x = x_centroid, y = y_centroid, color = nCount_RNA)) +
  # ggplot(aes(x = x_centroid, y = y_centroid, color = log10GenesPerUMI)) +
  scattermore::geom_scattermore() +
  scale_color_gradient(low = "gray", high = "red")

# log10GenesPerUMI
sc_6685n2[[]] %>% 
  rownames_to_column("cell_id") %>%
  right_join(coords_6685N2[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  dplyr::filter(
    x_centroid < 2000,
    y_centroid < 3000
  ) %>%
  
  dplyr::filter(
    log10GenesPerUMI > 0.7, 
    log10GenesPerUMI < 0.985
  ) %>%


  # ggplot(aes(x = x_centroid, y = y_centroid, color = nFeature_RNA)) +
  # ggplot(aes(x = x_centroid, y = y_centroid, color = nCount_RNA)) +
  ggplot(aes(x = x_centroid, y = y_centroid, color = log10GenesPerUMI)) +
  scattermore::geom_scattermore() +
  scale_color_gradient(low = "gray", high = "red")


# result
# log10GenesPerUMI > 0.8 -> log10GenesPerUMI > 0.7


# read ROI cells
roi_6685n2 <-
  read.table(
    "6685N2_coordinates.csv", 
    skip = 2, 
    comment.char = "#", 
    sep = "\t", 
    header = TRUE
  )


sc_6685n2[, roi_6685n2$Cell.ID] %>% .[[]]


all_data[, str_c("3", roi_6685n2$Cell.ID, sep = "_")]






