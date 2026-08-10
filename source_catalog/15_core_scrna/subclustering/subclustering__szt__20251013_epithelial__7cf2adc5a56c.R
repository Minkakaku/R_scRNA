library(Seurat)
library(tidyverse)

sp_epi_dim20 <- readRDS("sp_subcluster/sp_epi-pariedsample-dim20.rds")

epi_out_dir <- "20251013 epi"
dir.create(epi_out_dir)


# de analysis RNA_snn_res.0.6 ---------
sp_epi_dim20 %>% DimPlot(group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap")
# sp_epi_dim20 %>% DimPlot(group.by = "annotation_epi", label = TRUE, reduction = "harmonyumap")


all_markers_all <- 
  FindAllMarkers(
    sp_epi_dim20,
    group.by = "RNA_snn_res.0.6",
    return.thresh = 1.1
  )

write.csv(
  all_markers_all,
  file = file.path(epi_out_dir, 'all markers - RNA_snn_res.0.6.csv')
)

# 2.4 de analysis of 9 12 13 and others ---------
sp_epi_dim20$epi_annotation <- 
  ifelse(
    sp_epi_dim20$RNA_snn_res.0.6 %in% c('9', '12', '13'),
    "cluster_9_12_13",
    "Others"
  )
sp_epi_dim20$epi_annotation %>% unique

all_markers_all2 <- 
  FindAllMarkers(
    sp_epi_dim20,
    group.by = "epi_annotation",
    return.thresh = 1.1
  )

write.csv(
  all_markers_all2,
  file = file.path(epi_out_dir, 'all markers - epi_annotation.csv')
)


# 2.5 de analysis of 9 12 13 and other clusters ---------
sp_epi_dim20$epi_annotation2 <- 
  ifelse(
    sp_epi_dim20$RNA_snn_res.0.6 %in% c('9', '12', '13'),
    "cluster_9_12_13",
    as.character(sp_epi_dim20$RNA_snn_res.0.6)
  )
sp_epi_dim20$epi_annotation2 %>% unique

all_markers_all3 <- 
  FindAllMarkers(
    sp_epi_dim20,
    group.by = "epi_annotation2",
    return.thresh = 1.1
  )

write.csv(
  all_markers_all3,
  file = file.path(epi_out_dir, 'all markers - epi_annotation2.csv')
)











