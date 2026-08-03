library(Seurat)
library(tidyverse)

all_data_full <- readRDS("all_data_full.rds")

all_data_sc <- readRDS('output/singlecell/anno_epi.rds')
all_data <- readRDS('sp_subcluster/sp_epi_used.rds')


label_transfer_dir <- "epi_label_transfer"
dir.create(label_transfer_dir)

# all_data_sc 
# all_data_sc_used <- CreateSeuratObject(
#   CreateAssay5Object(
#     counts = all_data_sc@assays$RNA$counts,
#     data = all_data_sc@assays$RNA$data
#   ),
#   meta.data = all_data_sc@meta.data
# )
all_data_sc_used <- all_data_sc


all_data_sc_used <-
  all_data_sc_used %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "nCount_RNA") %>%
  RunPCA()

# all_data
# all_data_used <- CreateSeuratObject(
#   CreateAssay5Object(
#     counts = all_data@assays$RNA$counts,
#     data = all_data@assays$RNA$data
#   ),
#   meta.data = all_data@meta.data
# )
# all_data_used[['pca']] <- all_data[['pca']] 
all_data_used <- all_data 


# label transfer ------------------

ElbowPlot(all_data_sc_used, ndims = 50)
ElbowPlot(all_data_used, ndims = 50)


all_data.anchors <- 
  FindTransferAnchors(
    reference = all_data_sc_used,
    query = all_data_used,
    # normalization.method = "SCT",
    dims = 1:20,
    reference.reduction = "pca"
  )

predictions <- 
  TransferData(
    anchorset = all_data.anchors, 
    refdata = all_data_sc$annotation, 
    dims = 1:20
  )


all_data <- 
  AddMetaData(
    all_data, 
    metadata = predictions
  )

# saveRDS(all_data, 'sc_epi_label_transfer.rds')

pdf(
  file = file.path(label_transfer_dir, "label transfer using seurat.pdf"),
  height = 6,
  width = 8
)
all_data %>%
  DimPlot(
    reduction = "harmonyumap",
    group.by = "predicted.id",
    label = TRUE
  ) +
  ggsci::scale_color_igv()

# all_data %>%
#   DimPlot(
#     reduction = "harmonyumap",
#     group.by = "annotation_full",
#     label = TRUE
#   ) +
#   ggsci::scale_color_igv()
dev.off()


# pdf(
#   file = "label transfer using seurat - dotplot.pdf",
#   height = 6,
#   width = 10
# )
# DotPlot(
#   all_data_sc,
#   features = c(
#     'CD3D', 'CD3E',  "CD4", "IL7R", "CD8A", "FCGR3A", "NCAM1", "FOXP3",
#     'VCAN', "FCN1", "CD14", "CD68", "C1QA",
#     "EPCAM", "PECAM1", "COL1A1",
#     "CLEA9A", "IDO1", "CLEC10A", "CD1C",
#     "S100A8", "S100A9", "CSF3R", 
#     "TPSAB1",
#     "MS4A1", "JCHAIN"
#     
#   ),
#   group.by = "CellType"
# ) +
#   ggtitle("single cell") +
#   RotatedAxis()
# 
# 
# DotPlot(
#   all_data,
#   features = c(
#     'CD3D', 'CD3E', "CD4", "IL7R", "CD8A", "FCGR3A", "NCAM1", "FOXP3",
#     'VCAN', "FCN1", "CD14", "CD68", "C1QA",
#     "EPCAM", "PECAM1", "COL1A1",
#     "CLEA9A", "IDO1", "CLEC10A", "CD1C",
#     "S100A8", "S100A9", "CSF3R", 
#     "TPSAB1",
#     "MS4A1", "JCHAIN"
#     
#   ),
#   group.by = "predicted.id"
# ) +
#   ggtitle("label transfer") +
#   RotatedAxis()
# 
# dev.off()



# spatial plot of all celltypes - annotation --------------
cells_position_ls <-
  cbind(
    all_data_full[['sample']],
    Embeddings(all_data_full[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

annotation_col3 <- 'predicted.id'
plot_dat_sub3 <- 
  FetchData(
    all_data, 
    vars = c(annotation_col3)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls3 <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords # %>% left_join(plot_dat_sub, by = "cell_id")
    dat3 <- coords %>% inner_join(plot_dat_sub3, by = "cell_id")
    
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        height = height, 
        width = width,
        fill = "gray90"
      ) +
      geom_tile(
        data = dat3,
        aes(fill = .data[[annotation_col3]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(label_transfer_dir, "spatial plot of epi annotation in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls3 %>% walk(print)
dev.off()



