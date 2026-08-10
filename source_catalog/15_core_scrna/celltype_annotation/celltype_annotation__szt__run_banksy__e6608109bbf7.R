library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)


# HD示例数据：/mnt/NexData_file/home/ztt/test_data/Visium_HD/Colorectal_Cancer/output/data

out_dir <- "banksy"
dir.create(out_dir)


# meta_data_anno <- read.csv('adata_celltype_obs.csv')

all_data <- readRDS("sp_fibro_rctd_anno.rds")

all_data2 <- all_data[, all_data$project == "A02994A6"]
all_data2 <- CreateSeuratObject(
  CreateAssay5Object(
    counts = all_data2@assays$RNA$counts.A02994A6.1
  ),
  meta.data = all_data2@meta.data
)


# all_data$x <- Embeddings(all_data[['spatial']])[, 'SPATIAL_1'] 
# all_data$y <- Embeddings(all_data[['spatial']])[, 'SPATIAL_2'] 

# # coords 平移
# rng_coord <-
#   all_data[[]] %>%
#   group_by(sample) %>%
#   mutate(
#     xmin = min(x),
#     xmax = max(x),
#     rng = xmax - xmin
#   )
# 
# all_data$xx <- all_data$x - rng_coord$xmin + as.numeric(factor(all_data$sample)) * 1.2 * max(rng_coord$rng)
# 
# data.frame(
#   all_data$xx,
#   all_data$y
# ) %>%
#   slice_sample(n = 10000) %>%
#   plot

# banksy ----------
all_data2 <- 
  all_data2 %>% 
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) 


all_data2 <- 
  RunBanksy(all_data2,
            lambda = 0.8, 
            verbose = TRUE,
            assay = "RNA", # "Spatial.008um", 
            slot = "data",
            features = "variable",
            k_geom = 50,
            
            dimx = "x", # use coord 
            dimy = "y",
            
            group = "roi" # <-----
  )

DefaultAssay(all_data2) <- "BANKSY"

# data_own <- SeuratWrappers:::get_data(all_data2, assay = "RNA", slot = 'data', features = 'variable', verbose = F)
# locs <- SeuratWrappers:::get_locs(all_data2, dimx = "x", dimy = "y", dimz = NULL,
#                                   ndim = 2, data_own = data_own, group = 'sample', verbose = F)
# 
# 
# locs %>%
#   slice_sample(n = 10000) %>%
#   plot





# Run PCA and UMAP -----------
all_data2 <- RunPCA(
  all_data2, 
  assay = 'BANKSY', 
  features = rownames(all_data2[['BANKSY']]), 
  npcs = 30
)

# umap
ElbowPlot(all_data2, ndims = 30)
all_data2 <- RunUMAP(all_data2, dims = 1:10)


# harmony
all_data2 <- all_data2 %>%
  harmony::RunHarmony(
    dims.use = 1:10,
    group.by.vars = "roi",
    reduction.save = "harmony"
  )

all_data2 <- RunUMAP(
  all_data2,
  dims = 1:10,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)


DimPlot(
  all_data2,
  group.by = "roi"
)


# Clustering
all_data2 <- FindNeighbors(all_data2, dims = 1:10, reduction = "harmony")
all_data2 <- FindClusters(all_data2, resolution = c(0.2, 0.4, 0.6, 0.8))



DimPlot(
  all_data2,
  group.by = "BANKSY_snn_res.0.2"
)


cluster_lt_1000 <- all_data2$BANKSY_snn_res.0.2 %>% table %>% {names(.)[. < 1000]}

all_data22 <- all_data2[, !all_data2$BANKSY_snn_res.0.2 %in% cluster_lt_1000]



DimPlot(
  all_data22,
  group.by = "BANKSY_snn_res.0.2"
)

DimPlot(
  all_data22,
  group.by = "roi"
)


saveRDS(all_data2, 'all_data2_banksy.rds')


# plot  ---------------
all_data2_banksy <- readRDS('all_data2_banksy.rds')
# export pdf
pdf(
  file.path(out_dir, "umap plot of banksy.pdf"),
  height = 6,
  width = 8
)
DimPlot(
  all_data2_banksy,
  group.by = "BANKSY_snn_res.0.2"
)
DimPlot(
  all_data2_banksy,
  group.by = "sample"
)
dev.off()


# spatial plot
annotation_col <- c("BANKSY_snn_res.0.2")

cells_position_ls <-
  all_data2_banksy[[c("x", "y", "sample")]] %>%
  rownames_to_column('cell_id') %>%
  dplyr::rename(
    SPATIAL_1 = x,
    SPATIAL_2 = y
  ) %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

plot_dat <-
  FetchData(
    all_data2_banksy, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')


p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 15
    width <- 15
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col]]),
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
  file.path(out_dir, "spatial plot of celltype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()


