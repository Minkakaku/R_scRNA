library(tidyverse)
library(Seurat)


out_dir_260513 <- "20260513 dat_analysis"
dir.create(out_dir_260513)

# read data
all_data <- readRDS("all_data_full.rds")
sp_epi3 <- readRDS('sp_subcluster/sp_epi3.rds')

rctd_fibro <- readRDS("RCTD_fibro/rctd_result_ls.rds")
df_rctd_fibro <- rctd_fibro %>% map(~ .x@results$results_df) %>% bind_rows()

df_rctd_fibro$first_type %>% unique()

df_epi <- sp_epi3[['annotation_epi']] %>% rownames_to_column('cell_id')
df_fibro <- df_rctd_fibro %>% dplyr::select(annotation_fibro = first_type) %>% rownames_to_column('cell_id')

# merge annotation --------------
meta_anno_full <- 
  all_data[['annotation']] %>%
  rownames_to_column('cell_id') %>%
  left_join(df_epi, by = "cell_id") %>% 
  left_join(df_fibro, by = "cell_id")

meta_anno_full2 <-
  meta_anno_full %>% 
  map_dfc(as.character) %>%
  mutate(annotation3 = ifelse(
    is.na(annotation_epi),
    annotation,
    annotation_epi
  )) %>%
  mutate(
    annotation3 = ifelse(
      is.na(annotation_fibro),
      annotation3,
      annotation_fibro
    )
  ) %>% 
  mutate(annotation3 = str_replace(annotation3, "Others", "Stromal")) %>%
  dplyr::select(cell_id, annotation3) %>% 
  column_to_rownames('cell_id')


all_data <- all_data %>% AddMetaData(meta_anno_full2)

# saveRDS(all_data, 'all_data_full2.rds') ## not  run

# plot celltype ------------------
# spatial 
cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()
some_genes <- c()
annotation_col <- 'annotation3'
plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')


all_cells <- all_data[[annotation_col]] %>% unlist %>% unique
allcolour <- structure(
  names = all_cells,
  ggsci::pal_igv()(length(all_cells))
)
p_sp_ct_ls2 <- 
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width <- 5
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
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir_260513, "1. spatial plot of all celltypes in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls2 %>% walk(print)
dev.off()
