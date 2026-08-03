library(tidyverse)
library(Seurat)
library(spacexr)

rctd_out_dir <- 'RCTD_annotation'
dir.create(rctd_out_dir)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")


all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data@meta.data %>% head()

sc_data <- readRDS('RCTD/SC2510275964.rds')
sc_data@meta.data %>% head()
sc_data$cell_type %>% table()
sc_data$cell_type_sub %>% table()




# annotation ------------
sc_data_anno <- as.factor(sc_data$cell_type)
names(sc_data_anno) <- rownames(sc_data@meta.data)

## make sc ref -----------------------------------------------------------------
sc_data_ref <- Reference(
  GetAssayData(sc_data, layer = 'counts'), 
  sc_data_anno, 
  sc_data$total_counts
)

## prepare spatial data --------------------------------------------------------
spatial_rna <- 
  SpatialRNA(
    all_data@meta.data[c('x', 'y')], 
    LayerData(all_data, layer = 'counts', assay = 'RNA'), 
    colSums(LayerData(all_data, layer = 'counts', assay = 'RNA'))
  )

## run RCTD --------------------------------------------------------------------
rctd_obj <- 
  create.RCTD(
    spatial_rna, 
    sc_data_ref, 
    max_cores = 10
  )

rctd_result <- 
  run.RCTD(
    rctd_obj,
    doublet_mode = 'doublet'
  )

saveRDS(rctd_result@results$results_df, file.path(rctd_out_dir, 'rctd_result_df.rds'))




## visualization ---------------------------------------------------------------
names(allcolour) <- sc_data$cell_type %>% unique()

meta <- rctd_result@results$results_df['first_type']
meta$first_type_annotation <- meta$first_type_annotation
all_data <- AddMetaData(all_data, metadata = meta['first_type_annotation'])
all_data$first_type_annotation <- as.character(all_data$first_type_annotation)
all_data$first_type_annotation[is.na(all_data$first_type_annotation)] <- "Unknown"
all_data %>% head()


# spatial plot #
all_data_meta_ls <- 
  all_data@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  select(cell_id, x, y, first_type_annotation, sample, group) %>% 
  group_by(sample) %>% 
  nest() %>% 
  deframe()

all_data$first_type_annotation %>% table()
names(allcolour) <- all_data$first_type_annotation %>% unique()
allcolour['Unknown'] <- 'lightgrey'
allcolour 
  
spatial_dimplot_ls <- 
  all_data_meta_ls %>% 
  imap(function(dat, sample){ # browser()
    dat %>% 
      
      ggplot(aes(x = x, y = y, color = first_type_annotation)) +
      geom_point(size = 0.6, stroke = 0.1)  +
      scale_y_reverse() +
      coord_equal() + 
      scale_color_manual(values = allcolour, guide = guide_legend(override.aes = list(size = 3))) + 
      theme_void() +
      labs(fill = NULL, title = paste(sample, dat$group, sep = "_"))
  })
  
spatial_dimplot_ls %>% 
  iwalk(function(spatial_dimplot, sample){
    pdf(file.path(rctd_out_dir, str_glue('spatial dimplot - {sample} - first_type.pdf')), width = 15, height = 15)
    print(spatial_dimplot)
    dev.off()
  })




## save --------------------------------------------------------------------------
all_data_meta <- 
  all_data@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  select(cell_id, x, y, first_type_annotation, sample, group)

saveRDS(all_data_meta, file.path(rctd_out_dir, "all_data_meta_RCTD_annotation.rds"))



# export csv
all_data_meta %>%
  write.csv(file = file.path(rctd_out_dir, "all_data_meta_RCTD_annotation.csv"), row.names = FALSE)




## Annotation info output -------------------------------------------------------
all_data$first_type_annotation %>% table()
cellType_anno <- 
  all_data@meta.data[c('project', 'first_type_annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*?_', ''), 
    group = first_type_annotation
  ) %>% 
  dplyr::select(cell_id, group)

cellType_anno %>% 
  write.csv(str_glue('RCTD/cell type_annotation_RCTD.csv'), row.names = F)




## cellType distribution ---------------------------------------------------------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = first_type)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(rctd_out_dir, "cell_freq_project_RCTD_annotation.pdf"),
  height = 6,
  width  = 6
)


p_bar_group  <-
  all_data_meta %>%
  group_by( first_type,project) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  tidyr::complete(first_type, project, fill = list(count = 0)) %>%
  mutate(freq = count/sum(count)) %>% 
  
  
  ggplot(aes(x = first_type, y = freq, fill = project)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_bar_group,
  filename = file.path(rctd_out_dir, "boxplot_cell_freq_RCTD_annotation.pdf"),
  height = 8,
  width  = 15
)