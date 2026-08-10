library(Seurat)
library(tidyverse)


# 映射 -----------------
all_data <- readRDS('R_OUTPUT/all_data.rds')
all_data <- subset(JoinLayers(all_data), sample != 'D2')
Epithelial <- readRDS("Epithelial_subtype_remove_D2/Epithelial.rds")
Epithelial$subtype %>% table()
Fibroblast <- readRDS("Fibroblast_subtype_remove_D2/Fibroblast.rds")
Fibroblast$subtype %>% table()

library(tibble)
library(data.table)
subtype <- rbind(Epithelial[['subtype']] %>% rownames_to_column('cells'), 
                 Fibroblast[['subtype']] %>% rownames_to_column('cells'))
all_data$annotation %>% table()

meta.data.t <- 
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_full = ifelse((as.character(annotation) %in% c('Epithelial', 'Fibroblast')),
                             as.character(subtype), 
                             as.character(annotation)))%>%
  dplyr::select(cells, annotation_full) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation_full'])
all_data$annotation %>% table()
all_data$annotation_full %>% table()
all_data$annotation %>% table() %>% sum()
all_data$annotation_full %>% table() %>% sum()


cellType_all <- all_data$annotation_full %>% unique()
cellType_to_rm <- c(
  'low-nCounts', 
  'mix', 
  'noise'
)
cellType_to_keep <- cellType_all[!cellType_all %in% cellType_to_rm]
all_data <- subset(all_data, annotation_full %in% cellType_to_keep)
saveRDS(all_data, "all_data_annotation_removeD2_filtered.rds")




## spatial plot -----------------------------------------------------------------
spatial_theme <- 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
        #axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.title = element_text(hjust = 0.5, size = 20))


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')

cells_boundary_ls$D2 <- NULL


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_cellType.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_full, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_subType.pdf")), 
        height = 10, 
        width = 15)
    print(p)
    dev.off()
  })  




## Annotation info output -------------------------------------------------------
dir.create('anno_info')
all_data$annotation %>% table()
cellType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

cellType_anno_ls %>% 
  iwalk(function(cellType, project){
    write.csv(cellType, str_glue('anno_info/cell type {project}.csv'), row.names = F)
  })

all_data$annotation_full %>% table()
subType_anno_ls <- 
  all_data@meta.data[c('project', 'annotation_full')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation_full
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

subType_anno_ls %>% 
  iwalk(function(subType, project){
    write.csv(subType, str_glue('anno_info/sub cell type {project}.csv'), row.names = F)
  })
