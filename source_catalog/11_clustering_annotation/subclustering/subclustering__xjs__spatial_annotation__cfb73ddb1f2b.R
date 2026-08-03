library(Seurat)
library(tidyverse)


# 映射 -----------------
all_data <- readRDS('R_OUTPUT/all_data.rds')
All_T_SQ <- readRDS('T_SQ_subtype/All_T.rds')
All_T_SH <- readRDS('T_SH_subtype/All_T.rds')
Macrophage_SQ <- readRDS('Macrophage_SQ_subtype/Macrophage.rds')
Macrophage_SH <- readRDS('Macrophage_SH_subtype/Macrophage.rds')
Epithelial_SQ <- readRDS('Epithelial_SQ_subtype/Epithelial2.rds')
Epithelial_SH <- readRDS('Epithelial_SH_subtype/Epithelial.rds')

library(tibble)
library(data.table)
subtype <- rbind(All_T_SQ[['subtype']] %>% rownames_to_column('cells'), 
                 All_T_SH[['subtype']] %>% rownames_to_column('cells'), 
                 Macrophage_SQ[['subtype']] %>% rownames_to_column('cells'), 
                 Macrophage_SH[['subtype']] %>% rownames_to_column('cells'), 
                 Epithelial_SQ[['subtype']] %>% rownames_to_column('cells'), 
                 Epithelial_SH[['subtype']] %>% rownames_to_column('cells'))
all_data$annotation %>% table()

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_full = ifelse((as.character(annotation) %in% c('T cell', 'Macrophage', 'Epithelial/Tumor cell')),
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
  'CD4 T cell', 
  'low-nCounts', 
  'mix', 
  'Mix'
)
cellType_to_keep <- cellType_all[!cellType_all %in% cellType_to_rm]
all_data <- subset(all_data, annotation_full %in% cellType_to_keep)
saveRDS(all_data, "all_data_annotation_filtered.rds")




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
  })

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


all_data$annotation_full %>% table()
all_data$annotation4 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c('REST_tumor cell', 'Tumor cell')]
)
all_data$annotation4 %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


## highlight -----------
### Tumor ----------
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation4, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey","#DC143C","#40E0D0")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_tumor highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  


### CD8 Tem/ REST Tumor/ M1 ---------
all_data$annotation_full %>% table()
all_data$annotation5 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c('CD8 Tem', 'REST_tumor cell', 'M1 macrophage')]
)
all_data$annotation5 %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')

p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation5, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey","#228B22","#40E0D0","#DC143C")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_CD8 Tem - REST_tumor cell - M1 macrophage_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 


### CD8 Tem/ Tumor/ M1 ---------
all_data$annotation_full %>% table()
all_data$annotation6 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c('CD8 Tem', 'Tumor cell', 'M1 macrophage', 'M2 macrophage')]
)
all_data$annotation6 %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')

p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation6, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey","#228B22","#40E0D0","#FA8072","#DC143C")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_CD8 Tem - Tumor cell - M1 macrophage - M2 macrophage_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 


 



### target gene + celltype ------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

target_gene_ls <- c(
  'SLC6A4', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR1E', 'HTR2A', 'HTR2B', 'HTR2C', 'HTR3A', 'HTR4', 'HTR5A', 'HTR6', 'HTR7', 'HTRA1', 'HTRA2', 
  'SLC18A2', 'MAOA', 'MAOB', 'TPH1', 'TPH2'
)
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data)


#### REST Tumor -----------------
all_data$annotation_full %>% table()
all_data$annotation7 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c('REST_tumor cell', 'M1 macrophage', 'M2 macrophage')]
)
all_data$annotation7 %>% table()

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')


p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        # if (!x %in% colnames(all_data_meta)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = annotation7, group = cell_id), 
            size = 0.3,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", "#40E0D0","#FA8072","#DC143C")) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.5,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_target gene - REST_tumor cell - M1 - M2_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 


#### Tumor -----------------
all_data$annotation_full %>% table()
all_data$annotation8 <- fct_collapse(
  all_data$annotation_full, 
  'Others' = unique(all_data$annotation_full)[!unique(all_data$annotation_full) %in% c('Tumor cell', 'M1 macrophage', 'M2 macrophage')]
)
all_data$annotation8 %>% table()

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')


p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        # if (!x %in% colnames(all_data_meta)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = annotation8, group = cell_id), 
            size = 0.3,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", "#40E0D0","#FA8072","#DC143C")) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.5,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_target gene - Tumor cell - M1 - M2_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 




## Annotation info output -------------------------------------------------------
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
