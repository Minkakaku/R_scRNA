library(Seurat)
library(tidyverse)


# 映射 -----------------
all_data <- readRDS(file = file.path("R_OUTPUT", "all_data.rds"))
gut <- readRDS("gut_annotation/gut.rds")
heart <- readRDS("heart_annotation/heart.rds")
kidney <- readRDS("kidney_annotation/kidney.rds")
liver <- readRDS("liver_annotation/liver.rds")
lung <- readRDS("lung_annotation/lung.rds")

library(tibble)
library(data.table)
subtype <- rbind(gut[['subtype']] %>% rownames_to_column('cells'), 
                 heart[['subtype']] %>% rownames_to_column('cells'), 
                 kidney[['subtype']] %>% rownames_to_column('cells'), 
                 liver[['subtype']] %>% rownames_to_column('cells'), 
                 lung[['subtype']] %>% rownames_to_column('cells'))
subtype %>% dim()
all_data

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = subtype
    )%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
all_data$annotation %>% table()
all_data$annotation2 %>% table()
all_data$annotation %>% table() %>% sum()
all_data$annotation2 %>% table() %>% sum()


saveRDS(all_data, "all_data_annotation.rds")




# spatial plot -----------------------------------------------------------------
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
gut_meta <- gut[[]] %>% rownames_to_column('cell_id')
heart_meta <- heart[[]] %>% rownames_to_column('cell_id')
kidney_meta <- kidney[[]] %>% rownames_to_column('cell_id')
liver_meta <- liver[[]] %>% rownames_to_column('cell_id')
lung_meta <- lung[[]] %>% rownames_to_column('cell_id')

## heart --------
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(heart_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(heart_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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


## liver --------
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(liver_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(liver_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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


## gut --------
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(gut_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(gut_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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


## lung ----------
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(lung_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(lung_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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


## kidney ---------- 
p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(kidney_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(kidney_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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
    if (!(sample %in% unique(gut_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(gut_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = subtype, group = cell_id)) + 
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




# split spatial plot ------------
## heart ----------
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(heart$sample)){ # browser()
      anno_data   <- heart_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(heart_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$subtype %>% levels(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$subtype %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$subtype <- as.character(dat2$subtype)
          dat2$subtype[dat2$subtype != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = subtype, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path('3. spatial plot', str_glue("spatial dim plot {name}_cellType - split.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


## liver ----------
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(liver$sample)){ # browser()
      anno_data   <- liver_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(liver_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$subtype %>% levels(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$subtype %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$subtype <- as.character(dat2$subtype)
          dat2$subtype[dat2$subtype != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = subtype, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path('3. spatial plot', str_glue("spatial dim plot {name}_cellType - split.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


## gut ----------
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(gut$sample)){ # browser()
      anno_data   <- gut_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(gut_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$subtype %>% levels(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$subtype %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$subtype <- as.character(dat2$subtype)
          dat2$subtype[dat2$subtype != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = subtype, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path('3. spatial plot', str_glue("spatial dim plot {name}_cellType - split.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


## lung ----------
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(lung$sample)){ # browser()
      anno_data   <- lung_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(lung_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$subtype %>% levels(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$subtype %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$subtype <- as.character(dat2$subtype)
          dat2$subtype[dat2$subtype != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = subtype, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path('3. spatial plot', str_glue("spatial dim plot {name}_cellType - split.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


## kidney ----------
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(kidney$sample)){ # browser()
      anno_data   <- kidney_meta[c("sample", "group")] %>% distinct()
      group_info  <- anno_data$group[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(kidney_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$subtype %>% levels(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$subtype %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$subtype <- as.character(dat2$subtype)
          dat2$subtype[dat2$subtype != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = subtype, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path('3. spatial plot', str_glue("spatial dim plot {name}_cellType - split.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })
