setwd('~/')

library(Seurat)
library(tidyverse)

out_dir <- 'CISS_result'
dir.create(out_dir)

heart_obj <- readRDS("heart_annotation/heart.rds")
liver_obj <- readRDS("liver_annotation/liver.rds")
gut_obj   <- readRDS("gut_annotation/gut.rds")
lung_obj  <- readRDS("lung_annotation/lung.rds")
kidney_obj<- readRDS("kidney_annotation/kidney.rds")

seurat_list <- list(
  heart = heart_obj,
  liver = liver_obj,
  gut   = gut_obj,
  lung  = lung_obj,
  kidney= kidney_obj
)

# apoptosis score ----------
inflammation <- c('Hmgb1', 'Nlrp3', 'Il1b', 'Ccl3', 'Tnf', 'Il6')
gene_ls <- list(inflammation)
names(gene_ls) <- c('inflammation')

tissues <- c("heart", "liver", "gut", "lung", "kidney")

seurat_ls <- 
  seurat_list %>% 
  imap(function(seurat_obj, tissue){ # browser()
    # 取出Seurat对象
    seurat_obj <- seurat_list[[tissue]]
    
    # 提取HS组
    if(tissue == "gut"){
      seurat_hs <- subset(seurat_obj, group %in% c("HS0_gut", "HS6_gut", "HS12_gut", "HS_gut"))
    } else {
      hs_group <- paste0("HS_", tissue)
      seurat_hs <- subset(seurat_obj, group == hs_group)
    }
    
    # 打分
    seurat_hs <- AddModuleScore(
      object = seurat_hs,
      features = gene_ls, 
      name = names(gene_ls)
    ) 
    
    seurat_hs
  })




# CISS annotation ------------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb


seurat_ls2 <- 
  seurat_ls %>% 
  imap(function(seurat_hs, tissue){ # browser()
    meta_CISS <- 
      seurat_hs@meta.data %>% 
      rownames_to_column('cell_id') %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(meta){
        meta$CISS <- ifelse(meta$inflammation1 > quantile(meta$inflammation1, probs = 0.9), TRUE, FALSE)
        # print(meta$DSS %>% table())
        meta
      }) %>% 
      bind_rows() %>% 
      column_to_rownames('cell_id')
    
    seurat_hs <- AddMetaData(seurat_hs, meta_CISS['CISS'])    
    
    seurat_hs
  })

seurat_ls2

seurat_ls2$heart@meta.data %>% head()
seurat_list$heart@meta.data %>% head()

seurat_ls2$heart$CISS %>% table()
seurat_ls2$liver$CISS %>% table()
seurat_ls2$gut$CISS %>% table()
seurat_ls2$lung$CISS %>% table()
seurat_ls2$kidney$CISS %>% table()

saveRDS(seurat_ls2, 'R_OUTPUT/seurat_CISS_HSgroup_ls.rds')




# density plot ----------
out_dir <- 'CISS_density_plot'
dir.create(out_dir)


all_data_meta <- 
  seurat_ls2 %>% 
  map(function(seurat_obj){
    seurat_obj@meta.data %>% 
      rownames_to_column('cell_id')
  }) %>% 
  bind_rows()

all_data_meta$CISS %>% table()

cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')

contour_plot_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'CISS', 'inflammation1')], by = 'cell_id') %>% 
      dplyr::rename(inflammation1 = inflammation1)
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(CISS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = inflammation1, color = CISS)) + 
      # geom_point(aes(color = CISS)) +
      # geom_tile(aes(fill = CISS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = inflammation1, color = CISS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_gradient2(
        low = "white",             
        high = "red",              # 高表达
        midpoint = 0,  # 中间值
        na.value = "lightgrey"     # NA值的颜色
      ) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'lightgrey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      geom_density2d(
        data = ~ dplyr::filter(., CISS == "TRUE"), 
        aes(color = CISS),
        contour_var = "density",
        show.legend = FALSE,
        # n = 250
        # h = 500
      ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

contour_plot_ls2 %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("contour plot - combine - {sample} - score.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })




# boxplot ---------
## boxplot / annotation --------------
boxplot_CISS_prop_annotation <- 
  all_data_meta %>% 
  dplyr::select(CISS, annotation, sample, group) %>% 
  filter(CISS == 'TRUE') %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(annotation) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = annotation, y = freq, fill = annotation)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of CISS - {group}'))
  }) %>% 
  patchwork::wrap_plots(n = 3)

ggsave(boxplot_CISS_prop_annotation, filename = file.path(out_dir, 'boxplot_CISS_prop_by_annotation.pdf'), width = 15, height = 15)  


## boxplot / cell_type --------------
boxplot_CISS_prop_cell_type <- 
  all_data_meta %>% 
  dplyr::select(CISS, subtype, sample, group) %>% 
  dplyr::rename(cell_type = subtype) %>% 
  filter(CISS == 'TRUE') %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(cell_type) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = cell_type, y = freq, fill = cell_type)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of CISS - {group}'))
  }) %>% 
  patchwork::wrap_plots(n = 3)

ggsave(boxplot_CISS_prop_cell_type, filename = file.path(out_dir, 'boxplot_CISS_prop_by_cell_type.pdf'), width = 15, height = 15)  


## boxplot / subtype --------------
boxplot_CISS_prop_subtype <- 
  all_data_meta %>% 
  dplyr::select(CISS, subtype_re, sample, group) %>% 
  dplyr::rename(subtype = subtype_re) %>% 
  filter(CISS == 'TRUE') %>% 
  filter(!is.na(subtype)) %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){ # browser()
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(subtype) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = subtype, y = freq, fill = subtype)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of CISS - {group}'))
  })

ggsave(boxplot_CISS_prop_subtype$HS_liver, filename = file.path(out_dir, 'boxplot_CISS_prop_by_subtype.pdf'), width = 5, height = 5)  




# CISS spatial plot -----------
## spatial plot / cell_type -------------
CISS_spatial_plot_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'CISS', 'subtype')], by = 'cell_id') %>% 
      mutate(
        cell_type = ifelse(CISS == 'TRUE', as.character(subtype), 'non-CISS'), 
        cell_type = factor(cell_type, levels = c(unique(as.character(subtype)), 'non-CISS'))
      )
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(CISS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = cell_type, color = CISS)) + 
      # geom_point(aes(color = CISS)) +
      # geom_tile(aes(fill = CISS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = cell_type, color = CISS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_manual(values = c(allcolour[1:length(levels(dat2$cell_type)) - 1], 'lightgrey'), drop = FALSE) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'grey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      # geom_density2d(
      #   data = ~ dplyr::filter(., CISS == "TRUE"), 
      #   aes(color = CISS),
      #   contour_var = "density",
      #   show.legend = FALSE,
      #   # n = 250
      #   # h = 500
      # ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

CISS_spatial_plot_ls %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("CISS - {sample} - cell_type.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })


## spatial plot / subtype -------------
CISS_spatial_plot_subtype_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'CISS', 'subtype_re')], by = 'cell_id') %>% 
      filter(!is.na(subtype_re)) %>% 
      mutate(
        subtype = 
          ifelse(CISS == 'TRUE', as.character(subtype_re), 'non-CISS') %>% 
          str_replace_all(c('Netrophil' = 'Neutrophil')), 
        subtype = factor(subtype, levels = c('B cell', 'Cholangiocyte', 'Endothelial', 'Fibroblast', 'Macrophage', 
                                             'Neutrophil', 'Pericentral-hepatocyte', 'Periportal-hepatocyte', 
                                             'Pericyte', 'non-CISS'))
      )
    
    if (dim(dat)[1] == 0) return(NULL)
    
    # browser()
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(CISS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = subtype, color = CISS)) + 
      # geom_point(aes(color = CISS)) +
      # geom_tile(aes(fill = CISS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = subtype, color = CISS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_manual(values = c(allcolour[1:9], 'lightgrey'), drop = FALSE) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'grey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      # geom_density2d(
      #   data = ~ dplyr::filter(., CISS == "TRUE"), 
      #   aes(color = CISS),
      #   contour_var = "density",
      #   show.legend = FALSE,
      #   # n = 250
      #   # h = 500
      # ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

CISS_spatial_plot_subtype_ls %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("CISS - {sample} - subtype.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })
