
library(Seurat)
library(tidyverse)
library(ggsci)


all_data <- readRDS("all_data_ROI.rds")
cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


# 1.RCTD spatial plot --------------
results_outdir <- "20260630-modify"
spatial_outdir <- file.path(results_outdir, "1.spatial_plot")
dir.create(spatial_outdir, recursive = T)

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id")

cols <- pal_igv("default")(51)
celltype <- sort(unique(all_data$RCTD))

plan2cell <- cols[1:length(unique(all_data$RCTD))]
plan2cell <- setNames(plan2cell, celltype)


library(tibble)
cells_boundary_df <- bind_rows(
  cells_boundary_ls, .id = "sample"
)

pdf(file.path(spatial_outdir, "spatial plot of RCTD.pdf"), width = 20, height = 20)
cells_boundary_df %>% 
  right_join(all_data_meta, by = "cell_id") %>% 
  dplyr::filter(!is.na(first_type)) %>% 
  dplyr::filter(first_type != "Unknown") %>% 
  mutate(RCTD = first_type) %>% 
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(fill = RCTD, group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  scale_fill_manual(values = plan2cell) +
  theme_void()
dev.off()



major_marker_ls <- c(
  "B cell"        = "CD19",
  "B cell"        = "CD79A",
  "B cell"        = "MS4A1",
  'Endothelial'   = 'PECAM1', 
  # 'Endothelial' = 'VWF', 
  "Fibroblast"    = "PDGFRA", 
  "Fibroblast"    = "CCDC80",
  "Fibroblast"    = "SFRP4",
  # "Fibroblast"  = "COL1A1",
  # "Fibroblast"  = "COL1A2",
  # "Fibroblast"  = "DCN",
  # "Fibroblast"  = "ACTA2",
  "Monocyte"      = "CD14", 
  "Monocyte"      = "FCGR3A",
  "Macrophage"    = "CD68",
  "Macrophage"    = "CD163",
  # "Macrophage"  = "C1QA",
  # "Macrophage"  = "C1QC", 
  'Macrophage'    = 'CSF1R',
  # "Mast cell"   = "CPA3",
  "Mast cell"     = "CTSG",
  # "Mast cell"   = "TPSAB1",
  "Mast cell"     = "KIT",
  'Netrophil'     = "CSF3R",
  # 'Netrophil'   = "G0S2",
  # 'Netrophil'   = "RETNLG",
  # 'T cell'      = 'CD3D', 
  'T cell'        = 'CD3E', 
  'T cell'        = 'CD3G',
  # 'T cell'      = 'TRAC'
  'Epithelial'    = 'EPCAM'
  # 'Epithelial'  = 'KRT5',
  # 'Epithelial'  = 'KRT14',
)
plot_markers <- c("CD19", "CD79A", "MS4A1", 
                  "PECAM1", 
                  "PDGFRA", "CCDC80", "SFRP4", 
                  "CD14", "FCGR3A", 
                  "CD68", "CD163", "CSF1R", 
                  "CTSG", "KIT", 
                  "CSF3R", 
                  "CD3E", "CD3G", 
                  "EPCAM")
all_data_filt <- all_data[, all_data@meta.data$RCTD != "Unknown"]
Idents(all_data_filt) <- "RCTD"
cell_order <- c("B cell", "Endothelium", "Fibroblast", "Monocyte", "Macrophage", "Mast cell", "Neutrophil", "T cell", "Normal", "Tumor")
Idents(all_data_filt) <- 
  factor(all_data_filt@meta.data$RCTD, levels = cell_order)

pdf(file.path(spatial_outdir, "dotplot of rctd annotation.pdf"), width = 12, height = 6)
DotPlot(all_data_filt, features = plot_markers, dot.scale = 10, cols = c("lightgrey", "red"), col.min = 0) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()





# 2.cell ratio barplot by group ----------
cell_ratio_outdir <- file.path(results_outdir, "2.Cell_Ratio")
dir.create(cell_ratio_outdir)

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  filter(RCTD != "Unknown")

p_col_proj <- 
  all_data_meta %>% 
  group_by(group, RCTD) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  ggplot(aes(x = group, y = freq, fill = RCTD)) + 
  geom_col() + 
  # ggsci::scale_fill_igv() +
  scale_fill_manual(values = plan2cell) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <- 
  all_data_meta %>% 
  group_by(group, sample, RCTD) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  ggplot(aes(x = RCTD, y = freq, fill = group)) + 
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) + 
  ggsci::scale_fill_igv() +
  # scale_fill_manual(values = plan2cell) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj, 
  filename = file.path(cell_ratio_outdir, "cell_freq_project.pdf"), 
  height = 6, 
  width  = 5
)

ggsave(
  p_box_group, 
  filename = file.path(cell_ratio_outdir, "cell_freq_group.pdf"), 
  height = 6, 
  width  = 10
)





# 3.CN barplot by group --------
CN_outdir <- file.path(results_outdir, "3.CN analysis")
dir.create(CN_outdir)

All_type_meta_cn10 <- readRDS("RCTD_new/CN analysis/all_data_meta_cn10.rds")

All_type_cn10_count_df <- 
  All_type_meta_cn10 %>% 
  # group_by(sample, group, CN) %>% 
  group_by(group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))
All_type_cn10_count_df$CN <- factor(All_type_cn10_count_df$CN, levels = paste0("CN", 1:10))

p_cn10_props <- 
  All_type_cn10_count_df %>% 
  
  ggplot(aes(x = group, y = freq, fill = CN)) + 
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )


pdf(
  file.path(CN_outdir, "CN10_boxplot_and_barplot.pdf"), 
  height = 6, 
  width  = 5
)
p_cn10_props %>% print
# p_cn10_props_group %>% print
dev.off()



All_type_meta_cn15 <- readRDS("RCTD_new/CN analysis/all_data_meta_cn15.rds")
All_type_cn15_count_df <- 
  All_type_meta_cn15 %>% 
  # group_by(sample, group, CN) %>% 
  group_by(group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))
All_type_cn15_count_df $CN <- factor(All_type_cn15_count_df $CN, levels = paste0("CN", 1:15))

p_cn15_props <- 
  All_type_cn15_count_df %>% 
  
  ggplot(aes(x = group, y = freq, fill = CN)) + 
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )


pdf(
  file.path(CN_outdir, "CN15_boxplot_and_barplot.pdf"), 
  height = 6, 
  width  = 5
)
p_cn15_props %>% print
# p_cn15_props_group %>% print
dev.off()





# 4.tumor subcluster spatial plot ---------
Tumor_outdir <- file.path(results_outdir, "4.Tumor sub plot")
dir.create(Tumor_outdir)

colour = c('#C6DADF', '#E26EB8', '#918FBE', '#A2C5A6', '#EED784', '#D4B796', 
           '#F5BDC8', '#F0BDD8', '#77A1EC', '#F8D39C', '#F17591', '#CCC9E6')

Tumor <- readRDS("RCTD_new/Tumor/Tumor.rds")

pdf(file.path(Tumor_outdir, "DimPlot harmony UMAP.pdf"), height = 6, width = 7)
DimPlot(Tumor, group.by = "subtype", label = TRUE, reduction = "harmonyumap", repel = TRUE, cols = colour) + NoLegend()
dev.off()


Tumor_meta <- 
  Tumor[[]] %>% 
  rownames_to_column("cell_id")

p_col_proj <- 
  Tumor_meta %>% 
  group_by(group, subtype) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() + 
  # ggsci::scale_fill_igv() + 
  scale_fill_manual(values = colour) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj, 
  filename = file.path(Tumor_outdir, "cell_freq_project.pdf"), 
  height = 6, 
  width  = 5
)


## spatial plot 
celltype <- sort(unique(Tumor$subtype))
plan2cell <- colour[1:length(unique(Tumor_meta$subtype))]
plan2cell <- setNames(plan2cell, celltype)

library(tibble)
cells_boundary_df <- bind_rows(
  cells_boundary_ls, .id = "sample"
)

p_combined <- cells_boundary_df %>% 
  right_join(Tumor_meta, by = "cell_id") %>% 
  dplyr::filter(!is.na(subtype)) %>% 
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(fill = subtype, group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  scale_fill_manual(values = plan2cell) +
  theme_void()

ggsave(file.path(Tumor_outdir, "Tumor subcluster spatial plot combined.pdf"), p_combined, width = 20, height = 20)





# 5.Tumor c1 c4 CN plot ---------
CN_outdir <- file.path(results_outdir, "5.Tumor sub CN analysis")
dir.create(CN_outdir)

All_type_meta_cn10 <- readRDS("2025-11-04-Spatial analysis/CN analysis/all_data_meta_cn10.rds")

All_type_cn10_count_df <- 
  All_type_meta_cn10 %>% 
  # group_by(sample, group, CN) %>% 
  group_by(group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))
All_type_cn10_count_df$CN <- factor(All_type_cn10_count_df$CN, levels = paste0("CN", 1:10))

p_cn10_props <- 
  All_type_cn10_count_df %>% 
  
  ggplot(aes(x = group, y = freq, fill = CN)) + 
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )


pdf(
  file.path(CN_outdir, "CN10_boxplot_and_barplot.pdf"), 
  height = 6, 
  width  = 5
)
p_cn10_props %>% print
# p_cn10_props_group %>% print
dev.off()



All_type_meta_cn15 <- readRDS("2025-11-04-Spatial analysis/CN analysis/all_data_meta_cn15.rds")
All_type_cn15_count_df <- 
  All_type_meta_cn15 %>% 
  # group_by(sample, group, CN) %>% 
  group_by(group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))
All_type_cn15_count_df $CN <- factor(All_type_cn15_count_df $CN, levels = paste0("CN", 1:15))

p_cn15_props <- 
  All_type_cn15_count_df %>% 
  
  ggplot(aes(x = group, y = freq, fill = CN)) + 
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )


pdf(
  file.path(CN_outdir, "CN15_boxplot_and_barplot.pdf"), 
  height = 6, 
  width  = 5
)
p_cn15_props %>% print
# p_cn15_props_group %>% print
dev.off()





# 6.tumor infiltration ----------
tumor_infiltration_outdir <- file.path(results_outdir, "6.Tumor infiltration")
dir.create(tumor_infiltration_outdir)

library(sp)

sample_df <- read.csv("20260429-tumor infiltration/sample_group.csv")

reduceFindNN <- function(dat, k = 10, query = NULL, searchtype = "radius", radius = 100, mpp = 1, max_k = 2000){
  if(is.null(query)){
    query <- dat
  }
  
  stopifnot(all(c('x', 'y') %in% colnames(dat)))
  stopifnot(all(c('x', 'y') %in% colnames(query)))
  
  dat   <- dat[c('x', 'y')]
  query <- query[c('x', 'y')]
  
  nn_dat <- RANN::nn2(
    dat, 
    query = query,
    k = k,
    radius = radius / mpp # radius, um -> pixel
  )
  
  k_enough <- sum(apply(nn_dat$nn.idx == 0, 1, any))/nrow(nn_dat$nn.idx)
  cat(str_glue("k is {k} and is enough for {k_enough * 100}% cells.\n\n"))
  
  while(k < max_k && k_enough < 0.99){
    
    if(k < 700){
      k <- k + 100
    }else if(k < 2000){
      k <- k + 200
    }else if(k < 5000){
      k <- k + 500
    }else{
      k <- k + 1000
    }
    
    nn_dat <- RANN::nn2(
      dat, 
      query = query,
      searchtype = searchtype,
      k = k,
      radius = radius / mpp # radius, um -> pixel
    )
    
    k_enough <- sum(apply(nn_dat$nn.idx == 0, 1, any))/nrow(nn_dat$nn.idx)
    cat(str_glue("k is {k} and is enough for {k_enough * 100}% cells.\n\n"))
    
  }
  return(nn_dat) 
}

#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}

## coords ----------------
cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
# names(cells_position_ls) <- sample_df$group[match(names(cells_position_ls), sample_df$sample)]

# roi boundary
all_roi_file_ls <- list.files("20260429-tumor infiltration/tumor boundary/ROI", full.names = TRUE)  
names(all_roi_file_ls) <-  all_roi_file_ls %>% basename() %>% str_remove("\\..*$") %>% str_remove("_coordinates")

all_roi_dat_ls <-
  all_roi_file_ls %>%
  map(read.csv, skip = 2) %>%
  map(~ .x %>% group_by(Selection) %>% nest %>% deframe) %>%
  flatDoubleList()


all_tumor_area_ls <- 
  all_roi_dat_ls %>%
  imap(function(dat, name){
    name <- str_split_1(name, "---") %>% .[1]
    cells_position <- cells_position_ls[[name]]
    
    
    cell_idx <-
      sp::point.in.polygon(
        cells_position$x_centroid, # point
        cells_position$y_centroid, # point
        dat$X, # polygon
        dat$Y
      ) %>%
      as.logical  # only 0 is FALSE, all others is TRUE
    
    tumor_area_1 <- cells_position[cell_idx, ]  %>% mutate(x = x_centroid, y = y_centroid )
    tumor_area_2 <- cells_position[!cell_idx, ] %>% mutate(x = x_centroid, y = y_centroid )
    
    # return
    list(
      forward = list(
        tumor_area_1 = tumor_area_1,
        tumor_area_2 = tumor_area_2
      ),
      reverse = list(
        tumor_area_1 = tumor_area_2,
        tumor_area_2 = tumor_area_1
      )
    )
  }) %>%
  flatDoubleList()

## find tumor boundary ----------
all_tumor_boundary_ls <-
  all_tumor_area_ls %>%
  imap(function(dat, name){ # browser()
    
    roi_dat <- all_roi_dat_ls[[name]]
    
    # dat: tumor_area_1
    # query: tumor_area_2
    nn_left_right <- 
      RANN::nn2(
        dat$tumor_area_1[2:3],
        dat$tumor_area_2[2:3],
        k = 10,
        searchtype = 'radius',
        radius = 20
      )
    
    non_boundary_idx <- apply(nn_left_right$nn.idx == 0, 1, all)
    if(sum(non_boundary_idx) == nrow(dat$tumor_area_2)){
      print("radius is small, no boundary point.")
    }
    
    boundary_idx_1   <- nn_left_right$nn.idx[! non_boundary_idx, 1]
    boundary_1       <- dat$tumor_area_1[boundary_idx_1, ]
    boundary_idx_2   <- ! non_boundary_idx
    boundary_2       <- dat$tumor_area_2[boundary_idx_2, ]
    
    # comuted boudary coords
    boundary_coords_df <-
      data.frame(
        x = (boundary_1$x_centroid + boundary_2$x_centroid) / 2,
        y = (boundary_1$y_centroid + boundary_2$y_centroid) / 2
      )
    
    # return
    list(
      boundary_1 = boundary_1,
      boundary_2 = boundary_2,
      boundary_coords_df = boundary_coords_df
    )
  })

# validate tumor boundary
# 若边界未能找到，则调整radius的值
p_all_boundary_ls <- 
  all_tumor_boundary_ls %>%
  imap(function(boundary, name){ # browser()
    # name <- str_split_1(name, "---") %>% .[1]
    area <- all_tumor_area_ls[[name]]
    
    ggplot() +
      geom_point(
        # left area
        data = area$tumor_area_1 %>% dplyr::slice_sample(prop = 1), # 仅绘制10%的点，加快绘制速度
        aes(x = x_centroid, y = -y_centroid),
        color = "gray",
        size = 0.01
      ) +
      geom_point(
        # right area
        data = area$tumor_area_2 %>% dplyr::slice_sample(prop = 1), # 仅绘制10%的点，加快绘制速度
        aes(x = x_centroid, y = -y_centroid),
        color = "gray50",
        size = 0.01
      ) +
      
      geom_point(
        # comute left_boundary
        data = boundary$boundary_1,
        aes(x = x, y = -y),
        color = "#ea7070",
        size = 0.01
      ) + 
      geom_point(
        # compute right_boundary
        data = boundary$boundary_2,
        aes(x = x, y = -y),
        color = "#2694ab",
        size = 0.01
      ) +
      geom_point(
        # the final tumor boundary
        data = boundary$boundary_coords_df,
        aes(x = x, y = -y),
        color = "#600473",
        size = 0.01
      ) +
      theme_void()  +
      labs(
        caption = "tumor infiltration computed in light gray area"
      )
  })

# export
pdf(
  file.path(tumor_infiltration_outdir, "tumor boundary - plot.pdf"),
  width  = 10,
  height = 6
)
p_all_boundary_ls %>% imap(~ .x + ggtitle(.y)) %>% walk(print)
dev.off()


## compute tumor infiltration -------------------
# spline_points <-
#   data.frame(
#     x = boundary$boundary_coords_df$x,
#     y = c(0, diff(boundary$boundary_coords_df$x ))
#   ) %>%
#   mutate(yy = ifelse((y / 2) > 0, floor(y/2), ceiling(y/2))) %>%
#   mutate(yyy = yy %>% map(function(x){
#     if(x == 0) return(data.frame(yyy = 0))
#     if(x < 0) return(data.frame(yyy = seq(x, -1, by = 1) * 5))
#     data.frame(yyy = seq(1, x, by = 1) * 5)
#   })) %>%
#   dplyr::select(-y, -yy) %>%
#   unnest(yyy) %>%
#   mutate(res = x - yyy)
# 
# approx(
#   boundary$boundary_coords_df$x,
#   boundary$boundary_coords_df$y,
#   xout = spline_points$res
# ) %>% plot()


all_coords_tumor_infiltration_df_ls <-
  all_tumor_boundary_ls %>%
  imap(function(boundary, name){ # browser()
    area <- all_tumor_area_ls[[name]]
    
    # run
    radius_ls <- c(250, 500, 750) %>% set_names(., .)
    
    coords_tumor_infiltration_ls <- 
      radius_ls %>%
      map(function(radius){ # browser()
        dat_nn <-
          reduceFindNN(
            dat   = area$tumor_area_1, # %>% rbind(tumor_boundary_coords) %>% distinct(),
            query = boundary$boundary_coords_df, # tumor_boundary_coords,
            radius = radius,
            mpp    = 1,
            k      = 100,
            max_k  = nrow(area$tumor_area_1) # 2000
          )
        
        
        # all nn points is tumor infiltration points
        points_idx <- dat_nn$nn.idx %>% c() %>% .[. != 0] %>% unique
        
        # infiltration area
        coords    <- area$tumor_area_1[points_idx, ]
        coords
      })
    
    # 当前浸润带减去上一层的浸润带
    coords_tumor_infiltration_ls2 <- list()
    for(i in rev(seq_along(coords_tumor_infiltration_ls))){
      name <- names(coords_tumor_infiltration_ls)[i]
      if(i == 1){
        coords_tumor_infiltration_ls2[[name]] <- coords_tumor_infiltration_ls[[1]]
        next()
      }
      present_df <- coords_tumor_infiltration_ls[[i]]
      next_df    <- coords_tumor_infiltration_ls[[i-1]]
      coords_tumor_infiltration_ls2[[name]] <- present_df[!present_df$cell_id %in% next_df$cell_id, ]
    }
    
    # print(coords_tumor_infiltration_ls2)
    coords_tumor_infiltration_df <-
      coords_tumor_infiltration_ls2 %>%
      enframe(name = "radius") %>%
      unnest(value)
    
    # print(coords_tumor_infiltration_df)
    coords_tumor_infiltration_df
  })



# validate infiltration area ---
p_all_tumor_infiltration_ls <- 
  all_coords_tumor_infiltration_df_ls %>%
  imap(function(coords_tumor_infiltration_df, name){ # browser()
    name2 <- str_split_1(name, "---") %>% .[1]
    
    cells_position <- cells_position_ls[[name2]]
    boundary_coords_df <- all_tumor_boundary_ls[[name]]$boundary_coords_df
    
    cells_position %>%
      dplyr::slice_sample(prop = 1) %>%
      ggplot(aes(x = x_centroid , y = -y_centroid)) +
      geom_point(size = 0.01) +
      geom_point(
        data = boundary_coords_df,
        aes(x = x, y = -y),
        color = "red",
        size = 0.01
      ) + 
      geom_point(
        data = coords_tumor_infiltration_df,
        aes(x = x, y = -y, color = radius),
        # color = "blue",
        alpha = 0.5,
        size = 0.01
      ) +
      theme_void()
  })

# export
pdf(
  file.path(tumor_infiltration_outdir, "tumor infiltration - plot.pdf"),
  width  = 16,
  height =8
)
p_all_tumor_infiltration_ls %>% imap(~ .x + ggtitle(.y)) %>% walk(print)
dev.off()




# export tumor infiltration area
all_coords_tumor_infiltration_df_ls %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      file.path(tumor_infiltration_outdir, str_glue("{name} - tumor infiltration area.csv"))
    )
  })






# all_coords_tumor_infiltration_df_ls添加CellType
metadata <- all_data@meta.data %>% 
  rownames_to_column("cell_id") %>% 
  select(cell_id, RCTD)

all_coords_tumor_infiltration_df_ls_1 <- lapply(all_coords_tumor_infiltration_df_ls, function(df) {
  df <- left_join(df, metadata, by = "cell_id")
  return(df)
})

all_coords_tumor_infiltration_df_ls_1$`PC---PC---forward` %>% head()


# export tumor infiltration celltype
all_coords_tumor_infiltration_df_ls_1 %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      file.path(tumor_infiltration_outdir, str_glue("{name} - tumor infiltration celltype.csv"))
    )
  })




## cell ratio --------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

all_coords_tumor_infiltration_df <- bind_rows(all_coords_tumor_infiltration_df_ls_1, .id = "list_name")


all_coords_tumor_infiltration_df <- all_coords_tumor_infiltration_df %>% 
  mutate(
    direction = map_chr(str_split(list_name, "---"), ~ tail(.x, 1)), 
    Group = paste(direction, radius, sep = "_")
  )
head(all_coords_tumor_infiltration_df[, c("list_name", "radius", "direction", "RCTD", "Group")])


data <- all_coords_tumor_infiltration_df %>% filter(RCTD != "Unknown")
colnames(data)
data <- data[, c(3, 21, 19)]
colnames(data) <- c("Cells", "Group", "Annotation")
data_clean <- data %>%
  filter(!(Group %in% c("reverse_250", "reverse_500", "reverse_750") & Annotation == "Tumor"))


# mydata <- data %>% count(Sample, Annotation)
# mymtrix <- dcast(mydata, Sample~Annotation)
# mymtrix[is.na(mymtrix)] <- 0
# mymtrix_out <- t(mymtrix)
# write.table(as.data.frame(mymtrix_out), file = file.path(out_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")
# 
# rownames(mymtrix) <- mymtrix[, "Sample"]
# mymtrix <- mymtrix[, -1]
# mymtrix <- as.matrix(mymtrix)
# percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
#   mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
# }))
# rownames(percent) <- rownames(mymtrix)
# data4plot1 <- melt(percent)
# percent <- t(percent)
# DATA <- cbind(rownames(percent), percent)
# DATA <- rbind(c("Sample", colnames(percent)), DATA)
# write.table(DATA, file = file.path(out_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
# colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data_clean

write.table(final_data, file = file.path(tumor_infiltration_outdir, "Tumor infiltration.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# forward 
final_data_forward <- final_data %>% filter(Group %in% c("forward_250", "forward_500", "forward_750"))

mydata <- final_data_forward %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(tumor_infiltration_outdir, "Cell_type_stat_group_forward.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(tumor_infiltration_outdir, "Cell_type_stat_percent_group_forward.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
# data4plot2$Group <- factor(data4plot2$Group, levels = c("forward_250", "forward_500", "forward_750", "reverse_250", "reverse_500", "reverse_750"))

plot1 <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Group", y = "", title = "") + 
  scale_fill_manual(values = plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))


ggsave(plot1, file = file.path(tumor_infiltration_outdir, "ratio_barplot_forward.pdf"), width = 8, height = 8)
ggsave(plot1, file = file.path(tumor_infiltration_outdir, "ratio_barplot_forward.png"), width = 8, height = 8, dpi = 300)


# reverse
final_data_reverse <- final_data %>% filter(Group %in% c("reverse_250", "reverse_500", "reverse_750"))

mydata <- final_data_reverse %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(tumor_infiltration_outdir, "Cell_type_stat_group_reverse.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(tumor_infiltration_outdir, "Cell_type_stat_percent_group_reverse.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
# data4plot2$Group <- factor(data4plot2$Group, levels = c("forward_250", "forward_500", "forward_750", "reverse_250", "reverse_500", "reverse_750"))

plot2 <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Group", y = "", title = "") + 
  scale_fill_manual(values = plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))


ggsave(plot2, file = file.path(tumor_infiltration_outdir, "ratio_barplot_reverse.pdf"), width = 8, height = 8)
ggsave(plot2, file = file.path(tumor_infiltration_outdir, "ratio_barplot_reverse.png"), width = 8, height = 8, dpi = 300)








# 7.stat_density_2d -------------
plot_outdir <- file.path(results_outdir, "7.contour plot")
dir.create(plot_outdir)

color <- c('lightgrey', 'blue','seagreen2')

sp_genes <- c("PECAM1")
dat_expr <- FetchData(all_data_filt, vars = sp_genes) %>% rownames_to_column("cell_id")
cells_boundary_df <- bind_rows(cells_boundary_ls, .id = "sample")

dat <- cells_boundary_df %>% right_join(dat_expr, by = "cell_id")


pdf(file.path(plot_outdir, "cell contour plot PECAM1.pdf"), width = 20, height = 20)
dat %>% 
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(fill = .data[["PECAM1"]], group = cell_id)) +
  stat_density_2d(mapping = aes(x = vertex_x, y = vertex_y),
                  linemitre = 1000) + 
  scale_y_reverse() + 
  coord_equal() + 
  scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") + 
  theme_void() + 
  labs(
    title = "PECAM1 Expression contour plot"
  )
dev.off()


pdf(file.path(plot_outdir, "cell contour plot PECAM1.pdf"), width = 20, height = 20)
dat %>% 
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(fill = .data[["PECAM1"]], group = cell_id)) +
  stat_density_2d() + 
  scale_y_reverse() + 
  coord_equal() + 
  scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") + 
  theme_void() + 
  labs(
    title = "PECAM1 Expression contour plot"
  )
dev.off()


cells_boundary_ls %>% 
  .[3] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data_filt, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    dat <- 
      cells_boundary %>% 
      right_join(dat_expr, by = "cell_id")
    
    p_ls <- 
      sp_genes %>% 
      set_names(., .) %>% 
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>% 
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          # geom_polygon(aes(fill = .data[[gene]], group = cell_id)) + 
          # stat_density_2d(aes(alpha = after_stat(level)), color = "black", bins = 10) +
          stat_density_2d(aes(fill = .data[[gene]]), geom = "polygon") +
          # stat_density_2d(aes(alpha = after_stat(level)), geom = "polygon", bins = 10, fill = "mistyrose", color = NA) +
          # stat_density_2d(
          #   aes(weight = .data[[gene]], alpha = after_stat(level)),
          #   geom = "polygon",
          #   bins = 10,
          #   fill = "mistyrose",
          #   color = NA) +
          # stat_density_2d(
          #   # aes(weight = .data[[gene]], color = after_stat(level)),
          #   aes(weight = .data[[gene]], color = "black"),
          #   bins = 10,
          #   size = 0.2,
          #   alpha = 0.5) +
          scale_y_reverse() + 
          coord_equal() + 
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") + 
          theme_void() + 
          labs(
            title = project
          )
        p
      })
    p_ls
  })





# 8.ROI cell ratio ---------
roi_outdir <- file.path(results_outdir, "8.ROI cell ratio")
dir.create(roi_outdir)

cols <- pal_igv("default")(51)
celltype <- sort(unique(all_data$RCTD))
plan2cell <- cols[1:length(unique(all_data$RCTD))]
plan2cell <- setNames(plan2cell, celltype)


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  filter(RCTD != "Unknown") %>% 
  filter(ROI != "Other")
all_data_meta$ROI <- factor(all_data_meta$ROI, levels = c("Normal_cell", "adenoma", "Non_classic", "Tumor"))

p_col_proj <- 
  all_data_meta %>% 
  group_by(ROI, RCTD) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  ggplot(aes(x = ROI, y = freq, fill = RCTD)) + 
  geom_col() + 
  # ggsci::scale_fill_igv() +
  scale_fill_manual(values = plan2cell) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj, 
  filename = file.path(roi_outdir, "cell_freq_project.pdf"), 
  height = 6, 
  width  = 6
)






# 9.B cell subcluster ---------
result_dir <- file.path(results_outdir, "9.B cell subcluster")
dir.create(result_dir)

B_cell <- all_data_filt[, all_data_filt@meta.data$RCTD == "B cell"]

### PCA --------
DefaultAssay(B_cell) <- "RNA"
B_cell <- FindVariableFeatures(B_cell)
B_cell <- ScaleData(B_cell, vars.to.regress = 'nCount_RNA')
B_cell <- RunPCA(B_cell)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B_cell, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(B_cell, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(B_cell, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(B_cell, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(B_cell, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 16

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


B_cell <- RunUMAP(B_cell, dims = 1:dims,return.model = TRUE)



pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(B_cell, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



B_cell <- B_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

B_cell <- RunUMAP(
  B_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(B_cell, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
B_cell  <-
  FindNeighbors(
    B_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


B_cell <-
  FindClusters(
    B_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




### 去除1,4,6,7,8 ----------
B_cell_filt <- B_cell[, B_cell$RNA_snn_res.0.4 %in% c(0, 2, 3, 5, 9)]

## annotation -----------
cell_cluster <- c(
  "0" = "C1",
  "2" = "C2",
  "3" = "C3",
  "5" = "C4",
  "9" = "C5"
)

Idents(B_cell_filt) <- "RNA_snn_res.0.4"

B_cell_filt <- RenameIdents(B_cell_filt, cell_cluster)
B_cell_filt$subtype <- Idents(B_cell_filt)





## find markers -------
Idents(B_cell_filt) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(B_cell_filt), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    B_cell_filt, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 10
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(B_cell_filt, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.8), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(B_cell_filt, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- B_cell_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 57)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()

for (x in c(1:4)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[5]
plot_data_celltype <- data[which(data$Celltype == celltype[5]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[5]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 1)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 16, height = 8)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 16, height = 8, dpi = 300)











# 10.Monocyte subcluster ---------
result_dir <- file.path(results_outdir, "10.Monocyte subcluster")
dir.create(result_dir)

Monocyte <- all_data_filt[, all_data_filt@meta.data$RCTD == "Monocyte"]

### PCA --------
DefaultAssay(Monocyte) <- "RNA"
Monocyte <- FindVariableFeatures(Monocyte)
Monocyte <- ScaleData(Monocyte, vars.to.regress = 'nCount_RNA')
Monocyte <- RunPCA(Monocyte)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Monocyte, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Monocyte, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Monocyte, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Monocyte, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Monocyte, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 13

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Monocyte <- RunUMAP(Monocyte, dims = 1:dims,return.model = TRUE)



pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Monocyte, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Monocyte, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Monocyte <- Monocyte %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Monocyte <- RunUMAP(
  Monocyte, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Monocyte, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Monocyte, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Monocyte  <-
  FindNeighbors(
    Monocyte ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Monocyte <-
  FindClusters(
    Monocyte,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Monocyte, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Monocyte, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Monocyte, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Monocyte, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot
theme_45_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 10, height = 5)
DotPlot(Monocyte, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_45_angle
DotPlot(Monocyte, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme_45_angle
DotPlot(Monocyte, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme_45_angle
DotPlot(Monocyte, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme_45_angle
dev.off()


Mon_Markers <- c(
  "cMo" = "FCN1",
  "cMo" = "CD14",
  "cMo" = "CD36",
  "cMo" = "SELL",
  "cMo" = "S100A8",
  "cMo" = "S100A12",
  "cMo" = "CLEC12A",
  "cMo" = "MS4A6A",
  "cMo" = "CXCL14",
  # "ncMo" = "FCN1",
  "ncMo" = "FCGR3A",
  "ncMo" = "FCGR3B",
  "ncMo" = "LILRA1",
  "ncMo" = "LILRB1",
  "ncMo" = "CX3CR1",
  "ncMo" = "IFITM1",
  "ncMo" = "ICAM2",
  "ncMo" = "MTSS1",
  "ncMo" = "CDKN1C",
  "ncMo" = "CDH23",
  "ncMo" = "SLC44A2",
  # "intMo" = "FCN1",
  # "intMo" = "CD14",
  # "intMo" = "FCGR3A",
  "intMo" = "CD300E",
  "intMo" = "FAM65B"
)
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} Monocyte_marker.pdf")), width = 10, height = 5)
DotPlot(Monocyte, group.by = "RNA_snn_res.0.2", features = Mon_Markers) + theme_45_angle
dev.off()

Mon_Markers_list <- split(Mon_Markers, names(Mon_Markers))

all_genes <- rownames(Monocyte)
filtered_pathways <- lapply(Mon_Markers_list, function(genes) {
  intersect(genes, all_genes)
})
filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > 0]
pathway_names <- names(filtered_pathways)

Monocyte <- AddModuleScore(
  object = Monocyte,
  features = filtered_pathways,
  name = "Pathway_",
  ctrl = 50,
  seed = 42
)
colnames(Monocyte@meta.data)[grep("Pathway_", colnames(Monocyte@meta.data))] <- pathway_names

marker_meta <- Monocyte@meta.data %>%
  mutate(cluster = paste0("C", RNA_snn_res.0.2)) %>%
  dplyr::select(cluster, all_of(pathway_names))

pathway_means <- marker_meta %>%
  group_by(cluster) %>%
  summarise(across(all_of(pathway_names), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)

H <- pheatmap(
  mat = t(as.matrix(pathway_means)),
  color = my_palette,
  border_color = NA,
  cluster_rows = FALSE,  # 是否对行聚类
  cluster_cols = TRUE,  # 是否对列聚类
  clustering_method = "ward.D2",  # 聚类方法
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize_row = 8,      # 行名字体大小
  fontsize_col = 10,     # 列名字体大小
  angle_col = 45,        # 列名角度
  main = "Pathway Activity",  # 标题
  legend = TRUE,         # 显示图例
  scale = "column",       # 列标准化（使不同通路可比）
  #height = 15,
  #width = 5,
  cellwidth = 18,         # 单元格宽度
  cellheight = 12         # 单元格高度
)

ggsave(H, filename = file.path(result_dir, "Monocyte_subtype_heatmap.pdf"), height = 10, width = 10,)





### 去除2,  ----------
Monocyte_filt <- Monocyte[, Monocyte$RNA_snn_res.0.2 %in% c(0, 1, 3, 4, 5, 6)]

## annotation -----------
cell_cluster <- c(
  "0" = "cMo",
  "1" = "cMo",
  "3" = "intMo",
  "4" = "cMo",
  "5" = "ncMo",
  "6" = "ncMo"
)

Idents(Monocyte_filt) <- "RNA_snn_res.0.2"

Monocyte_filt <- RenameIdents(Monocyte_filt, cell_cluster)
Monocyte_filt$subtype <- Idents(Monocyte_filt)





## find markers -------
Idents(Monocyte_filt) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(Monocyte_filt), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    Monocyte_filt, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 10
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(Monocyte_filt, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.1, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(Monocyte_filt, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- Monocyte_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 57)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()

for (x in c(1:2)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[3]
plot_data_celltype <- data[which(data$Celltype == celltype[3]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[3]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 1)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 12, height = 8)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 12, height = 8, dpi = 300)











# 11.Macrophage subcluster ---------
result_dir <- file.path(results_outdir, "11.Macrophage subcluster")
dir.create(result_dir)

Macrophage <- all_data_filt[, all_data_filt@meta.data$RCTD == "Macrophage"]

### PCA --------
DefaultAssay(Macrophage) <- "RNA"
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 12

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)



pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Macrophage, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Macrophage, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot
theme_45_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 10, height = 5)
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_45_angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme_45_angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme_45_angle
DotPlot(Macrophage, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme_45_angle
dev.off()


major_marker <- c(
  # 'SPP1',
  M1 = 'NOS2',
  M1 = 'IL12',
  M1 = 'FCGR1A',
  M1 = 'FCGR1B',
  M1 = 'FCGR1C',
  M1 = 'CD80',
  M1 = 'CXCR10',
  M1 = 'IL23',
  M1 = 'CXCL9',
  M1 = 'CXCL10',
  M1 = 'CXCL11',
  M1 = 'CD86',
  M1 = 'IL1A',
  M1 = 'IL1B',
  M1 =  'IL6',
  M1 =  'TNF',
  M1 = 'MHCII',
  M1 = 'CCL5',
  M1 = 'IRF5',
  M1 = 'IRF1',
  M1 = 'CD40',
  M1 = 'KYNU',
  M1 = 'CCR7',
  M1 =  'VCAN',
  M1 =  'CD300E',
  M1 =  'INHBA',
  M1 =  'TNIP3',
  M2 = 'ARG1',
  M2 = 'IL10',
  M2 = 'CD32',
  M2 = 'CD163',
  M2 = 'FCER2',
  M2 = 'CD200R1',
  M2 = 'PDCD1LG2',
  M2 = 'CD274',
  M2 = 'MARCO',
  M2 = 'CSF1R',
  M2 = 'MRC1',
  M2 = 'IL1RN',
  M2 = 'Il1R2',
  M2 = 'IL4R',
  M2 = 'CCL4',
  M2 = 'CCL13',
  M2 = 'CCL20',
  M2 = 'CCL17',
  M2 = 'CCL18',
  M2 = 'CCL22',
  M2 = 'CCL24',
  M2 = 'LYVE1',
  M2 = 'VEGFA',
  M2 = 'VEGFB',
  M2 = 'VEGFC',
  M2 = 'VEGFD',
  M2 = 'EGF',
  M2 = 'CTSA',
  M2 = 'CTSB',
  M2 = 'CSTC',
  M2 = 'CTSD',
  M2 = 'TGFB1',
  M2 = 'CD209',
  M2 = 'CD163L1',
  M2 = 'FOLR2'
)
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} Macrophage_marker.pdf")), width = 10, height = 5)
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = major_marker) + theme_45_angle
dev.off()

pdf(file.path(result_dir, str_glue("vlnplot dims{dims} {reduction.use} Macrophage_marker.pdf")), width = 10, height = 6)
VlnPlot(Macrophage, features = c("M1", "M2"), pt.size = 0)
dev.off()



Mon_Markers_list <- split(major_marker, names(major_marker))

all_genes <- rownames(Macrophage)
filtered_pathways <- lapply(Mon_Markers_list, function(genes) {
  intersect(genes, all_genes)
})
filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > 0]
pathway_names <- names(filtered_pathways)

Macrophage <- AddModuleScore(
  object = Macrophage,
  features = filtered_pathways,
  name = "Pathway_",
  ctrl = 50,
  seed = 42
)
colnames(Macrophage@meta.data)[grep("Pathway_", colnames(Macrophage@meta.data))] <- pathway_names

marker_meta <- Macrophage@meta.data %>%
  mutate(cluster = paste0("C", RNA_snn_res.0.4)) %>%
  dplyr::select(cluster, all_of(pathway_names))

pathway_means <- marker_meta %>%
  group_by(cluster) %>%
  summarise(across(all_of(pathway_names), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)

H <- pheatmap(
  mat = t(as.matrix(pathway_means)),
  color = my_palette,
  border_color = NA,
  cluster_rows = FALSE,  # 是否对行聚类
  cluster_cols = TRUE,  # 是否对列聚类
  clustering_method = "ward.D2",  # 聚类方法
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize_row = 8,      # 行名字体大小
  fontsize_col = 10,     # 列名字体大小
  angle_col = 45,        # 列名角度
  main = "Pathway Activity",  # 标题
  legend = TRUE,         # 显示图例
  scale = "column",       # 列标准化（使不同通路可比）
  #height = 15,
  #width = 5,
  cellwidth = 18,         # 单元格宽度
  cellheight = 12         # 单元格高度
)

ggsave(H, filename = file.path(result_dir, "Macrophage_subtype_heatmap.pdf"), height = 10, width = 10,)





### 去除2,  ----------
# Macrophage_filt <- Macrophage[, Macrophage$RNA_snn_res.0.2 %in% c(0, 1, 3, 4, 5, 6)]

## annotation -----------
cell_cluster <- c(
  "0" = "Epithelial_like_Macrophage",
  "1" = "M1",
  "2" = "M2",
  "3" = "M1",
  "4" = "Cycling_M2",
  "5" = "Fibroblast_like_Macrophage",
  "6" = "Fibroblast_like_Macrophage",
  "7" = "M1"
)


Idents(Macrophage) <- "RNA_snn_res.0.4"
Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)




## find markers -------
Idents(Macrophage) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]
clean_top5 <- unique(all_top5_gene)

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    Macrophage, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 10
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(Macrophage, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 8), #图例字体大小和间距
    legend.key.height = unit(0.4, "cm"),
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.80, 0.88), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(Macrophage, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- Macrophage[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 59)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()

for (x in c(1:4)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 15, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[5]
plot_data_celltype <- data[which(data$Celltype == celltype[5]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[5]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 15, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 1)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 15, height = 8)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 15, height = 8, dpi = 300)













# 12.Mast subcluster ---------
result_dir <- file.path(results_outdir, "12.Mast cell subcluster")
dir.create(result_dir)

Mast <- all_data_filt[, all_data_filt@meta.data$RCTD == "Mast cell"]

### PCA --------
DefaultAssay(Mast) <- "RNA"
Mast <- FindVariableFeatures(Mast)
Mast <- ScaleData(Mast, vars.to.regress = 'nCount_RNA')
Mast <- RunPCA(Mast)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Mast, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Mast, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Mast, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Mast, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Mast, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 7

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Mast <- RunUMAP(Mast, dims = 1:dims,return.model = TRUE)



pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Mast, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Mast, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Mast <- Mast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Mast <- RunUMAP(
  Mast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Mast, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Mast, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Mast  <-
  FindNeighbors(
    Mast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Mast <-
  FindClusters(
    Mast,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Mast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Mast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot
theme_45_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 10, height = 5)
DotPlot(Mast, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_45_angle
DotPlot(Mast, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme_45_angle
DotPlot(Mast, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme_45_angle
DotPlot(Mast, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme_45_angle
dev.off()


# Mon_Markers <- c(
#   "cMo" = "FCN1",
#   "cMo" = "CD14",
#   "cMo" = "CD36",
#   "cMo" = "SELL",
#   "cMo" = "S100A8",
#   "cMo" = "S100A12",
#   "cMo" = "CLEC12A",
#   "cMo" = "MS4A6A",
#   "cMo" = "CXCL14",
#   # "ncMo" = "FCN1",
#   "ncMo" = "FCGR3A",
#   "ncMo" = "FCGR3B",
#   "ncMo" = "LILRA1",
#   "ncMo" = "LILRB1",
#   "ncMo" = "CX3CR1",
#   "ncMo" = "IFITM1",
#   "ncMo" = "ICAM2",
#   "ncMo" = "MTSS1",
#   "ncMo" = "CDKN1C",
#   "ncMo" = "CDH23",
#   "ncMo" = "SLC44A2",
#   # "intMo" = "FCN1",
#   # "intMo" = "CD14",
#   # "intMo" = "FCGR3A",
#   "intMo" = "CD300E",
#   "intMo" = "FAM65B"
# )
# pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} Mast_marker.pdf")), width = 10, height = 5)
# DotPlot(Mast, group.by = "RNA_snn_res.0.2", features = Mon_Markers) + theme_45_angle
# dev.off()
# 
# Mon_Markers_list <- split(Mon_Markers, names(Mon_Markers))
# 
# all_genes <- rownames(Mast)
# filtered_pathways <- lapply(Mon_Markers_list, function(genes) {
#   intersect(genes, all_genes)
# })
# filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > 0]
# pathway_names <- names(filtered_pathways)
# 
# Mast <- AddModuleScore(
#   object = Mast,
#   features = filtered_pathways,
#   name = "Pathway_",
#   ctrl = 50,
#   seed = 42
# )
# colnames(Mast@meta.data)[grep("Pathway_", colnames(Mast@meta.data))] <- pathway_names
# 
# marker_meta <- Mast@meta.data %>%
#   mutate(cluster = paste0("C", RNA_snn_res.0.2)) %>%
#   dplyr::select(cluster, all_of(pathway_names))
# 
# pathway_means <- marker_meta %>%
#   group_by(cluster) %>%
#   summarise(across(all_of(pathway_names), mean, na.rm = TRUE)) %>%
#   tibble::column_to_rownames("cluster") %>%
#   as.matrix()
# 
# my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
# 
# H <- pheatmap(
#   mat = t(as.matrix(pathway_means)),
#   color = my_palette,
#   border_color = NA,
#   cluster_rows = FALSE,  # 是否对行聚类
#   cluster_cols = TRUE,  # 是否对列聚类
#   clustering_method = "ward.D2",  # 聚类方法
#   show_rownames = TRUE,  # 显示行名
#   show_colnames = TRUE,  # 显示列名
#   fontsize_row = 8,      # 行名字体大小
#   fontsize_col = 10,     # 列名字体大小
#   angle_col = 45,        # 列名角度
#   main = "Pathway Activity",  # 标题
#   legend = TRUE,         # 显示图例
#   scale = "column",       # 列标准化（使不同通路可比）
#   #height = 15,
#   #width = 5,
#   cellwidth = 18,         # 单元格宽度
#   cellheight = 12         # 单元格高度
# )
# 
# ggsave(H, filename = file.path(result_dir, "Mast_subtype_heatmap.pdf"), height = 10, width = 10,)
# 




### 去除1,  ----------
Mast_filt <- Mast[, Mast$RNA_snn_res.0.2 %in% c(0, 2, 3)]

## annotation -----------
cell_cluster <- c(
  "0" = "C1",
  "2" = "C2",
  "3" = "C3"
)

Idents(Mast_filt) <- "RNA_snn_res.0.2"

Mast_filt <- RenameIdents(Mast_filt, cell_cluster)
Mast_filt$subtype <- Idents(Mast_filt)





## find markers -------
Idents(Mast_filt) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(Mast_filt), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    Mast_filt, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 10
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(Mast_filt, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.1, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(Mast_filt, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- Mast_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 57)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()

for (x in c(1:2)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[3]
plot_data_celltype <- data[which(data$Celltype == celltype[3]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[3]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 1)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 10, height = 8)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 10, height = 8, dpi = 300)










# 13.Neutrophil subcluster ---------
result_dir <- file.path(results_outdir, "13.Neutrophil subcluster")
dir.create(result_dir)

Neutrophil <- all_data_filt[, all_data_filt@meta.data$RCTD == "Neutrophil"]

### PCA --------
DefaultAssay(Neutrophil) <- "RNA"
Neutrophil <- FindVariableFeatures(Neutrophil)
Neutrophil <- ScaleData(Neutrophil, vars.to.regress = 'nCount_RNA')
Neutrophil <- RunPCA(Neutrophil)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Neutrophil, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Neutrophil, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Neutrophil, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Neutrophil, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Neutrophil, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 8

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Neutrophil <- RunUMAP(Neutrophil, dims = 1:dims,return.model = TRUE)

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Neutrophil, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Neutrophil, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Neutrophil <- Neutrophil %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Neutrophil <- RunUMAP(
  Neutrophil, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Neutrophil, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Neutrophil, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Neutrophil  <-
  FindNeighbors(
    Neutrophil ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Neutrophil <-
  FindClusters(
    Neutrophil,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Neutrophil, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neutrophil, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neutrophil, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Neutrophil, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot
theme_45_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 10, height = 5)
DotPlot(Neutrophil, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_45_angle
DotPlot(Neutrophil, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme_45_angle
DotPlot(Neutrophil, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme_45_angle
DotPlot(Neutrophil, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme_45_angle
dev.off()






### 去除0  ----------
Neutrophil_filt <- Neutrophil[, Neutrophil$RNA_snn_res.0.4 %in% c(1, 2, 3)]

## annotation -----------
cell_cluster <- c(
  "1" = "C1",
  "2" = "C2",
  "3" = "C3"
)

Idents(Neutrophil_filt) <- "RNA_snn_res.0.4"

Neutrophil_filt <- RenameIdents(Neutrophil_filt, cell_cluster)
Neutrophil_filt$subtype <- Idents(Neutrophil_filt)





## find markers -------
Idents(Neutrophil_filt) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(Neutrophil_filt), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    Neutrophil_filt, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 10
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(Neutrophil_filt, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.1, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(Neutrophil_filt, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- Neutrophil_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 57)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()

for (x in c(1:2)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[3]
plot_data_celltype <- data[which(data$Celltype == celltype[3]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[3]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 23, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 1)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 10, height = 8)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 10, height = 8, dpi = 300)











# 14.Tcell subcluster ---------
result_dir <- file.path(results_outdir, "14.Tcell subcluster")
dir.create(result_dir)

Tcell <- all_data_filt[, all_data_filt@meta.data$RCTD == "T cell"]

### PCA --------
DefaultAssay(Tcell) <- "RNA"
Tcell <- FindVariableFeatures(Tcell)
Tcell <- ScaleData(Tcell, vars.to.regress = 'nCount_RNA')
Tcell <- RunPCA(Tcell)


### pca inspector ------
pdf(file.path(result_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Tcell, ndims = 30)
dev.off()

pdf(file.path(result_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Tcell, dims = 1:5)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Tcell, dims = 6:10)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Tcell, dims = 11:15)
dev.off()
pdf(file.path(result_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Tcell, dims = 16:20)
dev.off()


### dimension reduction --------
dims <- 12

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Tcell <- RunUMAP(Tcell, dims = 1:dims,return.model = TRUE)

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Tcell, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - split.pdf")), height = 5, width = 12)
DimPlot(Tcell, group.by = "sample", reduction = "umap", split.by = "sample")
dev.off()



Tcell <- Tcell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Tcell <- RunUMAP(
  Tcell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Tcell, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} sample - harmony - split.pdf")), height = 5, width = 12)
DimPlot(Tcell, group.by = "sample", reduction = "harmonyumap", split.by = "sample")
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Tcell  <-
  FindNeighbors(
    Tcell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Tcell <-
  FindClusters(
    Tcell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(result_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Tcell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Tcell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot
theme_45_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf(file.path(result_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 10, height = 5)
DotPlot(Tcell, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme_45_angle
DotPlot(Tcell, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme_45_angle
DotPlot(Tcell, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme_45_angle
DotPlot(Tcell, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme_45_angle
dev.off()


### dotplot -------
major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  # 'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCL5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  # 'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

major_marker_ls_CD4 <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'T' ='CD3E',
  'CD4'='CD4',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2'
)


major_marker_ls_CD8 <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAYN',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Trm'='ZNF683',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)
# dotplot
pdf(file.path(result_dir, str_glue("dotplot dims{dims} major_marker_CD4_CD8.pdf")), width = 15, height = 6)
# DotPlot(All_T, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Tcell, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Tcell, group.by = "RNA_snn_res.0.4", features = major_marker_ls_CD4)+RotatedAxis()
DotPlot(Tcell, group.by = "RNA_snn_res.0.4", features = major_marker_ls_CD8)+RotatedAxis()
dev.off()



### 去除0  ----------
# Tcell_filt <- Tcell[, Tcell$RNA_snn_res.0.4 %in% c(1, 2, 3)]

## annotation -----------
cell_cluster <- c(
  "0" = "CD8_Tem",
  "1" = "Fibroblast_like_T",
  "2" = "CD4_Tcm",
  "3" = "Myeloid_like_T",
  "4" = "Epithelial_like_T",
  "5" = "CD4_naive",
  "6" = "Endothelial_like_T",
  "7" = "NK",
  "8" = "CD4_Treg",
  "9" = "CD8_Proliferating"
)

Idents(Tcell) <- "RNA_snn_res.0.4"

Tcell <- RenameIdents(Tcell, cell_cluster)
Tcell$subtype <- Idents(Tcell)





## find markers -------
Idents(Tcell) <- "subtype"
all_markers <-
  FindAllMarkers(
    JoinLayers(Tcell), 
    test.use        = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct         = 0.1, 
    only.pos        = TRUE
  )

# export marker gene
all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(result_dir, "marker.csv"), 
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    Tcell, 
    features = clean_top5,
    dot.scale = 10,
    col.min = -1,
    dot.min = 0.025
    # cols = c("blue", "white", "#e40000")
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("blue", "white", "#e40000"))

ggsave(
  p_dot_top5, 
  filename = file.path(result_dir, "top 5 gene dotplot.pdf"), 
  height = 5, 
  width  = 18
)






colour_py <- c(
  "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", 
  "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", 
  "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", 
  "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

pdf(file.path(result_dir, str_glue("Dimplot cluster.pdf")), height = 5.5, width = 6)
DimPlot(Tcell, group.by = "subtype", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 8), #图例字体大小和间距
    legend.key.height = unit(0.4, "cm"),
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.138, 0.80), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()


pdf(file.path(result_dir, str_glue("Dimplot group.pdf")), height = 5.5, width = 6)
DimPlot(Tcell, group.by = "group", label = F, reduction = 'harmonyumap', cols = colour_py) +
  theme_bw() +
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    legend.position = c(0.08, 0.85), #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Group") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
dev.off()




data <- Tcell[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 57)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(result_dir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(result_dir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")


library(ggpubr)
data <- read.table(file.path(result_dir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", check.names = F)
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", check.names = F)

data <- melt(data)
colnames(data) <- c("Celltype", "SampleID", "percent")
data <- left_join(data, sample_map)

orders <- read.table(file.path(result_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data$COND1 <- factor(data$COND1, levels = unique(sample_map[, 2]))
data$SampleID <- factor(data$SampleID, levels = unique(sample_map[, 1]))
data$Celltype <- factor(data$Celltype, levels = orders[, 1])
data <- data[which(data$Celltype != "Unclassified"), ]

# 设置比较组方案
my_comparisons <- list( c("PC", "PA") )

# 提取细胞类型
celltype <- sort(unique(data$Celltype))

# 绘图
plot_list <- list()
for (x in c(1:4, 6:9)){
  celltype_select <- celltype[x]
  plot_data_celltype <- data[which(data$Celltype == celltype[x]), ]
  max_num <- max(plot_data_celltype$percent)
  plot_list[[x]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
    geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
    geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
    stat_compare_means(size = 7, 
                       color = "red", 
                       label.y = c(max_num, max_num*1.05, max_num*1.15),  # 有几个比较组就要设置几个位置
                       comparisons = my_comparisons, 
                       label = "p.signif") + 
    ylim(c(0, max_num*1.3)) + 
    theme_test() + 
    ylab("percent(%)") + 
    theme(plot.title = element_text(size = 18, face = 'bold', vjust = 0.5, hjust = 0.5), 
          legend.text = element_text(colour = "black", size = 15), 
          legend.title = element_text(face = "bold", colour = "black", size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20, face = 'bold'), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 15, color = "black"), 
          legend.position = "right") + 
    labs(x = "", fill = "Group", title = celltype_select) + NoLegend() + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细
}

celltype_select <- celltype[5]
plot_data_celltype <- data[which(data$Celltype == celltype[5]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[5]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 18, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

celltype_select <- celltype[10]
plot_data_celltype <- data[which(data$Celltype == celltype[10]),]
max_num <- max(plot_data_celltype$percent)
plot_list[[10]] <- ggplot(plot_data_celltype, aes(COND1, percent, fill = COND1)) + 
  geom_boxplot(width = 0.7, cex = 1, notchwidth = 0.7, outlier.colour = NA) + 
  geom_point(aes(color = SampleID), position = position_dodge(width = 0), size = 5) + 
  stat_compare_means(size = 7, 
                     color = "red", 
                     label.y = c(max_num, max_num*1.05, max_num*1.15), 
                     comparisons = my_comparisons, 
                     label = "p.signif") + 
  ylim(c(0, max_num*1.3)) + 
  theme_test() + 
  ylab("percent(%)") + 
  theme(plot.title = element_text(size = 18, face = 'bold', vjust = 0.5, hjust = 0.5), 
        legend.box = "horizontal", 
        legend.text = element_text(colour = "black", size = 15), 
        legend.title = element_text(face = "bold", colour = "black", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"), 
        legend.position = "right") + 
  labs(x = "", fill = "Group", title = celltype_select) + 
  guides(color = guide_legend(override.aes = list(size = 5)), fill = FALSE) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) #边框粗细

plot_final <- wrap_plots(plot_list, nrow = 2)
ggsave(filename = file.path(result_dir, "celltype_boxplot.pdf"), plot_final, width = 15, height = 15)
ggsave(filename = file.path(result_dir, "celltype_boxplot.png"), plot_final, width = 15, height = 15, dpi = 300)













# 15.diff analysis PA vs PC ------------
diff_outdir <- file.path(results_outdir, "15.diff analysis PA vs PC")
dir.create(diff_outdir)

library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Hs.eg.db")
# library("org.Mm.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)

diffmarkers <-
  FindMarkers(
    JoinLayers(all_data_filt), 
    ident.1 = 'PA', # 处理组
    ident.2 = 'PC', # 比较组
    test.use = "wilcox", 
    logfc.threshold = 0.25, 
    min.pct = 0.1, 
    group.by = 'group' # 分组列
  )

if (nrow(diffmarkers) == 0) {
  cat("No DEGs found. Skipping enrichment analysis.\n")
}


diffmarkers_up <- diffmarkers %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
diffmarkers_down <- diffmarkers %>% filter(avg_log2FC < 0) %>% arrange(desc(avg_log2FC))


write.csv(diffmarkers, file.path(diff_outdir, 'group diffmarker.csv'))
write.csv(diffmarkers_up, file.path(diff_outdir, 'group diffmarker up.csv'))
write.csv(diffmarkers_down, file.path(diff_outdir, 'group diffmarker down.csv'))



enrichGO_up_ls <- enrichGO(
  rownames(diffmarkers_up), 
  OrgDb = org.Hs.eg.db, 
  ont='ALL', 
  pAdjustMethod = 'BH', 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = 'SYMBOL')
write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, 'Go up enriment table.csv'), row.names = FALSE)

enrichGO_up_dot <- enrichplot::dotplot(
  enrichGO_up_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichGO_up_dot, 
  filename = file.path(diff_outdir, 'Go up dotplot.pdf'), 
  height = 6, 
  width = 8)
enrichGO_up_bar <- barplot(
  enrichGO_up_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichGO_up_bar, 
  filename = file.path(diff_outdir, 'Go up barplot.pdf'), 
  height = 6, 
  width = 8)


# GO down
enrichGO_down_ls <- enrichGO(
  rownames(diffmarkers_down), 
  OrgDb = org.Hs.eg.db, 
  ont='ALL', 
  pAdjustMethod = 'BH', 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = 'SYMBOL')
write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, 'Go down enriment table.csv'), row.names = FALSE)

enrichGO_down_dot <- enrichplot::dotplot(
  enrichGO_down_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichGO_down_dot, 
  filename = file.path(diff_outdir, 'Go down dotplot.pdf'), 
  height = 6, 
  width = 8)
enrichGO_down_bar <- barplot(
  enrichGO_down_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichGO_down_bar, 
  filename = file.path(diff_outdir, 'Go down barplot.pdf'), 
  height = 6, 
  width = 8)

# KEGG up
de_up_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(diffmarkers_up), 
    OrgDb = "org.Hs.eg.db", 
    fromType = "SYMBOL", 
    toType = "ENTREZID"
  ) %>% 
  .$ENTREZID
enrichKEGG_up_ls <- enrichKEGG(
  de_up_genes_id_for_enrichmenet, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10, 
  maxGSSize = 500, 
  qvalueCutoff = 0.2, 
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, 'KEGG up enriment table.csv'), row.names = FALSE)

enrichKEGG_up_dot <- enrichplot::dotplot(
  enrichKEGG_up_ls, 
  showCategory = 10, 
  label_format = 50)
enrichKEGG_up_bar <- barplot(
  enrichKEGG_up_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichKEGG_up_dot, 
  filename = file.path(diff_outdir, 'KEGG up dotplot.pdf'), 
  height = 6, 
  width = 8)
ggsave(
  enrichKEGG_up_bar, 
  filename = file.path(diff_outdir, 'KEGG up barplot.pdf'), 
  height = 6, 
  width = 8)


# KEGG down
de_down_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(diffmarkers_down), 
    OrgDb = "org.Hs.eg.db", 
    fromType = "SYMBOL", 
    toType = "ENTREZID"
  ) %>% 
  .$ENTREZID
enrichKEGG_down_ls <- enrichKEGG(
  de_down_genes_id_for_enrichmenet, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10, 
  maxGSSize = 500, 
  qvalueCutoff = 0.2, 
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, 'KEGG down enriment table.csv'), row.names = FALSE)

enrichKEGG_down_dot <- enrichplot::dotplot(
  enrichKEGG_down_ls, 
  showCategory = 10, 
  label_format = 50)
enrichKEGG_down_bar <- barplot(
  enrichKEGG_down_ls, 
  showCategory = 10, 
  label_format = 50)
ggsave(
  enrichKEGG_down_dot, 
  filename = file.path(diff_outdir, 'KEGG down dotplot.pdf'), 
  height = 6, 
  width = 8)
ggsave(
  enrichKEGG_down_bar, 
  filename = file.path(diff_outdir, 'KEGG down barplot.pdf'), 
  height = 6, 
  width = 8)










