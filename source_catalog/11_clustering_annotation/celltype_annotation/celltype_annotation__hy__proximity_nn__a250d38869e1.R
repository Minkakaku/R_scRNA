library(tidyverse)
source("spatial analysis.r")


# setting variable -----------

cell_colors <- c(
  'Acinar cells'      =  "#1C959A",
  'Ductal cells'      =  "#2C76AB",
  'Endocrine cells'   =  "#1B518D",
  # 'Epithelial cells'  =  "#",
  'Malignant/epithelial cells'   =  "#5B62AA",
  
  'B cells'           =  "#DCCF7B",
  'CD4Ts'             =  "#70621F",
  'CD8Ts'             =  "#5D7748",
  'Tregs'             =  "#52B97A",
  'NKT cells'         =  "#6A6E54",
  
  'Granulocytes'      =  "#EB5582",
  'Myeloid cells'     =  "#DA6258",
  
  'Endothelial cells' =  "#68216E",
  'Fibroblasts'       =  "#BA3B94"
)

major_cell_colors <- c(
  'Epithelial cells'      =  "#1C959A",
  
  'B cells'           =  "#DCCF7B",
  'CD4Ts'             =  "#70621F",
  'CD8Ts'             =  "#5D7748",
  'Tregs'             =  "#52B97A",
  'T/NK cells'         =  "#6A6E54",
  
  'Granulocytes'      =  "#EB5582",
  'Myeloid cells'     =  "#DA6258",
  
  'Endothelial cells' =  "#68216E",
  'Fibroblasts'       =  "#BA3B94"
)

grp_colors <- c(
  P   = "#71C3A3",
  WT  = "#87CEEB",
  MUT = "#F8766C"
)

epi_cell_colors <- c(
  'Acinar cells'      =  "#1C959A",
  'Ductal cells'      =  "#2C76AB",
  'Endocrine cells'   =  "#1B518D",
  'Malignant cells'   =  "#5B62AA"
)



# read data -----------------
all_data <- readRDS(file = "all_data_fixed.rds")

out_dir <- "figure_v2"
sp_out_dir <- "figure_v2/spatial"
dir.create(sp_out_dir)


cells_position_ls <- readRDS("materials_for_analysis/cells_position_ls.rds")
cells_position <- bind_rows(cells_position_ls)


## total data clean
# ##annotation_final levels里有字符NA，会影响CN分析，此处去掉NA
# all_data_meta <-
#   all_data_meta %>%
#   mutate(annotation_final = final_annotation %>%
#            as.character() %>%
#            factor(., levels = unique(.))
#   )
# all_data_meta <- all_data[[]]


all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      all_data[[]] %>%
      dplyr::select(final_annotation, Group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = final_annotation,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")



# proximity ---------------
res_pa_ls <- list()
res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

# saveRDS(res_pa_ls,'res_pa_ls.rds')
# res_pa_ls <- readRDS('res_pa_ls.rds')

# epithelial proximity ----------
epi_idx <- all_data_meta_pos$Classification1 == "Malignant/epithelial cells"
res_pa_epi_ls <- res_pa_ls %>% map(~ .x[, epi_idx])
all_data_epi  <- all_data_meta_pos[epi_idx, ]


# distance diff
diff_idx_ls <- 
  seq_along(res_pa_epi_ls) %>%
  rev() %>%
  { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>%
  apply(1, c, simplify = FALSE)

diff_idx_ls


# proximity 
# concentric circles -> concentric annulus
res_pa_epi_diff_ls <- vector(length = 5, mode = "list")
for(idx in diff_idx_ls){
  res_pa_epi_diff_ls[[ idx[1] ]] <- res_pa_epi_ls[[ idx[1] ]] - res_pa_epi_ls[[ idx[2] ]]
}
res_pa_epi_diff_ls[[1]] <- res_pa_epi_ls[[1]]
names(res_pa_epi_diff_ls) <- names(res_pa_epi_ls)


# compute median grouped by slide_id
res_pa_epi_count <-
  res_pa_epi_diff_ls %>% 
  map(~ apply(.x, 1, function(x) tapply(x, all_data_epi$Slide_ID, median)))  


res_pa_epi_count <-
  res_pa_epi_count %>% 
  imap(function(x, y){
    x <- x %>%
      as.data.frame(check.names = FALSE) %>%
      rownames_to_column('Slide_ID')
    x$distance <- y
    x
  }) %>% 
  bind_rows() 


# add Group
all_data_group_epi <- all_data_epi[c('Group', 'Slide_ID')] %>% distinct()


res_pa_epi_count_df <-
  res_pa_epi_count %>%
  pivot_longer(
    -c(Slide_ID, distance), 
    names_to  = "Cell", 
    values_to = "count"
  ) %>%
  mutate(
    distance = str_c(distance, "um"),
    distance = distance %>% reorder(., as.numeric(str_remove(., "um$")))
  ) %>%
  left_join(
    all_data_group_epi, 
    by = "Slide_ID"
  ) 


# area 
area_ls <-
  res_pa_epi_count_df$distance %>% 
  unique %>% 
  str_remove("um") %>% 
  as.numeric() %>%
  sort() %>%
  map2(c(0, .[-length(.)]), ., function(x, y){
    print(c(x, y))
    pi * (y - x)^2
  })
  

## plot: facet by distance barplot ------
p_proximity_cell_count_distribution <-

  res_pa_epi_count_df %>%
  
  ggplot(aes(x = Cell, y = count, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ distance, scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median count'
  )

ggsave(
  p_proximity_cell_count_distribution, 
  height = 7,
  width  = 5 * 3,
  filename = file.path(out_dir, "2 proximity analysis - barlot - Group - combined.pdf")
)
 

# remove Group P
p_proximity_cell_count_distribution_omit_GroupP <-
  
  res_pa_epi_count_df %>%
  
  dplyr::filter(Group != 'P') %>%
  
  ggplot(aes(x = Cell, y = count, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ distance, scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median count'
  )

ggsave(
  p_proximity_cell_count_distribution_omit_GroupP, 
  height = 7,
  width  = 5 * 3,
  filename = file.path(out_dir, "2 proximity analysis - barlot - Group - remove GroupP - combined.pdf")
)

## plot: facet by celltype barplot ------
# 20250420 add
p_proximity_cell_count_distribution_by_celltype <-
  
  res_pa_epi_count_df %>%
  
  ggplot(aes(x = distance, y = count, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Cell, scales = "free") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median count'
  )

ggsave(
  p_proximity_cell_count_distribution_by_celltype, 
  height = 4 * 4,
  width  = 6 * 4,
  filename = file.path(out_dir, "2 proximity analysis - barlot - Group - combined - facet by celltype.pdf")
)


# remove Group P
p_proximity_cell_count_distribution_omit_GroupP_by_celltype <-
  
  res_pa_epi_count_df %>%
  
  dplyr::filter(Group != 'P') %>%
  
  ggplot(aes(x = distance, y = count, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Cell, scales = "free") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median count'
  )

ggsave(
  p_proximity_cell_count_distribution_omit_GroupP_by_celltype, 
  height = 4 * 4,
  width  = 6 * 4,
  filename = file.path(out_dir, "2 proximity analysis - barlot - Group - remove GroupP - combined - facet by celltype.pdf")
)

## smooth line plot ------------------
# res_pa_epi_count_raw <-
#   res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") %>%
#   mutate(distance = as.numeric(distance)) 
# 
# res_pa_epi_count_tidy <-
#   1:nrow(res_pa_epi_count_raw) %>%  
#   map(function(x){
#     count <- res_pa_epi_count_raw$count[x] 
#     if(count == 0) return(NULL)
#     xx <- 
#       res_pa_epi_count_raw[x, ] %>%
#       unlist(use.names = F) %>%
#       rep(count) %>%
#       matrix(nrow = count, byrow = T) %>%
#       as.data.frame()
#     
#     xx
#   }) %>%
#   purrr::discard(is.null) %>%
#   bind_rows() 
# 
# colnames(res_pa_epi_count_tidy) <- colnames(res_pa_epi_count_raw)

# head(res_pa_epi_count_tidy, 10)


p_proximity_smooth <-
  res_pa_epi_count_df %>%
  group_by(Cell, distance, Group) %>%
  summarise(count = median(count)) %>%
  
  mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%

  ggplot(aes(x = distance, y = count, color = Group)) +
  # geom_density(alpha = 0.6) +
  geom_smooth(method = 'loess', se = TRUE) +
  facet_wrap(~ Cell, scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = grp_colors) +
  labs(
    x = 'Distance(um)',
    y = 'Median count'
  )

# p_proximity_smooth


ggsave(
  p_proximity_smooth, 
  height = 3 * 3,
  width  = 4 * 3,
  filename = file.path(out_dir, "2 proximity analysis - smooth lineplot - Group - combined.pdf")
)

# p_proximity_smooth - remove GroupP
p_proximity_smooth2 <-
  res_pa_epi_count_df %>%
  dplyr::filter(Group != 'P') %>%
  group_by(Cell, distance, Group) %>%
  summarise(count = median(count)) %>%
  
  mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
  
  ggplot(aes(x = distance, y = count, color = Group)) +
  # geom_density(alpha = 0.6) +
  geom_smooth(method = 'loess') +
  facet_wrap(~ Cell, scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = grp_colors) +
  labs(
    x = 'Distance(um)',
    y = 'Median count'
  )

ggsave(
  p_proximity_smooth2, 
  height = 3 * 3,
  width  = 4 * 3,
  filename = file.path(out_dir, "2 proximity analysis - smooth lineplot - Group -remove GroupP - combined.pdf")
)


## radar plot ----------
p_proximity_radar <-
  res_pa_epi_count %>%
  pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
  left_join(all_data_group_epi, by = "Slide_ID")  %>%
  mutate(distance = as.numeric(distance)) %>%
  # mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
  group_by(Group, Cell, distance) %>%
  summarise(count = median(count)) %>%
  
  group_by(Cell) %>%
  mutate(count = count/max(count)) %>%
  # 
  ggplot(aes(x = Cell, y = distance, size = count, fill = Group)) +
  geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
  geom_hline(yintercept = c(10, 50, 100, 150, 200), color = "gray") +
  geom_point(
    shape = 21,
    stroke = 1.2,
    color = "gray30",
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0, #10,
      dodge.width = 0.75
    )) +
  geom_text(
    data = data.frame(x = 0, y = c(10, 50, 100, 150, 200)),
    aes(x = x, y = y, label = paste0(y, "um")),
    inherit.aes = FALSE,
    vjust = -1.
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(1,1,1,6, 'mm')
  ) +
  scale_size(range = c(3, 8)) +
  scale_y_continuous(limits = c(-100, 220)) +
  scale_fill_manual(values = grp_colors)+
  coord_polar(theta = 'x', clip = 'off') +
  labs(x = NULL, y = NULL)

# p_proximity_radar

ggsave(
  p_proximity_radar, 
  height = 8,
  width  = 12,
  filename = file.path(out_dir, "2 proximity analysis - radar plot - Group.pdf")
)


# 2024.9.29 IDO1+ Malignant cells proximity analysis ------------

norm_counts <- all_data@assays$RNA$data
IDO1_Malignant_cells_idx  <- norm_counts['IDO1', ] > 0 & all_data$final_annotation == "Malignant cells" 

all_data_meta_pos_ido <- all_data_meta_pos
all_data_meta_pos_ido$Classification1[IDO1_Malignant_cells_idx] <- "IDO1 Malignant cells" 

res_pa_ido_200 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_ido,
    query = all_data_meta_pos_ido,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )


res_pa_ido_200_freq <- sweep(res_pa_ido_200, 2, colSums(res_pa_ido_200), "/")

res_pa_ido_epi <- res_pa_ido_200[, IDO1_Malignant_cells_idx]
res_pa_ido_epi_freq <- res_pa_ido_200_freq[, IDO1_Malignant_cells_idx]
all_data_ido_epi  <- all_data_meta_pos[IDO1_Malignant_cells_idx, ]


res_pa_ido_epi_count <-
  apply(res_pa_ido_epi, 1, function(x) tapply(x, all_data_ido_epi$Slide_ID, median) )

res_pa_ido_epi_freq <-
  apply(res_pa_ido_epi_freq, 1, function(x) tapply(x, all_data_ido_epi$Slide_ID, median) )

# add Group

res_pa_ido_epi_count_df <-
  res_pa_ido_epi_count %>%
  as.data.frame() %>%
  rownames_to_column("Slide_ID") %>%
  pivot_longer(
    -Slide_ID, 
    names_to  = "Cell", 
    values_to = "count"
  ) %>%
  left_join(
    distinct(all_data_ido_epi[c('Group', 'Slide_ID')]), 
    by = "Slide_ID"
  ) 

res_pa_ido_epi_freq_df <-
  res_pa_ido_epi_freq %>%
  as.data.frame() %>%
  rownames_to_column("Slide_ID") %>%
  pivot_longer(
    -Slide_ID, 
    names_to  = "Cell", 
    values_to = "freq"
  ) %>%
  left_join(
    distinct(all_data_ido_epi[c('Group', 'Slide_ID')]), 
    by = "Slide_ID"
  ) 

# plot
p_pa_ido_epi <-
  res_pa_ido_epi_count_df %>%
  ggplot(aes(x = Cell, y = count, fill = Group)) +
  geom_boxplot() +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median count'
  )

p_pa_ido_epi_freq <-
  res_pa_ido_epi_freq_df %>%
  ggplot(aes(x = Cell, y = freq, fill = Group)) +
  geom_boxplot() +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = grp_colors) +
  labs(
    x = NULL, # 'Distance(um)',
    y = 'Median freq'
  )

ggsave(
  p_pa_ido_epi, 
  height = 6,
  width  = 9,
  filename = file.path(out_dir, "2 proximity analysis - boxplot plot - IDO cells.pdf")
)

ggsave(
  p_pa_ido_epi_freq, 
  height = 6,
  width  = 9,
  filename = file.path(out_dir, "2 proximity analysis - boxplot plot - IDO cells - freq.pdf")
)
